library("ProM")
library("parallel")
library("MASS")
library("dplyr")
library("ggplot2")
library("gridExtra")
library("Seurat")
library("ggrepel")
library("foreach")
install.packages("doMC")
library("doMC")
registerDoMC(30)
inputDir<-"/scratch/data/GSE210347/"
simDir<-"/scratch/ReconstructionCV_pancancer/"
library("ProM")
library("parallel")
library("MASS")
library("dplyr")
library("ggplot2")
library("gridExtra")
library("Seurat")
library("ggrepel")
library("foreach")
install.packages("doMC")
library("doMC")
registerDoMC(30)
inputDir<-"/scratch/data/GSE210347/"
simDir<-"/scratch/ReconstructionCV_pancancer/"
dir.create(simDir)


###################data
sc<-readRDS("/data/GSE210347/GSE210347_counts.Rds")
meta<-read.table("/data/GSE210347/GSE210347_meta.txt",header=TRUE,row.names=1)
sc<-CreateSeuratObject(sc,meta.data=meta)
sc<-sc[,sc$cluster!="undefined"]          
t<-mclapply(unique(sc$tissue),function(tis)
{
  print(tis)
  Obj<-sc[,sc$tissue==tis]
  clusterN<-table(Obj$cluster)
  Obj<-Obj[,Obj$cluster%in%names(clusterN)[clusterN>100]]
  ###filter non expressing genes
  Obj<-Obj[rowSums(GetAssayData(Obj)>0)>10,]
  ###
  Obj<-NormalizeData(Obj)
  save(Obj,file=paste0(simDir,"Obj_",tis,".rda"))
},mc.cores=length(unique(sc$tissue)))

meta<-read.table("/data/GSE210347/GSE210347_meta.txt",header=TRUE,row.names=1)
for(tis in setdiff(unique(meta$tissue),c("liver","bladder","colorectal")))
{
  print(tis)
 load(file=paste0(simDir,"Obj_",tis,".rda"))
 totalbulk<-sapply(unique(Obj$SampleID),function(ss)rowSums(Obj[,Obj$SampleID==ss]))
 totalbulk<-sweep(totalbulk,2,colSums(totalbulk),"/")*10^6
 sampleNames<-unique(Obj$SampleID)
 totalObj<-Obj
 
t=foreach(index = 1:length(sampleNames))%dopar%
{
  print(index)
#############
sample2remove<-sampleNames[index]
print(paste(sample2remove,"is LOO sample"))
Obj<-totalObj[,totalObj$SampleID!=sample2remove]
minorcellfreqM<-table(Obj$cluster,Obj$SampleID)
minorcellfreqM<-sweep(minorcellfreqM,2,colSums(minorcellfreqM),"/")

########################################################################################
##################prepare ProM input files
prepareInput(totalbulk,Obj,outfile=paste0(simDir,"ProM_input_CV_",tis,"_",index,".rda"),majorlabel="celltype",minorlabel="cluster",samplelabel="SampleID")
a<-load(paste0(simDir,"ProM_input_CV_",tis,"_",index,".rda"))
############################################
#####################################################################run ProM
#####some clusters may have no markers(we will ignore the deconvolution of those clusters, marker-based ProM will not work then)
proM<-TRUE
if(proM)
{
cluster2remove<-setdiff(colnames(minorsc.basis[[1]]),unique(markers$cluster))
if(length(cluster2remove)>0)
  {
  minorsc.basis[[1]]<-minorsc.basis[[1]][,!colnames(minorsc.basis[[1]])%in%cluster2remove]
  minorsc.basis[[2]]<-minorsc.basis[[2]][,!colnames(minorsc.basis[[2]])%in%cluster2remove]
  minorcellfreqM<-minorcellfreqM[!rownames(minorcellfreqM)%in%cluster2remove,]
}
t1<-ProM(bulk=bulk,sc.basis=minorsc.basis,markers=markers,priorcellfreqM=minorcellfreqM)
t2<-ProM(bulk=bulk,sc.basis=minorsc.basis,markers=markers,priorcellfreqM=minorcellfreqM,method="profileOnly")
t3<-ProM(bulk=bulk,sc.basis=minorsc.basis,markers=markers,priorcellfreqM=minorcellfreqM,method="markerOnly")
finalEstimateL<-list(ProM=t1,ProM_profileOnly=t2,ProM_markerOnly=t3)
save(finalEstimateL,file=paste0(simDir,"ProMresults_topN20/finalEstimateL_",tis,"_",index,".rda"))
}
######################################################################################################
##########################################################################other deconvolution Results

###########################################
other<-TRUE
    if(other)
        {
otherL<-vector("list",6)
names(otherL)<-c("DSA","nnls","MuSiC","BisqueRNA","CIBERSORT","DeconRNASeq")
for(m in names(otherL))
{
  print(m)
  otherL[[m]]<-otherDeconvolution(totalbulk,minorsc.basis,markers,Obj,method=m,ncores=1, clusterlabel = "cluster", 
                                  samplelabel = "SampleID")
}
save(otherL,file=paste0(simDir,"otherL_",tis,"_",index,".rda"))
}



######################################################################################################################################
#####################
ted<-TRUE
if(ted)
{
library("TED")
library("Seurat")
gg<-unique(markers$gene)
print(paste(length(gg),"genes included in the deconvolution"))

ref.dat<-as.matrix(expm1(GetAssayData(Obj[gg,]))*100)

ref.dat.filtered <- cleanup.genes(ref.dat= t(ref.dat),
                                  species="hs",
                                  gene.type=c("RB","chrM","chrX","chrY"),
                                  input.type="scRNA",
                                  exp.cells=5)
meta<-Obj@meta.data[,c("celltype","cluster")]

curalpha<-1
ted<-try(run.Ted (ref.dat= ref.dat.filtered,
                  X=t(totalbulk),
                  cell.type.labels= meta$celltype,
                  cell.subtype.labels= meta$cluster,
                  tum.key="Epithelium",
                  input.type="scRNA",
                  n.cores=9,
                  first.gibbs.only=FALSE,
                  seed=1,
                  if.vst=FALSE,
                  alpha=curalpha))
    if(length(ted)==1)
{
####sometimes, there is error in hclust since all zero? we will adjust alpha to try
  while(length(ted)==1&curalpha>0.9)
    {
    curalpha<-curalpha-0.01
    print(curalpha)
  ted<-try(run.Ted (ref.dat= ref.dat.filtered,
                    X=t(totalbulk),
                    cell.type.labels= meta$celltype,
                    cell.subtype.labels= meta$cluster,
                    tum.key="Epithelium",
                    input.type="scRNA",
                    n.cores=9,
                    first.gibbs.only=FALSE,
                    seed=1,
                    if.vst=FALSE,
                    alpha=curalpha))
  }
}

ted<-t(ted$res$first.gibbs$gibbs.theta)
save(ted,file=paste0(simDir,"ted2step_",tis,"_",index,".rda"))

}
}
                   }

                   
###########################################
#############################################################################################
###########major cell to minor cell map
minor2major<-function(minorfM,cellmap,mapcol="celltype")
{
  cellmap<-cellmap[rownames(minorfM),]
  majorfM<-sapply(unique(cellmap[,mapcol]),function(x){
    colSums(minorfM[rownames(cellmap)[cellmap[,mapcol]==x],,drop=FALSE])})
  colnames(majorfM)<-unique(cellmap[,mapcol])
  t(majorfM)
}


#####
meta<-read.table("/data/GSE210347/GSE210347_meta.txt",header=TRUE,row.names=1)
cellmap<-unique(meta[,c("cluster","celltype")])
cellmap<-subset(cellmap,cluster!="undefined")
rownames(cellmap)<-cellmap[,"cluster"]

############################################################################################################
################################obtain LOO results
#######################################

ddL<-lapply( unique(meta$tissue),function(tis)
{
  print(tis)
load(file=paste0(simDir,"Obj_",tis,".rda"))
sampleNames<-unique(Obj$SampleID)
minorcellfreqM<-table(Obj$cluster,Obj$SampleID)
minorcellfreqM<-sweep(minorcellfreqM,2,colSums(minorcellfreqM),"/")
    ######################load proM CV data 
sampleN<-length(sampleNames)
finalEstimateLCV<-lapply(c("ProM","ProM_profileOnly","ProM_markerOnly"),function(methodIndex)
{
  finalEstimate<-minorcellfreqM
  finalEstimate[,]<-0
  for(index in 1:sampleN)
  {
    sample2remove<-sampleNames[index]
    load(paste0(simDir,"/ProMresults_topN20/finalEstimateL_",tis,"_",index,".rda"))
    t<-finalEstimateL[[methodIndex]][,sample2remove]
    finalEstimate[names(t),sample2remove]<-t
  }
  finalEstimate
})
names(finalEstimateLCV)<-c("ProM","ProM_profileOnly","ProM_markerOnly")
 
######################################################add Ted results
finalEstimate<-minorcellfreqM
finalEstimate[,]<-0
for(index in 1:sampleN)
{
  sample2remove<-sampleNames[index]
  load(paste0(simDir,"/ted2step_",tis,"_",index,".rda"))
  finalEstimate[rownames(ted),sample2remove]<-ted[,sample2remove]
}
finalEstimateLCV[["BayesPrism"]]<-finalEstimate

#################CV results of other methods

otherLCV<-lapply(c("DSA",'CibersortX','CibersortX_batchB','CibersortX_batchS',"DeconRNASeq","nnls","MuSiC","BisqueRNA"),function(methodIndex)
{
  finalEstimate<-minorcellfreqM
  finalEstimate[,]<-0
  for(index in 1:sampleN)
  {
    sample2remove<-sampleNames[index]
    load(paste0(simDir,"/otherL_",tis,"_",index,".rda"))
    t<-otherL[[methodIndex]][,sample2remove]
    finalEstimate[names(t),sample2remove]<-t
  }
  finalEstimate
})
names(otherLCV)<-c("DSA",'CibersortX','CibersortX_batchB','CibersortX_batchS',"DeconRNASeq","nnls","MuSiC","BisqueRNA")
finalEstimateL<-c(finalEstimateLCV,otherLCV)
########################calculate correlation and MSE
####################################################################################################################accuracy analysis Random
getCorMeanPerCell<-function(x){t<-diag(cor(t(x),t(sminorcellfreqM[rownames(x),colnames(x)]),method="pearson"));t[is.na(t)]<-0;mean(t);}
getCorAll<-function(x){t<-cor(as.vector(x),as.vector(sminorcellfreqM[rownames(x),colnames(x)]),method="pearson");t;}
getMSE<-function(x){t<-mean((as.vector(x)-as.vector(sminorcellfreqM[rownames(x),colnames(x)]))^2);t;}
dd<-c()
ss<-sampleNames
    
####minor cluster    
sminorcellfreqM<-minorcellfreqM
d<-t(sapply(finalEstimateL,function(finalEstimate)c(MeanPerCell=getCorMeanPerCell(finalEstimate[,ss]),AllCell=getCorAll(finalEstimate[,ss]),MSE=getMSE(finalEstimate[,ss]))))
dd<-data.frame(d,method=names(finalEstimateL),celltype="minorcluster")
####major cluster
sminorcellfreqM<-minor2major(minorcellfreqM,cellmap)
majorfinalEstimateL<-lapply(finalEstimateL,function(x)minor2major(x,cellmap))
d<-t(sapply(majorfinalEstimateL,function(finalEstimate)c(MeanPerCell=getCorMeanPerCell(finalEstimate[,ss]),AllCell=getCorAll(finalEstimate[,ss]),MSE=getMSE(finalEstimate[,ss]))))
dd<-rbind(dd,data.frame(d,method=names(finalEstimateL),celltype="majorcluster"))

######
colnames(dd)<-c("mean_PCC_percell","overall_PCC","MSE","method","clustertype")
dd;
})

###############################################
####################################plot accuracy     
names(ddL)<-unique(meta$tissue)

dd<-do.call("rbind",lapply(names(ddL),function(tis)cbind(ddL[[tis]],tissue=tis)))
library(reshape2)
dd<-melt(dd,id.vars=c("method","tissue","clustertype"),measure.vars=c("mean_PCC_percell","overall_PCC","MSE"),variable.name="measurement")
dd<-subset(dd,method!="CibersortX_batchS")
####remove liver and ICC (<10 samples)
dd<-subset(dd,measurement!="overall_PCC"&!tissue%in%c("liver","ICC"))
                           
pdf(paste0(simDir,"accuracy_tedstep2_top20.pdf"),width=8,height=5)
ggplot(subset(dd,measurement!="MSE"),aes(x=reorder(method,value,median),y=value))+geom_boxplot()+geom_point(aes(color=tissue))+facet_wrap(~measurement+clustertype,ncol=2)+xlab("")+ylab("")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(subset(dd,measurement=="MSE"),aes(x=reorder(method,value,median),y=value))+geom_boxplot()+geom_point(aes(color=tissue))+facet_wrap(~clustertype,ncol=2)+xlab("")+ylab("")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(subset(dd,measurement!="MSE"),aes(x=reorder(tissue,value,median),y=value))+geom_boxplot()+geom_point(aes(color=method))+facet_wrap(~measurement+clustertype,ncol=2)+xlab("")+ylab("")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggplot(subset(dd,measurement=="MSE"),aes(x=reorder(tissue,value,median),y=value))+geom_boxplot()+geom_point(aes(color=method))+facet_wrap(~clustertype,ncol=2)+xlab("")+ylab("")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()    