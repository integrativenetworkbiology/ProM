

#############################
#' Prepare input files for ProM
#'
#' This function takes the input bulk expression and scRNAseq used for reference. It renormalizes the bulk and sc expression using common genes,generates the mean and variance reference matrix from sc data and computes the marker genes for each cell cluster.
#'
#' @param bulk Gene by sample TPM/count matrix for deconvolution
#' @param sc Singlecell data in Seurat format
#' @param outfile Rfile used to store R data output from this function 
#' @param cellN downsample cellN to speed up calculation of markers
#' @param majorlabel The column in Seurat object sc used to store major cluster annotation
#' @param minorlabel The column in Seurat object sc used to store minor cluster annotation
#' @param samplelabel The column in Seurat object sc used to store sample annotation
#' @param markerN The maxium number of markers per cell cluster
#' @return if no outfile supplied, it will return a list of four objects. majorsc.basis: reference mean and variance matrix for major cluster.minorsc.basis: reference mean and variance matrix for minor cluster.bulk:the bulk matrix.markers: cell markers for each major and minor cell clusters.
#' @export
prepareInput<-function(bulk,sc,outfile,cellN=2000,majorlabel="majorlabel",minorlabel="minorlabel",samplelabel="sampleID",markerN=200)
  {
    
    commongenes<-intersect(rownames(bulk),rownames(sc))
    print(paste("There are", length(commongenes),"common genes:renormaliation"))
    sc<-sc[commongenes,]
    bulk<-as.matrix(bulk[commongenes,])
 ####renormalize using common genes 
    sc<-Seurat::NormalizeData(sc)   
    bulk<-sweep(bulk,2,10^6/colSums(bulk),"*")

############obtain markers
    print("calculate major markers")
Seurat::Idents(sc)<-sc@meta.data[,majorlabel]
markers<-Seurat::FindAllMarkers(sc,assay="RNA",only.pos=TRUE,max.cells.per.ident=cellN)
markers<-markers[-c(grep("^RPL|^RPS",markers$gene)),]

    print("calculate minor markers")
minormarkerL<-vector("list",length(unique(sc@meta.data[,majorlabel])))
names(minormarkerL)<-unique(sc@meta.data[,majorlabel])
for(i in unique(sc@meta.data[,majorlabel]))
  {
    subset<-sc[,sc@meta.data[,majorlabel]==i]
    if(length(unique(subset@meta.data[,minorlabel]))>1)
      {
        print(i)
        Seurat::Idents(subset)<-subset@meta.data[,minorlabel]
        minormarkerL[[i]]<-Seurat::FindAllMarkers(subset,assay="RNA",only.pos=TRUE,max.cells.per.ident=cellN)
      }
  }

minormarkerL<-lapply(minormarkerL,function(m){if(length(grep("^RPL|^RPS",m$gene))>0)m<-m[-c(grep("^RPL|^RPS",m$gene)),];m;})

markers<-rbind(markers,do.call("rbind",minormarkerL))
markers<-gettopmarkers(markers,topN=markerN)
    
    print("calculate mean and variance of reference profile") 
 majorsc.basis<-getMeanVar(sc, meanMethod="poolmean", celltype=majorlabel, sample=samplelabel, verbose = TRUE)
  minorsc.basis<-getMeanVar(sc, meanMethod="poolmean", celltype=minorlabel, sample=samplelabel, verbose = TRUE)

######fill in some NA values in SD of minorsc.basis (subtype only exists in one sample)
  NAminor<-colnames(minorsc.basis[[2]])[colSums(is.na(minorsc.basis[[2]]))>0]
  if(length(NAminor)>0)
  {
  map<-unique(sc@meta.data[,c(majorlabel,minorlabel)])
  rownames(map)<-map[,2]
  for(i in NAminor)minorsc.basis[[2]][which(is.na(minorsc.basis[[2]][,i])),i]<-majorsc.basis[[2]][which(is.na(minorsc.basis[[2]][,i])),map[i,1]]
  }
 ##### 
    if(!missing(outfile))save(majorsc.basis,minorsc.basis,bulk,markers,file=outfile)else{return(list(majorsc.basis=majorsc.basis,minorsc.basis=minorsc.basis,bulk=bulk,markers=markers))}
    
}


#######################################################################
##################get the cell frequency distribution from single cell data
#############################
#' Calculate cell frequency per sample from Seurat Object
#'
#' This function takes Seurat Object as input, and output cell frequency per sample.
#'
#' @param sc Singlecell data in Seurat format
#' @param clusterlabel The column in Seurat object sc used to store cell cluster annotation
#' @param samplelabel The column in Seurat object sc used to store sample annotation
#' @return Cell frequency Matrix
#' @export
getcellfreqM<-function(sc,clusterlabel="minorlabel",samplelabel="sampleID")
  {
cellfreqM<-table(sc@meta.data[,clusterlabel],sc@meta.data[,samplelabel])
sweep(cellfreqM,2,colSums(cellfreqM),"/")
}




#######################################################################################
#############still slightly different from musci_basis in that we used the normalized data from single cell to calculate mean and var, while they use raw data counts.
getMeanVar<-function (sc, meanMethod=c("samplemean","poolmean"),selectmarkers = NULL, celltype, sample, verbose = TRUE,minN=0) 
{
   if (!is.null(selectmarkers)) {
     ids <- intersect(unlist(selectmarkers), rownames(sc))
     sc<-sc[ids,]
        }
  
    clusters <- as.character(sc@meta.data[,celltype])
   samples <- as.character(sc@meta.data[, sample])
      
   x<-as.matrix(exp(Seurat::GetAssayData(sc,slot="data"))-1)*100
    if(meanMethod=="samplemean")
      {
        M <- sapply(unique(clusters), function(ct) {
          ss<-table(samples[clusters%in%ct])
     ss<-names(ss)[ss>minN]
          y<-sapply(ss, function(sid)rowMeans(x[, clusters %in% ct & samples %in% sid, drop = FALSE]))
          apply(y,1,mean,na.rm=TRUE)
        })
      }
   
     if(meanMethod=="poolmean")
      {
        M <- sapply(unique(clusters), function(ct) rowMeans(x[, clusters %in% ct, drop = FALSE]))
      }
    if (verbose) {
        message("Creating Relative Abudance Matrix...")
    }
    
   Sigma <- sapply(unique(clusters), function(ct) {
     ss<-table(samples[clusters%in%ct])
     ss<-names(ss)[ss>minN]
          y<-sapply(ss, function(sid)rowMeans(x[, clusters %in% ct & samples %in% sid, drop = FALSE]))
          apply(y,1,var,na.rm=TRUE)
        })                     
        
        if (verbose) {
            message("Creating Variance Matrix...")
        }
    
        return(list(M = M,Sigma = Sigma))
    
 }


##############################
gettopmarkers<-function(markers,topN=200)
  {
    do.call("rbind",lapply(unique(markers$cluster),function(majori){
     m<-subset(markers,cluster==majori)
     m[1:min(nrow(m),topN),]
   }))
  }
