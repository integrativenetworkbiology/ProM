#########################################################################################
#############################
#' Wrapper function for other deconvolution methods
#'
#' This function takes bulk expression profile, sc reference profile, and marker genes to estimate cell frequency for each sample in the bulk expression.
#'
#' @param bulk Gene by sample TPM matrix, output from prepareInput() function
#' @param sc.basis Reference mean and variance gene profile for each cell type, output from prepareInput() function
#' @param markers Cell markers for each minor and major cell clusters, output from prepareInput() function
#' @param method One of the following deconvolution methods: c("CIBERSORT","DSA","DeconRNASeq","nnls","MuSiC","BisqueRNA") 
#' @return Estimated cell frequency matrix

#' @export

library("parallel")
otherDeconvolution<-function(bulk,sc.basis,markers,sc,method=c("CIBERSORT","DSA","DeconRNASeq","nnls","MuSiC","BisqueRNA"),clusterlabel="minorlabel",samplelabel="sampleID",ncores=10)
  {
    
    gg<-unique(markers$gene)
    print(paste(length(gg),"genes included in the deconvolution"))
        C = sc.basis[["M"]][gg,]
        T = bulk[gg,]
        
if(method=="CIBERSORT"){ #without QN. By default, CIBERSORT performed QN (only) on the mixture.
        
        index<-sample(1:ncores,ncol(T),replace=TRUE)
        RESULTS = do.call("rbind",mclapply(1:ncores,function(i)CIBERSORT(sig_matrix = C, mixture_file = T[,index==i,drop=FALSE], QN = FALSE),mc.cores=ncores)) 
        finalEstimate = t(RESULTS[,1:(ncol(RESULTS)-3)])
        #if there is one sample in one index, it will cause issue because of no colnames
        #print(colnames(finalEstimate))
        #print(colnames(T))
        finalEstimate=finalEstimate[,colnames(T)]

      }

if(method=="DSA")
  {
    mm<-markers[,c("gene","cluster")]
    mm<-subset(mm,cluster%in%colnames(sc.basis[[1]]))
             finalEstimate<-myDSA(T,mm,robust=TRUE)
     }
        

 if(method=="DeconRNASeq"){
       library(pcaMethods)
        finalEstimate = t(DeconRNASeq::DeconRNASeq(datasets = as.data.frame(T), signatures = as.data.frame(C),fig = FALSE)$out.all)
        colnames(finalEstimate) = colnames(T)
     }
 

if (method=="nnls"){

        RESULTS = do.call(cbind.data.frame,lapply(apply(T,2,function(x) nnls::nnls(as.matrix(C),x)), function(y) y$x))
        finalEstimate = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
        rownames(finalEstimate) <- colnames(C)

      }
    
        #############

    T.eset<-Biobase::ExpressionSet(assayData = as.matrix(bulk))

    phenoData<-data.frame(cellID=colnames(sc),cellType=sc@meta.data[,clusterlabel],SubjectName=sc@meta.data[,samplelabel])
    rownames(phenoData)<-phenoData$cellID
    C.eset <- Biobase::ExpressionSet(assayData = as.matrix(exp(GetAssayData(sc))-1)*100,phenoData = Biobase::AnnotatedDataFrame(phenoData))
    
    if(method=="BisqueRNA")finalEstimate <- try(BisqueRNA::ReferenceBasedDecomposition(T.eset, C.eset, markers=gg,use.overlap = FALSE)$bulk.props) 
    if(method=="MuSiC")
      {
        require(xbioc)
        finalEstimate  = try(t(MuSiC::music_prop(bulk.eset = T.eset, sc.eset = C.eset, clusters = 'cellType', markers = gg, samples = 'SubjectName', verbose = F)$Est.prop.weighted))
      }

    return(finalEstimate)

  }



