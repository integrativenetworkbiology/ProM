#########################################################################################
#############################
#' Wrapper function for ProM deconvolution
#'
#' This function takes bulk expression profile, sc reference profile, and marker genes to estimate cell frequency for each sample in the bulk expression.
#'
#' @param bulk Gene by sample TPM matrix, output from prepareInput() function
#' @param sc.basis Reference mean and variance gene profile for each cell type, output from prepareInput() function
#' @param markers Cell markers for each minor and major cell clusters, output from prepareInput() function
#' @param priorcellfreqM Prior cell frequency distribution, can be calculated using getcellfreqM() function
#' @param method One of the three ProM stratigies: 1)using both reference profile and cell marker based deconvolution(combined); 2)using only reference profile based deconvolution(profileOnly);3)using only cell marker based deconvolution(markerOnly)
#' @param usebulkSigma Indicating whether gene weight/variance that captures the difference between sc and bulk are considered
#' @param Zcutoff Cell clusters with Zfrequency (percentage of samples without this cell cluster) larger than Zcutoff need to use linear regression with intercept in marker-based deconvolution.
#' @param topmarkerN Number of top-ranking markers per cell cluster selected in marker-based deconvolution.
#' @param CCcutoff Correlation Coefficient cutoff used to select top-ranking markers per cell cluster in marker-based deconvolution.
#' @param transformation Whether bulk gene expression need to transform to reduce the possible platform difference between sc and bulk RNAseq. 1)none: no transformation; 2)profile:transformation only applies to reference profile-based deconvolution; 3)marker: transformation only applied to marker-based deconvolution; 4)both: transformation applied to both.
#' @return Estimated cell frequency matrix
#' @export
ProM<-function(bulk,sc.basis,markers,priorcellfreqM,method="combined",usebulkSigma=TRUE,Zcutoff=0.5,transformation="none",CCcutoff,topmarkerN=20,markerMethod="simulation")
  {
    ##########check the inputs
    if(any(sort(rownames(bulk))!=sort(rownames(sc.basis[[1]]))))stop("bulk and sc.basis have different genes. Please run prepareInput() to prepare the inputs")
    if(!transformation %in%c("none","profile","marker","both"))stop("transformation mode is not recognizable")
    if(!method %in% c("combined","profileOnly","markerOnly")) stop("method is not recognizable")
    if(missing(markers))stop("markers are needed")
    if(transformation!="none"&missing(priorcellfreqM)) print("priorcellfreqM is needed for transformation")
    if(method!="profileOnly"&missing(priorcellfreqM)&markerMethod=="simulation") print("priorcellfreqM is needed for marker-based deconvolution from simulation")
    if(!missing(priorcellfreqM)&any(colSums(priorcellfreqM>1)>0)) print("Values in priorcellfreqM can not be larger than 1")
    if(!missing(priorcellfreqM))priorcellfreqM[priorcellfreqM==1]<-0.99 ###cell frequency can't be 1, otherwize beta estimation will be error
    
      
      
    ###############run ProM with and without markerintercept for each cell type based on NZpercent
    finalEstimate<-promDeconvolution(bulk,sc.basis,markers,priorcellfreqM,method,usebulkSigma,markerintercept=FALSE,transformation,CCcutoff,topmarkerN=topmarkerN,markerMethod=markerMethod)
    
      Zpercent<-rowMeans(priorcellfreqM==0)
      ct<-names(Zpercent)[Zpercent>Zcutoff]
      if(length(ct)>0)
        {
      print(paste("The following celltype needs ProM with markerintercept",paste(ct,collapse=" ")))
      t<-promDeconvolution(bulk,sc.basis,markers,priorcellfreqM,method,usebulkSigma,markerintercept=TRUE,transformation,CCcutoff,topmarkerN=topmarkerN,markerMethod=markerMethod)
      finalEstimate[ct,]<-t[ct,]
      finalEstimate<-sweep(finalEstimate,2,colSums(finalEstimate),"/")
    }
    
    return(finalEstimate)
  }


   








###################################################################################################
promDeconvolution<-function(sbulk,minorsc.basis,minormarkers,calibrateminorcellfreqM, method=c("combined","profileOnly","markerOnly"),usebulkSigma=TRUE,markerintercept=FALSE,transformation=c("none","profile","marker","both"),CCcutoff,iter.max=100,simSeed=1,simN=1000,topmarkerN=20,markerMethod=c("simulation","simple"),normCellFreq=TRUE)
  {

############    
    sbulk<-sbulk[rownames(minorsc.basis[[1]]),]
    
    if(method %in% c("profileOnly","combined"))
      {
        print("estimation from profile")
        if(transformation%in%c("profile","both"))
      {
        print("transformation")
        minorreconstructM<-minorsc.basis[["M"]]%*%calibrateminorcellfreqM[colnames(minorsc.basis[["M"]]),]
        sbulk<-batchcorrection(minorreconstructM,sbulk)
      }

        profileEstimate<-estimatefromprofile(minorsc.basis,sbulk,subsetgenes=unique(minormarkers$gene),iter.max=iter.max)
        if(usebulkSigma)
          {
            print("estimation from profile with bulkSigma")
            bulksigma<-getbulksigma(minorsc.basis,profileEstimate,sbulk)
            profileEstimate<-estimatefromprofile(minorsc.basis,sbulk,bulkSigma=bulksigma,subsetgenes=unique(minormarkers$gene),iter.max=iter.max)
          }
      }
    
    if(method %in% c("markerOnly","combined"))
      {
        print("estimation from marker")
        print(paste("marker method:",markerMethod))
         if(transformation%in%c("marker","both"))
      {
        print("transformation")
        minorreconstructM<-minorsc.basis[["M"]]%*%calibrateminorcellfreqM[colnames(minorsc.basis[["M"]]),]
        sbulk<-batchcorrection(minorreconstructM,sbulk)
      }
        
        print("simulation without bulkSigma")
        simData<-simulateBulk(calibrateminorcellfreqM,minorsc.basis,N=simN,seed=simSeed)
       if(markerMethod=="simulation")markerEstimate<-predictfromsimulation_mean(simData$sbulk,simData$sminorcellfreqM,minormarkers,sbulk,cutoff=CCcutoff,topmarkerN=topmarkerN,intercept=markerintercept)
       if(markerMethod=="simple")  markerEstimate<-predictfrommarker_simple(simData$sbulk,simData$sminorcellfreqM,minormarkers,sbulk,topmarkerN=topmarkerN,intercept=markerintercept)
         if(usebulkSigma)
          {
            print("simulation with bulkSigma")
            bulksigma<-getbulksigma(minorsc.basis,markerEstimate,sbulk)
           simData<-simulateBulk(calibrateminorcellfreqM,minorsc.basis,N=simN,bulksigma=bulksigma,seed=simSeed)
         if(markerMethod=="simulation") markerEstimate<-predictfromsimulation_mean(simData$sbulk,simData$sminorcellfreqM,minormarkers,sbulk,cutoff=CCcutoff,topmarkerN=topmarkerN,intercept=markerintercept)
          if(markerMethod=="simple")  markerEstimate<-predictfrommarker_simple(simData$sbulk,simData$sminorcellfreqM,minormarkers,sbulk,topmarkerN=topmarkerN,intercept=markerintercept)  
          }
      }
        

    if(method %in% "combined")
      {
        print("combine estimation from profile and marker")
        finalEstimate<-(markerEstimate+profileEstimate[rownames(markerEstimate),])/2
      }
    if(method=="profileOnly")finalEstimate<-profileEstimate
    if(method=="markerOnly")finalEstimate<-markerEstimate
   
    norm<-function(transpM)sweep(transpM,2,colSums(transpM),"/")
    if(normCellFreq) {
      print("normalizing cell frequency")
      finalEstimate<-norm(finalEstimate)
      #### in some extreme cases, especially for marker-based deconvolution with interception, the cell freq could be all zero, hence NA is generated after normalization
      finalEstimate[is.na(finalEstimate)]<-0
    }
    return(finalEstimate)
  }




######################
batchcorrection <- function(Y.train, X.pred,rescale=TRUE,log=TRUE,percent=1) {
  
  if(log){
    Y.train<-log2(Y.train+1);
    X.pred<-log2(X.pred+1);
  }
  Y.center <- apply(Y.train,1,mean)
    n <- ncol(Y.train)
   shrink.scale <- apply(Y.train,1,sd)
  X.center<-apply(X.pred,1,mean)
  X.scale<-apply(X.pred,1,sd)
  Y.center<-(Y.center-X.center)*percent+X.center
  shrink.scale<-(shrink.scale-X.scale)*percent+X.scale
  
  if(!rescale)shrink.scale<-apply(X.pred,1,sd)

  X.pred.scaled<-t(apply(X.pred,1,scale))
  X.pred.scaled[is.na(X.pred.scaled)]<-0
   Y.pred<-sweep(X.pred.scaled,1, shrink.scale,"*")
   Y.pred<-sweep(Y.pred,1,Y.center,"+")
  colnames(Y.pred)<-colnames(X.pred)
  if(log)Y.pred<-2^Y.pred-1
  Y.pred[Y.pred<0]<-0
  return(Y.pred)
}



############################################################################################################################################################################################

estimatefromprofile<-function(sc.basis,bulk,bulkSigma,subsetgenes,iter.max=1000,nu = 1e-04,eps = 0.01,useresidual=TRUE,usescSigma=TRUE)
  {
    if(missing(subsetgenes))subsetgenes<-rownames(sc.basis[["M"]])
    X<-sc.basis[["M"]][subsetgenes,]
    scSigma<-sc.basis[["Sigma"]][subsetgenes,]
    if(!usescSigma)scSigma[,]<-0
    if(missing(bulkSigma)){
      bulkSigma<-rep(0,length(subsetgenes))
      names(bulkSigma)<-subsetgenes
    }else{
      bulkSigma<-bulkSigma[subsetgenes]
    }
    
    transpM<-sapply(1:ncol(bulk),function(i){
      #print(i)
       Y<-bulk[subsetgenes,i]
      weightnnls(Y, X, scSigma,bulkSigma=bulkSigma,iter.max=iter.max, nu = nu, eps = eps,residual=useresidual)
            })
    rownames(transpM)<-colnames(sc.basis[["M"]])
    colnames(transpM)<-colnames(bulk)
    transpM<-sweep(transpM,2,colSums(transpM),"/")
    return(transpM)
  }
                             
#########
                             
####################

 weightnnls<-function (Y, X, scSigma,bulkSigma=0,iter.max=1000, nu = 1e-04, eps = 0.01,trace=FALSE,residual=TRUE) 
{
  if(trace)pM<-c()
  lm.D = nnls::nnls(X, Y)
  r = resid(lm.D)
  p=lm.D$x
  for (iter in 1:iter.max) {
    if(!residual)r<-0
      weight.gene = 1/(nu +r^2+ bulkSigma + colSums(p^2 * t(scSigma)))
        Y.weight = Y * sqrt(weight.gene)
        D.weight = sweep(X, 1, sqrt(weight.gene), "*")
        lm.D.weight = nnls::nnls(D.weight, Y.weight)
         p.weight = lm.D.weight$x     
        if (sum(abs(p.weight - p)) < eps) break; 
      p = p.weight
      r = resid(lm.D.weight)
 if(trace)pM<-rbind(pM,p)
    }
  #print(paste("iterations",iter))
  if(trace) return(pM)
  p.weight;
}
 
##################
getbulksigma<-function(minorsc.basis,minorcellfreqM,bulk,transformation=FALSE)
  {
    minorcellfreqM<-minorcellfreqM[colnames(minorsc.basis[["M"]]),]
    minorreconstructM<-minorsc.basis[["M"]]%*%minorcellfreqM
    ss<-intersect(colnames(bulk),colnames(minorcellfreqM))
    print(paste(length(ss),"paired samples used to get bulksigma"))
    bulk<-bulk[rownames(minorreconstructM),]
    if(transformation)
      {
        minorreconstructM[,ss]<-batchcorrection(bulk[,ss],minorreconstructM[,ss])
      }   
    totalSigma<-apply((bulk[,ss]-minorreconstructM[,ss])^2,1,mean)
    return(totalSigma)
}


#############
predictfromsimulation_mean<-function(trainbulk,traincellfreqM,markers,testbulk,log=TRUE,cutoff,topmarkerN=20,intercept=FALSE)
  {
    usecutoff<-FALSE
    if(!missing(cutoff))usecutoff<-TRUE;
    precellfreqM<-t(sapply(rownames(traincellfreqM),function(minori){
      #print(minori)
      m<-subset(markers,cluster==minori)$gene
      cc<-cor(t(trainbulk[m,]),traincellfreqM[minori,],method="pearson");
        ################using old CC cutoff to select genes, default cutoff =0.3
     if(usecutoff){
      r<-rownames(cc)[cc>cutoff]
      curcutoff<-cutoff
      while(length(r)<20&curcutoff>=(-1))
        {
          curcutoff<-curcutoff-0.1
          #print(curcutoff)
          r<-rownames(cc)[cc>curcutoff]
        }
         }else{
         ######using topmarkerN to select genes
         r<-rownames(cc)[order(as.vector(cc),decreasing=TRUE)]
        # print(cc[r,])
         r<-r[1:min(topmarkerN,length(r))]        
     }
           
   
      if(log){ trainx<-exp(colMeans(log(trainbulk[r,,drop=FALSE]+1)))-1;
               testx<-exp(colMeans(log(testbulk[r,,drop=FALSE]+1)))-1}
      if(!log){ trainx<-colMeans(trainbulk[r,,drop=FALSE]);
               testx<-colMeans(testbulk[r,,drop=FALSE])}
      if(!intercept)m<-lm(cellfreq ~ exp-1, data=data.frame(cellfreq=traincellfreqM[minori,],exp=trainx))
      if(intercept)m<-lm(cellfreq ~ exp, data=data.frame(cellfreq=traincellfreqM[minori,],exp=trainx))
        
     predict(m,new=data.frame(cellfreq=NA,exp=testx)) 
    }))
    colnames(precellfreqM)<-colnames(testbulk)
    precellfreqM[precellfreqM<0]<-0
    return(precellfreqM);
  }



###########################################################################################################simulate data

simulateBulk<-function(minorcellfreqM,minorsc.basis,N=1000,outfile,bulksigma,targetbulk,seed=1,seed2=1)
  {
posRate<-apply(minorcellfreqM>0,1,mean)
fitR<-apply(minorcellfreqM,1,function(x){t<-x[x>0];if(length(t)==1)t<-c(t+0.001,t);fitdistrplus::fitdist(t,"beta")$estimate})

simulate<-function(N,posrate,shape1,shape2)
  { 
    s<-rbeta(N,shape1,shape2)
    s[runif(N)>posrate]<-0
    s
  }
set.seed(seed)
sminorcellfreqM<-sapply(rownames(minorcellfreqM),function(i)simulate(N,posRate[i],fitR[1,i],fitR[2,i]))
sminorcellfreqM<-sweep(sminorcellfreqM,1,rowSums(sminorcellfreqM),"/")
sminorcellfreqM<-t(sminorcellfreqM)


sprofileL<-lapply(1:ncol(minorsc.basis[[1]]),function(j)
                  {#print(j);
                    sapply(1:nrow(minorsc.basis[[1]]),function(i)rnorm(N,minorsc.basis[[1]][i,j],sqrt(minorsc.basis[[2]][i,j])))})
                names(sprofileL)<-colnames(minorsc.basis[[1]])
                for(j in names(sprofileL))
                  {#print(j)
                   sprofileL[[j]]<-sweep(sprofileL[[j]],1,sminorcellfreqM[j,],"*")
                   sprofileL[[j]][sprofileL[[j]]<0]<-0
                  }
sbulk<-t(Reduce("+",sprofileL))
rownames(sbulk)<-rownames(minorsc.basis[[1]])
colnames(sbulk)<-paste0("S",1:N)

if(!missing(bulksigma)){
  set.seed(seed2)
  bulksigma<-bulksigma[rownames(sbulk)];
  sbulk<-sbulk+t(sapply(1:length(bulksigma),function(i)rnorm(N,0,sqrt(bulksigma[i]))))
}
sbulk[sbulk<0]<-0  

if(!missing(targetbulk)) sbulk<-batchcorrection(targetbulk,sbulk)
sbulk<-sweep(sbulk,2,colSums(sbulk),"/")*10^6

if(!missing(outfile)){
  save(sbulk,sminorcellfreqM,file=outfile)}else{
    return(list(sbulk=sbulk,sminorcellfreqM=sminorcellfreqM))
  }

}


#############
predictfrommarker_simple<-function(trainbulk,traincellfreqM,markers,testbulk,log=TRUE,topmarkerN=20,intercept=FALSE)
  {
    precellfreqM<-t(sapply(rownames(traincellfreqM),function(minori){
      #print(minori)
      m<-subset(markers,cluster==minori)$gene
      r<-m[1:min(topmarkerN,length(m))]
      if(log){ trainx<-exp(colMeans(log(trainbulk[r,,drop=FALSE]+1)))-1;
               testx<-exp(colMeans(log(testbulk[r,,drop=FALSE]+1)))-1}
      if(!log){ trainx<-colMeans(trainbulk[r,,drop=FALSE]);
               testx<-colMeans(testbulk[r,,drop=FALSE])}
      if(!intercept)m<-lm(cellfreq ~ exp-1, data=data.frame(cellfreq=traincellfreqM[minori,],exp=trainx))
      if(intercept)m<-lm(cellfreq ~ exp, data=data.frame(cellfreq=traincellfreqM[minori,],exp=trainx))
        
     predict(m,new=data.frame(cellfreq=NA,exp=testx)) 
    }))
    colnames(precellfreqM)<-colnames(testbulk)
    precellfreqM[precellfreqM<0]<-0
    return(precellfreqM);
  }
