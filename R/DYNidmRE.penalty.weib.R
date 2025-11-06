### Code:
##' @title Penalised Illness-death model algorithm with weibull baseline risk
##' @param b  parameters not fixed
##' @param size_V number of parameters
##' @param fix0 indicators of fixed and unfixed parameters
##' @param ctime classification of subject according to their observations
##' @param N number of subjects
##' @param ve01 variables for transition 0 -->1 
##' @param ve02 variables for transition 0 -->2
##' @param ve12 variables for transition 1 -->2
##' @param dimnva01 number of variables for transition 0 -->1 
##' @param dimnva02 number of variables for transition 0 -->2
##' @param dimnva12 number of variables for transition 1 -->2
##' @param nvat01 number of variables for transition 0 -->1 
##' @param nvat02 number of variables for transition 0 -->2
##' @param nvat12 number of variables for transition 1 -->2
##' @param t0 time entry
##' @param t1 time L
##' @param t2 time R
##' @param t3 time of event/out
##' @param epsa control convergence parameter for beta 
##' @param epsb control convergence parameter for loglik
##' @param epsd control convergence for distance to minimum rdm
##' @param eps.eigen the power of convergence for eigen values of covariance matrix only
##' @param clustertype in which cluster to work
##' @param nproc number of cluster
##' @param maxiter Maximum number of iterations. The default is 200.
##' @param maxiter.pena Maximum number of iterations for penalised coefficients
##' @param troncature indicator if troncature or not
##' @param lambda01 Lambda on transition 0 --> 1
##' @param lambda02 Lambda on transition 0 --> 2
##' @param lambda12 Lambda on transition 1 --> 2
##' @param nlambda01 number of Lambda on transition 0 --> 1
##' @param nlambda02 number of Lambda on transition 0 --> 2
##' @param nlambda12 number of Lambda on transition 1 --> 2
##' @param alpha alpha on all transitions 
##' @param penalty which penalty to consider
##' @param penalty.factor which variable should be penalised
##' @param gausspoint number of points in gauss quadrature
##' @param partialH which hessian is computed for Newton-Raphson path of penalised parameters
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
#' @importFrom foreach "%do%"
#' @importFrom foreach "%dopar%"
#' @useDynLib HIDeM


DYNidmRE.penalty.weib<-function(b,fix0,size_V,
                   clustertype,epsa,epsb,epsd,eps.eigen,nproc,maxiter,maxiter.pena,
                   ctime,N,
                   ve01,ve02,ve12,dimnva01,dimnva02,dimnva12,nvat01,nvat02,nvat12,
                   t0,t1,t2,t3,troncature,
                   nlambda01,lambda01,nlambda02,lambda02,nlambda12,lambda12,
                   alpha,penalty.factor,penalty,partialH,
                   modelY,dataLongi,dataSurv,
                   Nsample,BLUP,seed,timeVar,id,formLong,
                   outcome01,outcome02,
                   outcome12,NtimePoints,
                   p01,p02,p12,
                   dimp01,dimp02,dimp12,scale.X){
  


  # create grid 3
  lambda<-expand.grid(lambda01,lambda02,lambda12)
  lambda<-unique(lambda)
  nlambda<-dim(lambda)[1]
  
  
  npm01<-ifelse(nvat01>0,sum(fix0[7:(7+nvat01-1)]==0),0)
  npm01Y<-ifelse(p01>0,sum(fix0[(7+nvat01+nvat02+nvat12):(6+nvat01+nvat02+nvat12+p01)]==0),0)
  npm02<-ifelse(nvat02>0,sum(fix0[(7+nvat01):(6+nvat01+nvat02)]==0),0)
  npm02Y<-ifelse(p02>0,sum(fix0[(7+nvat01+nvat02+nvat12+p01):(6+nvat01+nvat02+nvat12+p01+p02)]==0),0)
  npm12<-ifelse(nvat12>0,sum(fix0[(7+nvat01+nvat02):(7+nvat01+nvat02+nvat12)]==0),0)
  npm12Y<-ifelse(p12>0,sum(fix0[(7+nvat01+nvat02+nvat12+p01+p02):(6+nvat01+nvat02+nvat12+p01+p12+p02)]==0),0)
  
  nvat01<-npm01+npm01Y
  nvat02<-npm02+npm02Y
  nvat12<-npm12+npm12Y
  
  dimnva01<-ifelse(nvat01==0,1,nvat01)
  dimnva02<-ifelse(nvat02==0,1,nvat02)
  dimnva12<-ifelse(nvat12==0,1,nvat12)
  
  outputall<-list()
  length(outputall)<-Nsample
 
  browser()
  if(nlambda>1){
  
  for(idsample in 1:Nsample){
    
    # do prediction for the sample 
    
    browser()
    if(modelY$method=="INLA"){
      
      dataY<-INLAidmpredY(timeVar=timeVar,
                          truncated=troncature,
                          formLong=formLong,
                          dataSurv=dataSurv,
                          dataLongi=dataLongi,
                          id=id,
                          Nsample=1,
                          t0=t0,t1=t1,t2=t2,t3=t3,
                          ctime=ctime,
                          modelY=modelY,
                          seed=seed+idsample,
                          BLUP=BLUP,
                          nproc=1,
                          clustertype=clustertype)
    }else{
      
      dataY<-JMidmpredY(timeVar=timeVar,
                          truncated=troncature,
                          formLong=formLong,
                          dataSurv=dataSurv,
                          dataLongi=dataLongi,
                          id=id,
                          Nsample=1,
                          t0=t0,t1=t1,t2=t2,t3=t3,
                          ctime=ctime,
                          modelY=modelY,
                          seed=seed+idsample,
                          BLUP=BLUP)
    }
    
    for( m in unique(c(outcome01,outcome02,outcome12))){
      subdata<-dataY[dataY$Outcome==m,]
      x<-table(subdata[,colnames(subdata)%in%id])
      if(any(x!=NtimePoints)){stop("Prediction of marker ",m," could not be perform for each quadrature points, try Ypredmethod equi")}
      
    }
    
    dataY$Outcome<-as.character(dataY$Outcome)
    # attention if NtimePoints equidistant with INLA then NtimePoints takes 
    # need ID to be numeric -- then 
    dataY[,colnames(dataY)%in%id]<-as.numeric(dataY[,colnames(dataY)%in%id])
    # to keep tracks of time order for each individual 
    dataY$order<-as.numeric(ave(dataY[,colnames(dataY)%in%id], cumsum(c(TRUE, diff(dataY[,colnames(dataY)%in%id]) != 0)), FUN = seq_along))
    
    if(scale.X==T){
      
      # Compute group means and sds
      ym <- tapply(dataY[[4]], dataY$Outcome, mean)
      ys <- tapply(dataY[[4]], dataY$Outcome, sd)
      
      # Normalize (min-max) within each group
      dataY[[4]] <- ave(dataY[[4]], dataY$Outcome,
                        FUN = function(x) (x - min(x)) / (max(x) - min(x)))
      
    }
    
    
    if(length(outcome01)>=1){
      y01k<-dataY[dataY$Outcome%in%outcome01,]
      # order  by individual and timeline 
      y01k<-y01k[order(y01k[,colnames(y01k)%in%id],y01k$order),4]
      ve01k<-c(ve01,y01k)
      
    }else{
      ve01k<-ve01
    }
    
    if(length(outcome02)>=1){
      y02k<-dataY[dataY$Outcome%in%outcome02,]
      # order  by individual and timeline 
      y02k<-y02k[order(y02k[,colnames(y02k)%in%id],y02k$order),]
      ve02k<-c(ve02,y02k)
    }else{
      ve02k<-ve02
    }
    
    if(length(outcome12)>=1){
      y12k<-dataY[dataY$Outcome%in%outcome12,]
      # order  by individual and timeline 
      y12k<-y12k[order(y12k[,colnames(y12k)%in%id],y12k$order),4]
      ve12k<-c(ve12,y12k)
    }else{
      ve12k<-ve12
    }
    
    browser()

    
    outputall[[idsample]]<- idm.penalty.weib(b=b,
                                             fix0=fix0,
                                             size_V=size_V,
                                             clustertype=clustertype,
                                             epsa=epsa,
                                             epsb=epsb,
                                             epsd=epsd,
                                             eps.eigen=eps.eigen,
                                             nproc=nproc,
                                             maxiter=maxiter,
                                             maxiter.pena=maxiter.pena,
                                             ctime=ctime,
                                             N=N,
                                             ve01=ve01k,
                                             ve02=ve02k,
                                             ve12=ve12k,
                                             dimnva01=dimnva01,
                                             dimnva02=dimnva02,
                                             dimnva12=dimnva12,
                                             nvat01=nvat01,
                                             nvat02=nvat02,
                                             nvat12=nvat12,
                                             t0=t0,
                                             t1=t1,
                                             t2=t2,
                                             t3=t3,
                                             troncature=troncature,
                                             nlambda01=nlambda01,
                                             lambda01=lambda01,
                                             nlambda02=nlambda02,
                                             lambda02=lambda02,
                                             nlambda12=nlambda12,
                                             lambda12=lambda12,
                                             alpha=alpha,
                                             penalty.factor=penalty.factor,
                                             penalty=penalty,
                                             gausspoint=15,
                                             analytics=T,
                                             partialH=partialH)
    
  }
    
  }else{
    if(nproc==1){
      for(idsample in 1:Nsample){
        
        if(modelY$method=="INLA"){
          
          dataY<-INLAidmpredY(timeVar=timeVar,
                              truncated=troncature,
                              formLong=formLong,
                              dataSurv=dataSurv,
                              dataLongi=dataLongi,
                              id=id,
                              Nsample=1,
                              t0=t0,t1=t1,t2=t2,t3=t3,
                              ctime=ctime,
                              modelY=modelY,
                              seed=seed+idsample,
                              BLUP=BLUP,
                              nproc=1,
                              clustertype=clustertype)
        }else{
          
          dataY<-JMidmpredY(timeVar=timeVar,
                            truncated=troncature,
                            formLong=formLong,
                            dataSurv=dataSurv,
                            dataLongi=dataLongi,
                            id=id,
                            Nsample=1,
                            t0=t0,t1=t1,t2=t2,t3=t3,
                            ctime=ctime,
                            modelY=modelY,
                            seed=seed+idsample,
                            BLUP=BLUP)
        }
        
        for( m in unique(c(outcome01,outcome02,outcome12))){
          subdata<-dataY[dataY$Outcome==m,]
          x<-table(subdata[,colnames(subdata)%in%id])
          if(any(x!=NtimePoints)){stop("Prediction of marker ",m," could not be perform for each quadrature points, try Ypredmethod equi")}
          
        }
        
        dataY$Outcome<-as.character(dataY$Outcome)
        # attention if NtimePoints equidistant with INLA then NtimePoints takes 
        # need ID to be numeric -- then 
        dataY[,colnames(dataY)%in%id]<-as.numeric(dataY[,colnames(dataY)%in%id])
        # to keep tracks of time order for each individual 
        dataY$order<-as.numeric(ave(dataY[,colnames(dataY)%in%id], cumsum(c(TRUE, diff(dataY[,colnames(dataY)%in%id]) != 0)), FUN = seq_along))
        
        if(scale.X==T){
          
          # Compute group means and sds
          ym <- tapply(dataY[[4]], dataY$Outcome, mean)
          ys <- tapply(dataY[[4]], dataY$Outcome, sd)
          
          # Normalize (min-max) within each group
          dataY[[4]] <- ave(dataY[[4]], dataY$Outcome,
                            FUN = function(x) (x - min(x)) / (max(x) - min(x)))
          
        }
        
        
        if(length(outcome01)>=1){
          y01k<-dataY[dataY$Outcome%in%outcome01,]
          # order  by individual and timeline 
          y01k<-y01k[order(y01k[,colnames(y01k)%in%id],y01k$order),4]
          ve01k<-c(ve01,y01k)
          
        }else{
          ve01k<-ve01
        }
        
        if(length(outcome02)>=1){
          y02k<-dataY[dataY$Outcome%in%outcome02,]
          # order  by individual and timeline 
          y02k<-y02k[order(y02k[,colnames(y02k)%in%id],y02k$order),]
          ve02k<-c(ve02,y02k)
        }else{
          ve02k<-ve02
        }
        
        if(length(outcome12)>=1){
          y12k<-dataY[dataY$Outcome%in%outcome12,]
          # order  by individual and timeline 
          y12k<-y12k[order(y12k[,colnames(y12k)%in%id],y12k$order),4]
          ve12k<-c(ve12,y12k)
        }else{
          ve12k<-ve12
        }
        
        outputall[[idsample]]<-idm.penalty.weib(b=b,
                                                fix0=fix0,
                                                size_V=size_V,
                                                clustertype=clustertype,
                                                epsa=epsa,
                                                epsb=epsb,
                                                epsd=epsd,
                                                eps.eigen=eps.eigen,
                                                nproc=nproc,
                                                maxiter=maxiter,
                                                maxiter.pena=maxiter.pena,
                                                ctime=ctime,
                                                N=N,
                                                ve01=ve01k,
                                                ve02=ve02k,
                                                ve12=ve12k,
                                                dimnva01=dimnva01,
                                                dimnva02=dimnva02,
                                                dimnva12=dimnva12,
                                                nvat01=nvat01,
                                                nvat02=nvat02,
                                                nvat12=nvat12,
                                                t0=t0,
                                                t1=t1,
                                                t2=t2,
                                                t3=t3,
                                                troncature=troncature,
                                                nlambda01=nlambda01,
                                                lambda01=lambda01,
                                                nlambda02=nlambda02,
                                                lambda02=lambda02,
                                                nlambda12=nlambda12,
                                                lambda12=lambda12,
                                                alpha=alpha,
                                                penalty.factor=penalty.factor,
                                                penalty=penalty,
                                                gausspoint=15,
                                                analytics=T,
                                                partialH=partialH)
      }
      }else{
        
        if(is.null(clustertype)){
          clustpar <- parallel::makeCluster(nproc)#, outfile="")
        }
        else{
          clustpar <- parallel::makeCluster(nproc, type=clustertype)#, outfile="")
        }
        
        doParallel::registerDoParallel(clustpar)
        
        outputall<-foreach::foreach(idsample=1:Nsample,
                                 .combine = 'c',.multicombine=TRUE,
                                 .errorhandling = "remove")%dopar%{
                                   
                                   
                                   if(modelY$method=="INLA"){
                                     
                                     dataY<-INLAidmpredY(timeVar=timeVar,
                                                         truncated=troncature,
                                                         formLong=formLong,
                                                         dataSurv=dataSurv,
                                                         dataLongi=dataLongi,
                                                         id=id,
                                                         Nsample=1,
                                                         t0=t0,t1=t1,t2=t2,t3=t3,
                                                         ctime=ctime,
                                                         modelY=modelY,
                                                         seed=seed+idsample,
                                                         BLUP=BLUP,
                                                         nproc=1,
                                                         clustertype=clustertype)
                                   }else{
                                     
                                     dataY<-JMidmpredY(timeVar=timeVar,
                                                       truncated=troncature,
                                                       formLong=formLong,
                                                       dataSurv=dataSurv,
                                                       dataLongi=dataLongi,
                                                       id=id,
                                                       Nsample=1,
                                                       t0=t0,t1=t1,t2=t2,t3=t3,
                                                       ctime=ctime,
                                                       modelY=modelY,
                                                       seed=seed+idsample,
                                                       BLUP=BLUP)
                                   }
                                   
                                   for( m in unique(c(outcome01,outcome02,outcome12))){
                                     subdata<-dataY[dataY$Outcome==m,]
                                     x<-table(subdata[,colnames(subdata)%in%id])
                                     if(any(x!=NtimePoints)){stop("Prediction of marker ",m," could not be perform for each quadrature points, try Ypredmethod equi")}
                                     
                                   }
                                   
                                   dataY$Outcome<-as.character(dataY$Outcome)
                                   # attention if NtimePoints equidistant with INLA then NtimePoints takes 
                                   # need ID to be numeric -- then 
                                   dataY[,colnames(dataY)%in%id]<-as.numeric(dataY[,colnames(dataY)%in%id])
                                   # to keep tracks of time order for each individual 
                                   dataY$order<-as.numeric(ave(dataY[,colnames(dataY)%in%id], cumsum(c(TRUE, diff(dataY[,colnames(dataY)%in%id]) != 0)), FUN = seq_along))
                                   
                                   if(scale.X==T){
                                     
                                     # Compute group means and sds
                                     ym <- tapply(dataY[[4]], dataY$Outcome, mean)
                                     ys <- tapply(dataY[[4]], dataY$Outcome, sd)
                                     
                                     # Normalize (min-max) within each group
                                     dataY[[4]] <- ave(dataY[[4]], dataY$Outcome,
                                                       FUN = function(x) (x - min(x)) / (max(x) - min(x)))
                                     
                                   }
                                   
                                   
                                   if(length(outcome01)>=1){
                                     y01k<-dataY[dataY$Outcome%in%outcome01,]
                                     # order  by individual and timeline 
                                     y01k<-y01k[order(y01k[,colnames(y01k)%in%id],y01k$order),4]
                                     ve01k<-c(ve01,y01k)
                                     
                                   }else{
                                     ve01k<-ve01
                                   }
                                   
                                   if(length(outcome02)>=1){
                                     y02k<-dataY[dataY$Outcome%in%outcome02,]
                                     # order  by individual and timeline 
                                     y02k<-y02k[order(y02k[,colnames(y02k)%in%id],y02k$order),]
                                     ve02k<-c(ve02,y02k)
                                   }else{
                                     ve02k<-ve02
                                   }
                                   
                                   if(length(outcome12)>=1){
                                     y12k<-dataY[dataY$Outcome%in%outcome12,]
                                     # order  by individual and timeline 
                                     y12k<-y12k[order(y12k[,colnames(y12k)%in%id],y12k$order),4]
                                     ve12k<-c(ve12,y12k)
                                   }else{
                                     ve12k<-ve12
                                   }
                                   
                                   outputNsample<-idm.penalty.weib(b=b,
                                                                   fix0=fix0,
                                                                   size_V=size_V,
                                                                   clustertype=clustertype,
                                                                   epsa=epsa,
                                                                   epsb=epsb,
                                                                   epsd=epsd,
                                                                   eps.eigen=eps.eigen,
                                                                   nproc=1,
                                                                   maxiter=maxiter,
                                                                   maxiter.pena=maxiter.pena,
                                                                   ctime=ctime,
                                                                   N=N,
                                                                   ve01=ve01k,
                                                                   ve02=ve02k,
                                                                   ve12=ve12k,
                                                                   dimnva01=dimnva01,
                                                                   dimnva02=dimnva02,
                                                                   dimnva12=dimnva12,
                                                                   nvat01=nvat01,
                                                                   nvat02=nvat02,
                                                                   nvat12=nvat12,
                                                                   t0=t0,
                                                                   t1=t1,
                                                                   t2=t2,
                                                                   t3=t3,
                                                                   troncature=troncature,
                                                                   nlambda01=nlambda01,
                                                                   lambda01=lambda01,
                                                                   nlambda02=nlambda02,
                                                                   lambda02=lambda02,
                                                                   nlambda12=nlambda12,
                                                                   lambda12=lambda12,
                                                                   alpha=alpha,
                                                                   penalty.factor=penalty.factor,
                                                                   penalty=penalty,
                                                                   gausspoint=15,
                                                                   analytics=T,
                                                                   partialH=partialH)
                                    
                                   return(outputNsample)
                                 }
        
        parallel::stopCluster(clustpar)
      }
  }
  
  
  return(output=outputall)
}



