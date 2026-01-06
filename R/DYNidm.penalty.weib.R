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


DYNidm.penalty.weib<-function(b,fix0,size_V,
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
                   p01,p02,p12,assoc,
                   dimp01,dimp02,dimp12,scale.X,defpositive,warmstart){
  


  # need to keep original fix to calculate for beta 
  
  fix00<-fix0
  
  # create grid 3
  lambda<-expand.grid(lambda01,lambda02,lambda12)
  lambda<-unique(lambda)
  nlambda<-dim(lambda)[1]
  

  # need to check that same variable in each transition :
  
  # Initiate value of spline
  
   s.start<-b[1:6]
  
  # Initiate value of eta : all the same for each lambda
  
  beta.start<-b[(6+1):(size_V)]
  
  # fix0 will be used to calculate derivatives and second derivatives only
  # for Beta and not modelPar01,02,12
  fix0[1:6]<-rep(1,6)
  fix0.beta<-fix00
  fix0.beta[(6+1):size_V]<-rep(1,size_V-6)

  
  npm<-sum(fix0==0)
  npmweib<-sum(fix00[1:6]==0)
  
  npm01<-ifelse(nvat01>0,sum(fix0[7:(7+nvat01-1)]==0),0)
  npm01Y<-ifelse(p01>0,sum(fix0[(7+nvat01+nvat02+nvat12):(6+nvat01+nvat02+nvat12+p01)]==0),0)
  npm02<-ifelse(nvat02>0,sum(fix0[(7+nvat01):(6+nvat01+nvat02)]==0),0)
  npm02Y<-ifelse(p02>0,sum(fix0[(7+nvat01+nvat02+nvat12+p01):(6+nvat01+nvat02+nvat12+p01+p02)]==0),0)
  npm12<-ifelse(nvat12>0,sum(fix0[(7+nvat01+nvat02):(7+nvat01+nvat02+nvat12)]==0),0)
  npm12Y<-ifelse(p12>0,sum(fix0[(7+nvat01+nvat02+nvat12+p01+p02):(6+nvat01+nvat02+nvat12+p01+p12+p02)]==0),0)
  
  if(partialH==F){
    min<-(npm*(npm+1))/2
  }else{
    min<-npm
  }
  
  outputall<-list()
  length(outputall)<-Nsample
  browser()
  if(nlambda>1){
  
  for(idsample in 1:Nsample){
    
    # do prediction for the sample 
    
    if(modelY$method=="INLA"){
      
      dataY<-INLAidmpredY(timeVar=timeVar,
                          truncated=troncature,
                          formLong=formLong,
                          dataSurv=dataSurv,
                          dataLongi=dataLongi,
                          id=id,
                          t0=t0,t1=t1,t2=t2,t3=t3,
                          assoc=assoc,
                          ctime=ctime,
                          modelY=modelY,
                          seed=seed+idsample,
                          BLUP=BLUP,
                          scale.X=scale.X)
    
      
    }else{
      
      dataY<-JMidmpredY(timeVar=timeVar,
                          truncated=troncature,
                          formLong=formLong,
                          dataSurv=dataSurv,
                          dataLongi=dataLongi,
                          id=id,
                          t0=t0,t1=t1,t2=t2,t3=t3,
                          ctime=ctime,
                          modelY=modelY,
                          seed=seed+idsample,
                          BLUP=BLUP,
                        assoc=assoc,
                        scale.X=scale.X)
   
      
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

    
    if(length(outcome01)>=1){
      y01k<-dataY[dataY$Outcome%in%outcome01,]
      # order  by individual and timeline 
      y01k<-y01k[order(y01k[,colnames(y01k)%in%id],y01k$order),4]
      
      
    }else{
      y01k<-rep(0,N*NtimePoints)
    }
    
    if(length(outcome02)>=1){
      y02k<-dataY[dataY$Outcome%in%outcome02,]
      # order  by individual and timeline 
      y02k<-y02k[order(y02k[,colnames(y02k)%in%id],y02k$order),4]
      
    }else{
      y02k<-rep(0,N*NtimePoints)
    }
    
    if(length(outcome12)>=1){
      y12k<-dataY[dataY$Outcome%in%outcome12,]
      # order  by individual and timeline 
      y12k<-y12k[order(y12k[,colnames(y12k)%in%id],y12k$order),4]
      
    }else{
      y12k<-rep(0,N*NtimePoints)
    }
    
 


  if(nproc >1){
    
    outputNsample<-DYNidm.penalty.weib.nproc(beta.start=beta.start,
                                             s.start=s.start,
                                             npm01=npm01,
                                             npm02=npm02,
                                             npm12=npm12,
                                             npm01Y=npm01Y,
                                             npm02Y=npm02Y,
                                             npm12Y=npm12Y,
                                             npm=npm,
                                             fix0=fix0,
                                             fix00=fix00,
                                             fix0.beta=fix0.beta,
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
                                    ve01=ve01,
                                    ve02=ve02,
                                    ve12=ve12,
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
                                    nlambda=nlambda,
                                    lambda=lambda,
                                    alpha=alpha,
                                    penalty.factor=penalty.factor,
                                    penalty=penalty,
                                    partialH=partialH,
                                    Nsample=Nsample,
                                    NtimePoints=NtimePoints,
                                    p01=p01,
                                    p02=p02,
                                    p12=p12,
                                    dimp01=dimp01,
                                    dimp02=dimp02,
                                    dimp12=dimp12,
                                    defpositive=defpositive,
                                    warmstart=warmstart,
                                    y01k=y01k,
                                    y02k=y02k,
                                    y12k=y12k,
                                    min=min)
    
  }else{
    

      outputNsample<-DYNidm.penalty.weib.nonproc(beta.start=beta.start,
                                               s.start=s.start,
                                               npm01=npm01,
                                               npm02=npm02,
                                               npm12=npm12,
                                               npm01Y=npm01Y,
                                               npm02Y=npm02Y,
                                               npm12Y=npm12Y,
                                               npm=npm,
                                               fix0=fix0,
                                               fix00=fix00,
                                               fix0.beta=fix0.beta,
                                               size_V=size_V,
                                               epsa=epsa,
                                               epsb=epsb,
                                               epsd=epsd,
                                               eps.eigen=eps.eigen,
                                               maxiter=maxiter,
                                               maxiter.pena=maxiter.pena,
                                               ctime=ctime,
                                               N=N,
                                               ve01=ve01,
                                               ve02=ve02,
                                               ve12=ve12,
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
                                               nlambda=nlambda,
                                               lambda=lambda,
                                               alpha=alpha,
                                               penalty.factor=penalty.factor,
                                               penalty=penalty,
                                               partialH=partialH,
                                               Nsample=Nsample,
                                               NtimePoints=NtimePoints,
                                               p01=p01,
                                               p02=p02,
                                               p12=p12,
                                               dimp01=dimp01,
                                               dimp02=dimp02,
                                               dimp12=dimp12,
                                               defpositive=defpositive,
                                               warmstart=warmstart,
                                               y01k=y01k,
                                               y02k=y02k,
                                               y12k=y12k,
                                               min=min)
      
  
  }
    
    outputall[[idsample]]<-outputNsample
  }
    
  }else{
    if(nproc==1){
      
      outputall<-DYNidm.penalty.weib.nonproc.onepenalty(beta.start=beta.start,
                                                      s.start=s.start,
                                                      npm01=npm01,
                                                      npm02=npm02,
                                                      npm12=npm12,
                                                      npm01Y=npm01Y,
                                                      npm02Y=npm02Y,
                                                      npm12Y=npm12Y,
                                                      npm=npm,
                                                      fix0=fix0,
                                                      fix00=fix00,
                                                      fix0.beta=fix0.beta,
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
                                                      ve01=ve01,
                                                      ve02=ve02,
                                                      ve12=ve12,
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
                                                      nlambda=nlambda,
                                                      lambda=lambda,
                                                      alpha=alpha,
                                                      penalty.factor=penalty.factor,
                                                      penalty=penalty,
                                                      partialH=partialH,
                                                      Nsample=Nsample,
                                                      NtimePoints=NtimePoints,
                                                      p01=p01,
                                                      p02=p02,
                                                      p12=p12,
                                                      dimp01=dimp01,
                                                      dimp02=dimp02,
                                                      dimp12=dimp12,
                                                      defpositive=defpositive,
                                                      min=min,
                                                      modelY=modelY,
                                                      outcome01=outcome01,
                                                      outcome02=outcome02,
                                                      outcome12=outcome12,
                                                      timeVar=timeVar,
                                                      formLong=formLong,
                                                      dataSurv=dataSurv,
                                                      dataLongi=dataLongi,
                                                      id=id,
                                                      assoc=assoc,
                                                      seed=seed,
                                                      BLUP=BLUP)
      
      
      }else{
        
        
        outputall<-DYNidm.penalty.weib.nproc.onepenalty(beta.start=beta.start,
                                             s.start=s.start,
                                             npm01=npm01,
                                             npm02=npm02,
                                             npm12=npm12,
                                             npm01Y=npm01Y,
                                             npm02Y=npm02Y,
                                             npm12Y=npm12Y,
                                             npm=npm,
                                             fix0=fix0,
                                             fix00=fix00,
                                             fix0.beta=fix0.beta,
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
                                             ve01=ve01,
                                             ve02=ve02,
                                             ve12=ve12,
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
                                             nlambda=nlambda,
                                             lambda=lambda,
                                             alpha=alpha,
                                             penalty.factor=penalty.factor,
                                             penalty=penalty,
                                             partialH=partialH,
                                             Nsample=Nsample,
                                             NtimePoints=NtimePoints,
                                             p01=p01,
                                             p02=p02,
                                             p12=p12,
                                             dimp01=dimp01,
                                             dimp02=dimp02,
                                             dimp12=dimp12,
                                             defpositive=defpositive,
                                             min=min,
                                             modelY=modelY,
                                             outcome01=outcome01,
                                             outcome02=outcome02,
                                             outcome12=outcome12,
                                             timeVar=timeVar,
                                             formLong=formLong,
                                             dataSurv=dataSurv,
                                             dataLongi=dataLongi,
                                             id=id,
                                             assoc=assoc,
                                             seed=seed,
                                             BLUP=BLUP)
        
       
      }
  }
  
  
  return(output=outputall)
}



