### Code:
##' @title Penalised Illness-death model algorithm with M-spline baseline risk
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


DYNidm.penalty.splines<-function(b,fix0,size_V,size_spline,
                   clustertype,epsa,epsb,epsd,eps.eigen,nproc,maxiter,maxiter.pena,
                   ctime,N,
                   ve01,ve02,ve12,dimnva01,dimnva02,dimnva12,nvat01,nvat02,nvat12,
                   t0,t1,t2,t3,troncature,knots01,knots02,knots12,
                   nknots01,nknots02,nknots12,
                   nlambda01,lambda01,nlambda02,lambda02,nlambda12,lambda12,
                   alpha,penalty.factor,penalty,partialH,
                   modelY,dataLongi,dataSurv,
                   Nsample,BLUP,seed,timeVar,id,formLong,
                   outcome01,outcome02,
                   outcome12,NtimePoints,
                   p01,p02,p12, assoc,
                   dimp01,dimp02,dimp12,scale.X){
  


  # need to keep original fix to calculate for beta 

  V0<-NA
  fix00<-fix0
  
  # create grid 3
  lambda<-expand.grid(lambda01,lambda02,lambda12)
  lambda<-unique(lambda)
  nlambda<-dim(lambda)[1]
  
  # computation pbr 
  pbr_compu<-0
  
  # combine model 
  combine_lambda_mla<-function(x,newx){
    
    if(newx$combine==2){
      list(b=cbind(x$b,newx$b),
           V=cbind(x$V,newx$V),
           H=cbind(x$H,newx$H),
           fix=cbind(x$fix,newx$fix),
           lambda=cbind(x$lambda,newx$lambda),
           alpha=c(x$alpha,newx$alpha),
           fn.value=c(x$fn.value,newx$fn.value),
           fn.value.pena=c(x$fn.value.pena,newx$fn.value.pena),
           ni=c(x$ni,newx$ni),
           istop=c(x$istop,newx$istop),
           ca.beta=cbind(x$ca.beta,newx$ca.beta),
           ca.spline=cbind(x$ca.spline,newx$ca.spline),
           ca.validity=cbind(x$ca.validity,newx$ca.validity),
           cb=cbind(x$cb,newx$cb))
      
    }else{
      list(b=cbind(x$b,newx$b),
           V=cbind(x$V,newx$V),
           H=cbind(x$H,newx$H),
           fix=cbind(x$fix,newx$fix),
           lambda=cbind(x$lambda,newx$lambda),
           alpha=c(x$alpha,newx$alpha),
           fn.value=c(x$fn.value,newx$fn.value),
           fn.value.pena=c(x$fn.value.pena,newx$fn.value.pena),
           ni=c(x$ni,newx$ni),
           istop=c(x$istop,newx$istop),
           ca.beta=cbind(x$ca.beta,newx$ca.beta),
           ca.spline=cbind(x$ca.spline,newx$ca.spline),
           ca.validity=cbind(x$ca.validity,newx$ca.validity),
           cb=cbind(x$cb,newx$cb))}
    
  }
  
  
  # need to check that same variable in each transition :
  
  # Initiate value of spline
   s.start<-b[1:size_spline]
  
  # Initiate value of eta : all the same for each lambda
  
  beta.start<-b[(size_spline+1):(size_V)]
  
  
  combine<-0
  
  # fix0 will be used to calculate derivatives and second derivatives only
  # for Beta and not modelPar01,02,12
  fix0[1:size_spline]<-rep(1,size_spline)
  fix0.beta<-fix00
  fix0.beta[(size_spline+1):size_V]<-rep(1,size_V-size_spline)

  npm<-sum(fix0==0)
  npmspline<-sum(fix00[1:size_spline]==0)
  
  npm01<-ifelse(nvat01>0,sum(fix0[(size_spline+1):(size_spline+nvat01)]==0),0)
  npm01Y<-ifelse(p01>0,sum(fix0[(size_spline+1+nvat01+nvat02+nvat12):(size_spline+nvat01+nvat02+nvat12+p01)]==0),0)
  
  npm02<-ifelse(nvat02>0,sum(fix0[(size_spline+1+nvat01):(size_spline+nvat01+nvat02)]==0),0)
  npm02Y<-ifelse(p02>0,sum(fix0[(size_spline+1+nvat01+nvat02+nvat12+p01):(size_spline+nvat01+nvat02+nvat12+p01+p02)]==0),0)
  npm12<-ifelse(nvat12>0,sum(fix0[(size_spline+1+nvat01+nvat02):(size_spline+1+nvat01+nvat02+nvat12)]==0),0)
  npm12Y<-ifelse(p12>0,sum(fix0[(size_spline+1+nvat01+nvat02+nvat12+p01+p02):(size_spline+nvat01+nvat02+nvat12+p01+p12+p02)]==0),0)
  
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
    
    if(modelY$method=="INLA"){
      # so far cannot do nproc > 1 as not enough memory
      dataY<-INLAidmpredY(timeVar=timeVar,
                          truncated=troncature,
                          formLong=formLong,
                          dataSurv=dataSurv,
                          dataLongi=dataLongi,
                          id=id,
                          Nsample=1,
                          t0=t0,t1=t1,t2=t2,t3=t3,
                          assoc=assoc,
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
      
      
    }else{
      y01k<-rep(0,N*NtimePoints)
    }
    
    if(length(outcome02)>=1){
      y02k<-dataY[dataY$Outcome%in%outcome02,]
      # order  by individual and timeline 
      y02k<-y02k[order(y02k[,colnames(y02k)%in%id],y02k$order),]
      
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
    
    
    if(is.null(clustertype)){
      clustpar <- parallel::makeCluster(nproc)#, outfile="")
    }
    else{
      clustpar <- parallel::makeCluster(nproc, type=clustertype)#, outfile="")
    }
    
    doParallel::registerDoParallel(clustpar)
    
    
    
    id.lambda<-NULL # for cran check 
    if(partialH==F){
    outputNsample<-foreach::foreach(id.lambda=1:nlambda,
                             .combine = combine_lambda_mla,
                             .errorhandling = "remove")%dopar%{
                               
                               # computation pbr 
                               
                               pbr_compu<-0
                               
                               beta<-beta.start
                               s<-s.start
                               
                               if(idsample==1){
                                 beta<-beta.start
                                 s<-s.start
                               }else{
                                 s<-output[[1]]$b[1:size_spline,id.lambda]
                                 beta<-output[[1]]$b[(size_spline+1):size_V,id.lambda]
                              }
                               
                               converged<-F
                               ite<-0
                               # if beta not change do not need to recalculate weights 
                               H<-T
                               
                               eval.cv.spline<-rep(NA,maxiter+1)
                               eval.cv.beta<-rep(NA,maxiter+1)
                               eval.cv.loglik<-rep(NA,maxiter+1)
                               eval.loglik<-rep(NA,maxiter+1)
                               eval.validity<-rep(NA,maxiter+1)
                               
                               
                              
                               
                               while(converged==F & ite<=maxiter){
                                 
                                 
                                 b<-c(s,beta)
                                 bfix<-b[fix0==1]
                                 b<-b[fix0==0]
                                 # derivative of loglik
                                 
                                 output<-DYNderiva(
                                                    h=1e-8,b=b,
                                                    npm=npm,
                                                    npar=size_V,
                                                    bfix=bfix,
                                                    fix=fix0,
                                                    zi01=knots01,
                                                    zi02=knots02,
                                                    zi12=knots12,
                                                    ctime=ctime,
                                                    no=N,
                                                    nz01=nknots01,
                                                    nz02=nknots02,
                                                    nz12=nknots12,
                                                    ve01=ve01,
                                                    ve02=ve02,
                                                    ve12=ve12,
                                                    dimnva01=dimnva01,
                                                    dimnva02=dimnva02,
                                                    dimnva12=dimnva12,
                                                    nva01=nvat01,
                                                    nva02=nvat02,
                                                    nva12=nvat12,
                                                    t0=t0,
                                                    t1=t1,
                                                    t2=t2,
                                                    t3=t3,
                                                    troncature=troncature,
                                                    y01=y01k,
                                                    y02=y02k,
                                                    y12=y12k,
                                                    p01=p01,
                                                    p02=p02,
                                                    p12=p12,
                                                    dimp01=dimp01,
                                                    dimp02=dimp02,
                                                    dimp12=dimp12,
                                                    Ntime=NtimePoints)
                                 output<-output$v
                                
                                 if(ite==0){# loglik penalised
                                   fn.value<-gaussDYNidmlLikelihoodpena(b=b,
                                                                    npm=npm,
                                                                    npar=size_V,
                                                                    bfix=bfix,
                                                                    fix=fix0,
                                                                    zi01=knots01,
                                                                    zi02=knots02,
                                                                    zi12=knots12,
                                                                    ctime=ctime,
                                                                    no=N,
                                                                    nz01=nknots01,
                                                                    nz02=nknots02,
                                                                    nz12=nknots12,
                                                                    ve01=ve01,
                                                                    ve02=ve02,
                                                                    ve12=ve12,
                                                                    dimnva01=dimnva01,
                                                                    dimnva02=dimnva02,
                                                                    dimnva12=dimnva12,
                                                                    nva01=nvat01,
                                                                    nva02=nvat02,
                                                                    nva12=nvat12,
                                                                    t0=t0,
                                                                    t1=t1,
                                                                    t2=t2,
                                                                    t3=t3,
                                                                    troncature=troncature,
                                                                    y01=y01k,
                                                                    y02=y02k,
                                                                    y12=y12k,
                                                                    p01=p01,
                                                                    p02=p02,
                                                                    p12=p12,
                                                                    dimp01=dimp01,
                                                                    dimp02=dimp02,
                                                                    dimp12=dimp12,
                                                                    Ntime=NtimePoints,
                                                                   
                                                                    lambda=lambda[id.lambda,],
                                                                    alpha=alpha,
                                                                    penalty.factor=penalty.factor,
                                                                    penalty=penalty)
                                 }
                                 
                                 if(any(is.na(output))|any(output==Inf) |any(output==-Inf)){
                                   warning("Computational error for calculation of the hessian : division by 0 or Infinite value")
                                   if(ite==0){
                                     
                                    
                                     fu<-output[(min+1):length(output)]
                                     V<-matrix(0,nrow=npm,ncol=npm)
                                     V[upper.tri(V,diag=T)]<-output[1:min]
                                     V<-V+t(V)
                                     diag(V)<-diag(V)/2
                                     # deriva gives information matrix
                                     tr <- sum(diag(V))/npm
                                     V0<-V}
                                   ite<-ite+1
                                   pbr_compu<-1
                                   break
                                 }
                                 
                                
                                 fu<-output[(min+1):length(output)]
                                 V<-matrix(0,nrow=npm,ncol=npm)
                                 V[upper.tri(V,diag=T)]<-output[1:min]
                                 V<-V+t(V)
                                 diag(V)<-diag(V)/2
                                 
                                 # deriva gives information matrix
                                 tr <- sum(diag(V))/npm
                                 V0<-V
                                 
                                 eigen.values<-eigen(V,symmetric=T,only.values=T)$values
                                 
                                 idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                 
                                 
                                 idpos0<-idpos
                                 
                                 ncount<-da<-ga<-0
 
                                 while(idpos != 0){
                                   
                                   if(ncount==0){ 
                                     ga <- 0.01
                                     da <- 1E-2
                                   }else{
                                     if(((ncount <= 3) | (ga >= 1)) ){
                                       da <- da * 5
                                     }else{# if ncount > 10 only update ga 
                                       ga <- ga * 5
                                       # do not put ga at 1 as no countmax otherwise infinite while 
                                       if(ga > 1) ga <- 1
                                     }
                                   }
                                   
                                   ncount <- ncount + 1
                                   
                                   diagV <- diag(V)
                                   # put abs (1-ga) better than 1-ga cause ga can now be >1
                                   diagV<-ifelse(diagV!=0,diagV+da*(abs((1.e0-ga))*abs(diagV)+ga*tr),
                                                 da*ga*tr)
                                   
                                   diag(V)<-diagV
                                   # if we have a convex log-vraisemblance in eta then :
                                   # all eigen  values of the hessienne are >0.
                                   
                                   if(sum(V==Inf)>0|sum(V==-Inf)>0){break}
                                   eigen.values<-eigen(V,symmetric=T,only.values=T)$values
                                   # check if hessienne defined positive
                                   
                                   idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                   
                                   # if(def.positive==T){
                                   #   idpos<-ifelse(any(eigen.values<=0),1,0)
                                   # }else{idpos<-ifelse(any(abs(eigen.values)==0),1,0)}
                                   
                                 }
                                 
                                
                                 if(idpos!=0){
                                  
                                   warning("Hessian not defined positive")
                                   pbr_compu<-2
                                   ite<-ite+1
                                   break
                                 }
                                 
                                 
                                # update for beta 
                                 output.cv<-DYNcv.model(beta=beta,
                                                     nva01=npm01,
                                                     nva02=npm02,
                                                     nva12=npm12,
                                                     nva01Y=npm01Y,
                                                     nva02Y=npm02Y,
                                                     nva12Y=npm12Y,
                                                     fix=fix0[(size_spline+1):size_V],
                                                     penalty.factor=penalty.factor,
                                                     penalty=penalty,
                                                     v=V,
                                                     fu=fu,
                                                     lambda=lambda[id.lambda,],
                                                     alpha=alpha
                                 )
                                 
                                 # verify validity of parameters update 
                                 # and that we are better than previous estimates 
                                 
                                 b<-c(s,output.cv$b)
                                 
                                 betanew<-b[(size_spline+1):size_V]
                                 
                                 # penalised loglik see if inferior to previous
                                 res<-gaussDYNidmlLikelihoodpena(b=b,
                                                             npm=size_V,
                                                             npar=size_V,
                                                             bfix=1,
                                                             fix=rep(0,size_V),
                                                             zi01=knots01,
                                                             zi02=knots02,
                                                             zi12=knots12,
                                                             ctime=ctime,
                                                             no=N,
                                                             nz01=nknots01,
                                                             nz02=nknots02,
                                                             nz12=nknots12,
                                                             ve01=ve01,
                                                             ve02=ve02,
                                                             ve12=ve12,
                                                             dimnva01=dimnva01,
                                                             dimnva02=dimnva02,
                                                             dimnva12=dimnva12,
                                                             nva01=nvat01,
                                                             nva02=nvat02,
                                                             nva12=nvat12,
                                                             t0=t0,
                                                             t1=t1,
                                                             t2=t2,
                                                             t3=t3,
                                                             troncature=troncature,
                                                             y01=y01k,
                                                             y02=y02k,
                                                             y12=y12k,
                                                             p01=p01,
                                                             p02=p02,
                                                             p12=p12,
                                                             dimp01=dimp01,
                                                             dimp02=dimp02,
                                                             dimp12=dimp12,
                                                             Ntime=NtimePoints,
                                                             lambda=lambda[id.lambda,],
                                                             alpha=alpha,
                                                             penalty.factor=penalty.factor,
                                                             penalty=penalty)
                                 
                                 
                                 # we want to maximise the loglik thus : 
                                 # we have issue if res is NA or if not higher than previous one 
                                 # if not better or do not exist need to readjust
                                 # value of beta 
                                 if(res %in%c(-1e9,1e9) | res < fn.value){
                                   
                                   th<-1e-5
                                   step<-log(1.5)
                                   delta<-output.cv$b-c(beta)
                                   
                                   maxt <- max(abs(delta)) 
                                   
                                   if(maxt == 0){
                                     vw <- th
                                   }else{
                                     vw <- th/maxt
                                   }
                                   if(ite>0){
                                     res.out.error <- list("old.b"=round(c(s,beta)),
                                                           "old.rl"=round(fn.value),
                                                           "old.ca"=round(eval.cv.beta[ite]),
                                                           "old.cb"=round(eval.cv.loglik[ite]))
                                   }else{
                                     res.out.error <- list("old.b"=round(c(s,beta)),
                                                           "old.rl"=round(fn.value),
                                                           "old.ca"=round(1),
                                                           "old.cb"=round(1))
                                   }
                                   # from mla package 
                                   sears<-searpas(vw=vw,
                                                  step=step,
                                                  b=beta,
                                                  delta=delta,
                                                  funcpa=gaussDYNidmlLikelihoodpena,
                                                  res.out.error=res.out.error,
                                                  npm=npm,
                                                  npar=size_V,
                                                  bfix=s,
                                                  fix=fix0,
                                                  zi01=knots01,
                                                  zi02=knots02,
                                                  zi12=knots12,
                                                  ctime=ctime,
                                                  no=N,
                                                  nz01=nknots01,
                                                  nz02=nknots02,
                                                  nz12=nknots12,
                                                  ve01=ve01,
                                                  ve02=ve02,
                                                  ve12=ve12,
                                                  dimnva01=dimnva01,
                                                  dimnva02=dimnva02,
                                                  dimnva12=dimnva12,
                                                  nva01=nvat01,
                                                  nva02=nvat02,
                                                  nva12=nvat12,
                                                  t0=t0,
                                                  t1=t1,
                                                  t2=t2,
                                                  t3=t3,
                                                  troncature=troncature,
                                                  y01=y01k,
                                                  y02=y02k,
                                                  y12=y12k,
                                                  p01=p01,
                                                  p02=p02,
                                                  p12=p12,
                                                  dimp01=dimp01,
                                                  dimp02=dimp02,
                                                  dimp12=dimp12,
                                                  Ntime=NtimePoints,
                                                  lambda=lambda[id.lambda,],
                                                  alpha=alpha,
                                                  penalty.factor=penalty.factor,
                                                  penalty=penalty)
                                   
                                   
                                   betanew<-beta+delta*sears$vw
                                   b<-c(s,betanew)
                                   
                                   
                                   res<-gaussDYNidmlLikelihoodpena(b=b,
                                                               npm=size_V,
                                                               npar=size_V,
                                                               bfix=1,
                                                               fix=rep(0,size_V),
                                                               zi01=knots01,
                                                               zi02=knots02,
                                                               zi12=knots12,
                                                               ctime=ctime,
                                                               no=N,
                                                               nz01=nknots01,
                                                               nz02=nknots02,
                                                               nz12=nknots12,
                                                               ve01=ve01,
                                                               ve02=ve02,
                                                               ve12=ve12,
                                                               dimnva01=dimnva01,
                                                               dimnva02=dimnva02,
                                                               dimnva12=dimnva12,
                                                               nva01=nvat01,
                                                               nva02=nvat02,
                                                               nva12=nvat12,
                                                               t0=t0,
                                                               t1=t1,
                                                               t2=t2,
                                                               t3=t3,
                                                               troncature=troncature,
                                                               y01=y01k,
                                                               y02=y02k,
                                                               y12=y12k,
                                                               p01=p01,
                                                               p02=p02,
                                                               p12=p12,
                                                               dimp01=dimp01,
                                                               dimp02=dimp02,
                                                               dimp12=dimp12,
                                                               Ntime=NtimePoints,
                                                               lambda=lambda[id.lambda,],
                                                               alpha=alpha,
                                                               penalty.factor=penalty.factor,
                                                               penalty=penalty)
                                 }
                                 # if not better or do not exist need to readjust
                                 # value of beta 
                                 if(res %in%c(-1e9,1e9) | any(is.infinite(c(s,betanew)))){
                                   
                                   ite<-ite+1
                                   validity<-F
                                   eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                   eval.validity[ite]<-validity
                                   # save output over iterations 
                                   
                                   pbr_compu<-3
                                   break
                                 }else{validity<-T}
                                 
                                 
                                 # betanew already include s
                                 b<-c(s,betanew)
                                 
                                 bfix<-b[fix0.beta==1]
                                 b<-b[fix0.beta==0]
                                 # update for modelPar
                                 
                                 output.mla<- marqLevAlg::mla(b=b,
                                                  fn=gaussDYNidmlLikelihood,
                                                  epsa=epsa,
                                                  epsb=epsb,
                                                  epsd=epsd,
                                                  maxiter=maxiter.pena,
                                                  minimize=F,
                                                  npm= npmspline,
                                                  npar=size_V,
                                                  bfix=bfix,
                                                  fix=fix0.beta,
                                                  zi01=knots01,
                                                  zi02=knots02,
                                                  zi12=knots12,
                                                  ctime=ctime,
                                                  no=N,
                                                  nz01=nknots01,
                                                  nz02=nknots02,
                                                  nz12=nknots12,
                                                  ve01=ve01,
                                                  ve02=ve02,
                                                  ve12=ve12,
                                                  dimnva01=dimnva01,
                                                  dimnva02=dimnva02,
                                                  dimnva12=dimnva12,
                                                  nva01=nvat01,
                                                  nva02=nvat02,
                                                  nva12=nvat12,
                                                  t0=t0,
                                                  t1=t1,
                                                  t2=t2,
                                                  t3=t3,
                                                  troncature=troncature,
                                                  y01=y01k,
                                                  y02=y02k,
                                                  y12=y12k,
                                                  p01=p01,
                                                  p02=p02,
                                                  p12=p12,
                                                  dimp01=dimp01,
                                                  dimp02=dimp02,
                                                  dimp12=dimp12,
                                                  Ntime=NtimePoints)
                                 
                                 # look at convergence for each lambda :
                                 # mla output is loglik
                                 # new values for splines:
                                 snew<-s
                                 snew[fix00[1:size_spline]==0]<-output.mla$b
                               
                                 
                                 
                                 if(nvat01>0){
                                   b01<-betanew[1:nvat01][penalty.factor[1:nvat01]==1]
                                   if(p01>0){
                                     b01<-c(b01,betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1])
                                   }
                                 }else{
                                   if(p01>0){
                                     b01<-betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1]
                                   }else{
                                     b01<-0
                                   }
                                 }
                                 
                                 if(nvat02>0){
                                   b02<-betanew[(nvat01+1):(nvat01+nvat02)][penalty.factor[(nvat01+1):(nvat01+nvat02)]==1]
                                   if(p02>0){
                                     b02<-c(b02,betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1])
                                   }
                                 }else{
                                   if(p02>0){
                                     b02<-betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1]
                                   }else{b02<-0}
                                 }
                                 
                                 if(nvat12>0){
                                   b12<-betanew[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)][penalty.factor[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)]==1]
                                   if(p12>0){
                                     b12<-c(b12,betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1])
                                   }
                                 }else{
                                   if(p12>0){
                                     b12<-betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1]
                                   }else{b12<-0}
                                 }
                                 
                                 # calculate loglik pen 
                                 if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
                                   fn.valuenew<-output.mla$fn.value-lambda[id.lambda,1]*alpha*sum(abs(b01))-lambda[id.lambda,1]*(1-alpha)*sum(b01*b01)
                                   fn.valuenew<-fn.valuenew-lambda[id.lambda,2]*alpha*sum(abs(b02))-lambda[id.lambda,2]*(1-alpha)*sum(b02*b02)
                                   fn.valuenew<-fn.valuenew-lambda[id.lambda,3]*alpha*sum(abs(b12))-lambda[id.lambda,3]*(1-alpha)*sum(b12*b12)
                                 }
                                 
                                 
                                 if(penalty=="mcp"){
                                   
                                   p01<-rep(alpha*lambda[id.lambda,1]*lambda[id.lambda,1]/2,length(b01))
                                   idbeta<-which(b01<=alpha*lambda[id.lambda,1])
                                   p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])-((b01[idbeta]*b01[idbeta])/2*alpha)
                                   
                                   p02<-rep(alpha*lambda[id.lambda,2]*lambda[id.lambda,2]/2,length(b02))
                                   idbeta<-which(b02<=alpha*lambda[id.lambda,2])
                                   p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])-((b02[idbeta]*b02[idbeta])/2*alpha)
                                   
                                   p12<-rep(alpha*lambda[id.lambda,3]*lambda[id.lambda,3]/2,length(b12))
                                   idbeta<-which(b12<=alpha*lambda[id.lambda,3])
                                   p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])-((b12[idbeta]*b12[idbeta])/2*alpha)
                                   
                                   fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                   
                                 }
                                 
                                 if(penalty=="scad"){
                                   
                                   p01<-rep((lambda[id.lambda,1]^2)*(alpha+1)/2,length(b01))
                                   idbeta<-which(b01<=lambda[id.lambda,1])
                                   p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])
                                   idbeta<-which(abs(b01)<lambda[id.lambda,1]*alpha)
                                   p01[idbeta]<-(2*alpha*lambda[id.lambda,1]*abs(b01[idbeta])-b01[idbeta]^2-lambda[id.lambda,1]^2)/(2*(alpha-1))
                                   
                                   p02<-rep((lambda[id.lambda,2]^2)*(alpha+1)/2,length(b02))
                                   idbeta<-which(b02<=lambda[id.lambda,2])
                                   p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])
                                   idbeta<-which(abs(b02)<lambda[id.lambda,2]*alpha)
                                   p02[idbeta]<-(2*alpha*lambda[id.lambda,2]*abs(b02[idbeta])-b02[idbeta]^2-lambda[id.lambda,2]^2)/(2*(alpha-1))
                                   
                                   p12<-rep((lambda[id.lambda,3]^2)*(alpha+1)/2,length(b12))
                                   idbeta<-which(b12<=lambda[id.lambda,3])
                                   p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])
                                   idbeta<-which(abs(b12)<lambda[id.lambda,3]*alpha)
                                   p12[idbeta]<-(2*alpha*lambda[id.lambda,3]*abs(b12[idbeta])-b12[idbeta]^2-lambda[id.lambda,3]^2)/(2*(alpha-1))
                                   
                                   fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                   
                                 }
                                 
                                 ite<-ite+1
                                 
                                 #check cv 
                                 eval.cv.spline[ite]<-sum((snew-s)^2)
                                 eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                 eval.cv.loglik[ite]<-abs((fn.valuenew-fn.value)/fn.value)
                                 eval.loglik[ite]<-fn.valuenew
                                 eval.validity[ite]<-validity
                                 
                                 s<-snew
                                 beta<-betanew
                                 fn.value<-fn.valuenew
                                 
                                 # eval.cv beta valid only if validity.param=T
                                 if(((eval.cv.beta[ite] + eval.cv.spline[ite])<epsa) & eval.cv.loglik[ite]<epsb & validity==T){
                                   converged<-T}
                                 
                                 
                               }
                               
                               if(maxiter<=ite & converged==F){
                                 istop<-2
                               }else{
                                 if(ite<=maxiter & converged==T){
                                   istop<-1 
                                   
                                   # need to recalculate second derivatives 
                                   # if converged 
                                   
                                   b<-c(s,beta)
                                   bfix<-b[fix0==1]
                                   b<-b[fix0==0]
                                   
                                   output<-DYNderiva(
                                     h=1e-8,
                                                  b=b,
                                                  npm=npm,
                                                  npar=size_V,
                                                  bfix=bfix,
                                                  fix=fix0,
                                                  zi01=knots01,
                                                  zi02=knots02,
                                                  zi12=knots12,
                                                  ctime=ctime,
                                                  no=N,
                                                  nz01=nknots01,
                                                  nz02=nknots02,
                                                  nz12=nknots12,
                                                  ve01=ve01,
                                                  ve02=ve02,
                                                  ve12=ve12,
                                                  dimnva01=dimnva01,
                                                  dimnva02=dimnva02,
                                                  dimnva12=dimnva12,
                                                  nva01=nvat01,
                                                  nva02=nvat02,
                                                  nva12=nvat12,
                                                  t0=t0,
                                                  t1=t1,
                                                  t2=t2,
                                                  t3=t3,
                                                  troncature=troncature,
                                                  y01=y01k,
                                                  y02=y02k,
                                                  y12=y12k,
                                                  p01=p01,
                                                  p02=p02,
                                                  p12=p12,
                                                  dimp01=dimp01,
                                                  dimp02=dimp02,
                                                  dimp12=dimp12,
                                                  Ntime=NtimePoints)
                                
                                 output<-output$v
                                 
                                 fu<-output[(min+1):length(output)]
                                 V<-matrix(0,nrow=npm,ncol=npm)
                                 V[upper.tri(V,diag=T)]<-output[1:min]
                                 V<-V+t(V)
                                 diag(V)<-diag(V)/2
                                 
                                 # deriva gives information matrix
                                   V0<-V
                                 }else{
                                   if(pbr_compu==1){istop<-3}
                                   if(pbr_compu==2){istop<-4}
                                   if(pbr_compu==3){istop<-5}
                                 }
                               }
                               
                               # if stop==1 we can give matrix of second derivatives 
                               
                               
                               combine<-combine+1
                               return(list(b=c(s,beta),
                                           H=V0,
                                           lambda=as.double(lambda[id.lambda,]),
                                           alpha=alpha,
                                           fn.value=ifelse(!exists("output.mla"),NA,output.mla$fn.value), # loglik
                                           fn.value.pena=fn.value, # penalised loglik
                                           ni=ite,
                                           ca.beta=eval.cv.beta,
                                           ca.spline=eval.cv.spline,
                                           ca.validity=eval.validity,
                                           cb=eval.loglik,
                                           istop=istop,
                                           combine=combine))
                             }
    }else{
      outputNsample<-foreach::foreach(id.lambda=1:nlambda,
                               .combine = combine_lambda_mla,
                               .errorhandling = "remove")%dopar%{
                                 
                                 # computation pbr 
                                 
                                 pbr_compu<-0
                                 
                                 if(idsample==1){
                                   beta<-beta.start
                                   s<-s.start
                                 }else{
                                   s<-outputNsample[[1]]$b[1:size_spline,id.lambda]
                                   beta<-outputNsample[[1]]$b[(size_spline+1):size_V,id.lambda]
                                 }
                                 
                                 converged<-F
                                 ite<-0
                                 # if beta not change do not need to recalculate weights 
                                 H<-T
                                 
                                 eval.cv.spline<-rep(NA,maxiter+1)
                                 eval.cv.beta<-rep(NA,maxiter+1)
                                 eval.cv.loglik<-rep(NA,maxiter+1)
                                 eval.loglik<-rep(NA,maxiter+1)
                                 eval.validity<-rep(NA,maxiter+1)
                                 
                                 
                                 while(converged==F & ite<=maxiter){
                                   
                                   
                                   b<-c(s,beta)
                                   bfix<-b[fix0==1]
                                   b<-b[fix0==0]
                                   # derivative of loglik
                                   
                                     
                                     output<-DYNderivadiag(
                                                   h=1e-8,
                                                    b=b,
                                              npm=npm,
                                              npar=size_V,
                                              bfix=bfix,
                                              fix=fix0,
                                              zi01=knots01,
                                              zi02=knots02,
                                              zi12=knots12,
                                              ctime=ctime,
                                              no=N,
                                              nz01=nknots01,
                                              nz02=nknots02,
                                              nz12=nknots12,
                                              ve01=ve01,
                                              ve02=ve02,
                                              ve12=ve12,
                                              dimnva01=dimnva01,
                                              dimnva02=dimnva02,
                                              dimnva12=dimnva12,
                                              nva01=nvat01,
                                              nva02=nvat02,
                                              nva12=nvat12,
                                              t0=t0,
                                              t1=t1,
                                              t2=t2,
                                              t3=t3,
                                              troncature=troncature,
                                              y01=y01k,
                                              y02=y02k,
                                              y12=y12k,
                                              p01=p01,
                                              p02=p02,
                                              p12=p12,
                                              dimp01=dimp01,
                                              dimp02=dimp02,
                                              dimp12=dimp12,
                                              Ntime=NtimePoints)
                                     output<-output$v
                                   
                                     if(ite==0){# loglik penalised
                                       fn.value<-gaussDYNidmlLikelihoodpena(b=b,
                                                                           npm=npm,
                                                                           npar=size_V,
                                                                           bfix=bfix,
                                                                           fix=fix0,
                                                                           zi01=knots01,
                                                                           zi02=knots02,
                                                                           zi12=knots12,
                                                                           ctime=ctime,
                                                                           no=N,
                                                                           nz01=nknots01,
                                                                           nz02=nknots02,
                                                                           nz12=nknots12,
                                                                           ve01=ve01,
                                                                           ve02=ve02,
                                                                           ve12=ve12,
                                                                           dimnva01=dimnva01,
                                                                           dimnva02=dimnva02,
                                                                           dimnva12=dimnva12,
                                                                           nva01=nvat01,
                                                                           nva02=nvat02,
                                                                           nva12=nvat12,
                                                                           t0=t0,
                                                                           t1=t1,
                                                                           t2=t2,
                                                                           t3=t3,
                                                                           troncature=troncature,
                                                                           y01=y01k,
                                                                           y02=y02k,
                                                                           y12=y12k,
                                                                           p01=p01,
                                                                           p02=p02,
                                                                           p12=p12,
                                                                           dimp01=dimp01,
                                                                           dimp02=dimp02,
                                                                           dimp12=dimp12,
                                                                           Ntime=NtimePoints,
                                                                           
                                                                           lambda=lambda[id.lambda,],
                                                                           alpha=alpha,
                                                                           penalty.factor=penalty.factor,
                                                                           penalty=penalty)
                                     }
                                     
                                     if(any(is.na(output))|any(output==Inf) |any(output==-Inf)){
                                       warning("Computational error for calculation of the hessian : division by 0 or Infinite value")
                                       if(ite==0){
                                        
                                         fu<-output[(min+1):length(output)]
                                         V<-matrix(0,nrow=npm,ncol=npm)
                                         diag(V)<-output[1:min]
                                         # deriva gives information matrix
                                         tr <- sum(diag(V))/npm
                                         V0<-V}
                                       ite<-ite+1
                                       pbr_compu<-1
                                       break
                                     }
                                     
                                     
                                     
                                     fu<-output[(min+1):length(output)]
                                     V<-matrix(0,nrow=npm,ncol=npm)
                                     diag(V)<-output[1:min]
                                     # deriva gives information matrix
                                     tr <- sum(diag(V))/npm
                                     V0<-V
                                     
                                     eigen.values<-diag(V)
                                     
                                     idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                     
                                     
                                     idpos0<-idpos
                                     
                                     ncount<-da<-ga<-0
                                     
                                     while(idpos != 0){
                                       
                                       if(ncount==0){ 
                                         ga <- 0.01
                                         da <- 1E-2
                                       }else{
                                         if(((ncount <= 3) | (ga >= 1)) ){
                                           da <- da * 5
                                         }else{# if ncount > 10 only update ga 
                                           ga <- ga * 5
                                           # do not put ga at 1 as no countmax otherwise infinite while 
                                           if(ga > 1) ga <- 1
                                         }
                                       }
                                       
                                       ncount <- ncount + 1
                                       
                                       diagV <- diag(V)
                                       # put abs (1-ga) better than 1-ga cause ga can now be >1
                                       diagV<-ifelse(diagV!=0,diagV+da*(abs((1.e0-ga))*abs(diagV)+ga*tr),
                                                     da*ga*tr)
                                       
                                       diag(V)<-diagV
                                       # if we have a convex log-vraisemblance in eta then :
                                       # all eigen  values of the hessienne are >0.
                                       
                                       if(sum(V==Inf)>0|sum(V==-Inf)>0){break}
                                       eigen.values<-diag(V)
                                       # check if hessienne defined positive
                                       
                                       idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                       
                                       # if(def.positive==T){
                                       #   idpos<-ifelse(any(eigen.values<=0),1,0)
                                       # }else{idpos<-ifelse(any(abs(eigen.values)==0),1,0)}
                                       
                                     }
                                     
                                     if(idpos!=0){
                                       
                                       warning("Hessian not defined positive")
                                       pbr_compu<-2
                                       ite<-ite+1
                                       break
                                     }
                                     
                                     
                                   
                                   
                                
                                   output.cv<-DYNcv.model(beta=beta,
                                                          nva01=npm01,
                                                          nva02=npm02,
                                                          nva12=npm12,
                                                          nva01Y=npm01Y,
                                                          nva02Y=npm02Y,
                                                          nva12Y=npm12Y,
                                                          fix=fix0[(size_spline+1):size_V],
                                                          penalty.factor=penalty.factor,
                                                          penalty=penalty,
                                                          v=V,
                                                          fu=fu,
                                                          lambda=lambda[id.lambda,],
                                                          alpha=alpha
                                   )
                                   
                                   # verify validity of parameters update 
                                   # and that we are better than previous estimates 
                                   
                                   b<-c(s,output.cv$b)
                                   
                                   betanew<-b[(size_spline+1):size_V]
                                   
                                   # penalised loglik see if inferior to previous
                                   res<-gaussDYNidmlLikelihoodpena(b=b,
                                                                  npm=size_V,
                                                                  npar=size_V,
                                                                  bfix=1,
                                                                  fix=rep(0,size_V),
                                                                  zi01=knots01,
                                                                  zi02=knots02,
                                                                  zi12=knots12,
                                                                  ctime=ctime,
                                                                  no=N,
                                                                  nz01=nknots01,
                                                                  nz02=nknots02,
                                                                  nz12=nknots12,
                                                                  ve01=ve01,
                                                                  ve02=ve02,
                                                                  ve12=ve12,
                                                                  dimnva01=dimnva01,
                                                                  dimnva02=dimnva02,
                                                                  dimnva12=dimnva12,
                                                                  nva01=nvat01,
                                                                  nva02=nvat02,
                                                                  nva12=nvat12,
                                                                  t0=t0,
                                                                  t1=t1,
                                                                  t2=t2,
                                                                  t3=t3,
                                                                  troncature=troncature,
                                                                  y01=y01k,
                                                                  y02=y02k,
                                                                  y12=y12k,
                                                                  p01=p01,
                                                                  p02=p02,
                                                                  p12=p12,
                                                                  dimp01=dimp01,
                                                                  dimp02=dimp02,
                                                                  dimp12=dimp12,
                                                                  Ntime=NtimePoints,
                                                                  
                                                                  lambda=lambda[id.lambda,],
                                                                  alpha=alpha,
                                                                  penalty.factor=penalty.factor,
                                                                  penalty=penalty)
                                   
                                   
                                   # we want to maximise the loglik thus : 
                                   # we have issue if res is NA or if not higher than previous one 
                                   # if not better or do not exist need to readjust
                                   # value of beta 
                                   if(res %in%c(-1e9,1e9) | res < fn.value){
                                     
                                     th<-1e-5
                                     step<-log(1.5)
                                     delta<-output.cv$b-c(beta)
                                     
                                     maxt <- max(abs(delta)) 
                                     
                                     if(maxt == 0){
                                       vw <- th
                                     }else{
                                       vw <- th/maxt
                                     }
                                     if(ite>0){
                                       res.out.error <- list("old.b"=round(c(s,beta)),
                                                             "old.rl"=round(fn.value),
                                                             "old.ca"=round(eval.cv.beta[ite]),
                                                             "old.cb"=round(eval.cv.loglik[ite]))
                                     }else{
                                       res.out.error <- list("old.b"=round(c(s,beta)),
                                                             "old.rl"=round(fn.value),
                                                             "old.ca"=round(1),
                                                             "old.cb"=round(1))
                                     }
                                     # from mla package 
                                     sears<-searpas(vw=vw,
                                                    step=step,
                                                    b=beta,
                                                    delta=delta,
                                                    funcpa=gaussDYNidmlLikelihoodpena,
                                                    res.out.error=res.out.error,
                                                    npm=npm,
                                                    npar=size_V,
                                                    bfix=s,
                                                    fix=fix0,
                                                    zi01=knots01,
                                                    zi02=knots02,
                                                    zi12=knots12,
                                                    ctime=ctime,
                                                    no=N,
                                                    nz01=nknots01,
                                                    nz02=nknots02,
                                                    nz12=nknots12,
                                                    ve01=ve01,
                                                    ve02=ve02,
                                                    ve12=ve12,
                                                    dimnva01=dimnva01,
                                                    dimnva02=dimnva02,
                                                    dimnva12=dimnva12,
                                                    nva01=nvat01,
                                                    nva02=nvat02,
                                                    nva12=nvat12,
                                                    t0=t0,
                                                    t1=t1,
                                                    t2=t2,
                                                    t3=t3,
                                                    troncature=troncature,
                                                    y01=y01k,
                                                    y02=y02k,
                                                    y12=y12k,
                                                    p01=p01,
                                                    p02=p02,
                                                    p12=p12,
                                                    dimp01=dimp01,
                                                    dimp02=dimp02,
                                                    dimp12=dimp12,
                                                    Ntime=NtimePoints,
                                                    lambda=lambda[id.lambda,],
                                                    alpha=alpha,
                                                    penalty.factor=penalty.factor,
                                                    penalty=penalty)
                                     
                                     betanew<-beta+delta*sears$vw
                                     b<-c(s,betanew)
                                     
                                     
                                     res<-gaussDYNidmlLikelihoodpena(b=b,
                                                                    npm=size_V,
                                                                    npar=size_V,
                                                                    bfix=1,
                                                                    fix=rep(0,size_V),
                                                                    zi01=knots01,
                                                                    zi02=knots02,
                                                                    zi12=knots12,
                                                                    ctime=ctime,
                                                                    no=N,
                                                                    nz01=nknots01,
                                                                    nz02=nknots02,
                                                                    nz12=nknots12,
                                                                    ve01=ve01,
                                                                    ve02=ve02,
                                                                    ve12=ve12,
                                                                    dimnva01=dimnva01,
                                                                    dimnva02=dimnva02,
                                                                    dimnva12=dimnva12,
                                                                    nva01=nvat01,
                                                                    nva02=nvat02,
                                                                    nva12=nvat12,
                                                                    t0=t0,
                                                                    t1=t1,
                                                                    t2=t2,
                                                                    t3=t3,
                                                                    troncature=troncature,
                                                                    y01=y01k,
                                                                    y02=y02k,
                                                                    y12=y12k,
                                                                    p01=p01,
                                                                    p02=p02,
                                                                    p12=p12,
                                                                    dimp01=dimp01,
                                                                    dimp02=dimp02,
                                                                    dimp12=dimp12,
                                                                    Ntime=NtimePoints,
                                                                  
                                                                    lambda=lambda[id.lambda,],
                                                                    alpha=alpha,
                                                                    penalty.factor=penalty.factor,
                                                                    penalty=penalty)
                                   }
                                   # if not better or do not exist need to readjust
                                   # value of beta 
                                   if(res %in%c(-1e9,1e9) | any(is.infinite(c(s,betanew)))){
                                     
                                     ite<-ite+1
                                     validity<-F
                                     eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                     eval.validity[ite]<-validity
                                     # save output over iterations 
                                     
                                     pbr_compu<-3
                                     break
                                   }else{validity<-T}
                                   
                                   
                                   # betanew already include s
                                   b<-c(s,betanew)
                                   
                                   bfix<-b[fix0.beta==1]
                                   b<-b[fix0.beta==0]
                                   # update for modelPar
                                  
                                   output.mla<- marqLevAlg::mla(b=b,
                                                                fn=gaussDYNidmlLikelihood,
                                                                epsa=epsa,
                                                                epsb=epsb,
                                                                epsd=epsd,
                                                                maxiter=maxiter.pena,
                                                                minimize=F,
                                                                npm=npmspline,
                                                                npar=size_V,
                                                                bfix=bfix,
                                                                fix=fix0.beta,
                                                                zi01=knots01,
                                                                zi02=knots02,
                                                                zi12=knots12,
                                                                ctime=ctime,
                                                                no=N,
                                                                nz01=nknots01,
                                                                nz02=nknots02,
                                                                nz12=nknots12,
                                                                ve01=ve01,
                                                                ve02=ve02,
                                                                ve12=ve12,
                                                                dimnva01=dimnva01,
                                                                dimnva02=dimnva02,
                                                                dimnva12=dimnva12,
                                                                nva01=nvat01,
                                                                nva02=nvat02,
                                                                nva12=nvat12,
                                                                t0=t0,
                                                                t1=t1,
                                                                t2=t2,
                                                                t3=t3,
                                                                troncature=troncature,
                                                                y01=y01k,
                                                                y02=y02k,
                                                                y12=y12k,
                                                                p01=p01,
                                                                p02=p02,
                                                                p12=p12,
                                                                dimp01=dimp01,
                                                                dimp02=dimp02,
                                                                dimp12=dimp12,
                                                                Ntime=NtimePoints)

                                   # look at convergence for each lambda :
                                   # mla output is loglik
                                   # new values for splines:
                                   snew<-s
                                   snew[fix00[1:size_spline]==0]<-output.mla$b
                                   if(nvat01>0){
                                     b01<-betanew[1:nvat01][penalty.factor[1:nvat01]==1]
                                     if(p01>0){
                                       b01<-c(b01,betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1])
                                     }
                                   }else{
                                     if(p01>0){
                                       b01<-betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1]
                                     }else{
                                       b01<-0
                                     }
                                   }
                                   
                                   if(nvat02>0){
                                     b02<-betanew[(nvat01+1):(nvat01+nvat02)][penalty.factor[(nvat01+1):(nvat01+nvat02)]==1]
                                     if(p02>0){
                                       b02<-c(b02,betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1])
                                     }
                                   }else{
                                     if(p02>0){
                                       b02<-betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1]
                                     }else{b02<-0}
                                   }
                                   
                                   if(nvat12>0){
                                     b12<-betanew[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)][penalty.factor[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)]==1]
                                     if(p12>0){
                                       b12<-c(b12,betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1])
                                     }
                                   }else{
                                     if(p12>0){
                                       b12<-betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1]
                                     }else{b12<-0}
                                   }
                                   
                                   # calculate loglik pen 
                                   if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
                                     fn.valuenew<-output.mla$fn.value-lambda[id.lambda,1]*alpha*sum(abs(b01))-lambda[id.lambda,1]*(1-alpha)*sum(b01*b01)
                                     fn.valuenew<-fn.valuenew-lambda[id.lambda,2]*alpha*sum(abs(b02))-lambda[id.lambda,2]*(1-alpha)*sum(b02*b02)
                                     fn.valuenew<-fn.valuenew-lambda[id.lambda,3]*alpha*sum(abs(b12))-lambda[id.lambda,3]*(1-alpha)*sum(b12*b12)
                                   }
                                   
                                   
                                   if(penalty=="mcp"){
                                     
                                     p01<-rep(alpha*lambda[id.lambda,1]*lambda[id.lambda,1]/2,length(b01))
                                     idbeta<-which(b01<=alpha*lambda[id.lambda,1])
                                     p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])-((b01[idbeta]*b01[idbeta])/2*alpha)
                                     
                                     p02<-rep(alpha*lambda[id.lambda,2]*lambda[id.lambda,2]/2,length(b02))
                                     idbeta<-which(b02<=alpha*lambda[id.lambda,2])
                                     p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])-((b02[idbeta]*b02[idbeta])/2*alpha)
                                     
                                     p12<-rep(alpha*lambda[id.lambda,3]*lambda[id.lambda,3]/2,length(b12))
                                     idbeta<-which(b12<=alpha*lambda[id.lambda,3])
                                     p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])-((b12[idbeta]*b12[idbeta])/2*alpha)
                                     
                                     fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                     
                                   }
                                   
                                   if(penalty=="scad"){
                                     
                                     p01<-rep((lambda[id.lambda,1]^2)*(alpha+1)/2,length(b01))
                                     idbeta<-which(b01<=lambda[id.lambda,1])
                                     p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])
                                     idbeta<-which(abs(b01)<lambda[id.lambda,1]*alpha)
                                     p01[idbeta]<-(2*alpha*lambda[id.lambda,1]*abs(b01[idbeta])-b01[idbeta]^2-lambda[id.lambda,1]^2)/(2*(alpha-1))
                                     
                                     p02<-rep((lambda[id.lambda,2]^2)*(alpha+1)/2,length(b02))
                                     idbeta<-which(b02<=lambda[id.lambda,2])
                                     p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])
                                     idbeta<-which(abs(b02)<lambda[id.lambda,2]*alpha)
                                     p02[idbeta]<-(2*alpha*lambda[id.lambda,2]*abs(b02[idbeta])-b02[idbeta]^2-lambda[id.lambda,2]^2)/(2*(alpha-1))
                                     
                                     p12<-rep((lambda[id.lambda,3]^2)*(alpha+1)/2,length(b12))
                                     idbeta<-which(b12<=lambda[id.lambda,3])
                                     p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])
                                     idbeta<-which(abs(b12)<lambda[id.lambda,3]*alpha)
                                     p12[idbeta]<-(2*alpha*lambda[id.lambda,3]*abs(b12[idbeta])-b12[idbeta]^2-lambda[id.lambda,3]^2)/(2*(alpha-1))
                                     
                                     fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                     
                                   }
                                   
                                   ite<-ite+1
                                   
                                   #check cv 
                                   eval.cv.spline[ite]<-sum((snew-s)^2)
                                   eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                   eval.cv.loglik[ite]<-abs((fn.valuenew-fn.value)/fn.value)
                                   eval.loglik[ite]<-fn.valuenew
                                   eval.validity[ite]<-validity
                                   
                                   s<-snew
                                   beta<-betanew
                                   fn.value<-fn.valuenew
                                   
                                   # eval.cv beta valid only if validity.param=T
                                   if(((eval.cv.beta[ite] + eval.cv.spline[ite])<epsa) & eval.cv.loglik[ite]<epsb & validity==T){
                                     converged<-T}
                                   
                                   
                                 }
                                 
                                 if(maxiter<=ite & converged==F){
                                   istop<-2
                                 }else{
                                   if(ite<=maxiter & converged==T){
                                     istop<-1 
                                     
                                     # need to recalculate second derivatives 
                                     # if converged 
                                     
                                     b<-c(s,beta)
                                     bfix<-b[fix0==1]
                                     b<-b[fix0==0]
                                     
                                     output<-DYNderivadiag(
                                                         h=1e-8,
                                                        b=b,
                                                        npm=npm,
                                                        npar=size_V,
                                                        bfix=bfix,
                                                        fix=fix0,
                                                        zi01=knots01,
                                                        zi02=knots02,
                                                        zi12=knots12,
                                                        ctime=ctime,
                                                        no=N,
                                                        nz01=nknots01,
                                                        nz02=nknots02,
                                                        nz12=nknots12,
                                                        ve01=ve01,
                                                        ve02=ve02,
                                                        ve12=ve12,
                                                        dimnva01=dimnva01,
                                                        dimnva02=dimnva02,
                                                        dimnva12=dimnva12,
                                                        nva01=nvat01,
                                                        nva02=nvat02,
                                                        nva12=nvat12,
                                                        t0=t0,
                                                        t1=t1,
                                                        t2=t2,
                                                        t3=t3,
                                                        troncature=troncature,
                                                        y01=y01k,
                                                        y02=y02k,
                                                        y12=y12k,
                                                        p01=p01,
                                                        p02=p02,
                                                        p12=p12,
                                                        dimp01=dimp01,
                                                        dimp02=dimp02,
                                                        dimp12=dimp12,
                                                        Ntime=NtimePoints)
                                     
                                      output<-output$v
                                      fu<-output[(npm+1):length(output)]
                                      V<-matrix(0,nrow=npm,ncol=npm)
                                      diag(V)<-output[1:npm]
                                      
                                     # deriva gives information matrix
                                     V0<-V
                                   }else{
                                     if(pbr_compu==1){istop<-3}
                                     if(pbr_compu==2){istop<-4}
                                     if(pbr_compu==3){istop<-5}
                                   }
                                 }
                                 
                                 # if stop==1 we can give matrix of second derivatives 
                                 
                                 
                                 combine<-combine+1
                                 return(list(b=c(s,beta),
                                             H=V0,
                                             lambda=as.double(lambda[id.lambda,]),
                                             alpha=alpha,
                                             fn.value=ifelse(!exists("output.mla"),NA,output.mla$fn.value), # loglik
                                             fn.value.pena=fn.value, # penalised loglik
                                             ni=ite,
                                             ca.beta=eval.cv.beta,
                                             ca.spline=eval.cv.spline,
                                             ca.validity=eval.validity,
                                             cb=eval.loglik,
                                             istop=istop,
                                             combine=combine))
                               }
    }
    
    parallel::stopCluster(clustpar)
    
  }else{
    
  
    id.lambda<-NULL # for cran check 
    if(partialH==F){
      
    outputNsample<-foreach::foreach(id.lambda=1:nlambda,
                             .combine = combine_lambda_mla,
                             .errorhandling = "remove")%do%{
                               
                               
                               # computation pbr 
                               
                               pbr_compu<-0
                               
                               if(idsample==1){
                                 beta<-beta.start
                                 s<-s.start
                               }else{
                                 s<-output[[1]]$b[1:size_spline,id.lambda]
                                 beta<-output[[1]]$b[(size_spline+1):size_V,id.lambda]
                               }
                               
                               converged<-F
                               ite<-0
                               # if beta not change do not need to recalculate weights 
                               H<-T
                               
                               eval.cv.spline<-rep(NA,maxiter+1)
                               eval.cv.beta<-rep(NA,maxiter+1)
                               eval.cv.loglik<-rep(NA,maxiter+1)
                               eval.loglik<-rep(NA,maxiter+1)
                               eval.validity<-rep(NA,maxiter+1)
                               
                              
                               while(converged==F & ite<=maxiter){
                               
                                
                                 b<-c(s,beta)
                                 bfix<-b[fix0==1]
                                 b<-b[fix0==0]
                               
                                 
                                 output<-DYNderiva(
                                                h=1e-8,
                                                b=b,
                                                npm=npm,
                                                npar=size_V,
                                                bfix=bfix,
                                                fix=fix0,
                                                zi01=knots01,
                                                zi02=knots02,
                                                zi12=knots12,
                                                ctime=ctime,
                                                no=N,
                                                nz01=nknots01,
                                                nz02=nknots02,
                                                nz12=nknots12,
                                                        ve01=ve01,
                                                        ve02=ve02,
                                                        ve12=ve12,
                                                        dimnva01=dimnva01,
                                                        dimnva02=dimnva02,
                                                        dimnva12=dimnva12,
                                                        nva01=nvat01,
                                                        nva02=nvat02,
                                                        nva12=nvat12,
                                                        t0=t0,
                                                        t1=t1,
                                                        t2=t2,
                                                        t3=t3,
                                                        troncature=troncature,
                                                        y01=y01k,
                                                        y02=y02k,
                                                        y12=y12k,
                                                        p01=p01,
                                                        p02=p02,
                                                        p12=p12,
                                                        dimp01=dimp01,
                                                        dimp02=dimp02,
                                                        dimp12=dimp12,
                                                        Ntime=NtimePoints)
                                 
                                 
                                
                                 #sparseHessianFD
                                 # h=1e-16 --> problem in hessian evaluation
                                 # h=1e-8 --> no issue
                                 output<-output$v
                                
                                 
                               if(ite==0){
                                   fn.value<-gaussDYNidmlLikelihoodpena(b=b,
                                                                       npm=npm,
                                                                       npar=size_V,
                                                                       bfix=bfix,
                                                                       fix=fix0,
                                                                       zi01=knots01,
                                                                       zi02=knots02,
                                                                       zi12=knots12,
                                                                       ctime=ctime,
                                                                       no=N,
                                                                       nz01=nknots01,
                                                                       nz02=nknots02,
                                                                       nz12=nknots12,
                                                                       ve01=ve01,
                                                                       ve02=ve02,
                                                                       ve12=ve12,
                                                                       dimnva01=dimnva01,
                                                                       dimnva02=dimnva02,
                                                                       dimnva12=dimnva12,
                                                                       nva01=nvat01,
                                                                       nva02=nvat02,
                                                                       nva12=nvat12,
                                                                       t0=t0,
                                                                       t1=t1,
                                                                       t2=t2,
                                                                       t3=t3,
                                                                       troncature=troncature,
                                                                       y01=y01k,
                                                                       y02=y02k,
                                                                       y12=y12k,
                                                                       p01=p01,
                                                                       p02=p02,
                                                                       p12=p12,
                                                                       dimp01=dimp01,
                                                                       dimp02=dimp02,
                                                                       dimp12=dimp12,
                                                                       Ntime=NtimePoints,
                                                                       
                                                                       lambda=lambda[id.lambda,],
                                                                       alpha=alpha,
                                                                       penalty.factor=penalty.factor,
                                                                       penalty=penalty)
                               }
                                 
                                
                                   if(any(is.na(output))|any(output==Inf) |any(output==-Inf)){
                                     warning("Computational error for calculation of the hessian : division by 0 or Infinite value")
                                     if(ite==0){
                                       
                                       fu<-output[(min+1):length(output)]
                                       V<-matrix(0,nrow=npm,ncol=npm)
                                       V[upper.tri(V,diag=T)]<-output[1:min]
                                       V<-V+t(V)
                                       diag(V)<-diag(V)/2
                                       # deriva gives information matrix
                                       tr <- sum(diag(V))/npm
                                       V0<-V}
                                     ite<-ite+1
                                     pbr_compu<-1
                                     break
                                   }
                                  
                                 
                                 fu<-output[(min+1):length(output)]
                                 V<-matrix(0,nrow=npm,ncol=npm)
                                 V[upper.tri(V,diag=T)]<-output[1:min]
                                 V<-V+t(V)
                                 diag(V)<-diag(V)/2
                                 
                                 # projection V into space of eigen vectors
                                 # tol <- 1e-12
                                 # eig <- eigen(V, symmetric = TRUE)
                                 # vals <- eig$values
                                 # vecs <- eig$vectors
                                 # 
                                 # # tronquer les valeurs propres trop petites
                                 # vals_clipped <- pmax(vals, tol)
                                 # 
                                 # # reconstruire la matrice projetée
                                 # H_proj <- vecs %*% diag(vals_clipped) %*% t(vecs)
                                 # 
                                 # 
                                 # pp<-matrix(0,nrow=npm,ncol=npm)
                                 # diag(pp)<-diag(V)
                                   # deriva gives information matrix
                                   tr <- sum(diag(V))/npm
                                   V0<-V
                                
                                 
                                 eigen.values<-eigen(V,symmetric=T,only.values=T)$values
                                
                                 idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                 
                                 
                                 idpos0<-idpos
                                 
                                 ncount<-da<-ga<-0
                              
                                 while(idpos != 0){
                                   
                                   if(ncount==0){ 
                                     ga <- 0.01
                                     da <- 1E-2
                                   }else{
                                     if(((ncount <= 3) | (ga >= 1)) ){
                                       da <- da * 5
                                     }else{# if ncount > 10 only update ga 
                                       ga <- ga * 5
                                       # do not put ga at 1 as no countmax otherwise infinite while 
                                       if(ga > 1) ga <- 1
                                     }
                                   }
                                   
                                   ncount <- ncount + 1
                                   
                                   diagV <- diag(V)
                                   # put abs (1-ga) better than 1-ga cause ga can now be >1
                                   diagV<-ifelse(diagV!=0,diagV+da*(abs((1.e0-ga))*abs(diagV)+ga*tr),
                                                 da*ga*tr)
                                   
                                   diag(V)<-diagV
                                   # if we have a convex log-vraisemblance in eta then :
                                   # all eigen  values of the hessienne are >0.
                                   
                                   if(sum(V==Inf)>0|sum(V==-Inf)>0){break}
                                   eigen.values<-eigen(V,symmetric=T,only.values=T)$values
                                   # check if hessienne defined positive
                                   idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                   
                                   
                                   
                                 }
                                 
                               
                                 if(idpos!=0){
                                  
                                   warning("Hessian not defined positive")
                                   pbr_compu<-2
                                   ite<-ite+1
                                   break
                                 }
                                 
                            
                                 # update beta
                                
                                 output.cv<-DYNcv.model(beta=beta,
                                                        nva01=npm01,
                                                        nva02=npm02,
                                                        nva12=npm12,
                                                        nva01Y=npm01Y,
                                                        nva02Y=npm02Y,
                                                        nva12Y=npm12Y,
                                                        fix=fix0[(size_spline+1):size_V],
                                                        penalty.factor=penalty.factor,
                                                        penalty=penalty,
                                                        v=V,
                                                        fu=fu,
                                                        lambda=lambda[id.lambda,],
                                                        alpha=alpha
                                 )
                                 
                                 # verify validity of parameters update 
                                 # and that we are better than previous estimates 
                               
                                 b<-c(s,output.cv$b)
                                 
                                 betanew<-b[(size_spline+1):size_V]
                                 
                                 # penalised loglik see if inferior to previous
                                 res<-gaussDYNidmlLikelihoodpena(b=b,
                                                                npm=size_V,
                                                                npar=size_V,
                                                                bfix=1,
                                                                fix=rep(0,size_V),
                                                                zi01=knots01,
                                                                zi02=knots02,
                                                                zi12=knots12,
                                                                ctime=ctime,
                                                                no=N,
                                                                nz01=nknots01,
                                                                nz02=nknots02,
                                                                nz12=nknots12,
                                                                ve01=ve01,
                                                                ve02=ve02,
                                                                ve12=ve12,
                                                                dimnva01=dimnva01,
                                                                dimnva02=dimnva02,
                                                                dimnva12=dimnva12,
                                                                nva01=nvat01,
                                                                nva02=nvat02,
                                                                nva12=nvat12,
                                                                t0=t0,
                                                                t1=t1,
                                                                t2=t2,
                                                                t3=t3,
                                                                troncature=troncature,
                                                                y01=y01k,
                                                                y02=y02k,
                                                                y12=y12k,
                                                                p01=p01,
                                                                p02=p02,
                                                                p12=p12,
                                                                dimp01=dimp01,
                                                                dimp02=dimp02,
                                                                dimp12=dimp12,
                                                                Ntime=NtimePoints,
                                                                lambda=lambda[id.lambda,],
                                                                alpha=alpha,
                                                                penalty.factor=penalty.factor,
                                                                penalty=penalty)
                                 
                                 
                                 # we want to maximise the loglik thus : 
                                 # we have issue if res is NA or if not higher than previous one 
                                 # if not better or do not exist need to readjust
                                 # value of beta 
                                
                              
                                if(res %in%c(-1e9,1e9) | res < fn.value){
                                  
                                  print(paste0("needed update at ite :",ite))
                                   th<-1e-5
                                   step<-log(1.5)
                                   delta<-output.cv$b-c(beta)
                                   
                                   maxt <- max(abs(delta)) 
                                   
                                   if(maxt == 0){
                                     vw <- th
                                   }else{
                                     vw <- th/maxt
                                   }
                                   if(ite>0){
                                   res.out.error <- list("old.b"=round(c(s,beta)),
                                                         "old.rl"=round(fn.value),
                                                         "old.ca"=round(eval.cv.beta[ite]),
                                                         "old.cb"=round(eval.cv.loglik[ite]))
                                   }else{
                                     res.out.error <- list("old.b"=round(c(s,beta)),
                                                           "old.rl"=round(fn.value),
                                                           "old.ca"=round(1),
                                                           "old.cb"=round(1))
                                   }
                                   m<-1
                                   betanew<-beta
                                   # from mla package
                                   
                                   sears<-searpas(vw=vw,
                                                  step=step,
                                                  b=beta,
                                                  delta=delta,
                                                  funcpa=gaussDYNidmlLikelihoodpena,
                                                  res.out.error=res.out.error,
                                                  npm=npm,
                                                  npar=size_V,
                                                  bfix=s,
                                                  fix=fix0,
                                                  zi01=knots01,
                                                  zi02=knots02,
                                                  zi12=knots12,
                                                  ctime=ctime,
                                                  no=N,
                                                  nz01=nknots01,
                                                  nz02=nknots02,
                                                  nz12=nknots12,
                                                  ve01=ve01,
                                                  ve02=ve02,
                                                  ve12=ve12,
                                                  dimnva01=dimnva01,
                                                  dimnva02=dimnva02,
                                                  dimnva12=dimnva12,
                                                  nva01=nvat01,
                                                  nva02=nvat02,
                                                  nva12=nvat12,
                                                  t0=t0,
                                                  t1=t1,
                                                  t2=t2,
                                                  t3=t3,
                                                  troncature=troncature,
                                                  y01=y01k,
                                                  y02=y02k,
                                                  y12=y12k,
                                                  p01=p01,
                                                  p02=p02,
                                                  p12=p12,
                                                  dimp01=dimp01,
                                                  dimp02=dimp02,
                                                  dimp12=dimp12,
                                                  Ntime=NtimePoints,
                                                  lambda=lambda[id.lambda,],
                                                  alpha=alpha,
                                                  penalty.factor=penalty.factor,
                                                  penalty=penalty)
                                   
                                   
                                   betanew<-beta+delta*sears$vw
                                   b<-c(s,betanew)
                                   
                                   
                                   res<-gaussDYNidmlLikelihoodpena(b=b,
                                                                  npm=size_V,
                                                                  npar=size_V,
                                                                  bfix=1,
                                                                  fix=rep(0,size_V),
                                                                  zi01=knots01,
                                                                  zi02=knots02,
                                                                  zi12=knots12,
                                                                  ctime=ctime,
                                                                  no=N,
                                                                  nz01=nknots01,
                                                                  nz02=nknots02,
                                                                  nz12=nknots12,
                                                                  ve01=ve01,
                                                                  ve02=ve02,
                                                                  ve12=ve12,
                                                                  dimnva01=dimnva01,
                                                                  dimnva02=dimnva02,
                                                                  dimnva12=dimnva12,
                                                                  nva01=nvat01,
                                                                  nva02=nvat02,
                                                                  nva12=nvat12,
                                                                  t0=t0,
                                                                  t1=t1,
                                                                  t2=t2,
                                                                  t3=t3,
                                                                  troncature=troncature,
                                                                  y01=y01k,
                                                                  y02=y02k,
                                                                  y12=y12k,
                                                                  p01=p01,
                                                                  p02=p02,
                                                                  p12=p12,
                                                                  dimp01=dimp01,
                                                                  dimp02=dimp02,
                                                                  dimp12=dimp12,
                                                                  Ntime=NtimePoints,
                                                                  lambda=lambda[id.lambda,],
                                                                  alpha=alpha,
                                                                  penalty.factor=penalty.factor,
                                                                  penalty=penalty)
                                   
                                 }
                                 # if not better or do not exist need to readjust
                                 # value of beta 
                                if(res %in%c(-1e9,1e9) | any(is.infinite(c(s,betanew)))){
                                  
                                     ite<-ite+1
                                     validity<-F
                                     eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                     eval.validity[ite]<-validity
                                     # save output over iterations 
                                     pbr_compu<-3
                                     break
                                   }else{validity<-T}
                                   
                                   
                                 # betanew already include s
                                 b<-c(s,betanew)
                                 
                                 bfix<-b[fix0.beta==1]
                                 b<-b[fix0.beta==0]
                                 # update modelPar
                              
                                
                                 output.mla<- marqLevAlg::mla(b=b,
                                                              fn=gaussDYNidmlLikelihood,
                                                              epsa=epsa,
                                                              epsb=epsb,
                                                              epsd=epsd,
                                                              maxiter=maxiter.pena,
                                                              minimize=F,
                                                              npm=npmspline,
                                                              npar=size_V,
                                                              bfix=bfix,
                                                              fix=fix0.beta,
                                                              zi01=knots01,
                                                              zi02=knots02,
                                                              zi12=knots12,
                                                              ctime=ctime,
                                                              no=N,
                                                              nz01=nknots01,
                                                              nz02=nknots02,
                                                              nz12=nknots12,
                                                              ve01=ve01,
                                                              ve02=ve02,
                                                              ve12=ve12,
                                                              dimnva01=dimnva01,
                                                              dimnva02=dimnva02,
                                                              dimnva12=dimnva12,
                                                              nva01=nvat01,
                                                              nva02=nvat02,
                                                              nva12=nvat12,
                                                              t0=t0,
                                                              t1=t1,
                                                              t2=t2,
                                                              t3=t3,
                                                              troncature=troncature,
                                                              y01=y01k,
                                                              y02=y02k,
                                                              y12=y12k,
                                                              p01=p01,
                                                              p02=p02,
                                                              p12=p12,
                                                              dimp01=dimp01,
                                                              dimp02=dimp02,
                                                              dimp12=dimp12,
                                                              Ntime=NtimePoints)
                                 
                                 # look at convergence for each lambda :
                                 
                                 # new values for splines:
                                 snew<-s
                                 snew[fix00[1:size_spline]==0]<-output.mla$b
                                 if(nvat01>0){
                                   b01<-betanew[1:nvat01][penalty.factor[1:nvat01]==1]
                                   if(p01>0){
                                     b01<-c(b01,betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1])
                                   }
                                 }else{
                                   if(p01>0){
                                     b01<-betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1]
                                   }else{
                                     b01<-0
                                   }
                                 }
                                 
                                 if(nvat02>0){
                                   b02<-betanew[(nvat01+1):(nvat01+nvat02)][penalty.factor[(nvat01+1):(nvat01+nvat02)]==1]
                                   if(p02>0){
                                     b02<-c(b02,betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1])
                                   }
                                 }else{
                                   if(p02>0){
                                     b02<-betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1]
                                   }else{b02<-0}
                                 }
                                 
                                 if(nvat12>0){
                                   b12<-betanew[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)][penalty.factor[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)]==1]
                                   if(p12>0){
                                     b12<-c(b12,betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1])
                                   }
                                 }else{
                                   if(p12>0){
                                     b12<-betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1]
                                   }else{b12<-0}
                                 }
                                 # maximisation issue : lpen =l - pen
                                 if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
                                   fn.valuenew<-output.mla$fn.value-lambda[id.lambda,1]*alpha*sum(abs(b01))-lambda[id.lambda,1]*(1-alpha)*sum(b01*b01)
                                   fn.valuenew<-fn.valuenew-lambda[id.lambda,2]*alpha*sum(abs(b02))-lambda[id.lambda,2]*(1-alpha)*sum(b02*b02)
                                   fn.valuenew<-fn.valuenew-lambda[id.lambda,3]*alpha*sum(abs(b12))-lambda[id.lambda,3]*(1-alpha)*sum(b12*b12)

                                   }
                                 
                                 if(penalty=="mcp"){
                                   
                                   p01<-rep(alpha*lambda[id.lambda,1]*lambda[id.lambda,1]/2,length(b01))
                                   idbeta<-which(b01<=alpha*lambda[id.lambda,1])
                                   p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])-((b01[idbeta]*b01[idbeta])/2*alpha)
                                   
                                   p02<-rep(alpha*lambda[id.lambda,2]*lambda[id.lambda,2]/2,length(b02))
                                   idbeta<-which(b02<=alpha*lambda[id.lambda,2])
                                   p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])-((b02[idbeta]*b02[idbeta])/2*alpha)
                                   
                                   p12<-rep(alpha*lambda[id.lambda,3]*lambda[id.lambda,3]/2,length(b12))
                                   idbeta<-which(b12<=alpha*lambda[id.lambda,3])
                                   p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])-((b12[idbeta]*b12[idbeta])/2*alpha)
                                   
                                   fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                   
                                 }
                                 
                                 if(penalty=="scad"){
                                   
                                   p01<-rep((lambda[id.lambda,1]^2)*(alpha+1)/2,length(b01))
                                   idbeta<-which(b01<=lambda[id.lambda,1])
                                   p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])
                                   idbeta<-which(abs(b01)<lambda[id.lambda,1]*alpha)
                                   p01[idbeta]<-(2*alpha*lambda[id.lambda,1]*abs(b01[idbeta])-b01[idbeta]^2-lambda[id.lambda,1]^2)/(2*(alpha-1))
                                   
                                   p02<-rep((lambda[id.lambda,2]^2)*(alpha+1)/2,length(b02))
                                   idbeta<-which(b02<=lambda[id.lambda,2])
                                   p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])
                                   idbeta<-which(abs(b02)<lambda[id.lambda,2]*alpha)
                                   p02[idbeta]<-(2*alpha*lambda[id.lambda,2]*abs(b02[idbeta])-b02[idbeta]^2-lambda[id.lambda,2]^2)/(2*(alpha-1))
                                   
                                   p12<-rep((lambda[id.lambda,3]^2)*(alpha+1)/2,length(b12))
                                   idbeta<-which(b12<=lambda[id.lambda,3])
                                   p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])
                                   idbeta<-which(abs(b12)<lambda[id.lambda,3]*alpha)
                                   p12[idbeta]<-(2*alpha*lambda[id.lambda,3]*abs(b12[idbeta])-b12[idbeta]^2-lambda[id.lambda,3]^2)/(2*(alpha-1))
                                   
                                   fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                   
                                 }
                                 
                                 ite<-ite+1
                                 
                                 # check cv criterias 
                               
                                 eval.cv.spline[ite]<-sum((snew-s)^2)
                                 eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                 eval.cv.loglik[ite]<-abs((fn.valuenew-fn.value)/fn.value)
                                 eval.loglik[ite]<-fn.valuenew
                                 eval.validity[ite]<-validity
                                 
                                 s<-snew
                                 beta<-betanew
                                 fn.value<-fn.valuenew
                               
                                 # eval.cv beta valid only if validity.param=T
                                 # add ite >1 so that do not stop at first ite ?
                                 if(((eval.cv.beta[ite]+eval.cv.spline[ite])<epsa) & eval.cv.loglik[ite]<epsb & validity==T){
                                   converged<-T}
                                 
                                 
                               }
                               
                               if(maxiter<=ite & converged==F){
                                 istop<-2
                               }else{
                                 if(ite<=maxiter & converged==T){
                                   istop<-1
                                   # need to recalculate second derivatives 
                                   
                                   b<-c(s,beta)
                                   bfix<-b[fix0==1]
                                   b<-b[fix0==0]
                                   
                                   output<-DYNderiva(
                                                   h=1e-8,
                                                  b=b,
                                                  npm=npm,
                                                  npar=size_V,
                                                  bfix=bfix,
                                                  fix=fix0,
                                                  zi01=knots01,
                                                  zi02=knots02,
                                                  zi12=knots12,
                                                  ctime=ctime,
                                                  no=N,
                                                  nz01=nknots01,
                                                  nz02=nknots02,
                                                  nz12=nknots12,
                                                  ve01=ve01,
                                                  ve02=ve02,
                                                  ve12=ve12,
                                                  dimnva01=dimnva01,
                                                  dimnva02=dimnva02,
                                                  dimnva12=dimnva12,
                                                  nva01=nvat01,
                                                  nva02=nvat02,
                                                  nva12=nvat12,
                                                  t0=t0,
                                                  t1=t1,
                                                  t2=t2,
                                                  t3=t3,
                                                  troncature=troncature,
                                                  y01=y01k,
                                                  y02=y02k,
                                                  y12=y12k,
                                                  p01=p01,
                                                  p02=p02,
                                                  p12=p12,
                                                  dimp01=dimp01,
                                                  dimp02=dimp02,
                                                  dimp12=dimp12,
                                                  Ntime=NtimePoints)
                                   
                                   output<-output$v
                                   fu<-output[(min+1):length(output)]
                                   V<-matrix(0,nrow=npm,ncol=npm)
                                   V[upper.tri(V,diag=T)]<-output[1:min]
                                   V<-V+t(V)
                                   diag(V)<-diag(V)/2
                                   # deriva gives information matrix
                                   V0<-V
                                 }else{
                                   if(pbr_compu==1){istop<-3}
                                   if(pbr_compu==2){istop<-4}
                                   if(pbr_compu==3){istop<-5}
                                 }
                               }
                               
                               # if stop==1 we can give matrix of second derivatives 
                             
                               combine<-combine+1
                               return(list(b=c(s,beta),
                                           H=V0,
                                           lambda=as.double(lambda[id.lambda,]),
                                           alpha=alpha,
                                           fn.value=ifelse(!exists("output.mla"),NA,output.mla$fn.value),
                                           fn.value.pena=fn.value,
                                           ni=ite,
                                           ca.beta=eval.cv.beta,
                                           ca.spline=eval.cv.spline,
                                           ca.validity=eval.validity,
                                           cb=eval.loglik,
                                           istop=istop,
                                           combine=combine))
                             }
    }else{outputNsample<-foreach::foreach(id.lambda=1:nlambda,
                                   .combine = combine_lambda_mla,
                                   .errorhandling = "remove")%do%{
                                     
                                     
                                     # computation pbr 
                                     
                                     pbr_compu<-0
                                     
                                     if(idsample==1){
                                       beta<-beta.start
                                       s<-s.start
                                     }else{
                                       s<-outputNsample[[1]]$b[1:size_spline,id.lambda]
                                       beta<-outputNsample[[1]]$b[(size_spline+1):size_V,id.lambda]
                                     }
                                     
                                     converged<-F
                                     ite<-0
                                     # if beta not change do not need to recalculate weights 
                                     H<-T
                                     
                                     eval.cv.spline<-rep(NA,maxiter+1)
                                     eval.cv.beta<-rep(NA,maxiter+1)
                                     eval.cv.loglik<-rep(NA,maxiter+1)
                                     eval.loglik<-rep(NA,maxiter+1)
                                     eval.validity<-rep(NA,maxiter+1)
                                     
                                     
                                     while(converged==F & ite<=maxiter){
                                      
                                       b<-c(s,beta)
                                       bfix<-b[fix0==1]
                                       b<-b[fix0==0]
                                       
                                 
                                         
                                         output<-DYNderivadiag(
                                                        h=1e-8,
                                                        b=b,
                                                                   npm=npm,
                                                                   npar=size_V,
                                                                   bfix=bfix,
                                                                   fix=fix0,
                                                        zi01=knots01,
                                                        zi02=knots02,
                                                        zi12=knots12,
                                                        ctime=ctime,
                                                        no=N,
                                                        nz01=nknots01,
                                                        nz02=nknots02,
                                                        nz12=nknots12,
                                                                   ve01=ve01,
                                                                   ve02=ve02,
                                                                   ve12=ve12,
                                                                   dimnva01=dimnva01,
                                                                   dimnva02=dimnva02,
                                                                   dimnva12=dimnva12,
                                                                   nva01=nvat01,
                                                                   nva02=nvat02,
                                                                   nva12=nvat12,
                                                                   t0=t0,
                                                                   t1=t1,
                                                                   t2=t2,
                                                                   t3=t3,
                                                                   troncature=troncature,
                                                                   y01=y01k,
                                                                   y02=y02k,
                                                                   y12=y12k,
                                                                   p01=p01,
                                                                   p02=p02,
                                                                   p12=p12,
                                                                   dimp01=dimp01,
                                                                   dimp02=dimp02,
                                                                   dimp12=dimp12,
                                                                   Ntime=NtimePoints)
                                         
                                       
                                         output<-output$v
                                         if(ite==0){
                                           fn.value<-gaussDYNidmlLikelihoodpena(b=b,
                                                                               npm=npm,
                                                                               npar=size_V,
                                                                               bfix=bfix,
                                                                               fix=fix0,
                                                                               zi01=knots01,
                                                                               zi02=knots02,
                                                                               zi12=knots12,
                                                                               ctime=ctime,
                                                                               no=N,
                                                                               nz01=nknots01,
                                                                               nz02=nknots02,
                                                                               nz12=nknots12,
                                                                               ve01=ve01,
                                                                               ve02=ve02,
                                                                               ve12=ve12,
                                                                               dimnva01=dimnva01,
                                                                               dimnva02=dimnva02,
                                                                               dimnva12=dimnva12,
                                                                               nva01=nvat01,
                                                                               nva02=nvat02,
                                                                               nva12=nvat12,
                                                                               t0=t0,
                                                                               t1=t1,
                                                                               t2=t2,
                                                                               t3=t3,
                                                                               troncature=troncature,
                                                                               y01=y01k,
                                                                               y02=y02k,
                                                                               y12=y12k,
                                                                               p01=p01,
                                                                               p02=p02,
                                                                               p12=p12,
                                                                               dimp01=dimp01,
                                                                               dimp02=dimp02,
                                                                               dimp12=dimp12,
                                                                               Ntime=NtimePoints,
                                                                               
                                                                               lambda=lambda[id.lambda,],
                                                                               alpha=alpha,
                                                                               penalty.factor=penalty.factor,
                                                                               penalty=penalty)
                                         }
                                         
                                         if(any(is.na(output$V))|any(output$V==Inf) |any(output$V==-Inf)|any(is.na(output$fu))|any(output$fu==Inf) |any(output$fu==-Inf)){
                                           warning("Computational error for calculation of the hessian : division by 0 or Infinite value")
                                           if(ite==0){
                                             
                                             fu<-output[(min+1):(2*min)]
                                             V<-matrix(0,nrow=npm,ncol=npm)
                                             diag(V)<-output[1:min]
                                             # deriva gives information matrix
                                             tr <- sum(diag(V))/npm
                                             V0<-V}
                                           ite<-ite+1
                                           pbr_compu<-1
                                           break
                                         }
                                         
                                         
                                         fu<-output[(min+1):(2*min)]
                                         V<-matrix(0,nrow=npm,ncol=npm)
                                         diag(V)<-output[1:min]
                                         # deriva gives information matrix
                                         tr <- sum(diag(V))/npm
                                         V0<-V
                                         
                                         eigen.values<-diag(V)
                                         
                                         idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                         
                                         
                                         idpos0<-idpos
                                         
                                         ncount<-da<-ga<-0
                                         
                                         
                                         while(idpos != 0){
                                           
                                           if(ncount==0){ 
                                             ga <- 0.01
                                             da <- 1E-2
                                           }else{
                                             if(((ncount <= 3) | (ga >= 1)) ){
                                               da <- da * 5
                                             }else{# if ncount > 10 only update ga 
                                               ga <- ga * 5
                                               # do not put ga at 1 as no countmax otherwise infinite while 
                                               if(ga > 1) ga <- 1
                                             }
                                           }
                                           
                                           ncount <- ncount + 1
                                           
                                           diagV <- diag(V)
                                           # put abs (1-ga) better than 1-ga cause ga can now be >1
                                           diagV<-ifelse(diagV!=0,diagV+da*(abs((1.e0-ga))*abs(diagV)+ga*tr),
                                                         da*ga*tr)
                                           
                                           diag(V)<-diagV
                                           # if we have a convex log-vraisemblance in eta then :
                                           # all eigen  values of the hessienne are >0.
                                           
                                           if(sum(V==Inf)>0|sum(V==-Inf)>0){break}
                                           eigen.values<-diag(V)
                                           # check if hessienne defined positive
                                           idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                           
                                           
                                           
                                         }
                                         
                                         if(idpos!=0){
                                           
                                           warning("Hessian not defined positive")
                                           pbr_compu<-2
                                           ite<-ite+1
                                           break
                                         }
                                         
                                       
                                       # update beta
                                       output.cv<-DYNcv.model(beta=beta,
                                                              nva01=npm01,
                                                              nva02=npm02,
                                                              nva12=npm12,
                                                              nva01Y=npm01Y,
                                                              nva02Y=npm02Y,
                                                              nva12Y=npm12Y,
                                                              fix=fix0[(size_spline+1):size_V],
                                                              penalty.factor=penalty.factor,
                                                              penalty=penalty,
                                                              v=V,
                                                              fu=fu,
                                                              lambda=lambda[id.lambda,],
                                                              alpha=alpha
                                       )
                                       
                                       # verify validity of parameters update 
                                       # and that we are better than previous estimates 
                                       
                                       b<-c(s,output.cv$b)
                                       
                                       betanew<-b[(size_spline+1):size_V]
                                       
                                       # penalised loglik see if inferior to previous
                                       res<-gaussDYNidmlLikelihoodpena(b=b,
                                                                      npm=size_V,
                                                                      npar=size_V,
                                                                      bfix=1,
                                                                      fix=rep(0,size_V),
                                                                      zi01=knots01,
                                                                      zi02=knots02,
                                                                      zi12=knots12,
                                                                      ctime=ctime,
                                                                      no=N,
                                                                      nz01=nknots01,
                                                                      nz02=nknots02,
                                                                      nz12=nknots12,
                                                                      ve01=ve01,
                                                                      ve02=ve02,
                                                                      ve12=ve12,
                                                                      dimnva01=dimnva01,
                                                                      dimnva02=dimnva02,
                                                                      dimnva12=dimnva12,
                                                                      nva01=nvat01,
                                                                      nva02=nvat02,
                                                                      nva12=nvat12,
                                                                      t0=t0,
                                                                      t1=t1,
                                                                      t2=t2,
                                                                      t3=t3,
                                                                      troncature=troncature,
                                                                      y01=y01k,
                                                                      y02=y02k,
                                                                      y12=y12k,
                                                                      p01=p01,
                                                                      p02=p02,
                                                                      p12=p12,
                                                                      dimp01=dimp01,
                                                                      dimp02=dimp02,
                                                                      dimp12=dimp12,
                                                                      Ntime=NtimePoints,
                                                                     
                                                                      lambda=lambda[id.lambda,],
                                                                      alpha=alpha,
                                                                      penalty.factor=penalty.factor,
                                                                      penalty=penalty)
                                       
                                       # we want to maximise the loglik thus : 
                                       # we have issue if res is NA or if not higher than previous one 
                                       # if not better or do not exist need to readjust
                                       # value of beta 
                                    
                                       if(res %in%c(-1e9,1e9) | res < fn.value){
                                         
                                         th<-1e-5
                                         step<-log(1.5)
                                         delta<-output.cv$b-c(beta)
                                         
                                         maxt <- max(abs(delta)) 
                                         
                                         if(maxt == 0){
                                           vw <- th
                                         }else{
                                           vw <- th/maxt
                                         }
                                         if(ite>0){
                                           res.out.error <- list("old.b"=round(c(s,beta)),
                                                                 "old.rl"=round(fn.value),
                                                                 "old.ca"=round(eval.cv.beta[ite]),
                                                                 "old.cb"=round(eval.cv.loglik[ite]))
                                         }else{
                                           res.out.error <- list("old.b"=round(c(s,beta)),
                                                                 "old.rl"=round(fn.value),
                                                                 "old.ca"=round(1),
                                                                 "old.cb"=round(1))
                                         }
                                         m<-1
                                         betanew<-beta
                                         # from mla package
                                         
                                         sears<-searpas(vw=vw,
                                                        step=step,
                                                        b=beta,
                                                        delta=delta,
                                                        funcpa=gaussDYNidmlLikelihoodpena,
                                                        res.out.error=res.out.error,
                                                        npm=npm,
                                                        npar=size_V,
                                                        bfix=s,
                                                        fix=fix0,
                                                        zi01=knots01,
                                                        zi02=knots02,
                                                        zi12=knots12,
                                                        ctime=ctime,
                                                        no=N,
                                                        nz01=nknots01,
                                                        nz02=nknots02,
                                                        nz12=nknots12,
                                                        ve01=ve01,
                                                        ve02=ve02,
                                                        ve12=ve12,
                                                        dimnva01=dimnva01,
                                                        dimnva02=dimnva02,
                                                        dimnva12=dimnva12,
                                                        nva01=nvat01,
                                                        nva02=nvat02,
                                                        nva12=nvat12,
                                                        t0=t0,
                                                        t1=t1,
                                                        t2=t2,
                                                        t3=t3,
                                                        troncature=troncature,
                                                        y01=y01k,
                                                        y02=y02k,
                                                        y12=y12k,
                                                        p01=p01,
                                                        p02=p02,
                                                        p12=p12,
                                                        dimp01=dimp01,
                                                        dimp02=dimp02,
                                                        dimp12=dimp12,
                                                        Ntime=NtimePoints,
                                                        lambda=lambda[id.lambda,],
                                                        alpha=alpha,
                                                        penalty.factor=penalty.factor,
                                                        penalty=penalty)
                                         
                                         
                                         betanew<-beta+delta*sears$vw
                                         b<-c(s,betanew)
                                         
                                         
                                         res<-gaussDYNidmlLikelihoodpena(b=b,
                                                                        npm=size_V,
                                                                        npar=size_V,
                                                                        bfix=1,
                                                                        fix=rep(0,size_V),
                                                                        zi01=knots01,
                                                                        zi02=knots02,
                                                                        zi12=knots12,
                                                                        ctime=ctime,
                                                                        no=N,
                                                                        nz01=nknots01,
                                                                        nz02=nknots02,
                                                                        nz12=nknots12,
                                                                        ve01=ve01,
                                                                        ve02=ve02,
                                                                        ve12=ve12,
                                                                        dimnva01=dimnva01,
                                                                        dimnva02=dimnva02,
                                                                        dimnva12=dimnva12,
                                                                        nva01=nvat01,
                                                                        nva02=nvat02,
                                                                        nva12=nvat12,
                                                                        t0=t0,
                                                                        t1=t1,
                                                                        t2=t2,
                                                                        t3=t3,
                                                                        troncature=troncature,
                                                                        y01=y01k,
                                                                        y02=y02k,
                                                                        y12=y12k,
                                                                        p01=p01,
                                                                        p02=p02,
                                                                        p12=p12,
                                                                        dimp01=dimp01,
                                                                        dimp02=dimp02,
                                                                        dimp12=dimp12,
                                                                        Ntime=NtimePoints,
                                                                       
                                                                        lambda=lambda[id.lambda,],
                                                                        alpha=alpha,
                                                                        penalty.factor=penalty.factor,
                                                                        penalty=penalty)
                                         
                                       }
                                       # if not better or do not exist need to readjust
                                       # value of beta 
                                       if(res %in%c(-1e9,1e9) | any(is.infinite(c(s,betanew)))){
                                         
                                         ite<-ite+1
                                         validity<-F
                                         eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                         eval.validity[ite]<-validity
                                         # save output over iterations 
                                         pbr_compu<-3
                                         break
                                       }else{validity<-T}
                                       
                                       
                                       # betanew already include s
                                       b<-c(s,betanew)
                                       
                                       bfix<-b[fix0.beta==1]
                                       b<-b[fix0.beta==0]
                                       # update modelPar
                                       
                                       
                                       output.mla<- marqLevAlg::mla(b=b,
                                                                    fn=gaussDYNidmlLikelihood,
                                                                    epsa=epsa,
                                                                    epsb=epsb,
                                                                    epsd=epsd,
                                                                    maxiter=maxiter.pena,
                                                                    minimize=F,
                                                                    npm=npmspline,
                                                                    npar=size_V,
                                                                    bfix=bfix,
                                                                    fix=fix0.beta,
                                                                    zi01=knots01,
                                                                    zi02=knots02,
                                                                    zi12=knots12,
                                                                    ctime=ctime,
                                                                    no=N,
                                                                    nz01=nknots01,
                                                                    nz02=nknots02,
                                                                    nz12=nknots12,
                                                                    ve01=ve01,
                                                                    ve02=ve02,
                                                                    ve12=ve12,
                                                                    dimnva01=dimnva01,
                                                                    dimnva02=dimnva02,
                                                                    dimnva12=dimnva12,
                                                                    nva01=nvat01,
                                                                    nva02=nvat02,
                                                                    nva12=nvat12,
                                                                    t0=t0,
                                                                    t1=t1,
                                                                    t2=t2,
                                                                    t3=t3,
                                                                    troncature=troncature,
                                                                    y01=y01k,
                                                                    y02=y02k,
                                                                    y12=y12k,
                                                                    p01=p01,
                                                                    p02=p02,
                                                                    p12=p12,
                                                                    dimp01=dimp01,
                                                                    dimp02=dimp02,
                                                                    dimp12=dimp12,
                                                                    Ntime=NtimePoints)
                                       
                                       # look at convergence for each lambda :
                                       
                                       # new values for splines:
                                       snew<-s
                                       snew[fix00[1:size_spline]==0]<-output.mla$b
                                       if(nvat01>0){
                                         b01<-betanew[1:nvat01][penalty.factor[1:nvat01]==1]
                                         if(p01>0){
                                           b01<-c(b01,betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1])
                                         }
                                       }else{
                                         if(p01>0){
                                           b01<-betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1]
                                         }else{
                                           b01<-0
                                         }
                                       }
                                       
                                       if(nvat02>0){
                                         b02<-betanew[(nvat01+1):(nvat01+nvat02)][penalty.factor[(nvat01+1):(nvat01+nvat02)]==1]
                                         if(p02>0){
                                           b02<-c(b02,betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1])
                                         }
                                       }else{
                                         if(p02>0){
                                           b02<-betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1]
                                         }else{b02<-0}
                                       }
                                       
                                       if(nvat12>0){
                                         b12<-betanew[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)][penalty.factor[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)]==1]
                                         if(p12>0){
                                           b12<-c(b12,betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1])
                                         }
                                       }else{
                                         if(p12>0){
                                           b12<-betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1]
                                         }else{b12<-0}
                                       }
                                       # maximisation issue : lpen =l - pen
                                       if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
                                         fn.valuenew<-output.mla$fn.value-lambda[id.lambda,1]*alpha*sum(abs(b01))-lambda[id.lambda,1]*(1-alpha)*sum(b01*b01)
                                         fn.valuenew<-fn.valuenew-lambda[id.lambda,2]*alpha*sum(abs(b02))-lambda[id.lambda,2]*(1-alpha)*sum(b02*b02)
                                         fn.valuenew<-fn.valuenew-lambda[id.lambda,3]*alpha*sum(abs(b12))-lambda[id.lambda,3]*(1-alpha)*sum(b12*b12)
                                         
                                       }
                                       
                                       if(penalty=="mcp"){
                                         
                                         p01<-rep(alpha*lambda[id.lambda,1]*lambda[id.lambda,1]/2,length(b01))
                                         idbeta<-which(b01<=alpha*lambda[id.lambda,1])
                                         p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])-((b01[idbeta]*b01[idbeta])/2*alpha)
                                         
                                         p02<-rep(alpha*lambda[id.lambda,2]*lambda[id.lambda,2]/2,length(b02))
                                         idbeta<-which(b02<=alpha*lambda[id.lambda,2])
                                         p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])-((b02[idbeta]*b02[idbeta])/2*alpha)
                                         
                                         p12<-rep(alpha*lambda[id.lambda,3]*lambda[id.lambda,3]/2,length(b12))
                                         idbeta<-which(b12<=alpha*lambda[id.lambda,3])
                                         p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])-((b12[idbeta]*b12[idbeta])/2*alpha)
                                         
                                         fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                         
                                       }
                                       
                                       if(penalty=="scad"){
                                         
                                         p01<-rep((lambda[id.lambda,1]^2)*(alpha+1)/2,length(b01))
                                         idbeta<-which(b01<=lambda[id.lambda,1])
                                         p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])
                                         idbeta<-which(abs(b01)<lambda[id.lambda,1]*alpha)
                                         p01[idbeta]<-(2*alpha*lambda[id.lambda,1]*abs(b01[idbeta])-b01[idbeta]^2-lambda[id.lambda,1]^2)/(2*(alpha-1))
                                         
                                         p02<-rep((lambda[id.lambda,2]^2)*(alpha+1)/2,length(b02))
                                         idbeta<-which(b02<=lambda[id.lambda,2])
                                         p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])
                                         idbeta<-which(abs(b02)<lambda[id.lambda,2]*alpha)
                                         p02[idbeta]<-(2*alpha*lambda[id.lambda,2]*abs(b02[idbeta])-b02[idbeta]^2-lambda[id.lambda,2]^2)/(2*(alpha-1))
                                         
                                         p12<-rep((lambda[id.lambda,3]^2)*(alpha+1)/2,length(b12))
                                         idbeta<-which(b12<=lambda[id.lambda,3])
                                         p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])
                                         idbeta<-which(abs(b12)<lambda[id.lambda,3]*alpha)
                                         p12[idbeta]<-(2*alpha*lambda[id.lambda,3]*abs(b12[idbeta])-b12[idbeta]^2-lambda[id.lambda,3]^2)/(2*(alpha-1))
                                         
                                         fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                         
                                       }
                                       
                                       ite<-ite+1
                                       
                                       # check cv criterias 
                                       
                                       eval.cv.spline[ite]<-sum((snew-s)^2)
                                       eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                       eval.cv.loglik[ite]<-abs((fn.valuenew-fn.value)/fn.value)
                                       eval.loglik[ite]<-fn.valuenew
                                       eval.validity[ite]<-validity
                                       
                                       s<-snew
                                       beta<-betanew
                                       fn.value<-fn.valuenew
                                       
                                       # eval.cv beta valid only if validity.param=T
                                       if(((eval.cv.beta[ite]+ eval.cv.spline[ite])<epsa) & eval.cv.loglik[ite]<epsb & validity==T){
                                         converged<-T}
                                       
                                       
                                     }
                                     if(maxiter<=ite & converged==F){
                                       istop<-2
                                     }else{
                                       if(ite<=maxiter & converged==T){
                                         istop<-1
                                         # need to recalculate second derivatives 
                                         
                                         b<-c(s,beta)
                                         bfix<-b[fix0==1]
                                         b<-b[fix0==0]
                                         
                                         output<-DYNderivadiag(
                                                        h=1e-8,
                                                        b=b,
                                                        npm=npm,
                                                        npar=size_V,
                                                        bfix=bfix,
                                                        fix=fix0,
                                                        zi01=knots01,
                                                        zi02=knots02,
                                                        zi12=knots12,
                                                        ctime=ctime,
                                                        no=N,
                                                        nz01=nknots01,
                                                        nz02=nknots02,
                                                        nz12=nknots12,
                                                        ve01=ve01,
                                                        ve02=ve02,
                                                        ve12=ve12,
                                                        dimnva01=dimnva01,
                                                        dimnva02=dimnva02,
                                                        dimnva12=dimnva12,
                                                        nva01=nvat01,
                                                        nva02=nvat02,
                                                        nva12=nvat12,
                                                        t0=t0,
                                                        t1=t1,
                                                        t2=t2,
                                                        t3=t3,
                                                        troncature=troncature,
                                                        y01=y01k,
                                                        y02=y02k,
                                                        y12=y12k,
                                                        p01=p01,
                                                        p02=p02,
                                                        p12=p12,
                                                        dimp01=dimp01,
                                                        dimp02=dimp02,
                                                        dimp12=dimp12,
                                                        Ntime=NtimePoints)
                                        
                                         output<-output$v
                                         
                                         
                                         fu<-output[(min+1):(2*min)]
                                         V<-matrix(0,nrow=npm,ncol=npm)
                                         diag(V)<-output[1:min]
                                         # deriva gives information matrix
                                         V0<-V
                                       }else{
                                         if(pbr_compu==1){istop<-3}
                                         if(pbr_compu==2){istop<-4}
                                         if(pbr_compu==3){istop<-5}
                                       }
                                     }
                                     
                                     # if stop==1 we can give matrix of second derivatives 
                                     
                                     combine<-combine+1
                                    
                                     return(list(b=c(s,beta),
                                                 H=V0,
                                                 lambda=as.double(lambda[id.lambda,]),
                                                 alpha=alpha,
                                                 fn.value=ifelse(!exists("output.mla"),NA,output.mla$fn.value),
                                                 fn.value.pena=fn.value,
                                                 ni=ite,
                                                 ca.beta=eval.cv.beta,
                                                 ca.spline=eval.cv.spline,
                                                 ca.validity=eval.validity,
                                                 cb=eval.loglik,
                                                 istop=istop,
                                                 combine=combine))
                                   }}
  }
    
    outputall[[idsample]]<-outputNsample
  }
    
    
  }else{ #can do parallelisation on the sample 
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
                              assoc=assoc,
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
          
          
        }else{
          y01k<-rep(0,N*NtimePoints)
        }
        
        if(length(outcome02)>=1){
          y02k<-dataY[dataY$Outcome%in%outcome02,]
          # order  by individual and timeline 
          y02k<-y02k[order(y02k[,colnames(y02k)%in%id],y02k$order),]
          
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
        
        
          id.lambda<-NULL # for cran check 
          if(partialH==F){
            
            outputNsample<-foreach::foreach(id.lambda=1:nlambda,
                                     .combine = combine_lambda_mla,
                                     .errorhandling = "remove")%do%{
                                       
                                       
                                       # computation pbr 
                                       
                                       pbr_compu<-0
                                       
                                       if(idsample==1){
                                         beta<-beta.start
                                         s<-s.start
                                       }else{
                                         s<-output[[1]]$b[1:size_spline,id.lambda]
                                         beta<-output[[1]]$b[(size_spline+1):size_V,id.lambda]
                                       }
                                       
                                       converged<-F
                                       ite<-0
                                       # if beta not change do not need to recalculate weights 
                                       H<-T
                                       
                                       eval.cv.spline<-rep(NA,maxiter+1)
                                       eval.cv.beta<-rep(NA,maxiter+1)
                                       eval.cv.loglik<-rep(NA,maxiter+1)
                                       eval.loglik<-rep(NA,maxiter+1)
                                       eval.validity<-rep(NA,maxiter+1)
                                       
                                       
                                       while(converged==F & ite<=maxiter){
                                         
                                         
                                         b<-c(s,beta)
                                         bfix<-b[fix0==1]
                                         b<-b[fix0==0]
                                         
                                      
                                         output<-DYNderiva(
                                           h=1e-8,
                                           b=b,
                                           npm=npm,
                                           npar=size_V,
                                           bfix=bfix,
                                           fix=fix0,
                                           zi01=knots01,
                                           zi02=knots02,
                                           zi12=knots12,
                                           ctime=ctime,
                                           no=N,
                                           nz01=nknots01,
                                           nz02=nknots02,
                                           nz12=nknots12,
                                           ve01=ve01,
                                           ve02=ve02,
                                           ve12=ve12,
                                           dimnva01=dimnva01,
                                           dimnva02=dimnva02,
                                           dimnva12=dimnva12,
                                           nva01=nvat01,
                                           nva02=nvat02,
                                           nva12=nvat12,
                                           t0=t0,
                                           t1=t1,
                                           t2=t2,
                                           t3=t3,
                                           troncature=troncature,
                                           y01=y01k,
                                           y02=y02k,
                                           y12=y12k,
                                           p01=p01,
                                           p02=p02,
                                           p12=p12,
                                           dimp01=dimp01,
                                           dimp02=dimp02,
                                           dimp12=dimp12,
                                           Ntime=NtimePoints)
                                         
                                         
                                         
                                         #sparseHessianFD
                                         # h=1e-16 --> problem in hessian evaluation
                                         # h=1e-8 --> no issue
                                         output<-output$v
                                         
                                         
                                         if(ite==0){
                                           fn.value<-gaussDYNidmlLikelihoodpena(b=b,
                                                                                npm=npm,
                                                                                npar=size_V,
                                                                                bfix=bfix,
                                                                                fix=fix0,
                                                                                zi01=knots01,
                                                                                zi02=knots02,
                                                                                zi12=knots12,
                                                                                ctime=ctime,
                                                                                no=N,
                                                                                nz01=nknots01,
                                                                                nz02=nknots02,
                                                                                nz12=nknots12,
                                                                                ve01=ve01,
                                                                                ve02=ve02,
                                                                                ve12=ve12,
                                                                                dimnva01=dimnva01,
                                                                                dimnva02=dimnva02,
                                                                                dimnva12=dimnva12,
                                                                                nva01=nvat01,
                                                                                nva02=nvat02,
                                                                                nva12=nvat12,
                                                                                t0=t0,
                                                                                t1=t1,
                                                                                t2=t2,
                                                                                t3=t3,
                                                                                troncature=troncature,
                                                                                y01=y01k,
                                                                                y02=y02k,
                                                                                y12=y12k,
                                                                                p01=p01,
                                                                                p02=p02,
                                                                                p12=p12,
                                                                                dimp01=dimp01,
                                                                                dimp02=dimp02,
                                                                                dimp12=dimp12,
                                                                                Ntime=NtimePoints,
                                                                                
                                                                                lambda=lambda[id.lambda,],
                                                                                alpha=alpha,
                                                                                penalty.factor=penalty.factor,
                                                                                penalty=penalty)
                                         }
                                         
                                         
                                         if(any(is.na(output))|any(output==Inf) |any(output==-Inf)){
                                           warning("Computational error for calculation of the hessian : division by 0 or Infinite value")
                                           if(ite==0){
                                             
                                             fu<-output[(min+1):length(output)]
                                             V<-matrix(0,nrow=npm,ncol=npm)
                                             V[upper.tri(V,diag=T)]<-output[1:min]
                                             V<-V+t(V)
                                             diag(V)<-diag(V)/2
                                             # deriva gives information matrix
                                             tr <- sum(diag(V))/npm
                                             V0<-V}
                                           ite<-ite+1
                                           pbr_compu<-1
                                           break
                                         }
                                         
                                         
                                         fu<-output[(min+1):length(output)]
                                         V<-matrix(0,nrow=npm,ncol=npm)
                                         V[upper.tri(V,diag=T)]<-output[1:min]
                                         V<-V+t(V)
                                         diag(V)<-diag(V)/2
                                         
                                         # projection V into space of eigen vectors
                                         # tol <- 1e-12
                                         # eig <- eigen(V, symmetric = TRUE)
                                         # vals <- eig$values
                                         # vecs <- eig$vectors
                                         # 
                                         # # tronquer les valeurs propres trop petites
                                         # vals_clipped <- pmax(vals, tol)
                                         # 
                                         # # reconstruire la matrice projetée
                                         # H_proj <- vecs %*% diag(vals_clipped) %*% t(vecs)
                                         # 
                                         # 
                                         # pp<-matrix(0,nrow=npm,ncol=npm)
                                         # diag(pp)<-diag(V)
                                         # deriva gives information matrix
                                         tr <- sum(diag(V))/npm
                                         V0<-V
                                         
                                         
                                         eigen.values<-eigen(V,symmetric=T,only.values=T)$values
                                         
                                         idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                         
                                         
                                         idpos0<-idpos
                                         
                                         ncount<-da<-ga<-0
                                         
                                         while(idpos != 0){
                                           
                                           if(ncount==0){ 
                                             ga <- 0.01
                                             da <- 1E-2
                                           }else{
                                             if(((ncount <= 3) | (ga >= 1)) ){
                                               da <- da * 5
                                             }else{# if ncount > 10 only update ga 
                                               ga <- ga * 5
                                               # do not put ga at 1 as no countmax otherwise infinite while 
                                               if(ga > 1) ga <- 1
                                             }
                                           }
                                           
                                           ncount <- ncount + 1
                                           
                                           diagV <- diag(V)
                                           # put abs (1-ga) better than 1-ga cause ga can now be >1
                                           diagV<-ifelse(diagV!=0,diagV+da*(abs((1.e0-ga))*abs(diagV)+ga*tr),
                                                         da*ga*tr)
                                           
                                           diag(V)<-diagV
                                           # if we have a convex log-vraisemblance in eta then :
                                           # all eigen  values of the hessienne are >0.
                                           
                                           if(sum(V==Inf)>0|sum(V==-Inf)>0){break}
                                           eigen.values<-eigen(V,symmetric=T,only.values=T)$values
                                           # check if hessienne defined positive
                                           idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                           
                                           
                                           
                                         }
                                         
                                         
                                         if(idpos!=0){
                                           
                                           warning("Hessian not defined positive")
                                           pbr_compu<-2
                                           ite<-ite+1
                                           break
                                         }
                                         
                                         
                                         # update beta
                                         
                                         output.cv<-DYNcv.model(beta=beta,
                                                                nva01=npm01,
                                                                nva02=npm02,
                                                                nva12=npm12,
                                                                nva01Y=npm01Y,
                                                                nva02Y=npm02Y,
                                                                nva12Y=npm12Y,
                                                                fix=fix0[(size_spline+1):size_V],
                                                                penalty.factor=penalty.factor,
                                                                penalty=penalty,
                                                                v=V,
                                                                fu=fu,
                                                                lambda=lambda[id.lambda,],
                                                                alpha=alpha
                                         )
                                         
                                         # verify validity of parameters update 
                                         # and that we are better than previous estimates 
                                         
                                         b<-c(s,output.cv$b)
                                         
                                         betanew<-b[(size_spline+1):size_V]
                                         
                                         # penalised loglik see if inferior to previous
                                         res<-gaussDYNidmlLikelihoodpena(b=b,
                                                                         npm=size_V,
                                                                         npar=size_V,
                                                                         bfix=1,
                                                                         fix=rep(0,size_V),
                                                                         zi01=knots01,
                                                                         zi02=knots02,
                                                                         zi12=knots12,
                                                                         ctime=ctime,
                                                                         no=N,
                                                                         nz01=nknots01,
                                                                         nz02=nknots02,
                                                                         nz12=nknots12,
                                                                         ve01=ve01,
                                                                         ve02=ve02,
                                                                         ve12=ve12,
                                                                         dimnva01=dimnva01,
                                                                         dimnva02=dimnva02,
                                                                         dimnva12=dimnva12,
                                                                         nva01=nvat01,
                                                                         nva02=nvat02,
                                                                         nva12=nvat12,
                                                                         t0=t0,
                                                                         t1=t1,
                                                                         t2=t2,
                                                                         t3=t3,
                                                                         troncature=troncature,
                                                                         y01=y01k,
                                                                         y02=y02k,
                                                                         y12=y12k,
                                                                         p01=p01,
                                                                         p02=p02,
                                                                         p12=p12,
                                                                         dimp01=dimp01,
                                                                         dimp02=dimp02,
                                                                         dimp12=dimp12,
                                                                         Ntime=NtimePoints,
                                                                         lambda=lambda[id.lambda,],
                                                                         alpha=alpha,
                                                                         penalty.factor=penalty.factor,
                                                                         penalty=penalty)
                                         
                                         
                                         # we want to maximise the loglik thus : 
                                         # we have issue if res is NA or if not higher than previous one 
                                         # if not better or do not exist need to readjust
                                         # value of beta 
                                         
                                         
                                         if(res %in%c(-1e9,1e9) | res < fn.value){
                                           
                                           print(paste0("needed update at ite :",ite))
                                           th<-1e-5
                                           step<-log(1.5)
                                           delta<-output.cv$b-c(beta)
                                           
                                           maxt <- max(abs(delta)) 
                                           
                                           if(maxt == 0){
                                             vw <- th
                                           }else{
                                             vw <- th/maxt
                                           }
                                           if(ite>0){
                                             res.out.error <- list("old.b"=round(c(s,beta)),
                                                                   "old.rl"=round(fn.value),
                                                                   "old.ca"=round(eval.cv.beta[ite]),
                                                                   "old.cb"=round(eval.cv.loglik[ite]))
                                           }else{
                                             res.out.error <- list("old.b"=round(c(s,beta)),
                                                                   "old.rl"=round(fn.value),
                                                                   "old.ca"=round(1),
                                                                   "old.cb"=round(1))
                                           }
                                           m<-1
                                           betanew<-beta
                                           # from mla package
                                           
                                           sears<-searpas(vw=vw,
                                                          step=step,
                                                          b=beta,
                                                          delta=delta,
                                                          funcpa=gaussDYNidmlLikelihoodpena,
                                                          res.out.error=res.out.error,
                                                          npm=npm,
                                                          npar=size_V,
                                                          bfix=s,
                                                          fix=fix0,
                                                          zi01=knots01,
                                                          zi02=knots02,
                                                          zi12=knots12,
                                                          ctime=ctime,
                                                          no=N,
                                                          nz01=nknots01,
                                                          nz02=nknots02,
                                                          nz12=nknots12,
                                                          ve01=ve01,
                                                          ve02=ve02,
                                                          ve12=ve12,
                                                          dimnva01=dimnva01,
                                                          dimnva02=dimnva02,
                                                          dimnva12=dimnva12,
                                                          nva01=nvat01,
                                                          nva02=nvat02,
                                                          nva12=nvat12,
                                                          t0=t0,
                                                          t1=t1,
                                                          t2=t2,
                                                          t3=t3,
                                                          troncature=troncature,
                                                          y01=y01k,
                                                          y02=y02k,
                                                          y12=y12k,
                                                          p01=p01,
                                                          p02=p02,
                                                          p12=p12,
                                                          dimp01=dimp01,
                                                          dimp02=dimp02,
                                                          dimp12=dimp12,
                                                          Ntime=NtimePoints,
                                                          lambda=lambda[id.lambda,],
                                                          alpha=alpha,
                                                          penalty.factor=penalty.factor,
                                                          penalty=penalty)
                                           
                                           
                                           betanew<-beta+delta*sears$vw
                                           b<-c(s,betanew)
                                           
                                           
                                           res<-gaussDYNidmlLikelihoodpena(b=b,
                                                                           npm=size_V,
                                                                           npar=size_V,
                                                                           bfix=1,
                                                                           fix=rep(0,size_V),
                                                                           zi01=knots01,
                                                                           zi02=knots02,
                                                                           zi12=knots12,
                                                                           ctime=ctime,
                                                                           no=N,
                                                                           nz01=nknots01,
                                                                           nz02=nknots02,
                                                                           nz12=nknots12,
                                                                           ve01=ve01,
                                                                           ve02=ve02,
                                                                           ve12=ve12,
                                                                           dimnva01=dimnva01,
                                                                           dimnva02=dimnva02,
                                                                           dimnva12=dimnva12,
                                                                           nva01=nvat01,
                                                                           nva02=nvat02,
                                                                           nva12=nvat12,
                                                                           t0=t0,
                                                                           t1=t1,
                                                                           t2=t2,
                                                                           t3=t3,
                                                                           troncature=troncature,
                                                                           y01=y01k,
                                                                           y02=y02k,
                                                                           y12=y12k,
                                                                           p01=p01,
                                                                           p02=p02,
                                                                           p12=p12,
                                                                           dimp01=dimp01,
                                                                           dimp02=dimp02,
                                                                           dimp12=dimp12,
                                                                           Ntime=NtimePoints,
                                                                           lambda=lambda[id.lambda,],
                                                                           alpha=alpha,
                                                                           penalty.factor=penalty.factor,
                                                                           penalty=penalty)
                                           
                                         }
                                         # if not better or do not exist need to readjust
                                         # value of beta 
                                         if(res %in%c(-1e9,1e9) | any(is.infinite(c(s,betanew)))){
                                           
                                           ite<-ite+1
                                           validity<-F
                                           eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                           eval.validity[ite]<-validity
                                           # save output over iterations 
                                           pbr_compu<-3
                                           break
                                         }else{validity<-T}
                                         
                                         
                                         # betanew already include s
                                         b<-c(s,betanew)
                                         
                                         bfix<-b[fix0.beta==1]
                                         b<-b[fix0.beta==0]
                                         # update modelPar
                                         
                                         
                                         output.mla<- marqLevAlg::mla(b=b,
                                                                      fn=gaussDYNidmlLikelihood,
                                                                      epsa=epsa,
                                                                      epsb=epsb,
                                                                      epsd=epsd,
                                                                      maxiter=maxiter.pena,
                                                                      minimize=F,
                                                                      npm=npmspline,
                                                                      npar=size_V,
                                                                      bfix=bfix,
                                                                      fix=fix0.beta,
                                                                      zi01=knots01,
                                                                      zi02=knots02,
                                                                      zi12=knots12,
                                                                      ctime=ctime,
                                                                      no=N,
                                                                      nz01=nknots01,
                                                                      nz02=nknots02,
                                                                      nz12=nknots12,
                                                                      ve01=ve01,
                                                                      ve02=ve02,
                                                                      ve12=ve12,
                                                                      dimnva01=dimnva01,
                                                                      dimnva02=dimnva02,
                                                                      dimnva12=dimnva12,
                                                                      nva01=nvat01,
                                                                      nva02=nvat02,
                                                                      nva12=nvat12,
                                                                      t0=t0,
                                                                      t1=t1,
                                                                      t2=t2,
                                                                      t3=t3,
                                                                      troncature=troncature,
                                                                      y01=y01k,
                                                                      y02=y02k,
                                                                      y12=y12k,
                                                                      p01=p01,
                                                                      p02=p02,
                                                                      p12=p12,
                                                                      dimp01=dimp01,
                                                                      dimp02=dimp02,
                                                                      dimp12=dimp12,
                                                                      Ntime=NtimePoints)
                                         
                                         # look at convergence for each lambda :
                                         
                                         # new values for splines:
                                         snew<-s
                                         snew[fix00[1:size_spline]==0]<-output.mla$b
                                         if(nvat01>0){
                                           b01<-betanew[1:nvat01][penalty.factor[1:nvat01]==1]
                                           if(p01>0){
                                             b01<-c(b01,betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1])
                                           }
                                         }else{
                                           if(p01>0){
                                             b01<-betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1]
                                           }else{
                                             b01<-0
                                           }
                                         }
                                         
                                         if(nvat02>0){
                                           b02<-betanew[(nvat01+1):(nvat01+nvat02)][penalty.factor[(nvat01+1):(nvat01+nvat02)]==1]
                                           if(p02>0){
                                             b02<-c(b02,betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1])
                                           }
                                         }else{
                                           if(p02>0){
                                             b02<-betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1]
                                           }else{b02<-0}
                                         }
                                         
                                         if(nvat12>0){
                                           b12<-betanew[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)][penalty.factor[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)]==1]
                                           if(p12>0){
                                             b12<-c(b12,betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1])
                                           }
                                         }else{
                                           if(p12>0){
                                             b12<-betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1]
                                           }else{b12<-0}
                                         }
                                         # maximisation issue : lpen =l - pen
                                         if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
                                           fn.valuenew<-output.mla$fn.value-lambda[id.lambda,1]*alpha*sum(abs(b01))-lambda[id.lambda,1]*(1-alpha)*sum(b01*b01)
                                           fn.valuenew<-fn.valuenew-lambda[id.lambda,2]*alpha*sum(abs(b02))-lambda[id.lambda,2]*(1-alpha)*sum(b02*b02)
                                           fn.valuenew<-fn.valuenew-lambda[id.lambda,3]*alpha*sum(abs(b12))-lambda[id.lambda,3]*(1-alpha)*sum(b12*b12)
                                           
                                         }
                                         
                                         if(penalty=="mcp"){
                                           
                                           p01<-rep(alpha*lambda[id.lambda,1]*lambda[id.lambda,1]/2,length(b01))
                                           idbeta<-which(b01<=alpha*lambda[id.lambda,1])
                                           p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])-((b01[idbeta]*b01[idbeta])/2*alpha)
                                           
                                           p02<-rep(alpha*lambda[id.lambda,2]*lambda[id.lambda,2]/2,length(b02))
                                           idbeta<-which(b02<=alpha*lambda[id.lambda,2])
                                           p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])-((b02[idbeta]*b02[idbeta])/2*alpha)
                                           
                                           p12<-rep(alpha*lambda[id.lambda,3]*lambda[id.lambda,3]/2,length(b12))
                                           idbeta<-which(b12<=alpha*lambda[id.lambda,3])
                                           p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])-((b12[idbeta]*b12[idbeta])/2*alpha)
                                           
                                           fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                           
                                         }
                                         
                                         if(penalty=="scad"){
                                           
                                           p01<-rep((lambda[id.lambda,1]^2)*(alpha+1)/2,length(b01))
                                           idbeta<-which(b01<=lambda[id.lambda,1])
                                           p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])
                                           idbeta<-which(abs(b01)<lambda[id.lambda,1]*alpha)
                                           p01[idbeta]<-(2*alpha*lambda[id.lambda,1]*abs(b01[idbeta])-b01[idbeta]^2-lambda[id.lambda,1]^2)/(2*(alpha-1))
                                           
                                           p02<-rep((lambda[id.lambda,2]^2)*(alpha+1)/2,length(b02))
                                           idbeta<-which(b02<=lambda[id.lambda,2])
                                           p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])
                                           idbeta<-which(abs(b02)<lambda[id.lambda,2]*alpha)
                                           p02[idbeta]<-(2*alpha*lambda[id.lambda,2]*abs(b02[idbeta])-b02[idbeta]^2-lambda[id.lambda,2]^2)/(2*(alpha-1))
                                           
                                           p12<-rep((lambda[id.lambda,3]^2)*(alpha+1)/2,length(b12))
                                           idbeta<-which(b12<=lambda[id.lambda,3])
                                           p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])
                                           idbeta<-which(abs(b12)<lambda[id.lambda,3]*alpha)
                                           p12[idbeta]<-(2*alpha*lambda[id.lambda,3]*abs(b12[idbeta])-b12[idbeta]^2-lambda[id.lambda,3]^2)/(2*(alpha-1))
                                           
                                           fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                           
                                         }
                                         
                                         ite<-ite+1
                                         
                                         # check cv criterias 
                                         
                                         eval.cv.spline[ite]<-sum((snew-s)^2)
                                         eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                         eval.cv.loglik[ite]<-abs((fn.valuenew-fn.value)/fn.value)
                                         eval.loglik[ite]<-fn.valuenew
                                         eval.validity[ite]<-validity
                                         
                                         s<-snew
                                         beta<-betanew
                                         fn.value<-fn.valuenew
                                         
                                         # eval.cv beta valid only if validity.param=T
                                         # add ite >1 so that do not stop at first ite ?
                                         if(((eval.cv.beta[ite]+eval.cv.spline[ite])<epsa) & eval.cv.loglik[ite]<epsb & validity==T){
                                           converged<-T}
                                         
                                         
                                       }
                                       
                                       if(maxiter<=ite & converged==F){
                                         istop<-2
                                       }else{
                                         if(ite<=maxiter & converged==T){
                                           istop<-1
                                           # need to recalculate second derivatives 
                                           
                                           b<-c(s,beta)
                                           bfix<-b[fix0==1]
                                           b<-b[fix0==0]
                                           
                                           output<-DYNderiva(
                                             h=1e-8,
                                             b=b,
                                             npm=npm,
                                             npar=size_V,
                                             bfix=bfix,
                                             fix=fix0,
                                             zi01=knots01,
                                             zi02=knots02,
                                             zi12=knots12,
                                             ctime=ctime,
                                             no=N,
                                             nz01=nknots01,
                                             nz02=nknots02,
                                             nz12=nknots12,
                                             ve01=ve01,
                                             ve02=ve02,
                                             ve12=ve12,
                                             dimnva01=dimnva01,
                                             dimnva02=dimnva02,
                                             dimnva12=dimnva12,
                                             nva01=nvat01,
                                             nva02=nvat02,
                                             nva12=nvat12,
                                             t0=t0,
                                             t1=t1,
                                             t2=t2,
                                             t3=t3,
                                             troncature=troncature,
                                             y01=y01k,
                                             y02=y02k,
                                             y12=y12k,
                                             p01=p01,
                                             p02=p02,
                                             p12=p12,
                                             dimp01=dimp01,
                                             dimp02=dimp02,
                                             dimp12=dimp12,
                                             Ntime=NtimePoints)
                                           
                                           output<-output$v
                                           fu<-output[(min+1):length(output)]
                                           V<-matrix(0,nrow=npm,ncol=npm)
                                           V[upper.tri(V,diag=T)]<-output[1:min]
                                           V<-V+t(V)
                                           diag(V)<-diag(V)/2
                                           # deriva gives information matrix
                                           V0<-V
                                         }else{
                                           if(pbr_compu==1){istop<-3}
                                           if(pbr_compu==2){istop<-4}
                                           if(pbr_compu==3){istop<-5}
                                         }
                                       }
                                       
                                       # if stop==1 we can give matrix of second derivatives 
                                       
                                       combine<-combine+1
                                       return(list(b=c(s,beta),
                                                   H=V0,
                                                   lambda=as.double(lambda[id.lambda,]),
                                                   alpha=alpha,
                                                   fn.value=ifelse(!exists("output.mla"),NA,output.mla$fn.value),
                                                   fn.value.pena=fn.value,
                                                   ni=ite,
                                                   ca.beta=eval.cv.beta,
                                                   ca.spline=eval.cv.spline,
                                                   ca.validity=eval.validity,
                                                   cb=eval.loglik,
                                                   istop=istop,
                                                   combine=combine))
                                     }
          }else{outputNsample<-foreach::foreach(id.lambda=1:nlambda,
                                         .combine = combine_lambda_mla,
                                         .errorhandling = "remove")%do%{
                                           
                                           
                                           # computation pbr 
                                           
                                           pbr_compu<-0
                                           
                                           if(idsample==1){
                                             beta<-beta.start
                                             s<-s.start
                                           }else{
                                             s<-output[[1]]$b[1:size_spline,id.lambda]
                                             beta<-output[[1]]$b[(size_spline+1):size_V,id.lambda]
                                           }
                                           
                                           converged<-F
                                           ite<-0
                                           # if beta not change do not need to recalculate weights 
                                           H<-T
                                           
                                           eval.cv.spline<-rep(NA,maxiter+1)
                                           eval.cv.beta<-rep(NA,maxiter+1)
                                           eval.cv.loglik<-rep(NA,maxiter+1)
                                           eval.loglik<-rep(NA,maxiter+1)
                                           eval.validity<-rep(NA,maxiter+1)
                                           
                                           
                                           while(converged==F & ite<=maxiter){
                                             
                                             b<-c(s,beta)
                                             bfix<-b[fix0==1]
                                             b<-b[fix0==0]
                                             
                                             
                                             
                                             output<-DYNderivadiag(
                                               h=1e-8,
                                               b=b,
                                               npm=npm,
                                               npar=size_V,
                                               bfix=bfix,
                                               fix=fix0,
                                               zi01=knots01,
                                               zi02=knots02,
                                               zi12=knots12,
                                               ctime=ctime,
                                               no=N,
                                               nz01=nknots01,
                                               nz02=nknots02,
                                               nz12=nknots12,
                                               ve01=ve01,
                                               ve02=ve02,
                                               ve12=ve12,
                                               dimnva01=dimnva01,
                                               dimnva02=dimnva02,
                                               dimnva12=dimnva12,
                                               nva01=nvat01,
                                               nva02=nvat02,
                                               nva12=nvat12,
                                               t0=t0,
                                               t1=t1,
                                               t2=t2,
                                               t3=t3,
                                               troncature=troncature,
                                               y01=y01k,
                                               y02=y02k,
                                               y12=y12k,
                                               p01=p01,
                                               p02=p02,
                                               p12=p12,
                                               dimp01=dimp01,
                                               dimp02=dimp02,
                                               dimp12=dimp12,
                                               Ntime=NtimePoints)
                                             
                                             
                                             output<-output$v
                                             if(ite==0){
                                               fn.value<-gaussDYNidmlLikelihoodpena(b=b,
                                                                                    npm=npm,
                                                                                    npar=size_V,
                                                                                    bfix=bfix,
                                                                                    fix=fix0,
                                                                                    zi01=knots01,
                                                                                    zi02=knots02,
                                                                                    zi12=knots12,
                                                                                    ctime=ctime,
                                                                                    no=N,
                                                                                    nz01=nknots01,
                                                                                    nz02=nknots02,
                                                                                    nz12=nknots12,
                                                                                    ve01=ve01,
                                                                                    ve02=ve02,
                                                                                    ve12=ve12,
                                                                                    dimnva01=dimnva01,
                                                                                    dimnva02=dimnva02,
                                                                                    dimnva12=dimnva12,
                                                                                    nva01=nvat01,
                                                                                    nva02=nvat02,
                                                                                    nva12=nvat12,
                                                                                    t0=t0,
                                                                                    t1=t1,
                                                                                    t2=t2,
                                                                                    t3=t3,
                                                                                    troncature=troncature,
                                                                                    y01=y01k,
                                                                                    y02=y02k,
                                                                                    y12=y12k,
                                                                                    p01=p01,
                                                                                    p02=p02,
                                                                                    p12=p12,
                                                                                    dimp01=dimp01,
                                                                                    dimp02=dimp02,
                                                                                    dimp12=dimp12,
                                                                                    Ntime=NtimePoints,
                                                                                    
                                                                                    lambda=lambda[id.lambda,],
                                                                                    alpha=alpha,
                                                                                    penalty.factor=penalty.factor,
                                                                                    penalty=penalty)
                                             }
                                             
                                             if(any(is.na(output$V))|any(output$V==Inf) |any(output$V==-Inf)|any(is.na(output$fu))|any(output$fu==Inf) |any(output$fu==-Inf)){
                                               warning("Computational error for calculation of the hessian : division by 0 or Infinite value")
                                               if(ite==0){
                                                 
                                                 fu<-output[(min+1):(2*min)]
                                                 V<-matrix(0,nrow=npm,ncol=npm)
                                                 diag(V)<-output[1:min]
                                                 # deriva gives information matrix
                                                 tr <- sum(diag(V))/npm
                                                 V0<-V}
                                               ite<-ite+1
                                               pbr_compu<-1
                                               break
                                             }
                                             
                                             
                                             fu<-output[(min+1):(2*min)]
                                             V<-matrix(0,nrow=npm,ncol=npm)
                                             diag(V)<-output[1:min]
                                             # deriva gives information matrix
                                             tr <- sum(diag(V))/npm
                                             V0<-V
                                             
                                             eigen.values<-diag(V)
                                             
                                             idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                             
                                             
                                             idpos0<-idpos
                                             
                                             ncount<-da<-ga<-0
                                             
                                             
                                             while(idpos != 0){
                                               
                                               if(ncount==0){ 
                                                 ga <- 0.01
                                                 da <- 1E-2
                                               }else{
                                                 if(((ncount <= 3) | (ga >= 1)) ){
                                                   da <- da * 5
                                                 }else{# if ncount > 10 only update ga 
                                                   ga <- ga * 5
                                                   # do not put ga at 1 as no countmax otherwise infinite while 
                                                   if(ga > 1) ga <- 1
                                                 }
                                               }
                                               
                                               ncount <- ncount + 1
                                               
                                               diagV <- diag(V)
                                               # put abs (1-ga) better than 1-ga cause ga can now be >1
                                               diagV<-ifelse(diagV!=0,diagV+da*(abs((1.e0-ga))*abs(diagV)+ga*tr),
                                                             da*ga*tr)
                                               
                                               diag(V)<-diagV
                                               # if we have a convex log-vraisemblance in eta then :
                                               # all eigen  values of the hessienne are >0.
                                               
                                               if(sum(V==Inf)>0|sum(V==-Inf)>0){break}
                                               eigen.values<-diag(V)
                                               # check if hessienne defined positive
                                               idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                               
                                               
                                               
                                             }
                                             
                                             if(idpos!=0){
                                               
                                               warning("Hessian not defined positive")
                                               pbr_compu<-2
                                               ite<-ite+1
                                               break
                                             }
                                             
                                             
                                             # update beta
                                             output.cv<-DYNcv.model(beta=beta,
                                                                    nva01=npm01,
                                                                    nva02=npm02,
                                                                    nva12=npm12,
                                                                    nva01Y=npm01Y,
                                                                    nva02Y=npm02Y,
                                                                    nva12Y=npm12Y,
                                                                    fix=fix0[(size_spline+1):size_V],
                                                                    penalty.factor=penalty.factor,
                                                                    penalty=penalty,
                                                                    v=V,
                                                                    fu=fu,
                                                                    lambda=lambda[id.lambda,],
                                                                    alpha=alpha
                                             )
                                             
                                             # verify validity of parameters update 
                                             # and that we are better than previous estimates 
                                             
                                             b<-c(s,output.cv$b)
                                             
                                             betanew<-b[(size_spline+1):size_V]
                                             
                                             # penalised loglik see if inferior to previous
                                             res<-gaussDYNidmlLikelihoodpena(b=b,
                                                                             npm=size_V,
                                                                             npar=size_V,
                                                                             bfix=1,
                                                                             fix=rep(0,size_V),
                                                                             zi01=knots01,
                                                                             zi02=knots02,
                                                                             zi12=knots12,
                                                                             ctime=ctime,
                                                                             no=N,
                                                                             nz01=nknots01,
                                                                             nz02=nknots02,
                                                                             nz12=nknots12,
                                                                             ve01=ve01,
                                                                             ve02=ve02,
                                                                             ve12=ve12,
                                                                             dimnva01=dimnva01,
                                                                             dimnva02=dimnva02,
                                                                             dimnva12=dimnva12,
                                                                             nva01=nvat01,
                                                                             nva02=nvat02,
                                                                             nva12=nvat12,
                                                                             t0=t0,
                                                                             t1=t1,
                                                                             t2=t2,
                                                                             t3=t3,
                                                                             troncature=troncature,
                                                                             y01=y01k,
                                                                             y02=y02k,
                                                                             y12=y12k,
                                                                             p01=p01,
                                                                             p02=p02,
                                                                             p12=p12,
                                                                             dimp01=dimp01,
                                                                             dimp02=dimp02,
                                                                             dimp12=dimp12,
                                                                             Ntime=NtimePoints,
                                                                             
                                                                             lambda=lambda[id.lambda,],
                                                                             alpha=alpha,
                                                                             penalty.factor=penalty.factor,
                                                                             penalty=penalty)
                                             
                                             # we want to maximise the loglik thus : 
                                             # we have issue if res is NA or if not higher than previous one 
                                             # if not better or do not exist need to readjust
                                             # value of beta 
                                             
                                             if(res %in%c(-1e9,1e9) | res < fn.value){
                                               
                                               th<-1e-5
                                               step<-log(1.5)
                                               delta<-output.cv$b-c(beta)
                                               
                                               maxt <- max(abs(delta)) 
                                               
                                               if(maxt == 0){
                                                 vw <- th
                                               }else{
                                                 vw <- th/maxt
                                               }
                                               if(ite>0){
                                                 res.out.error <- list("old.b"=round(c(s,beta)),
                                                                       "old.rl"=round(fn.value),
                                                                       "old.ca"=round(eval.cv.beta[ite]),
                                                                       "old.cb"=round(eval.cv.loglik[ite]))
                                               }else{
                                                 res.out.error <- list("old.b"=round(c(s,beta)),
                                                                       "old.rl"=round(fn.value),
                                                                       "old.ca"=round(1),
                                                                       "old.cb"=round(1))
                                               }
                                               m<-1
                                               betanew<-beta
                                               # from mla package
                                               
                                               sears<-searpas(vw=vw,
                                                              step=step,
                                                              b=beta,
                                                              delta=delta,
                                                              funcpa=gaussDYNidmlLikelihoodpena,
                                                              res.out.error=res.out.error,
                                                              npm=npm,
                                                              npar=size_V,
                                                              bfix=s,
                                                              fix=fix0,
                                                              zi01=knots01,
                                                              zi02=knots02,
                                                              zi12=knots12,
                                                              ctime=ctime,
                                                              no=N,
                                                              nz01=nknots01,
                                                              nz02=nknots02,
                                                              nz12=nknots12,
                                                              ve01=ve01,
                                                              ve02=ve02,
                                                              ve12=ve12,
                                                              dimnva01=dimnva01,
                                                              dimnva02=dimnva02,
                                                              dimnva12=dimnva12,
                                                              nva01=nvat01,
                                                              nva02=nvat02,
                                                              nva12=nvat12,
                                                              t0=t0,
                                                              t1=t1,
                                                              t2=t2,
                                                              t3=t3,
                                                              troncature=troncature,
                                                              y01=y01k,
                                                              y02=y02k,
                                                              y12=y12k,
                                                              p01=p01,
                                                              p02=p02,
                                                              p12=p12,
                                                              dimp01=dimp01,
                                                              dimp02=dimp02,
                                                              dimp12=dimp12,
                                                              Ntime=NtimePoints,
                                                              lambda=lambda[id.lambda,],
                                                              alpha=alpha,
                                                              penalty.factor=penalty.factor,
                                                              penalty=penalty)
                                               
                                               
                                               betanew<-beta+delta*sears$vw
                                               b<-c(s,betanew)
                                               
                                               
                                               res<-gaussDYNidmlLikelihoodpena(b=b,
                                                                               npm=size_V,
                                                                               npar=size_V,
                                                                               bfix=1,
                                                                               fix=rep(0,size_V),
                                                                               zi01=knots01,
                                                                               zi02=knots02,
                                                                               zi12=knots12,
                                                                               ctime=ctime,
                                                                               no=N,
                                                                               nz01=nknots01,
                                                                               nz02=nknots02,
                                                                               nz12=nknots12,
                                                                               ve01=ve01,
                                                                               ve02=ve02,
                                                                               ve12=ve12,
                                                                               dimnva01=dimnva01,
                                                                               dimnva02=dimnva02,
                                                                               dimnva12=dimnva12,
                                                                               nva01=nvat01,
                                                                               nva02=nvat02,
                                                                               nva12=nvat12,
                                                                               t0=t0,
                                                                               t1=t1,
                                                                               t2=t2,
                                                                               t3=t3,
                                                                               troncature=troncature,
                                                                               y01=y01k,
                                                                               y02=y02k,
                                                                               y12=y12k,
                                                                               p01=p01,
                                                                               p02=p02,
                                                                               p12=p12,
                                                                               dimp01=dimp01,
                                                                               dimp02=dimp02,
                                                                               dimp12=dimp12,
                                                                               Ntime=NtimePoints,
                                                                               
                                                                               lambda=lambda[id.lambda,],
                                                                               alpha=alpha,
                                                                               penalty.factor=penalty.factor,
                                                                               penalty=penalty)
                                               
                                             }
                                             # if not better or do not exist need to readjust
                                             # value of beta 
                                             if(res %in%c(-1e9,1e9) | any(is.infinite(c(s,betanew)))){
                                               
                                               ite<-ite+1
                                               validity<-F
                                               eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                               eval.validity[ite]<-validity
                                               # save output over iterations 
                                               pbr_compu<-3
                                               break
                                             }else{validity<-T}
                                             
                                             
                                             # betanew already include s
                                             b<-c(s,betanew)
                                             
                                             bfix<-b[fix0.beta==1]
                                             b<-b[fix0.beta==0]
                                             # update modelPar
                                             
                                             
                                             output.mla<- marqLevAlg::mla(b=b,
                                                                          fn=gaussDYNidmlLikelihood,
                                                                          epsa=epsa,
                                                                          epsb=epsb,
                                                                          epsd=epsd,
                                                                          maxiter=maxiter.pena,
                                                                          minimize=F,
                                                                          npm=npmspline,
                                                                          npar=size_V,
                                                                          bfix=bfix,
                                                                          fix=fix0.beta,
                                                                          zi01=knots01,
                                                                          zi02=knots02,
                                                                          zi12=knots12,
                                                                          ctime=ctime,
                                                                          no=N,
                                                                          nz01=nknots01,
                                                                          nz02=nknots02,
                                                                          nz12=nknots12,
                                                                          ve01=ve01,
                                                                          ve02=ve02,
                                                                          ve12=ve12,
                                                                          dimnva01=dimnva01,
                                                                          dimnva02=dimnva02,
                                                                          dimnva12=dimnva12,
                                                                          nva01=nvat01,
                                                                          nva02=nvat02,
                                                                          nva12=nvat12,
                                                                          t0=t0,
                                                                          t1=t1,
                                                                          t2=t2,
                                                                          t3=t3,
                                                                          troncature=troncature,
                                                                          y01=y01k,
                                                                          y02=y02k,
                                                                          y12=y12k,
                                                                          p01=p01,
                                                                          p02=p02,
                                                                          p12=p12,
                                                                          dimp01=dimp01,
                                                                          dimp02=dimp02,
                                                                          dimp12=dimp12,
                                                                          Ntime=NtimePoints)
                                             
                                             # look at convergence for each lambda :
                                             
                                             # new values for splines:
                                             snew<-s
                                             snew[fix00[1:size_spline]==0]<-output.mla$b
                                             if(nvat01>0){
                                               b01<-betanew[1:nvat01][penalty.factor[1:nvat01]==1]
                                               if(p01>0){
                                                 b01<-c(b01,betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1])
                                               }
                                             }else{
                                               if(p01>0){
                                                 b01<-betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1]
                                               }else{
                                                 b01<-0
                                               }
                                             }
                                             
                                             if(nvat02>0){
                                               b02<-betanew[(nvat01+1):(nvat01+nvat02)][penalty.factor[(nvat01+1):(nvat01+nvat02)]==1]
                                               if(p02>0){
                                                 b02<-c(b02,betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1])
                                               }
                                             }else{
                                               if(p02>0){
                                                 b02<-betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1]
                                               }else{b02<-0}
                                             }
                                             
                                             if(nvat12>0){
                                               b12<-betanew[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)][penalty.factor[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)]==1]
                                               if(p12>0){
                                                 b12<-c(b12,betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1])
                                               }
                                             }else{
                                               if(p12>0){
                                                 b12<-betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1]
                                               }else{b12<-0}
                                             }
                                             # maximisation issue : lpen =l - pen
                                             if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
                                               fn.valuenew<-output.mla$fn.value-lambda[id.lambda,1]*alpha*sum(abs(b01))-lambda[id.lambda,1]*(1-alpha)*sum(b01*b01)
                                               fn.valuenew<-fn.valuenew-lambda[id.lambda,2]*alpha*sum(abs(b02))-lambda[id.lambda,2]*(1-alpha)*sum(b02*b02)
                                               fn.valuenew<-fn.valuenew-lambda[id.lambda,3]*alpha*sum(abs(b12))-lambda[id.lambda,3]*(1-alpha)*sum(b12*b12)
                                               
                                             }
                                             
                                             if(penalty=="mcp"){
                                               
                                               p01<-rep(alpha*lambda[id.lambda,1]*lambda[id.lambda,1]/2,length(b01))
                                               idbeta<-which(b01<=alpha*lambda[id.lambda,1])
                                               p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])-((b01[idbeta]*b01[idbeta])/2*alpha)
                                               
                                               p02<-rep(alpha*lambda[id.lambda,2]*lambda[id.lambda,2]/2,length(b02))
                                               idbeta<-which(b02<=alpha*lambda[id.lambda,2])
                                               p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])-((b02[idbeta]*b02[idbeta])/2*alpha)
                                               
                                               p12<-rep(alpha*lambda[id.lambda,3]*lambda[id.lambda,3]/2,length(b12))
                                               idbeta<-which(b12<=alpha*lambda[id.lambda,3])
                                               p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])-((b12[idbeta]*b12[idbeta])/2*alpha)
                                               
                                               fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                               
                                             }
                                             
                                             if(penalty=="scad"){
                                               
                                               p01<-rep((lambda[id.lambda,1]^2)*(alpha+1)/2,length(b01))
                                               idbeta<-which(b01<=lambda[id.lambda,1])
                                               p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])
                                               idbeta<-which(abs(b01)<lambda[id.lambda,1]*alpha)
                                               p01[idbeta]<-(2*alpha*lambda[id.lambda,1]*abs(b01[idbeta])-b01[idbeta]^2-lambda[id.lambda,1]^2)/(2*(alpha-1))
                                               
                                               p02<-rep((lambda[id.lambda,2]^2)*(alpha+1)/2,length(b02))
                                               idbeta<-which(b02<=lambda[id.lambda,2])
                                               p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])
                                               idbeta<-which(abs(b02)<lambda[id.lambda,2]*alpha)
                                               p02[idbeta]<-(2*alpha*lambda[id.lambda,2]*abs(b02[idbeta])-b02[idbeta]^2-lambda[id.lambda,2]^2)/(2*(alpha-1))
                                               
                                               p12<-rep((lambda[id.lambda,3]^2)*(alpha+1)/2,length(b12))
                                               idbeta<-which(b12<=lambda[id.lambda,3])
                                               p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])
                                               idbeta<-which(abs(b12)<lambda[id.lambda,3]*alpha)
                                               p12[idbeta]<-(2*alpha*lambda[id.lambda,3]*abs(b12[idbeta])-b12[idbeta]^2-lambda[id.lambda,3]^2)/(2*(alpha-1))
                                               
                                               fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                               
                                             }
                                             
                                             ite<-ite+1
                                             
                                             # check cv criterias 
                                             
                                             eval.cv.spline[ite]<-sum((snew-s)^2)
                                             eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                             eval.cv.loglik[ite]<-abs((fn.valuenew-fn.value)/fn.value)
                                             eval.loglik[ite]<-fn.valuenew
                                             eval.validity[ite]<-validity
                                             
                                             s<-snew
                                             beta<-betanew
                                             fn.value<-fn.valuenew
                                             
                                             # eval.cv beta valid only if validity.param=T
                                             if(((eval.cv.beta[ite]+ eval.cv.spline[ite])<epsa) & eval.cv.loglik[ite]<epsb & validity==T){
                                               converged<-T}
                                             
                                             
                                           }
                                           if(maxiter<=ite & converged==F){
                                             istop<-2
                                           }else{
                                             if(ite<=maxiter & converged==T){
                                               istop<-1
                                               # need to recalculate second derivatives 
                                               
                                               b<-c(s,beta)
                                               bfix<-b[fix0==1]
                                               b<-b[fix0==0]
                                               
                                               output<-DYNderivadiag(
                                                 h=1e-8,
                                                 b=b,
                                                 npm=npm,
                                                 npar=size_V,
                                                 bfix=bfix,
                                                 fix=fix0,
                                                 zi01=knots01,
                                                 zi02=knots02,
                                                 zi12=knots12,
                                                 ctime=ctime,
                                                 no=N,
                                                 nz01=nknots01,
                                                 nz02=nknots02,
                                                 nz12=nknots12,
                                                 ve01=ve01,
                                                 ve02=ve02,
                                                 ve12=ve12,
                                                 dimnva01=dimnva01,
                                                 dimnva02=dimnva02,
                                                 dimnva12=dimnva12,
                                                 nva01=nvat01,
                                                 nva02=nvat02,
                                                 nva12=nvat12,
                                                 t0=t0,
                                                 t1=t1,
                                                 t2=t2,
                                                 t3=t3,
                                                 troncature=troncature,
                                                 y01=y01k,
                                                 y02=y02k,
                                                 y12=y12k,
                                                 p01=p01,
                                                 p02=p02,
                                                 p12=p12,
                                                 dimp01=dimp01,
                                                 dimp02=dimp02,
                                                 dimp12=dimp12,
                                                 Ntime=NtimePoints)
                                               
                                               output<-output$v
                                               
                                               
                                               fu<-output[(min+1):(2*min)]
                                               V<-matrix(0,nrow=npm,ncol=npm)
                                               diag(V)<-output[1:min]
                                               # deriva gives information matrix
                                               V0<-V
                                             }else{
                                               if(pbr_compu==1){istop<-3}
                                               if(pbr_compu==2){istop<-4}
                                               if(pbr_compu==3){istop<-5}
                                             }
                                           }
                                           
                                           # if stop==1 we can give matrix of second derivatives 
                                           
                                           combine<-combine+1
                                           
                                           return(list(b=c(s,beta),
                                                       H=V0,
                                                       lambda=as.double(lambda[id.lambda,]),
                                                       alpha=alpha,
                                                       fn.value=ifelse(!exists("output.mla"),NA,output.mla$fn.value),
                                                       fn.value.pena=fn.value,
                                                       ni=ite,
                                                       ca.beta=eval.cv.beta,
                                                       ca.spline=eval.cv.spline,
                                                       ca.validity=eval.validity,
                                                       cb=eval.loglik,
                                                       istop=istop,
                                                       combine=combine))
                                         }}
        
        
        outputall[[idsample]]<-outputNsample
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
                                                       assoc=assoc,
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
                                   
                                   
                                 }else{
                                   y01k<-rep(0,N*NtimePoints)
                                 }
                                 
                                 if(length(outcome02)>=1){
                                   y02k<-dataY[dataY$Outcome%in%outcome02,]
                                   # order  by individual and timeline 
                                   y02k<-y02k[order(y02k[,colnames(y02k)%in%id],y02k$order),]
                                   
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
        
        
        
          
          
          id.lambda<-NULL # for cran check 
          if(partialH==F){
            
            outputNsample<-foreach::foreach(id.lambda=1:nlambda,
                                     .combine = combine_lambda_mla,
                                     .errorhandling = "remove")%do%{
                                       
                                       
                                       # computation pbr 
                                       
                                       pbr_compu<-0
                                       
                                      
                                         beta<-beta.start
                                         s<-s.start
                                       
                                       
                                       converged<-F
                                       ite<-0
                                       # if beta not change do not need to recalculate weights 
                                       H<-T
                                       
                                       eval.cv.spline<-rep(NA,maxiter+1)
                                       eval.cv.beta<-rep(NA,maxiter+1)
                                       eval.cv.loglik<-rep(NA,maxiter+1)
                                       eval.loglik<-rep(NA,maxiter+1)
                                       eval.validity<-rep(NA,maxiter+1)
                                       
                                       
                                       while(converged==F & ite<=maxiter){
                                         
                                         
                                         b<-c(s,beta)
                                         bfix<-b[fix0==1]
                                         b<-b[fix0==0]
                                         
                                         
                                         output<-DYNderiva(
                                           h=1e-8,
                                           b=b,
                                           npm=npm,
                                           npar=size_V,
                                           bfix=bfix,
                                           fix=fix0,
                                           zi01=knots01,
                                           zi02=knots02,
                                           zi12=knots12,
                                           ctime=ctime,
                                           no=N,
                                           nz01=nknots01,
                                           nz02=nknots02,
                                           nz12=nknots12,
                                           ve01=ve01,
                                           ve02=ve02,
                                           ve12=ve12,
                                           dimnva01=dimnva01,
                                           dimnva02=dimnva02,
                                           dimnva12=dimnva12,
                                           nva01=nvat01,
                                           nva02=nvat02,
                                           nva12=nvat12,
                                           t0=t0,
                                           t1=t1,
                                           t2=t2,
                                           t3=t3,
                                           troncature=troncature,
                                           y01=y01k,
                                           y02=y02k,
                                           y12=y12k,
                                           p01=p01,
                                           p02=p02,
                                           p12=p12,
                                           dimp01=dimp01,
                                           dimp02=dimp02,
                                           dimp12=dimp12,
                                           Ntime=NtimePoints)
                                         
                                         
                                         
                                         #sparseHessianFD
                                         # h=1e-16 --> problem in hessian evaluation
                                         # h=1e-8 --> no issue
                                         output<-output$v
                                         
                                         
                                         if(ite==0){
                                           fn.value<-gaussDYNidmlLikelihoodpena(b=b,
                                                                                npm=npm,
                                                                                npar=size_V,
                                                                                bfix=bfix,
                                                                                fix=fix0,
                                                                                zi01=knots01,
                                                                                zi02=knots02,
                                                                                zi12=knots12,
                                                                                ctime=ctime,
                                                                                no=N,
                                                                                nz01=nknots01,
                                                                                nz02=nknots02,
                                                                                nz12=nknots12,
                                                                                ve01=ve01,
                                                                                ve02=ve02,
                                                                                ve12=ve12,
                                                                                dimnva01=dimnva01,
                                                                                dimnva02=dimnva02,
                                                                                dimnva12=dimnva12,
                                                                                nva01=nvat01,
                                                                                nva02=nvat02,
                                                                                nva12=nvat12,
                                                                                t0=t0,
                                                                                t1=t1,
                                                                                t2=t2,
                                                                                t3=t3,
                                                                                troncature=troncature,
                                                                                y01=y01k,
                                                                                y02=y02k,
                                                                                y12=y12k,
                                                                                p01=p01,
                                                                                p02=p02,
                                                                                p12=p12,
                                                                                dimp01=dimp01,
                                                                                dimp02=dimp02,
                                                                                dimp12=dimp12,
                                                                                Ntime=NtimePoints,
                                                                                
                                                                                lambda=lambda[id.lambda,],
                                                                                alpha=alpha,
                                                                                penalty.factor=penalty.factor,
                                                                                penalty=penalty)
                                         }
                                         
                                         
                                         if(any(is.na(output))|any(output==Inf) |any(output==-Inf)){
                                           warning("Computational error for calculation of the hessian : division by 0 or Infinite value")
                                           if(ite==0){
                                             
                                             fu<-output[(min+1):length(output)]
                                             V<-matrix(0,nrow=npm,ncol=npm)
                                             V[upper.tri(V,diag=T)]<-output[1:min]
                                             V<-V+t(V)
                                             diag(V)<-diag(V)/2
                                             # deriva gives information matrix
                                             tr <- sum(diag(V))/npm
                                             V0<-V}
                                           ite<-ite+1
                                           pbr_compu<-1
                                           break
                                         }
                                         
                                         
                                         fu<-output[(min+1):length(output)]
                                         V<-matrix(0,nrow=npm,ncol=npm)
                                         V[upper.tri(V,diag=T)]<-output[1:min]
                                         V<-V+t(V)
                                         diag(V)<-diag(V)/2
                                         
                                         # projection V into space of eigen vectors
                                         # tol <- 1e-12
                                         # eig <- eigen(V, symmetric = TRUE)
                                         # vals <- eig$values
                                         # vecs <- eig$vectors
                                         # 
                                         # # tronquer les valeurs propres trop petites
                                         # vals_clipped <- pmax(vals, tol)
                                         # 
                                         # # reconstruire la matrice projetée
                                         # H_proj <- vecs %*% diag(vals_clipped) %*% t(vecs)
                                         # 
                                         # 
                                         # pp<-matrix(0,nrow=npm,ncol=npm)
                                         # diag(pp)<-diag(V)
                                         # deriva gives information matrix
                                         tr <- sum(diag(V))/npm
                                         V0<-V
                                         
                                         
                                         eigen.values<-eigen(V,symmetric=T,only.values=T)$values
                                         
                                         idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                         
                                         
                                         idpos0<-idpos
                                         
                                         ncount<-da<-ga<-0
                                         
                                         while(idpos != 0){
                                           
                                           if(ncount==0){ 
                                             ga <- 0.01
                                             da <- 1E-2
                                           }else{
                                             if(((ncount <= 3) | (ga >= 1)) ){
                                               da <- da * 5
                                             }else{# if ncount > 10 only update ga 
                                               ga <- ga * 5
                                               # do not put ga at 1 as no countmax otherwise infinite while 
                                               if(ga > 1) ga <- 1
                                             }
                                           }
                                           
                                           ncount <- ncount + 1
                                           
                                           diagV <- diag(V)
                                           # put abs (1-ga) better than 1-ga cause ga can now be >1
                                           diagV<-ifelse(diagV!=0,diagV+da*(abs((1.e0-ga))*abs(diagV)+ga*tr),
                                                         da*ga*tr)
                                           
                                           diag(V)<-diagV
                                           # if we have a convex log-vraisemblance in eta then :
                                           # all eigen  values of the hessienne are >0.
                                           
                                           if(sum(V==Inf)>0|sum(V==-Inf)>0){break}
                                           eigen.values<-eigen(V,symmetric=T,only.values=T)$values
                                           # check if hessienne defined positive
                                           idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                           
                                           
                                           
                                         }
                                         
                                         
                                         if(idpos!=0){
                                           
                                           warning("Hessian not defined positive")
                                           pbr_compu<-2
                                           ite<-ite+1
                                           break
                                         }
                                         
                                         
                                         # update beta
                                         
                                         output.cv<-DYNcv.model(beta=beta,
                                                                nva01=npm01,
                                                                nva02=npm02,
                                                                nva12=npm12,
                                                                nva01Y=npm01Y,
                                                                nva02Y=npm02Y,
                                                                nva12Y=npm12Y,
                                                                fix=fix0[(size_spline+1):size_V],
                                                                penalty.factor=penalty.factor,
                                                                penalty=penalty,
                                                                v=V,
                                                                fu=fu,
                                                                lambda=lambda[id.lambda,],
                                                                alpha=alpha
                                         )
                                         
                                         # verify validity of parameters update 
                                         # and that we are better than previous estimates 
                                         
                                         b<-c(s,output.cv$b)
                                         
                                         betanew<-b[(size_spline+1):size_V]
                                         
                                         # penalised loglik see if inferior to previous
                                         res<-gaussDYNidmlLikelihoodpena(b=b,
                                                                         npm=size_V,
                                                                         npar=size_V,
                                                                         bfix=1,
                                                                         fix=rep(0,size_V),
                                                                         zi01=knots01,
                                                                         zi02=knots02,
                                                                         zi12=knots12,
                                                                         ctime=ctime,
                                                                         no=N,
                                                                         nz01=nknots01,
                                                                         nz02=nknots02,
                                                                         nz12=nknots12,
                                                                         ve01=ve01,
                                                                         ve02=ve02,
                                                                         ve12=ve12,
                                                                         dimnva01=dimnva01,
                                                                         dimnva02=dimnva02,
                                                                         dimnva12=dimnva12,
                                                                         nva01=nvat01,
                                                                         nva02=nvat02,
                                                                         nva12=nvat12,
                                                                         t0=t0,
                                                                         t1=t1,
                                                                         t2=t2,
                                                                         t3=t3,
                                                                         troncature=troncature,
                                                                         y01=y01k,
                                                                         y02=y02k,
                                                                         y12=y12k,
                                                                         p01=p01,
                                                                         p02=p02,
                                                                         p12=p12,
                                                                         dimp01=dimp01,
                                                                         dimp02=dimp02,
                                                                         dimp12=dimp12,
                                                                         Ntime=NtimePoints,
                                                                         lambda=lambda[id.lambda,],
                                                                         alpha=alpha,
                                                                         penalty.factor=penalty.factor,
                                                                         penalty=penalty)
                                         
                                         
                                         # we want to maximise the loglik thus : 
                                         # we have issue if res is NA or if not higher than previous one 
                                         # if not better or do not exist need to readjust
                                         # value of beta 
                                         
                                         
                                         if(res %in%c(-1e9,1e9) | res < fn.value){
                                           
                                           print(paste0("needed update at ite :",ite))
                                           th<-1e-5
                                           step<-log(1.5)
                                           delta<-output.cv$b-c(beta)
                                           
                                           maxt <- max(abs(delta)) 
                                           
                                           if(maxt == 0){
                                             vw <- th
                                           }else{
                                             vw <- th/maxt
                                           }
                                           if(ite>0){
                                             res.out.error <- list("old.b"=round(c(s,beta)),
                                                                   "old.rl"=round(fn.value),
                                                                   "old.ca"=round(eval.cv.beta[ite]),
                                                                   "old.cb"=round(eval.cv.loglik[ite]))
                                           }else{
                                             res.out.error <- list("old.b"=round(c(s,beta)),
                                                                   "old.rl"=round(fn.value),
                                                                   "old.ca"=round(1),
                                                                   "old.cb"=round(1))
                                           }
                                           m<-1
                                           betanew<-beta
                                           # from mla package
                                           
                                           sears<-searpas(vw=vw,
                                                          step=step,
                                                          b=beta,
                                                          delta=delta,
                                                          funcpa=gaussDYNidmlLikelihoodpena,
                                                          res.out.error=res.out.error,
                                                          npm=npm,
                                                          npar=size_V,
                                                          bfix=s,
                                                          fix=fix0,
                                                          zi01=knots01,
                                                          zi02=knots02,
                                                          zi12=knots12,
                                                          ctime=ctime,
                                                          no=N,
                                                          nz01=nknots01,
                                                          nz02=nknots02,
                                                          nz12=nknots12,
                                                          ve01=ve01,
                                                          ve02=ve02,
                                                          ve12=ve12,
                                                          dimnva01=dimnva01,
                                                          dimnva02=dimnva02,
                                                          dimnva12=dimnva12,
                                                          nva01=nvat01,
                                                          nva02=nvat02,
                                                          nva12=nvat12,
                                                          t0=t0,
                                                          t1=t1,
                                                          t2=t2,
                                                          t3=t3,
                                                          troncature=troncature,
                                                          y01=y01k,
                                                          y02=y02k,
                                                          y12=y12k,
                                                          p01=p01,
                                                          p02=p02,
                                                          p12=p12,
                                                          dimp01=dimp01,
                                                          dimp02=dimp02,
                                                          dimp12=dimp12,
                                                          Ntime=NtimePoints,
                                                          lambda=lambda[id.lambda,],
                                                          alpha=alpha,
                                                          penalty.factor=penalty.factor,
                                                          penalty=penalty)
                                           
                                           
                                           betanew<-beta+delta*sears$vw
                                           b<-c(s,betanew)
                                           
                                           
                                           res<-gaussDYNidmlLikelihoodpena(b=b,
                                                                           npm=size_V,
                                                                           npar=size_V,
                                                                           bfix=1,
                                                                           fix=rep(0,size_V),
                                                                           zi01=knots01,
                                                                           zi02=knots02,
                                                                           zi12=knots12,
                                                                           ctime=ctime,
                                                                           no=N,
                                                                           nz01=nknots01,
                                                                           nz02=nknots02,
                                                                           nz12=nknots12,
                                                                           ve01=ve01,
                                                                           ve02=ve02,
                                                                           ve12=ve12,
                                                                           dimnva01=dimnva01,
                                                                           dimnva02=dimnva02,
                                                                           dimnva12=dimnva12,
                                                                           nva01=nvat01,
                                                                           nva02=nvat02,
                                                                           nva12=nvat12,
                                                                           t0=t0,
                                                                           t1=t1,
                                                                           t2=t2,
                                                                           t3=t3,
                                                                           troncature=troncature,
                                                                           y01=y01k,
                                                                           y02=y02k,
                                                                           y12=y12k,
                                                                           p01=p01,
                                                                           p02=p02,
                                                                           p12=p12,
                                                                           dimp01=dimp01,
                                                                           dimp02=dimp02,
                                                                           dimp12=dimp12,
                                                                           Ntime=NtimePoints,
                                                                           lambda=lambda[id.lambda,],
                                                                           alpha=alpha,
                                                                           penalty.factor=penalty.factor,
                                                                           penalty=penalty)
                                           
                                         }
                                         # if not better or do not exist need to readjust
                                         # value of beta 
                                         if(res %in%c(-1e9,1e9) | any(is.infinite(c(s,betanew)))){
                                           
                                           ite<-ite+1
                                           validity<-F
                                           eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                           eval.validity[ite]<-validity
                                           # save output over iterations 
                                           pbr_compu<-3
                                           break
                                         }else{validity<-T}
                                         
                                         
                                         # betanew already include s
                                         b<-c(s,betanew)
                                         
                                         bfix<-b[fix0.beta==1]
                                         b<-b[fix0.beta==0]
                                         # update modelPar
                                         
                                         
                                         output.mla<- marqLevAlg::mla(b=b,
                                                                      fn=gaussDYNidmlLikelihood,
                                                                      epsa=epsa,
                                                                      epsb=epsb,
                                                                      epsd=epsd,
                                                                      maxiter=maxiter.pena,
                                                                      minimize=F,
                                                                      npm=npmspline,
                                                                      npar=size_V,
                                                                      bfix=bfix,
                                                                      fix=fix0.beta,
                                                                      zi01=knots01,
                                                                      zi02=knots02,
                                                                      zi12=knots12,
                                                                      ctime=ctime,
                                                                      no=N,
                                                                      nz01=nknots01,
                                                                      nz02=nknots02,
                                                                      nz12=nknots12,
                                                                      ve01=ve01,
                                                                      ve02=ve02,
                                                                      ve12=ve12,
                                                                      dimnva01=dimnva01,
                                                                      dimnva02=dimnva02,
                                                                      dimnva12=dimnva12,
                                                                      nva01=nvat01,
                                                                      nva02=nvat02,
                                                                      nva12=nvat12,
                                                                      t0=t0,
                                                                      t1=t1,
                                                                      t2=t2,
                                                                      t3=t3,
                                                                      troncature=troncature,
                                                                      y01=y01k,
                                                                      y02=y02k,
                                                                      y12=y12k,
                                                                      p01=p01,
                                                                      p02=p02,
                                                                      p12=p12,
                                                                      dimp01=dimp01,
                                                                      dimp02=dimp02,
                                                                      dimp12=dimp12,
                                                                      Ntime=NtimePoints)
                                         
                                         # look at convergence for each lambda :
                                         
                                         # new values for splines:
                                         snew<-s
                                         snew[fix00[1:size_spline]==0]<-output.mla$b
                                         if(nvat01>0){
                                           b01<-betanew[1:nvat01][penalty.factor[1:nvat01]==1]
                                           if(p01>0){
                                             b01<-c(b01,betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1])
                                           }
                                         }else{
                                           if(p01>0){
                                             b01<-betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1]
                                           }else{
                                             b01<-0
                                           }
                                         }
                                         
                                         if(nvat02>0){
                                           b02<-betanew[(nvat01+1):(nvat01+nvat02)][penalty.factor[(nvat01+1):(nvat01+nvat02)]==1]
                                           if(p02>0){
                                             b02<-c(b02,betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1])
                                           }
                                         }else{
                                           if(p02>0){
                                             b02<-betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1]
                                           }else{b02<-0}
                                         }
                                         
                                         if(nvat12>0){
                                           b12<-betanew[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)][penalty.factor[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)]==1]
                                           if(p12>0){
                                             b12<-c(b12,betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1])
                                           }
                                         }else{
                                           if(p12>0){
                                             b12<-betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1]
                                           }else{b12<-0}
                                         }
                                         # maximisation issue : lpen =l - pen
                                         if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
                                           fn.valuenew<-output.mla$fn.value-lambda[id.lambda,1]*alpha*sum(abs(b01))-lambda[id.lambda,1]*(1-alpha)*sum(b01*b01)
                                           fn.valuenew<-fn.valuenew-lambda[id.lambda,2]*alpha*sum(abs(b02))-lambda[id.lambda,2]*(1-alpha)*sum(b02*b02)
                                           fn.valuenew<-fn.valuenew-lambda[id.lambda,3]*alpha*sum(abs(b12))-lambda[id.lambda,3]*(1-alpha)*sum(b12*b12)
                                           
                                         }
                                         
                                         if(penalty=="mcp"){
                                           
                                           p01<-rep(alpha*lambda[id.lambda,1]*lambda[id.lambda,1]/2,length(b01))
                                           idbeta<-which(b01<=alpha*lambda[id.lambda,1])
                                           p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])-((b01[idbeta]*b01[idbeta])/2*alpha)
                                           
                                           p02<-rep(alpha*lambda[id.lambda,2]*lambda[id.lambda,2]/2,length(b02))
                                           idbeta<-which(b02<=alpha*lambda[id.lambda,2])
                                           p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])-((b02[idbeta]*b02[idbeta])/2*alpha)
                                           
                                           p12<-rep(alpha*lambda[id.lambda,3]*lambda[id.lambda,3]/2,length(b12))
                                           idbeta<-which(b12<=alpha*lambda[id.lambda,3])
                                           p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])-((b12[idbeta]*b12[idbeta])/2*alpha)
                                           
                                           fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                           
                                         }
                                         
                                         if(penalty=="scad"){
                                           
                                           p01<-rep((lambda[id.lambda,1]^2)*(alpha+1)/2,length(b01))
                                           idbeta<-which(b01<=lambda[id.lambda,1])
                                           p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])
                                           idbeta<-which(abs(b01)<lambda[id.lambda,1]*alpha)
                                           p01[idbeta]<-(2*alpha*lambda[id.lambda,1]*abs(b01[idbeta])-b01[idbeta]^2-lambda[id.lambda,1]^2)/(2*(alpha-1))
                                           
                                           p02<-rep((lambda[id.lambda,2]^2)*(alpha+1)/2,length(b02))
                                           idbeta<-which(b02<=lambda[id.lambda,2])
                                           p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])
                                           idbeta<-which(abs(b02)<lambda[id.lambda,2]*alpha)
                                           p02[idbeta]<-(2*alpha*lambda[id.lambda,2]*abs(b02[idbeta])-b02[idbeta]^2-lambda[id.lambda,2]^2)/(2*(alpha-1))
                                           
                                           p12<-rep((lambda[id.lambda,3]^2)*(alpha+1)/2,length(b12))
                                           idbeta<-which(b12<=lambda[id.lambda,3])
                                           p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])
                                           idbeta<-which(abs(b12)<lambda[id.lambda,3]*alpha)
                                           p12[idbeta]<-(2*alpha*lambda[id.lambda,3]*abs(b12[idbeta])-b12[idbeta]^2-lambda[id.lambda,3]^2)/(2*(alpha-1))
                                           
                                           fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                           
                                         }
                                         
                                         ite<-ite+1
                                         
                                         # check cv criterias 
                                         
                                         eval.cv.spline[ite]<-sum((snew-s)^2)
                                         eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                         eval.cv.loglik[ite]<-abs((fn.valuenew-fn.value)/fn.value)
                                         eval.loglik[ite]<-fn.valuenew
                                         eval.validity[ite]<-validity
                                         
                                         s<-snew
                                         beta<-betanew
                                         fn.value<-fn.valuenew
                                         
                                         # eval.cv beta valid only if validity.param=T
                                         # add ite >1 so that do not stop at first ite ?
                                         if(((eval.cv.beta[ite]+eval.cv.spline[ite])<epsa) & eval.cv.loglik[ite]<epsb & validity==T){
                                           converged<-T}
                                         
                                         
                                       }
                                       
                                       if(maxiter<=ite & converged==F){
                                         istop<-2
                                       }else{
                                         if(ite<=maxiter & converged==T){
                                           istop<-1
                                           # need to recalculate second derivatives 
                                           
                                           b<-c(s,beta)
                                           bfix<-b[fix0==1]
                                           b<-b[fix0==0]
                                           
                                           output<-DYNderiva(
                                             h=1e-8,
                                             b=b,
                                             npm=npm,
                                             npar=size_V,
                                             bfix=bfix,
                                             fix=fix0,
                                             zi01=knots01,
                                             zi02=knots02,
                                             zi12=knots12,
                                             ctime=ctime,
                                             no=N,
                                             nz01=nknots01,
                                             nz02=nknots02,
                                             nz12=nknots12,
                                             ve01=ve01,
                                             ve02=ve02,
                                             ve12=ve12,
                                             dimnva01=dimnva01,
                                             dimnva02=dimnva02,
                                             dimnva12=dimnva12,
                                             nva01=nvat01,
                                             nva02=nvat02,
                                             nva12=nvat12,
                                             t0=t0,
                                             t1=t1,
                                             t2=t2,
                                             t3=t3,
                                             troncature=troncature,
                                             y01=y01k,
                                             y02=y02k,
                                             y12=y12k,
                                             p01=p01,
                                             p02=p02,
                                             p12=p12,
                                             dimp01=dimp01,
                                             dimp02=dimp02,
                                             dimp12=dimp12,
                                             Ntime=NtimePoints)
                                           
                                           output<-output$v
                                           fu<-output[(min+1):length(output)]
                                           V<-matrix(0,nrow=npm,ncol=npm)
                                           V[upper.tri(V,diag=T)]<-output[1:min]
                                           V<-V+t(V)
                                           diag(V)<-diag(V)/2
                                           # deriva gives information matrix
                                           V0<-V
                                         }else{
                                           if(pbr_compu==1){istop<-3}
                                           if(pbr_compu==2){istop<-4}
                                           if(pbr_compu==3){istop<-5}
                                         }
                                       }
                                       
                                       # if stop==1 we can give matrix of second derivatives 
                                       
                                       combine<-combine+1
                                       return(list(b=c(s,beta),
                                                   H=V0,
                                                   lambda=as.double(lambda[id.lambda,]),
                                                   alpha=alpha,
                                                   fn.value=ifelse(!exists("output.mla"),NA,output.mla$fn.value),
                                                   fn.value.pena=fn.value,
                                                   ni=ite,
                                                   ca.beta=eval.cv.beta,
                                                   ca.spline=eval.cv.spline,
                                                   ca.validity=eval.validity,
                                                   cb=eval.loglik,
                                                   istop=istop,
                                                   combine=combine))
                                     }
          }else{outputNsample<-foreach::foreach(id.lambda=1:nlambda,
                                         .combine = combine_lambda_mla,
                                         .errorhandling = "remove")%do%{
                                           
                                           
                                           # computation pbr 
                                           
                                           pbr_compu<-0
                                           
                                           
                                             beta<-beta.start
                                             s<-s.start
                                           
                                           
                                           converged<-F
                                           ite<-0
                                           # if beta not change do not need to recalculate weights 
                                           H<-T
                                           
                                           eval.cv.spline<-rep(NA,maxiter+1)
                                           eval.cv.beta<-rep(NA,maxiter+1)
                                           eval.cv.loglik<-rep(NA,maxiter+1)
                                           eval.loglik<-rep(NA,maxiter+1)
                                           eval.validity<-rep(NA,maxiter+1)
                                           
                                           
                                           while(converged==F & ite<=maxiter){
                                             
                                             b<-c(s,beta)
                                             bfix<-b[fix0==1]
                                             b<-b[fix0==0]
                                             
                                             
                                             
                                             output<-DYNderivadiag(
                                               h=1e-8,
                                               b=b,
                                               npm=npm,
                                               npar=size_V,
                                               bfix=bfix,
                                               fix=fix0,
                                               zi01=knots01,
                                               zi02=knots02,
                                               zi12=knots12,
                                               ctime=ctime,
                                               no=N,
                                               nz01=nknots01,
                                               nz02=nknots02,
                                               nz12=nknots12,
                                               ve01=ve01,
                                               ve02=ve02,
                                               ve12=ve12,
                                               dimnva01=dimnva01,
                                               dimnva02=dimnva02,
                                               dimnva12=dimnva12,
                                               nva01=nvat01,
                                               nva02=nvat02,
                                               nva12=nvat12,
                                               t0=t0,
                                               t1=t1,
                                               t2=t2,
                                               t3=t3,
                                               troncature=troncature,
                                               y01=y01k,
                                               y02=y02k,
                                               y12=y12k,
                                               p01=p01,
                                               p02=p02,
                                               p12=p12,
                                               dimp01=dimp01,
                                               dimp02=dimp02,
                                               dimp12=dimp12,
                                               Ntime=NtimePoints)
                                             
                                             
                                             output<-output$v
                                             if(ite==0){
                                               fn.value<-gaussDYNidmlLikelihoodpena(b=b,
                                                                                    npm=npm,
                                                                                    npar=size_V,
                                                                                    bfix=bfix,
                                                                                    fix=fix0,
                                                                                    zi01=knots01,
                                                                                    zi02=knots02,
                                                                                    zi12=knots12,
                                                                                    ctime=ctime,
                                                                                    no=N,
                                                                                    nz01=nknots01,
                                                                                    nz02=nknots02,
                                                                                    nz12=nknots12,
                                                                                    ve01=ve01,
                                                                                    ve02=ve02,
                                                                                    ve12=ve12,
                                                                                    dimnva01=dimnva01,
                                                                                    dimnva02=dimnva02,
                                                                                    dimnva12=dimnva12,
                                                                                    nva01=nvat01,
                                                                                    nva02=nvat02,
                                                                                    nva12=nvat12,
                                                                                    t0=t0,
                                                                                    t1=t1,
                                                                                    t2=t2,
                                                                                    t3=t3,
                                                                                    troncature=troncature,
                                                                                    y01=y01k,
                                                                                    y02=y02k,
                                                                                    y12=y12k,
                                                                                    p01=p01,
                                                                                    p02=p02,
                                                                                    p12=p12,
                                                                                    dimp01=dimp01,
                                                                                    dimp02=dimp02,
                                                                                    dimp12=dimp12,
                                                                                    Ntime=NtimePoints,
                                                                                    
                                                                                    lambda=lambda[id.lambda,],
                                                                                    alpha=alpha,
                                                                                    penalty.factor=penalty.factor,
                                                                                    penalty=penalty)
                                             }
                                             
                                             if(any(is.na(output$V))|any(output$V==Inf) |any(output$V==-Inf)|any(is.na(output$fu))|any(output$fu==Inf) |any(output$fu==-Inf)){
                                               warning("Computational error for calculation of the hessian : division by 0 or Infinite value")
                                               if(ite==0){
                                                 
                                                 fu<-output[(min+1):(2*min)]
                                                 V<-matrix(0,nrow=npm,ncol=npm)
                                                 diag(V)<-output[1:min]
                                                 # deriva gives information matrix
                                                 tr <- sum(diag(V))/npm
                                                 V0<-V}
                                               ite<-ite+1
                                               pbr_compu<-1
                                               break
                                             }
                                             
                                             
                                             fu<-output[(min+1):(2*min)]
                                             V<-matrix(0,nrow=npm,ncol=npm)
                                             diag(V)<-output[1:min]
                                             # deriva gives information matrix
                                             tr <- sum(diag(V))/npm
                                             V0<-V
                                             
                                             eigen.values<-diag(V)
                                             
                                             idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                             
                                             
                                             idpos0<-idpos
                                             
                                             ncount<-da<-ga<-0
                                             
                                             
                                             while(idpos != 0){
                                               
                                               if(ncount==0){ 
                                                 ga <- 0.01
                                                 da <- 1E-2
                                               }else{
                                                 if(((ncount <= 3) | (ga >= 1)) ){
                                                   da <- da * 5
                                                 }else{# if ncount > 10 only update ga 
                                                   ga <- ga * 5
                                                   # do not put ga at 1 as no countmax otherwise infinite while 
                                                   if(ga > 1) ga <- 1
                                                 }
                                               }
                                               
                                               ncount <- ncount + 1
                                               
                                               diagV <- diag(V)
                                               # put abs (1-ga) better than 1-ga cause ga can now be >1
                                               diagV<-ifelse(diagV!=0,diagV+da*(abs((1.e0-ga))*abs(diagV)+ga*tr),
                                                             da*ga*tr)
                                               
                                               diag(V)<-diagV
                                               # if we have a convex log-vraisemblance in eta then :
                                               # all eigen  values of the hessienne are >0.
                                               
                                               if(sum(V==Inf)>0|sum(V==-Inf)>0){break}
                                               eigen.values<-diag(V)
                                               # check if hessienne defined positive
                                               idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                               
                                               
                                               
                                             }
                                             
                                             if(idpos!=0){
                                               
                                               warning("Hessian not defined positive")
                                               pbr_compu<-2
                                               ite<-ite+1
                                               break
                                             }
                                             
                                             
                                             # update beta
                                             output.cv<-DYNcv.model(beta=beta,
                                                                    nva01=npm01,
                                                                    nva02=npm02,
                                                                    nva12=npm12,
                                                                    nva01Y=npm01Y,
                                                                    nva02Y=npm02Y,
                                                                    nva12Y=npm12Y,
                                                                    fix=fix0[(size_spline+1):size_V],
                                                                    penalty.factor=penalty.factor,
                                                                    penalty=penalty,
                                                                    v=V,
                                                                    fu=fu,
                                                                    lambda=lambda[id.lambda,],
                                                                    alpha=alpha
                                             )
                                             
                                             # verify validity of parameters update 
                                             # and that we are better than previous estimates 
                                             
                                             b<-c(s,output.cv$b)
                                             
                                             betanew<-b[(size_spline+1):size_V]
                                             
                                             # penalised loglik see if inferior to previous
                                             res<-gaussDYNidmlLikelihoodpena(b=b,
                                                                             npm=size_V,
                                                                             npar=size_V,
                                                                             bfix=1,
                                                                             fix=rep(0,size_V),
                                                                             zi01=knots01,
                                                                             zi02=knots02,
                                                                             zi12=knots12,
                                                                             ctime=ctime,
                                                                             no=N,
                                                                             nz01=nknots01,
                                                                             nz02=nknots02,
                                                                             nz12=nknots12,
                                                                             ve01=ve01,
                                                                             ve02=ve02,
                                                                             ve12=ve12,
                                                                             dimnva01=dimnva01,
                                                                             dimnva02=dimnva02,
                                                                             dimnva12=dimnva12,
                                                                             nva01=nvat01,
                                                                             nva02=nvat02,
                                                                             nva12=nvat12,
                                                                             t0=t0,
                                                                             t1=t1,
                                                                             t2=t2,
                                                                             t3=t3,
                                                                             troncature=troncature,
                                                                             y01=y01k,
                                                                             y02=y02k,
                                                                             y12=y12k,
                                                                             p01=p01,
                                                                             p02=p02,
                                                                             p12=p12,
                                                                             dimp01=dimp01,
                                                                             dimp02=dimp02,
                                                                             dimp12=dimp12,
                                                                             Ntime=NtimePoints,
                                                                             
                                                                             lambda=lambda[id.lambda,],
                                                                             alpha=alpha,
                                                                             penalty.factor=penalty.factor,
                                                                             penalty=penalty)
                                             
                                             # we want to maximise the loglik thus : 
                                             # we have issue if res is NA or if not higher than previous one 
                                             # if not better or do not exist need to readjust
                                             # value of beta 
                                             
                                             if(res %in%c(-1e9,1e9) | res < fn.value){
                                               
                                               th<-1e-5
                                               step<-log(1.5)
                                               delta<-output.cv$b-c(beta)
                                               
                                               maxt <- max(abs(delta)) 
                                               
                                               if(maxt == 0){
                                                 vw <- th
                                               }else{
                                                 vw <- th/maxt
                                               }
                                               if(ite>0){
                                                 res.out.error <- list("old.b"=round(c(s,beta)),
                                                                       "old.rl"=round(fn.value),
                                                                       "old.ca"=round(eval.cv.beta[ite]),
                                                                       "old.cb"=round(eval.cv.loglik[ite]))
                                               }else{
                                                 res.out.error <- list("old.b"=round(c(s,beta)),
                                                                       "old.rl"=round(fn.value),
                                                                       "old.ca"=round(1),
                                                                       "old.cb"=round(1))
                                               }
                                               m<-1
                                               betanew<-beta
                                               # from mla package
                                               
                                               sears<-searpas(vw=vw,
                                                              step=step,
                                                              b=beta,
                                                              delta=delta,
                                                              funcpa=gaussDYNidmlLikelihoodpena,
                                                              res.out.error=res.out.error,
                                                              npm=npm,
                                                              npar=size_V,
                                                              bfix=s,
                                                              fix=fix0,
                                                              zi01=knots01,
                                                              zi02=knots02,
                                                              zi12=knots12,
                                                              ctime=ctime,
                                                              no=N,
                                                              nz01=nknots01,
                                                              nz02=nknots02,
                                                              nz12=nknots12,
                                                              ve01=ve01,
                                                              ve02=ve02,
                                                              ve12=ve12,
                                                              dimnva01=dimnva01,
                                                              dimnva02=dimnva02,
                                                              dimnva12=dimnva12,
                                                              nva01=nvat01,
                                                              nva02=nvat02,
                                                              nva12=nvat12,
                                                              t0=t0,
                                                              t1=t1,
                                                              t2=t2,
                                                              t3=t3,
                                                              troncature=troncature,
                                                              y01=y01k,
                                                              y02=y02k,
                                                              y12=y12k,
                                                              p01=p01,
                                                              p02=p02,
                                                              p12=p12,
                                                              dimp01=dimp01,
                                                              dimp02=dimp02,
                                                              dimp12=dimp12,
                                                              Ntime=NtimePoints,
                                                              lambda=lambda[id.lambda,],
                                                              alpha=alpha,
                                                              penalty.factor=penalty.factor,
                                                              penalty=penalty)
                                               
                                               
                                               betanew<-beta+delta*sears$vw
                                               b<-c(s,betanew)
                                               
                                               
                                               res<-gaussDYNidmlLikelihoodpena(b=b,
                                                                               npm=size_V,
                                                                               npar=size_V,
                                                                               bfix=1,
                                                                               fix=rep(0,size_V),
                                                                               zi01=knots01,
                                                                               zi02=knots02,
                                                                               zi12=knots12,
                                                                               ctime=ctime,
                                                                               no=N,
                                                                               nz01=nknots01,
                                                                               nz02=nknots02,
                                                                               nz12=nknots12,
                                                                               ve01=ve01,
                                                                               ve02=ve02,
                                                                               ve12=ve12,
                                                                               dimnva01=dimnva01,
                                                                               dimnva02=dimnva02,
                                                                               dimnva12=dimnva12,
                                                                               nva01=nvat01,
                                                                               nva02=nvat02,
                                                                               nva12=nvat12,
                                                                               t0=t0,
                                                                               t1=t1,
                                                                               t2=t2,
                                                                               t3=t3,
                                                                               troncature=troncature,
                                                                               y01=y01k,
                                                                               y02=y02k,
                                                                               y12=y12k,
                                                                               p01=p01,
                                                                               p02=p02,
                                                                               p12=p12,
                                                                               dimp01=dimp01,
                                                                               dimp02=dimp02,
                                                                               dimp12=dimp12,
                                                                               Ntime=NtimePoints,
                                                                               
                                                                               lambda=lambda[id.lambda,],
                                                                               alpha=alpha,
                                                                               penalty.factor=penalty.factor,
                                                                               penalty=penalty)
                                               
                                             }
                                             # if not better or do not exist need to readjust
                                             # value of beta 
                                             if(res %in%c(-1e9,1e9) | any(is.infinite(c(s,betanew)))){
                                               
                                               ite<-ite+1
                                               validity<-F
                                               eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                               eval.validity[ite]<-validity
                                               # save output over iterations 
                                               pbr_compu<-3
                                               break
                                             }else{validity<-T}
                                             
                                             
                                             # betanew already include s
                                             b<-c(s,betanew)
                                             
                                             bfix<-b[fix0.beta==1]
                                             b<-b[fix0.beta==0]
                                             # update modelPar
                                             
                                             
                                             output.mla<- marqLevAlg::mla(b=b,
                                                                          fn=gaussDYNidmlLikelihood,
                                                                          epsa=epsa,
                                                                          epsb=epsb,
                                                                          epsd=epsd,
                                                                          maxiter=maxiter.pena,
                                                                          minimize=F,
                                                                          npm=npmspline,
                                                                          npar=size_V,
                                                                          bfix=bfix,
                                                                          fix=fix0.beta,
                                                                          zi01=knots01,
                                                                          zi02=knots02,
                                                                          zi12=knots12,
                                                                          ctime=ctime,
                                                                          no=N,
                                                                          nz01=nknots01,
                                                                          nz02=nknots02,
                                                                          nz12=nknots12,
                                                                          ve01=ve01,
                                                                          ve02=ve02,
                                                                          ve12=ve12,
                                                                          dimnva01=dimnva01,
                                                                          dimnva02=dimnva02,
                                                                          dimnva12=dimnva12,
                                                                          nva01=nvat01,
                                                                          nva02=nvat02,
                                                                          nva12=nvat12,
                                                                          t0=t0,
                                                                          t1=t1,
                                                                          t2=t2,
                                                                          t3=t3,
                                                                          troncature=troncature,
                                                                          y01=y01k,
                                                                          y02=y02k,
                                                                          y12=y12k,
                                                                          p01=p01,
                                                                          p02=p02,
                                                                          p12=p12,
                                                                          dimp01=dimp01,
                                                                          dimp02=dimp02,
                                                                          dimp12=dimp12,
                                                                          Ntime=NtimePoints)
                                             
                                             # look at convergence for each lambda :
                                             
                                             # new values for splines:
                                             snew<-s
                                             snew[fix00[1:size_spline]==0]<-output.mla$b
                                             if(nvat01>0){
                                               b01<-betanew[1:nvat01][penalty.factor[1:nvat01]==1]
                                               if(p01>0){
                                                 b01<-c(b01,betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1])
                                               }
                                             }else{
                                               if(p01>0){
                                                 b01<-betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1]
                                               }else{
                                                 b01<-0
                                               }
                                             }
                                             
                                             if(nvat02>0){
                                               b02<-betanew[(nvat01+1):(nvat01+nvat02)][penalty.factor[(nvat01+1):(nvat01+nvat02)]==1]
                                               if(p02>0){
                                                 b02<-c(b02,betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1])
                                               }
                                             }else{
                                               if(p02>0){
                                                 b02<-betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1]
                                               }else{b02<-0}
                                             }
                                             
                                             if(nvat12>0){
                                               b12<-betanew[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)][penalty.factor[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)]==1]
                                               if(p12>0){
                                                 b12<-c(b12,betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1])
                                               }
                                             }else{
                                               if(p12>0){
                                                 b12<-betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1]
                                               }else{b12<-0}
                                             }
                                             # maximisation issue : lpen =l - pen
                                             if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
                                               fn.valuenew<-output.mla$fn.value-lambda[id.lambda,1]*alpha*sum(abs(b01))-lambda[id.lambda,1]*(1-alpha)*sum(b01*b01)
                                               fn.valuenew<-fn.valuenew-lambda[id.lambda,2]*alpha*sum(abs(b02))-lambda[id.lambda,2]*(1-alpha)*sum(b02*b02)
                                               fn.valuenew<-fn.valuenew-lambda[id.lambda,3]*alpha*sum(abs(b12))-lambda[id.lambda,3]*(1-alpha)*sum(b12*b12)
                                               
                                             }
                                             
                                             if(penalty=="mcp"){
                                               
                                               p01<-rep(alpha*lambda[id.lambda,1]*lambda[id.lambda,1]/2,length(b01))
                                               idbeta<-which(b01<=alpha*lambda[id.lambda,1])
                                               p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])-((b01[idbeta]*b01[idbeta])/2*alpha)
                                               
                                               p02<-rep(alpha*lambda[id.lambda,2]*lambda[id.lambda,2]/2,length(b02))
                                               idbeta<-which(b02<=alpha*lambda[id.lambda,2])
                                               p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])-((b02[idbeta]*b02[idbeta])/2*alpha)
                                               
                                               p12<-rep(alpha*lambda[id.lambda,3]*lambda[id.lambda,3]/2,length(b12))
                                               idbeta<-which(b12<=alpha*lambda[id.lambda,3])
                                               p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])-((b12[idbeta]*b12[idbeta])/2*alpha)
                                               
                                               fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                               
                                             }
                                             
                                             if(penalty=="scad"){
                                               
                                               p01<-rep((lambda[id.lambda,1]^2)*(alpha+1)/2,length(b01))
                                               idbeta<-which(b01<=lambda[id.lambda,1])
                                               p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])
                                               idbeta<-which(abs(b01)<lambda[id.lambda,1]*alpha)
                                               p01[idbeta]<-(2*alpha*lambda[id.lambda,1]*abs(b01[idbeta])-b01[idbeta]^2-lambda[id.lambda,1]^2)/(2*(alpha-1))
                                               
                                               p02<-rep((lambda[id.lambda,2]^2)*(alpha+1)/2,length(b02))
                                               idbeta<-which(b02<=lambda[id.lambda,2])
                                               p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])
                                               idbeta<-which(abs(b02)<lambda[id.lambda,2]*alpha)
                                               p02[idbeta]<-(2*alpha*lambda[id.lambda,2]*abs(b02[idbeta])-b02[idbeta]^2-lambda[id.lambda,2]^2)/(2*(alpha-1))
                                               
                                               p12<-rep((lambda[id.lambda,3]^2)*(alpha+1)/2,length(b12))
                                               idbeta<-which(b12<=lambda[id.lambda,3])
                                               p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])
                                               idbeta<-which(abs(b12)<lambda[id.lambda,3]*alpha)
                                               p12[idbeta]<-(2*alpha*lambda[id.lambda,3]*abs(b12[idbeta])-b12[idbeta]^2-lambda[id.lambda,3]^2)/(2*(alpha-1))
                                               
                                               fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                               
                                             }
                                             
                                             ite<-ite+1
                                             
                                             # check cv criterias 
                                             
                                             eval.cv.spline[ite]<-sum((snew-s)^2)
                                             eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                             eval.cv.loglik[ite]<-abs((fn.valuenew-fn.value)/fn.value)
                                             eval.loglik[ite]<-fn.valuenew
                                             eval.validity[ite]<-validity
                                             
                                             s<-snew
                                             beta<-betanew
                                             fn.value<-fn.valuenew
                                             
                                             # eval.cv beta valid only if validity.param=T
                                             if(((eval.cv.beta[ite]+ eval.cv.spline[ite])<epsa) & eval.cv.loglik[ite]<epsb & validity==T){
                                               converged<-T}
                                             
                                             
                                           }
                                           if(maxiter<=ite & converged==F){
                                             istop<-2
                                           }else{
                                             if(ite<=maxiter & converged==T){
                                               istop<-1
                                               # need to recalculate second derivatives 
                                               
                                               b<-c(s,beta)
                                               bfix<-b[fix0==1]
                                               b<-b[fix0==0]
                                               
                                               output<-DYNderivadiag(
                                                 h=1e-8,
                                                 b=b,
                                                 npm=npm,
                                                 npar=size_V,
                                                 bfix=bfix,
                                                 fix=fix0,
                                                 zi01=knots01,
                                                 zi02=knots02,
                                                 zi12=knots12,
                                                 ctime=ctime,
                                                 no=N,
                                                 nz01=nknots01,
                                                 nz02=nknots02,
                                                 nz12=nknots12,
                                                 ve01=ve01,
                                                 ve02=ve02,
                                                 ve12=ve12,
                                                 dimnva01=dimnva01,
                                                 dimnva02=dimnva02,
                                                 dimnva12=dimnva12,
                                                 nva01=nvat01,
                                                 nva02=nvat02,
                                                 nva12=nvat12,
                                                 t0=t0,
                                                 t1=t1,
                                                 t2=t2,
                                                 t3=t3,
                                                 troncature=troncature,
                                                 y01=y01k,
                                                 y02=y02k,
                                                 y12=y12k,
                                                 p01=p01,
                                                 p02=p02,
                                                 p12=p12,
                                                 dimp01=dimp01,
                                                 dimp02=dimp02,
                                                 dimp12=dimp12,
                                                 Ntime=NtimePoints)
                                               
                                               output<-output$v
                                               
                                               
                                               fu<-output[(min+1):(2*min)]
                                               V<-matrix(0,nrow=npm,ncol=npm)
                                               diag(V)<-output[1:min]
                                               # deriva gives information matrix
                                               V0<-V
                                             }else{
                                               if(pbr_compu==1){istop<-3}
                                               if(pbr_compu==2){istop<-4}
                                               if(pbr_compu==3){istop<-5}
                                             }
                                           }
                                           
                                           # if stop==1 we can give matrix of second derivatives 
                                           
                                           combine<-combine+1
                                           
                                           return(list(b=c(s,beta),
                                                       H=V0,
                                                       lambda=as.double(lambda[id.lambda,]),
                                                       alpha=alpha,
                                                       fn.value=ifelse(!exists("output.mla"),NA,output.mla$fn.value),
                                                       fn.value.pena=fn.value,
                                                       ni=ite,
                                                       ca.beta=eval.cv.beta,
                                                       ca.spline=eval.cv.spline,
                                                       ca.validity=eval.validity,
                                                       cb=eval.loglik,
                                                       istop=istop,
                                                       combine=combine))
                                         }}
        
        
        return(outputNsample)
                               }
      parallel::stopCluster(clustpar)
    }
    
    
    }
  
  
  return(output=outputall)
}



