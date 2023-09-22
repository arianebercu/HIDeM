##' programmation for penalised idm with weibull as risk bases 
##' Code:
##' @title idm penalty weib
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
##' @param troncature0 indicator if troncature or not
##' @param epsa control convergence parameter for beta 
##' @param epsb control convergence parameter for loglik
##' @param epsd control convergence for distance to minimum rdm
##' @param eps.eigen the power of convergence for eigen values of covariance matrix only
##' @param print.info shloud we print info during mla convergence
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
##' @param idd indicator of death
##' @param idm indicator of illness
##' @param ts time of death or last news
##' @param troncature0 indicator if troncature or not
#' @importFrom foreach "%do%"
#' @importFrom foreach "%dopar%"
#' @useDynLib SmoothHazardoptim9


idm.penalty.weib<-function(b,fix0,size_V,
                   clustertype,epsa,epsb,epsd,eps.eigen,print.info,nproc,maxiter,maxiter.pena,
                   ctime,N,
                   ve01,ve02,ve12,dimnva01,dimnva02,dimnva12,nvat01,nvat02,nvat12,
                   t0,t1,t2,t3,idd,idm,ts,troncature,
                   nlambda01,lambda01,nlambda02,lambda02,nlambda12,lambda12,
                   alpha,penalty.factor,penalty){
  
  
  # need to keep original fix to calculate for beta 
  
  V0<-NA
  fix00<-fix0
  
  # create grid 3
  lambda<-expand.grid(lambda01,lambda02,lambda12)
  lambda<-unique(lambda)
  nlambda<-dim(lambda)[1]
  
  # computation pbr 
  pbr_compu<-0
  
  combine_lambda<-function(x,newx){
    
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
  
   s.start<-b[1:6]
  
  # Initiate value of eta : all the same for each lambda
  
  beta.start<-b[(6+1):(size_V)]
  
  
  combine<-0
  
  # fix0 will be used to calculate derivatives and second derivatives only
  # for Beta and not theta01,02,12
  fix0[1:6]<-rep(1,6)
  fix0.beta<-fix00
  fix0.beta[(6+1):size_V]<-rep(1,size_V-6)
  
                

  
  if(nproc >1){
    
    
    if(is.null(clustertype)){
      clustpar <- parallel::makeCluster(nproc)#, outfile="")
    }
    else{
      clustpar <- parallel::makeCluster(nproc, type=clustertype)#, outfile="")
    }
    
    doParallel::registerDoParallel(clustpar)
    
    
    
    
    output<-foreach::foreach(id.lambda=1:nlambda,
                             .combine = combine_lambda,
                             .errorhandling = "remove")%dopar%{
                               
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
                               
                               npm<-sum(fix0==0)
                               
                               npm01<-ifelse(nvat01>0,sum(fix0[7:(7+nvat01-1)]==0),0)
                               npm02<-ifelse(nvat02>0,sum(fix0[(7+nvat01):(6+nvat01+nvat02)]==0),0)
                               npm12<-ifelse(nvat12>0,sum(fix0[(7+nvat01+nvat02):size_V]==0),0)
                               
                               
                               while(converged==F & ite<=maxiter){
                                 
                                 
                                 b<-c(s,beta)
                                 bfix<-b[fix0==1]
                                 b<-b[fix0==0]
                                 output<-derivaweib(b0=b,
                                                    np0=length(b),
                                                    npar0=size_V,
                                                    bfix0=bfix,
                                                    fix0=fix0,
                                                    c0=ctime,
                                                    no0=N,
                                                    ve010=ve01,
                                                    ve020=ve02,
                                                    ve120=ve12,
                                                    dimnva01=dimnva01,
                                                    dimnva02=dimnva02,
                                                    dimnva12=dimnva12,
                                                    nva01=nvat01,
                                                    nva02=nvat02,
                                                    nva12=nvat12,
                                                    t00=t0,
                                                    t10=t1,
                                                    t20=t2,
                                                    t30=t3,
                                                    troncature0=troncature)
                                 
                                 if(ite==0){
                                   fn.value<-idmlLikelihoodweibpena(b=b,
                                                                    np0=length(b),
                                                                    npar0=size_V,
                                                                    bfix0=bfix,
                                                                    fix0=fix0,
                                                                    c0=ctime,
                                                                    no0=N,
                                                                    ve010=ve01,
                                                                    ve020=ve02,
                                                                    ve120=ve12,
                                                                    dimnva01=dimnva01,
                                                                    dimnva02=dimnva02,
                                                                    dimnva12=dimnva12,
                                                                    nva01=nvat01,
                                                                    nva02=nvat02,
                                                                    nva12=nvat12,
                                                                    t00=t0,
                                                                    t10=t1,
                                                                    t20=t2,
                                                                    t30=t3,
                                                                    troncature0=troncature,
                                                                    lambda=lambda[id.lambda,],
                                                                    alpha=alpha,
                                                                    penalty.factor=penalty.factor,
                                                                    penalty=penalty)
                                 }
                                 
                                 if(any(is.na(output))|any(output==Inf) |any(output==-Inf)){
                                   warning("Computational error for calculation of the hessian : division by 0 or Infinite value")
                                   if(ite==0){
                                     fu <- output[1:npm]
                                     min<-npm
                                     V01<- matrix(0,npm01,npm01)
                                     V01[lower.tri(V01,diag=TRUE)] <- output[(min+1):(min+npm01*(npm01+1)/2)]
                                     
                                     
                                     min<-min+(npm01*(npm01+1)/2)
                                     if(npm01>0&npm02>0){
                                       V0102<- matrix(data=output[(min+1):(min+npm02*npm01)],
                                                      nrow=npm02,ncol=npm01)
                                     }
                                     
                                     min<-min+npm02*npm01
                                     
                                     if(npm01>0&npm12>0){
                                       V0112<- matrix(data=output[(min+1):(min+npm12*npm01)],
                                                      nrow=npm12,ncol=npm01)
                                     }
                                     
                                     
                                     min<-min+npm12*npm01
                                     V02<- matrix(0,npm02,npm02)
                                     V02[lower.tri(V02,diag=TRUE)] <- output[(min+1):(min+npm02*(npm02+1)/2)]
                                     
                                     
                                     min<-min+(npm02*(npm02+1)/2)
                                     
                                     if(npm02>0&npm12>0){
                                       V0212<- matrix(data=output[(min+1):(min+npm12*npm02)],
                                                      nrow=npm12,ncol=npm02)}
                                     
                                     
                                     min<-min+npm12*npm02
                                     V12<- matrix(0,npm12,npm12)
                                     V12[lower.tri(V12,diag=TRUE)] <- output[(min+1):length(output)]
                                     
                                     
                                     V<- matrix(0,npm,npm)
                                     if(npm01>0){
                                       V[1:npm01,1:npm01]<-V01
                                       if(npm02>0){V[(npm01+1):(npm01+npm02),1:npm01]<-V0102}
                                       if(npm12>0){V[(npm01+npm02+1):npm,1:npm01]<-V0112}
                                     }
                                     if(npm02>0){
                                       V[(npm01+1):(npm01+npm02),(npm01+1):(npm01+npm02)]<-V02
                                       if(npm12>0){V[(npm01+npm02+1):npm,(npm01+1):(npm01+npm02)]<-V0212}
                                     }
                                     
                                     if(npm12>0){
                                       V[(npm01+npm02+1):npm,(npm01+npm02+1):(npm)]<-V12}
                                     
                                     V<-V+t(V)
                                     diag(V)<-diag(V)/2
                                     # hessian is - second derivatives 
                                     V<--V
                                     
                                     tr <- sum(diag(V))/npm
                                     V0<-V}
                                   ite<-ite+1
                                   pbr_compu<-1
                                   break
                                 }
                                 
                                 fu <- output[1:npm]
                                 
                                 min<-npm
                                 V01<- matrix(0,npm01,npm01)
                                 V01[lower.tri(V01,diag=TRUE)] <- output[(min+1):(min+npm01*(npm01+1)/2)]
                                 
                                 
                                 min<-min+(npm01*(npm01+1)/2)
                                 if(npm01>0&npm02>0){
                                   V0102<- matrix(data=output[(min+1):(min+npm02*npm01)],
                                                  nrow=npm02,ncol=npm01)
                                 }
                                 
                                 min<-min+npm02*npm01
                                 
                                 if(npm01>0&npm12>0){
                                   V0112<- matrix(data=output[(min+1):(min+npm12*npm01)],
                                                  nrow=npm12,ncol=npm01)
                                 }
                                 
                                 
                                 min<-min+npm12*npm01
                                 V02<- matrix(0,npm02,npm02)
                                 V02[lower.tri(V02,diag=TRUE)] <- output[(min+1):(min+npm02*(npm02+1)/2)]
                                 
                                 
                                 min<-min+(npm02*(npm02+1)/2)
                                 
                                 if(npm02>0&npm12>0){
                                   V0212<- matrix(data=output[(min+1):(min+npm12*npm02)],
                                                  nrow=npm12,ncol=npm02)}
                                 
                                 
                                 min<-min+npm12*npm02
                                 V12<- matrix(0,npm12,npm12)
                                 V12[lower.tri(V12,diag=TRUE)] <- output[(min+1):length(output)]
                                 
                                 
                                 V<- matrix(0,npm,npm)
                                 if(npm01>0){
                                   V[1:npm01,1:npm01]<-V01
                                   if(npm02>0){V[(npm01+1):(npm01+npm02),1:npm01]<-V0102}
                                   if(npm12>0){V[(npm01+npm02+1):npm,1:npm01]<-V0112}
                                 }
                                 if(npm02>0){
                                   V[(npm01+1):(npm01+npm02),(npm01+1):(npm01+npm02)]<-V02
                                   if(npm12>0){V[(npm01+npm02+1):npm,(npm01+1):(npm01+npm02)]<-V0212}
                                 }
                                 
                                 if(npm12>0){
                                   V[(npm01+npm02+1):npm,(npm01+npm02+1):(npm)]<-V12}
                                 
                                 V<-V+t(V)
                                 diag(V)<-diag(V)/2
                                 # hessian is - second derivatives 
                                 V<--V
                                 
                                 tr <- sum(diag(V))/npm
                                 V0<-V
                                 
                                 eigen.values<-eigen(V,symmetric=T,only.values=T)$values
                                 
                                 idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                 
                                 
                                 idpos0<-idpos
                                 
                                 ncount<-0
                                 
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
                                 
                                 
                                 
                                 # attention -V as second derivative = -H
                                 output.cv<-cv.model(beta=beta,
                                                     nva01=npm01,
                                                     nva02=npm02,
                                                     nva12=npm12,
                                                     fix=fix0[7:size_V],
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
                                 
                                 betanew<-b[(6+1):size_V]
                                 
                                 res<-idmlLikelihoodweibpena(b0=b,
                                                             np0=length(b),
                                                             npar0=size_V,
                                                             bfix0=1,
                                                             fix0=rep(0,size_V),
                                                             c0=ctime,
                                                             no0=N,
                                                             ve010=ve01,
                                                             ve020=ve02,
                                                             ve120=ve12,
                                                             dimnva01=dimnva01,
                                                             dimnva02=dimnva02,
                                                             dimnva12=dimnva12,
                                                             nva01=nvat01,
                                                             nva02=nvat02,
                                                             nva12=nvat12,
                                                             t00=t0,
                                                             t10=t1,
                                                             t20=t2,
                                                             t30=t3,
                                                             troncature0=troncature,
                                                             lambda=lambda[id.lambda,],
                                                             alpha=alpha,
                                                             penalty.factor=penalty.factor,
                                                             penalty=penalty)
                                 
                                 
                                 # we want to maximise the loglik thus : 
                                 # we have issue if res is NA or if not higher than previous one 
                                 
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
                                   
                                   sears<-searpas(vw=vw,
                                                  step=step,
                                                  b=beta,
                                                  delta=delta,
                                                  funcpa=idmlLikelihoodweibpena,
                                                  res.out.error=res.out.error,
                                                  np0=length(beta),
                                                  npar0=size_V,
                                                  bfix0=s,
                                                  fix0=fix0,
                                                  c0=ctime,
                                                  no0=N,
                                                  ve010=ve01,
                                                  ve020=ve02,
                                                  ve120=ve12,
                                                  dimnva01=dimnva01,
                                                  dimnva02=dimnva02,
                                                  dimnva12=dimnva12,
                                                  nva01=nvat01,
                                                  nva02=nvat02,
                                                  nva12=nvat12,
                                                  t00=t0,
                                                  t10=t1,
                                                  t20=t2,
                                                  t30=t3,
                                                  troncature0=troncature,
                                                  lambda=lambda[id.lambda,],
                                                  alpha=alpha,
                                                  penalty.factor=penalty.factor,
                                                  penalty=penalty)
                                   
                                   
                                   betanew<-beta+delta*sears$vw
                                   b<-c(s,betanew)
                                   
                                   
                                   res<-idmlLikelihoodweibpena(b0=b,
                                                               np0=length(b),
                                                               npar0=size_V,
                                                               bfix0=1,
                                                               fix0=rep(0,size_V),
                                                               c0=ctime,
                                                               no0=N,
                                                               ve010=ve01,
                                                               ve020=ve02,
                                                               ve120=ve12,
                                                               dimnva01=dimnva01,
                                                               dimnva02=dimnva02,
                                                               dimnva12=dimnva12,
                                                               nva01=nvat01,
                                                               nva02=nvat02,
                                                               nva12=nvat12,
                                                               t00=t0,
                                                               t10=t1,
                                                               t20=t2,
                                                               t30=t3,
                                                               troncature0=troncature,
                                                               lambda=lambda[id.lambda,],
                                                               alpha=alpha,
                                                               penalty.factor=penalty.factor,
                                                               penalty=penalty)
                                 }
                                 
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
                                 
                                 output.mla<- marqLevAlg::mla(b=b,
                                                  fn=idmlLikelihoodweib,
                                                  epsa=epsa,
                                                  epsb=epsb,
                                                  epsd=epsd,
                                                  print.info = print.info,
                                                  maxiter=maxiter.pena,
                                                  minimize=F,
                                                  np0=length(b),
                                                  npar0=size_V,
                                                  bfix0=bfix,
                                                  fix0=fix0.beta,
                                                  c0=ctime,
                                                  no0=N,
                                                  ve010=ve01,
                                                  ve020=ve02,
                                                  ve120=ve12,
                                                  dimnva01=dimnva01,
                                                  dimnva02=dimnva02,
                                                  dimnva12=dimnva12,
                                                  nva01=nvat01,
                                                  nva02=nvat02,
                                                  nva12=nvat12,
                                                  t00=t0,
                                                  t10=t1,
                                                  t20=t2,
                                                  t30=t3,
                                                  troncature0=troncature)
                                 
                                 # look at convergence for each lambda :
                                
                                 # new values for splines:
                                 snew<-s
                                 snew[fix00[1:6]==0]<-output.mla$b
                                 if(nvat01>0){
                                 b01<-betanew[1:nvat01][penalty.factor[1:nvat01]==1]
                                 }else{b01<-0}
                                 if(nvat02>0){
                                 b02<-betanew[(nvat01+1):(nvat01+nvat02)][penalty.factor[(nvat01+1):(nvat01+nvat02)]==1]
                                 }else{b02<-0}
                                 if(nvat12>0){
                                 b12<-betanew[(nvat01+nvat02+1):length(betanew)][penalty.factor[(nvat01+nvat02+1):length(betanew)]==1]
                                 }else{b12<-0}
                                 if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
                                   fn.valuenew<-output.mla$fn.value+lambda[id.lambda,1]*alpha*sum(abs(b01))+lambda[id.lambda,1]*(1-alpha)*sum(b01*b01)
                                   fn.valuenew<-fn.valuenew+lambda[id.lambda,2]*alpha*sum(abs(b02))+lambda[id.lambda,2]*(1-alpha)*sum(b02*b02)
                                   fn.valuenew<-fn.valuenew+lambda[id.lambda,3]*alpha*sum(abs(b12))+lambda[id.lambda,3]*(1-alpha)*sum(b12*b12)
                                 }
                                 
                                 
                                 if(penalty=="mcp"){
                                   p01<-unlist(sapply(1:length(b01),FUN=function(x){
                                     if(b01[x]<=alpha*lambda[id.lambda,1]){
                                       return(lambda[id.lambda,1]*abs(b01[x])-((b01[x]*b01[x])/2*alpha))
                                     }else{return(alpha*lambda[id.lambda,1]*lambda[id.lambda,1]/2)}
                                   }))
                                   
                                   p02<-unlist(sapply(1:length(b02),FUN=function(x){
                                     if(b02[x]<=alpha*lambda[id.lambda,2]){
                                       return(lambda[id.lambda,2]*abs(b02[x])-((b02[x]*b02[x])/2*alpha))
                                     }else{return(alpha*lambda[id.lambda,2]*lambda[id.lambda,2]/2)}
                                   }))
                                   p12<-unlist(sapply(1:length(b12),FUN=function(x){
                                     if(b12[x]<=alpha*lambda[id.lambda,3]){
                                       return(lambda[id.lambda,3]*abs(b12[x])-((b12[x]*b12[x])/2*alpha))
                                     }else{return(alpha*lambda[id.lambda,3]*lambda[id.lambda,3]/2)}
                                   }))
                                   fn.valuenew<-output.mla$fn.value+sum(p01)+sum(p02)+sum(p12)
                                   
                                 }
                                 
                                 if(penalty=="scad"){
                                   
                                   p01<-unlist(sapply(1:length(b01),FUN=function(x){
                                     if(b01[x]<=lambda[id.lambda,1]){
                                       return(lambda[id.lambda,1]*abs(b01[x]))
                                     }else{
                                       if(abs(b01[x])<lambda[id.lambda,1]*alpha){
                                         return((2*alpha*lambda[id.lambda,1]*abs(b01[x])-b01[x]^2-lambda[id.lambda,1]^2)/(2*(alpha-1)))
                                       }else{
                                         return((lambda[id.lambda,1]^2)*(alpha+1)/2)
                                       }
                                     }
                                   }))
                                   
                                   p02<-unlist(sapply(1:length(b02),FUN=function(x){
                                     if(b02[x]<=lambda[id.lambda,2]){
                                       return(lambda[id.lambda,2]*abs(b02[x]))
                                     }else{
                                       if(abs(b02[x])<lambda[id.lambda,2]*alpha){
                                         return((2*alpha*lambda[id.lambda,2]*abs(b02[x])-b02[x]^2-lambda[id.lambda,2]^2)/(2*(alpha-1)))
                                       }else{
                                         return((lambda[id.lambda,2]^2)*(alpha+1)/2)
                                       }
                                     }
                                   }))
                                   
                                   p12<-unlist(sapply(1:length(b12),FUN=function(x){
                                     if(b12[x]<=lambda[id.lambda,3]){
                                       return(lambda[id.lambda,3]*abs(b12[x]))
                                     }else{
                                       if(abs(b12[x])<lambda[id.lambda,3]*alpha){
                                         return((2*alpha*lambda[id.lambda,3]*abs(b12[x])-b12[x]^2-lambda[id.lambda,3]^2)/(2*(alpha-1)))
                                       }else{
                                         return((lambda[id.lambda,3]^2)*(alpha+1)/2)
                                       }
                                     }
                                   }))
                                   fn.valuenew<-output.mla$fn.value+sum(p01)+sum(p02)+sum(p12)
                                   
                                 }
                                 
                                   
                                 ite<-ite+1
                                 
                                 eval.cv.spline[ite]<-sum((snew-s)^2)
                                 eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                 eval.cv.loglik[ite]<-abs(fn.valuenew-fn.value)
                                 eval.loglik[ite]<-fn.valuenew
                                 eval.validity[ite]<-validity
                                 
                                 s<-snew
                                 beta<-betanew
                                 fn.value<-fn.valuenew
                                 
                                 # eval.cv beta valid only if validity.param=T
                                 if(eval.cv.beta[ite]<epsa & eval.cv.spline[ite]<epsa & eval.cv.loglik[ite]<epsb & validity==T){
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
                                   
                                   output<-derivaweib(b0=b,
                                                      np0=length(b),
                                                      npar0=size_V,
                                                      bfix0=bfix,
                                                      fix0=fix0,
                                                      c0=ctime,
                                                      no0=N,
                                                      ve010=ve01,
                                                      ve020=ve02,
                                                      ve120=ve12,
                                                      dimnva01=dimnva01,
                                                      dimnva02=dimnva02,
                                                      dimnva12=dimnva12,
                                                      nva01=nvat01,
                                                      nva02=nvat02,
                                                      nva12=nvat12,
                                                      t00=t0,
                                                      t10=t1,
                                                      t20=t2,
                                                      t30=t3,
                                                      troncature0=troncature)
                                   min<-npm
                                   V01<- matrix(0,npm01,npm01)
                                   V01[lower.tri(V01,diag=TRUE)] <- output[(min+1):(min+npm01*(npm01+1)/2)]
                                   
                                   
                                   min<-min+(npm01*(npm01+1)/2)
                                   if(npm01>0&npm02>0){
                                     V0102<- matrix(data=output[(min+1):(min+npm02*npm01)],
                                                    nrow=npm02,ncol=npm01)
                                   }
                                   
                                   min<-min+npm02*npm01
                                   
                                   if(npm01>0&npm12>0){
                                     V0112<- matrix(data=output[(min+1):(min+npm12*npm01)],
                                                    nrow=npm12,ncol=npm01)
                                   }
                                   
                                   
                                   min<-min+npm12*npm01
                                   V02<- matrix(0,npm02,npm02)
                                   V02[lower.tri(V02,diag=TRUE)] <- output[(min+1):(min+npm02*(npm02+1)/2)]
                                   
                                   
                                   min<-min+(npm02*(npm02+1)/2)
                                   
                                   if(npm02>0&npm12>0){
                                     V0212<- matrix(data=output[(min+1):(min+npm12*npm02)],
                                                    nrow=npm12,ncol=npm02)}
                                   
                                   
                                   min<-min+npm12*npm02
                                   V12<- matrix(0,npm12,npm12)
                                   V12[lower.tri(V12,diag=TRUE)] <- output[(min+1):length(output)]
                                   
                                   
                                   V<- matrix(0,npm,npm)
                                   if(npm01>0){
                                     V[1:npm01,1:npm01]<-V01
                                     if(npm02>0){V[(npm01+1):(npm01+npm02),1:npm01]<-V0102}
                                     if(npm12>0){V[(npm01+npm02+1):npm,1:npm01]<-V0112}
                                   }
                                   if(npm02>0){
                                     V[(npm01+1):(npm01+npm02),(npm01+1):(npm01+npm02)]<-V02
                                     if(npm12>0){V[(npm01+npm02+1):npm,(npm01+1):(npm01+npm02)]<-V0212}
                                   }
                                   
                                   if(npm12>0){
                                     V[(npm01+npm02+1):npm,(npm01+npm02+1):(npm)]<-V12}
                                   
                                   V<-V+t(V)
                                   diag(V)<-diag(V)/2
                                   # hessian is - second derivatives 
                                   V<--V
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
    parallel::stopCluster(clustpar)
    
  }else{
    

    output<-foreach::foreach(id.lambda=1:nlambda,
                             .combine = combine_lambda,
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
                               
                               npm<-sum(fix0==0)
                               
                               npm01<-ifelse(nvat01>0,sum(fix0[7:(7+nvat01-1)]==0),0)
                               npm02<-ifelse(nvat02>0,sum(fix0[(7+nvat01):(6+nvat01+nvat02)]==0),0)
                               npm12<-ifelse(nvat12>0,sum(fix0[(7+nvat01+nvat02):size_V]==0),0)
                               
                               
                               while(converged==F & ite<=maxiter){
                               
                                 b<-c(s,beta)
                                 bfix<-b[fix0==1]
                                 b<-b[fix0==0]
                                 
                                 output<-derivaweib(b0=b,
                                                    np0=length(b),
                                                    npar0=size_V,
                                                    bfix0=bfix,
                                                    fix0=fix0,
                                                    c0=ctime,
                                                    no0=N,
                                                    ve010=ve01,
                                                    ve020=ve02,
                                                    ve120=ve12,
                                                    dimnva01=dimnva01,
                                                    dimnva02=dimnva02,
                                                    dimnva12=dimnva12,
                                                    nva01=nvat01,
                                                    nva02=nvat02,
                                                    nva12=nvat12,
                                                    t00=t0,
                                                    t10=t1,
                                                    t20=t2,
                                                    t30=t3,
                                                    troncature0=troncature)
                                 
                               if(ite==0){
                                   fn.value<-idmlLikelihoodweibpena(b=b,
                                                                    np0=length(b),
                                                                    npar0=size_V,
                                                                    bfix0=bfix,
                                                                    fix0=fix0,
                                                                    c0=ctime,
                                                                    no0=N,
                                                                    ve010=ve01,
                                                                    ve020=ve02,
                                                                    ve120=ve12,
                                                                    dimnva01=dimnva01,
                                                                    dimnva02=dimnva02,
                                                                    dimnva12=dimnva12,
                                                                    nva01=nvat01,
                                                                    nva02=nvat02,
                                                                    nva12=nvat12,
                                                                    t00=t0,
                                                                    t10=t1,
                                                                    t20=t2,
                                                                    t30=t3,
                                                                    troncature0=troncature,
                                                                    lambda=lambda[id.lambda,],
                                                                    alpha=alpha,
                                                                    penalty.factor=penalty.factor,
                                                                    penalty=penalty)
                               }
                                 
                                 if(any(is.na(output))|any(output==Inf) |any(output==-Inf)){
                                   warning("Computational error for calculation of the hessian : division by 0 or Infinite value")
                                   if(ite==0){
                                     fu <- output[1:npm]
                                     min<-npm
                                     V01<- matrix(0,npm01,npm01)
                                     V01[lower.tri(V01,diag=TRUE)] <- output[(min+1):(min+npm01*(npm01+1)/2)]
                                     
                                     
                                     min<-min+(npm01*(npm01+1)/2)
                                     if(npm01>0&npm02>0){
                                       V0102<- matrix(data=output[(min+1):(min+npm02*npm01)],
                                                      nrow=npm02,ncol=npm01)
                                     }
                                     
                                     min<-min+npm02*npm01
                                     
                                     if(npm01>0&npm12>0){
                                       V0112<- matrix(data=output[(min+1):(min+npm12*npm01)],
                                                      nrow=npm12,ncol=npm01)
                                     }
                                     
                                     
                                     min<-min+npm12*npm01
                                     V02<- matrix(0,npm02,npm02)
                                     V02[lower.tri(V02,diag=TRUE)] <- output[(min+1):(min+npm02*(npm02+1)/2)]
                                     
                                     
                                     min<-min+(npm02*(npm02+1)/2)
                                     
                                     if(npm02>0&npm12>0){
                                       V0212<- matrix(data=output[(min+1):(min+npm12*npm02)],
                                                      nrow=npm12,ncol=npm02)}
                                     
                                     
                                     min<-min+npm12*npm02
                                     V12<- matrix(0,npm12,npm12)
                                     V12[lower.tri(V12,diag=TRUE)] <- output[(min+1):length(output)]
                                     
                                     
                                     V<- matrix(0,npm,npm)
                                     if(npm01>0){
                                       V[1:npm01,1:npm01]<-V01
                                       if(npm02>0){V[(npm01+1):(npm01+npm02),1:npm01]<-V0102}
                                       if(npm12>0){V[(npm01+npm02+1):npm,1:npm01]<-V0112}
                                     }
                                     if(npm02>0){
                                       V[(npm01+1):(npm01+npm02),(npm01+1):(npm01+npm02)]<-V02
                                       if(npm12>0){V[(npm01+npm02+1):npm,(npm01+1):(npm01+npm02)]<-V0212}
                                     }
                                     
                                     if(npm12>0){
                                       V[(npm01+npm02+1):npm,(npm01+npm02+1):(npm)]<-V12}
                                     
                                     V<-V+t(V)
                                     diag(V)<-diag(V)/2
                                     # hessian is - second derivatives 
                                     V<--V
                                     
                                     tr <- sum(diag(V))/npm
                                     V0<-V}
                                   ite<-ite+1
                                   pbr_compu<-1
                                   break
                                 }
                                 
                                 fu <- output[1:npm]
                                 
                                 min<-npm
                                 V01<- matrix(0,npm01,npm01)
                                 V01[lower.tri(V01,diag=TRUE)] <- output[(min+1):(min+npm01*(npm01+1)/2)]
                                 
                                 
                                 min<-min+(npm01*(npm01+1)/2)
                                 if(npm01>0&npm02>0){
                                 V0102<- matrix(data=output[(min+1):(min+npm02*npm01)],
                                                nrow=npm02,ncol=npm01)
                                 }
                                 
                                 min<-min+npm02*npm01
                                 
                                 if(npm01>0&npm12>0){
                                 V0112<- matrix(data=output[(min+1):(min+npm12*npm01)],
                                                nrow=npm12,ncol=npm01)
                                 }
                                 
                                 
                                 min<-min+npm12*npm01
                                 V02<- matrix(0,npm02,npm02)
                                 V02[lower.tri(V02,diag=TRUE)] <- output[(min+1):(min+npm02*(npm02+1)/2)]
                                 
                                 
                                 min<-min+(npm02*(npm02+1)/2)
                                 
                                 if(npm02>0&npm12>0){
                                 V0212<- matrix(data=output[(min+1):(min+npm12*npm02)],
                                                nrow=npm12,ncol=npm02)}
                                 
                                 
                                 min<-min+npm12*npm02
                                 V12<- matrix(0,npm12,npm12)
                                 V12[lower.tri(V12,diag=TRUE)] <- output[(min+1):length(output)]
                                 
                                 
                                 V<- matrix(0,npm,npm)
                                 if(npm01>0){
                                   V[1:npm01,1:npm01]<-V01
                                   if(npm02>0){V[(npm01+1):(npm01+npm02),1:npm01]<-V0102}
                                   if(npm12>0){V[(npm01+npm02+1):npm,1:npm01]<-V0112}
                                 }
                                 if(npm02>0){
                                 V[(npm01+1):(npm01+npm02),(npm01+1):(npm01+npm02)]<-V02
                                 if(npm12>0){V[(npm01+npm02+1):npm,(npm01+1):(npm01+npm02)]<-V0212}
                                 }
                                 
                                 if(npm12>0){
                                 V[(npm01+npm02+1):npm,(npm01+npm02+1):(npm)]<-V12}
                                 
                                 V<-V+t(V)
                                 diag(V)<-diag(V)/2
                                # hessian is - second derivatives 
                                 V<--V
                                 
                                 tr <- sum(diag(V))/npm
                                 V0<-V
                                 
                                 eigen.values<-eigen(V,symmetric=T,only.values=T)$values
                                
                                 idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                 
                                 
                                 idpos0<-idpos
                                 
                                 ncount<-0
                                 
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
                                 
                                 # attention -V as second derivative = -H
                                 output.cv<-cv.model(beta=beta,
                                                    nva01=npm01,
                                                    nva02=npm02,
                                                    nva12=npm12,
                                                    fix=fix0[7:size_V],
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
                                 
                                 betanew<-b[(6+1):size_V]
                                 
                                 res<-idmlLikelihoodweibpena(b0=b,
                                                         np0=length(b),
                                                         npar0=size_V,
                                                         bfix0=1,
                                                         fix0=rep(0,size_V),
                                                         c0=ctime,
                                                         no0=N,
                                                         ve010=ve01,
                                                         ve020=ve02,
                                                         ve120=ve12,
                                                         dimnva01=dimnva01,
                                                         dimnva02=dimnva02,
                                                         dimnva12=dimnva12,
                                                         nva01=nvat01,
                                                         nva02=nvat02,
                                                         nva12=nvat12,
                                                         t00=t0,
                                                         t10=t1,
                                                         t20=t2,
                                                         t30=t3,
                                                         troncature0=troncature,
                                                         lambda=lambda[id.lambda,],
                                                         alpha=alpha,
                                                         penalty.factor=penalty.factor,
                                                         penalty=penalty)
                                 
                                 
                                 # we want to maximise the loglik thus : 
                                 # we have issue if res is NA or if not higher than previous one 
                                 
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
                                   
                                   sears<-searpas(vw=vw,
                                                   step=step,
                                                   b=beta,
                                                   delta=delta,
                                                   funcpa=idmlLikelihoodweibpena,
                                                   res.out.error=res.out.error,
                                                   np0=length(beta),
                                                   npar0=size_V,
                                                   bfix0=s,
                                                   fix0=fix0,
                                                   c0=ctime,
                                                   no0=N,
                                                   ve010=ve01,
                                                   ve020=ve02,
                                                   ve120=ve12,
                                                   dimnva01=dimnva01,
                                                   dimnva02=dimnva02,
                                                   dimnva12=dimnva12,
                                                   nva01=nvat01,
                                                   nva02=nvat02,
                                                   nva12=nvat12,
                                                   t00=t0,
                                                   t10=t1,
                                                   t20=t2,
                                                   t30=t3,
                                                   troncature0=troncature,
                                                   lambda=lambda[id.lambda,],
                                                   alpha=alpha,
                                                  penalty.factor=penalty.factor,
                                                  penalty=penalty)
                                   
                                   
                                   betanew<-beta+delta*sears$vw
                                   b<-c(s,betanew)
                                   
                                   
                                   res<-idmlLikelihoodweibpena(b0=b,
                                                               np0=length(b),
                                                               npar0=size_V,
                                                               bfix0=1,
                                                               fix0=rep(0,size_V),
                                                               c0=ctime,
                                                               no0=N,
                                                               ve010=ve01,
                                                               ve020=ve02,
                                                               ve120=ve12,
                                                               dimnva01=dimnva01,
                                                               dimnva02=dimnva02,
                                                               dimnva12=dimnva12,
                                                               nva01=nvat01,
                                                               nva02=nvat02,
                                                               nva12=nvat12,
                                                               t00=t0,
                                                               t10=t1,
                                                               t20=t2,
                                                               t30=t3,
                                                               troncature0=troncature,
                                                               lambda=lambda[id.lambda,],
                                                               alpha=alpha,
                                                               penalty.factor=penalty.factor,
                                                               penalty=penalty)
                                 }
                                   
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
                                 
                                 output.mla<- marqLevAlg::mla(b=b,
                                                  fn=idmlLikelihoodweib,
                                                  epsa=epsa,
                                                  epsb=epsb,
                                                  epsd=epsd,
                                                  print.info = print.info,
                                                  maxiter=maxiter.pena,
                                                  minimize=F,
                                                  np0=length(b),
                                                  npar0=size_V,
                                                  bfix0=bfix,
                                                  fix0=fix0.beta,
                                                  c0=ctime,
                                                  no0=N,
                                                  ve010=ve01,
                                                  ve020=ve02,
                                                  ve120=ve12,
                                                  dimnva01=dimnva01,
                                                  dimnva02=dimnva02,
                                                  dimnva12=dimnva12,
                                                  nva01=nvat01,
                                                  nva02=nvat02,
                                                  nva12=nvat12,
                                                  t00=t0,
                                                  t10=t1,
                                                  t20=t2,
                                                  t30=t3,
                                                  troncature0=troncature)
                                
                                 # look at convergence for each lambda :
                                 
                                 # new values for splines:
                                 snew<-s
                                 snew[fix00[1:6]==0]<-output.mla$b
                                 if(nvat01>0){
                                 b01<-betanew[1:nvat01][penalty.factor[1:nvat01]==1]
                                 }else{b01<-0}
                                 if(nvat02>0){
                                 b02<-betanew[(nvat01+1):(nvat01+nvat02)][penalty.factor[(nvat01+1):(nvat01+nvat02)]==1]
                                 }else{b02<-0}
                                 if(nvat12>0){
                                 b12<-betanew[(nvat01+nvat02+1):length(betanew)][penalty.factor[(nvat01+nvat02+1):length(betanew)]==1]
                                 }else{b12<-0}
                                 
                                 if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
                                   fn.valuenew<-output.mla$fn.value+lambda[id.lambda,1]*alpha*sum(abs(b01))+lambda[id.lambda,1]*(1-alpha)*sum(b01*b01)
                                   fn.valuenew<-fn.valuenew+lambda[id.lambda,2]*alpha*sum(abs(b02))+lambda[id.lambda,2]*(1-alpha)*sum(b02*b02)
                                   fn.valuenew<-fn.valuenew+lambda[id.lambda,3]*alpha*sum(abs(b12))+lambda[id.lambda,3]*(1-alpha)*sum(b12*b12)
                                 }
                                 
                                 if(penalty=="mcp"){
                                   p01<-unlist(sapply(1:length(b01),FUN=function(x){
                                     if(b01[x]<=alpha*lambda[id.lambda,1]){
                                       return(lambda[id.lambda,1]*abs(b01[x])-((b01[x]*b01[x])/2*alpha))
                                     }else{return(alpha*lambda[id.lambda,1]*lambda[id.lambda,1]/2)}
                                   }))
                                   
                                   p02<-unlist(sapply(1:length(b02),FUN=function(x){
                                     if(b02[x]<=alpha*lambda[id.lambda,2]){
                                       return(lambda[id.lambda,2]*abs(b02[x])-((b02[x]*b02[x])/2*alpha))
                                     }else{return(alpha*lambda[id.lambda,2]*lambda[id.lambda,2]/2)}
                                   }))
                                   p12<-unlist(sapply(1:length(b12),FUN=function(x){
                                     if(b12[x]<=alpha*lambda[id.lambda,3]){
                                       return(lambda[id.lambda,3]*abs(b12[x])-((b12[x]*b12[x])/2*alpha))
                                     }else{return(alpha*lambda[id.lambda,3]*lambda[id.lambda,3]/2)}
                                   }))
                                   fn.valuenew<-output.mla$fn.value+sum(p01)+sum(p02)+sum(p12)
                                   
                                 }
                                 
                                 if(penalty=="scad"){
                                   
                                   p01<-unlist(sapply(1:length(b01),FUN=function(x){
                                     if(b01[x]<=lambda[id.lambda,1]){
                                       return(lambda[id.lambda,1]*abs(b01[x]))
                                     }else{
                                       if(abs(b01[x])<lambda[id.lambda,1]*alpha){
                                         return((2*alpha*lambda[id.lambda,1]*abs(b01[x])-b01[x]^2-lambda[id.lambda,1]^2)/(2*(alpha-1)))
                                       }else{
                                         return((lambda[id.lambda,1]^2)*(alpha+1)/2)
                                       }
                                     }
                                   }))
                                   
                                   p02<-unlist(sapply(1:length(b02),FUN=function(x){
                                     if(b02[x]<=lambda[id.lambda,2]){
                                       return(lambda[id.lambda,2]*abs(b02[x]))
                                     }else{
                                       if(abs(b02[x])<lambda[id.lambda,2]*alpha){
                                         return((2*alpha*lambda[id.lambda,2]*abs(b02[x])-b02[x]^2-lambda[id.lambda,2]^2)/(2*(alpha-1)))
                                       }else{
                                         return((lambda[id.lambda,2]^2)*(alpha+1)/2)
                                       }
                                     }
                                   }))
                                   
                                   p12<-unlist(sapply(1:length(b12),FUN=function(x){
                                     if(b12[x]<=lambda[id.lambda,3]){
                                       return(lambda[id.lambda,3]*abs(b12[x]))
                                     }else{
                                       if(abs(b12[x])<lambda[id.lambda,3]*alpha){
                                         return((2*alpha*lambda[id.lambda,3]*abs(b12[x])-b12[x]^2-lambda[id.lambda,3]^2)/(2*(alpha-1)))
                                       }else{
                                         return((lambda[id.lambda,3]^2)*(alpha+1)/2)
                                       }
                                     }
                                   }))
                                   fn.valuenew<-output.mla$fn.value+sum(p01)+sum(p02)+sum(p12)
                                   
                                 }
                                 
                                   
                                 ite<-ite+1
                                 
                                 eval.cv.spline[ite]<-sum((snew-s)^2)
                                 eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                 eval.cv.loglik[ite]<-abs(fn.valuenew-fn.value)
                                 eval.loglik[ite]<-fn.valuenew
                                 eval.validity[ite]<-validity
                                 
                                 s<-snew
                                 beta<-betanew
                                 fn.value<-fn.valuenew
                                 
                                 # eval.cv beta valid only if validity.param=T
                                 if(eval.cv.beta[ite]<epsa & eval.cv.spline[ite]<epsa & eval.cv.loglik[ite]<epsb & validity==T){
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
                                   
                                   output<-derivaweib(b0=b,
                                                      np0=length(b),
                                                      npar0=size_V,
                                                      bfix0=bfix,
                                                      fix0=fix0,
                                                      c0=ctime,
                                                      no0=N,
                                                      ve010=ve01,
                                                      ve020=ve02,
                                                      ve120=ve12,
                                                      dimnva01=dimnva01,
                                                      dimnva02=dimnva02,
                                                      dimnva12=dimnva12,
                                                      nva01=nvat01,
                                                      nva02=nvat02,
                                                      nva12=nvat12,
                                                      t00=t0,
                                                      t10=t1,
                                                      t20=t2,
                                                      t30=t3,
                                                      troncature0=troncature)
                                   min<-npm
                                   V01<- matrix(0,npm01,npm01)
                                   V01[lower.tri(V01,diag=TRUE)] <- output[(min+1):(min+npm01*(npm01+1)/2)]
                                   
                                   
                                   min<-min+(npm01*(npm01+1)/2)
                                   if(npm01>0&npm02>0){
                                     V0102<- matrix(data=output[(min+1):(min+npm02*npm01)],
                                                    nrow=npm02,ncol=npm01)
                                   }
                                   
                                   min<-min+npm02*npm01
                                   
                                   if(npm01>0&npm12>0){
                                     V0112<- matrix(data=output[(min+1):(min+npm12*npm01)],
                                                    nrow=npm12,ncol=npm01)
                                   }
                                   
                                   
                                   min<-min+npm12*npm01
                                   V02<- matrix(0,npm02,npm02)
                                   V02[lower.tri(V02,diag=TRUE)] <- output[(min+1):(min+npm02*(npm02+1)/2)]
                                   
                                   
                                   min<-min+(npm02*(npm02+1)/2)
                                   
                                   if(npm02>0&npm12>0){
                                     V0212<- matrix(data=output[(min+1):(min+npm12*npm02)],
                                                    nrow=npm12,ncol=npm02)}
                                   
                                   
                                   min<-min+npm12*npm02
                                   V12<- matrix(0,npm12,npm12)
                                   V12[lower.tri(V12,diag=TRUE)] <- output[(min+1):length(output)]
                                   
                                   
                                   V<- matrix(0,npm,npm)
                                   if(npm01>0){
                                     V[1:npm01,1:npm01]<-V01
                                     if(npm02>0){V[(npm01+1):(npm01+npm02),1:npm01]<-V0102}
                                     if(npm12>0){V[(npm01+npm02+1):npm,1:npm01]<-V0112}
                                   }
                                   if(npm02>0){
                                     V[(npm01+1):(npm01+npm02),(npm01+1):(npm01+npm02)]<-V02
                                     if(npm12>0){V[(npm01+npm02+1):npm,(npm01+1):(npm01+npm02)]<-V0212}
                                   }
                                   
                                   if(npm12>0){
                                     V[(npm01+npm02+1):npm,(npm01+npm02+1):(npm)]<-V12}
                                   
                                   V<-V+t(V)
                                   diag(V)<-diag(V)/2
                                   # hessian is - second derivatives 
                                   V<--V
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
  }
  
  
  
  
  
  
  return(output=output)
}
