### Code:
##' @title Illness-death model algorithm with weibull baseline risk
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
##' @param clustertype in which cluster to work
##' @param nproc number of cluster
##' @param maxiter Maximum number of iterations. The default is 200.
##' @param troncature indicator if troncature or not
##' @param idd number of subjects that died
##' @param idm number of subjects that had illness
##' @param ts delay in the study
##' @param gausspoint number of points in gauss quadrature
##' @param methodCV methodolgy to optimise penalised parameters 
##' @param weib the form of the weibull parameters 
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
#' @useDynLib SmoothHazardoptim9

idm.weib<-function(b,fix0,size_V,
                   clustertype,epsa,epsb,epsd,nproc,maxiter,
                   ctime,N,
                   ve01,ve02,ve12,dimnva01,dimnva02,dimnva12,nvat01,nvat02,nvat12,
                   t0,t1,t2,t3,idd,idm,ts,troncature,gausspoint,weib,methodCV){

  bfix<-b[fix0==1]
  b<-b[fix0==0]

  if(methodCV%in%c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                 "Brent")){
    

  out<- optim(par=b,
              method=methodCV,
              hessian = T,
              control=list(maxit=maxiter,
                           fnscale=-1, #to maximise function
                           abstol=epsb,
                           parscale=rep(1,length(b))),
              
                        fn=idmlLikelihoodweib,
                        npm=length(b),
                        npar=size_V,
                        bfix=bfix,
                        fix=fix0,
                        ctime=ctime,
                        no=N,
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
                        gausspoint=gausspoint,weib=weib)
  V<-solve(-out$hessian)
  
  return(list(b=out$par,
              fn.value=out$value,
              ni=out$counts[1],
              istop=out$convergence,
              v=V,
              grad=NULL,
              ca=NULL,
              cb=NULL,
              rdm=NULL,
              bfix=bfix,
              fix0=fix0))
  }else{

  # maximise loglik 
    #browser()
  out<- marqLevAlg::mla(b=b,
                        print.info=T,
                    fn=idmlLikelihoodweib,
                    gr=grmlaweib,
                    hess = hessianmlaweib,
                    epsa=epsa,
                    epsb=epsb,
                    epsd=epsd,
                    nproc=nproc,
                    clustertype=clustertype,
                    maxiter=maxiter,
                    minimize=F,
                    npm=length(b),
                    npar=size_V,
                    bfix=bfix,
                    fix=fix0,
                    ctime=ctime,
                    no=N,
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
                    gausspoint=gausspoint,weib=weib)
  
  
  if(out$istop==4){
    stop("Problem in the loglikelihood computation.")
  }
  return(list(b=out$b,
              fn.value=out$fn.value,
              ni=out$ni,
              istop=out$istop,
              v=out$v,
              grad=out$grad,
              ca=out$ca,
              cb=out$cb,
              rdm=out$rdm,
              bfix=bfix,
              fix0=fix0))
  }
  


  
}
   
   
  
  
     
