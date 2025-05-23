causal.idm.weib<-function(b,fix0,size_V,
                   clustertype,epsa,epsb,epsd,eps.eigen,nproc,maxiter,
                   ctime,N,
                   ve01,ve02,ve12,dimnva01,dimnva02,dimnva12,nvat01,nvat02,nvat12,nvat12dep,
                   semiMark,t0,t1,t2,t3,idd,idm,ts,troncature){

  # initiate parameters values 
  if(!is.null(b)){
  
   s.start<-b[1:6]
  
  if(size_V>6){
  beta.start<-b[(6+1):(size_V)]
  }else{beta.start<-NULL}}
  

  # if no intial values for beta do some iterations with mla 
  if(is.null(b)){
      
      s.start<-c(1,sqrt(sum(idm)/ts),1,sqrt(sum(idd)/ts),1,sqrt(sum(idd)/ts))
      beta.start<-rep(0,size_V-6)
      # initialise parameters of weib without any variables of interest
      start.weib<-c(s.start)
      bfix.weib<-c(start.weib[fix0[1:6]==1],beta.start)
      start.weib<-start.weib[fix0[1:6]==0]
      fix.weib<-fix0
      if(size_V>6){
        fix.weib[(6+1):size_V]<-1}

      output.mla<- marqLevAlg::mla(b=start.weib,
                       fn=causalidmlLikelihoodweib,
                       epsa=epsa,
                       epsb=epsb,
                       epsd=epsd,
                       nproc=nproc,
                       clustertype = clustertype,
                       maxiter=maxiter,
                       minimize=F,
                       npm=length(start.weib),
                       npar=size_V,
                       bfix=bfix.weib,
                       fix=fix.weib,
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
                       nva12dep=nvat12dep,
                       semiMark=semiMark,
                       t0=t0,
                       t1=t1,
                       t2=t2,
                       t3=t3,
                       troncature=troncature)

      if(output.mla$istop==1){
        s.start<-output.mla$b}
    }                           
  beta<-beta.start
  s<-s.start
  b<-c(s,beta)
  bfix<-b[fix0==1]
  b<-b[fix0==0]
  
  # maximise loglik 
  out<- marqLevAlg::mla(b=b,
                    fn=causalidmlLikelihoodweib,
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
                    nva12dep=nvat12dep,
                    semiMark=semiMark,
                    t0=t0,
                    t1=t1,
                    t2=t2,
                    t3=t3,
                    troncature=troncature)

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
   
   
  
  
     
