idm.no.step.sequential<-function(b,clustertype,epsa,epsb,epsd,print.info,nproc,
                                 maxiter,size_V,bfix,fix0,knots01,knots02,knots12,
                                 ctime,N,nknots01,nknots02,nknots12,ve01,ve02,ve12,
                                 dimnva01,dimnva02,dimnva12,nvat01,nvat02,nvat12,t0,t1,
                                 t2,t3,troncature,gauss.point,penalty,lambda,gamma){

  posfix<-which(fix0==1)
  npm<-length(b)
  
  
  if(length(lambda)==1){
    
    out<-marqLevAlg::mla(b=b,
             m=npm,
             fn=idmlLikelihoodpena,
             clustertype=clustertype,
             epsa=epsa,
             epsb=epsb,
             epsd=epsd,
             print.info = print.info,
             nproc=nproc,
             maxiter=maxiter,
             minimize=F,
             np0=npm,
             npar0=size_V,
             bfix0=bfix,
             fix0=fix0,
             zi010=knots01,
             zi020=knots02,
             zi120=knots12,
             c0=ctime,
             no0=N,
             nz01=nknots01,
             nz02=nknots02,
             nz12=nknots12,
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
             gausspoint0=gauss.point,
             penalty=penalty,
             lambda=lambda,
             gamma=gamma
    )
    
    
    beta<-rep(NA,size_V)
    beta[which(fix0==0)]<-out$b
    beta[posfix]<-bfix
    

    Vr <- matrix(0,npm,npm)
    Vr[upper.tri(Vr,diag=TRUE)] <-out$v[1:(npm*(npm+1)/2)]
    Vr <- t(Vr)
    Vr[upper.tri(Vr,diag=TRUE)] <- out$v[1:(npm*(npm+1)/2)]
    V <- matrix(0,size_V,size_V)
    V[setdiff(1:size_V,posfix),setdiff(1:size_V,posfix)] <- Vr
    
    
    if(out$istop==1){
      
      Hr<--solve(Vr)
      H<- matrix(0,size_V,size_V)
      H[setdiff(1:size_V,posfix),setdiff(1:size_V,posfix)] <- Hr
      
      
     
      deriv <- -deriva(1,out$b,funcpa=idmlLikelihood,
                      np0=npm,
                      npar0=size_V,
                      bfix0=bfix,
                      fix0=fix0,
                      zi010=knots01,
                      zi020=knots02,
                      zi120=knots12,
                      c0=ctime,
                      no0=N,
                      nz01=nknots01,
                      nz02=nknots02,
                      nz12=nknots12,
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
                      gausspoint0=gauss.point)$v
      
      Hnopena <- matrix(0,npm,npm)
      Hnopena[upper.tri(Hnopena,diag=TRUE)] <-deriv[1:(npm*(npm+1)/2)]
      Hnopena <- t(Hnopena)
      Hnopena[upper.tri(Hnopena,diag=TRUE)] <- deriv[1:(npm*(npm+1)/2)]
      H.nopena <- matrix(0,size_V,size_V)
      H.nopena[setdiff(1:size_V,posfix),setdiff(1:size_V,posfix)] <- Hnopena

      
      theta<-beta
      size_spline<-nknots01+nknots02+nknots12+6
      theta[size_spline:length(theta)]<-rep(0,size_spline+length(theta)-1)
      theta<-theta[which(fix0==0)]
      browser()
      V.theta<-4*theta*Vr%*%diag(theta)
      H.nopena.theta<-4*theta*Hnopena%*%diag(theta)
      browser()
      trLCV.theta<-lava::tr(-V.theta%*%H.nopena.theta)
      trLCV<-lava::tr(-V%*%H.nopena)
      
      
    }else{
      H<-matrix(0,size_V,size_V)
      trLCV<-0
      trLCV.theta

    }
    
    
    
    return(list(b=as.matrix(beta),
                lambda=lambda,
                gamma=gamma,
                v=V,
                H=H,
                trLCV=trLCV,
                trLCV.theta=trLCV.theta,
                fn.value=out$fn.value,
                ni=out$ni,
                istop=out$istop,
                ca=out$ca,
                cb=out$cb,
                rdm=out$rdm,
                fix0=fix0))
    
    
  }else{
    
    if(nproc==1){
      warning("When having multiple lambda, nproc must be at least 2, by default it will be put at 2")
      nproc<-2
    }
    
    if(is.null(clustertype)){
        clustpar <- parallel::makeCluster(nproc)#, outfile="")
      }
      else{
        clustpar <- parallel::makeCluster(nproc, type=clustertype)#, outfile="")
      }
      
    doParallel::registerDoParallel(clustpar)
    
    
    i<-NULL # for cran check 
    res<-foreach::foreach(i = 1:length(lambda),.combine=cbind) %dopar% {
      
      
      out<-marqLevAlg::mla(b=b,
               m=npm,
               fn=idmlLikelihoodpena,
               clustertype=clustertype,
               epsa=epsa,
               epsb=epsb,
               epsd=epsd,
               print.info = print.info,
               nproc=1,
               maxiter=maxiter,
               minimize=F,
               np0=npm,
               npar0=size_V,
               bfix0=bfix,
               fix0=fix0,
               zi010=knots01,
               zi020=knots02,
               zi120=knots12,
               c0=ctime,
               no0=N,
               nz01=nknots01,
               nz02=nknots02,
               nz12=nknots12,
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
               gausspoint0=gauss.point,
               penalty=penalty,
               lambda=lambda[i],
               gamma=gamma
      )
      
      
      c(lambda[i],out$fn.value,
        out$ni,
        out$istop,
        out$ca,
        out$cb,
        out$rdm,
        out$b,
        out$v)
      
    }
    
    parallel::stopCluster(clustpar)
    
    lambda<-res[1,]
    fn.value<-res[2,]
    ni<-res[3,]
    istop<-res[4,]
    ca<-res[5,]
    cb<-res[6,]
    rdm<-res[7,]
    
    
    beta<-matrix(NA,nrow=size_V,ncol=length(lambda))
    beta[which(fix0==0),]<-res[8:(7+npm),]
    beta[posfix,]<-bfix
    
    
    v<-matrix(NA,nrow=(npm*(npm+1)/2),ncol=length(lambda))
    v<-res[(8+npm):(7+npm+(npm*(npm+1)/2)),]
    
  
    combine_var<-function(x,newx){
      list(V=list(x[[1]],newx[[1]]),
           H=list(x[[2]],newx[[2]]),
           trLCV=list(x[[3]],newx[[3]]))
    }
    
    if(is.null(clustertype)){
      clustpar <- parallel::makeCluster(nproc)#, outfile="")
    }
    else{
      clustpar <- parallel::makeCluster(nproc, type=clustertype)#, outfile="")
    }
    
    doParallel::registerDoParallel(clustpar)
    
    
    i<-NULL # for cran check 
    res<-foreach::foreach(i=1:length(lambda),.combine = combine_var)%dopar%{
      
      Vr <- matrix(0,npm,npm)
      Vr[upper.tri(Vr,diag=TRUE)] <-v[,i]
      Vr <- t(Vr)
      Vr[upper.tri(Vr,diag=TRUE)] <- v[,i]
      V <- matrix(0,size_V,size_V)
      V[setdiff(1:size_V,posfix),setdiff(1:size_V,posfix)] <- Vr
      
      
      if(istop[i]==1){
        
        Hr<--solve(Vr)
        H<- matrix(0,size_V,size_V)
        H[setdiff(1:size_V,posfix),setdiff(1:size_V,posfix)] <- Hr
        
        x<-as.vector(beta[,i])[which(fix0==0)]

     
        
        deriv <- -deriva(1,x,funcpa=idmlLikelihood,
                         np0=length(x),
                         npar0=size_V,
                         bfix0=bfix,
                         fix0=fix0,
                         zi010=knots01,
                         zi020=knots02,
                         zi120=knots12,
                         c0=ctime,
                         no0=N,
                         nz01=nknots01,
                         nz02=nknots02,
                         nz12=nknots12,
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
                         gausspoint0=gauss.point)$v
        
        Hnopena <- matrix(0,npm,npm)
        Hnopena[upper.tri(Hnopena,diag=TRUE)] <-deriv[1:(npm*(npm+1)/2)]
        Hnopena <- t(Hnopena)
        Hnopena[upper.tri(Hnopena,diag=TRUE)] <- deriv[1:(npm*(npm+1)/2)]
        H.nopena <- matrix(0,size_V,size_V)
        H.nopena[setdiff(1:size_V,posfix),setdiff(1:size_V,posfix)] <- Hnopena
        
        theta<-beta[,i]
        size_spline<-nknots01+nknots02+nknots12+6
        theta[size_spline:length(theta)]<-rep(0,size_spline+length(theta)-1)
        theta<-theta[which(fix0==0)]
        browser()
        V.theta<-4*theta*Vr%*%diag(theta)
        H.nopena.theta<-4*theta*Hnopena%*%diag(theta)
        browser()
        
        trLCV.theta<-lava::tr(-V.theta%*%H.nopena.theta)
        trLCV<-lava::tr(-V%*%H.nopena)
        
        
        
        
      }else{
        H<-matrix(0,size_V,size_V)
        Hnopena<-matrix(0,size_V,size_V)
        trLCV<-0
        trLCV.theta<-0
      }
      
      list(V,H,trLCV,trLCV.theta)
    }
    
    parallel::stopCluster(clustpar)
    
    
    return(list(b=beta,
                lambda=lambda,
                gamma=gamma,
                v=res$V,
                H=res$H,
                trLCV=unlist(res$trLCV),
                trLCV.theta=unlist(res$trLCV.theta),
                fn.value=fn.value,
                ni=ni,
                istop=istop,
                ca=ca,
                cb=cb,
                rdm=rdm,
                fix0=fix0))
  }
}
