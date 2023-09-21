idm.step.sequential<-function(b,clustertype,epsa,epsb,epsd,print.info,nproc,
                                 maxiter,size_V,size_spline,noVar,bfix,fix0,knots01,knots02,knots12,
                                 ctime,N,nknots01,nknots02,nknots12,ve01,ve02,ve12,
                                 dimnva01,dimnva02,dimnva12,nvat01,nvat02,nvat12,t0,t1,
                                 t2,t3,troncature,gauss.point,option.sequential,penalty,lambda,gamma,
                              ca, cb,rdm){

  # First step : 
  # fix B.spline parameters at 0 or its values with first lambda :  
  fix00<-fix0
  maxstep<-0
  
  if(ca<=epsa & cb<=epsb & rdm>epsd){
  # put fix at 1 if necessary 
  btot<-rep(NA,size_V)
  btot[which(fix0==0)]<-b
  btot[which(fix0==1)]<-bfix
  
  nparspline<-sum(fix0[1:size_spline]==0) 
  
  # if we still have parameters to fix possibly we evaluate it 
  if(nparspline>0){
    
    fix0<-ifelse(fix0==1 | ((abs(btot)<=option.sequential$cutoff) & fix0==0),1,0)
    
    if(min(noVar)==0){ # if we have explanatory variable we look at their fix 
      fix0[(size_spline+1):length(fix0)]<-fix00[(size_spline+1):length(fix0)]
      posfix<-which(fix0==1)
      posfix<-posfix[which(posfix>size_spline)]
      bfix<-c(rep(0,sum(fix0[1:size_spline])),btot[posfix])
    }else{
      bfix<-rep(0,sum(fix0[1:size_spline]))
    }
  }
  
  b<-btot[which(fix0==0)]
  }
  
    
    while(maxstep<maxiter){
      

      
      out<-marqLevAlg::mla(b=b,
               m=length(b),
               fn=idmlLikelihoodpena,
               clustertype=clustertype,
               epsa=epsa,
               epsb=epsb,
               epsd=epsd,
               print.info = print.info,
               nproc=nproc,
               maxiter=option.sequential$step,
               minimize=F,
               np0=length(b),
               npar0=size_V,
               bfix0=bfix,
               fix0=fix0,
               zi010=knots01,
               zi020=knots02,
               zi120=knots12,
               c0=ctime,
               no0=N,
               nz010=nknots01,
               nz020=nknots02,
               nz120=nknots12,
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
      

      if(out$istop==4){
        stop("Problem in the loglikelihood computation.")
      }
      
      maxstep<-maxstep+out$ni
      b<-out$b
      
      
      
      while(out$ca<=epsa & out$cb<=epsb & out$rdm>epsd){
        
        
        if(maxstep>=maxiter){break}
        
        btot<-rep(NA,size_V)
        btot[which(fix0==0)]<-b
        btot[which(fix0==1)]<-bfix
        
        nparspline<-sum(fix0[1:size_spline]==0) 
        
        # if we still have parameters to fix possibly we evaluate it 
        if(nparspline>0){
          
          fix0<-ifelse(fix0==1 | ((abs(btot)<=option.sequential$cutoff) & fix0==0),1,0)
          
          if(min(noVar)==0){ # if we have explanatory variable we look at their fix 
            fix0[(size_spline+1):length(fix0)]<-fix00[(size_spline+1):length(fix0)]
            posfix<-which(fix0==1)
            posfix<-posfix[which(posfix>size_spline)]
            bfix<-c(rep(0,sum(fix0[1:size_spline])),btot[posfix])
          }else{
              bfix<-rep(0,sum(fix0[1:size_spline]))
            }
          }
          
          b<-btot[which(fix0==0)]
          
        
        
        out<-marqLevAlg::mla(b=b,
                 m=length(b),
                 fn=idmlLikelihoodpena,
                 clustertype=clustertype,
                 epsa=epsa,
                 epsb=epsb,
                 epsd=epsd,
                 print.info = print.info,
                 nproc=nproc,
                 maxiter=option.sequential$step,
                 minimize=F,
                 np0=length(b),
                 npar0=size_V,
                 bfix0=bfix,
                 fix0=fix0,
                 zi010=knots01,
                 zi020=knots02,
                 zi120=knots12,
                 c0=ctime,
                 no0=N,
                 nz010=nknots01,
                 nz020=nknots02,
                 nz120=nknots12,
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
        
        
        if(out$istop==4){
          stop("Problem in the loglikelihood computation.")
        }
        
        
        
        maxstep<-maxstep+out$ni
        b<-out$b
        
      }
      
      if(out$ca<=epsa & out$cb<=epsb & out$rdm<=epsd){
        
        break
      }
      
    }
    
    
    beta<-rep(NA,size_V)
    beta[which(fix0==0)]<-out$b
    beta[which(fix0==1)]<-bfix
    
    
    npm<-length(out$b)
    posfix<-which(fix0==1)
    
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
                       nz010=nknots01,
                       nz020=nknots02,
                       nz120=nknots12,
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
      trLCV.theta<-0
    }
    

    return(list(b=as.matrix(beta),
                lambda=lambda,
                gamma=gamma,
                v=V,
                H=H,
                trLCV=trLCV,
                trLCV.theta=trLCV.theta,
                fn.value=out$fn.value,
                ni=maxstep,
                istop=out$istop,
                grad=out$grad,
                ca=out$ca,
                cb=out$cb,
                rdm=out$rdm,
                fix0=fix0))
}