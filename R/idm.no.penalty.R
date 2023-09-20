idm.no.penalty<-function(b,clustertype,epsa,epsb,epsd,print.info,nproc,maxiter,size_V,size_spline,noVar,bfix,
                         fix0,knots01,knots02,knots12,ctime,N,nknots01,nknots02,nknots12,
                         ve01,ve02,ve12,dimnva01,dimnva02,dimnva12,nvat01,nvat02,nvat12,
                         t0,t1,t2,t3,troncature,gauss.point,step.sequential,option.sequential){
  fix00<-fix0
  
  if(step.sequential==F){
    
    
    out<-mla(b=b,
             m=length(b),
             fn=idmlLikelihood,
             clustertype=clustertype,
             epsa=epsa,
             epsb=epsb,
             epsd=epsd,
             print.info = print.info,
             nproc=nproc,
             maxiter=maxiter,
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
             gausspoint0=gauss.point
    )
    
    
    if(out$istop==4){
      stop("Problem in the loglikelihood computation.")
    }
    
    ni<-out$ni
    #simplifier modele ou two step
    
  }else{
    
    if(option.sequential$min>0){
      
      out<-mla(b=b,
               m=length(b),
               fn=idmlLikelihood,
               clustertype=clustertype,
               epsa=epsa,
               epsb=epsb,
               epsd=epsd,
               print.info = print.info,
               nproc=nproc,
               maxiter=option.sequential$min,
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
               gausspoint0=gauss.point
      )
      
      
      if(out$istop==4){
        stop("Problem in the loglikelihood computation.")
      }
      
      maxstep<-out$ni
      ni<-maxstep
      b<-out$b
      out_first_mla<-out$istop
      # change condition from : out$ca<=epsa & out$cb<=epsb & out$rdm>epsd
      # to : 
      
      if(out$rdm>epsd){
        
        
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
    }else{
      
      maxstep<-0
      ni<-maxstep
      out_first_mla<-2
    }
    
    
    while(maxstep<maxiter & out_first_mla!=1){
      
      
      out<-mla(b=b,
               m=length(b),
               fn=idmlLikelihood,
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
               gausspoint0=gauss.point
      )
      
      #browser()
      
      if(out$istop==4){
        stop("Problem in the loglikelihood computation.")
      }
      
      maxstep<-maxstep+out$ni
      b<-out$b
      
      if(maxstep>=maxiter){break}
      
      # replace while(out$ca<=epsa & out$cb<=epsb & out$rdm>epsd) by : 
      if(out$ca<=epsa & out$cb<=epsb & out$rdm>epsd){
        
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
      
      #ni<-ni+out$ni
      
      
      
      if(out$ca<=epsa & out$cb<=epsb & out$rdm<=epsd){break}
      
    }
    ni<-maxstep
  }
  
  return(list(b=out$b,
              fn.value=out$fn.value,
              ni=ni,
              istop=out$istop,
              v=out$v,
              grad=out$grad,
              ca=out$ca,
              cb=out$cb,
              rdm=out$rdm,
              bfix=bfix,
              fix0=fix0))
}