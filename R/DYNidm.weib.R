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
##' @param weib the form of the weibull parameters 
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
#' @useDynLib HIDeM

DYNidm.weib<-function(b,fix0,size_V,
                      clustertype,epsa,epsb,epsd,nproc,maxiter,
                      ctime,N,
                      ve01,ve02,ve12,dimnva01,dimnva02,dimnva12,nvat01,nvat02,nvat12,
                      t0,t1,t2,t3,idd,idm,ts,troncature, dataY,
                      Longitransition,
                      NtimePoints,
                      Ypredmethod,
                      timeVar,ynames,id,
                      outcome,outcome01,outcome02,outcome12,
                      p01,p02,p12,
                      dimp01,dimp02,dimp12){
  
  bfix<-b[fix0==1]
  b<-b[fix0==0]
  npm<-size_V-sum(fix0)
  
  
  Nsample<-dim(dataY)[2]-3

  
  if(Ypredmethod=="equi"){
  
  # if prediction did not work 
  # check which time to keep : 
  time<-dataY[dataY$Outcome%in%ynames[1] & dataY[,colnames(dataY)%in%id]==unique(dataY[,colnames(dataY)%in%id])[1],colnames(dataY)%in%timeVar]
  valid<-unlist(lapply(time,FUN=function(x){
    dd<-dataY[dataY[,colnames(dataY)%in%timeVar]==x,]
    if(dim(dd)[1]!=N*length(ynames)){return(F)}else{return(T)}
  }))
  
  time<-time[valid]
  NtimePoints<-length(time)
  
  dataY<-dataY[dataY[,colnames(dataY)%in%timeVar]%in%time,]
  dataY$Outcome<-as.character(dataY$Outcome)
  # attention if NtimePoints equidistant with INLA then NtimePoints takes 
  # need ID to be numeric -- then 
  dataY[,colnames(dataY)%in%id]<-as.numeric(dataY[,colnames(dataY)%in%id])
  # to keep tracks of time order for each individual 
  dataY$order<-as.numeric(ave(dataY[,colnames(dataY)%in%id], cumsum(c(TRUE, diff(dataY[,colnames(dataY)%in%id]) != 0)), FUN = seq_along))
  
  if(length(outcome01)>=1){
  y01<-dataY[dataY$Outcome%in%outcome01,]
  # order  by individual and timeline 
  y01<-y01[order(y01[,colnames(y01)%in%id],y01$order),]
  
  
  }else{
    y01<-as.double(rep(0,N*NtimePoints))
  }
  
  if(length(outcome02)>=1){
    y02<-dataY[dataY$Outcome%in%outcome02,]
    y02<-y02[order(y02$ID,y02$order),]
    
  }else{
    y02<-as.double(rep(0,N*NtimePoints))
  }
  
  if(length(outcome12)>=1){
    y12<-dataY[dataY$Outcome%in%outcome12,]
    # order  by individual and timeline 
    y12<-y12[order(y12$ID,y12$order),]
    
  }else{
    y12<-as.double(rep(0,N*NtimePoints))
  }

 
 out<-list()
 length(out)<-Nsample
 
 # # error comes from ctime = 1 ? 
 # 
 # idd<-which(ctime==1)
 # ctime<-ctime[-idd]
 # t0<-t0[-idd]
 # t1<-t1[-idd]
 # t2<-t2[-idd]
 # t3<-t3[-idd]
 # N<-N-length(idd)
 # y01<-y01[!y01$ID%in%idd,]
 # y02<-y02[!y02$ID%in%idd,]
 # y12<-y12[!y12$ID%in%idd,]
 
 #cl <- makeCluster(2)
 #registerDoParallel(cl)
 #stopCluster(cl)
 
 #browser()
  for(k in 1:Nsample){
  
  out[[k]]<- tryCatch({ marqLevAlg::mla(b=b,
                        fn=DYNidmlLikelihoodweib,
                        epsa=epsa,
                        epsb=epsb,
                        epsd=epsd,
                        nproc=nproc,
                        maxiter=maxiter,
                        print.info=T,
                        minimize=F,
                        npm=npm,
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
                        y01=y01[,colnames(y01)%in%paste0("Sample_",k)],
                        y02=y02[,colnames(y02)%in%paste0("Sample_",k)],
                        y12=y12[,colnames(y12)%in%paste0("Sample_",k)],
                        p01=p01,
                        p02=p02,
                        p12=p12,
                        dimp01=dimp01,
                        dimp02=dimp02,
                        dimp12=dimp12,
                        Ntime=NtimePoints,
                        time=time)
  }, error = function(e) {
    # Return NULL on error to skip this patient
    NULL
  })
  if(out[[k]]$istop==1){
    b<-out[[k]]$b
  }
 
  }
  }else{
    
    # check if predictions could be performed 
    #browser()
    if(troncature==T){
      NtimePoints<-271
    }else{
      NtimePoints<-256
    }
    
    for( k in outcome){
      subdata<-dataY[dataY$Outcome==k,]
      x<-table(subdata[,colnames(subdata)%in%id])
      if(any(x!=NtimePoints)){stop("Prediction of marker ",k," could not be perform for each quadrature points, try Ypredmethod equi")}
      
    }
    
    dataY$Outcome<-as.character(dataY$Outcome)
    # attention if NtimePoints equidistant with INLA then NtimePoints takes 
    # need ID to be numeric -- then 
    dataY[,colnames(dataY)%in%id]<-as.numeric(dataY[,colnames(dataY)%in%id])
    # to keep tracks of time order for each individual 
    dataY$order<-as.numeric(ave(dataY[,colnames(dataY)%in%id], cumsum(c(TRUE, diff(dataY[,colnames(dataY)%in%id]) != 0)), FUN = seq_along))
    

    if(length(outcome01)>=1){
      y01<-dataY[dataY$Outcome%in%outcome01,]
      # order  by individual and timeline 
      y01<-y01[order(y01[,colnames(y01)%in%id],y01$order),]
      
      
    }else{
      y01<-as.double(rep(0,N*NtimePoints))
    }
    
    if(length(outcome02)>=1){
      y02<-dataY[dataY$Outcome%in%outcome02,]
      # order  by individual and timeline 
      y02<-y02[order(y02$ID,y02$order),]
      
    }else{
      y02<-as.double(rep(0,N*NtimePoints))
    }
    
    if(length(outcome12)>=1){
      y12<-dataY[dataY$Outcome%in%outcome12,]
      # order  by individual and timeline 
      y12<-y12[order(y12$ID,y12$order),]
      
    }else{
      y12<-as.double(rep(0,N*NtimePoints))
    }
    
    
    out<-list()
    length(out)<-Nsample
#browser()

    for(k in 1:Nsample){
      print(paste0("Estimating illness-death model on sample ",k))
      out[[k]]<- tryCatch({ marqLevAlg::mla(b=b,
                                            fn=gaussDYNidmlLikelihoodweib,
                                            epsa=epsa,
                                            epsb=epsb,
                                            epsd=epsd,
                                            nproc=nproc,
                                            clustertype=clustertype,
                                            maxiter=maxiter,
                                            minimize=F,
                                            print.info = T,
                                            npm=npm,
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
                                            y01=y01[,colnames(y01)%in%paste0("Sample_",k)],
                                            y02=y02[,colnames(y02)%in%paste0("Sample_",k)],
                                            y12=y12[,colnames(y12)%in%paste0("Sample_",k)],
                                            p01=p01,
                                            p02=p02,
                                            p12=p12,
                                            dimp01=dimp01,
                                            dimp02=dimp02,
                                            dimp12=dimp12,
                                            Ntime=NtimePoints)
      }, error = function(e) {
        # Return NULL on error to skip this patient
        NULL
      })
      if(out[[k]]$istop==1){
        b<-out[[k]]$b
      }
      
    }
    
  }
  
  
  
  
  return(out)
  
  
  
  
  
}

