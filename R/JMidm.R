### Code:
##' @title Calculate predictions for time-depend covariates using INLA
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  


JMidm<-function(timeVar,
                functional_forms, 
                truncated,
                formLong,
                formSurv,
                dataSurv,
                dataLongi,
                id,
                Nsample,
                n_iter,
                n_burnin,
                n_thin,
                n_chain,
                nproc,t0,t1,t2,t3,idm,idd,
                clustertype,
                Ypredmethod,
                NtimesPoints){
  
  # define timePoints of prediction : 
  

  idsubjects<-unique(dataSurv[,colnames(dataSurv)%in%id])
 
  
  timePointsdata<-do.call(rbind, lapply(idsubjects,FUN=function(x){
    index<-which(dataSurv[seq(1,dim(dataSurv)[1],by=2),colnames(dataSurv)%in%id]==x)
    max.int<-ifelse(idm[index]==1,t2[index],t3[index])
    min.int<-t1[index]
    
    timePoints<-gauss_kronrod_points(lower.intdouble=min.int,
                                     upper.intdouble=max.int,
                                     end.time=t3[index],
                                     truncated=truncated,
                                     entry.time=t0[index])
    return(data.frame(timePoints=timePoints,
                      index=x))}))
  
  
  print("Start running joint univarite models")
  
  Yall<-list()
  length(Yall)<-length(formLong)
  
  #browser()
  for(indice in 1:length(formLong)){

    # need to have all elements of joint
    # global variables otherwise error in predict

    JMmodel<-JMbayes2::jm(formSurv, list(formLong[[indice]]), time_var = timeVar, 
                             functional_forms = functional_forms[[indice]],
                             n_iter =n_iter, n_burnin = n_burnin, n_thin =n_thin,
                             n_chains=n_chain, data_Surv = dataSurv,
                             cores=nproc,save_random_effects=T)
  
    
    #if(is.null(JMmodel)){stop("The JMbayes2 model for your marker could not be run, see above warnings.")}
    #browser()

    if(Ypredmethod=="gauss"){
   
    if(nproc==1){
    PredYx<-do.call(rbind,lapply(idsubjects,FUN=function(x){
      print(x)
      indexLongi<-which(dataLongi[,colnames(dataLongi)%in%id]==x)
      datai<-dataLongi[indexLongi,]
      
      datai<-merge(x=datai,y=dataSurv,by=id,all.x=T)
      timeJM<-timePointsdata[timePointsdata$index==x,colnames(timePointsdata)=="timePoints"]
      
      ND<-list(newdataL=dataLongi[indexLongi,],newdataE=datai)
      
      # we have no prediction if superior to visit times
      if(any(timeJM>max(dataLongi$visit))){print("timeJM superior")}
      
      timeJM<-ifelse(timeJM>max(dataLongi$visit),max(dataLongi$visit),
                     timeJM)
      
      predY<-predict(JMmodel, 
                     newdata = ND,process="longitudinal",
                     times=timeJM,
                     #return_newdata=T,
                     type_pred = "link",type="subject_specific",
                     control=list(n_samples=Nsample,
                                  all_times=T,
                                  return_mcmc=T))
        #$preds$y is the mean value over the n_samples
        
        Y<-data.frame(predY$newdata2$mcmc)
        colnames(Y)<-paste0("Nsample",c(1:Nsample))
        
        Y$time<-timeJM
        Y$ID<-x
        
       
      return(Y)
      
    }))
    }else{
      
        if(is.null(clustertype)){
          clustpar <- parallel::makeCluster(nproc)#, outfile="")
        }
        else{
          clustpar <- parallel::makeCluster(nproc, type=clustertype)#, outfile="")
        }
  
      PredYx<-do.call(rbind,parLapply(cl=clustpar,idsubjects,fun=function(x){
        print(x)
        indexLongi<-which(dataLongi[,colnames(dataLongi)%in%id]==x)
        datai<-dataLongi[indexLongi,]
        
        datai<-merge(x=datai,y=dataSurv,by=id,all.x=T)
        timeJM<-timePointsdata[timePointsdata$index==x,colnames(timePointsdata)=="timePoints"]
        
        ND<-list(newdataL=dataLongi[indexLongi,],newdataE=datai)
        
        # we have no prediction if superior to visit times
        if(any(timeJM>max(dataLongi$visit))){print("timeJM superior")}
        
        timeJM<-ifelse(timeJM>max(dataLongi$visit),max(dataLongi$visit),
                       timeJM)
        
        predY<-predict(JMmodel, 
                       newdata = ND,process="longitudinal",
                       times=timeJM,
                       #return_newdata=T,
                       #type_pred = "link",
                       type="subject_specific",
                       control=list(n_samples=Nsample,
                                    all_times=T,
                                    return_mcmc=T))
        #$preds$y is the mean value over the n_samples
        
        Y<-data.frame(predY$newdata2$mcmc)
        colnames(Y)<-paste0("Nsample",c(1:Nsample))
        
        Y$time<-timeJM
        Y$ID<-x
        
        
        return(Y)
        
      }))
     
    }
    
   
    PredYx$outcome<-names(functional_forms)[[indice]]
    colnames(PredYx)[colnames(PredYx)%in%"ID"]<-id
    colnames(PredYx)[colnames(PredYx)%in%"time"]<-timeVar
    
    
    }else{
      
      # maxmimum need to not exceed time observed
      
      times<-seq(min(t0),max(dataLongi[,colnames(dataLongi)%in%timeVar]),length.out=NtimesPoints)
      if(nproc==1){
      PredYx<-do.call(rbind,lapply(idsubjects,FUN=function(x){
    
        indexLongi<-which(dataLongi[,colnames(dataLongi)%in%id]==x)
        datai<-dataLongi[indexLongi,]
        
        datai<-merge(x=datai,y=dataSurv,by=id,all.x=T)
        
        ND<-list(newdataL=dataLongi[indexLongi,],newdataE=datai)
        
        predY<-predict(JMmodel, 
                       newdata = ND,
                       process="longitudinal",
                       times=times, #cannot put return_data=T with return mcm = T will create error
                       #type_pred = "link",
                       type="subject_specific",
                       
                       control=list(n_samples=Nsample,
                                    all_times=T, # do not forget to compute for all values
                                    return_mcmc=T))
        #$preds$y is the mean value over the n_samples
        
        Y<-data.frame(predY$newdata2$mcmc)
        colnames(Y)<-paste0("Nsample",c(1:Nsample))
        Y$ID<-x
        return(Y)
        
      }))
      }else{
        if(is.null(clustertype)){
          clustpar <- parallel::makeCluster(nproc)#, outfile="")
        }
        else{
          clustpar <- parallel::makeCluster(nproc, type=clustertype)#, outfile="")
        }
        
        PredYx<-do.call(rbind,parLapply(cl=clustpar,idsubjects,fun=function(x){
          
          indexLongi<-which(dataLongi[,colnames(dataLongi)%in%id]==x)
          datai<-dataLongi[indexLongi,]
          
          datai<-merge(x=datai,y=dataSurv,by=id,all.x=T)
          
          ND<-list(newdataL=dataLongi[indexLongi,],newdataE=datai)
          
          predY<-predict(JMmodel, 
                         newdata = ND,
                         process="longitudinal",
                         times=times, #cannot put return_data=T with return mcm = T will create error
                         #type_pred = "link",
                         type="subject_specific",
                         
                         control=list(n_samples=Nsample,
                                      all_times=T, # do not forget to compute for all values
                                      return_mcmc=T))
          #$preds$y is the mean value over the n_samples
          
          Y<-data.frame(predY$newdata2$mcmc)
          colnames(Y)<-paste0("Nsample",c(1:Nsample))
          Y$ID<-x
          return(Y)
          
        }))
        
      }
      
      colnames(PredYx)[colnames(PredYx)%in%"ID"]<-id
      
      PredYx$outcome<-names(functional_forms)[[indice]]
      
      PredYx$time<-rep(times,length(idsubjects))
      colnames(PredYx)[colnames(PredYx)%in%"time"]<-timeVar
    }
    
    Yall[[indice]]<- PredYx
  }
  
  Yall<-do.call(rbind,Yall)
  print("End of running joint univarite models")
  return(Yall[,colnames(Yall)%in%c(id,timeVar,"outcome",paste0("Nsample",c(1:Nsample)))])
  
}