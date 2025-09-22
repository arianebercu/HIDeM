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
                NtimesPoints,seed){
  
  # define timePoints of prediction : 
  
  set.seed(seed = seed)
  

  if(Ypredmethod=="gauss"){
    
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
    return(data.frame(index=x,timePoints=timePoints
                      ))}))
  
  }
  
  print("Start running joint univarite models")
  
  Yall<-list()
  length(Yall)<-length(formLong)
  

  for(indice in 1:length(formLong)){

    # need to have all elements of joint
    # global variables otherwise error in predict
    print(paste0("For marker: ",names(functional_forms)[[indice]]))

    JMmodel<-JMbayes2::jm(formSurv, list(formLong[[indice]]), time_var = timeVar, 
                             functional_forms = functional_forms[[indice]],
                             n_iter =n_iter, n_burnin = n_burnin, n_thin =n_thin,
                             n_chains=n_chain, data_Surv = dataSurv,
                             cores=nproc,save_random_effects=T)
  
    
    #if(is.null(JMmodel)){stop("The JMbayes2 model for your marker could not be run, see above warnings.")}
  
  
    if(Ypredmethod=="gauss"){
   
    # if(nproc==1){
    # PredYx<-do.call(rbind,lapply(idsubjects,FUN=function(x){
    #   print(x)
    #   indexLongi<-which(dataLongi[,colnames(dataLongi)%in%id]==x)
    #   datai<-dataLongi[indexLongi,]
    #   
    #   datai<-merge(x=datai,y=dataSurv,by=id,all.x=T)
    #   timeJM<-timePointsdata[timePointsdata$index==x,colnames(timePointsdata)=="timePoints"]
    #   
    #   ND<-list(newdataL=dataLongi[indexLongi,],newdataE=datai)
    #   
    #   # we have no prediction if superior to visit times
    #   if(any(timeJM>max(dataLongi$visit))){print("timeJM superior")}
    #   
    #   timeJM<-ifelse(timeJM>max(dataLongi$visit),max(dataLongi$visit),
    #                  timeJM)
    #   
    #   predY<-predict(JMmodel, 
    #                  newdata = ND,process="longitudinal",
    #                  times=timeJM,
    #                  #return_newdata=T,
    #                  type_pred = "link",type="subject_specific",
    #                  control=list(n_samples=Nsample,
    #                               all_times=T,
    #                               return_mcmc=T))
    #     #$preds$y is the mean value over the n_samples
    #     
    #     Y<-data.frame(predY$newdata2$mcmc)
    #     colnames(Y)<-paste0("Nsample",c(1:Nsample))
    #     
    #     Y$time<-timeJM
    #     Y$ID<-x
    #     
    #    
    #   return(Y)
    #   
    # }))
    # }else{
    #   
    #     if(is.null(clustertype)){
    #       clustpar <- parallel::makeCluster(nproc)#, outfile="")
    #     }
    #     else{
    #       clustpar <- parallel::makeCluster(nproc, type=clustertype)#, outfile="")
    #     }
    # 
    #   PredYx<-do.call(rbind,parLapply(cl=clustpar,idsubjects,fun=function(x){
    #     print(x)
    #     indexLongi<-which(dataLongi[,colnames(dataLongi)%in%id]==x)
    #     datai<-dataLongi[indexLongi,]
    #     
    #     datai<-merge(x=datai,y=dataSurv,by=id,all.x=T)
    #     timeJM<-timePointsdata[timePointsdata$index==x,colnames(timePointsdata)=="timePoints"]
    #     
    #     ND<-list(newdataL=dataLongi[indexLongi,],newdataE=datai)
    #     
    #     # we have no prediction if superior to visit times
    #     if(any(timeJM>max(dataLongi$visit))){print("timeJM superior")}
    #     
    #     timeJM<-ifelse(timeJM>max(dataLongi$visit),max(dataLongi$visit),
    #                    timeJM)
    #     
    #     predY<-predict(JMmodel, 
    #                    newdata = ND,process="longitudinal",
    #                    times=timeJM,
    #                    #return_newdata=T,
    #                    #type_pred = "link",
    #                    type="subject_specific",
    #                    control=list(n_samples=Nsample,
    #                                 all_times=T,
    #                                 return_mcmc=T))
    #     #$preds$y is the mean value over the n_samples
    #     
    #     Y<-data.frame(predY$newdata2$mcmc)
    #     colnames(Y)<-paste0("Nsample",c(1:Nsample))
    #     
    #     Y$time<-timeJM
    #     Y$ID<-x
    #     
    #     
    #     return(Y)
    #     
    #   }))
    #  
    # }
    # 
    # 
    # PredYx$Outcome<-names(functional_forms)[[indice]]
    # colnames(PredYx)[colnames(PredYx)%in%"ID"]<-id
    # colnames(PredYx)[colnames(PredYx)%in%"time"]<-timeVar
    # 
      
     
      
      N<-length(unique(dataLongi[,id]))
      newdataLongi<-timePointsdata
      colnames(newdataLongi)<-c(id,timeVar)
      
      terms_FE<-JMmodel$model_info$terms$terms_FE[[1]]
      terms_RE<-JMmodel$model_info$terms$terms_RE[[1]]
      
      X<-model.matrix(reformulate(attr(terms(terms_FE), "term.labels")),data=newdataLongi)
      Z <- model.matrix(reformulate(attr(terms(terms_RE), "term.labels")),data=newdataLongi)
      
      betas<-do.call(rbind,JMmodel$mcmc$betas1)
      # ind_RE <- JMmodel$model_data$ind_RE
      # idd<-dataLongi[[JMmodel$model_info$var_names$idVar]]
      
      d1<-dim(JMmodel$mcmc$b[[1]])[3]
      d2<-JMmodel$control$n_chains
      b_mat <- array(0.0, dim = c(dim(JMmodel$mcmc$b[[1]])[1:2], d1*d2))
      b_mat[, , seq(1, d1)]<-JMmodel$mcmc$b[[1]]
      for( i in 2:JMmodel$control$n_chains){
        b_mat[, , seq(d1*(i-1)+1, d1*i)] <- JMmodel$mcmc$b[[i]]
      }
      
      # K <- length(ind_RE)
      # M <- dim(b_mat)[3L]
      M<-dim(betas)[1]
      NtimesPoints<-ifelse(truncated==F,256,271)
      nn<-rep(NtimesPoints,N)
      ends   <- cumsum(nn)
      starts <- ends - nn + 1
      
      if(Nsample>1){
        
      idNsample<-sample(x=c(1:M),size=Nsample-1)

      Fixed <- X %*% t(betas[idNsample, , drop = FALSE])   # (n x m)
      
      # Terme aléatoire, on remplit directement un tableau vide
      Random_all <- matrix(0, nrow = nrow(Z), ncol = Nsample-1)
      Random_mean <- matrix(0,nrow=nrow(Z),ncol=1)
      
      for (j in seq_len(N)) {
        rows <- starts[j]:ends[j]
        Zj   <- Z[rows, , drop = FALSE]                        # (nn[j] x q)
        Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
        Random_all[rows, ] <- Zj %*% Bj[1,,]
        Random_mean[rows, ] <- Zj %*% JMmodel$statistics$Mean$b[j,]
      }
      
      PredYx <- Fixed + Random_all
      
      }else{
        PredYx<-NULL
        
        # Terme aléatoire, on remplit directement un tableau vide
        Random_mean <- matrix(0,nrow=nrow(Z),ncol=1)
        
        for (j in seq_len(N)) {
          rows <- starts[j]:ends[j]
          Zj   <- Z[rows, , drop = FALSE]
          Random_mean[rows, ] <- Zj %*% JMmodel$statistics$Mean$b[j,]
        }
        
        
      }
      
      PredYmean<-X%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean
      PredYx<-cbind(PredYmean,PredYx)
      # PredYxx <- do.call(cbind, lapply(idNsample, function(i) {
      #     fixed <- X %*% betas[i,]
      # 
      #     # concaténer les contributions en une seule fois
      #     random <- do.call(rbind, lapply(seq_along(starts), function(j) {
      #       Z[starts[j]:ends[j], , drop = FALSE] %*% b_mat[j, , i]
      #     }))
      # 
      #     fixed + random
      #   }))

      Outcome<-names(functional_forms)[[indice]]
      PredYx<-cbind(newdataLongi,Outcome,PredYx)
      colnames(PredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
      
    
    }else{
      
      # maxmimum need to not exceed time observed
      
      # times<-seq(min(t0),max(dataLongi[,colnames(dataLongi)%in%timeVar]),length.out=NtimesPoints)
      # if(nproc==1){
      # PredYx<-do.call(rbind,lapply(idsubjects,FUN=function(x){
      # 
      #   indexLongi<-which(dataLongi[,colnames(dataLongi)%in%id]==x)
      #   datai<-dataLongi[indexLongi,]
      #   
      #   datai<-merge(x=datai,y=dataSurv,by=id,all.x=T)
      #   
      #   ND<-list(newdataL=dataLongi[indexLongi,],newdataE=datai)
      #   
      #   predY<-predict(JMmodel, 
      #                  newdata = ND,
      #                  process="longitudinal",
      #                  times=times, #cannot put return_data=T with return mcm = T will create error
      #                  #type_pred = "link",
      #                  type="subject_specific",
      #                  
      #                  control=list(n_samples=Nsample,
      #                               all_times=T, # do not forget to compute for all values
      #                               return_mcmc=T))
      #   #$preds$y is the mean value over the n_samples
      #   
      #   Y<-data.frame(predY$newdata2$mcmc)
      #   colnames(Y)<-paste0("Nsample",c(1:Nsample))
      #   Y$ID<-x
      #   return(Y)
      #   
      # }))
      # }else{
      #   if(is.null(clustertype)){
      #     clustpar <- parallel::makeCluster(nproc)#, outfile="")
      #   }
      #   else{
      #     clustpar <- parallel::makeCluster(nproc, type=clustertype)#, outfile="")
      #   }
      #   
      #   PredYx<-do.call(rbind,parLapply(cl=clustpar,idsubjects,fun=function(x){
      #     
      #     indexLongi<-which(dataLongi[,colnames(dataLongi)%in%id]==x)
      #     datai<-dataLongi[indexLongi,]
      #     
      #     datai<-merge(x=datai,y=dataSurv,by=id,all.x=T)
      #     
      #     ND<-list(newdataL=dataLongi[indexLongi,],newdataE=datai)
      #     
      #     predY<-predict(JMmodel, 
      #                    newdata = ND,
      #                    process="longitudinal",
      #                    times=times, #cannot put return_data=T with return mcm = T will create error
      #                    #type_pred = "link",
      #                    type="subject_specific",
      #                    
      #                    control=list(n_samples=Nsample,
      #                                 all_times=T, # do not forget to compute for all values
      #                                 return_mcmc=T))
      #     #$preds$y is the mean value over the n_samples
      #     
      #     Y<-data.frame(predY$newdata2$mcmc)
      #     colnames(Y)<-paste0("Nsample",c(1:Nsample))
      #     Y$ID<-x
      #     return(Y)
      #     
      #   }))
      
      times<-seq(min(t0),max(dataLongi[,colnames(dataLongi)%in%timeVar]),length.out=NtimesPoints)
      
      N<-length(unique(dataLongi[,id]))
      
      newdataLongi<-do.call(rbind,lapply(unique(dataLongi[,id]),FUN=function(x){
        data.frame(id=rep(x,NtimesPoints),
                   time=times)}))
      colnames(newdataLongi)<-c(id,timeVar)
      
      terms_FE<-JMmodel$model_info$terms$terms_FE[[1]]
      terms_RE<-JMmodel$model_info$terms$terms_RE[[1]]
      
      X<-model.matrix(reformulate(attr(terms(terms_FE), "term.labels")),data=newdataLongi)
      Z <- model.matrix(reformulate(attr(terms(terms_RE), "term.labels")),data=newdataLongi)
      
      betas<-do.call(rbind,JMmodel$mcmc$betas1)
      # ind_RE <- JMmodel$model_data$ind_RE
      # idd<-dataLongi[[JMmodel$model_info$var_names$idVar]]
      
      d1<-dim(JMmodel$mcmc$b[[1]])[3]
      d2<-JMmodel$control$n_chains
      b_mat <- array(0.0, dim = c(dim(JMmodel$mcmc$b[[1]])[1:2], d1*d2))
      b_mat[, , seq(1, d1)]<-JMmodel$mcmc$b[[1]]
      for( i in 2:JMmodel$control$n_chains){
        b_mat[, , seq(d1*(i-1)+1, d1*i)] <- JMmodel$mcmc$b[[i]]
      }
      
      # K <- length(ind_RE)
      # M <- dim(b_mat)[3L]
      M<-dim(betas)[1]
      nn<-rep(NtimesPoints,N)
      ends   <- cumsum(nn)
      starts <- ends - nn + 1
      
      if(Nsample>1){
      idNsample<-sample(x=c(1:M),size=Nsample-1)
      
      Fixed <- X %*% t(betas[idNsample, , drop = FALSE])   # (n x m)
    
      # Terme aléatoire, on remplit directement un tableau vide
      Random_all <- matrix(0, nrow = nrow(Z), ncol = Nsample-1)
      Random_mean <- matrix(0,nrow=nrow(Z),ncol=1)
      for (j in seq_len(N)) {
        rows <- starts[j]:ends[j]
        Zj   <- Z[rows, , drop = FALSE]                        # (nn[j] x q)
        Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
        Random_all[rows, ] <- Zj %*% Bj[1,,]
        Random_mean[rows, ] <- Zj %*% JMmodel$statistics$Mean$b[j,]
      }
     
      PredYx <- Fixed + Random_all
      }else{
        PredYx<-NULL
        
        Random_mean <- matrix(0,nrow=nrow(Z),ncol=1)
        for (j in seq_len(N)) {
          rows <- starts[j]:ends[j]
          Zj   <- Z[rows, , drop = FALSE]
          Random_mean[rows, ] <- Zj %*% JMmodel$statistics$Mean$b[j,]
        }
        }
     
      PredYmean<-X%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean
      PredYx<-cbind(PredYmean,PredYx)
      
      #add mean value of fixed and random effects 
      
      
      # PredYxx <- do.call(cbind, lapply(idNsample, function(i) {
      #     fixed <- X %*% betas[i,]
      # 
      #     # concaténer les contributions en une seule fois
      #     random <- do.call(rbind, lapply(seq_along(starts), function(j) {
      #       Z[starts[j]:ends[j], , drop = FALSE] %*% b_mat[j, , i]
      #     }))
      # 
      #     fixed + random
      #   }))

      Outcome<-names(functional_forms)[[indice]]
      PredYx<-cbind(newdataLongi,Outcome,PredYx)
      colnames(PredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
      
     
      }
      
      
    
    Yall[[indice]]<- PredYx
  }
  
  Yall<-do.call(rbind,Yall)
  print("End of running joint univarite models")
  return(Yall[,colnames(Yall)%in%c(id,timeVar,"Outcome",paste0("Sample_",c(1:Nsample)))])
  
}






