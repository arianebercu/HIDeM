### Code:
##' @title Calculate predictions for time-depend covariates using INLA
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  


INLAidm<-function(timeVar,family,basRisk,assoc,
                  truncated,formLong,formSurv,dataSurv,dataLongi,id,
                  Nsample,nproc,t0,t1,t2,t3,
                  idm,idd,clustertype,Ypredmethod,NtimesPoints){
  
  # define timePoints of prediction : 

browser()
  idsubjects<-unique(dataSurv[,colnames(dataSurv)%in%id])
  
  if(Ypredmethod=="gauss"){
  
  timePointsdata<-do.call(rbind, lapply(idsubjects,FUN=function(x){
    index<-which(dataSurv[,colnames(dataSurv)%in%id]==x)
    max.int<-ifelse(idm[index]==1,t2[index],t3[index])
    min.int<-t1[index]
  
    
    timePoints<-gauss_kronrod_points(lower.intdouble=min.int,
                                     upper.intdouble=max.int,
                                     end.time=t3[index],
                                     truncated=truncated,
                                     entry.time=t0[index])
    return(data.frame(index=x,timePoints=timePoints))}))
  
  
  }else{
    
  
    timePointsdata<-do.call(rbind, lapply(idsubjects,FUN=function(x){
      return(data.frame(index=x,timePoints=seq(min(t0),max(dataLongi[,colnames(dataLongi)%in%timeVar]))))}))
  }

  browser()
  
  ## augmentation of the data 
  colnames(timePointsdata)<-c(id,timeVar)
  dataLongi_augmented<-merge(dataLongi,timePointsdata,by=c(id,timeVar),all.x=T,all.y=T)
  rownames(dataLongi_augmented)<-NULL
  dataLongi_augmented<-dataLongi_augmented[order(dataLongi_augmented[,colnames(dataLongi_augmented)%in%id],
                                                 dataLongi_augmented[,colnames(dataLongi_augmented)%in%timeVar]),]
  
  print("Start running joint univarite models")

  Yall<-list()
  length(Yall)<-length(formLong)
  dataSurv<-dataSurv[order(dataSurv[,colnames(dataSurv)%in%id]),]
  
 # formSurvinla<<-formSurv
#  dataLongiinla<<-dataLongi_augmented
#  dataSurvinla<<-dataSurv
#  idinla<<-id
#  timeVarinla<<-timeVar
  

    for(indice in 1:length(formLong)){

      # need to have all elements of joint
      # global variables otherwise error in predict
     # formLonginla<<-formLong[[indice]]
    #  familyinla<<-family[indice]
    #  basRiskinla<<-basRisk[indice]
    #  associnla<<-assoc[[indice]]
    
      # cannot have lightmode as need BLUP

      INLAmodel<-INLAjoint::joint(formSurv = formSurv,
                                       formLong = formLong[[indice]],
                                       dataLong = dataLongi_augmented, dataSurv=dataSurv, id = id, timeVar = timeVar,
                                       family = family[indice],
                                       basRisk = basRisk[indice], NbasRisk = 15, assoc = assoc[[indice]],
                                       control=list(int.strategy="eb"))
    

      browser()
      
      if(is.null(INLAmodel)){stop("The inla model for your marker could not be run, see above warnings.")}
      
      
      if(Nsample>1){
        
        #samples 
        SMP <- inla.posterior.sample(Nsample-1, INLAmodel)
        linPred <- sapply(SMP, function(x) x$latent) 

        # keep only indice we want: 
        # Collapse each row into a string
        key1 <- do.call(paste, c(dataLongi_augmented[,colnames(dataLongi_augmented)%in%c(id,timeVar)], sep = "\r"))
        key2 <- do.call(paste, c(timePointsdata, sep = "\r"))
        
        # Find indices of rows from data1 that are in data2
        indices <- which(key1 %in% key2)
        PredYx<-linPred[indices,]
        
        #add BLUP first column
        PredYx<-cbind(INLAmodel$summary.linear.predictor$mean[indices],PredYx)
        #add informations 
        Outcome<-all.vars(terms(formLong[[indice]]))[1]
        PredYx<-cbind(timePointsdata,Outcome,PredYx)
        colnames(PredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
        
        
        
        
        
      }else{
        
        # keep only indice we want: 
        # Collapse each row into a string
        key1 <- do.call(paste, c(dataLongi_augmented[,colnames(dataLongi_augmented)%in%c(id,timeVar)], sep = "\r"))
        key2 <- do.call(paste, c(timePointsdata, sep = "\r"))
        
        # Find indices of rows from data1 that are in data2
        indices <- which(key1 %in% key2)
        #add BLUP first column
        #add informations 
        Outcome<-all.vars(terms(formLong[[indice]]))[1]
        PredYx<-cbind(timePointsdata,Outcome,INLAmodel$summary.linear.predictor$mean[indices])
        colnames(PredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
        
      }

      Yall[[indice]]<-PredYx
      
  
      }

  Yall<-do.call(rbind,Yall)
  print("End of running joint univarite models")
  return(Yall[,colnames(Yall)%in%c(id,timeVar,"Outcome",paste0("Sample_",c(1:Nsample)))])
  
}
