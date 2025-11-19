### Code:
##' @title Calculate predictions for time-depend covariates using INLA
#' @useDynLib HIDeM
#' @importFrom INLAjoint joint
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  


INLAidm<-function(timeVar,family,basRisk,assoc,
                  truncated,formLong,formSurv,dataSurv,dataLongi,id,
                  nproc,t0,t1,t2,t3,
                  idm,idd,clustertype,lightmode){
  
  # define timePoints of prediction : 

  idsubjects<-unique(dataSurv[,colnames(dataSurv)%in%id])
  

  
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
  
  
  

  
  ## augmentation of the data 
  colnames(timePointsdata)<-c(id,timeVar)
  dataLongi_augmented<-merge(timePointsdata,dataLongi,by=c(id,timeVar),all.x=T,all.y=T)
  rownames(dataLongi_augmented)<-NULL
  dataLongi_augmented<-dataLongi_augmented[order(dataLongi_augmented[,colnames(dataLongi_augmented)%in%id],
                                                 dataLongi_augmented[,colnames(dataLongi_augmented)%in%timeVar]),]
  
 
  print("Start running joint univarite models")
  
  modelY<-list()
  length(modelY)<-length(formLong)
  
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
      
      
      #start with run = False -- have structure for slopes
    
      
      print(paste0("For marker: ",formLong[[indice]][[2]]))
     
     
      #int.strategy="eb" in previous version later 4.5.1
      INLAmodel<-INLAjoint::joint(formSurv = formSurv,
                                       formLong = formLong[[indice]],
                                       dataLong = dataLongi_augmented, dataSurv=dataSurv, id = id, timeVar = timeVar,
                                       family = family[indice],
                                       basRisk = basRisk[indice], NbasRisk = 15, assoc = assoc[[indice]],
                                       control=list(int.strategy="eb"))
      
      #browser()
      if(lightmode==T){
        erase<-c(".args","marginals.random","dic","waic","mode","residuals",
                 "logfile","selection","internal.marginals.hyperpar",
                 "marginals.hyperpar","size.random","priors_used","marginals.fixed",
                 "internal.summary.hyperpar","joint.hyper","formLong",
                 "summary.hyperpar","formSurv","SurvInfo","cpu.intern","call",
                 "size.linear.predictor","version","mlik","gcpo","cpo","cpu.used",
                 "model.random","basRisk","survOutcome","REstruc","assoc","names.fixed",
                 "assoc_Names","famLongi","fixRE","timeVar","id","longOutcome","cureVar",
                 "dataSurv","dataLong","corRE","mat_k","range","control.link","corLong",
                 "NbasRisk","variant","run","nhyper","survFacChar","lonFacChar",
                 "marginals.fitted.values","marginals.linear.predictor",
                 "summary.fitted.values","po","all.hyper")
        # for BLUP need : "summary.random","summary.fixed"
        nameINLA<-names(INLAmodel)
        classINLA<-class(INLAmodel)
        INLAmodel<-INLAmodel[-which(nameINLA%in%erase)]
        class(INLAmodel)<-classINLA
        
      }
      

      modelY[[indice]]<-INLAmodel
    }
  print("End of running univarite models")

  return(modelY)
  }
  
      
