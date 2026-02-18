### Code:
##' @title Calculate predictions for time-depend covariates using INLA
#' @importFrom foreach "%do%"
#' @importFrom foreach "%dopar%"
#' @importFrom Deriv "Deriv"
#' @importFrom INLA "inla.posterior.sample"
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  


DYNpredIDM_INLA<-function(object,newdata,s,horizon,event, nb.draws=1000,
                     modelY,seed,BLUP,scale.X,id,timeVar,formLong){

  
  
  # define timePoints of prediction : 

  set.seed(seed)
  
  idsubjects<-unique(newdata[,colnames(newdata)%in%id])
  
  # force only observation before s or at s 
  newdata<-newdata[newdata[,colnames(newdata)%in%timeVar]<=s,]
  

  #create gauss kronrod time point between s and s+t
  timePointsdata<-do.call(rbind, lapply(idsubjects,FUN=function(x){
    timePoints<-gauss_kronrod_points(lower.intdouble=s,
                                     upper.intdouble=horizon,
                                     end.time=horizon,
                                     truncated=T,
                                     entry.time=s)
    return(data.frame(index=x,timePoints=timePoints))}))
  

  
  ## augmentation of the newdata 
  colnames(timePointsdata)<-c(id,timeVar)
  dataLongi_augmented<-merge(timePointsdata,newdata,by=c(id,timeVar),all.x=T,all.y=T)
  rownames(dataLongi_augmented)<-NULL
  dataLongi_augmented<-dataLongi_augmented[order(dataLongi_augmented[,colnames(dataLongi_augmented)%in%id],
                                                 dataLongi_augmented[,colnames(dataLongi_augmented)%in%timeVar]),]
  
  Yall<-list()
  length(Yall)<-length(formLong)
  
  
  
  formSurv<-list(inla.surv(time=TimeCR, event=idm) ~ -1,
                 inla.surv(time=TimeCR, event=iddCR) ~ -1)
  
  newdata$TimeCR<-NA
  newdata$idm<-NA
  newdata$iddCR<-NA
  #
#  if(scale.X==T){
#    tcenter<-ifelse(truncated==T,0,min(t0))
#    dataCenter<-data.frame(ID=idsubjects,time=tcenter)
#    colnames(dataCenter)<-c(id,timeVar)
#  }
  
  timepred<-unique(timePointsdata[,colnames(timePointsdata)%in%timeVar])
  NtimePoints<-length(unique(timepred))
  # should center by median or mean ? 
  
    
    for(indice in 1:length(formLong)){

      INLAmodel<-object$modelY[[indice]]
      
      INLAmodel$run<-F
      
      pp<-joint.run(model = INLAmodel,data= newdata,
                  formSurv = formSurv,
              id = id, timeVar = timeVar,
              family = family[indice],
              basRisk = basRisk[indice], NbasRisk = 15, assoc = assoc[[indice]],
              control=list(int.strategy="eb"),run=F)
      
      # data structure
      ct <- INLAmodel$misc$configs$contents
      
      if(is.null(INLAmodel)){stop("The inla model for your marker could not be run, see above warnings.")}
      choiceY<-na.omit(unlist(assoc[[indice]]))
      
      ids<-str_detect(rownames(INLAmodel$summary.hyperpar),"Theta")
      inla.hyperpar.sample(n=2,INLAmodel)
      tt<-inla.hyperpar(result=INLAmodel)
      
      
      Vparam<-INLAmodel$summary.hyperpar[ids,"mode"]
      p<-(-1+sqrt(1+8*length(Vparam)))/2
      matCov<-matrix(0,ncol=p,nrow=p)
      matCov[upper.tri(matCov,diag=T)]<-Vparam
      
      
      if(BLUP==F){
        
        
        #samples seed=seed cannot do parallel estimation on it 
        SMP <- INLA::inla.posterior.sample(1, INLAmodel,seed=seed)
        
        res<-NULL
        key1 <- do.call(paste, c(dataLongi_augmented[,colnames(dataLongi_augmented)%in%c(id,timeVar)], sep = "\r"))
        key2 <- do.call(paste, c(timePointsdata, sep = "\r"))
        # keep only indice we want: 
        # Collapse each row into a string
        
        # Find indices of rows from data1 that are in data2
        # while keeping order of key2
        indices <- match(key2, key1)
        
        if("value"%in% choiceY){
        
        Y<-as.matrix(make_XINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=SMP[[1]]))
        Y<-Y[indices,]
        Outcome<-all.vars(terms(formLong[[indice]]))[1]
        PredYx<-cbind(timePointsdata,Outcome=Outcome,Y)
        if(scale.X==T){
          Ycenter<-make_XINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=SMP[[1]])
          PredYx$Y<-(PredYx$Y-mean(Ycenter))/sd(Ycenter)
        }
        colnames(PredYx)[4]<-"Sample_1"
        res<-rbind(res,PredYx)
        }
        
        if("slope"%in% choiceY){
        dY<-as.matrix(make_dXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=SMP[[1]]))
        dY<-dY[indices,]
        slopeOutcome<-paste0("slope_",Outcome)
        slopePredYx<-cbind(timePointsdata,Outcome=slopeOutcome,dY)
        if(scale.X==T){
          dYcenter<-make_dXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=SMP[[1]])
          slopePredYx$dY<-(slopePredYx$dY-mean(dYcenter))/sd(dYcenter)
        }
        colnames(slopePredYx)[4]<-"Sample_1"
        res<-rbind(res,slopePredYx)
        }
        
        if("RE"%in% choiceY){
          
        REY<-as.matrix(make_REXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=SMP[[1]]))
       REY<-REY[indices,]
       
       if(scale.X==T){
         REYcenter<-as.matrix(make_REXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=SMP[[1]]))
         REY<-do.call(cbind,lapply(c(1:dim(REY)[2]),FUN=function(x){
           (REY[,x]-mean(REYcenter[,x]))/sd(REYcenter[,x])
         }))
       }
       namesREY<-unlist(lapply(1:dim(REY)[2],FUN=function(x){
         rep(paste0("RE_",colnames(REY)[x],"_",Outcome),dim(REY)[1])
       }))
       dataREY<- do.call(rbind, replicate(dim(REY)[2], timePointsdata, simplify = FALSE))
       REPredYx<-cbind(dataREY,Outcome=namesREY,as.vector(REY))
       
       
       colnames(REPredYx)[4]<-"Sample_1"
       res<-rbind(res,REPredYx)
        }
        
        Yall[[indice]]<- res
        
        
        
      }else{
        
        
        res<-NULL
        key1 <- do.call(paste, c(dataLongi_augmented[,colnames(dataLongi_augmented)%in%c(id,timeVar)], sep = "\r"))
        key2 <- do.call(paste, c(timePointsdata, sep = "\r"))
        # keep only indice we want: 
        # Collapse each row into a string
        
        # Find indices of rows from data1 that are in data2
        # while keeping order of key2
        indices <- match(key2, key1)
        
        if("value"%in% choiceY){
          Y<-as.matrix(make_XINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=INLAmodel))
          Y<-Y[indices,]
          Outcome<-all.vars(terms(formLong[[indice]]))[1]
          PredYx<-cbind(timePointsdata,Outcome=Outcome,Y)
          
          if(scale.X==T){
            Ycenter<-make_XINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=INLAmodel)
            PredYx$Y<-(PredYx$Y-mean(Ycenter))/sd(Ycenter)
          }
          colnames(PredYx)[4]<-"Sample_1"
          res<-rbind(res,PredYx)
        }
        
        if("slope"%in% choiceY){
        dY<-as.matrix(make_dXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=INLAmodel))
        dY<-dY[indices,]
        slopeOutcome<-paste0("slope_",Outcome)
        slopePredYx<-cbind(timePointsdata,Outcome=slopeOutcome,dY)
        
        if(scale.X==T){
          dYcenter<-make_dXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=INLAmodel)
          slopePredYx$dY<-(slopePredYx$dY-mean(dYcenter))/sd(dYcenter)
        }
        colnames(slopePredYx)[4]<-"Sample_1"
        res<-rbind(res,slopePredYx)
        }
        
        if("RE"%in% choiceY){
        REY<-as.matrix(make_REXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=INLAmodel))
        REY<-REY[indices,]
        
        if(scale.X==T){
          REYcenter<-as.matrix(make_REXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=INLAmodel))
          REY<-do.call(cbind,lapply(c(1:dim(REY)[2]),FUN=function(x){
            (REY[,x]-mean(REYcenter[,x]))/sd(REYcenter[,x])
          }))
        }
        namesREY<-unlist(lapply(1:dim(REY)[2],FUN=function(x){
          rep(paste0("RE_",colnames(REY)[x],"_",Outcome),dim(REY)[1])
        }))
        dataREY<- do.call(rbind, replicate(dim(REY)[2], timePointsdata, simplify = FALSE))
        REPredYx<-cbind(dataREY,Outcome=namesREY,as.vector(REY))
        colnames(REPredYx)[4]<-"Sample_1"
        res<-rbind(res,REPredYx)
        }
        
        Yall[[indice]]<- res
        
        
        
      }

      
  
    }
      
      
      


  Yall<-do.call(rbind,Yall)
  
  
  return(Yall[,colnames(Yall)%in%c(id,timeVar,"Outcome","Sample_1")])
  

  
}


