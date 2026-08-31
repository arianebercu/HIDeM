### Code:
##' @title Calculate predictions for time-depend covariates using INLA
##' @param pred A vector containing the prediction for the probability of having the illness between s and horizon time
##' @param objectSurvival A idm object from HIDeM package containing the illness-death model estimation with time-fixed covariates 
##' @param newdata The newdata on which we want to calculate the BS, note : it must be the same as the one calculated for the predictions and in the same order as the predictions 
##' @param s entry time 
##' @param horizon horizon time 
##' @param envir working environment 
#' @importFrom foreach "%do%"
#' @importFrom foreach "%dopar%"
#' @importFrom Deriv "Deriv"
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  
#' @export


BS<-function(pred,
             objectSurvival,
             s,
             horizon,
             newdata,
             envir=parent.frame())
  {


  call <- match.call()
  ptm <- proc.time()

  if(missing(pred)){stop("Need to specify a vector of predictions in pred")}
  if(!(class(pred)%in%"numeric")){stop("Need to specify a vector of numeric predictions in pred")}
  if(missing(objectSurvival)){stop("Need to specify objectSurvival as a idm object")}
  if(!inherits(objectSurvival,"idm")){stop("Need to specify objectSurvival as a idm object")}

  
  if(missing(newdata)){stop("Need to provide newdata as a data.frame")}
  if(!inherits(newdata,"data.frame")){stop("Need to provide newdata as a data.frame")}
  
  if(sum(is.na(newdata))>0)stop("Need a new data frame with no missing data.")
  if(missing(s)|missing(horizon))stop("Need to specify s and horizon")
  
  if(!inherits(s,c("numeric","integer")))stop("s need to be an integer or numeric")
    if(length(s)!=1)stop("Length of s need to be 1")
    if((s < 0))stop("s need to be numeric superior or equal to 0 with s < horizon")
  
  
  if(!inherits(horizon,c("numeric","integer")))stop("horizon need to be an integer or numeric")
  if(length(horizon)!=1)stop("Length of horizon need to be 1")
  if(is.null(s))stop("landmark time need to be provided")
  if((s < 0) | (horizon < 0) | (s >= horizon))stop("s and horizon need to be numeric superior or equal to 0 with s < horizon")
  

  if(dim(newdata)[1]!=length(pred)){stop("The number of lines in newdata needs to be the size of pred")}
  #erase subjects having the event before s 
  N<-length(unique(newdata[,colnames(newdata)%in%objectSurvival$id]))
  
  
  #################################################################################
  ####################  prepare censored event times  #############################
  #################################################################################
  
  newdataS <- newdata[!duplicated(newdata[,colnames(newdata)%in%objectSurvival$id], fromLast = TRUE), ]
  m01 <- model.frame(objectSurvival$formula01, data = newdataS)
  m02 <- model.frame(objectSurvival$formula02, data = newdataS)
  m12 <- model.frame(objectSurvival$formula12, data = newdataS)

  responseTrans <- stats::model.response(m01)
  responseAbs <- stats::model.response(m02)

  isIntervalCensored <- attr(responseTrans,"cens.type")=="intervalCensored"
  truncated <- nchar(attr(responseAbs,"entry.type"))>1
  abstime <- as.double(responseAbs[,"time"])
  ## It may happen that the illness time is observed exactly, in which case
  ## the status is 1, thus we need two criteria to declare illness status:
  ## 1. exact observations with illness status ==1
  ## 2. interval censored with any illness status. FIXME: check the corresponding likelihood


  idm <- responseTrans[,"status"]==(as.integer(isIntervalCensored)+1)
  if (isIntervalCensored)
    idm[(responseTrans[,"status"]==1 & (responseTrans[,"L"]==responseTrans[,"R"]))] <- 1
  ## exit status
  idd <- responseAbs[,"status"]==1


  #N <- length(abstime)
  if (truncated==0){
    entrytime <- as.double(NULL)
  }else{
    entrytime <- as.double(responseAbs[,"entry"])
  }
  if (isIntervalCensored){
    Ltime <- as.double(responseTrans[,"L",drop=TRUE])
    Rtime <- as.double(responseTrans[,"R",drop=TRUE])
    ## if (any(Rtime<abstime & idm ==0))
    ## warning(paste("For ",
    ## sum(Rtime<abstime & idm ==0),
    ## " cases where the ill status is not observed\n and the last inspection time (R) is smaller than the right censored time (T)\n the time R is set to T."))
  }else{# exactly observed transition times
    Ltime <- as.double(responseTrans[,"time",drop=TRUE])
    Rtime <- as.double(responseTrans[,"time",drop=TRUE])
    Ltime[idm==0] <- abstime[idm==0]
    Rtime[idm==0] <- abstime[idm==0]
  }
  ## find time boundaries
  if (length(entrytime)>0){
    alltimes <- sort(unique(c(Ltime, Rtime,entrytime,abstime)))
    amax <- max(alltimes)
    amin <- min(alltimes)
  }
  else{
    alltimes <- sort(unique(c(Ltime, Rtime,abstime)))
    amax <- max(alltimes)
    amin <- 0
  }


  t0<-rep(s,N)
  t1<-Ltime
  t2<-Rtime
  t3<-abstime
  t4<-rep(horizon,N)
  ctime<-rep(NA,N)
  if (isIntervalCensored){

  ctime<-ifelse(t1>horizon,0,NA)
  ctime<-ifelse((t1>s) & (t2<=horizon) & (idm==1),1,ctime)
  ctime<-ifelse((t1<=horizon) & (t1>s) & (t3>horizon) & (t2>horizon),2,ctime)
  ctime<-ifelse((idm==1) & (t1<s) & (t2>s) & (t2<=horizon) & (t3>horizon),3,ctime)
  ctime<-ifelse((idm==1) & (t1<s) & (t2>s) & (t2<=horizon) & (idd==1) & (t3<=horizon) & (t3>s),4,ctime)
  ctime<-ifelse((idd==1) & (t1<=horizon) & (t1>s) & (idm==0) & (t3<=horizon),5,ctime)
  ctime<-ifelse((idd==1) & (t1<=s) & (idm==0) & (t3<=horizon) & (t3>s),6,ctime)
  ctime<-ifelse((t1<=s) & (t2>=horizon) & (t3>=horizon),6,ctime)

  if(sum(is.na(ctime))>0){
    id<-which(is.na(ctime))
    id<-paste0(id,collapse=";")
    stop(paste0("For subject(s) : ",id," error of classification ctime, please check if the subjects are at risk between s and horizon"))
  }

  }else{ #t1=t2
    
    ctime<-ifelse(t1>horizon,0,NA)
    ctime<-ifelse((t1>s) & (t1<=horizon) & (idm==0) & (t3<horizon) & (idd==1),1,ctime)
    ctime<-ifelse((t1>s) & (t1<=horizon) & (idm==0) & (t3>horizon),2,ctime)
    ctime<-ifelse((t1>s) & (t1<=horizon) & (idm==1),3,ctime)
    ctime<-ifelse((t1<s) & (idm==0) & (idd==1) & (t3<horizon)& (t3>s),4,ctime)
    ctime<-ifelse((t1<s) & (idm==0) & (t3>horizon),5,ctime)
    
    if(sum(is.na(ctime))>0){
      id<-which(is.na(ctime))
      id<-paste0(id,collapse=";")
      stop(paste0("For subject(s) : ",id," error of classification ctime, please check if the subjects are at risk between s and horizon"))
    }
    
  }
  browser()
  ####################### calculate weights over subjects ######################
  if (isIntervalCensored){
  w1<-unlist(lapply(c(1:N),function(x){
    
    if(ctime[x]==0){return(0)}
    if(ctime[x]==1){return(1)}
    if(ctime[x]%in%c(2,7)){
      predtiming<-predict(objectSurvival,s=t1[x],t=horizon,conf.int=F,nsim=1)
      p01<-predtiming$p01
      
      predtiming<-predict(objectSurvival,s=0,t=horizon,conf.int=F,nsim=1)
      p00<-predtiming$p00
      
      if(k==2){
        return(p01/(p00+p01))
      }else{
        predtiming<-predict(objectSurvival,s=s,t=horizon,conf.int=F,nsim=1)
        p01_1<-predtiming$p01
        return(p01_1/(p00+p01))
      }
    }
    
    if(ctime[x]%in%c(3,4)){
      
      predtiming<-predict(objectSurvival,s=s,t=t2[x],conf.int=F,nsim=1)
      p01_1<-predtiming$p01
     
      
      predtiming<-predict(objectSurvival,s=t1[x],t=t2[x],conf.int=F,nsim=1)
      p01_2<-predtiming$p01
      
      return(p01_1/p01_2)
      
    }
    
    if(ctime[x]%in%c(5,6)){
      
      predtiming<-predict(objectSurvival,s=t1[x],t=t3[x],conf.int=F,nsim=1)
      p01<-predtiming$p01
      p02<-exp(-predtiming$cumulativeintensity01-predtiming$cumulativeintensity02)*predtiming$intensity02
   
      if(k==5){
        return((p01*predtiming$intensity12)/(p01*predtiming$intensity12+p02))
      }else{
        predtiming<-predict(objectSurvival,s=s,t=t3[x],conf.int=F,nsim=1)
        p01_1<-predtiming$p01
        return((p01_1*predtiming$intensity12)/(p01*predtiming$intensity12+p02))
      }
      
    }
    
  }))
  
  
  w0<-unlist(lapply(c(1:N),function(x){
    
    if(ctime[x]==0){return(1)}
    if(ctime[x]%in%c(1,3,4)){return(0)}
    if(ctime[x]%in%c(2,5)){
      return(1-w1[x])
    }
    
    if(ctime[x]==6){
      
      predtiming<-predict(objectSurvival,s=0,t=t3[x],conf.int=F,nsim=1)
      p02<-predtiming$p02_0
      
      
      predtiming<-predict(objectSurvival,s=t1[x],t=t3[x],conf.int=F,nsim=1)
      p01_2<-predtiming$p01*predtiming$intensity02
      
      return(p02/(p01_2+p02))
      
    }
    
    if(ctime[x]==7){
      
      predtiming<-predict(objectSurvival,s=0,t=horizon,conf.int=F,nsim=1)
      p00<-predtiming$p00
      
      
      predtiming<-predict(objectSurvival,s=t1[x],t=horizon,conf.int=F,nsim=1)
      p01<-predtiming$p01*predtiming$intensity02
      
      return(p00/(p01+p00))
             
    }
    

    
  }))

  }else{
    
    w1<-unlist(lapply(c(1:N),function(x){
      
      if(ctime[x]==0){return(0)}
      if(ctime[x]==3){return(1)}
      if(ctime[x]%in%c(2,5)){
        predtiming<-predict(objectSurvival,s=t1[x],t=horizon,conf.int=F,nsim=1)
        p01<-predtiming$p01
        
        predtiming<-predict(objectSurvival,s=0,t=horizon,conf.int=F,nsim=1)
        p00<-predtiming$p00
        
        if(k==2){
          return(p01/(p00+p01))
        }else{
          predtiming<-predict(objectSurvival,s=s,t=horizon,conf.int=F,nsim=1)
          p01_1<-predtiming$p01
          return(p01_1/(p00+p01))
        }
      }
      
      
      if(ctime[x]%in%c(1,4)){
        
        predtiming<-predict(objectSurvival,s=t1[x],t=t3[x],conf.int=F,nsim=1)
        p01<-predtiming$p01
        p02<-exp(-predtiming$cumulativeintensity01-predtiming$cumulativeintensity02)*predtiming$intensity02
        
        if(k==1){
          return((p01*predtiming$intensity12)/(p01*predtiming$intensity12+p02))
        }else{
          predtiming<-predict(objectSurvival,s=s,t=t3[x],conf.int=F,nsim=1)
          p01_1<-predtiming$p01
          return((p01_1*predtiming$intensity12)/(p01*predtiming$intensity12+p02))
        }
        
      }
      
    }))
    
    
    w0<-unlist(lapply(c(1:N),function(x){
      
      if(ctime[x]==0){return(1)}
      if(ctime[x]==3){return(0)}
      if(ctime[x]%in%c(1,2)){
        return(1-w1[x])
      }
      
      
      if(ctime[x]==4){
        
        predtiming<-predict(objectSurvival,s=t1[x],t=t3[x],conf.int=F,nsim=1)
        p01<-predtiming$p01
        p02<-exp(-predtiming$cumulativeintensity01-predtiming$cumulativeintensity02)*predtiming$intensity02
        
          return(p02/(p01*predtiming$intensity12+p02))
        
      }
      
      if(ctime[x]==5){
        
        predtiming<-predict(objectSurvival,s=t1[x],t=horizon,conf.int=F,nsim=1)
        p01<-predtiming$p01
        p00<-exp(-predtiming$cumulativeintensity01-predtiming$cumulativeintensity02)
        
        return(p00/(p01+p00))
        
      }
      
    }))
    
  }
  if(any(is.na(w1)))stop("Weights missing for illness cases")
  if(any(is.na(w0)))stop("Weights missing for illness-free cases")
  # calculate brier score 
    BS<-sum(w1*(1-pred)^2+w0*(0-pred)^2)
  res<-list(BS=BS,
            s=s,
            horizon=horizon)
  return(res)
 
}


