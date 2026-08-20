### Code:
##' @title Calculate predictions for time-depend covariates using INLA
##' @param objectY  A predYidm object from HIDeM package containing the univariate joint model with competing risk from step 1
##' @param objectSurvival A DYNidm object from HIDeM package containing the illness-death model estimation with time-dependent covariates 
##' @param newdata The newdata on which we want to do survival predictions 
##' @param s entry time 
##' @param horizon horizon time 
##' @param envir working environment 
#' @importFrom foreach "%do%"
#' @importFrom foreach "%dopar%"
#' @importFrom Deriv "Deriv"
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  
#' @export


DYNpredIDM<-function(objectY,
                          objectSurvival,
                          newdata,s=NULL,
                     horizon=NULL,
                          envir=parent.frame(),
                         predicted.newdata=NULL,
                     control=list(NsampleHY=1,NsampleFE=1,NsampleRE=1,return.data=F)){

  call <- match.call()
  ptm <- proc.time()

  if(missing(objectY)){stop("Need to specify objectY as a predYidm object")}
  if(!inherits(objectY,"predYidm")){stop("Need to specify objectY as a predYidm object")}
  if(missing(objectSurvival)){stop("Need to specify objectSurvival as a regDYNidm or DYNidm object")}
  if(!inherits(objectSurvival,"DYNidm")){stop("Need to specify objectSurvival as a DYNidm object")}

  
  if(missing(newdata)){stop("Need to provide newData as a data.frame")}
  if(!inherits(newdata,"data.frame")){stop("Need to provide newData as a data.frame")}
  
  if(sum(is.na(newdata))>0)stop("Need a new data frame with no missing data.")
  if(missing(s)|missing(horizon))stop("Need to specify s and horizon")
  
  if(!inherits(s,c("numeric","integer")))stop("s need to be an integer or numeric")
    if(length(s)!=1)stop("Length of s need to be 1")
    if((s < 0))stop("s need to be numeric superior or equal to 0 with s < horizon")
  
  
  if(!inherits(horizon,c("numeric","integer")))stop("horizon need to be an integer or numeric")
  if(length(horizon)!=1)stop("Length of horizon need to be 1")
  if(is.null(s))stop("landmark time need to be provided")
  if((s < 0) | (horizon < 0) | (s >= horizon))stop("s and horizon need to be numeric superior or equal to 0 with s < horizon")
  
  
  timeVar<-objectY$timeVar
  id<-objectY$id
  
  variables<-c(id,
               timeVar,
               objectSurvival$Xnames01,
               objectSurvival$Xnames02,
               objectSurvival$Xnames12,
               objectSurvival$timedepXnames01,
               objectSurvival$timedepXnames02,
               objectSurvival$timedepXnames12)
  variables<-unique(variables)
  
  if(any((variables%in%colnames(newdata))==F)){stop(paste0("Need a new data frame with no missing data on the variables : ",paste0(variables,collapse = ", ")))}
  
  if(objectY$method!="INLA"){stop("Prediction not available with JMBayes 2 method")}
  
  ############################## supress subjects censored or having the event before s ###
  #iderase<-newdata[which(newdata[,colnames(newdata])]
  ############################## supress info after s #########################
  
  newdata<-newdata[newdata[,colnames(newdata)%in%timeVar]<=s,]
  newdata<-newdata[,colnames(newdata)%in%variables]
  newdata<-na.omit(newdata)
  if(dim(newdata)[1]==0){stop("No follow-up after s, should proceed with smaller values for s")}
  #erase subjects having the event before s 
  N<-length(unique(newdata[,colnames(newdata)%in%objectY$id]))
  
  
  ######################## prepare for survival ################################
  
  #################################################################################
  ####################  prepare censored event times  #############################
  #################################################################################
  

  horizon1<-s1<-rep(NA,N)
    
  m01 <- model.frame(objectSurvival$formula01, data = newdata)
  m02 <- model.frame(objectSurvival$formula01, data = newdata)
 
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
  
  
  N <- length(abstime)
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
  t4<-rep(NA,N)
  ctime<-rep(NA,N)
  
  
  ctime<-ifelse(idm==0 & idd==0 & t1==t3,1,NA)
  ctime<-ifelse(idm==1 & idd==0 & t1<t2,2,ctime)
  ctime<-ifelse(idm==1 & idd==0 & t1==t2,3,ctime)
  ctime<-ifelse(idm==1 & idd==1 & t1<t2,4,ctime)
  ctime<-ifelse(idm==1 & idd==1 & t1==t2,5,ctime)
  ctime<-ifelse(idm==0 & idd==0 & t1<t3,6,ctime)
  ctime<-ifelse(idm==0 & idd==1,7,ctime)
  
  if(sum(is.na(ctime))>0){stop("For subject with no event, time for event 01 cannot be equal to time for 12 and 02")}
  
  t2<-ifelse(ctime==1 | ctime==3 | ctime==5,t1,
             ifelse(ctime==2 | ctime==4,t2,
                    ifelse(ctime==6 | ctime==7, t3,NA)
             )
  )
  
  t3<-ifelse(ctime==1, t1,
             ifelse(ctime==2 | ctime==3 | ctime==4 | ctime==5 | ctime==6 | ctime==7,t3,NA))
  
  t4<-ifelse(ctime==1 | ctime==2 | ctime==3 | ctime==4 | ctime==5, t1,
             ifelse( ctime==6 | ctime==7,t3,NA))
  
  ######################### association check RE ###############################
  assoc<-objectSurvival$assoc
  verifRE<-0
  
  for(k in 1:length(assoc)){
      vv<-na.omit(unlist(assoc[[k]]))
      if(sum(vv%in%"RE")==length(vv)){verifRE<-verifRE+1}
    }
    
  ####################### time fixed var #######################################
  x01<-newdata[,colnames(newdata)%in%c(objectSurvival$Xnames01,id)]
  x01<-unique(newdata)
  x01<-x01[,colnames(x01)%in%c(objectSurvival$Xnames01)]
  NC01<-NCOL(x01)
  
  x02<-newdata[,colnames(newdata)%in%c(objectSurvival$Xnames02,id)]
  x02<-unique(newdata)
  x02<-x02[,colnames(x02)%in%c(objectSurvival$Xnames02)]
  NC02<-NCOL(x02)
  
  
  x12<-newdata[,colnames(newdata)%in%c(objectSurvival$Xnames12,id)]
  x12<-unique(newdata)
  x12<-x12[,colnames(x12)%in%c(objectSurvival$Xnames12)]
  NC12<-NCOL(x12)
  
  size1 <- NC01+NC02+NC12
  
  noVar<-c(ifelse(as.integer(NC01)>0,0,1),
           ifelse(as.integer(NC02)>0,0,1),
           ifelse(as.integer(NC12)>0,0,1))
  
  nvat01 <- ifelse(noVar[1]==1,0,NC01)
  nvat02 <- ifelse(noVar[2]==1,0,NC02)
  nvat12 <- ifelse(noVar[3]==1,0,NC12)
  #nvat12<-length(unique(objectSurvival$Xnames12))
  
  dimnva01<-ifelse(nvat01==0,1,nvat01)
  dimnva02<-ifelse(nvat02==0,1,nvat02)
  dimnva12<-ifelse(nvat12==0,1,nvat12)
  
  NC<-c(NC01,NC02,NC12)
  
  if(noVar[1]==1){ve01<-as.double(rep(0,N))}else{ve01<-as.double(x01)}
  if(noVar[2]==1){ve02<-as.double(rep(0,N))}else{ve02<-as.double(x02)}
  if(noVar[3]==1){ve12<-as.double(rep(0,N))}else{ve12<-as.double(x12)}
  
  t0<-rep(s,N)
  t1<-rep(horizon,N)
  
  ########################## time dep var ######################################
  outcome01<-objectSurvival$linktimedepXnames01
  outcome02<-objectSurvival$linktimedepXnames02
  outcome12<-objectSurvival$linktimedepXnames12
  
  p01<-length(outcome01)
  dimp01<-ifelse(length(outcome01)>0,length(outcome01),1)
  p02<-length(outcome02)
  dimp02<-ifelse(length(outcome02)>0,length(outcome02),1)
  p12<-length(outcome12)
  dimp12<-ifelse(length(outcome12)>0,length(outcome12),1)
  if(p01==0 & p02==0){stop("No time dependent variable please refer to HIDeM::intensity for calculation of intensities with only time fixed covariates")}
  
 
  if(objectSurvival$method=="splines"){
    
   nknots01 <- length(unique(objectSurvival$knots01))
   nknots02<- length(unique(objectSurvival$knots02))
   amin<-min(newdata[,colnames(newdata)%in%timeVar])
   amax<-max(newdata[,colnames(newdata)%in%timeVar])
   
   if(knots01[1]> amin) stop(paste("Transition 0->1: Smallest time point should not appear before the starting knots at :",knots01[1]))
   if (knots01[length(knots01)] < amax) stop(paste("Transition 0->1: Highest time point should not appear after the last knots at :",knots01[length(knots01)]))
   
   if(knots02[1]> amin) stop(paste("Transition 0->2: Smallest time point should not appear before the starting knots at :",knots02[1]))
   if (knots02[length(knots02)] < amax) stop(paste("Transition 0->2: Highest time point should not appear after the last knots at :",knots02[length(knots02)]))
   if(knots12[1]> amin) stop(paste("Transition 1->2: Smallest time point should not appear before the starting knots at :",knots12[1]))
   if (knots12[length(knots12)] < amax) stop(paste("Transition 1->2: Highest time point should not appear after the last knots at :",knots12[length(knots12)]))
   
      
    size_spline<-nknots01+nknots02+nknots12+6
    
    
    # index_size1<-c(1:(nknots01+nknots02+nknots12+6))
    # start<-nknots01+nknots02+nknots12+6
    # if(nvat01>0|nvat02>0){
    # index_size1<-c(index_size1,((start+1):(nvat01+nvat02+start)))
    # }
    # start<-nvat01+nvat02+nvat12+start

  }else{ 
    
    # index_size1<-c(1:6)
    # if(!is.null(objectSurvival$posfix)){
    # nstart<-sum(objectSurvival$posfix%in%c(1:6))
    # }else{
    #   nstart<-7
    # }
    
    # if(nvat01>0|nvat02>0){
    # index_size1<-c(index_size1,(nstart:(nvat01+nvat02+nstart-1)))
    # }
    # start<-nvat01+nvat02+nvat12+6
    
    size_spline<-6 
    }
  

  size1<-size1+size_spline
  
  ### index to keep with binit
  # jump x12 and end after p02
  #index<-c(index_size1,(start+1):(start+p01+p02))
  
  size_V<-size1
  ############################################################################
  #################### defines initiate values ###############################
  ############################################################################
  
  # if length of regDYNidm or DYNidm superior to 1 need to summaries over replicates 
  if(length(objectSurvival$DYNidm)!=1){
    istop<-lapply(objectSurvival$DYNidm,FUN = function(x){
      if(x$istop%in%c(1,3)){return(T)}else{return(F)}
    })
    istop<-do.call(c,istop)
    if(sum(istop==F)==length(istop)){stop("All the survival models did not converged")}
    #binit<-prepareData(object=objectSurvival,istop=istop,index=index)
    binit<-prepareData(object=objectSurvival,istop=istop)
  }else{
    if(!objectSurvival$DYNidm[[1]]$istop%in%c(1,3)){stop("The survival model did not converged")}
    
    b_all<-rep(NA,length(objectSurvival$DYNidm[[1]]$b)+length(objectSurvival$posfix))
    if(!is.null(objectSurvival$posfix)){
    b_all[objectSurvival$posfix]<-objectSurvival$bfix
    b_all[-objectSurvival$posfix]<-objectSurvival$DYNidm[[1]]$b
    }else{
      b_all<-objectSurvival$DYNidm[[1]]$b
    }
    #binit<-b_all[index]
    binit<-b_all
  }

  ############################ keep index of variables selected ################
  
  varY<-unique(c(objectSurvival$timedepXnames01,
           objectSurvival$timedepXnames02,
           objectSurvival$timedepXnames12))
  Yindex<-lapply(objectY$formLong,FUN=function(x){as.character(x[[2]])})
  Yindex<-do.call(c,Yindex)
  index<-which(Yindex%in%varY)

  ############################ perform prediction of Y #########################
  
  if(is.null(predicted.newdata)){
  dataY<-DYNINLAidmpredY(object=objectY$modelY,
                         newdata=newdata,
                         s=s,
                         horizon=horizon,
                         scale.X=objectSurvival$scale.X,
                         assoc=objectY$assoc,
                         assocSurv=objectSurvival$assoc,
                         id=id,
                         timeVar=timeVar,
                         formLong=objectY$formLong,
                         formSurv=objectY$formSurv,
                         basRisk=objectY$basRisk,
                         index=index,
                         family=objectY$family,
                         envir=envir,
                         NsampleRE=control$NsampleRE,
                         NsampleHY=control$NsampleHY,
                         NsampleFE=control$NsampleFE)
  
  dataY$Outcome<-as.character(dataY$Outcome)
  # attention if NtimePoints equidistant with INLA then NtimePoints takes 
  # need ID to be numeric -- then 
  dataY[,colnames(dataY)%in%id]<-as.numeric(dataY[,colnames(dataY)%in%id])
  # to keep tracks of time order for each individual 
  dataY$order<-as.numeric(ave(dataY[,colnames(dataY)%in%id], cumsum(c(TRUE, diff(dataY[,colnames(dataY)%in%id]) != 0)), FUN = seq_along))
  
  
  
  }else{
    dataY<-predicted.newdata
    dataY<-data.frame(dataY)
  }
  
  ########################## check prediction ##################################
  NtimePoints<-271
  
  for( m in unique(c(outcome01,outcome02,outcome12))){
    subdata<-dataY[dataY$Outcome==m,]
    x<-table(subdata[,colnames(subdata)%in%id])
    if(any(x!=NtimePoints)){stop("Prediction of marker ",m," could not be perform for each quadrature points")}
    
  }
  

  
  
 
 
  if(length(outcome01)>=1){
    y01<-dataY[dataY$Outcome%in%outcome01,]
    # order  by individual and timeline 
    y01<-y01[order(y01[,colnames(y01)%in%id],y01$order),4]
    
    
  }else{
    y01<-rep(0,N*NtimePoints)
  }
  
  if(length(outcome02)>=1){
    y02<-dataY[dataY$Outcome%in%outcome02,]
    # order  by individual and timeline 
    y02<-y02[order(y02[,colnames(y02)%in%id],y02$order),4]
    
  }else{
    y02<-rep(0,N*NtimePoints)
  }
  
  if(length(outcome12)>=1){
    y12<-dataY[dataY$Outcome%in%outcome12,]
    # order  by individual and timeline 
    y12<-y12[order(y12[,colnames(y12)%in%id],y12$order),4]
    
  }else{
    y12<-rep(0,N*NtimePoints)
  }

  if(objectSurvival$method=="splines"){
    res<-rep(0,N)
    out1<- tryCatch({  .Fortran("citimedep",
                             ## input
                             as.double(binit),
                             as.integer(size_V),
                             as.double(knots01),
                             as.double(knots02),
                             as.double(knots12),
                             as.integer(N),
                             as.integer(nknots01),
                             as.integer(nknots02),
                             as.integer(nknots12),
                             as.double(ve01),
                             as.double(ve02),
                             as.double(ve12),
                             as.double(y01),
                             as.double(y02),
                             as.double(y12),
                             as.integer(p01),
                             as.integer(p02),
                             as.integer(p12),
                             as.integer(dimp01),
                             as.integer(dimp02),
                             as.integer(dimp12),
                             as.integer(NtimePoints),
                             as.integer(dimnva01),
                             as.integer(dimnva02),
                             as.integer(dimnva12),
                             as.integer(nvat01),
                             as.integer(nvat02),
                             as.integer(nvat12),
                             as.double(t0),
                             as.double(t1),
                             likelihood_res=as.double(res),
                             PACKAGE="HIDeM")$likelihood_res
  }, error = function(e) {
    # Return NULL on error to skip this patient
    NULL
  })
  }else{
    res<-rep(0,N)
    out1<- tryCatch({   .Fortran("ciweibtimedep",
                                ## input
                                as.double(binit),
                                as.integer(size_V),
                                as.integer(N),
                                as.double(ve01),
                                as.double(ve02),
                                as.double(ve12),
                                as.double(y01),
                                as.double(y02),
                                as.double(y12),
                                as.integer(p01),
                                as.integer(p02),
                                as.integer(p12),
                                as.integer(dimp01),
                                as.integer(dimp02),
                                as.integer(dimp12),
                                as.integer(NtimePoints),
                                as.integer(dimnva01),
                                as.integer(dimnva02),
                                as.integer(dimnva12),
                                as.integer(nvat01),
                                as.integer(nvat02),
                                as.integer(nvat12),
                                as.double(t0),
                                as.double(t1),
                                likelihood_res=as.double(res),
                                PACKAGE="HIDeM")$likelihood_res
    }, error = function(e) {
      # Return NULL on error to skip this patient
      NULL
    })
  }
  
  if(control$return.data==T){
    predicted.newdata<-dataY
  }else{predicted.newdata<-NULL}
  
  res<-list(CIF01=ifelse(!is.null(out1),out1[1:N],NULL),
            p00=ifelse(!is.null(out1),out1[(N+1):(2*N)],NULL),
            p02_0=ifelse(!is.null(out1),out1[(2*N+1):(3*N)],NULL),
            p12=ifelse(!is.null(out1),out1[(2*N+1):(3*N)],NULL),
            s=s,
            horizon=horizon,
            predicted.newdata= predicted.newdata)
  return(res)
 
}


# prepareData<-function(object,istop,index){
#   
# 
#   posfix<-object$posfix
#   bfix<-object$bfix
#   object<-object$DYNidm
#   Nrep<-sum(istop==T)
#   
#   b<-lapply(object,FUN=function(x){
#     if(x$istop%in%c(1,3)){
#       b_all<-rep(NA,length(x$b)+length(posfix))
#       if(!is.null(posfix)){
#       b_all[posfix]<-bfix
#       b_all[-posfix]<-x$b
#       }else{ b_all<-x$b}
#       return(b_all[index])}else{return(rep(NA,length(b_all[index])))}})
#   b<-do.call(rbind,b)
#   b<-na.omit(b)
#   b<-colSums(b)/Nrep
#   return(b)
# }


prepareData<-function(object,istop){
  
  
  posfix<-object$posfix
  bfix<-object$bfix
  object<-object$DYNidm
  Nrep<-sum(istop==T)
  
  b<-lapply(object,FUN=function(x){
    if(x$istop%in%c(1,3)){
      b_all<-rep(NA,length(x$b)+length(posfix))
      if(!is.null(posfix)){
        b_all[posfix]<-bfix
        b_all[-posfix]<-x$b
      }else{ b_all<-x$b}
      return(b_all)}else{return(rep(NA,length(b_all)))}})
  b<-do.call(rbind,b)
  b<-na.omit(b)
  b<-colSums(b)/Nrep
  return(b)
}
