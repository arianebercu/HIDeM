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
                          newdata,
                     s=NULL,
                     horizon=NULL,
                     times,
                          envir=parent.frame(),
                         predicted.newdata=NULL,
                     control=list(NsampleHY=1,NsampleFE=1,NsampleRE=1,NidLoop="auto",return.data=F)){


  call <- match.call()
  ptm <- proc.time()

  if(missing(objectY)){stop("Need to specify objectY as a predYidm object")}
  if(!inherits(objectY,"predYidm")){stop("Need to specify objectY as a predYidm object")}
  if(missing(objectSurvival)){stop("Need to specify objectSurvival as a regDYNidm or DYNidm object")}
  if(!inherits(objectSurvival,"DYNidm")){stop("Need to specify objectSurvival as a DYNidm object")}

  
  if(missing(newdata)){stop("Need to provide newdata as a data.frame")}
  if(!inherits(newdata,"data.frame")){stop("Need to provide newdata as a data.frame")}
  if(!missing(times)){
    if(!inherits(times,"list")){stop("times must be a list of 9 elements")}
    if(length(times)!=9){stop("times must be a list of 9 elements")}
    
    expected <- c("start", "end", "id", "p00", "p01", "p02_0","p02_1","p11","p12")
    if (!identical(names(times), expected)) {
      stop(
        "Invalid names in list.\n",
        "Expected: ", paste(expected, collapse = ", "), "\n",
        "Found:    ", paste(names(times), collapse = ", ")
      )
    }
    if(!inherits(times$p00,"logical")){stop("times$p00 must be a logical")}
    if(!inherits(times$p01,"logical")){stop("times$p01 must be a logical")}
    if(!inherits(times$p02_0,"logical")){stop("times$p00 must be a logical")}
    if(!inherits(times$p02_1,"logical")){stop("times$p00 must be a logical")}
    if(!inherits(times$p11,"logical")){stop("times$p11 must be a logical")}
    if(!inherits(times$p12,"logical")){stop("times$p12 must be a logical")}
    
    if(any((times$id%in%newdata[,objectY$id])==F)){stop("All id in times must be in newdata")}
    
    lens <- lengths(times)
    if (length(unique(lens)) != 1) {
      stop(
        "All elements must have the same length.\n",
        "Lengths found: ", paste(names(times), lens, sep = " = ", collapse = ", ")
      )
    }
    
    if(!inherits(times$start,c("numeric","integer")))stop("times$start need to be an integer or numeric")
    if(any(times$start < 0))stop("times$start need to be numeric superior or equal to 0 with s < horizon")
    if(!inherits(times$end,c("numeric","integer")))stop("times$end need to be an integer or numeric")
    if(any(times$end < 0))stop("times$end need to be numeric superior or equal to 0 with s < horizon")
    
  }else{
    times<-NULL
  }
  
  if(sum(is.na(newdata))>0)stop("Need a new data frame with no missing data.")
  if(missing(s)|missing(horizon))stop("Need to specify s and horizon")
  
  if(!inherits(s,c("numeric","integer")))stop("s need to be an integer or numeric")
    if(length(s)!=1)stop("Length of s need to be 1")
    if((s < 0))stop("s need to be numeric superior or equal to 0 with s < horizon")
  
  
  if(!inherits(horizon,c("numeric","integer")))stop("horizon need to be an integer or numeric")
  
  if(!inherits(control$NsampleHY,c("numeric","integer"))){
    if(round(control$NsampleHY)!=control$NsampleHY)stop("control$NsampleHY need to be an integer")}
  if(!inherits(control$NsampleFE,c("numeric","integer"))){
    if(round(control$NsampleFE)!=control$NsampleFE)stop("control$NsampleFE need to be an integer")}
  if(!inherits(control$NsampleRE,c("numeric","integer"))){
    if(round(control$NsampleRE)!=control$NsampleRE)stop("control$NsampleHYneed to be an integer")}

  if(!inherits(control$NidLoop,c("numeric","integer","character")))stop("control$NidLoop to be an integer or character")
if(length(control$NsampleHY)!=1)stop("Length of control$NsampleHY need to be 1")
  if(length(control$NsampleFE)!=1)stop("Length of control$NsampleFE need to be 1")
  if(length(control$NsampleRE)!=1)stop("Length of control$NsampleRE need to be 1")
  if(length(control$NidLoop)!=1)stop("Length of control$NidLoop need to be 1")
  if(inherits(control$NidLoop,"character")){
    if(control$NidLoop!="auto"){
    stop("control$NidLoop to be auto")}
  }else{
    if(round(control$NidLoop)!=control$NidLoop)stop("control$NidLoop to be an integer")
    if(control$NidLoop<=0)stop("control$NidLoop needs to be greater than 0")
  }
  if(control$NsampleFE<=0)stop("control$NsampleFE needs to be greater than 0")
  if(control$NsampleRE<=0)stop("control$NsampleRE needs to be greater than 0")
  if(control$NsampleHY<=0)stop("control$NsampleHY needs to be greater than 0")
  
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
  #newdata<-newdata[,colnames(newdata)%in%variables]
  newdata<-na.omit(newdata)
  if(dim(newdata)[1]==0){stop("No follow-up after s, should proceed with smaller values for s")}
  #erase subjects having the event before s 
  N<-length(unique(newdata[,colnames(newdata)%in%objectY$id]))
  
  if(inherits(control$NidLoop,c("numeric","integer"))){
  if(control$NidLoop>N)stop(paste0("control$NidLoop needs to less or equal to ",N))}
  ######################## prepare for survival ################################
  
  #################################################################################
  ####################  prepare censored event times  #############################
  #################################################################################
  
  #browser()
  # 31/08/2026 : do not estimate weights with our model 
  # horizon1<-s1<-rep(NA,N)
  #   
  # newdataS <- newdata[!duplicated(newdata[,colnames(newdata)%in%objectY$id], fromLast = TRUE), ]
  # m01 <- model.frame(objectSurvival$formula01, data = newdataS)
  # m02 <- model.frame(objectSurvival$formula02, data = newdataS)
  # 
  # responseTrans <- stats::model.response(m01)
  # responseAbs <- stats::model.response(m02)
  # 
  # isIntervalCensored <- attr(responseTrans,"cens.type")=="intervalCensored"
  # truncated <- nchar(attr(responseAbs,"entry.type"))>1
  # abstime <- as.double(responseAbs[,"time"])
  # ## It may happen that the illness time is observed exactly, in which case
  # ## the status is 1, thus we need two criteria to declare illness status:
  # ## 1. exact observations with illness status ==1
  # ## 2. interval censored with any illness status. FIXME: check the corresponding likelihood
  # 
  # 
  # idm <- responseTrans[,"status"]==(as.integer(isIntervalCensored)+1)
  # if (isIntervalCensored)
  #   idm[(responseTrans[,"status"]==1 & (responseTrans[,"L"]==responseTrans[,"R"]))] <- 1
  # ## exit status
  # idd <- responseAbs[,"status"]==1
  # 
  # 
  # #N <- length(abstime)
  # if (truncated==0){
  #   entrytime <- as.double(NULL)
  # }else{
  #   entrytime <- as.double(responseAbs[,"entry"])
  # }
  # if (isIntervalCensored){
  #   Ltime <- as.double(responseTrans[,"L",drop=TRUE])
  #   Rtime <- as.double(responseTrans[,"R",drop=TRUE])
  #   ## if (any(Rtime<abstime & idm ==0))
  #   ## warning(paste("For ",
  #   ## sum(Rtime<abstime & idm ==0),
  #   ## " cases where the ill status is not observed\n and the last inspection time (R) is smaller than the right censored time (T)\n the time R is set to T."))
  # }else{# exactly observed transition times
  #   Ltime <- as.double(responseTrans[,"time",drop=TRUE])
  #   Rtime <- as.double(responseTrans[,"time",drop=TRUE])
  #   Ltime[idm==0] <- abstime[idm==0]
  #   Rtime[idm==0] <- abstime[idm==0]
  # }
  # ## find time boundaries 
  # if (length(entrytime)>0){
  #   alltimes <- sort(unique(c(Ltime, Rtime,entrytime,abstime)))
  #   amax <- max(alltimes)
  #   amin <- min(alltimes)
  # }
  # else{
  #   alltimes <- sort(unique(c(Ltime, Rtime,abstime)))
  #   amax <- max(alltimes)
  #   amin <- 0
  # }
  # 
  # 
  # t0<-rep(s,N)
  # t1<-Ltime
  # t2<-Rtime
  # t3<-abstime
  # t4<-rep(horizon,N)
  # ctime<-rep(NA,N)
  # if (isIntervalCensored){
  # 
  # 
  # 
  # ctime<-ifelse(t1>horizon,0,NA)
  # ctime<-ifelse((t1>s) & (t2<=horizon) & (idm==1),1,ctime)
  # ctime<-ifelse((t1<=horizon) & (t1>s) & (t3>horizon),2,ctime)
  # ctime<-ifelse((idm==1) & (t1<s) & (t2>s) & (t2<=horizon),3,ctime)
  # ctime<-ifelse((idd==1) & (t1<=horizon) & (t1>s) & (idm==0) & (t3<=horizon),4,ctime)
  # ctime<-ifelse((idd==1) & (t1<=s) & (idm==0) & (t3<=horizon) & (t3>s),5,ctime)
  # ctime<-ifelse((t1<=s) & (idm==1) & (t2>=horizon) & (t3>=horizon),6,ctime)
  # 
  # # recode times l and r 
  # # case 0
  # t1[which(ctime==0)]<-s
  # t2[which(ctime==0)]<-horizon
  # #case 1 no changes 
  # #case 2 
  # t2[which(ctime==2)]<-horizon
  # # case 3
  # t1[which(ctime==3)]<-s
  # #case 4
  # t2[which(ctime==4)]<-t3[which(ctime==4)]
  # #case 5 
  # t1[which(ctime==5)]<-s
  # t2[which(ctime==5)]<-t3[which(ctime==5)]
  # 
  # #case 6
  # t1[which(ctime==6)]<-s
  # t2[which(ctime==6)]<-horizon
  # 
  # if(sum(is.na(ctime))>0){stop("For subject with no event, time for event 01 cannot be equal to time for 12 and 02")}
  # 
  # }
  
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
  
  if(!is.null(times)){
    
    get_ids <- function(col) if (!is.null(times)) times$id[times[[col]] %in% TRUE] else integer(0)
    ids00   <- get_ids("p00")
    ids01   <- get_ids("p01")
    ids02_0 <- get_ids("p02_0")
    ids02_1 <- get_ids("p02_1")
    ids11 <- get_ids("p11")
    ids12 <- get_ids("p12")
    
    
  }

  if(noVar[1]==1){
    ve01<-as.double(rep(0,N))
    if(!is.null(times)){
    ve01_00<-as.double(rep(0,length(ids00)))
    ve01_02_0<-as.double(rep(0,length(ids02_0)))
    ve01_02_1<-as.double(rep(0,length(ids02_1)))
    ve01_01<-as.double(rep(0,length(ids01)))
    ve01_11<-as.double(rep(0,length(ids11)))
    ve01_12<-as.double(rep(0,length(ids12)))
    }
  }else{
      ve01<-as.double(x01)
      if (!is.null(times)) {
        Xn01 <- objectSurvival$Xnames01
        ve01_00   <- get_ve(newdata, id, Xn, ids00)
        ve01_01   <- get_ve(newdata, id, Xn, ids01)
        ve01_02_0 <- get_ve(newdata, id, Xn, ids02_0)
        ve01_02_1 <- get_ve(newdata, id, Xn, ids02_1)
        ve01_11 <- get_ve(newdata, id, Xn, ids11)
        ve01_12 <- get_ve(newdata, id, Xn, ids12)
      }
      
      }
  if(noVar[2]==1){
    ve02<-as.double(rep(0,N))
    if(!is.null(times)){
      ve02_00<-as.double(rep(0,length(ids00)))
      ve02_02_0<-as.double(rep(0,length(ids02_0)))
      ve02_02_1<-as.double(rep(0,length(ids02_1)))
      ve02_01<-as.double(rep(0,length(ids01)))
      ve02_11<-as.double(rep(0,length(ids11)))
      ve02_12<-as.double(rep(0,length(ids12)))
    }
  }else{
      ve02<-as.double(x02)
      if (!is.null(times)) {
        Xn <- objectSurvival$Xnames02
        ve02_00   <- get_ve(newdata, id, Xn, ids00)
        ve02_01   <- get_ve(newdata, id, Xn, ids01)
        ve02_02_0 <- get_ve(newdata, id, Xn, ids02_0)
        ve02_02_1 <- get_ve(newdata, id, Xn, ids02_1)
        ve02_11 <- get_ve(newdata, id, Xn, ids11)
        ve02_12 <- get_ve(newdata, id, Xn, ids12)
      }
      
      }
  if(noVar[3]==1){
    ve12<-as.double(rep(0,N))
    if(!is.null(times)){
      ve12_00<-as.double(rep(0,length(ids00)))
      ve12_02_0<-as.double(rep(0,length(ids02_0)))
      ve12_02_1<-as.double(rep(0,length(ids02_1)))
      ve12_01<-as.double(rep(0,length(ids01)))
      ve12_11<-as.double(rep(0,length(ids11)))
      ve12_12<-as.double(rep(0,length(ids12)))
    }
  }else{
      ve12<-as.double(x12)
      
      if (!is.null(times)) {
        
        Xn <- objectSurvival$Xnames12
        ve12_00   <- get_ve(newdata, id, Xn, ids00)
        ve12_01   <- get_ve(newdata, id, Xn, ids01)
        ve12_02_0 <- get_ve(newdata, id, Xn, ids02_0)
        ve12_02_1 <- get_ve(newdata, id, Xn, ids02_1)
        ve12_11 <- get_ve(newdata, id, Xn, ids11)
        ve12_12 <- get_ve(newdata, id, Xn, ids12)
      }
      }

  
  t0<-rep(s,N)
  t1<-rep(horizon,N)
  
  ########################## time dep var ######################################
  outcome01<-objectSurvival$linktimedepXnames01
  outcome02<-objectSurvival$linktimedepXnames02
  outcome12<-objectSurvival$linktimedepXnames12
  outcome01<-outcome01[!is.na(outcome01)]
  outcome02<-outcome02[!is.na(outcome02)]
  outcome12<-outcome12[!is.na(outcome12)]
  
  p01<-length(outcome01)
  dimp01<-ifelse(p01>0,p01,1)
  p02<-length(outcome02)
  dimp02<-ifelse(p02>0,p02,1)
  p12<-length(outcome12)
  dimp12<-ifelse(p12>0,p12,1)
  
  if(p01==0 & p02==0){stop("No time dependent variable please refer to HIDeM::intensity for calculation of intensities with only time fixed covariates")}
  
 
  size1 <- size1+p01+p02+p12
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
  
  if(!("posfix"%in%names(objectSurvival))){
    objectSurvival$posfix<-NULL
  }
  if(length(objectSurvival$DYNidm)!=1){
    istop<-lapply(objectSurvival$DYNidm,FUN = function(x){
      if(x$istop%in%c(1,3)){return(T)}else{return(F)}
    })
    istop<-do.call(c,istop)
    if(sum(istop==F)!=length(istop)){warning("All the survival models did not converged")}
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
                         times=times,
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
                         NsampleFE=control$NsampleFE,
                         NidLoop=control$NidLoop)
                         # 31/08/2026 : do not estimate weights with our model 
                         # t0=t0,
                         # t1=t1,
                         # t2=t2,
                         # t3=t3,
                         # t4=t4,
                        # casei=ctime,
                        # isIntervalCensored=isIntervalCensored)
  
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
  
  #10/09/2026 : change integral from A(s,u) instead of A(0,u)
  #NtimePoints_0<-255
  NtimePoints_0<-240
  # 31/08/2026 : do not estimate weights with our model 
  # if(isIntervalCensored){
  # NtimePoints_1<-sum(ctime>=2)*(255*3)+sum(ctime==0)*30+sum(ctime==1)*255
  # NtimePoints_check<-sum(ctime>=2 & ctime<=5)*(255*3)+sum(ctime==0)*30+sum(ctime==1)*255+sum(ctime==6)*(255*2)
  # Ntime01<-NtimePoints_1*p01
  # Ntime02<-NtimePoints_1*p02
  # Ntime12<-NtimePoints_1*p12
  # NtimePoints_1max<-ifelse(sum(ctime==2)>0,255*3,
  #                          ifelse(sum(ctime==1)>0,255,30))
  # }


    
  for( m in unique(c(outcome01,outcome02,outcome12))){
    subdata<-dataY[dataY$Outcome==m,]
    
    # 31/08/2026 : do not estimate weights with our model 
    #if(dim(subdata)[1]!=(N*NtimePoints_0+NtimePoints_check)){stop("Prediction of marker ",m," could not be perform for each quadrature points")}
    if(is.null(times)){
    if(dim(subdata)[1]!=(N*NtimePoints_0)){stop("Prediction of marker ",m," could not be perform for each quadrature points")}
    }
    
  }
  
  #names(ctime)<-unique(newdata[,colnames(newdata)%in%id])
 
  # number of times each row's subject appears in a given ids vector (0 if absent)
  count_in <- function(y,ids) {
    n <- as.integer(table(ids)[as.character(y[[id]])])
    n[is.na(n)] <- 0L
    n
  }
  
  select_block <- function(n, off, size) {
    after & n > 0 & rank_after > off & rank_after <= off + n * size
  }

  if(length(outcome01)>=1){
    y01<-dataY[dataY$Outcome%in%outcome01,]
    #y01_0<-y01[order(y01[,colnames(y01)%in%id],y01$order),"Sample_1"]
    # order  by individual and timeline 
    # 31/08/2026 : do not estimate weights with our model 
     y01<-y01[order(y01[,colnames(y01)%in%id],y01$order),]
     y01_0<-y01[y01$order<=NtimePoints_0,"Sample_1"]
     
     if(!is.null(times)){
       
       n00   <- count_in(y01,ids00)
       n02_0 <- count_in(y01,ids02_0)
       n11   <- count_in(y01,ids11)
       n12   <- count_in(y01,ids12)
       n01   <- count_in(y01,ids01)
       n02_1 <- count_in(y01,ids02_1)
       
       # position among the rows after NtimePoints_0, within each subject
       # attention order restart at 1 for each subject and each Outcome Y 
       after      <- y01$order > NtimePoints_0
       grp        <- interaction(y01[[id]], y01$Outcome, drop = TRUE)
       rank_after <- ave(as.integer(after), grp, FUN = cumsum)
       
       # offsets, stacked in the order 00, 02_0, 01, 02_1
       off00   <- 0
       off02_0 <- off00   + 15  * n00
       off11   <- off02_0 + 16  * n02_0
       off12   <- off11 + 15  * n11
       off01   <- off12 + 16  * n12
       off02_1 <- off01   + 255 * n01
     
     
       if (length(ids00) > 0) {
         
         y01_00   <- y01[select_block(n00,   off00,   15),  "Sample_1"]
       }
       
       if (length(ids02_0) > 0) {
         
         y01_02_0   <- y01[select_block(n02_0, off02_0, 16),  "Sample_1"]
       }
       
       if (length(ids11) > 0) {
         
         y01_11   <- y01[select_block(n11, off11, 15),  "Sample_1"]
       }
       
       if (length(ids12) > 0) {
         
         y01_12   <- y01[select_block(n12, off12, 16),  "Sample_1"]
       }
       
       if (length(ids01) > 0) {
        
         y01_01   <- y01[select_block(n01,   off01,   255), "Sample_1"]
       }
       
       if (length(ids02_1) > 0) {
        
         y01_02_1 <- y01[select_block(n02_1, off02_1, 255), "Sample_1"]
       }
     }
     
    
  }else{
    y01_0<-rep(0,N*NtimePoints_0)
    if(!is.null(times)){
      y01_00<-rep(0,15*length(ids00))
      y01_02_0<-rep(0,16*length(ids02_0))
      y01_11<-rep(0,15*length(ids11))
      y01_12<-rep(0,16*length(ids12))
      y01_01<-rep(0,255*length(ids01))
      y01_02_1<-rep(0,255*length(ids02_1))
    }
   # y01_1<-rep(0,N*NtimePoints_1)
  }
 

  
  
  if(length(outcome02)>=1){
    y02<-dataY[dataY$Outcome%in%outcome02,]
    #y02_0<-y02[order(y02[,colnames(y02)%in%id],y02$order),"Sample_1"]
    # order  by individual and timeline 
    # 31/08/2026 : do not estimate weights with our model 
    y02<-y02[order(y02[,colnames(y02)%in%id],y02$order),]
    y02_0<-y02[y02$order<=NtimePoints_0,"Sample_1"]
    
    if(!is.null(times)){
      
      n00   <- count_in(y02,ids00)
      n02_0 <- count_in(y02,ids02_0)
      n11   <- count_in(y02,ids11)
      n12   <- count_in(y02,ids12)
      n01   <- count_in(y02,ids01)
      n02_1 <- count_in(y02,ids02_1)
      
      # position among the rows after NtimePoints_0, within each subject
      # attention order restart at 1 for each subject and each Outcome Y 
      after      <- y02$order > NtimePoints_0
      grp        <- interaction(y02[[id]], y02$Outcome, drop = TRUE)
      rank_after <- ave(as.integer(after), grp, FUN = cumsum)
      
      # offsets, stacked in the order 00, 02_0, 01, 02_1
      off00   <- 0
      off02_0 <- off00   + 15  * n00
      off11   <- off02_0 + 16  * n02_0
      off12   <- off11 + 15  * n11
      off01   <- off12 + 16  * n12
      off02_1 <- off01   + 255 * n01
      
      
      if (length(ids00) > 0) {
        
        y02_00   <- y02[select_block(n00,   off00,   15),  "Sample_1"]
      }
      
      if (length(ids02_0) > 0) {
        
        y02_02_0   <- y02[select_block(n02_0, off02_0, 16),  "Sample_1"]
      }
      
      if (length(ids11) > 0) {
        
        y02_11   <- y02[select_block(n11, off11, 15),  "Sample_1"]
      }
      
      if (length(ids12) > 0) {
        
        y02_12   <- y02[select_block(n12, off12, 16),  "Sample_1"]
      }
      
      if (length(ids01) > 0) {
        
        y02_01   <- y02[select_block(n01,   off01,   255), "Sample_1"]
      }
      
      if (length(ids02_1) > 0) {
        
        y02_02_1 <- y02[select_block(n02_1, off02_1, 255), "Sample_1"]
      }
    }
    
    
  }else{
    y02_0<-rep(0,N*NtimePoints_0)
    if(!is.null(times)){
      y02_00<-rep(0,15*length(ids00))
      y02_02_0<-rep(0,16*length(ids02_0))
      y02_11<-rep(0,15*length(ids11))
      y02_12<-rep(0,16*length(ids12))
      y02_01<-rep(0,255*length(ids01))
      y02_02_1<-rep(0,255*length(ids02_1))
    }
    # y01_1<-rep(0,N*NtimePoints_1)
  }
  
  
  if(length(outcome12)>=1){
    y12<-dataY[dataY$Outcome%in%outcome12,]
    #y12_0<-y12[order(y12[,colnames(y12)%in%id],y12$order),"Sample_1"]
    # 31/08/2026 : do not estimate weights with our model 
    # order  by individual and timeline 
     y12<-y12[order(y12[,colnames(y12)%in%id],y12$order),]
     y12_0<-y12[y12$order<=NtimePoints_0,"Sample_1"]
    # idcol <- colnames(y12)[colnames(y12) %in% id]
    # grpidx <- split(seq_len(nrow(y12)), y12[[idcol]])
    # 
    # y12_1 <- unlist(Map(function(idx, ct) y12[[4]][idx][if (ct == 6) TRUE else y12$order[idx] > NtimePoints_0],
    #                     grpidx, ctime[names(grpidx)]), use.names = FALSE)
    # 
    # 
     
     if(!is.null(times)){
       
       n00   <- count_in(y12,ids00)
       n02_0 <- count_in(y12,ids02_0)
       n11   <- count_in(y12,ids11)
       n12   <- count_in(y12,ids12)
       n01   <- count_in(y12,ids01)
       n02_1 <- count_in(y12,ids02_1)
       
       # position among the rows after NtimePoints_0, within each subject
       after      <- y12$order > NtimePoints_0
       grp        <- interaction(y12[[id]], y12$Outcome, drop = TRUE)
       rank_after <- ave(as.integer(after), grp, FUN = cumsum)
       
       # offsets, stacked in the order 00, 02_0, 01, 02_1
       off00   <- 0
       off02_0 <- off00   + 15  * n00
       off11   <- off02_0 + 16  * n02_0
       off12   <- off11 + 15  * n11
       off01   <- off12 + 16  * n12
       off02_1 <- off01   + 225 * n01
       
       
       if (length(ids00) > 0) {
         
         y12_00   <- y12[select_block(n00,   off00,   15),  "Sample_1"]
       }
       
       if (length(ids02_0) > 0) {
         
         y12_02_0   <- y12[select_block(n02_0, off02_0, 16),  "Sample_1"]
       }
       
       
       if (length(ids11) > 0) {
         
         y12_11   <- y12[select_block(n11,   off11,   15),  "Sample_1"]
       }
       
       if (length(ids12) > 0) {
         
         y12_12   <- y12[select_block(n12,   off12,   16),  "Sample_1"]
       }
       
       if (length(ids01) > 0) {
         
         y12_01   <- y12[select_block(n01,   off01,   225), "Sample_1"]
       }
       
       if (length(ids02_1) > 0) {
         
         y12_02_1 <- y12[select_block(n02_1, off02_1, 226), "Sample_1"]
       }
     }
  }else{
    y12_0<-rep(0,N*NtimePoints_0)
    #y12_1<-rep(0,N*NtimePoints_1)
    if(!is.null(times)){
      y12_00<-rep(0,15*length(ids00))
      y12_02_0<-rep(0,16*length(ids02_0))
      y12_11<-rep(0,15*length(ids11))
      y12_12<-rep(0,16*length(ids12))
      y12_01<-rep(0,225*length(ids01))
      y12_02_1<-rep(0,226*length(ids02_1))
    }
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
                             as.double(y01_0),
                             as.double(y02_0),
                             as.double(y12_0),
                             as.integer(p01),
                             as.integer(p02),
                             as.integer(p12),
                             as.integer(dimp01),
                             as.integer(dimp02),
                             as.integer(dimp12),
                             as.integer(NtimePoints_0),
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
    if(!is.null(times)){
      if (length(ids00) > 0){
        res<-rep(0,length(ids00))
        out00<- tryCatch({  .Fortran("p00timedep",
                                    ## input
                                    as.double(binit),
                                    as.integer(size_V),
                                    as.double(knots01),
                                    as.double(knots02),
                                    as.double(knots12),
                                    as.integer(length(ids00)),
                                    as.integer(nknots01),
                                    as.integer(nknots02),
                                    as.integer(nknots12),
                                    as.double(ve01_00),
                                    as.double(ve02_00),
                                    as.double(ve12_00),
                                    as.double(y01_00),
                                    as.double(y02_00),
                                    as.double(y12_00),
                                    as.integer(p01),
                                    as.integer(p02),
                                    as.integer(p12),
                                    as.integer(dimp01),
                                    as.integer(dimp02),
                                    as.integer(dimp12),
                                    as.integer(15),
                                    as.integer(dimnva01),
                                    as.integer(dimnva02),
                                    as.integer(dimnva12),
                                    as.integer(nvat01),
                                    as.integer(nvat02),
                                    as.integer(nvat12),
                                    as.double(times$start[ids00]),
                                    as.double(times$end[ids00]),
                                    likelihood_res=as.double(res),
                                    PACKAGE="HIDeM")$likelihood_res
        }, error = function(e) {
          # Return NULL on error to skip this patient
          NULL
        })
      }else{out00<-NULL}
      
      if (length(ids02_0) > 0){
        res<-rep(0,length(ids02_0))
        out02_0<- tryCatch({  .Fortran("p020timedep",
                                     ## input
                                     as.double(binit),
                                     as.integer(size_V),
                                     as.double(knots01),
                                     as.double(knots02),
                                     as.double(knots12),
                                     as.integer(length(ids02_0)),
                                     as.integer(nknots01),
                                     as.integer(nknots02),
                                     as.integer(nknots12),
                                     as.double(ve01_02_0),
                                     as.double(ve02_02_0),
                                     as.double(ve12_02_0),
                                     as.double(y01_02_0),
                                     as.double(y02_02_0),
                                     as.double(y12_02_0),
                                     as.integer(p01),
                                     as.integer(p02),
                                     as.integer(p12),
                                     as.integer(dimp01),
                                     as.integer(dimp02),
                                     as.integer(dimp12),
                                     as.integer(16),
                                     as.integer(dimnva01),
                                     as.integer(dimnva02),
                                     as.integer(dimnva12),
                                     as.integer(nvat01),
                                     as.integer(nvat02),
                                     as.integer(nvat12),
                                     as.double(times$start[ids02_0]),
                                     as.double(times$end[ids02_0]),
                                     likelihood_res=as.double(res),
                                     PACKAGE="HIDeM")$likelihood_res
        }, error = function(e) {
          # Return NULL on error to skip this patient
          NULL
        })
      }else{out02_0<-NULL}
      
      
      if (length(ids11) > 0){
        res<-rep(0,length(ids11))
        out11<- tryCatch({  .Fortran("p11timedep",
                                     ## input
                                     as.double(binit),
                                     as.integer(size_V),
                                     as.double(knots01),
                                     as.double(knots02),
                                     as.double(knots12),
                                     as.integer(length(ids11)),
                                     as.integer(nknots01),
                                     as.integer(nknots02),
                                     as.integer(nknots12),
                                     as.double(ve01_11),
                                     as.double(ve02_11),
                                     as.double(ve12_11),
                                     as.double(y01_11),
                                     as.double(y02_11),
                                     as.double(y12_11),
                                     as.integer(p01),
                                     as.integer(p02),
                                     as.integer(p12),
                                     as.integer(dimp01),
                                     as.integer(dimp02),
                                     as.integer(dimp12),
                                     as.integer(15),
                                     as.integer(dimnva01),
                                     as.integer(dimnva02),
                                     as.integer(dimnva12),
                                     as.integer(nvat01),
                                     as.integer(nvat02),
                                     as.integer(nvat12),
                                     as.double(times$start[ids11]),
                                     as.double(times$end[ids11]),
                                     likelihood_res=as.double(res),
                                     PACKAGE="HIDeM")$likelihood_res
        }, error = function(e) {
          # Return NULL on error to skip this patient
          NULL
        })
      }else{out11<-NULL}
      
      
      if (length(ids12) > 0){
        res<-rep(0,length(ids12))
        out12<- tryCatch({  .Fortran("p12timedep",
                                     ## input
                                     as.double(binit),
                                     as.integer(size_V),
                                     as.double(knots01),
                                     as.double(knots02),
                                     as.double(knots12),
                                     as.integer(length(ids12)),
                                     as.integer(nknots01),
                                     as.integer(nknots02),
                                     as.integer(nknots12),
                                     as.double(ve01_12),
                                     as.double(ve02_12),
                                     as.double(ve12_12),
                                     as.double(y01_12),
                                     as.double(y02_12),
                                     as.double(y12_12),
                                     as.integer(p01),
                                     as.integer(p02),
                                     as.integer(p12),
                                     as.integer(dimp01),
                                     as.integer(dimp02),
                                     as.integer(dimp12),
                                     as.integer(16),
                                     as.integer(dimnva01),
                                     as.integer(dimnva02),
                                     as.integer(dimnva12),
                                     as.integer(nvat01),
                                     as.integer(nvat02),
                                     as.integer(nvat12),
                                     as.double(times$start[ids12]),
                                     as.double(times$end[ids12]),
                                     likelihood_res=as.double(res),
                                     PACKAGE="HIDeM")$likelihood_res
        }, error = function(e) {
          # Return NULL on error to skip this patient
          NULL
        })
      }else{out12<-NULL}
      
      
      
      if (length(ids02_1) > 0){
        res<-rep(0,length(ids02_1))
        out02_1<- tryCatch({  .Fortran("p021timedep",
                                     ## input
                                     as.double(binit),
                                     as.integer(size_V),
                                     as.double(knots01),
                                     as.double(knots02),
                                     as.double(knots12),
                                     as.integer(length(ids02_1)),
                                     as.integer(nknots01),
                                     as.integer(nknots02),
                                     as.integer(nknots12),
                                     as.double(ve01_02_1),
                                     as.double(ve02_02_1),
                                     as.double(ve12_02_1),
                                     as.double(y01_02_1),
                                     as.double(y02_02_1),
                                     as.double(y12_02_1),
                                     as.integer(p01),
                                     as.integer(p02),
                                     as.integer(p12),
                                     as.integer(dimp01),
                                     as.integer(dimp02),
                                     as.integer(dimp12),
                                     as.integer(255),
                                     as.integer(226),
                                     as.integer(dimnva01),
                                     as.integer(dimnva02),
                                     as.integer(dimnva12),
                                     as.integer(nvat01),
                                     as.integer(nvat02),
                                     as.integer(nvat12),
                                     as.double(times$start[ids02_1]),
                                     as.double(times$end[ids02_1]),
                                     likelihood_res=as.double(res),
                                     PACKAGE="HIDeM")$likelihood_res
        }, error = function(e) {
          # Return NULL on error to skip this patient
          NULL
        })
      }else{out02_1<-NULL}
      
      
      if (length(ids01) > 0){
        res<-rep(0,length(ids01))
        out01<- tryCatch({  .Fortran("p01timedep",
                                     ## input
                                     as.double(binit),
                                     as.integer(size_V),
                                     as.double(knots01),
                                     as.double(knots02),
                                     as.double(knots12),
                                     as.integer(length(ids01)),
                                     as.integer(nknots01),
                                     as.integer(nknots02),
                                     as.integer(nknots12),
                                     as.double(ve01_01),
                                     as.double(ve02_01),
                                     as.double(ve12_01),
                                     as.double(y01_01),
                                     as.double(y02_01),
                                     as.double(y12_01),
                                     as.integer(p01),
                                     as.integer(p02),
                                     as.integer(p12),
                                     as.integer(dimp01),
                                     as.integer(dimp02),
                                     as.integer(dimp12),
                                     as.integer(255),
                                     as.integer(225),
                                     as.integer(dimnva01),
                                     as.integer(dimnva02),
                                     as.integer(dimnva12),
                                     as.integer(nvat01),
                                     as.integer(nvat02),
                                     as.integer(nvat12),
                                     as.double(times$start[ids01]),
                                     as.double(times$end[ids01]),
                                     likelihood_res=as.double(res),
                                     PACKAGE="HIDeM")$likelihood_res
        }, error = function(e) {
          # Return NULL on error to skip this patient
          NULL
        })
      }else{out01<-NULL}
    }
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
                                as.double(y01_0),
                                as.double(y02_0),
                                as.double(y12_0),
                                as.integer(p01),
                                as.integer(p02),
                                as.integer(p12),
                                as.integer(dimp01),
                                as.integer(dimp02),
                                as.integer(dimp12),
                                as.integer(NtimePoints_0),
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
    
    if(!is.null(times)){
      if (length(ids00) > 0) {
        res<-rep(0,length(ids00))
        
        out00<- tryCatch({   .Fortran("p00weibtimedep",
                                     ## input
                                     as.double(binit),
                                     as.integer(size_V),
                                     as.integer(length(ids00)),
                                     as.double(ve01_00),
                                     as.double(ve02_00),
                                     as.double(ve12_00),
                                     as.double(y01_00),
                                     as.double(y02_00),
                                     as.double(y12_00),
                                     as.integer(p01),
                                     as.integer(p02),
                                     as.integer(p12),
                                     as.integer(dimp01),
                                     as.integer(dimp02),
                                     as.integer(dimp12),
                                     as.integer(15),
                                     as.integer(dimnva01),
                                     as.integer(dimnva02),
                                     as.integer(dimnva12),
                                     as.integer(nvat01),
                                     as.integer(nvat02),
                                     as.integer(nvat12),
                                     as.double(times$start[ids00]),
                                     as.double(times$end[ids00]),
                                     likelihood_res=as.double(res),
                                     PACKAGE="HIDeM")$likelihood_res
        }, error = function(e) {
          # Return NULL on error to skip this patient
          NULL
        })
      }else{out00<-NULL}
      
      if (length(ids02_0) > 0) {
        res<-rep(0,length(ids02_0))
        out02_0<- tryCatch({   .Fortran("p020weibtimedep",
                                      ## input
                                      as.double(binit),
                                      as.integer(size_V),
                                      as.integer(length(ids02_0)),
                                      as.double(ve01_02_0),
                                      as.double(ve02_02_0),
                                      as.double(ve12_02_0),
                                      as.double(y01_02_0),
                                      as.double(y02_02_0),
                                      as.double(y12_02_0),
                                      as.integer(p01),
                                      as.integer(p02),
                                      as.integer(p12),
                                      as.integer(dimp01),
                                      as.integer(dimp02),
                                      as.integer(dimp12),
                                      as.integer(16),
                                      as.integer(dimnva01),
                                      as.integer(dimnva02),
                                      as.integer(dimnva12),
                                      as.integer(nvat01),
                                      as.integer(nvat02),
                                      as.integer(nvat12),
                                      as.double(times$start[ids02_0]),
                                      as.double(times$end[ids02_0]),
                                      likelihood_res=as.double(res),
                                      PACKAGE="HIDeM")$likelihood_res
        }, error = function(e) {
          # Return NULL on error to skip this patient
          NULL
        })
      }else{out02_0<-NULL}
      
      
      if (length(ids11) > 0) {
        res<-rep(0,length(ids11))
        
        out11<- tryCatch({   .Fortran("p11weibtimedep",
                                      ## input
                                      as.double(binit),
                                      as.integer(size_V),
                                      as.integer(length(ids11)),
                                      as.double(ve01_11),
                                      as.double(ve02_11),
                                      as.double(ve12_11),
                                      as.double(y01_11),
                                      as.double(y02_11),
                                      as.double(y12_11),
                                      as.integer(p01),
                                      as.integer(p02),
                                      as.integer(p12),
                                      as.integer(dimp01),
                                      as.integer(dimp02),
                                      as.integer(dimp12),
                                      as.integer(15),
                                      as.integer(dimnva01),
                                      as.integer(dimnva02),
                                      as.integer(dimnva12),
                                      as.integer(nvat01),
                                      as.integer(nvat02),
                                      as.integer(nvat12),
                                      as.double(times$start[ids11]),
                                      as.double(times$end[ids11]),
                                      likelihood_res=as.double(res),
                                      PACKAGE="HIDeM")$likelihood_res
        }, error = function(e) {
          # Return NULL on error to skip this patient
          NULL
        })
      }else{out11<-NULL}
      
      
      
      
      if (length(ids12) > 0) {
        res<-rep(0,length(ids12))
        
        out12<- tryCatch({   .Fortran("p12weibtimedep",
                                      ## input
                                      as.double(binit),
                                      as.integer(size_V),
                                      as.integer(length(ids12)),
                                      as.double(ve01_12),
                                      as.double(ve02_12),
                                      as.double(ve12_12),
                                      as.double(y01_12),
                                      as.double(y02_12),
                                      as.double(y12_12),
                                      as.integer(p01),
                                      as.integer(p02),
                                      as.integer(p12),
                                      as.integer(dimp01),
                                      as.integer(dimp02),
                                      as.integer(dimp12),
                                      as.integer(16),
                                      as.integer(dimnva01),
                                      as.integer(dimnva02),
                                      as.integer(dimnva12),
                                      as.integer(nvat01),
                                      as.integer(nvat02),
                                      as.integer(nvat12),
                                      as.double(times$start[ids12]),
                                      as.double(times$end[ids12]),
                                      likelihood_res=as.double(res),
                                      PACKAGE="HIDeM")$likelihood_res
        }, error = function(e) {
          # Return NULL on error to skip this patient
          NULL
        })
      }else{out12<-NULL}
      
      if (length(ids01) > 0) {
        res<-rep(0,length(ids01))
        out01<- tryCatch({   .Fortran("p01weibtimedep",
                                        ## input
                                        as.double(binit),
                                        as.integer(size_V),
                                        as.integer(length(ids01)),
                                        as.double(ve01_01),
                                        as.double(ve02_01),
                                        as.double(ve12_01),
                                        as.double(y01_01),
                                        as.double(y02_01),
                                        as.double(y12_01),
                                        as.integer(p01),
                                        as.integer(p02),
                                        as.integer(p12),
                                        as.integer(dimp01),
                                        as.integer(dimp02),
                                        as.integer(dimp12),
                                        as.integer(255),
                                        as.integer(225),
                                        as.integer(dimnva01),
                                        as.integer(dimnva02),
                                        as.integer(dimnva12),
                                        as.integer(nvat01),
                                        as.integer(nvat02),
                                        as.integer(nvat12),
                                        as.double(times$start[ids01]),
                                        as.double(times$end[ids01]),
                                        likelihood_res=as.double(res),
                                        PACKAGE="HIDeM")$likelihood_res
        }, error = function(e) {
          # Return NULL on error to skip this patient
          NULL
        })
      }else{out01<-NULL}
      
      if (length(ids02_1) > 0) {
        res<-rep(0,length(ids02_1))
        out02_1<- tryCatch({   .Fortran("p021weibtimedep",
                                        ## input
                                        as.double(binit),
                                        as.integer(size_V),
                                        as.integer(length(ids02_1)),
                                        as.double(ve01_02_1),
                                        as.double(ve02_02_1),
                                        as.double(ve12_02_1),
                                        as.double(y01_02_1),
                                        as.double(y02_02_1),
                                        as.double(y12_02_1),
                                        as.integer(p01),
                                        as.integer(p02),
                                        as.integer(p12),
                                        as.integer(dimp01),
                                        as.integer(dimp02),
                                        as.integer(dimp12),
                                        as.integer(255),
                                        as.integer(226),
                                        as.integer(dimnva01),
                                        as.integer(dimnva02),
                                        as.integer(dimnva12),
                                        as.integer(nvat01),
                                        as.integer(nvat02),
                                        as.integer(nvat12),
                                        as.double(times$start[ids02_1]),
                                        as.double(times$end[ids02_1]),
                                        likelihood_res=as.double(res),
                                        PACKAGE="HIDeM")$likelihood_res
        }, error = function(e) {
          # Return NULL on error to skip this patient
          NULL
        })
      }else{out02_1<-NULL}
    }
    
    # 31/08/2026 : do not estimate weights with our model 
    # if(isIntervalCensored){
    # res<-rep(0,4*N)
    # out2<- tryCatch({   .Fortran("idmciweibtimedep",
    #                              ## input
    #                              as.double(binit),
    #                              as.integer(size_V),
    #                              as.integer(N),
    #                              as.double(ve01),
    #                              as.double(ve02),
    #                              as.double(ve12),
    #                              as.double(y01_1),
    #                              as.double(y02_1),
    #                              as.double(y12_1),
    #                              as.integer(p01),
    #                              as.integer(p02),
    #                              as.integer(p12),
    #                              as.integer(dimp01),
    #                              as.integer(dimp02),
    #                              as.integer(dimp12),
    #                              as.integer(NtimePoints_1max),
    #                              as.integer(dimnva01),
    #                              as.integer(dimnva02),
    #                              as.integer(dimnva12),
    #                              as.integer(nvat01),
    #                              as.integer(nvat02),
    #                              as.integer(nvat12),
    #                              as.double(t1),
    #                              as.double(t2),
    #                              as.integer(ctime),
    #                              as.integer(Ntime01),
    #                              as.integer(Ntime02),
    #                              as.integer(Ntime12),
    #                              likelihood_res=as.double(res),
    #                              PACKAGE="HIDeM")$likelihood_res
    # }, error = function(e) {
    #   # Return NULL on error to skip this patient
    #   NULL
    # })
    # }else{out2<-NULL}
  }
  
  if(control$return.data==T){
    predicted.newdata<-dataY
  }else{predicted.newdata<-NULL}
  
  if(!is.null(out1)){
    CIF01<-out1[1:N]
    id_above1<-which(out1[1:N]>1)
    if(length(id_above1)>0){warnings("Attention some individual have predicted probabilities above 1")}
  }else{
    CIF01<-NULL
    id_above1<-NULL
  }
  
 
  if(!is.null(times)){
    
    times$prediction_p00<-rep(NA,length(times$id))
    if(sum(times$p00==T)>0){
      times$prediction_p00[which(times$p00==T)]<-out00
    }
    
    times$prediction_p01<-rep(NA,length(times$id))
    if(sum(times$p01==T)>0){
    times$prediction_p01[which(times$p01==T)]<-out01
    }
    
    times$prediction_p02_0<-rep(NA,length(times$id))
    if(sum(times$p02_0==T)>0){
      times$prediction_p02_0[which(times$p02_0==T)]<-out02_0
    }
    
    times$prediction_p02_1<-rep(NA,length(times$id))
    if(sum(times$p02_1==T)>0){
      times$prediction_p02_1[which(times$p02_1==T)]<-out02_1
    }
    
    times$prediction_p11<-rep(NA,length(times$id))
    if(sum(times$p11==T)>0){
      times$prediction_p11[which(times$p11==T)]<-out11
    }
    
    times$prediction_p12<-rep(NA,length(times$id))
    if(sum(times$p12==T)>0){
      times$prediction_p12[which(times$p12==T)]<-out12
    }
  }
  
  # 31/08/2026 : do not estimate weights with our model 
  # if(!is.null(out2)){
  #   BS<-sum(out2[1:N])/N
  #   p01<-out2[(N+1):(2*N)]
  #   p00<-out2[(2*N+1):(3*N)]
  #   weights<-out2[(3*N+1):(4*N)]
  # }else{
  #   BS<-p01<-p00<-weights<-NULL
  # }
  # cases<-ctime
  # Ltime<-t1
  # Rtime<-t2
  # res<-list(CIF01=CIF01,
  #           BS=BS,
  #           p01=p01,
  #           p00=p00,
  #           weights=weights,
  #           cases=cases,
  #           s=s,
  #           horizon=horizon,
  #           Ltime=Ltime,
  #           Rtime=Rtime,
  #           predicted.newdata= predicted.newdata)
  
  res<-list(CIF01=CIF01,
            times=times,
            s=s,
            id_above1=id_above1,
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
    b_all<-rep(NA,length(x$b)+length(posfix))
    if(x$istop%in%c(1,3)){
      if(!is.null(posfix)){
        b_all[posfix]<-bfix
        b_all[-posfix]<-x$b
      }else{ b_all<-x$b}
    }
      return(b_all)})
  b<-do.call(rbind,b)
  b<-na.omit(b)
  b<-colSums(b)/Nrep
  return(b)
}


get_ve <- function(newdata, id, Xn, ids, by_row = FALSE) {
  
  # nothing to extract
  if (length(Xn) == 0 || length(ids) == 0) return(numeric(0))
  
  # all needed columns must exist
  miss <- setdiff(c(id, Xn), colnames(newdata))
  if (length(miss) > 0)
    stop("Columns not found in newdata: ", paste(miss, collapse = ", "))
  
  # one row per subject (first row), used as a lookup table
  cols <- unique(c(id, Xn))
  x <- newdata[!duplicated(newdata[[id]]), cols, drop = FALSE]
  
  # one position per entry of ids: repeated ids give repeated rows
  pos <- match(ids, x[[id]])
  if (anyNA(pos))
    stop("Ids not found in newdata: ",
         paste(unique(ids[is.na(pos)]), collapse = ", "))
  
  m <- as.matrix(x[pos, Xn, drop = FALSE])
  
  # make sure everything is numeric (logical is fine, factor/character is not)
  if (!is.numeric(m)) {
    m2 <- suppressWarnings(matrix(as.numeric(m), nrow = nrow(m)))
    if (anyNA(m2) && !anyNA(m))
      stop("Non-numeric covariate(s) in: ", paste(Xn, collapse = ", "),
           ". Build a numeric design matrix first (model.matrix).")
    m <- m2
  }
  
  # by_row = FALSE: column-major (all subjects for covariate 1, then covariate 2, ...)
  # by_row = TRUE : one subject after another
  as.double(if (by_row) t(m) else m)
}
