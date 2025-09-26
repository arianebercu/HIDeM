### Code:
#' @title Fit an dynamic illness-death model with limited and numerous explanatory variables, including time dependent covariates
#' @description
#' Fit an illness-death model with a limited number of explanatory variables
#'  using semi-parametric approach with a linear combination of M-splines or 
#'  Weibull distribution on the transition intensities. 
#' Fit an illness-death model with a large number of explanatory variables 
#' using regularization approach of the log-likelihood, i.e. lasso, elastic-net,
#'  MCP or SCAD penalty.
#' Left-truncated, right-censored, and interval-censored data are allowed.
#' State 0 corresponds to the initial state, state 1 to the transient one,
#' state 2 to the absorbant one. The allowed transitions are: 0 --> 1, 0 --> 2
#' and 1 --> 2.
#'
#' The estimated parameters are obtained using the robust Marquardt algorithm
#' (Marquardt, 1963) or a proximal gradient algorithm. 
#'
#' @param formula01 A formula specifying a regression model for the
#' \code{0 --> 1} transition from the initial state to the transient
#' state of the illness-death model.  The right hand side of the
#' formula specifies the covariate terms, and the left hand side must
#' be an event history object as returned by the function \code{Hist}.
#' @param formula02 A formula specifying a regression model for the
#' \code{0 --> 2} transition from the initial state to the absorbing
#' state. The left hand side must be equal to the left hand side of
#' \code{formula01}. If missing it is set to \code{formula01}.
#' @param formula12 A formula specifying a regression model for the
#' \code{1 --> 2} transition from the transient state to the absorbing
#' state.  operator is not required. If missing it is set to
#' \code{formula01}.
#' @param data A survival data frame in which to interpret the variables of
#' \code{formula01}, \code{formula02} and \code{formula12}.
#' @param dataLongi A longitudinal data frame for the time dependent covariates
#' This argument is only used for models with M-splines.
#' @param na.action how NAs are treated. The default is first, any
#' na.action attribute of data, second a na.action setting of options,
#' and third 'na.fail' if that is unset. The 'factory-fresh' default
#' is na.omit. Another possible value is NULL.
#' @param clustertype in which cluster to work
#' @param nproc number of cluster
#' @param envir working environment 
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
#' @references D. Marquardt (1963). An algorithm for least-squares estimation
#' of nonlinear parameters.  \emph{SIAM Journal of Applied Mathematics},
#' 431-441.
#' @keywords illness-death
##'@importFrom grDevices col2rgb
##'@importFrom graphics lines
##'@importFrom graphics par
##'@importFrom graphics polygon segments
##'@importFrom stats as.formula formula integrate model.frame model.matrix na.fail na.omit pchisq pweibull qnorm quantile terms update.formula
#' @useDynLib HIDeM
#' @export
DYNpredY <- function(formula01,
                   formula02,
                   formula12,
                   data,
                   formLong,
                   Longitransition,
                   threshold,
                   dataLongi,
                   timeVar,
                   id,
                   Nsample,
                   methodJM=list(functional_forms,
                                 n_iter=3500, 
                                 n_burnin=500,
                                 n_thin=1,
                                 n_chains=3,
                                 seed=1234),
                   methodINLA=list(family="gaussian",
                                   basRisk= "rw1",
                                   assoc=NULL),
                   nproc=1,
                   clustertype="FORK",
                   Ypredmethod="gauss",
                   NtimePoints=250,
                   envir=parent.frame()){
  
  
  
  # {{{ check formula
  #################################################################################
  ################################ check what has been given ##################
  #################################################################################
  
  call <- match.call()
  ptm <- proc.time()
  

    if(missing(methodJM) & missing(methodINLA)){
      stop("Need to specify the method used for time-depend covariates by giving methodINLA or methodJM")
    }
    if(!missing(methodINLA)){
      methodlongi<-"INLA"
      print("Time-dependent covariates estimation will be performed using INLA")
    }else{
      methodlongi<-"JM"
      print("Time-dependent covariates estimation will be performed using JMBayes2")
    }
  if(missing(formula01))stop("Argument formula01 is missing.")
  if(missing(formula02))stop("Argument formula02 is missing.")
  if(!inherits(formula01,"formula"))stop("The argument formula01 must be a class 'formula'.")
  if(!inherits(formula02,"formula"))stop("The argument formula02 must be a class 'formula'.")
  
  if(!Ypredmethod%in%c("gauss","equi"))stop("The argument Ypredmethod needs to be either gauss or equi")
  
  ## if(missing(formula02)) formula02 <- formula01
  if(missing(formula12)) formula12 <- formula02
  # }}}
  # {{{ evaluate formula in data
  if(missing(data)) stop("Need a data frame.")
  if(sum(is.na(data))>0)stop("Need a survival data frame with no missing data.")
  
  if(!inherits(data,"data.frame"))stop("Argument 'data' must be a data.frame")
  
 
    if(missing(dataLongi)) stop("Need a dataLongi frame.")
    if(sum(is.na(dataLongi))>0)stop("Need a longitudinal data frame with no missing data.")
    
    if(!inherits(dataLongi,"data.frame"))stop("Argument 'dataLongi' must be a data.frame")
  
  ############################################################################
  ############################### get database defined by formulas ###########
  ############################################################################
  m <- match.call()
  m01 <- m02 <- m12 <- m[match(c("","data","na.action"),names(m),nomatch=0)]
  m01$formula <- formula01
  m02$formula <- formula02
  m12$formula <- formula12
  m01[[1]] <- m02[[1]] <- m12[[1]] <- as.name("model.frame")
  
  #################################################################################
  #################### dealing with missing data ##################################
  #################################################################################
 
  if(anyNA(data)){
    variables=unique(c(all.vars(formula01),all.vars(formula02),all.vars(formula12)))
    data=data[,variables]
    data=na.omit(data)
    m01[[2]] <- m02[[2]] <- m12[[2]] <- data
  }
  
  
    if(anyNA(dataLongi)){
      variables=unique(c(all.vars(formula01),all.vars(formula02),all.vars(formula12)))
      dataLongi=dataLongi[,variables]
      dataLongi=na.omit(dataLongi)
      # m01[[2]] <- m02[[2]] <- m12[[2]] <- data
    }
  
  m01 <- eval(m01,parent.frame())
  m02 <- eval(m02,parent.frame())
  m12 <- eval(m12,parent.frame())
  
  
  responseTrans <- stats::model.response(m01)
  responseAbs <- stats::model.response(m02)
  
  
  if(is.null(Longitransition)){
    Longitransition<-list()
    length(Longitransition)<-length(formLong)
    Longitransition<-lapply(Longitransition, function(x) c("01", "02", "12"))
  }else{
    
    if (!inherits(Longitransition, "list")){
      cli_abort(c(
        "{.var Longitransition} must be a list object",
        "x" = "You've supplied a {.cls {class(Longitransition)}} object"
      ))
    }
    if(length(Longitransition)!=length(formLong)){
      stop(paste0("Longitransition must be a list object of length ",length(formLong)))
    }
    
    for(k in 1:length(Longitransition)){
      if(!inherits(Longitransition[[k]], "character")){
        stop("Each element of the list Longitransition must be a character")
      }
      if(any(!Longitransition[[k]]%in%c("01","02","12"))){
        stop("Each element of the list Longitransition must be a character taking values in : 01,02 or/and 12")
      }
    }
    
  }
  
  

  
  #################################################################################
  ####################  prepare censored event times  #############################
  #################################################################################
  
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
  
  
  if (attr(responseAbs,"cens.type")=="intervalCensored") stop("No method available when the transtion to the absorbing state is interval censored.")
  if (isIntervalCensored && any(Rtime<Ltime)) stop("Misspecified transitition times:\nSome left interval limits are greater than the corresponding right limits.")
  
  #################################################################################
  ####################### check entry parameters ##################################
  #################################################################################
  
  
  if(!inherits(NtimePoints,c("numeric","integer"))|(NtimePoints!=floor(NtimePoints)))stop("NtimePoints has to be an integer greater than 100.")
  if(NtimePoints<100)stop("NtimePoints has to be an integer greater than 100.")
  
  if(!inherits(nproc,c("numeric","integer"))|(nproc!=floor(nproc)))stop("nproc has to be an integer.")

  troncature<-ifelse(truncated==T,1,0)
  
  #################################################################################
  ####################### define profile of subjects ##############################
  #################################################################################
  
  if(truncated==1){
    t0<-entrytime
  }else{t0<-rep(0,N)}
  
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
  
  ############################################################################
  #################### defines knots placements for splines ##################
  #####################         and initiate values         ##################
  ############################################################################
  
  if(methodlongi=="INLA"){
    
    
    if(!inherits(formLong,"list")){stop("formLong should be a list of equation")}
    type<-lapply(formLong,FUN=function(x)!inherits(x,"formula"))
    type<-unlist(type)
    if(any(type)){stop("formLong should be a list of equation")}
    type<-lapply(formLong,FUN=function(x){
      match(as.character(x[[2]]),colnames(dataLongi),nomatch=0)
    })
    type<-unlist(type)
    if(any(type==0)){stop("Each outcome in formLong should be in dataLongi")}
    
    # keep names of Y 
    ynames<-unlist(lapply(formLong,FUN=function(x){as.character(x[[2]])}))
    
  }else{
    
    if(!inherits(formLong,"list")){stop("formLong should be a list of lme models")}
    type<-lapply(formLong,FUN=function(x)!inherits(x,c("MixMod","lme")))
    type<-unlist(type)
    if(any(type)){stop("formLong should be a list of lme or MixMod models")}
    type<-lapply(formLong,FUN=function(x){
      if(class(x)=="lme"){
        match(as.character(x$terms[[2]]),colnames(dataLongi),nomatch=0)
      }else{
        match(as.character(x$Terms$termsX[[2]]),colnames(dataLongi),nomatch=0)
      }
    })
    type<-unlist(type)
    if(any(type==0)){stop("Each outcome in formLong should be in dataLongi")}
    
    
    CR_forms<-list()
    length(CR_forms)<-length(formLong)
    ynames<-lapply(formLong,FUN=function(x){
      if(class(x)=="lme"){
        as.character(x$terms[[2]])
      }else{
        as.character(x$Terms$termsX[[2]])
      }
    })
    ynames<-unlist(ynames)
    
    for(m in 1:length(formLong)){
      if("value"%in%methodJM$functional_forms[[m]] & "slope"%in%methodJM$functional_forms[[m]]){
        CR_forms[[m]]<-as.formula(paste0("~value(",ynames[m],"):CR + slope(",ynames[m],"):CR"))
      }else{
        if("value"%in%methodJM$functional_forms[[m]]){
          CR_forms[[m]]<-as.formula(paste0("~value(",ynames[m],"):CR"))
        }else{
          CR_forms[[m]]<-as.formula(paste0("~slope(",ynames[m],"):CR"))
        }
      }
    }
    names(CR_forms)<-ynames
    
    
  }
  
  
  outcome<-ynames
  outcome01<-  ynames[unlist(lapply(Longitransition,FUN=function(x){
    if("01"%in%x){
      return(T)}else{return(F)}
  }))]
  outcome02<-  ynames[unlist(lapply(Longitransition,FUN=function(x){
    if("02"%in%x){
      return(T)}else{return(F)}
  }))]
  
  outcome12<-  ynames[unlist(lapply(Longitransition,FUN=function(x){
    if("12"%in%x){
      return(T)}else{return(F)}
  }))]
  
  if(length(outcome01)>=1){
    p01<-length(outcome01)
    dimp01<-length(outcome01)
  }else{
    p01<-0
    dimp01<-1
  }
  
  if(length(outcome02)>=1){
    p02<-length(outcome02)
    dimp02<-length(outcome02)
  }else{
    p02<-0
    dimp02<-1
  }
  
  if(length(outcome12)>=1){
    p12<-length(outcome12)
    dimp12<-length(outcome12)
    
  }else{
    p12<-0
    dimp12<-1
  }
  

  ################################################################################
  ############################## RUN INLA or JM MODEL ################################## ################################################################################
  if(missing(threshold)){
    threshold<-max(t2-t1)
  }else{
    if(!inherits(threshold,c("numeric","integer")))stop("threshold has to be a numeric or integer.")
    if(threshold<=0)stop("threshold has be strictly positive")
    
  }
  
  
  
  if(missing(id)){
    stop("Need to specify id variable")
  }else{
    if(!inherits(id,"character")){stop("id need specify the column name")}
    if(class(dataLongi[,colnames(dataLongi)%in%id])!="integer"){stop("ID for subject needs to be an integer ")}
    if(length(id)!=1){stop("id need to be a character")}
    if(!id%in%colnames(data)|!id%in%colnames(dataLongi)){stop("id need to be in data and dataLongi")}}
  
  
  # competing risk definition 
  iddCR<-ifelse(idm==1,0,
                ifelse(idd==1 & t3<=t2+threshold,1,0))
  TimeCR<-ifelse(idm==1,(t1+t2)/2,
                 ifelse(idd==1 & t3<=t2+threshold,t3,t2))
  
    if(methodlongi=="INLA"){
      
      
      
      if(truncated==1){
        formSurv<-list(inla.surv(time=TimeCR, event=idm,truncation=t0) ~ -1,
                       inla.surv(time=TimeCR, event=iddCR,truncation=t0) ~ -1)
      }else{
        formSurv<-list(inla.surv(time=TimeCR, event=idm) ~ 1,
                       inla.surv(time=TimeCR, event=iddCR) ~ 1)
      }
      
      # set right environement for formulas 
      formSurv <- lapply(formSurv, function(f) {
        environment(f) <-envir
        f
      })
      
      formLong <- lapply(formLong, function(f) {
        environment(f) <- envir
        f
      })
      
      dataINLA<-data.frame(ID=data[,colnames(data)%in%id],
                           idm=idm,
                           iddCR=iddCR,
                           TimeCR=TimeCR)
      
      if(missing(Nsample)){
        Nsample<-100}else{
          if(!inherits(Nsample,c("integer","numeric")))stop("Nsample need to be a numeric or integer")
          if(length(Nsample)!=1)stop("Nsample need to be a numeric or integer")
        }
      
      dataY<-INLAidm(timeVar = timeVar,
                     family = methodINLA$family,
                     basRisk = methodINLA$basRisk,
                     assoc = methodINLA$assoc,
                     truncated=truncated,
                     formLong=formLong,
                     formSurv=formSurv,
                     dataSurv=dataINLA,
                     dataLongi=dataLongi,
                     id=id,
                     Nsample=Nsample,
                     nproc=nproc,
                     t0=t0,
                     t1=t1, 
                     t2=t2, 
                     t3=t3,
                     idm=idm, 
                     idd=idd,
                     clustertype=clustertype,
                     Ypredmethod=Ypredmethod,
                     NtimePoints=NtimePoints)
      
      
    }else{
      
      
      
      dataJM<-data.frame(ID=data[,colnames(data)%in%id],
                         status=ifelse(idm==1,1,
                                       ifelse(iddCR==1,2,0)),
                         TimeCR=TimeCR)
      colnames(dataJM)[1]<-id
      
      
      dataJM$status<-as.factor(dataJM$status)
      dataJM<-JMbayes2::crisk_setup(dataJM, statusVar = "status", censLevel = "0", 
                                    nameStrata = "CR")
      
      
      
      CoxFit_CR <- coxph(Surv(TimeCR, status2) ~ strata(CR),
                         data = dataJM)
      
      
      
      if(missing(Nsample)){
        Nsample<-100}else{
          if(!inherits(Nsample,c("integer","numeric")))stop("Nsample need to be a numeric or integer")
          if(length(Nsample)!=1)stop("Nsample need to be a numeric or integer")
        }
      
      dataY<-JMidm(timeVar = timeVar,
                   functional_forms = CR_forms, 
                   truncated=truncated,
                   formLong=formLong,
                   formSurv=CoxFit_CR ,
                   dataSurv=dataJM,
                   dataLongi=dataLongi,
                   id=id,
                   Nsample=Nsample,
                   n_iter=methodJM$n_iter,
                   n_burnin=methodJM$n_burnin,
                   n_thin=methodJM$n_thin,
                   n_chain=methodJM$n_chain,
                   nproc=nproc,
                   t0=t0,
                   t1=t1, 
                   t2=t2, 
                   t3=t3,
                   idm=idm, 
                   idd=idd,
                   clustertype=clustertype,
                   Ypredmethod=Ypredmethod,
                   NtimePoints=NtimePoints,
                   seed=methodJM$seed)
      
    }
  res<-list(Y=dataY,
       method=methodlongi,
       runtime=proc.time()-ptm)
  

  class(res) <- "predYidm"
  return(res)
  }