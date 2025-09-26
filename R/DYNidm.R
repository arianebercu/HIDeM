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
#' @param maxiter Maximum number of iterations. The default is 200.
#' @param maxiter.pena Maximum number of iterations for penalised coefficients
#' @param eps A vector of 3 integers >0 used to define the power of
#' three convergence criteria: 1. for the regression parameters,
#' 2. for the likelihood, 3. for the second derivatives. The default
#' is \code{c(5,5,3)} which is translated into convergence if the
#' respective values change less then \eqn{10^{-5}} (for regression
#' parameters and likelihood) and \eqn{10^{-3}} for the second
#' derivatives between two iterations.
#' @param n.knots For \code{method="splines"} only, a vector of length
#' 3 specifing the number of knots, one for each transition, for the
#' M-splines estimate of the baseline intensities in the order \code{0
#' --> 1}, \code{0 --> 2}, \code{1 --> 2}. The default is c(7,7,7). When \code{knots}
#' are specified as a list this argument is ignored.
#' The algorithm needs least 5 knots and at most 20 knots.
#' @param knots Argument only active for the penalized likelihood approach \code{method="Splines"}.
#' There are three ways to control the placement of the knots between the smallest and the largest
#' of all time points:
#' \describe{
#'  \item{\code{knots="equidistant"}}{Knots are placed with same distance on the time scale.}
#'  \item{\code{knots="quantile"}}{Knots are placed such that the number of observations is roughly the same between knots.}
#' \item{knots=list()}{List of 1 or 2 or three vectors. The list elements are the actual placements
#' (timepoints) of the knots for the M-spline. The list may contain
#' one vector of placements for each transition in the order \code{0 --> 1}, \code{0 --> 2}, \code{1 --> 2}.
#' If only vector is specifified the knots are used for all transitions. If only 2 vectors are specifified, the
#' knots for the \code{0 --> 1} transition are also used for the \code{1 --> 2} transition.}
#' }
#' The algorithm needs at least 3 knots in spline and allows no more than 20 knots.
#' @param type.quantile Argument only active for the likelihood approach \code{method="splines"}.
#' There are three ways to control the placement of the knots  according to the time considered between states :=
#' \describe{
#'  \item{\code{type.quantile=1}}{Time for \code{0 --> 1} is the imputed to the middle of the interval left and right for demence . Time for \code{0 --> 2}
#'  and \code{1 --> 2} is the same t, time of news. }
#'  \item{\code{type.quantile=2}}{Time for \code{0 --> 1} is the imputed to the middle of the interval left and right. Time for \code{0 --> 2}
#'  and \code{1 --> 2} is the same t, time of news. }
#' \item{\code{type.quantile=3}}{Time for \code{0 --> 1} is the imputed to the middle of the interval left and right. Time for \code{0 --> 2}
#'  is time of death for non demented sujects only. Time for \code{1 --> 2} is time of death for suject diagnose with dementia. }
#' \item{\code{type.quantile=4}}{Time for \code{0 --> 1} is left and right. Time for \code{0 --> 2}
#'  is time of death for non demented sujects only. Time for \code{1 --> 2} is time of death for suject diagnose with dementia. }
#' }
#' @param B vector of size the number of parameters, in the following order, first the parameters of splines \code{0 --> 1}, \code{0 --> 2}, \code{1 --> 2},
#' second the parameters of explanatory variables in order  \code{0 --> 1}, \code{0 --> 2}, \code{1 --> 2}.
#' This argument is only used for models with M-splines.
#' @param method type of estimation method: "splines" for a likelihood approach with approximation of the transition
#' intensities by M-splines, "Weib" for a parametric approach with a
#' Weibull distribution on the transition intensities. Default is
#' "Weib".
#' @param na.action how NAs are treated. The default is first, any
#' na.action attribute of data, second a na.action setting of options,
#' and third 'na.fail' if that is unset. The 'factory-fresh' default
#' is na.omit. Another possible value is NULL.
#' @param scale.X do you want to center and reduce your explanatory variables
#' @param posfix index of fixed parameters 
#' @param lambda01 Lambda on transition 0 --> 1
#' @param lambda02 Lambda on transition 0 --> 2
#' @param lambda12 Lambda on transition 1 --> 2
#' @param nlambda01 number of Lambda on transition 0 --> 1
#' @param nlambda02 number of Lambda on transition 0 --> 2
#' @param nlambda12 number of Lambda on transition 1 --> 2
#' @param alpha alpha on all transitions 
#' @param penalty which penalty to consider
#' @param penalty.factor which variable should be penalised
#' computed for the Newton-Raphson path for the penalised regression parameter. If FALSE, the 
#' complete hessian is computed. 
#' @param partialH True if should use only diagonal terms of the derivatives of beta
#' @param clustertype in which cluster to work
#' @param nproc number of cluster
#' @param envir working environment 
#' @return
#' \item{call}{the call that produced the result.} \item{coef}{regression
#' parameters.} \item{loglik}{vector containing the log-likelihood and
#' the penalised log-likelihood} \item{cv}{vector containing the convergence criteria based on 
#' stability of parameters, log-likelihood (or penalised) and relative distance to minimum/maximum}
#' \item{niter}{number of iterations.} \item{converged}{integer equal to 1 when
#' the model converged, 2, 3 or 4 otherwise.} \item{modelPar}{Weibull
#' parameters.} \item{N}{number of subjects.} \item{events1}{number of events 0
#' --> 1.} \item{events2}{number of events 0 --> 2 or 0 --> 1 --> 2.}
#' \item{NC}{vector containing the number of covariates on transitions 0 --> 1,
#' 0 --> 2, 1 --> 2.} \item{responseTrans}{model response for the 0 --> 1
#' transition. \code{Hist} or \code{Surv} object.} \item{responseAbs}{model
#' response for the 0 --> 2 transition. \code{Hist} or \code{Surv} object.}
#' \item{time}{times for which transition intensities have been evaluated for
#' plotting.} \item{maxtime}{times of last follow-up or event} \item{mintime}{times of first follow-up or event}
#'  \item{HR}{vector of hazard risks.}
#' \item{V}{variance-covariance matrix derived from the Hessian of the log-likelihood or penalised
#' log-likelihood}\item{Xnames01}{names of covariates on 0 --> 1.}
#' \item{Xnames02}{names of covariates on 0 --> 2.} \item{Xnames12}{names of
#' covariates on 1 --> 2.} \item{knots01}{knots to approximate by M-splines the
#' intensity of the 0 --> 1 transition.} \item{knots02}{knots to approximate by
#' M-splines the intensity of the 0 --> 2 transition.} \item{knots12}{knots to
#' approximate by M-splines the intensity of the 1 --> 2 transition.}
#' \item{nknots01}{number of knots on transition 0 --> 1.}
#' \item{nknots02}{number of knots on transition 0 --> 2.}
#' \item{nknots12}{number of knots on transition 1 --> 2.}
#' \item{theta01}{square root of splines coefficients for transition 0 --> 1.}
#' \item{theta02}{square root of splines coefficients for transition 0 --> 2.}
#' \item{theta12}{square root of splines coefficients for transition 1 --> 2.}
#' \item{alpha}{penalty parameter alpha on all transitions}
#' \item{lambda}{matrix of lambda penalty parameters, first line for 0 --> 1, second 
#' line for 0 --> 2 and third line for 1 --> 2}
#' \item{BIC}{Bayesian Information Criterion} 
#' \item{GCV}{Generalised Cross-Validation approximation}
#' \item{levels}{a list containing the type of the variable on all transitions and its level, 
#' useful for prediction on new-dataset.}
#' \item{runtime}{running time in seconds of the function.}
#' @seealso \code{\link{print.idm}}
#' \code{\link{summary.idm}}
#' \code{\link{predict.idm}}
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
DYNidm <- function(formula01,
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
                
                method="Weib",
                scale.X=T,
                maxiter=100,
                maxiter.pena=10,
                eps=c(5,5,3),
                n.knots=NULL,
                knots="equidistant",
                type.quantile=1,
                na.action = na.fail,
                B=NULL,
                posfix=NULL,
                gausspoint=10,
                lambda01=NULL,
                lambda02=NULL,
                lambda12=NULL,
                nlambda01=50,
                nlambda02=50,
                nlambda12=50,
                penalty=NULL,
                penalty.factor=NULL,
                alpha=ifelse(penalty=="scad",3.7,
                             ifelse(penalty=="mcp",3,
                                    ifelse(penalty%in%c("elasticnet"),0.5,1))),
                nproc=1,
                analytics=T,
                partialH=F,
                clustertype="FORK",
                Ypredmethod="gauss",
                NtimePoints=250,
                dataY=NULL,
                envir=parent.frame()){
  


    # {{{ check formula
  #################################################################################
    ################################ check what has been given ##################
  #################################################################################
  
    call <- match.call()
    ptm <- proc.time()
    if(is.null(dataY)){
    if(missing(methodJM) & missing(methodINLA)){
      stop("Need to specify the method used for time-depend covariates by giving methodINLA or methodJM")
    }
    if(!missing(methodINLA)){
      methodlongi<-"INLA"
      print("Time-dependent covariates estimation will be performed using INLA")
    }else{
      methodlongi<-"JM"
      print("Time-dependent covariates estimation will be performed using JMBayes2")
    }}else{
      if(!inherits(dataY,"predYidm"))stop("The argument dataY must be a class 'predYidm' obtained using DYNpredY.")
      methodlongi<-dataY$method
      dataY<-as.data.frame(dataY$Y)
    }
    if(missing(formula01))stop("Argument formula01 is missing.")
    if(missing(formula02))stop("Argument formula02 is missing.")
    if(!inherits(formula01,"formula"))stop("The argument formula01 must be a class 'formula'.")
    if(!inherits(formula02,"formula"))stop("The argument formula02 must be a class 'formula'.")
    
    if(!method%in%c("Weib","splines"))stop("The argument method needs to be either splines or Weib")
    if(!Ypredmethod%in%c("gauss","equi"))stop("The argument Ypredmethod needs to be either gauss or equi")
    
    if(length(analytics)!=1)stop("The argument analytics must be either T or F")
    if(!analytics%in%c(T,F))stop("The argument analytics must be either T or F")
    ## if(missing(formula02)) formula02 <- formula01
    if(missing(formula12)) formula12 <- formula02
    # }}}
    # {{{ evaluate formula in data
    if(missing(data)) stop("Need a data frame.")
    if(sum(is.na(data))>0)stop("Need a survival data frame with no missing data.")

    if(!inherits(data,"data.frame"))stop("Argument 'data' must be a data.frame")
    
    if(is.null(dataY)){
    if(missing(dataLongi)) stop("Need a dataLongi frame.")
    if(sum(is.na(dataLongi))>0)stop("Need a longitudinal data frame with no missing data.")
    
    if(!inherits(dataLongi,"data.frame"))stop("Argument 'dataLongi' must be a data.frame")}
 
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

    if(is.null(dataY)){
    if(anyNA(dataLongi)){
      variables=unique(c(all.vars(formula01),all.vars(formula02),all.vars(formula12)))
      dataLongi=dataLongi[,variables]
      dataLongi=na.omit(dataLongi)
     # m01[[2]] <- m02[[2]] <- m12[[2]] <- data
    }
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
    #####################   extract covariates   ####################################
    #################################################################################
    
    x01 <- model.matrix(formula01,data=m01)[, -1, drop = FALSE]
    NC01 <- NCOL(x01)

    if (NC01>0)
        Xnames01 <- colnames(x01)
    else
        Xnames01 <- NULL
    ## formula02
    x02 <- model.matrix(formula02,data=m02)[, -1, drop = FALSE]
    NC02 <- NCOL(x02)


    if (NC02>0)
        Xnames02 <- colnames(x02)
    else
        Xnames02 <- NULL
    
    x12 <- model.matrix(formula12,data=m12)[, -1, drop = FALSE]
    NC12 <- NCOL(x12)


    if (NC12>0){
      Xnames12 <- colnames(x12) 
    }else{
     
        Xnames12 <- NULL}


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
    

      
      if(!inherits(maxiter,c("numeric","integer"))|(maxiter!=floor(maxiter)))stop("Maxiter has to be an integer.")
    if(maxiter<0)stop("Maxiter has to be an integer greater or equal to 0.")

    
    if(!inherits(NtimePoints,c("numeric","integer"))|(NtimePoints!=floor(NtimePoints)))stop("NtimePoints has to be an integer greater than 100.")
    if(NtimePoints<100)stop("NtimePoints has to be an integer greater than 100.")
    
      if(!inherits(nproc,c("numeric","integer"))|(nproc!=floor(nproc)))stop("nproc has to be an integer.")
    
      # nbr of quadrature points for estimating integral in idm without penalisation
      if(!gausspoint%in%c(10,15,21,31,41,51,61))stop("Argument type.quantile has to a numeric : 10, 15, 21, 31, 51 or 61.")
      


    #################################################################################
    ########################### Define number of parameters per transitions #########
    #################################################################################
    
    size1 <- NC01 + NC02 + NC12
    

    noVar<-c(ifelse(as.integer(NC01)>0,0,1),
             ifelse(as.integer(NC02)>0,0,1),
             ifelse(as.integer(NC12)>0,0,1))


    nvat01 <- ifelse(noVar[1]==1,0,NC01)
    nvat02 <- ifelse(noVar[2]==1,0,NC02)
    nvat12 <- ifelse(noVar[3]==1,0,NC12)
    
    dimnva01<-ifelse(nvat01==0,1,nvat01)
    dimnva02<-ifelse(nvat02==0,1,nvat02)
    dimnva12<-ifelse(nvat12==0,1,nvat12)
    
    
    NC<-c(NC01,NC02,NC12)
    
    if(noVar[1]==1){ve01<-as.double(rep(0,N))}else{ve01<-as.double(x01)}
    if(noVar[2]==1){ve02<-as.double(rep(0,N))}else{ve02<-as.double(x02)}
    if(noVar[3]==1){ve12<-as.double(rep(0,N))}else{ve12<-as.double(x12)}


    #################################################################################
    ############################# Convergence criterias #############################
    #################################################################################
    epsa<-0.1^eps[1]
    epsb<-0.1^eps[2]
    epsd<-0.1^eps[3]
    eps.eigen<-0.1^2
    
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
    

    size1<-size1+p01+p02+p12
    
  
    if(method=="splines"){
      
      if (is.character(knots)){
        
        
        if(is.null(n.knots)) n.knots<-c(3,3,3)
        if ((length(n.knots)>3) || (length(n.knots)<1)) stop("Argument n.knots has to be a vector of at least one positive integer and at most 3 positive integers.")
        if (length(n.knots)==1) n.knots <- c(n.knots,n.knots,n.knots)
        if (length(n.knots)==2) n.knots <- c(n.knots,n.knots[1])
        nknots01 <- n.knots[1]
        nknots02 <- n.knots[2]
        nknots12 <- n.knots[3]
        if((!is.numeric(n.knots) && !is.integer(n.knots)) || (any(n.knots <3)) || (any(n.knots >20))) #AB : relax this condition to 2 if non penalized spline || (any(n.knots < 5))
          stop("Each element of n.knots has to be an integer between 3 and 20. See help(idm).")
        # AB : PBR codage quantile / equidistant only equidistant in Github
        
        if (!knots%in%c("quantile","equidistant"))stop("Knots need to be either 'equidistant', 'quantile' or directly its values")
        if(!type.quantile%in%c(1,2,3,4))stop("Argument type.quantile has to a numeric : 1, 2, 3 or 4.")
        
        if (knots=="quantile" & type.quantile==1){
          
          approx.illtimes <- (Rtime[idm==1]+Ltime[idm==1])/2
          #approx.illtimes <- Rtime[idm==1]
          knots01 <- quantile(approx.illtimes,seq(0,1,1/(nknots01-1)))
          
          death.time<-responseAbs[responseAbs[,"status"]%in%c(1,2),"time"]
          # Look only at time of events of death when already diagnose
          knots02 <- quantile(death.time,seq(0,1,1/(nknots02-1)))
          knots12 <- quantile(death.time,seq(0,1,1/(nknots12-1)))
        }
        
        if (knots=="quantile" & type.quantile==2){
          approx.illtimes <- (Rtime[idm==1]+Ltime[idm==1])/2
          #approx.illtimes <- Rtime[idm==1]
          knots01 <- quantile(approx.illtimes,seq(0,1,1/(nknots01-1)))
          
          # Look only at time of events of death when already diagnose
          knots02 <- quantile(abstime,seq(0,1,1/(nknots02-1)))
          knots12 <- quantile(abstime,seq(0,1,1/(nknots12-1)))
        }
        
        if (knots=="quantile" & type.quantile==3){
          approx.illtimes <- (Rtime[idm==1] + Ltime[idm==1])/2
          #approx.illtimes <- Rtime[idm==1]
          knots01 <- quantile(approx.illtimes,seq(0,1,1/(nknots01-1)))
          
          # Look only at time of events of death when already diagnose
          # responseTrans = data frame with statut =1 or 2 when dementia
          # responseAbs = data frame with statut =1 or 2 when death
          illdeathtimes <- responseAbs[responseTrans[,"status"]%in%c(1,2) & responseAbs[,"status"]%in%c(1,2),"time"]
          knots12 <- quantile(illdeathtimes,seq(0,1,1/(nknots12-1)))
          
          # Look only at time of events of death when not diagnose
          deathtimes <- responseAbs[responseTrans[,"status"]==0 & responseAbs[,"status"]%in%c(1,2),"time"]
          knots02 <- quantile(deathtimes,seq(0,1,1/(nknots02-1)))
        }
        
        if (knots=="quantile" & type.quantile==4){
          approx.illtimes <- c(Rtime[idm==1],Ltime[idm==1])
          #approx.illtimes <- Rtime[idm==1]
          knots01 <- quantile(approx.illtimes,seq(0,1,1/(nknots01-1)))
          
          # Look only at time of events of death when already diagnose
          illdeathtimes <- responseAbs[responseTrans[,"status"]%in%c(1,2) & responseAbs[,"status"]%in%c(1,2),"time"]
          knots12 <- quantile(illdeathtimes,seq(0,1,1/(nknots12-1)))
          
          # Look only at time of events of death when not diagnose
          deathtimes <- responseAbs[responseTrans[,"status"]==0 & responseAbs[,"status"]%in%c(1,2),"time"]
          knots02 <- quantile(deathtimes,seq(0,1,1/(nknots02-1)))
        }
        if(knots=="equidistant"){
          warning("Unknown specification of knots. Fall back to equidistant.")
          knots01 <- seq(amin,amax,(amax-amin)/(nknots01-1))
          knots02 <- seq(amin,amax,(amax-amin)/(nknots02-1))
          knots12 <- seq(amin,amax,(amax-amin)/(nknots12-1))}
        
        
      }else{## user specified knots
        if (!is.list(knots) || length(knots)==1)
          knots <- list(knots,knots,knots)
        if (length(knots)==2) ## re-use knots from 0->1 for 1->2
          knots <- c(knots,knots[1])
        if (!all(sapply(knots,is.numeric)))
          stop("Incorrect form of argument knots. See help(idm).")
        knots01 <- sort(knots[[1]])
        knots02 <- sort(knots[[2]])
        knots12 <- sort(knots[[3]])
        if (!is.null(knots01)){if(knots01[1]< amin) stop(paste("Transition 0->1: Smallest knot should not be smaller than the time point:",amin))}
        if (!is.null(knots01)){if (knots01[length(knots01)]> amax) stop(paste("Transition 0->1: Largest knot should not be larger than the time point:",amax))}
        if (!is.null(knots02)){if (knots02[1]< amin) stop(paste("Transition 0->2: Smallest knot should not be smaller than the time point:",amin))}
        if (!is.null(knots02)){if (knots02[length(knots02)]> amax) stop(paste("Transition 0->2: Largest knot should not be larger than the time point:",amax))}
        if (!is.null(knots12)){if (knots12[1]< amin) stop(paste("Transition 1->2: Smallest knot should not be smaller than the time point:",amin))}
        if (!is.null(knots12)){if (knots12[length(knots12)]> amax) stop(paste("Transition 1->2: Largest knot should not be larger than the time point:",amax))}
        ## FIXME: check if knots within amin, amax
        ## if (knots01[[1]] < amin) stop("Smallest knot ")
        nknots01 <- length(knots01)
        nknots02 <- length(knots02)
        nknots12 <- length(knots12)
        
        
      }

      if (any(c(nknots01,nknots02,nknots12)>20)){
        stop("Cannot handle more than 20 knots.")
      }
      
      ## AB : add min and max
      ## 0 -- 1
      if (min(knots01)>amin) knots01[1] <- amin
      if (max(knots01)<amax) knots01[length(knots01)] <- amax
      ## 0 -- 2
      if (min(knots02)>amin) knots02[1] <- amin
      if (max(knots02)<amax) knots02[length(knots02)] <- amax
      ## 1 -- 2
      if (min(knots12)>amin) knots12[1] <- amin
      if (max(knots12)<amax) knots12[length(knots12)] <- amax
      
      
      ## make fake knots needed for M-splines
      knots01 <- c(rep(knots01[1],3),knots01,rep(knots01[length(knots01)],3))
      knots02 <- c(rep(knots02[1],3),knots02,rep(knots02[length(knots02)],3))
      knots12 <- c(rep(knots12[1],3),knots12,rep(knots12[length(knots12)],3))
      size_spline<-nknots01+nknots02+nknots12 + 6
      size_V <- size1 + size_spline
      
      # if B is step put b at B otherwise 0.5 for splines and 0 for parameters
      # initiate parameters values for method=splines
      if (is.null(B)){
        
        b<-c(rep(0.5,size_spline),rep(0,size1))
        # define splines and initialize them
      }else{
        if(!inherits(B,"numeric")){stop(paste0("B need to be a numeric"))}
        # if(any(B[1:size_spline]<0)){stop(paste0("B need to be positive for spline parameters"))}
        if(length(B)!=(size_V)){stop(paste0("The length of the initialization must be : ",size_V))}
        b<-B}
      
      
    }

   
    ############################################################################
    #################### defines initiate values with weibull ##################
    ############################################################################
    
    if(method=="Weib"){
      
      size2 <- size1^2
      size_V <- size1 + 6
      # if B is step put b at B otherwise defined by events in data for
      # weibull parameters and 0 for beta
      ts<-sum(abstime-t0)
      
      if (!is.null(B)){
        
        if(!inherits(B,"numeric")){stop(paste0("B need to be a numeric"))}
        # if(any(B[1:size_spline]<0)){stop(paste0("B need to be positive for spline parameters"))}
        if(length(B)!=(size_V)){stop(paste0("The length of the initialization must be : ",size_V))}
        b<-B}else{
        
          b<-c(1,sqrt(sum(idm)/ts),1,sqrt(sum(idd)/ts),1,sqrt(sum(idd)/ts),rep(0,size_V-6))
          
        }
    }
    
    
    ############################################################################
    #################### check if we have penalty parameters  ##################
    ############################################################################
    
    fix0<-rep(0,size_V)
    if(is.null(penalty)){
      penalty<-"none"}
    if(!penalty%in%c("none","lasso","ridge","elasticnet","mcp","scad")){
      stop(paste0("Parameter penalty must be either : lasso, ridge, elasticnet, mcp or scad"))}
    
    
    if(!is.null(posfix)){
      
      if(!inherits(posfix,c("numeric","integer"))){stop(paste0("Posfix need to be a numeric"))}
      posfix<-na.omit(posfix)
      if(min(posfix)<=0){stop(paste0("The indexation of posfix need to start at 1"))}
      if(max(posfix)>(size_V)){stop(paste0("The indexation of posfix cannot exceed the number of parameters : ",size_V))}
      if(length(posfix)==(size_V)){stop(paste0("At least one parameter need to be non-fixed"))}
      
      # adapt to penalised setting :
      #if(any(posfix>size_spline) & penalty==T){stop(paste0("Fixed parameters can only be on spline when penalised model is request "))}
      if(length(posfix)==6 & penalty%in%c("lasso","ridge","elasticnet","mcp","scad") & !any(posfix>6)){
        stop(paste0("All weibull parameters cannot be fixed when penalised model is request "))}
      
      fix0[posfix]<-1
      
      
    }
    

    fit <- NULL
    ############################################################################
    ######################## Start algorithm to maximise #######################
    ########################       log-likelihodd        #######################
    ############################################################################
    
    
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
    
if(is.null(dataY)){
    if(methodlongi=="INLA"){
      
      
     
      
      
      if(truncated==1){
        formSurv<-list(inla.surv(time=TimeCR, event=idm,truncation=t0) ~ -1,
                       inla.surv(time=TimeCR, event=iddCR,truncation=t0) ~ -1)
      }else{
        formSurv<-list(inla.surv(time=TimeCR, event=idm) ~ -1,
                       inla.surv(time=TimeCR, event=iddCR) ~ -1)
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
}
   
                   
    
    ############################# scale explanatory variables #####################
    if(scale.X==T){
      # to know which variable to center and reduces : 
      
      if(nvat01>0){
        names<-as.vector(unlist(lapply(m01, class)))[-1]
        
        # standardize variable : 
        # we need to center and reduce variables when penalty is done 
        namesx01<-as.vector(unlist(lapply(m01[,colnames(m01)%in%Xnames01], class)))[-1]
        scaledx01<-scale(x01[,which(namesx01=="numeric")])
        # keep center and reduces parameters 
        xm01<-attr(scaledx01,"scaled:center")
        xs01<-attr(scaledx01,"scaled:scale")
        x01[,which(namesx01=="numeric")]<-scaledx01
        ve01<-as.double(x01)
      }
      
      if(nvat02>0){
        
        # to know wich variable to center and reduces : 
        names<-as.vector(unlist(lapply(m02, class)))[-1]
        
        # standardize variable : 
        
        # we need to center and reduce variables when penalty is done 
        namesx02<-as.vector(unlist(lapply(m02[,colnames(m02)%in%Xnames02], class)))[-1]
        scaledx02<-scale(x02[,which(namesx02=="numeric")])
        # keep center and reduces parameters 
        xm02<-attr(scaledx02,"scaled:center")
        xs02<-attr(scaledx02,"scaled:scale")
        x02[,which(namesx02=="numeric")]<-scaledx02
        ve02<-as.double(x02)
      }
      
      if(nvat12>0){
        
        # to know wich variable to center and reduces : 
        # no need [-1] as no element before ~
        names<-as.vector(unlist(lapply(m12, class)))
        
        # standardize variable : 
        
        # we need to center and reduce variables when penalty is done 
        namesx12<-as.vector(unlist(lapply(m12[,colnames(m12)%in%Xnames12], class)))
        scaledx12<-scale(x12[,which(namesx12=="numeric")])
        # keep center and reduces parameters 
        xm12<-attr(scaledx12,"scaled:center")
        xs12<-attr(scaledx12,"scaled:scale")
        x12[,which(namesx12=="numeric")]<-scaledx12
        ve12<-as.double(x12)
      }
      
      
     ym <- do.call(rbind, lapply(split(dataY[, 4:dim(dataY)[2]], dataY$Outcome), function(x) {
        apply(x, 2, mean)
      }))
      
      # sds per group
      ys <- do.call(rbind, lapply(split(dataY[, 4:dim(dataY)[2]], dataY$Outcome), function(x) {
        apply(x, 2, sd)
      }))
      
      dataY[, 4:dim(dataY)[2]] <- do.call(rbind,
                                  lapply(split( dataY[, 4:dim(dataY)[2]], dataY$Outcome), function(x) {
                                    apply(x, 2, function(col) (col - min(col)) / (max(col) - min(col)))
                                  })
      )
      
      
    }
    
    if(penalty=="none"){
 
######################### with M-spline baseline risk  #########################
      if(method=="splines"){
        
        if(sum(fix0)>0){
          bfix<-b[fix0==1]
          b<-b[fix0==0]
        }else{bfix<-1}
        
    
        out<-DYNidm.splines(b=b,
                         clustertype=clustertype,
                         epsa=epsa,
                         epsb=epsb,
                         epsd=epsd,
                         nproc=nproc,
                         maxiter=maxiter,
                         size_V=size_V,
                         size_spline=size_spline,
                         noVar=noVar,
                         bfix=bfix,
                         fix0=fix0,
                         knots01=knots01,
                         knots02=knots02,
                         knots12=knots12,
                         ctime=ctime,
                         N=N,
                         nknots01=nknots01,
                         nknots02=nknots02,
                         nknots12=nknots12,
                         ve01=ve01,
                         ve02=ve02,
                         ve12=ve12,
                         dimnva01=dimnva01,
                         dimnva02=dimnva02,
                         dimnva12=dimnva12,
                         nvat01=nvat01,
                         nvat02=nvat02,
                         nvat12=nvat12,
                         t0=t0,
                         t1=t1,
                         t2=t2,
                         t3=t3,
                         troncature=troncature,
                         
                         dataY=dataY,
                         Longitransition=Longitransition,
                         NtimePoints=NtimePoints,
                         Ypredmethod=Ypredmethod,
                         timeVar=timeVar,
                         ynames=ynames,
                         id=id,
                         outcome=outcome,
                         outcome01=outcome01,
                         outcome02=outcome02,
                         outcome12=outcome12,
                         p01=p01,p02=p02,p12=p12,
                         dimp01=dimp01,
                         dimp02=dimp02,
                         dimp12=dimp12)
      }
  
      

######################### with weibull baseline risk ###########################
      
      if(method=="Weib"){

   
        #save(ctime,file="testctime.RData")
        
        out <- DYNidm.weib(b=b,
                                    fix0=fix0,
                                    size_V=size_V,
                                    clustertype=clustertype,
                                    epsa=epsa,
                                    epsb=epsb,
                                    epsd=epsd,
                                    nproc=nproc,
                                    maxiter=maxiter,
                                    ctime=ctime,
                                    N=N,
                                    ve01=ve01,
                                    ve02=ve02,
                                    ve12=ve12,
                                    dimnva01=dimnva01,
                                    dimnva02=dimnva02,
                                    dimnva12=dimnva12,
                                    nvat01=nvat01,
                                    nvat02=nvat02,
                                    nvat12=nvat12,
                                    t0=t0,
                                    t1=t1,
                                    t2=t2,
                                    t3=t3,
                                    idm=idm,
                                    idd=idd,
                                    ts=ts,
                                    troncature=troncature,
                           
                        dataY=dataY,
                        Longitransition=Longitransition,
                        NtimePoints=NtimePoints,
                        Ypredmethod=Ypredmethod,
                        timeVar=timeVar,
                        ynames=ynames,
                        id=id,
                        outcome=outcome,
                        outcome01=outcome01,
                        outcome02=outcome02,
                        outcome12=outcome12,
                        p01=p01,p02=p02,p12=p12,
                        dimp01=dimp01,
                        dimp02=dimp02,
                        dimp12=dimp12)
        
        
            
        
        
       
      }
      
      
    }else{

      ############################################################################
      ######################## Start algorithm to maximise #######################
      ########################  penalised log-likelihodd    ######################
      ############################################################################

#################################### penalty check #############################
          if(nvat01==0 & nvat02==0 & nvat12==0 & p01==0 & p02==0 & p12==0)stop("To perform penalisation you need explanatory variables in each transition")
          
          # permits to not penalise on some parameters
          if(is.null(penalty.factor)){
            penalty.factor<-rep(1, nvat01+nvat02+nvat12+p01+p02+p12)
          }else{
            if(any(min(penalty.factor)<0) | any(max(penalty.factor)>1) | any(round(penalty.factor)!=penalty.factor) | length(penalty.factor)!=(nvat01+nvat02+nvat12+p01+p02+p12)){
              stop(paste0("Penalty.factor need to be a vector of 0 and 1 of length : ",nvat01+nvat02+nvat12+p01+p02+p12))
            }
          }


           
############################ set value of penalty parameters ###################
          if(penalty=="lasso"){alpha<-1}
          if(penalty=="ridge"){alpha<-0}
          if(length(alpha)>1)stop("Can only specify one value for alpha")
          if(penalty=="mcp"){
            if(!inherits(alpha,c("numeric","integer"))  | alpha<=1)stop("Alpha need to be a numeric and superior to 1")
            
          }
      
          if(penalty=="scad"){
            if(!inherits(alpha,c("numeric","integer")) | alpha<=2)stop("Alpha need to be a numeric and superior to 2")
            
          }
          if(penalty%in%c("elasticnet")){
           if(!inherits(alpha,c("numeric","integer"))  | alpha>1 | alpha <0)stop("Alpha need to be a numeric between 0 and 1")
          }
          if(!inherits(nlambda01,c("numeric","integer")) | round(nlambda01)!=nlambda01 | nlambda01<1)stop("Nlambda01 need to be an integer superior or equal to 1")
          if(!inherits(nlambda02,c("numeric","integer")) | round(nlambda02)!=nlambda02 | nlambda02<1)stop("Nlambda02 need to be an integer superior or equal to 1")
          if(!inherits(nlambda12,c("numeric","integer")) | round(nlambda12)!=nlambda12 | nlambda12<1)stop("Nlambda12 need to be an integer superior or equal to 1")
          
          pace.lambda<-ifelse(N<size_V,0.05,0.0001)
          
################################################################################
########################## perform penalty algorithm ###########################
##########################   with M-splines baseline risk ######################
################################################################################
          
          Nsample<-dim(dataY)[2]-3
          browser()
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
            
            
          }else{
            
            # check if predictions could be performed 
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
            
          }
          
          if(method=="splines"){
            # if user did not specified the lambda values 
            if(is.null(lambda01)|is.null(lambda02)|is.null(lambda12)){
              
              if(nproc>1){
                if(is.null(clustertype)){
                  clustpar <- parallel::makeCluster(nproc)#, outfile="")
                }
                else{
                  clustpar <- parallel::makeCluster(nproc, type=clustertype)#, outfile="")
                }
                
                doParallel::registerDoParallel(clustpar)
                
              }
             # calculate derivatives when all beta=0 to 
             # estimate range of values for lambda
              output<-deriva.gradient(b=c(b[1:size_spline],rep(0,size_V-size_spline)),
                                  nproc=nproc,
                                  funcpa=idmlLikelihood,
                                  npm=size_V,
                                  npar=size_V,
                                  bfix=1,
                                  fix=rep(0,size_V),
                                  zi01=knots01,
                                  zi02=knots02,
                                  zi12=knots12,
                                  ctime=ctime,
                                  no=N,
                                  nz01=nknots01,
                                  nz02=nknots02,
                                  nz12=nknots12,
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
                                  gausspoint=gausspoint)
              
              if(nproc>1){parallel::stopCluster(clustpar)}
              ## what should we do if max(output$v) == 0
              # take maximum of absolute derivatives
              if(penalty%in%c("ridge","mcp","scad")){
                #lambda.max<-ifelse(max(output$v)==0,0.001,max(output$v))
                lambda.max<-ifelse(max(abs(output$v))==0,0.001,max(abs(output$v)))
              }else{
                lambda.max<-ifelse(max(abs(output$v))==0,0.001,max(abs(output$v))/alpha)
              }
            }
            
            if(nvat01>0){
            if(!is.null(lambda01)){
              nlambda01<-length(lambda01)
              if(length(lambda01)<1)stop("Penalisation can be performed for at least one lambda01 ")
              if(min(lambda01)<=0)stop("Lambda01 must be composed of strictly positive values ")
            }else{
              
              lambda01<-lambda.max*((pace.lambda)^(c(1:nlambda01)/nlambda01))
            }
            }else{lambda01<-0.0001}
            
            if(nvat02>0){
            if(!is.null(lambda02)){
              nlambda02<-length(lambda02)
              if(length(lambda02)<1)stop("Penalisation can be performed for at least one lambda02 ")
              if(min(lambda02)<=0)stop("Lambda02 must be composed of strictly positive values ")
            }else{
              lambda02<-lambda.max*((pace.lambda)^(c(1:nlambda02)/nlambda02))
            }
            }else{lambda02<-0.0001}
            
            if(nvat12>0){
            if(!is.null(lambda12)){
              nlambda12<-length(lambda12)
              if(length(lambda12)<1)stop("Penalisation can be performed for at least one lambda12 ")
              if(min(lambda12)<=0)stop("Lambda12 must be composed of strictly positive values ")
            }else{
              lambda12<-lambda.max*((pace.lambda)^(c(1:nlambda12)/nlambda12))
            }
            }else{lambda12<-0.0001}
            
################################################################################
########################## perform penalty algorithm ###########################
##########################   with M-splines baseline risk ######################
################################################################################
            
            out<-idm.penalty.splines(b=b,
                             fix0=fix0,
                             size_V=size_V,
                             size_spline=size_spline,
                             clustertype=clustertype,
                             epsa=epsa,
                             epsb=epsb,
                             epsd=epsd,
                             eps.eigen=eps.eigen,
                             nproc=nproc,
                             maxiter=maxiter,
                             maxiter.pena=maxiter.pena,
                             knots01=knots01,
                             knots02=knots02,
                             knots12=knots12,
                             ctime=ctime,
                             N=N,
                             nknots01=nknots01,
                             nknots02=nknots02,
                             nknots12=nknots12,
                             ve01=ve01,
                             ve02=ve02,
                             ve12=ve12,
                             dimnva01=dimnva01,
                             dimnva02=dimnva02,
                             dimnva12=dimnva12,
                             nvat01=nvat01,
                             nvat02=nvat02,
                             nvat12=nvat12,
                             t0=t0,
                             t1=t1,
                             t2=t2,
                             t3=t3,
                             troncature=troncature,
                             gausspoint=gausspoint,
                             nlambda01=nlambda01,
                             lambda01=lambda01,
                             nlambda02=nlambda02,
                             lambda02=lambda02,
                             nlambda12=nlambda12,
                             lambda12=lambda12,
                             alpha=alpha,
                             penalty.factor=penalty.factor,
                             penalty=penalty,
                                    analytics=analytics,
                             partialH=partialH)
            
############################## Output   ########################################
############################## on beta and HR   ################################
            
           
            lambda<-out$lambda
            alpha<-out$alpha
            
            H<-out$H
            
            beta<-as.matrix(out$b)
            fix<-as.matrix(out$fix)
            lambda<-as.matrix(lambda)
            
            theta_names <- cbind(c(rep("theta01",(nknots01+2)),rep("theta02",(nknots02+2)),rep("theta12",(nknots12+2))),c((1:(nknots01+2)),(1:(nknots02+2)),(1:(nknots12+2))))
            theta_names <- as.vector(apply(theta_names,1,paste,collapse=" "))
            rownames(beta) <-c(theta_names,c(Xnames01,Xnames02,Xnames12))
            rownames(fix)<-c(theta_names,c(Xnames01,Xnames02,Xnames12))
            rownames(lambda)<-c("lambda01","lambda02","lambda12")
            
            theta01<-beta[1:(nknots01+2),]
            theta02<-beta[(nknots01+3):(nknots01+nknots02+4),]
            theta12<-beta[(nknots01+nknots02+5):(nknots01+nknots12+nknots02+6),]
            betaCoef <- beta[(size_spline+1):size_V,]
            betaCoef<-as.matrix(betaCoef)
            fit$coef <- betaCoef
            fit$HR <- exp(betaCoef)
            
########################## on  BIC and GCV #####################################
            if(dim(beta)[2]>1){
              fit$BIC<--2*out$fn.value+log(N)*colSums(beta[(size_spline+1):size_V,]!=0)
            }else{fit$BIC<--2*out$fn.value+log(N)*sum(beta[(size_spline+1):size_V,]!=0)}
            
            
            fit$GCV<-rep(NA,dim(lambda)[2])
            
########################## on  covariance matrix ###############################           
            npm<-sum(fix0[(size_spline+1):size_V]==0)
            npm01<-ifelse(nvat01>0,sum(fix0[(size_spline+1):(size_spline+nvat01)]==0),0)
            npm02<-ifelse(nvat02>0,sum(fix0[(size_spline+nvat01+1):(size_spline+nvat01+nvat02)]==0),0)
            npm12<-ifelse(nvat12>0,sum(fix0[(size_spline+nvat01+nvat02+1):size_V]==0),0)
            
            V<-matrix(0,ncol=dim(lambda)[2]*npm,nrow=npm)
          
            for(i in 1:dim(lambda)[2]){
              id.keep<-which(betaCoef[fix0[(size_spline+1):size_V]==0,i]!=0)
              if(length(id.keep)==0){
                fit$GCV[i]<--1/N*out$fn.value[i]
              }else{
                lambda.matrix<-matrix(0,npm,npm)
                if(nvat01>0){
                  
                  if(penalty%in%c("none","lasso","ridge","elasticnet")){
                    diag(lambda.matrix)[1:npm01]<-rep(2*(1-alpha[i])*lambda[1,i],npm01)}
                  
                  if(penalty=="mcp"){
                    idbeta<-which(abs(betaCoef[1:npm01,i])<=alpha[i]*lambda[1,i])
                    diag(lambda.matrix)[1:npm01]<-rep(0,npm01)
                    diag(lambda.matrix)[idbeta]<--1/alpha[i]
                    
                  }
                  if(penalty=="scad"){
                    idbeta<-which((abs(betaCoef[1:npm01,i])>lambda[1,i])&(abs(betaCoef[1:npm01,i])<=lambda[1,i]*alpha[i]))
                    diag(lambda.matrix)[1:npm01]<-rep(0,npm01)
                    diag(lambda.matrix)[idbeta]<--1/(alpha[i]-1)
                    
                  }
                }
                if(nvat02>0){
                  
                  if(penalty%in%c("none","lasso","ridge","elasticnet")){
                    diag(lambda.matrix)[(npm01+1):(npm01+npm02)]<-rep(2*(1-alpha[i])*lambda[2,i],npm02)}
                  
                  if(penalty=="mcp"){
                    
                    idbeta<-which(abs(betaCoef[(npm01+1):(npm01+npm02),i])<=alpha[i]*lambda[2,i])
                    diag(lambda.matrix)[(npm01+1):(npm01+npm02)]<-rep(0,npm02)
                    diag(lambda.matrix)[idbeta+npm01]<--1/alpha[i]
                    
                    
                  }
                  if(penalty=="scad"){
                    
                    idbeta<-which((abs(betaCoef[(npm01+1):(npm01+npm02):npm01,i])>lambda[2,i])&(abs(betaCoef[(npm01+1):(npm01+npm02),i])<=lambda[2,i]*alpha[i]))
                    diag(lambda.matrix)[(npm01+1):(npm01+npm02)]<-rep(0,npm02)
                    diag(lambda.matrix)[idbeta+npm01]<--1/(alpha[i]-1)
                    
                  }
                }
                if(nvat12>0){
                  if(penalty%in%c("none","lasso","ridge","elasticnet")){
                    diag(lambda.matrix)[(npm01+npm02+1):npm]<-rep(2*(1-alpha[i])*lambda[3,i],npm12)}
                  
                  
                  if(penalty=="mcp"){
                    idbeta<-which(abs(betaCoef[(npm01+npm02+1):npm,i])<=alpha[i]*lambda[3,i])
                    diag(lambda.matrix)[(npm01+npm02+1):npm]<-rep(0,npm12)
                    diag(lambda.matrix)[idbeta+npm01+npm02]<--1/alpha[i]
                    
                  }
                  
                  if(penalty=="scad"){
                    idbeta<-which((abs(betaCoef[(npm01+npm02+1):npm,i])>lambda[3,i])&(abs(betaCoef[(npm01+npm02+1):npm,i])<=lambda[3,i]*alpha[i]))
                    diag(lambda.matrix)[(npm01+npm02+1):npm]<-rep(0,npm12)
                    diag(lambda.matrix)[idbeta+npm01+npm02]<--1/(alpha[i]-1)
                    
                  }
                }
                # only beta different from 0 
                lambda.matrix<-lambda.matrix[id.keep,id.keep,drop=FALSE]
                
                if(out$istop[i]==1){
                  H_spec<-out$H[,(1+(i-1)*npm):(i*npm),drop=FALSE]
                  # only beta different from 0 
                  H_spec<-H_spec[id.keep,id.keep,drop=FALSE]
                  if(!(any(is.na(H_spec))|any(H_spec==Inf) |any(H_spec==-Inf))){
                  if(abs(det(H_spec))>1e-10 & kappa(H_spec)<1e10){
                    id.keepV<-(1+(i-1)*npm):(i*npm)
                    id.keepV<-id.keepV[id.keepV%in%((i-1)*npm+id.keep)]
                    V[id.keep,id.keepV]<-solve(H_spec)
                    # maximisation issue thus : 10/04/24
                    #H_pl<-H_spec+lambda.matrix
                    # as mla return v=-second derivatives when maximisation^-1
                    H_pl<-H_spec+lambda.matrix
                    trace_model<-lava::tr(solve(H_pl)%*%H_spec)
                    fit$GCV[i]<--1/N*(out$fn.value[i]-trace_model)}
                  }
                }
                
              }
            }
            cv<-list(ca.beta=out$ca.beta,ca.spline=out$ca.spline,ca.validity=out$ca.validity,
                     cb=out$cb,rdm=NULL)
            
            fit$knots01 <- knots01
            fit$knots02 <- knots02
            fit$knots12 <- knots12
            fit$nknots01 <- nknots01
            fit$nknots02 <- nknots02
            fit$nknots12 <- nknots12
            fit$theta01 <- as.matrix(theta01)
            fit$theta02 <- as.matrix(theta02)
            fit$theta12 <- as.matrix(theta12)

############# on times to do prediction on #################################
            fit$time <- matrix(NA,ncol=3,nrow=100)
            fit$time[,1]<-seq(from=knots01[1],to=knots01[length(knots01)],length.out=100)
            fit$time[,2]<-seq(from=knots02[1],to=knots02[length(knots02)],length.out=100)
            fit$time[,3]<-seq(from=knots12[1],to=knots12[length(knots12)],length.out=100)
            
          }
          
          ################################################################################
          ########################## perform penalty algorithm ###########################
          ##########################   with weibull baseline risk ######################
          ################################################################################
          
          if(method=="Weib"){
          
           
            
            if(is.null(lambda01)|is.null(lambda02)|is.null(lambda12)){
              if(nproc>1){
                if(is.null(clustertype)){
                  clustpar <- parallel::makeCluster(nproc)#, outfile="")
                }
                else{
                  clustpar <- parallel::makeCluster(nproc, type=clustertype)#, outfile="")
                }
                
                doParallel::registerDoParallel(clustpar)
                
              }
              # calculate derivatives at beta=0, to have a range of 
              # value for lambda 
              output<-deriva.gradient(b=c(b[1:6],rep(0,size_V-6)),
                                  nproc=nproc,
                                  funcpa=gaussDYNidmlLikelihoodweib,
                                  npm=size_V,
                                  npar=size_V,
                                  bfix=1,
                                  fix=rep(0,size_V),
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
                                  y01=y01[,4],
                                  y02=y02[,4],
                                  y12=y12[,4],
                                  p01=p01,
                                  p02=p02,
                                  p12=p12,
                                  dimp01=dimp01,
                                  dimp02=dimp02,
                                  dimp12=dimp12,
                                  Ntime=NtimePoints)
              
             
              if(nproc>1){parallel::stopCluster(clustpar)}
              
              
              ## what should we do if max(output$v) == 0
              if(penalty%in%c("ridge","mcp","scad")){
                #lambda.max<-ifelse(max(output$v)==0,0.001,max(output$v))
                lambda.max<-ifelse(max(abs(output$v))==0,0.001,max(abs(output$v)))
              }else{
                lambda.max<-ifelse(max(abs(output$v))==0,0.001,max(abs(output$v))/alpha)
              }
              
            }
            
            
            if(nvat01>0 | p01>0){
            if(!is.null(lambda01)){
              nlambda01<-length(lambda01)
              if(length(lambda01)<1)stop("Penalisation can be performed for at least one lambda01 ")
              if(min(lambda01)<=0)stop("Lambda01 must be composed of strictly positive values ")
            }else{
              
              lambda01<-lambda.max*((pace.lambda)^(c(1:nlambda01)/nlambda01))
            }
            }else{lambda01<-0.0001}
            
            if(nvat02>0 | p02 >0){
            if(!is.null(lambda02)){
              nlambda02<-length(lambda02)
              if(length(lambda02)<1)stop("Penalisation can be performed for at least one lambda02 ")
              if(min(lambda02)<=0)stop("Lambda02 must be composed of strictly positive values ")
            }else{
              
              lambda02<-lambda.max*((pace.lambda)^(c(1:nlambda02)/nlambda02))
            }
            }else{lambda02<-0.0001}
            
            if(nvat12>0 | p12 >0){
            if(!is.null(lambda12)){
              nlambda12<-length(lambda12)
              if(length(lambda12)<1)stop("Penalisation can be performed for at least one lambda12 ")
              if(min(lambda12)<=0)stop("Lambda12 must be composed of strictly positive values ")
            }else{
              
              lambda12<-lambda.max*((pace.lambda)^(c(1:nlambda12)/nlambda12))
            }
            }else{lambda12<-0.0001}
            
            
            
           
              out <- DYNidm.penalty.weib(b=b,
                               fix0=fix0,
                               size_V=size_V,
                               clustertype=clustertype,
                               epsa=epsa,
                               epsb=epsb,
                               epsd=epsd,
                               eps.eigen=eps.eigen,
                               nproc=nproc,
                               maxiter=maxiter,
                               maxiter.pena=maxiter.pena,
                               ctime=ctime,
                               N=N,
                               ve01=ve01,
                               ve02=ve02,
                               ve12=ve12,
                               dimnva01=dimnva01,
                               dimnva02=dimnva02,
                               dimnva12=dimnva12,
                               nvat01=nvat01,
                               nvat02=nvat02,
                               nvat12=nvat12,
                               t0=t0,
                               t1=t1,
                               t2=t2,
                               t3=t3,
                               troncature=troncature,
                               nlambda01=nlambda01,
                               lambda01=lambda01,
                               nlambda02=nlambda02,
                               lambda02=lambda02,
                               nlambda12=nlambda12,
                               lambda12=lambda12,
                               alpha=alpha,
                               penalty.factor=penalty.factor,
                               penalty=penalty,
                               partialH=partialH,
                               y01=y01,
                               y02=y02,
                               y12=y12,
                               NtimePoints=NtimePoints,
                               p01=p01,
                               p02=p02,
                               p12=p12,
                               dimp01=dimp01,
                               dimp02=dimp02,
                               dimp12=dimp12,
                               Nsample=Nsample)
            
              
              
######################### Output ###############################################
######################## on beta and HR ########################################
              
             

        }

}
        return(out)
}
