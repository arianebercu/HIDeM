#' Fit an illness-death model
#'
#' Fit an illness-death model using either a semi-parametric approach
#' (penalized likelihood with an approximation of the transition intensity
#' functions by linear combination of M-splines) or a parametric approach
#' (specifying Weibull distributions on the transition intensities).
#' Left-truncated, right-censored, and interval-censored data are allowed.
#' State 0 corresponds to the initial state, state 1 to the transient one,
#' state 2 to the absorbant one. The allowed transitions are: 0 --> 1, 0 --> 2
#' and 1 --> 2.
#'
#' The estimated parameters are obtained using the robust Marquardt algorithm
#' (Marquardt, 1963) which is a combination between a Newton-Raphson algorithm
#' and a steepest descent algorithm.
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
#' @param data A data frame in which to interpret the variables of
#' \code{formula01}, \code{formula02} and \code{formula12}.
#' @param maxiter Maximum number of iterations. The default is 200.
#' @param maxiter.pena Maximum number of iterations for penalised coefficients
#' @param eps A vector of 3 integers >0 used to define the power of
#' three convergence criteria: 1. for the regression parameters,
#' 2. for the likelihood, 3. for the second derivatives. The default
#' is \code{c(5,5,3)} which is translated into convergence if the
#' respective values change less then \eqn{10^{-5}} (for regression
#' parameters and likelihood) and \eqn{10^{-3}} for the second
#' derivatives between two iterations.
#' @param eps.spline the power of convergence for splines parameters only
#' @param eps.eigen the power of convergence for eigen values of covariance matrix only
#' @param n.knots For \code{method="splines"} only, a vector of length
#' 3 specifing the number of knots, one for each transition, for the
#' M-splines estimate of the baseline intensities in the order \code{0
#' --> 1}, \code{0 --> 2}, \code{1 --> 2}. The default is c(7,7,7). When \code{knots}
#' are specified as a list this argument is ignored.
#' The algorithm needs least 5 knots and at most 20 knots.
#' @param knots Argument only active for the penalized likelihood approach \code{method="Splines"}.
#' There are three ways to control the placement of the knots between the smallest and the largest
#' of all time points:
#' \itemize{
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
#' \itemize{
#' #'  \item{\code{type.quantile=1}}{Time for \code{0 --> 1} is the imputed to the middle of the interval left and right for demence . Time for \code{0 --> 2}
#'  and \code{1 --> 2} is the same t, time of news. }
#'  \item{\code{type.quantile=2}}{Time for \code{0 --> 1} is the imputed to the middle of the interval left and right. Time for \code{0 --> 2}
#'  and \code{1 --> 2} is the same t, time of news. }
#' \item{\code{type.quantile=3}}{Time for \code{0 --> 1} is the imputed to the middle of the interval left and right. Time for \code{0 --> 2}
#'  is time of death for non demented sujects only. Time for \code{1 --> 2} is time of death for suject diagnose with dementia. }
#'  #' \item{\code{type.quantile=4}}{Time for \code{0 --> 1} is left and right. Time for \code{0 --> 2}
#'  is time of death for non demented sujects only. Time for \code{1 --> 2} is time of death for suject diagnose with dementia. }
#' }
#' @param B vector of size the number of parameters, in the following order, first the parameters of splines \code{0 --> 1}, \code{0 --> 2}, \code{1 --> 2},
#' second the parameters of explanatory variables in order  \code{0 --> 1}, \code{0 --> 2}, \code{1 --> 2}.
#' This argument is only used for models with M-splines.
#' @param method type of estimation method: "splines" for a likelihood approach with approximation of the transition
#' intensities by M-splines, "Weib" for a parametric approach with a
#' Weibull distribution on the transition intensities. Default is
#' "Weib".
#' @param subset expression indicating the subset of the rows of data
#' to be used in the fit. All observations are included by default.
#' @param na.action how NAs are treated. The default is first, any
#' na.action attribute of data, second a na.action setting of options,
#' and third 'na.fail' if that is unset. The 'factory-fresh' default
#' is na.omit. Another possible value is NULL.
#' @param scale.X do you want to center and reduce your explanatory variables
#' @param posfix index of fixed parameters 
#' @param gauss.point gauss quadrature points in the approximation of integrals
#' @param lambda01 Lambda on transition 0 --> 1
#' @param lambda02 Lambda on transition 0 --> 2
#' @param lambda12 Lambda on transition 1 --> 2
#' @param nlambda01 number of Lambda on transition 0 --> 1
#' @param nlambda02 number of Lambda on transition 0 --> 2
#' @param nlambda12 number of Lambda on transition 1 --> 2
#' @param alpha alpha on all transitions 
#' @param penalty which penalty to consider
#' @param penalty.factor which variable should be penalised
#' @param step.sequential should we use the optimisation version to fix splines 
#' @param clustertype in which cluster to work
#' @param nproc number of cluster
#' @param option.sequential parameters to give if you want to do the optimisation version to
#'  fix splines
#' @return
#'
#' \item{call}{the call that produced the result.} \item{coef}{regression
#' parameters.} \item{loglik}{vector containing the log-likelihood without and
#' with covariate.} \item{cv}{vector containing the convergence criteria.}
#' \item{niter}{number of iterations.} \item{converged}{integer equal to 1 when
#' the model converged, 2, 3 or 4 otherwise.} \item{modelPar}{Weibull
#' parameters.} \item{N}{number of subjects.} \item{events1}{number of events 0
#' --> 1.} \item{events2}{number of events 0 --> 2 or 0 --> 1 --> 2.}
#' \item{NC}{vector containing the number of covariates on transitions 0 --> 1,
#' 0 --> 2, 1 --> 2.} \item{responseTrans}{model response for the 0 --> 1
#' transition. \code{Hist} or \code{Surv} object.} \item{responseAbs}{model
#' response for the 0 --> 2 transition. \code{Hist} or \code{Surv} object.}
#' \item{time}{times for which transition intensities have been evaluated for
#' plotting. Vector in the Weibull approach. Matrix in the penalized likelihhod
#' approach for which the colums corresponds to the transitions 0 --> 1, 1 -->
#' 2, 0 --> 2.} \item{RR}{vector of relative risks.}
#' \item{V}{variance-covariance matrix derived from the Hessian of the log-likelihood
#' if using method="Weib" or, from the Hessian of the penalized log-likelihood
#' if using method="Splines".}
#' \item{se}{standart errors of the
#' regression parameters.} \item{Xnames01}{names of covariates on 0 --> 1.}
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
#' \item{CV}{a binary variable equals to 1 when search of the smoothing
#' parameters \link{kappa} by approximated cross-validation, 1 otherwise. The
#' default is 0.} \item{kappa}{vector containing the smoothing parameters for
#' transition 0 --> 1, 0 --> 2, 1 --> 2 used to estimate the model by the
#' penalized likelihood approach.} \item{CVcrit}{cross validation criteria.}
#' \item{DoF}{degrees of freedom of the model.} \item{na.action}{observations
#' deleted if missing values.}
#' @author R: Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr> Fortran:
#' Pierre Joly <Pierre.Joly@@isped.u-bordeaux2.fr>
#' @seealso \code{\link{print.idm}}
#' \code{\link{summary.idm}}
#' \code{\link{predict.idm}}
#' @references D. Marquardt (1963). An algorithm for least-squares estimation
#' of nonlinear parameters.  \emph{SIAM Journal of Applied Mathematics},
#' 431-441.
#' @keywords illness-death
#'
##' @examples
##' library(lava)
##' library(prodlim)
##' set.seed(17)
##' d <- simulateIDM(100)$data
##' # right censored data
##' fitRC <- idm(formula01=Hist(time=observed.illtime,event=seen.ill)~X1+X2,
##'              formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
##'              formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d)
##' fitRC
##'
##' \dontrun{
##' set.seed(17)
##' d <- simulateIDM(300)$data
##' fitRC.splines <- idm(formula01=Hist(time=observed.illtime,event=seen.ill)~X1+X2,
##'              formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
##'              formula12=Hist(time=observed.lifetime,event=seen.exit)~1,data=d,
##'              method="splines")
##' }
##' # interval censored data
##' fitIC <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2,
##'              formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
##'              formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d)
##' fitIC
##'
##' \dontrun{
##'
##'     data(Paq1000)
##'
##'     # Illness-death model with certif on the 3 transitions
##'     # Weibull parametrization and likelihood maximization
##'
##'     fit.weib <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
##'                     formula01=Hist(time=list(l,r),event=dementia)~certif,
##'                     data=Paq1000)
##'
##'     # Illness-death model with certif on transitions 01 and 02
##'     # Splines parametrization and penalized likelihood maximization
##'     fit.splines <-  idm(formula02=Hist(time=t,event=death,entry=e)~certif,
##'                         formula01=Hist(time=list(l,r),event=dementia)~certif,
##'                         formula12=~1,
##'                         method="splines",
##'                         data=Paq1000)
##'     fit.weib
##'     summary(fit.splines)
##' }
##'
#' @importFrom prodlim Hist
#' @useDynLib SmoothHazardoptim9
#' @export
idm <- function(formula01,
                formula02,
                formula12,
                data,

                method="Weib",
                scale.X=T,
                maxiter=100,
                maxiter.pena=10,

                eps=c(5,5,3),
                eps.spline=3,
                eps.eigen=2,
                n.knots=NULL,
                knots="equidistant",
                type.quantile=1,
                subset=NULL,
                na.action = na.fail,

                B=NULL,
                posfix=NULL,

                gauss.point=10,

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
                step.sequential=F,
                option.sequential=list(cutoff=10^-3,
                                    min=20,
                                    step=10),

                nproc=1,
                clustertype="FORK",
                envir=parent.frame()){


    # {{{ check formula
    ################################ check what has been given ##################
    call <- match.call()
    ptm <- proc.time()
    if(missing(formula01))stop("Argument formula01 is missing.")
    if(missing(formula02))stop("Argument formula02 is missing.")
    if(!inherits(formula01,"formula"))stop("The argument formula01 must be a class 'formula'.")
    if(!inherits(formula02,"formula"))stop("The argument formula02 must be a class 'formula'.")
    
    if(!method%in%c("Weib","splines"))stop("The argument method needs to be either splines or Weib")
    ## if(missing(formula02)) formula02 <- formula01
    if(missing(formula12)) formula12 <- formula02
    # }}}
    # {{{ evaluate formula in data
    if(missing(data)) stop("Need a data frame.")
    if(sum(is.na(data))>0)stop("Need a data frame with no missing data.")

    if(!inherits(data,"data.frame"))stop("Argument 'data' must be a data.frame")
 
    ############################### get database defined by formulas ###########
    m <- match.call()
    m01 <- m02 <- m12 <- m[match(c("","data","subset","na.action"),names(m),nomatch=0)]
    m01$formula <- formula01
    m02$formula <- formula02
    m12$formula <- formula12
    m01[[1]] <- m02[[1]] <- m12[[1]] <- as.name("model.frame")

    ## dealing with missing data if no covariate on transition 0->2
    if(anyNA(data)){
      variables=unique(c(all.vars(formula01),all.vars(formula02),all.vars(formula12)))
      data=data[,variables]
      data=na.omit(data)
      m01[[2]] <- m02[[2]] <- m12[[2]] <- data
    }

    m01 <- eval(m01,parent.frame())
    m02 <- eval(m02,parent.frame())
    m12 <- eval(m12,parent.frame())


    responseTrans <- stats::model.response(m01)
    responseAbs <- stats::model.response(m02)
    # }}}
    # {{{ extract covariates
    ## formula01
   
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
    ## formula12
    x12 <- model.matrix(formula12,data=m12)[, -1, drop = FALSE]
    NC12 <- NCOL(x12)


    if (NC12>0)
        Xnames12 <- colnames(x12)
    else
        Xnames12 <- NULL


    # }}}
    # {{{ prepare censored event times
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
    # }}}
    # {{{ check data for integrity
    if (attr(responseAbs,"cens.type")=="intervalCensored") stop("No method available when the transtion to the absorbing state is interval censored.")
    if (isIntervalCensored && any(Rtime<Ltime)) stop("Misspecified transitition times:\nSome left interval limits are greater than the corresponding right limits.")
    # }}}
    # {{{ call Fortran function weib and collect results
    ## ===R provides===
    ## Variable name| Explanation|Dimension|Storage mode|Remark
    ## entrytime| truncation time|length N|double|if is_truncated = 0 then length 0
    ## l| right censored  time or left border of interval censored obs for illness|length N|double|
    ## r| right border of interval censored obs for illness|length N|double|without interval censoring r=l
    ## d| right censored  time or obs time for death|length N|double|
    ## idm| illness status (1= event, 0=no event)|length N|integer|
    ## idd| death status (1= event, 0=no event)|length N|integer|
    ## x01| covariate matrix in vector form for 0---> 1: c(X_11,X_12,...,X_1P,X_21,X22,...,X_P*N)|length N*P|double|dummy variables
    ## x02| covariate matrix in vector form for 0---> 2: c(X_11,X_12,...,X_1P,X_21,X22,...,X_P*N)|length N*P|double|
    ## x12| covariate matrix in vector form for 1---> 2: c(X_11,X_12,...,X_1P,X_21,X22,...,X_P*N)|length N*P|double|
    ## cluster| NOT YET |--|--|
    ## strata| NOT YET |--|--|
    ## N| number of subjects|length 1|integer|
    ## P01| number of covariates for 0---> 1|length 1|integer|
    ## P02| number of covariates for 0---> 2|length 1|integer|
    ## P12| number of covariates for 1---> 2|length 1|integer|
    ## is_truncated|0=no, 1=yes|length 1|integer|
    ## eps|convergence criteria: 1:likelihood,2:parameter est,3:gradient parameter est |length 3|integer|example eps=c(7,4,5) then use 10^-7,10^-4,10^-5. Defaults to c(5,5,3)
    ## maxiter| maximum number of iteration | length 1 | integer | > 0 default to 200


      if(!inherits(option.sequential$cutoff,c("numeric","integer")))stop("The cutoff for spline has to a numeric or integer.")
      if(!inherits(maxiter,c("numeric","integer"))|(maxiter!=floor(maxiter)))stop("Maxiter has to be an integer.")
      if(!inherits(nproc,c("numeric","integer"))|(nproc!=floor(nproc)))stop("nproc has to be an integer.")
    
      # nbr of quadrature points for estimating integral in idm without penalisation
      if(!gauss.point%in%c(10,15,21,31,41,51,61))stop("Argument type.quantile has to a numeric : 10, 15, 21, 31, 51 or 61.")
      if(!inherits(step.sequential,"logical"))stop("Argument step.sequential has to be TRUE or FALSE.")

    #################### if the user want to fix splines parameters put #######
    #################### (step.sequential==T) #################################
      if(step.sequential==T){

        if(!inherits(option.sequential$cutoff,c("numeric","integer")))stop("The cutoff for spline has to be a numeric or integer.")
        if(!(inherits(option.sequential$step,c("numeric","integer")))|(option.sequential$step!=floor(option.sequential$step)))stop("The step has to be a integer.")
        if(!(inherits(option.sequential$min,c("numeric","integer")))|(option.sequential$min!=floor(option.sequential$min)))stop("The min has to be a integer.")
        if(option.sequential$step<=0)stop("Steps need to be at least of 1.")
        if(option.sequential$min<=0)stop("The minimum of iteration must be at least 0")
        if(option.sequential$cutoff<=0)stop("The cutoff for spline need to be positive.")
        if(!(option.sequential$min<=maxiter))stop(paste0("The algorithm must do at least one iteration before reaching ",maxiter," , min inferior or equal to maxiter."))
        if(!(option.sequential$step<=maxiter))stop(paste0("The step must be inferior or equal to ",maxiter))

      }

      #  	cat("------ Program Splines ------ \n")
      ## check knots

    size1 <- NC01 + NC02 + NC12
    

    noVar<-c(ifelse(as.integer(NC01)>0,0,1),
             ifelse(as.integer(NC02)>0,0,1),
             ifelse(as.integer(NC12)>0,0,1))


    nvat01 <- ifelse(noVar[1]==1,0,NC01)
    nvat02 <- ifelse(noVar[2]==1,0,NC02)
    nvat12 <- ifelse(noVar[3]==1,0,NC12)


    #cv criterias
    epsa<-0.1^eps[1]
    epsb<-0.1^eps[2]
    epsd<-0.1^eps[3]
    eps.spline<-0.1^eps.spline[1]
    eps.eigen<-0.1^eps.eigen[1]
    
    troncature<-ifelse(truncated==T,1,0)
    
 ####################### define profile of subjects ############################
    
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
    
    #################### defines knots placements for splines ##################
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
    
    
###### check if we have penalty and if some parameters are fixed by user #######
    fix0<-rep(0,size_V)
    if(is.null(penalty)){
      penalty<-"none"}
    if(!penalty%in%c("none","lasso","ridge","elasticnet","mcp","scad")){
      stop(paste0("Parameter penalty must be either : lasso, ridge, elastic.net, mcp or scad"))}
    
    
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
    

    dimnva01<-ifelse(nvat01==0,1,nvat01)
    dimnva02<-ifelse(nvat02==0,1,nvat02)
    dimnva12<-ifelse(nvat12==0,1,nvat12)

    fit <- NULL

    ni<-0
    NC<-c(NC01,NC02,NC12)
    
    if(noVar[1]==1){ve01<-as.double(rep(0,N))}else{ve01<-as.double(x01)}
    if(noVar[2]==1){ve02<-as.double(rep(0,N))}else{ve02<-as.double(x02)}
    if(noVar[3]==1){ve12<-as.double(rep(0,N))}else{ve12<-as.double(x12)}
    
######################### algorithm if no penalty ##############################
    if(penalty=="none"){
 
######################### call for splines #####################################
      if(method=="splines"){
        
        if(sum(fix0)>0){
          bfix<-b[fix0==1]
          b<-b[fix0==0]
        }else{bfix<-1}
        
      
        out<-idm.no.penalty(b,clustertype,epsa,epsb,epsd,nproc,maxiter,size_V,size_spline,noVar,bfix,
                            fix0,knots01,knots02,knots12,ctime,N,nknots01,nknots02,nknots12,
                            ve01,ve02,ve12,dimnva01,dimnva02,dimnva12,nvat01,nvat02,nvat12,
                            t0,t1,t2,t3,troncature,gauss.point,step.sequential,option.sequential)
  
      
        fix<-out$fix0
        
        beta<-rep(NA,size_V)
        beta[fix==0]<-out$b
        beta[fix==1]<-out$bfix
        npm<-sum(fix==0)
### if CV is true keep V and solve V to have inverse for confidence intervals###
        if(out$istop==1){
          
          
          Vr <- matrix(0,npm,npm)
          Vr[upper.tri(Vr,diag=TRUE)] <- out$v[1:(npm*(npm+1)/2)]
          Vr <- t(Vr)
          Vr[upper.tri(Vr,diag=TRUE)] <- out$v[1:(npm*(npm+1)/2)]
          V <- matrix(0,size_V,size_V)
          
          posfix<-which(fix==1)
          V[setdiff(1:size_V,posfix),setdiff(1:size_V,posfix)] <- Vr
          
          Hr<--solve(Vr)
          H<- matrix(0,size_V,size_V)
          H[setdiff(1:size_V,posfix),setdiff(1:size_V,posfix)] <- Hr
          
        }else{ #otherwise put 0
          Hr <- matrix(0,npm,npm)
          Hr[upper.tri(Hr,diag=TRUE)] <- out$v[1:(npm*(npm+1)/2)]
          Hr <- t(Hr)
          Hr[upper.tri(Hr,diag=TRUE)] <- out$v[1:(npm*(npm+1)/2)]
          H <- matrix(0,size_V,size_V)
          
          posfix<-which(fix==1)
          H[setdiff(1:size_V,posfix),setdiff(1:size_V,posfix)] <- Hr
          
          V<- matrix(0,size_V,size_V)
          
          
        }
        
        # output of beta and HR 
        if (sum(NC)>0){  
          
          # if at least one covariate
          
          betaCoef <- beta[(size_spline+1):size_V]
          names(betaCoef) <- c(Xnames01,Xnames02,Xnames12)
          fit$coef <- betaCoef
          fit$HR <- exp(betaCoef)
          
        }
        
        
        # splines parameters values
        theta_names <- cbind(c(rep("theta01",(nknots01+2)),rep("theta02",(nknots02+2)),rep("theta12",(nknots12+2))),c((1:(nknots01+2)),(1:(nknots02+2)),(1:(nknots12+2))))
        theta_names <- as.vector(apply(theta_names,1,paste,collapse=" "))
        
        
        names(beta)<- c(theta_names,c(Xnames01,Xnames02,Xnames12))
        names(fix)<-c(theta_names,c(Xnames01,Xnames02,Xnames12))
        colnames(V) <- c(theta_names,c(Xnames01,Xnames02,Xnames12))
        rownames(V) <- c(theta_names,c(Xnames01,Xnames02,Xnames12))
        colnames(H) <- c(theta_names,c(Xnames01,Xnames02,Xnames12))
        rownames(H) <- c(theta_names,c(Xnames01,Xnames02,Xnames12))
        
        theta01<-beta[1:(nknots01+2)]
        theta02<-beta[(nknots01+3):(nknots01+nknots02+4)]
        theta12<-beta[(nknots01+nknots02+5):(nknots01+nknots12+nknots02+6)]
        
        # No value given by BIC as no penalty
        fit$BIC<-NULL
        
        cv<-list(ca=out$ca,cb=out$cb,rdm=out$rdm)
        
        fit$knots01 <- knots01
        fit$knots02 <- knots02
        fit$knots12 <- knots12
        fit$nknots01 <- nknots01
        fit$nknots02 <- nknots02
        fit$nknots12 <- nknots12
        fit$theta01 <- theta01
        fit$theta02 <- theta02
        fit$theta12 <- theta12
        
        
############# output times to do prediction on #################################
        fit$time <- matrix(NA,ncol=3,nrow=100)
        fit$time[,1]<-seq(from=knots01[1],to=knots01[length(knots01)],length.out=100)
        fit$time[,2]<-seq(from=knots02[1],to=knots02[length(knots02)],length.out=100)
        fit$time[,3]<-seq(from=knots12[1],to=knots12[length(knots12)],length.out=100)
        
        
      }
      ######################### call for weibull #####################################
      
      if(method=="Weib"){
        
        out <- idm.weib(b=b,
                                fix0=fix0,
                                size_V=size_V,
                                clustertype=clustertype,
                                epsa=epsa,
                                epsb=epsb,
                                epsd=epsd,
                                eps.eigen=eps.eigen,
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
                                troncature=troncature)
        
        
        
        fix<-out$fix0
        beta<-rep(NA,size_V)
        beta[fix==0]<-out$b
        beta[fix==1]<-out$bfix
        npm<-sum(fix==0)
        
### if CV is true keep V and solve V to have inverse for confidence intervals###
        
        if(out$istop==1){
          
          
          Vr <- matrix(0,npm,npm)
          Vr[upper.tri(Vr,diag=TRUE)] <- out$v[1:(npm*(npm+1)/2)]
          Vr <- t(Vr)
          Vr[upper.tri(Vr,diag=TRUE)] <- out$v[1:(npm*(npm+1)/2)]
          V <- matrix(0,size_V,size_V)
          
          posfix<-which(fix==1)
          V[setdiff(1:size_V,posfix),setdiff(1:size_V,posfix)] <- Vr
          
          Hr<--solve(Vr)
          H<- matrix(0,size_V,size_V)
          H[setdiff(1:size_V,posfix),setdiff(1:size_V,posfix)] <- Hr
          
        }else{
          Hr <- matrix(0,npm,npm)
          Hr[upper.tri(Hr,diag=TRUE)] <- out$v[1:(npm*(npm+1)/2)]
          Hr <- t(Hr)
          Hr[upper.tri(Hr,diag=TRUE)] <- out$v[1:(npm*(npm+1)/2)]
          H <- matrix(0,size_V,size_V)
          
          posfix<-which(fix==1)
          H[setdiff(1:size_V,posfix),setdiff(1:size_V,posfix)] <- Hr
          
          V<- matrix(0,size_V,size_V)
          
          
        }
        # output of beta and HR 
        if (sum(NC)>0){  
          
          # if at least one covariate
          
          betaCoef <- beta[(6+1):size_V]
          names(betaCoef) <- c(Xnames01,Xnames02,Xnames12)
          fit$coef <- betaCoef
          fit$HR <- exp(betaCoef)
          
        }
        
        # weibull parameters
        
        theta_names <- c("modelPar1 01",
                         "modelPar2 01",
                         "modelPar1 02",
                         "modelPar2 02",
                         "modelPar1 12",
                         "modelPar2 12")
        
        names(beta)<- c(theta_names,c(Xnames01,Xnames02,Xnames12))
        names(fix)<- c(theta_names,c(Xnames01,Xnames02,Xnames12))
        colnames(V) <- c(theta_names,c(Xnames01,Xnames02,Xnames12))
        rownames(V) <- c(theta_names,c(Xnames01,Xnames02,Xnames12))
        colnames(H) <- c(theta_names,c(Xnames01,Xnames02,Xnames12))
        rownames(H) <- c(theta_names,c(Xnames01,Xnames02,Xnames12))
        
        modelPar<-beta[1:6]
        
        # no value for BIC as no penalty
        fit$BIC<-NULL
        
        cv<-list(ca=out$ca,cb=out$cb,rdm=out$rdm)
        
        
        fit$modelPar <- modelPar
        
        
      }
      
      
      fit$loglik <- c(out$fn.value,NULL)
      lambda<-NULL
      alpha<-NULL
      
    }else{

############################# Penalty algorithm ################################

          if(nvat01==0 & nvat02==0 & nvat12==0)stop("To perform penalisation you need explanatory variables in each transition")
          
          # permits to not penalise on some parameters
          if(is.null(penalty.factor)){
            penalty.factor<-rep(1, nvat01+nvat02+nvat12)
          }else{
            if(min(penalty.factor)<0 | max(penalty.factor)>1 | round(penalty.factor)!=penalty.factor | length(penalty.factor)!=(nvat01+nvat02+nvat12)){
              stop(paste0("Penalty.factor need to be a vector of 0 and 1 of length : ",nvat01+nvat02+nvat12))
            }
          }


          if(scale.X==T){
            # to know which variable to center and reduces : 
            names<-as.vector(unlist(lapply(m01, class)))[-1]
            
            # standardize variable : 
            
            # we need to center and reduce variables when penalty is done 
            namesx01<-as.vector(unlist(lapply(m01, class)))[-1]
            scaledx01<-scale(x01[,which(namesx01=="numeric")])
            # keep center and reduces parameters 
            xm01<-attr(scaledx01,"scaled:center")
            xs01<-attr(scaledx01,"scaled:scale")
            x01[,which(namesx01=="numeric")]<-scaledx01
            
            
            
            # to know wich variable to center and reduces : 
            names<-as.vector(unlist(lapply(m02, class)))[-1]
            
            # standardize variable : 
            
            # we need to center and reduce variables when penalty is done 
            namesx02<-as.vector(unlist(lapply(m02, class)))[-1]
            scaledx02<-scale(x02[,which(namesx02=="numeric")])
            # keep center and reduces parameters 
            xm02<-attr(scaledx02,"scaled:center")
            xs02<-attr(scaledx02,"scaled:scale")
            x02[,which(namesx02=="numeric")]<-scaledx02
            
            
            # to know wich variable to center and reduces : 
            # no need [-1] as no element before ~
            names<-as.vector(unlist(lapply(m12, class)))
            
            # standardize variable : 
            
            # we need to center and reduce variables when penalty is done 
            namesx12<-as.vector(unlist(lapply(m12, class)))
            scaledx12<-scale(x12[,which(namesx12=="numeric")])
            # keep center and reduces parameters 
            xm12<-attr(scaledx12,"scaled:center")
            xs12<-attr(scaledx12,"scaled:scale")
            x12[,which(namesx12=="numeric")]<-scaledx12
            
            if(noVar[1]==1){ve01<-as.double(rep(0,N))}else{ve01<-as.double(x01)}
            if(noVar[2]==1){ve02<-as.double(rep(0,N))}else{ve02<-as.double(x02)}
            if(noVar[3]==1){ve12<-as.double(rep(0,N))}else{ve12<-as.double(x12)}
            
            
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
              output<-deriva_gradient(b=c(b[1:size_spline],rep(0,size_V-size_spline)),
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
                                  gausspoint=gauss.point)
              
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
            
            out<-idm.penalty(b=b,
                             fix0=fix0,
                             size_V=size_V,
                             size_spline=size_spline,
                             clustertype=clustertype,
                             epsa=epsa,
                             epsb=epsb,
                             epsd=epsd,
                             eps.spline=eps.spline,
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
                             gauss.point=gauss.point,
                             nlambda01=nlambda01,
                             lambda01=lambda01,
                             nlambda02=nlambda02,
                             lambda02=lambda02,
                             nlambda12=nlambda12,
                             lambda12=lambda12,
                             alpha=alpha,
                             penalty.factor=penalty.factor,
                             step.sequential=step.sequential,
                             option.sequential=option.sequential,
                             penalty=penalty)
            
            
            
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
            
########################## define BIC and GCV ##################################
            if(dim(beta)[2]>1){
              fit$BIC<--2*out$fn.value+log(N)*colSums(beta[(size_spline+1):size_V,]!=0)
            }else{fit$BIC<--2*out$fn.value+log(N)*sum(beta[(size_spline+1):size_V,]!=0)}
            
            
            fit$GCV<-rep(NA,dim(lambda)[2])
            
            
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
            fit$theta01 <- theta01
            fit$theta02 <- theta02
            fit$theta12 <- theta12

            fit$time <- matrix(NA,ncol=3,nrow=100)
            fit$time[,1]<-seq(from=knots01[1],to=knots01[length(knots01)],length.out=100)
            fit$time[,2]<-seq(from=knots02[1],to=knots02[length(knots02)],length.out=100)
            fit$time[,3]<-seq(from=knots12[1],to=knots12[length(knots12)],length.out=100)
            
          }
          if(method=="Weib"){
              #	cat("------ Program Weibull ------ \n")
            # some initial steps to have values for weibull parameters
            output.mla<- marqLevAlg::mla(b=b[1:6],
                             fn=idmlLikelihoodweib,
                             epsa=epsa,
                             epsb=epsb,
                             epsd=epsd,
                             nproc=nproc,
                             clustertype = clustertype,
                             maxiter=maxiter,
                             minimize=F,
                             npm=6,
                             npar=size_V,
                             bfix=b[7:(size_V-6)],
                             fix=c(rep(0,6),rep(1,size_V-6)),
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
                             troncature=troncature)
            
            # take thoses values if converged only otherwise thoses
            # by default or by the user
            if(output.mla$istop==1){
              b<-c(output.mla$b,b[7:size_V])}
                   
            
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
              output<-deriva_gradient(b=c(b[1:6],rep(0,size_V-6)),
                                  nproc=nproc,
                                  funcpa=idmlLikelihoodweib,
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
                                  troncature=troncature)
              
              if(nproc>1){parallel::stopCluster(clustpar)}
              
              
              ## what should we do if max(output$v) == 0
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
           
              out <- idm.penalty.weib(b=b,
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
                               penalty=penalty)
              
              
              
              
              beta<-out$b
              fix<-fix0
              lambda<-out$lambda
              alpha<-out$alpha
              
              H<-out$H
              
              beta<-as.matrix(beta)
              lambda<-as.matrix(lambda)
              
              theta_names<-c("modelPar1 01",
              "modelPar2 01",
              "modelPar1 02",
              "modelPar2 02",
              "modelPar1 12",
              "modelPar2 12")
              rownames(beta) <-c(theta_names,Xnames01,Xnames02,Xnames12)
              names(fix)<-c(theta_names,c(Xnames01,Xnames02,Xnames12))
              rownames(lambda)<-c("lambda01","lambda02","lambda12")
              
              modelPar<-beta[1:6,]
              betaCoef <- beta[7:size_V,]
              betaCoef<-as.matrix(betaCoef)
              fit$coef <- betaCoef
              fit$HR <- exp(betaCoef)
              
              
####################   calculate BIC    #######################################
              if(dim(beta)[2]>1){
                fit$BIC<--2*out$fn.value+log(N)*colSums(beta[7:size_V,]!=0)
              }else{fit$BIC<--2*out$fn.value+log(N)*sum(beta[7:size_V,]!=0)}
              
              
###################### calculate GCV ###########################################
              fit$GCV<-rep(NA,dim(lambda)[2])
              
              npm<-sum(fix0[7:size_V]==0)
              npm01<-ifelse(nvat01>0,sum(fix0[7:(6+nvat01)]==0),0)
              npm02<-ifelse(nvat02>0,sum(fix0[(7+nvat01):(6+nvat01+nvat02)]==0),0)
              npm12<-ifelse(nvat12>0,sum(fix0[(7+nvat01+nvat02):size_V]==0),0)
              
              
              V<-matrix(0,ncol=dim(lambda)[2]*npm,nrow=npm)

              
              
              for(i in 1:dim(lambda)[2]){
                id.keep<-which(betaCoef[fix0[7:size_V]==0,i]!=0)
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
                      #V[,(1+(i-1)*npm):(i*npm)]<-solve(out$H[,(1+(i-1)*npm):(i*npm),drop=FALSE])
                      id.keepV<-(1+(i-1)*npm):(i*npm)
                      id.keepV<-id.keepV[id.keepV%in%((i-1)*npm+id.keep)]
                      
                      V[id.keep,id.keepV]<-solve(H_spec)
                      # maximisation issue thus : 10/04/24
                      # as mla in maximisation return -second derivatives 
                      H_pl<-H_spec+lambda.matrix
                      trace_model<-lava::tr(solve(H_pl)%*%H_spec)
                      fit$GCV[i]<--1/N*(out$fn.value[i]-trace_model)}
                    }
                  }
                  
                }
              }
                
              cv<-list(ca.beta=out$ca.beta,ca.spline=out$ca.spline,ca.validity=out$ca.validity,
                       cb=out$cb,rdm=NULL)
              
             
              fit$modelPar<-modelPar
          }
          
          fit$loglik <- c(out$fn.value,out$fn.value.pena)
          
          



        }


        fit$call <- call
        fit$terms <- list("Formula01"=terms(formula01),
                          "Formula02"=terms(formula02),
                          "Formula12"=terms(formula12))
        
        fit$cv <- cv
        fit$niter <- out$ni

        fit$converged <- out$istop

        fit$N <- N
        fit$events1 <- sum(idm)
        fit$events2 <- sum(idd)
        fit$NC <- NC
        fit$responseAbs <- responseAbs
        fit$responseTrans <- responseTrans


        fit$V <- V
        fit$H <- H
        fit$fix<-fix
        fit$lambda<-lambda
        fit$alpha<-alpha
        fit$penalty<-penalty


        if(NC01>0) fit$Xnames01 <- Xnames01
        if(NC02>0) fit$Xnames02 <- Xnames02
        if(NC12>0) fit$Xnames12 <- Xnames12

       ################### need to have levels for values binary/qualitative ###
        ################## useful for prediction purpose #######################
        
        var.level<-var.class<-list(NULL)
        N_var<-c(all.vars(formula01)[all.vars(formula01)%in%labels(terms(formula01))],
          all.vars(formula02)[all.vars(formula02)%in%labels(terms(formula02))],
          all.vars(formula12)[all.vars(formula12)%in%labels(terms(formula12))])
        length(var.level)<-length(var.class)<-length(N_var)
        m<-1
        if(NC01>0){
        for(k in 2:dim(m01)[2]){
          var.class[[m]]<-class(m01[,k])
          if(class(m01[,k])[1]%in%c("character","factor")){
          var.level[[m]]<-names(table(m01[,k]))}
          m<-m+1
        }
        }
        if(NC02>0){
          for(k in 2:dim(m02)[2]){
            var.class[[m]]<-class(m02[,k])
            if(class(m02[,k])[1]%in%c("character","factor")){
              var.level[[m]]<-names(table(m02[,k]))}
            m<-m+1
          }
        }
        
        if(NC12>0){
          for(k in 1:dim(m12)[2]){
            var.class[[m]]<-class(m12[,k])
            if(class(m12[,k])[1]%in%c("character","factor")){
              var.level[[m]]<-names(table(m12[,k]))}
            m<-m+1
          }
        }
        
        fit$levels<-list(values=var.level,
                         class=var.class)
        

        fit$na.action <- "na.fail"


        fit$maxtime <- amax
        fit$mintime <- amin
        class(fit) <- "idm"
        fit$runtime <- proc.time()-ptm
        return(fit)
}
