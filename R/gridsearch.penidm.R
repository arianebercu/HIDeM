
###' Code : 
#' @param gridmethod On which indicator the grid should be based, either BIC or GCV
#' @param sizegrid The number of values that we want for each lambda, its a vector 
#' of three elements, \code{0 --> 1}, \code{0 --> 2} and \code{1 --> 2}. The size of each
#' element should not exceed respectively nlambda01, nlambda02 and nlambda12
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
#' @param na.action how NAs are treated. The default is first, any
#' na.action attribute of data, second a na.action setting of options,
#' and third 'na.fail' if that is unset. The 'factory-fresh' default
#' is na.omit. Another possible value is NULL.
#' @param scale.X do you want to center and reduce your explanatory variables
#' @param posfix index of fixed parameters 
#' @param gausspoint gauss quadrature points in the approximation of integrals
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
#' \item{m01}{Model estimated on 0 --> 1} \item{m02}{Model estimated on 0 --> 2} 
#' \item{m12}{Model estimated on 1 --> 2} \item{lambda01}{vector of lambda penalty parameters
#' on transition 0 --> 1 minimising the BIC or GCV in model m01} \item{lambda02}{vector of lambda penalty parameters
#' on transition 0 --> 2 minimising the BIC or GCV in model m02} \item{lambda12}{vector of lambda penalty parameters
#' on transition 1 --> 2 minimising the BIC or GCV in model m12}
#' \item{alpha}{the penalty parameter alpha} 
#' \item{sizegrid}{the size of lambda penalty parameters for each transition 0 -->1, 0 -->2 and 1 -->2}
#' \item{gridmethod}{On which indicator the grid should be based, either BIC or GCV}
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
#' @seealso \code{\link{print.idm}}
#' @examples
#' 
#' \dontrun{
#' library(lava)
#' library(prodlim)
#' set.seed(17)
#' d <- simulateIDM(n=1000,beta01=c(1,1,0,0.5,0.5,rep(0,5)),
#' beta02=c(1,0,0,0,0.5,rep(0,5)),beta12=c(1,0,0,0,0.5,rep(0,5)))$data
#' fitgrid<- gridsearch.penidm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
#' formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
#' formula12=~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10, nlambda01=20, nlambda02=20, nlambda12=20,
#' data=d,penalty="lasso")
#' }
#' \code{\link{summary.idm}}
#' \code{\link{predict.idm}}
#' @references D. Marquardt (1963). An algorithm for least-squares estimation
#' of nonlinear parameters.  \emph{SIAM Journal of Applied Mathematics},
#' 431-441.
#' @importFrom prodlim Hist
#' @useDynLib SmoothHazardoptim9
#' @export
gridsearch.penidm <- function(
                gridmethod="BIC",
                sizegrid=c(3,3,3),
                formula01,
                formula02,
                formula12,
                data,

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
                penalty="lasso",
                penalty.factor=NULL,
                alpha=ifelse(penalty=="scad",3.7,
                             ifelse(penalty=="mcp",3,
                                    ifelse(penalty%in%c("elasticnet","corrected.elasticnet"),0.5,1))),
                nproc=1,
                clustertype="FORK"){
  
  call <- match.call()

  if(!inherits(data,"data.frame"))stop("Argument 'data' must be a data.frame")
  
  if(!penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet","mcp","scad")){
    stop(paste0("Parameter penalty must be either : lasso, ridge, elastic.net, corrected.elasticnet, mcp or scad"))}
  # check formula 
  
  if(missing(formula01))stop("Argument formula01 is missing.")
  if(missing(formula02))stop("Argument formula02 is missing.")
  if(!inherits(formula01,"formula"))stop("The argument formula01 must be a class 'formula'.")
  if(!inherits(formula02,"formula"))stop("The argument formula02 must be a class 'formula'.")
  if(missing(formula12)) formula12 <- formula02
  
  if(length(sizegrid)!=3|(!class(sizegrid)%in%c("numeric","integer"))){
    stop("Sizegrid should be a numeric vector of 3 elements")
  }
  # do check on gridmethod
  
  if(!gridmethod%in%c("BIC","GCV")){
    stop("gridmethod should be BIC or GCV")
  }
  # do check on sizegrid 
  
  if(sizegrid[1]>nlambda01){sizegrid[1]<-nlambda01}
  if(sizegrid[2]>nlambda02){sizegrid[1]<-nlambda02}
  if(sizegrid[3]>nlambda12){sizegrid[3]<-nlambda12}
  
  if(sizegrid[1]<=1){sizegrid[1]<-1}
  if(sizegrid[2]<=1){sizegrid[1]<-1}
  if(sizegrid[3]<=1){sizegrid[3]<-1}
  

  
    if(method=="Weib"){
    posfix01<-c(posfix,3:6)
    posfix02<-c(posfix,1,2,5,6)
    posfix12<-c(posfix,1:4)
    }else{
      #define the nknots
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
      
      start<-nknots01+3
      end<-nknots01+nknots02+nknots12+6
      # fix 02,12 the knots at
      posfix01<-c(posfix,start:end)
      
      start<-nknots01+nknots02+5
      # fix 01, 12 the knots at
      posfix02<-c(posfix,1:(nknots01+2),start:end)
      
      end<-nknots01+nknots02+4
      # fix 01, 02 the knots at
      posfix12<-c(posfix,1:end)
    }
  
  ############################### 01 #####################################
 # be careful to formulas environnement 
  
  f12<-as.formula("~1")
  f02<-as.character(formula02)[2]
  f02<-as.formula(paste0(f02,"~1"))
  f01<-formula01
  environment(f01)<-environment(f02)<-environment(f12)
  
  
  m01 <-SmoothHazardoptim9::idm(formula02 = f02,
                                  formula01 = f01,
                                  formula12 = f12,
                                  data=data,
                                  nproc=nproc, 
                                  penalty=penalty,
                                  eps=eps,
                                  method=method,
                                  scale.X=scale.X,
                                  nlambda01 = nlambda01,
                                  lambda02=0.0001,
                                  lambda12=0.0001,
                                  posfix=posfix01,
                                  alpha=alpha,
                                  maxiter=maxiter,
                                  maxiter.pena = maxiter.pena,
                                  n.knots=n.knots,
                                  knots=knots,
                                  type.quantile=type.quantile,
                                  na.action =na.action,
                                  B=B,
                                  gausspoint=gausspoint,
                                  clustertype=clustertype)
  
  
  ################################ 02 ############################################
  
    
    f12<-as.formula("~1")
    f01<-as.character(formula01)[2]
    f01<-as.formula(paste0(f01,"~1"))
    f02<-formula02
    
    environment(f02)<-environment(f01)<-environment(f12)
    
    
    m02 <-SmoothHazardoptim9::idm(formula02 = f02,
                                  formula01 = f01,
                                  formula12 = f12,
                                  data=data,
                                  nproc=nproc, 
                                  penalty=penalty,
                                  eps=eps,
                                  method=method,
                                  scale.X=scale.X,
                                  lambda02=lambda02,
                                  nlambda02 = nlambda02,
                                  lambda01=0.0001,
                                  lambda12=0.0001,
                                  posfix=posfix02,
                                  alpha=alpha,
                                  maxiter=maxiter,
                                  maxiter.pena = maxiter.pena,
                                  n.knots=n.knots,
                                  knots=knots,
                                  type.quantile=type.quantile,
                                  na.action =na.action,
                                  B=B,
                                  gausspoint=gausspoint,
                                  clustertype=clustertype)
 
  
  ######################################## 12 ###################################
  
  
  f01<-as.character(formula01)[2]
  f01<-as.formula(paste0(f01,"~1"))
  f02<-as.character(formula02)[2]
  f02<-as.formula(paste0(f02,"~1"))
  f12<-formula12
  
  environment(f12)<-environment(f02)<-environment(f01)
  
  
  m12 <-SmoothHazardoptim9::idm(formula02 = f02,
                                  formula01 = f01,
                                  formula12 = f12,
                                  data=data,
                                  nproc=nproc, 
                                  penalty=penalty,
                                  eps=eps,
                                  method=method,
                                  scale.X=scale.X,
                                  lambda12=lambda12,
                                  nlambda12 = nlambda12,
                                  lambda01=0.0001,
                                  lambda02=0.0001,
                                  posfix=posfix12,
                                  alpha=alpha,
                                  maxiter=maxiter,
                                  maxiter.pena = maxiter.pena,
                                  n.knots=n.knots,
                                  knots=knots,
                                  type.quantile=type.quantile,
                                  na.action =na.action,
                                  B=B,
                                  gausspoint=gausspoint,
                                  clustertype=clustertype)
  
  
  if(gridmethod=="BIC"){
    
  # sort BIC 
   BIC12<-sort(m12$BIC[m12$converged==1])
   BIC02<-sort(m02$BIC[m02$converged==1])
   BIC01<-sort(m01$BIC[m01$converged==1])
   BIC12<-na.omit(unique(BIC12))
   BIC02<-na.omit(unique(BIC02))
   BIC01<-na.omit(unique(BIC01))
   
   if(length(BIC12)<sizegrid[3]){sizegrid[3]<-length(BIC12)}
   if(length(BIC02)<sizegrid[2]){sizegrid[2]<-length(BIC02)}
   if(length(BIC01)<sizegrid[1]){sizegrid[1]<-length(BIC01)}
   
   # keep only best values in size we want 
   
   BIC12<-BIC12[1:sizegrid[3]]
   BIC02<-BIC02[1:sizegrid[2]]
   BIC01<-BIC01[1:sizegrid[1]]
   
  # identify the lambda associated
  
  id12<-which(m12$BIC%in%BIC12)
  if(length(id12)==0){
    warning("No model converged for transition 1-->2, values without convergence will be given.")
    BIC12<-sort(m12$BIC)
    BIC12<-na.omit(unique(BIC12))
    sizegrid[3]<-ifelse(length(BIC12)<sizegrid[3],length(BIC12),sizegrid[3])
    BIC12<-BIC12[1:sizegrid[3]]
  }
  id12<-which(m12$BIC%in%BIC12)
  lambda12<-m12$lambda[3,id12]
  lambda12<-na.omit(lambda12)
  # verify that good size, if two BIC equal need to keep only one 
  lambda12<-lambda12[1:sizegrid[3]]
  
  id02<-which(m02$BIC%in%BIC02)
  if(length(id02)==0){
    warning("No model converged for transition 0-->2, values without convergence will be give.")
    BIC02<-sort(m02$BIC)[1:sizegrid[2]]
    BIC02<-na.omit(unique(BIC02))
    sizegrid[2]<-ifelse(length(BIC02)<sizegrid[2],length(BIC02),sizegrid[2])
    BIC02<-BIC02[1:sizegrid[2]]}
  id02<-which(m02$BIC%in%BIC02)
  lambda02<-m02$lambda[2,id02]
  lambda02<-na.omit(lambda02)
  # verify that good size, if two BIC equal need to keep only one 
  lambda02<-lambda02[1:sizegrid[2]]
  
  
  id01<-which(m01$BIC%in%BIC01)
  if(length(id01)==0){
    warning("No model converged for transition 0-->1, values without convergence will be give.")
    BIC01<-sort(m01$BIC)[1:sizegrid[1]]
    BIC01<-na.omit(unique(BIC01))
    sizegrid[1]<-ifelse(length(BIC01)<sizegrid[1],length(BIC01),sizegrid[1])
    BIC01<-BIC01[1:sizegrid[1]]}
  id01<-which(m01$BIC%in%BIC01)
  lambda01<-m01$lambda[1,id01]
  lambda01<-na.omit(lambda01)
  # verify that good size, if two BIC equal need to keep only one 
  lambda01<-lambda01[1:sizegrid[1]]
  
  }else{
    # Do on GCV
    
    # sort GCV 
    GCV12<-sort(m12$GCV[m12$converged==1])
    GCV02<-sort(m02$GCV[m02$converged==1])
    GCV01<-sort(m01$GCV[m01$converged==1])
    GCV12<-na.omit(unique(GCV12))
    GCV02<-na.omit(unique(GCV02))
    GCV01<-na.omit(unique(GCV01))
    
    if(length(GCV12)<sizegrid[3]){sizegrid[3]<-length(GCV12)}
    if(length(GCV02)<sizegrid[2]){sizegrid[2]<-length(GCV02)}
    if(length(GCV01)<sizegrid[1]){sizegrid[1]<-length(GCV01)}
    
    # keep only best values in size we want 
    
    GCV12<-GCV12[1:sizegrid[3]]
    GCV02<-GCV02[1:sizegrid[2]]
    GCV01<-GCV01[1:sizegrid[1]]
    
    # identify the lambda associated
    
    id12<-which(m12$GCV%in%GCV12)
    if(length(id12)==0){
      warning("No model converged for transition 1-->2, no values of GCV could be computed.")
      lambda12<-NULL
    }else{
      lambda12<-m12$lambda[3,id12]
      lambda12<-na.omit(lambda12)
      # verify that good size, if two GCV equal need to keep only one 
      lambda12<-lambda12[1:sizegrid[3]]
    }
    
    id02<-which(m02$GCV%in%GCV02)
    if(length(id02)==0){
      warning("No model converged for transition 0-->2, no values of GCV could be computed.")
      lambda02<-NULL
      }else{
        lambda02<-m02$lambda[2,id02]
        lambda02<-na.omit(lambda02)
        # verify that good size, if two GCV equal need to keep only one 
        lambda02<-lambda02[1:sizegrid[2]]
      }
    
    
    id01<-which(m01$GCV%in%GCV01)
    if(length(id01)==0){
      warning("No model converged for transition 0-->1, no values of GCV could be computed.")
      lambda01<-NULL
      }else{
          id01<-which(m01$GCV%in%GCV01)
          lambda01<-m01$lambda[1,id01]
          lambda01<-na.omit(lambda01)
          # verify that good size, if two GCV equal need to keep only one 
          lambda01<-lambda01[1:sizegrid[1]]
      }
  }
    
  fit<-NULL
  fit$m01<-m01
  fit$m02<-m02
  fit$m12<-m12
  fit$lambda12<-lambda12
  fit$lambda01<-lambda01
  fit$lambda02<-lambda02
  fit$alpha<-alpha
  fit$sizegrid<-sizegrid
  fit$gridmethod<-gridmethod
  return(fit)
    
}
