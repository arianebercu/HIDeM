#' Fit an illness-death model using sampling or boostrap
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
#' @param K number of repetitions of sampling or boostrap, by default values 100 
#' @param tau Percentage of observation to use for resampling can go u
#' @param seed value of the seed to initialise the random number generator and ensure reproducibility of the results
#' resampling approach. Possible values are: "subsampling" for sampling without replacement of a proportion tau of the observations, or "bootstrap" for sampling with replacement generating a resampled dataset with as many observations as in the full sample. Alternatively, this argument can be a function to use for resampling. This function must use arguments named data and tau and return the IDs of observations to be included in the resampled dataset.
#' @param calscore True if we want to calculate the calibration score otherwise BIC and GCV only will be given
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
#' 2, 0 --> 2.} \item{intensity01}{matched values of the intensities for
#' transition 0 --> 1.} \item{lowerIntensity01}{lower confidence intervals for
#' the values of the intensities for transition 0 --> 1.}
#' \item{upperIntensity01}{upper confidence intervals for the values of the
#' intensities for transition 0 --> 1.} \item{intensity02}{matched values of
#' the intensities for transition 0 --> 2.} \item{lowerIntensity02}{lower
#' confidence intervals for the values of the intensities for transition 0 -->
#' 2.} \item{upperIntensity02}{upper confidence intervals for the values of the
#' intensities for transition 0 --> 2.} \item{intensity12}{matched values of
#' the intensities for transition 1 --> 2.} \item{lowerIntensity12}{lower
#' confidence intervals for the values of the intensities for transition 1 -->
#' 2.} \item{upperIntensity12}{upper confidence intervals for the values of the
#' intensities for transition 1 --> 2.} \item{RR}{vector of relative risks.}
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
##'              formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d,
##'              conf.int=FALSE)
##' fitRC
##'
##' \dontrun{
##' set.seed(17)
##' d <- simulateIDM(300)$data
##' fitRC.splines <- idm(formula01=Hist(time=observed.illtime,event=seen.ill)~X1+X2,
##'              formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
##'              formula12=Hist(time=observed.lifetime,event=seen.exit)~1,data=d,
##'              conf.int=FALSE,method="splines")
##' }
##' # interval censored data
##' fitIC <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2,
##'              formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
##'              formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d,
##'              conf.int=FALSE)
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
calibrate.penidm <- function(
                K=100, # number of sample,
                tau = 0.5, #% in each fold
                seed = 1, # seed
                calscore=T,
                resampling = "subsampling",
                formula01,
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
                penalty="lasso",
                penalty.factor=NULL,
                alpha=ifelse(penalty=="scad",3.7,
                             ifelse(penalty=="mcp",3,
                                    ifelse(penalty%in%c("elasticnet","corrected.elasticnet"),0.5,1))),
                nproc=1,
                clustertype="FORK"){
  
  if((!class(K)%in%c("integer","numeric"))|K<=1){stop("K need to be an integer higher than 1")}
  K<-round(K) # in case if not integer
  if(missing(data)) stop("Need a data frame.")
  if(sum(is.na(data))>0)stop("Need a data frame with no missing data.")
  
  if(!inherits(data,"data.frame"))stop("Argument 'data' must be a data.frame")
  
  if(!penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet","mcp","scad")){
    stop(paste0("Parameter penalty must be either : lasso, ridge, elastic.net, corrected.elasticnet, mcp or scad"))}
  # resampling check value 
  
  if(!resampling%in%c("subsampling","bootstrap")){
    stop("resampling need to be either bootstrap or subsampling")
  }
  # if we do resampling and calculate calibration we need to have same 
  # lambda for each K thus : 
  
  if(is.null(lambda01)|is.null(lambda02)|is.null(lambda12)){
    stop("You need to specify lambda on all transition")}
  
  # keep in list
  model<-list()
  length(model)<-K
  
  # be careful to environment of formula 
  environment(formula01)<-environment(formula02)<-environment(formula12)<-environment()
  
  for(k in 1:K){
    set.seed(seed+k)
    if(resampling == "subsampling"){
      # sample 
      id<-sample(1:dim(data)[1],size=dim(data)[1]*tau)
      subdata<-data[id,]
      #run idm 
      # if start is specified or first ite 
      if(k>1 & is.null(B)){
        if(method=="Weib"){
          B<-c(model.idm$modelPar[,1],model.idm$coef[,1])
        }else{B<-c(model.idm$theta01[,1],model.idm$theta02[,1],model.idm$theta12[,1],model.idm$coef[,1])}
      }
      model.idm<-SmoothHazardoptim9::idm(formula01=formula01,
                 formula02=formula02,
                 formula12=formula12,
                 data=subdata,
                 method=method,
                 scale.X=scale.X,
                 maxiter=maxiter,
                 maxiter.pena=maxiter.pena,
                 eps=eps,
                 eps.spline=eps.spline,
                 eps.eigen=eps.eigen,
                 n.knots=n.knots,
                 knots=knots,
                 type.quantile=type.quantile,
                 subset=subset,
                 na.action =na.action,
                 B=B,
                 posfix=posfix,
                 gauss.point=gauss.point,
                 lambda01=lambda01,
                 lambda02=lambda02,
                 lambda12=lambda12,
                 nlambda01=nlambda01,
                 nlambda02=nlambda02,
                 nlambda12=nlambda12,
                 penalty=penalty,
                 penalty.factor=penalty.factor,
                 alpha=alpha,
                 step.sequential=F,
                 option.sequential=list(cutoff=10^-3,
                                             min=20,
                                             step=10),
                      
                 nproc=nproc,
                 clustertype=clustertype)
      
    }else{
        
      # sample boostrap
      id<-sample(1:dim(data)[1],size=dim(data)[1],replace = T)
      subdata<-data[id,]
      #run idm 
      # if start is specified or first ite 
      if(k>1 & is.null(B)){
        if(method=="Weib"){
          B<-c(model.idm$modelPar[,1],model.idm$coef[,1])
        }else{B<-c(model.idm$theta01[,1],model.idm$theta02[,1],model.idm$theta12[,1],model.idm$coef[,1])}
      }
      model.idm<-SmoothHazardoptim9::idm(formula01=formula01,
                     formula02=formula02,
                     formula12=formula12,
                     data=subdata,
                     method=method,
                     scale.X=scale.X,
                     maxiter=maxiter,
                     maxiter.pena=maxiter.pena,
                     eps=eps,
                     eps.spline=eps.spline,
                     eps.eigen=eps.eigen,
                     n.knots=n.knots,
                     knots=knots,
                     type.quantile=type.quantile,
                     subset=subset,
                     na.action =na.action,
                     B=B,
                     posfix=posfix,
                     gauss.point=gauss.point,
                     lambda01=lambda01,
                     lambda02=lambda02,
                     lambda12=lambda12,
                     nlambda01=nlambda01,
                     nlambda02=nlambda02,
                     nlambda12=nlambda12,
                     penalty=penalty,
                     penalty.factor=penalty.factor,
                     alpha=alpha,
                     step.sequential=F,
                     option.sequential=list(cutoff=10^-3,
                                            min=20,
                                            step=10),
                     
                     nproc=nproc,
                     clustertype=clustertype)
    }
    model[[k]]<-model.idm
    
      
  }
  
  
  BIC<-GCV<-CV<-matrix(NA,ncol=K,nrow=length(model[[1]]$BIC))
  pi<-score<-rep(NA,length(model[[1]]$BIC))
  
  for(k in 1:K){
    BIC[,k]<-model[[k]]$BIC
    GCV[,k]<-model[[k]]$GCV
    CV[,k]<-model[[k]]$converged
  }
  
  # same lambda for all
  lambda01<-model[[1]]$lambda[1,]
  lambda02<-model[[1]]$lambda[2,]
  lambda12<-model[[1]]$lambda[3,]
  
  meanBIC<-sapply(c(1:dim(BIC)[1]),FUN=function(x){
    y<-BIC[x,]
    y<-y[CV[x,]==1]
    mean(y)
  })
  
  meanGCV<-sapply(c(1:dim(GCV)[1]),FUN=function(x){
    y<-GCV[x,]
    y<-y[CV[x,]==1]
    mean(y)
  })
  
  optBIC<-list(lambda01=lambda01[which.min(meanBIC)],
               lambda02=lambda02[which.min(meanBIC)],
               lambda12=lambda12[which.min(meanBIC)],
               BIC=min(meanBIC))
  optGCV<-list(lambda01=lambda01[which.min(meanGCV)],
               lambda02=lambda02[which.min(meanGCV)],
               lambda12=lambda12[which.min(meanGCV)],
               GCV=min(meanGCV))
  
  init<-0
  pi_list <-seq(0.6, 0.9, by = 0.01)
  if(calscore==T){
    
    #selprop proportions de fois où le coef est sélectionné
    # do a loop on lambda
    for (m in 1:length(lambda01)[1]){
      coef<-rep(0,length(model[[1]]$coef[,m]))
      names(coef)<-rownames(model[[1]]$coef)
      for(k in 1:K){
        coef<-coef+ifelse(model[[k]]$coef[,m]!=0,1,0)
      }
      selprop<-coef/K
    
      if(any(selprop!=1)){
      value<-StabilityScore(selprop, 
                               pi_list = pi_list,
                               K=K, 
                               n_cat = 3)
      
      id_NA<-which(!is.na(value))
      value<-value[id_NA]
      pi_seq<-pi_list[id_NA]
      pi[m]<-pi_seq[which.min(value)]
      score[m]<-min(value)
      if(init==0){
        opt<-score[m]
        prop<-selprop
        init<-1
        }else{
          if(score[m]<opt){
            opt<-score[m]
            prop<-selprop
          }
        }
      
      }else{score[m]<-pi[m]<-NA}
      
  
        
    }
      
    
    optCal<-list(score=score,
                 pi=pi,
                 selprop=prop,
                 lambda01=lambda01,
                 lambda02=lambda02,
                 lambda12=lambda12
                 )
    
  }else{
    optCal<-NULL
  }
  

  
  return(list(optGCV=optGCV,
              optBIC=optBIC,
              optCal=optCal,
              model=model))
    # 
}
