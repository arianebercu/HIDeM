### Code:
#' @title Fit an illness-death model using sampling or boostrap
#' @param K number of repetitions of sampling or boostrap, by default values 100 
#' @param tau Percentage of observation to use for resampling can go u
#' @param seed value of the seed to initialise the random number generator and ensure reproducibility of the results
#' @param resampling approach. Possible values are: "subsampling" for sampling without replacement of a proportion tau of the observations, or "bootstrap" for sampling with replacement generating a resampled dataset with as many observations as in the full sample. Alternatively, this argument can be a function to use for resampling. This function must use arguments named data and tau and return the IDs of observations to be included in the resampled dataset.
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
#'
#' The boostrap returns a list of K elements of type idm, see idm for more details 
#' @seealso \code{\link{idm}}
#' @examples
#' 
#' \dontrun{
#'  library(lava)
#'  library(prodlim)
#'  set.seed(17)
#'  d <- simulateIDM(n=1000)$data
#'  fitweib <- bootstrap.idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2,
#'  formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
#'  formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d,penalty="none",K=3)
#'  print(fitweib)
#' 
#' }
##' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
#' @keywords illness-death
#' @importFrom prodlim Hist
#' @useDynLib SmoothHazardoptim9
#' @export
bootstrap.penidm <- function(
                K=100, # number of sample,
                tau = 0.5, #% in each fold
                seed = 1, # seed
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
                step.sequential=F,
                option.sequential=list(cutoff=10^-3,
                                       min=20,
                                       step=10),
                
                nproc=1,
                clustertype="FORK"){
  
  if((!class(K)%in%c("integer","numeric"))|K<=1){stop("K need to be an integer higher than 1")}
  K<-round(K) # in case if not integer
  if(missing(data)) stop("Need a data frame.")
  if(sum(is.na(data))>0)stop("Need a data frame with no missing data.")
  
  if(!inherits(data,"data.frame"))stop("Argument 'data' must be a data.frame")
  
  if(!penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet","mcp","scad","none")){
    stop(paste0("Parameter penalty must be either : none, lasso, ridge, elastic.net, corrected.elasticnet, mcp or scad"))}
  # resampling check value 
  
  if(!resampling%in%c("subsampling","bootstrap")){
    stop("resampling need to be either bootstrap or subsampling")
  }
  # if we do resampling and calculate calibration we need to have same 
  # lambda for each K thus : 
  
  if(penalty!="none" & (is.null(lambda01)|is.null(lambda02)|is.null(lambda12))){
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
        if(penalty!="none"){
          if(method=="Weib"){
            B<-ifelse(dim(model.idm$coef)[2]==1,
                      c(model.idm$modelPar,model.idm$coef[,1]),
                          c(model.idm$modelPar[,1],model.idm$coef[,1]))
          }else{B<-ifelse(dim(model.idm$coef)[2]==1,
                          c(model.idm$theta01,model.idm$theta02,model.idm$theta12,model.idm$coef[,1]),
                          c(model.idm$theta01[,1],model.idm$theta02[,1],model.idm$theta12[,1],model.idm$coef[,1]))}
        }else{
          if(method=="Weib"){
            B<-c(model.idm$modelPar,model.idm$coef)
          }else{B<-c(model.idm$theta01,model.idm$theta02,model.idm$theta12,model.idm$coef)}
        }
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
                 n.knots=n.knots,
                 knots=knots,
                 type.quantile=type.quantile,
                 na.action =na.action,
                 B=B,
                 posfix=posfix,
                 gausspoint=gausspoint,
                 lambda01=lambda01,
                 lambda02=lambda02,
                 lambda12=lambda12,
                 nlambda01=nlambda01,
                 nlambda02=nlambda02,
                 nlambda12=nlambda12,
                 penalty=penalty,
                 penalty.factor=penalty.factor,
                 alpha=alpha,
                 step.sequential=step.sequential,
                 option.sequential=option.sequential,
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
          B<-ifelse(dim(model.idm$coef)[2]==1,
                    c(model.idm$modelPar,model.idm$coef[,1]),
                    c(model.idm$modelPar[,1],model.idm$coef[,1]))
        }else{B<-ifelse(dim(model.idm$coef)[2]==1,
                        c(model.idm$theta01,model.idm$theta02,model.idm$theta12,model.idm$coef[,1]),
                        c(model.idm$theta01[,1],model.idm$theta02[,1],model.idm$theta12[,1],model.idm$coef[,1]))}
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
                     n.knots=n.knots,
                     knots=knots,
                     type.quantile=type.quantile,
                     na.action =na.action,
                     B=B,
                     posfix=posfix,
                     gausspoint=gausspoint,
                     lambda01=lambda01,
                     lambda02=lambda02,
                     lambda12=lambda12,
                     nlambda01=nlambda01,
                     nlambda02=nlambda02,
                     nlambda12=nlambda12,
                     penalty=penalty,
                     penalty.factor=penalty.factor,
                     alpha=alpha,
                     step.sequential=step.sequential,
                     option.sequential=option.sequential,
                     nproc=nproc,
                     clustertype=clustertype)
    }
    model[[k]]<-model.idm
    
      
  }
  
  
  return(model)
    # 
}
