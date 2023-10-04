### intensity.R ---
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Feb  6 2016 (08:47)
## Version:
## last-updated: Feb 25 2016 (13:17)
##           By: Thomas Alexander Gerds
##     Update #: 27
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
##'  M-spline estimate of the transition intensity function
##' and the cumulative transition intensity function
##' for survival and illness-death models
##'
##' The estimate of the transition intensity function is a linear
##' combination of M-splines and the estimate of the cumulative transition
##' intensity function is a linear combination of I-splines (the integral of a
##' M-spline is called I-spline). The coefficients \code{theta} are the same for
##' the M-splines and I-splines.
##'
##' Important: the theta parameters returned by \code{idm} are in fact
##' the square root of the splines coefficients. See examples.
##'
##' This function is a R-translation of a corresponding Fortran function called \code{susp}. \code{susp} is
##' used internally by \code{idm}.
##'
##' @title M-spline estimate of the transition intensity function
##' @param times Time points at which to estimate the intensity function
##' @param knots Knots for the M-spline
##' @param number.knots Number of knots for the M-splines (and I-splines see details)
##' @param theta The coefficients for the linear combination of M-splines (and I-splines see details), theta from the model
##' @param linear.predictor Linear predictor beta*Z. When it is non-zero,
##' transition and cumulative transition are multiplied by \code{exp(linear.predictor)}. Default is zero.
##' @param V matrix of variance co-variance
##' @param fix indicators if parameters are fixed 
##' @param converged indicator of convergence of the models
##' @param conf.int confidence intervals, 1 - alpha
##' @param method the methodology, splines or weib
##' @return
##' \item{times}{The time points at which the following estimates are evaluated.}
##' \item{intensity}{The transition intensity function evaluated at \code{times}.}
##' \item{cumulative.intensity}{The cumulative transition intensity function evaluated at \code{times}}
##' \item{survival}{The "survival" function, i.e., exp(-cumulative.intensity)}
##'
##' @examples
##'
##'
##' \dontrun{
##'            
##'   data(Paq1000)
##'   fit.idm <-  idm(formula02 = Hist(time = t, event = death, entry = e) ~ certif,
##'                   formula01 = Hist(time = list(l,r), event = dementia) ~ certif,
##'                   formula12 = ~ certif, method = "Splines", data = Paq1000)
##'   # Probability of survival in state 0 at age 80 for a subject with no cep given
##'     that he is in state 0 at 70
##'   su0 <- (intensity(times = 80, knots = fit.idm$knots01,
##'                    number.knots = fit.idm$nknots01,
##'                    theta = fit.idm$theta01^2)$survival
##'          *intensity(times = 80, knots = fit.idm$knots02,
##'                    number.knots = fit.idm$nknots02,
##'                    theta = fit.idm$theta02^2)$survival)/
##'         (intensity(times = 70, knots = fit.idm$knots01,
##'                    number.knots = fit.idm$nknots01,
##'                    theta = fit.idm$theta01^2)$survival
##'         *intensity(times = 70, knots = fit.idm$knots02,
##'                    number.knots = fit.idm$nknots02,
##'                    theta = fit.idm$theta02^2)$survival)
##'   # Same result as:
##'   predict(fit.idm, s = 70, t = 80, conf.int = FALSE) # see first element
##' }
#' @useDynLib SmoothHazardoptim9
##' @export
#' @author R: Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr> and Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' Fortran: Pierre Joly <Pierre.Joly@@isped.u-bordeaux2.fr>
#' 
intensity <- function(times,knots,number.knots,theta,linear.predictor=0,V=NULL,
                      fix=NULL,converged=NULL,conf.int = 0.95,
                      method="weib") {
  
  # need to compute intensity for weibull and splines output

  if(method=="weib"){
    if(length(theta)!=2) stop(paste0("The number of theta coefficients must be 2 but the length providded is",
                              length(theta)))
    if(sum(is.na(theta))>0 | sum(is.na(times))>0)stop("No missing data is allowed in theta and times")
    
  }else{
    if (length(theta)!= number.knots+2) stop(paste0("For ",
                  number.knots,
                  " knots we need ",
                  number.knots+2,
                  " coefficients. But, length of argument theta as provided is ",
                  length(theta)))
    if(sum(is.na(theta))>0 | sum(is.na(times))>0 | is.na(number.knots))stop("No missing data is allowed in theta, times, knots or number.knots")
    if(number.knots<3)stop("Number of knots need to be superior or equal to 3")
  }
  
    if(is.null(converged)){
      #warning(paste("Converged was not provided, it suppose that the model did not converged"))
      interval.calcul<-0
      converged<-0
    }
  
  if(converged!=1){
    #warning(paste("As the convergence was not obtained confidence intervals will not be provided"))
    interval.calcul<-0
  }else{interval.calcul<-1}
  
    if(is.null(V)){
    #warning(paste("The covariance matrix is not provided, the calculus of confidence intervals is not possible"))
      interval.calcul<-0
    }
  
   if(is.null(fix)){
    #warning(paste("The indicator on fixed splines parameters was not provided, it suppose that no parameter was fixed"))
    fix<-rep(0,length(theta))
   }
  
  if(length(fix)!=length(theta)){
    stop(paste0("Fixe index must have the same length as theta thus : ",length(theta)))
    
  }
  if(sum(fix)==length(theta)){
    #warning(paste("All parameters are fixed so no interval confidence will be provided"))
    interval.calcul<-0
    
  }
    cumulative.intensity=rep(0,length(times))   # risque cumule
    intensity=rep(0,length(times))  # risque
    survival=rep(0,length(times))   # survie

    theta.square<-theta^2
    
    if(method=="weib"){
      
      
      intensity<-theta[1]*(theta[2]^theta[1])*times^(theta[1]-1)
      cumulative.intensity<-(theta[2]*times)^theta[1]
      
      e = exp(linear.predictor)
      intensity=intensity*e
      cumulative.intensity=cumulative.intensity*e
      survival = exp(-cumulative.intensity)
      
      lowerintensity<-rep(NA,length(times))
      upperintensity<-rep(NA,length(times))
      lowercumulative.intensity<-rep(NA,length(times))
      uppercumulative.intensity<-rep(NA,length(times))
      
      if(interval.calcul==1){
        
        V<-V[fix==0,fix==0]
        
        if(dim(V)[1]!=sum(fix==0))stop("The number of parameters estimated is not equal to sum(fix==1)), need to change V or fix.spline")
        
        for( j in 1:length(times)){
          
          
          if(fix[1]==0){
            #deriv1int<-2*theta.square[1]*(times[j]^(theta.square[1]-1))*(theta.square[2]^theta.square[1])*
            # ((theta.square[1]^2)*(log(times[j])+log(theta.square[2]^2))+1)
            deriv1int<-(times[j]^(theta[1]-1))*(theta[2]^theta[1])*(theta[1]*log(times[j])+theta[1]*log(theta[2])+1)
            deriv1int<-deriv1int*(2*sqrt(theta[1]))
            
            #deriv1cumu<-2*theta[1]*((times[j]*theta.square[2])^theta.square[1])*log(times[j]*theta.square[2])
            deriv1cumu<-((times[j]*theta[2])^theta[1])*log(times[j]*theta[2])
            deriv1cumu<-deriv1cumu*(2*sqrt(theta[1]))
          }else{deriv1int<-deriv1cumu<-NULL}
          if(fix[2]==0){
            #deriv2int<-2*(theta.square[1]^2)*(times[j]^(theta.square[1]-1))*(theta.square[2]^theta.square[1])/theta[2]
            deriv2int<-(theta[1]^2)*(times[j]^(theta[1]-1))*(theta[2]^(theta[1]-1))
            deriv2int<-deriv2int*(2*sqrt(theta[2]))
            
            #deriv2cumu<-2*theta.square[1]*((times[j]*theta.square[2])^theta.square[1])/theta[2]
            deriv2cumu<-(theta[1]*((times[j]*theta[2])^theta[1]))/theta[2]
            deriv2cumu<-deriv2cumu*(2*sqrt(theta[2]))
          }else{deriv2int<-deriv2cumu<-NULL}

          derivint<-c(deriv1int,deriv2int)
          
          derivcumu<-c(deriv1cumu,deriv2cumu)
          
          Vthetaint<-t(derivint)%*%V%*%derivint
          Vthetacumu<-t(derivcumu)%*%V%*%derivcumu
          lowerintensity[j]<-intensity[j]+qnorm((1-conf.int)/2)*sqrt(Vthetaint)
          upperintensity[j]<-intensity[j]-qnorm((1-conf.int)/2)*sqrt(Vthetaint)
          lowercumulative.intensity[j]<-cumulative.intensity[j]+qnorm((1-conf.int)/2)*sqrt(Vthetacumu)
          uppercumulative.intensity[j]<-cumulative.intensity[j]-qnorm((1-conf.int)/2)*sqrt(Vthetacumu)
        }
        }
      }else{
      knots.unique<-unique(knots)
      knots.bound<-knots.unique[c(1,length(knots.unique))]
      
  
      knots.int<-knots.unique[-c(1,length(knots.unique))]
      #browser()
      if(is.null(times)){
        times<-seq(from=knots.bound[1],to=knots.bound[2],length.out=100)}
      
      msplines<-splines2::mSpline(x=times,knots=knots.int,Boundary.knots=knots.bound,intercept = T)
      isplines<-splines2::iSpline(x=times,knots=knots.int,Boundary.knots=knots.bound,intercept = T)
      
      intensity<-msplines%*%theta.square
      cumulative.intensity<-isplines%*%theta.square
      
      
      e = exp(linear.predictor)
      intensity=intensity*e
      cumulative.intensity=cumulative.intensity*e
      survival = exp(-cumulative.intensity)
      
      if(length(times)==1){
        msplines<-matrix(msplines[,fix==0],nrow=1,ncol=sum(fix==0))
        isplines<-matrix(isplines[,fix==0],nrow=1,ncol=sum(fix==0))
      }else{
        msplines<-msplines[,fix==0]
        isplines<-isplines[,fix==0]
      }
      
      
      
      lowerintensity<-rep(NA,dim(msplines)[1])
      upperintensity<-rep(NA,dim(msplines)[1])
      lowercumulative.intensity<-rep(NA,dim(isplines)[1])
      uppercumulative.intensity<-rep(NA,dim(isplines)[1])
  
      if(interval.calcul==1){
        
        V<-V[fix==0,fix==0]
        theta<-theta[fix==0]
        if(dim(V)[1]!=sum(fix==0))stop("The number of parameters estimated is not equal to sum(fix==1)), need to change V or fix.spline")
        Vtheta<-4*theta*V%*%diag(theta)
        
        for( j in 1:dim(msplines)[1]){
          lowerintensity[j]<-intensity[j]+qnorm((1-conf.int)/2)*sqrt(msplines[j,]%*%Vtheta%*%(msplines[j,]))
          upperintensity[j]<-intensity[j]-qnorm((1-conf.int)/2)*sqrt(msplines[j,]%*%Vtheta%*%(msplines[j,]))
          lowercumulative.intensity[j]<-cumulative.intensity[j]+qnorm((1-conf.int)/2)*sqrt(isplines[j,]%*%Vtheta%*%(isplines[j,]))
          uppercumulative.intensity[j]<-cumulative.intensity[j]-qnorm((1-conf.int)/2)*sqrt(isplines[j,]%*%Vtheta%*%(isplines[j,]))
          }
      }
    }

    return(list(times=times,intensity=as.vector(intensity),
                lowerintensity=as.vector(lowerintensity),
                upperintensity=as.vector(upperintensity),
                cumulative.intensity=as.vector(cumulative.intensity),
                lowercumulative.intensity=as.vector(lowercumulative.intensity),
                uppercumulative.intensity=as.vector(uppercumulative.intensity),
                survival=as.vector(survival)))
}

#----------------------------------------------------------------------
### intensity.R ends here
# # 
# library(SmoothHazard)
# data(Paq1000)
# 
# # Illness-death model with certif on the 3 transitions
# # Weibull parametrization and likelihood maximization
# 
# fit.weib <- SmoothHazard::idm(formula02=Hist(time=t,event=death,entry=e)~certif,
#                 formula01=Hist(time=list(l,r),event=dementia)~certif,
#                 data=Paq1000)
# 
# 
# fit.weib.8 <- SmoothHazardoptim8::idm(formula02=Hist(time=t,event=death,entry=e)~certif,
#                 formula01=Hist(time=list(l,r),event=dementia)~certif,
#                 data=Paq1000)
# test<-intensity(times=fit.weib$time,
#           theta=fit.weib.8$modelPar[1:2]^2,linear.predictor=0,V=fit.weib.8$V[1:2,1:2],
#           fix=rep(0,2),converged=1,conf.int = 0.95,
#           method="Weib")
# 
# round(test$intensity,4)==round(fit.weib$intensity01,4)
# test$lowerintensity
# fit.weib$lowerIntensity01
# 
# test$upperintensity
# fit.weib$upperIntensity01
# 
# fit.s<- SmoothHazardoptim8::idm(formula02=Hist(time=t,event=death,entry=e)~certif,
#                                       formula01=Hist(time=list(l,r),event=dementia)~certif,
#                                       data=Paq1000,method = "splines",n.knots = c(3,3,3))
# 
# 
# 
# test<-intensity(times=fit.s$time[,1],
#                 theta=fit.s$theta01,
#                 knots=fit.s$knots01,
#                 number.knots=fit.s$nknots01,linear.predictor=0,V=fit.s$V[1:5,1:5],
#                 fix=rep(0,5),converged=1,conf.int = 0.95,
#                 method="splines")
# 
# fit.idm <-  SmoothHazard::idm(formula02 = Hist(time = t, event = death, entry = e) ~ certif,
#                 formula01 = Hist(time = list(l,r), event = dementia) ~ certif,
#                 formula12 = ~ certif, method = "Splines", data = Paq1000)
# 
# test<-intensity(times=fit.idm$time[,1],
#                 theta=fit.idm$theta01,
#                 knots=fit.idm$knots01,
#                 number.knots=fit.idm$nknots01,linear.predictor=0,V=fit.idm$V[1:9,1:9],
#                 fix=rep(0,9),converged=1,conf.int = 0.95,
#                 method="splines")
# test$intensity
# fit.idm$intensity01
