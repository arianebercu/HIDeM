### Code:
##' @title M-spline or Weibull estimate of the transition intensity function
##' @description 
##' M-spline or Weibull estimate of the transition intensity function
##' and the cumulative transition intensity function
##' for survival and illness-death models
##' Important: the theta parameters returned by \code{idm} are in fact
##' the square root of the splines coefficients. See examples.
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
##'   library(lava)
##' library(prodlim)
##' set.seed(17)
##' d <- simulateIDM(n=1000)$data
##' fitweib <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2,
##'                formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
##'                formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d)
##' intensity(times=fitweib$time[,1],
##'                theta = fitweib$modelPar[1:2]^2)
##' fitsplines <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2,
##'                   formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
##'                   formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d,
##'                   method="splines")
##' intensity(times=fitsplines$time[,1],
##'                theta = fitsplines$theta01^2,
##'                knots = fitsplines$knots01,
##'                number.knots = fitsplines$nknots01,
##'                method = "splines")
##' }
##' @importFrom pracma gauss_kronrod
#' @useDynLib HIDeM
##' @export
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
#' 
intensity <- function(times,knots,number.knots=NULL,theta,linear.predictor=0,V=NULL,
                      fix=NULL,converged=NULL,conf.int = 0.95,
                      method="weib") {
  
  # need to compute intensity for weibull and splines output

  if(method=="weib"){
    if(length(theta)!=2) stop(paste0("The number of theta coefficients must be 2 but the length providded is",
                              length(theta)))
    if(sum(is.na(theta))>0 | sum(is.na(times))>0)stop("No missing data is allowed in theta and times")
    
  }else{
    if(is.null(number.knots))stop("Need to specify the number of knots")
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
            deriv1int<-(times[j]^(theta[1]-1))*(theta[2]^theta[1])*(theta[1]*log(times[j])+theta[1]*log(theta[2])+1)
            deriv1int<-deriv1int*(2*sqrt(theta[1]))
            
            deriv1cumu<-((times[j]*theta[2])^theta[1])*log(times[j]*theta[2])
            deriv1cumu<-deriv1cumu*(2*sqrt(theta[1]))
          }else{deriv1int<-deriv1cumu<-NULL}
          if(fix[2]==0){
            deriv2int<-(theta[1]^2)*(times[j]^(theta[1]-1))*(theta[2]^(theta[1]-1))
            deriv2int<-deriv2int*(2*sqrt(theta[2]))
            
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

      if(is.null(times)){
        times<-seq(from=knots.bound[1],to=knots.bound[2],length.out=100)
      }else{
        times<-times[times>=knots.bound[1] & times<=knots.bound[2]]
        }
      # need to sort knots.int  to have same as msplines
      # do not need to sort times not necessary same output 
      msplines<-splinesMI(x=times,knots=sort(knots.int),Boundary.knots=knots.bound)$Mspline
      isplines<-splinesMI(x=times,knots=sort(knots.int),Boundary.knots=knots.bound)$Ispline
      #splines2::iSpline(x=times,knots=knots.int,Boundary.knots=knots.bound,intercept=T)
      #browser()
      
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
#  library(HIDeM)
#  data(Paq1000)
#  
# #  Illness-death model with certif on the 3 transitions
# #  Weibull parametrization and likelihood maximization
#  
#  fit.weib <- HIDeM::idm(formula02=Hist(time=t,event=death,entry=e)~certif,
#                  formula01=Hist(time=list(l,r),event=dementia)~certif,
#                  data=Paq1000)
#  
#  test<-intensity(times=fit.weib$time,
#            theta=fit.weib$modelPar[1:2]^2,linear.predictor=0,V=fit.weib$V[1:2,1:2],
#            fix=rep(0,2),converged=1,conf.int = 0.95,
#            method="weib")
#  
#  round(test$intensity,4)==round(fit.weib$intensity01,4)
#  test$lowerintensity
#  fit.weib$lowerIntensity01
#  
#  test$upperintensity
#  fit.weib$upperIntensity01
#  
#  fit.s<- HIDeM::idm(formula02=Hist(time=t,event=death,entry=e)~certif,
#                                        formula01=Hist(time=list(l,r),event=dementia)~certif,
#                                        data=Paq1000,method = "splines",n.knots = c(3,3,3))
#  
#  
#  
#  test<-intensity(times=fit.s$time[,1],
#                  theta=fit.s$theta01,
#                  knots=fit.s$knots01,
#                  number.knots=fit.s$nknots01,linear.predictor=0,V=fit.s$V[1:5,1:5],
#                  fix=rep(0,5),converged=1,conf.int = 0.95,
#                  method="splines")
#  
#  fit.idm <-  HIDeM::idm(formula02 = Hist(time = t, event = death, entry = e) ~ certif,
#                  formula01 = Hist(time = list(l,r), event = dementia) ~ certif,
#                  formula12 = ~ certif, method = "Splines", data = Paq1000)
#  
#  test<-intensity(times=fit.idm$time[,1],
#                  theta=fit.idm$theta01,
#                  knots=fit.idm$knots01,
#                  number.knots=fit.idm$nknots01,linear.predictor=0,V=fit.idm$V[1:9,1:9],
#                  fix=rep(0,9),converged=1,conf.int = 0.95,
#                  method="splines")
#  test$intensity
#  fit.idm$intensity01
# 
#  same function as isplines and msplines from packages splines2 
#  --> rewrite to use this as internal 
# 
