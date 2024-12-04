### Code:
##' @title First and second order derivatives of the log-likelihood with weibull baseline risk
##' @param b  parameters on explanatory variables not fixed
##' @param bfix  parameters on explanatory variables fixed
##' @param npm  number of parameters not fixed, thus length of b
##' @param npar  number of parameters
##' @param fix indicator of length npar, values 1 if parameter fixed
##' @param ctime profile of patients from 1 to 7
##' @param no number of subjects 
##' @param ve01 variables for transition 0 -->1 
##' @param ve02 variables for transition 0 -->2
##' @param ve12 variables for transition 1 -->2
##' @param dimnva01 number of variables for transition 0 -->1, if not variables value 1
##' @param dimnva02 number of variables for transition 0 -->2, if not variables value 1
##' @param dimnva12 number of variables for transition 1 -->2, if not variables value 1
##' @param nva01 number of variables for transition 0 -->1, if not variables value 0
##' @param nva02 number of variables for transition 0 -->2, if not variables value 0
##' @param nva12 number of variables for transition 1 -->2, if not variables value 0
##' @param t0 time of entry
##' @param t1 time of last visit or last visit without diagnose of illness
##' @param t2 time of last visit or time diagnose of illness
##' @param t3 time of last visit or death
##' @param troncature indicator of troncature, value 1 if there is troncature otherwise 0.
#' @useDynLib SmoothHazardoptim9
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 

derivaweib<-function(b,npm,npar,bfix,fix,ctime,no,ve01,ve02,ve12,
                     dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                     t0,t1,t2,t3,troncature,weib){
  res<-rep(0,(npm*(npm+1)/2)+npm)
  # return first and second derivatives of the loglik
  .Fortran("derivaweib",
           ## input
           as.double(b),
           as.integer(npm),
           as.integer(npar),
           as.double(bfix),
           as.integer(fix),
           as.integer(ctime),
           as.integer(no),
           as.double(ve01),
           as.double(ve12),
           as.double(ve02),
           as.integer(dimnva01),
           as.integer(dimnva12),
           as.integer(dimnva02),
           as.integer(nva01),
           as.integer(nva12),
           as.integer(nva02),
           as.double(t0),
           as.double(t1),
           as.double(t2),
           as.double(t3),
           as.integer(troncature),
           as.integer(weib),
           likelihood_deriv=as.double(res),
           PACKAGE="SmoothHazardoptim9")$likelihood_deriv
}

derivaweibdiag<-function(b,npm,npar,bfix,fix,ctime,no,ve01,ve02,ve12,
                     dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                     t0,t1,t2,t3,troncature,weib){
  res<-rep(0,2*npm)
  # return first and second derivatives of the loglik
 .Fortran("derivaweibdiag",
           ## input
           as.double(b),
           as.integer(npm),
           as.integer(npar),
           as.double(bfix),
           as.integer(fix),
           as.integer(ctime),
           as.integer(no),
           as.double(ve01),
           as.double(ve12),
           as.double(ve02),
           as.integer(dimnva01),
           as.integer(dimnva12),
           as.integer(dimnva02),
           as.integer(nva01),
           as.integer(nva12),
           as.integer(nva02),
           as.double(t0),
           as.double(t1),
           as.double(t2),
           as.double(t3),
           as.integer(troncature),
           as.integer(weib),
           likelihood_deriv=as.double(res),
           PACKAGE="SmoothHazardoptim9")$likelihood_deriv
}


