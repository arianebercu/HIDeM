### idmlLikelihood.R ---
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
##' @title idm log likelihood
##' @param b  parameters not fixed
##' @param np  number of parameters not fixed
##' @param npar number of parameters
##' @param bfix parameters fixed
##' @param fix indicators of fixed and unfixed parameters
##' @param zi01 knots of transition 0 --> 1
##' @param zi02 knots of transition 0 --> 2
##' @param zi12 knots of transition 1 --> 2
##' @param ctime classification of subject according to their observations
##' @param no number of subjects
##' @param nz01 number of knots for transition 0 -->1 
##' @param nz02 number of knots for transition 0 -->2
##' @param nz12 number of knots for transition 1 -->2
##' @param ve01 variables for transition 0 -->1 
##' @param ve02 variables for transition 0 -->2
##' @param ve12 variables for transition 1 -->2
##' @param dimnva01 number of variables for transition 0 -->1 
##' @param dimnva02 number of variables for transition 0 -->2
##' @param dimnva12 number of variables for transition 1 -->2
##' @param nva01 number of variables for transition 0 -->1 
##' @param nva02 number of variables for transition 0 -->2
##' @param nva12 number of variables for transition 1 -->2
##' @param t0 time entry
##' @param t1 time L
##' @param t2 time R
##' @param t3 time of event/out
##' @param troncature indicator if troncature or not
##' @param gausspoint number of gausspoint quadrature
#' @useDynLib SmoothHazardoptim9
##' @export
#' @author R: Ariane Bercu, Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr> and Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' Fortran: Pierre Joly <Pierre.Joly@@isped.u-bordeaux2.fr>
#' 

idmlLikelihood<-function(b,np,npar,bfix,fix,zi01,zi02,zi12,ctime,no,nz01,nz02,nz12,ve01,ve02,ve12,
                         dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                         t0,t1,t2,t3,troncature,gausspoint){
  res<-0
  #browser()
  .Fortran("idmlikelihood",
           ## input
           as.double(b),
           as.integer(np),
           as.integer(npar),
           as.double(bfix),
           as.integer(fix),
           as.double(zi01),
           as.double(zi12),
           as.double(zi02),
           as.integer(ctime),
           as.integer(no),
           as.integer(nz01),
           as.integer(nz12),
           as.integer(nz02),
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
           as.integer(gausspoint),
           likelihood_res=as.double(res),
           PACKAGE="SmoothHazardoptim9")$likelihood_res
}


