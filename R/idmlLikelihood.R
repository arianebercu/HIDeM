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
##' @param b0  parameters not fixed
##' @param np0  number of parameters not fixed
##' @param npar0 number of parameters
##' @param bfix0 parameters fixed
##' @param fix0 indicators of fixed and unfixed parameters
##' @param zi010 knots of transition 0 --> 1
##' @param zi020 knots of transition 0 --> 2
##' @param zi120 knots of transition 1 --> 2
##' @param c0 classification of subject according to their observations
##' @param no0 number of subjects
##' @param nz01 number of knots for transition 0 -->1 
##' @param nz02 number of knots for transition 0 -->2
##' @param nz12 number of knots for transition 1 -->2
##' @param ve010 variables for transition 0 -->1 
##' @param ve020 variables for transition 0 -->2
##' @param ve120 variables for transition 1 -->2
##' @param dimnva01 number of variables for transition 0 -->1 
##' @param dimnva02 number of variables for transition 0 -->2
##' @param dimnva12 number of variables for transition 1 -->2
##' @param nva01 number of variables for transition 0 -->1 
##' @param nva02 number of variables for transition 0 -->2
##' @param nva12 number of variables for transition 1 -->2
##' @param t00 time entry
##' @param t10 time L
##' @param t20 time R
##' @param t30 time of event/out
##' @param troncature0 indicator if troncature or not
##' @param gausspoint0 number of gausspoint quadrature
#' @useDynLib SmoothHazardoptim9
##' @export
#' @author R: Ariane Bercu, Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr> and Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' Fortran: Pierre Joly <Pierre.Joly@@isped.u-bordeaux2.fr>
#' 

idmlLikelihood<-function(b0,np0,npar0,bfix0,fix0,zi010,zi020,zi120,c0,no0,nz01,nz02,nz12,ve010,ve020,ve120,
                         dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                         t00,t10,t20,t30,troncature0,gausspoint0){
  res<-0
  #browser()
  .Fortran("idmlikelihood",
           ## input
           as.double(b0),
           as.integer(np0),
           as.integer(npar0),
           as.double(bfix0),
           as.integer(fix0),
           as.double(zi010),
           as.double(zi120),
           as.double(zi020),
           as.integer(c0),
           as.integer(no0),
           as.integer(nz01),
           as.integer(nz12),
           as.integer(nz02),
           as.double(ve010),
           as.double(ve120),
           as.double(ve020),
           as.integer(dimnva01),
           as.integer(dimnva12),
           as.integer(dimnva02),
           as.integer(nva01),
           as.integer(nva12),
           as.integer(nva02),
           as.double(t00),
           as.double(t10),
           as.double(t20),
           as.double(t30),
           as.integer(troncature0),
           as.integer(gausspoint0),
           likelihood_res=as.double(res),
           PACKAGE="SmoothHazardoptim9")$likelihood_res
}


