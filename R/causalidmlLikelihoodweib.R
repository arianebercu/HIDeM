### idmlLikelihoodweib.R ---
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
##' @param npm  number of parameters not fixed
##' @param npar number of parameters
##' @param bfix parameters fixed
##' @param fix indicators of fixed and unfixed parameters
##' @param ctime classification of subject according to their observations
##' @param no number of subjects
##' @param ve01 variables for transition 0 -->1 
##' @param ve02 variables for transition 0 -->2
##' @param ve12 variables for transition 1 -->2
##' @param dimnva01 number of variables for transition 0 -->1 
##' @param dimnva02 number of variables for transition 0 -->2
##' @param dimnva12 number of variables for transition 1 -->2
##' @param nva01 number of variables for transition 0 -->1 
##' @param nva02 number of variables for transition 0 -->2
##' @param nva12 number of variables for transition 1 -->2
##' @param nva12dep time dep transition 12
##' @param semiMark indicator if semiMarkov or Markov
##' @param t0 time entry
##' @param t1 time L
##' @param t2 time R
##' @param t3 time of event/out
##' @param troncature indicator if troncature or not
#' @useDynLib HIDeM
##' @export
#' @author R: Ariane Bercu, Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr> and Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' Fortran: Pierre Joly <Pierre.Joly@@isped.u-bordeaux2.fr>
#' 

causalidmlLikelihoodweib<-function(b,npm,npar,bfix,fix,ctime,no,ve01,ve02,ve12,
                         dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,nva12dep,
                         semiMark,t0,t1,t2,t3,troncature){
  res<-0

  .Fortran("causalidmlikelihoodweib",
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
           as.integer(nva12dep),
           as.integer(semiMark),
           as.double(t0),
           as.double(t1),
           as.double(t2),
           as.double(t3),
           as.integer(troncature),
           likelihood_res=as.double(res),
           PACKAGE="HIDeM")$likelihood_res
}


