### Code:
##' @title Log likelihood with M-splines baseline risk
##' @param b  parameters not fixed
##' @param npm  number of parameters not fixed
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
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
DYNderiva<-function(h,b,npm,npar,bfix,fix,zi01,zi02,zi12,ctime,no,nz01,nz02,nz12,ve01,ve02,ve12,
                         dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                         t0,t1,t2,t3,troncature,y01,y02,y12,
                         p01,p02,p12,dimp01,dimp02,dimp12,Ntime){
  res<-rep(0,npm)
  fu<-.Fortran("firstderivaidmlikelihoodsplinetimedep",
           ## input
           as.double(b),
           as.integer(npm),
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
           as.double(y01),
           as.double(y02),
           as.double(y12),
           as.integer(p01),
           as.integer(p02),
           as.integer(p12),
           as.integer(dimp01),
           as.integer(dimp02),
           as.integer(dimp12),
           as.integer(Ntime),
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
           likelihood_res=as.double(res),
           PACKAGE="HIDeM")$likelihood_res
  npm <- length(b)
  bh2 <- bh <- rep(0,npm)
  v <- matrix(0,npm,npm)
  fcith <- fcith2 <- rep(0,npm)
  ## function 
  eps<-h^(1/2)
  h <- sapply(b,function(x){max(eps,(eps*abs(x)))})
  
  l<-0
  ## gradient null
  for(i in 1:npm){
    l <- l+1
    bh <- bh2 <- b
    bh[i] <- bh[i] + h[i]
    bh2[i] <- bh2[i] - h[i]
    
    fcith <- .Fortran("firstderivaidmlikelihoodsplinetimedep",
                      ## input
                      as.double(bh),
                      as.integer(npm),
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
                      as.double(y01),
                      as.double(y02),
                      as.double(y12),
                      as.integer(p01),
                      as.integer(p02),
                      as.integer(p12),
                      as.integer(dimp01),
                      as.integer(dimp02),
                      as.integer(dimp12),
                      as.integer(Ntime),
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
                      likelihood_res=as.double(res),
                      PACKAGE="HIDeM")$likelihood_res
    
    fcith2 <- .Fortran("firstderivaidmlikelihoodsplinetimedep",
                       ## input
                       as.double(bh2),
                       as.integer(npm),
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
                       as.double(y01),
                       as.double(y02),
                       as.double(y12),
                       as.integer(p01),
                       as.integer(p02),
                       as.integer(p12),
                       as.integer(dimp01),
                       as.integer(dimp02),
                       as.integer(dimp12),
                       as.integer(Ntime),
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
                       likelihood_res=as.double(res),
                       PACKAGE="HIDeM")$likelihood_res
    
    v[i,] <- -(fcith-fcith2)/(2*h[i])
  }
  
  v<-v[upper.tri(v,diag=T)]
  
  
  
  result <- list(v=c(v,fu))
  return(result)
}


DYNderivadiag<-function(h,b,npm,npar,bfix,fix,zi01,zi02,zi12,ctime,no,nz01,nz02,nz12,ve01,ve02,ve12,
                    dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                    t0,t1,t2,t3,troncature,y01,y02,y12,
                    p01,p02,p12,dimp01,dimp02,dimp12,Ntime){
  res<-rep(0,npm)
  fu<-.Fortran("firstderivaidmlikelihoodsplinetimedep",
               ## input
               as.double(b),
               as.integer(npm),
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
               as.double(y01),
               as.double(y02),
               as.double(y12),
               as.integer(p01),
               as.integer(p02),
               as.integer(p12),
               as.integer(dimp01),
               as.integer(dimp02),
               as.integer(dimp12),
               as.integer(Ntime),
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
               likelihood_res=as.double(res),
               PACKAGE="HIDeM")$likelihood_res
  npm <- length(b)
  bh2 <- bh <- rep(0,npm)
  v <- rep(0,npm)
  fcith <- fcith2 <- rep(0,npm)
  ## function 
  eps<-h^(1/2)
  h <- sapply(b,function(x){max(eps,(eps*abs(x)))})
  
  l<-0
  ## gradient null
  for(i in 1:npm){
    l <- l+1
    bh <- bh2 <- b
    bh[i] <- bh[i] + h[i]
    bh2[i] <- bh2[i] - h[i]
    
    fcith <- .Fortran("firstderivaidmlikelihoodsplinetimedep",
                      ## input
                      as.double(bh),
                      as.integer(npm),
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
                      as.double(y01),
                      as.double(y02),
                      as.double(y12),
                      as.integer(p01),
                      as.integer(p02),
                      as.integer(p12),
                      as.integer(dimp01),
                      as.integer(dimp02),
                      as.integer(dimp12),
                      as.integer(Ntime),
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
                      likelihood_res=as.double(res),
                      PACKAGE="HIDeM")$likelihood_res
    
    fcith2 <- .Fortran("firstderivaidmlikelihoodsplinetimedep",
                       ## input
                       as.double(bh2),
                       as.integer(npm),
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
                       as.double(y01),
                       as.double(y02),
                       as.double(y12),
                       as.integer(p01),
                       as.integer(p02),
                       as.integer(p12),
                       as.integer(dimp01),
                       as.integer(dimp02),
                       as.integer(dimp12),
                       as.integer(Ntime),
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
                       likelihood_res=as.double(res),
                       PACKAGE="HIDeM")$likelihood_res
    
    v[i] <- -(fcith[i]-fcith2[i])/(2*h[i])
  }
  
  
  result <- list(v=c(v,fu))
  return(result)
}





