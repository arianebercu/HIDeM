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
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
#' @export

DYNderivaweib<-function(h,b,npm,npar,bfix,fix,ctime,no,ve01,ve02,ve12,
                         dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                         t0,t1,t2,t3,troncature,
                         y01,y02,y12,p01,p02,p12,
                         dimp01,dimp02,dimp12,Ntime){
  
  res<-rep(0,npm)
  fu<-.Fortran("firstderivaidmlikelihoodweibtimedep",
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
      
      fcith <- .Fortran("firstderivaidmlikelihoodweibtimedep",
                           ## input
                           as.double(bh),
                           as.integer(npm),
                           as.integer(npar),
                           as.double(bfix),
                           as.integer(fix),
                           as.integer(ctime),
                           as.integer(no),
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
      fcith2 <- .Fortran("firstderivaidmlikelihoodweibtimedep",
                            ## input
                            as.double(bh2),
                            as.integer(npm),
                            as.integer(npar),
                            as.double(bfix),
                            as.integer(fix),
                            as.integer(ctime),
                            as.integer(no),
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




DYNderivaweibdiag<-function(h,b,npm,npar,bfix,fix,ctime,no,ve01,ve02,ve12,
                        dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                        t0,t1,t2,t3,troncature,
                        y01,y02,y12,p01,p02,p12,
                        dimp01,dimp02,dimp12,Ntime){
  
  res<-rep(0,npm)
  fu<-.Fortran("firstderivaidmlikelihoodweibtimedep",
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
    
    fcith <- .Fortran("firstderivaidmlikelihoodweibtimedep",
                      ## input
                      as.double(bh),
                      as.integer(npm),
                      as.integer(npar),
                      as.double(bfix),
                      as.integer(fix),
                      as.integer(ctime),
                      as.integer(no),
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
    fcith2 <- .Fortran("firstderivaidmlikelihoodweibtimedep",
                       ## input
                       as.double(bh2),
                       as.integer(npm),
                       as.integer(npar),
                       as.double(bfix),
                       as.integer(fix),
                       as.integer(ctime),
                       as.integer(no),
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




