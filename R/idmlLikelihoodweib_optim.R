
### Code:
##' @title Log likelihood with weibull baseline risk
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
##' @param t0 time entry
##' @param t1 time L
##' @param t2 time R
##' @param t3 time of event/out
##' @param troncature indicator if troncature or not
##' @param gausspoint number of gausspoint quadrature
#' @useDynLib SmoothHazardoptim9
##' @export
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
#' 
idmlLikelihoodweiboptim<-function(b,npm,npar,bfix,fix,ctime,no,ve01,ve02,ve12,
                                   dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                                   t0,t1,t2,t3,troncature,lambda,alpha,penalty.factor,penalty,gausspoint,weib){
  

  start<-sum(fix[1:6]==1)
  if(start==6){
  b<-b[1:npm]-b[(1+npm):(2*npm)]
  }else{
    b<-c(b[1:(6-start)],b[(6-start+1):(npm)]-b[(npm+1):(2*(npm-6+start)+6-start)])
  }
  b0<-b

  res<-idmlLikelihoodweib(b,npm,npar,bfix,fix,ctime,no,ve01,ve02,ve12,
                          dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                          t0,t1,t2,t3,troncature,gausspoint,weib)
  
  b<-rep(NA,npar)
  b[fix==0]<-b0
  b[fix==1]<-bfix
  
  if(nva01>0){
    b01<-b[(6+1):(6+nva01)][penalty.factor[1:nva01]==1]
  }else{b01<-0}
  
  if(nva02>0){
    b02<-b[(6+1+nva01):(6+nva01+nva02)][penalty.factor[(nva01+1):(nva01+nva02)]==1]
  }else{b02<-0}
  
  if(nva12>0){
    b12<-b[(6+1+nva01+nva02):npar][penalty.factor[(nva01+nva02+1):(nva01+nva02+nva12)]==1]
  }else{b12<-0}
  # lpen = l-pen
    res<-res-lambda[,1]*(1-alpha)*sum(b01*b01)
    res<-res-lambda[,2]*(1-alpha)*sum(b02*b02)
    res<-res-lambda[,3]*(1-alpha)*sum(b12*b12)
  
  
  return(as.double(-res))
}


