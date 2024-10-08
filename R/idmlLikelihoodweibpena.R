
### Code:
##' @title idm log likelihood with weibull baseline risk
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
##' @param lambda matrix of penalty lambda with column one for transiton from 
##' 0 --> 1, column two from 0 --> 2 and column 3 for 1 --> 2
##' @param alpha alpha on all transitions 
##' @param penalty which penalty to consider
##' @param penalty.factor which variable should be penalised
##' @param gausspoint number of points in gauss quadrature
##' @param gausspoint number of gausspoint quadrature
#' @useDynLib SmoothHazardoptim9
##' @export
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
#' 


idmlLikelihoodweibpena<-function(b,npm,npar,bfix,fix,ctime,no,ve01,ve02,ve12,
                         dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                         t0,t1,t2,t3,troncature,lambda,alpha,penalty.factor,penalty,gausspoint){
  res<-0
  b0<-b
  #browser()
  res<-.Fortran("idmlikelihoodweib",
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
           as.integer(gausspoint),
           likelihood_res=as.double(res),
           PACKAGE="SmoothHazardoptim9")$likelihood_res

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
  if(penalty%in%c("lasso","ridge","elasticnet")){
    res<-res-lambda[,1]*alpha*sum(abs(b01))-lambda[,1]*(1-alpha)*sum(b01*b01)
    res<-res-lambda[,2]*alpha*sum(abs(b02))-lambda[,2]*(1-alpha)*sum(b02*b02)
    res<-res-lambda[,3]*alpha*sum(abs(b12))-lambda[,3]*(1-alpha)*sum(b12*b12)
  }

   if(penalty=="mcp"){
  
     p01<-rep(alpha*lambda[,1]*lambda[,1]/2,length(b01))
     idbeta<-which(b01<=alpha*lambda[,1])
     p01[idbeta]<-lambda[,1]*abs(b01[idbeta])-((b01[idbeta]*b01[idbeta])/2*alpha)
  
     p02<-rep(alpha*lambda[,2]*lambda[,2]/2,length(b02))
     idbeta<-which(b02<=alpha*lambda[,2])
     p02[idbeta]<-lambda[,2]*abs(b02[idbeta])-((b02[idbeta]*b02[idbeta])/2*alpha)
  
     p12<-rep(alpha*lambda[,3]*lambda[,3]/2,length(b12))
     idbeta<-which(b12<=alpha*lambda[,3])
     p12[idbeta]<-lambda[,3]*abs(b12[idbeta])-((b12[idbeta]*b12[idbeta])/2*alpha)
  
     res<-res-sum(p01)-sum(p02)-sum(p12)
  
   }

   if(penalty=="scad"){
  
     p01<-rep((lambda[,1]^2)*(alpha+1)/2,length(b01))
     idbeta<-which(b01<=lambda[,1])
     p01[idbeta]<-lambda[,1]*abs(b01[idbeta])
     idbeta<-which(abs(b01)<lambda[,1]*alpha)
     p01[idbeta]<-(2*alpha*lambda[,1]*abs(b01[idbeta])-b01[idbeta]^2-lambda[,1]^2)/(2*(alpha-1))
  
     p02<-rep((lambda[,2]^2)*(alpha+1)/2,length(b02))
     idbeta<-which(b02<=lambda[,2])
     p02[idbeta]<-lambda[,2]*abs(b02[idbeta])
     idbeta<-which(abs(b02)<lambda[,2]*alpha)
     p02[idbeta]<-(2*alpha*lambda[,2]*abs(b02[idbeta])-b02[idbeta]^2-lambda[,2]^2)/(2*(alpha-1))
  
     p12<-rep((lambda[,3]^2)*(alpha+1)/2,length(b12))
     idbeta<-which(b12<=lambda[,3])
     p12[idbeta]<-lambda[,3]*abs(b12[idbeta])
     idbeta<-which(abs(b12)<lambda[,3]*alpha)
     p12[idbeta]<-(2*alpha*lambda[,3]*abs(b12[idbeta])-b12[idbeta]^2-lambda[,3]^2)/(2*(alpha-1))
  
     res<-res-sum(p01)-sum(p02)-sum(p12)
   }


  return(as.double(res))
}

