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
#' @export

derivaweib<-function(b,npm,npar,bfix,fix,ctime,no,ve01,ve02,ve12,
                     dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                     t0,t1,t2,t3,troncature){
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
           likelihood_deriv=as.double(res),
           PACKAGE="SmoothHazardoptim9")$likelihood_deriv
}

derivaweibdiag<-function(b,npm,npar,bfix,fix,ctime,no,ve01,ve02,ve12,
                     dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                     t0,t1,t2,t3,troncature){
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
           likelihood_deriv=as.double(res),
           PACKAGE="SmoothHazardoptim9")$likelihood_deriv
}

derivaweibdiagana<-function(b,npm,npar,bfix,fix,ctime,no,ve01,ve02,ve12,
                            dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                            t0,t1,t2,t3,troncature,gausspoint){
  
  res<-rep(0,npm*2)
  output<-.Fortran("derivaweiballparadiag",
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
                   likelihood_deriv=as.double(res),
                   PACKAGE="SmoothHazardoptim9")$likelihood_deriv
  
  if(any(output==Inf)| any(output==-Inf) | any(is.na(output)) | any(is.nan(output))){
    
    output[any(output==Inf)|any(is.na(output)) | any(is.nan(output))]<-.Machine$double.eps
    output[any(output==-Inf)]<--.Machine$double.eps
    
  }

  nweib<-sum(fix[1:6]==0)
  
  if(nweib>0){
    
    output[(npm+1):(npm+nweib)]<-4*b[1:nweib]*output[(npm+1):(npm+nweib)]*b[1:nweib]+output[1:nweib]*2
    output[1:nweib]<-2*output[1:nweib]*b[1:nweib]
    
  }
  return(output)
  
}


derivaweibana<-function(b,npm,npar,bfix,fix,ctime,no,ve01,ve02,ve12,
                        dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                        t0,t1,t2,t3,troncature,gausspoint){

  res<-rep(0,npm+npm*(npm+1)/2)
  output<-.Fortran("derivaweiballpara",
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
                   likelihood_deriv=as.double(res),
                   PACKAGE="SmoothHazardoptim9")$likelihood_deriv
  
  
  if(any(output==Inf)| any(output==-Inf) | any(is.na(output)) | any(is.nan(output))){
    
    output[any(output==Inf)|any(is.na(output)) | any(is.nan(output))]<-.Machine$double.eps
    output[any(output==-Inf)]<--.Machine$double.eps
    
  }
  nweib<-sum(fix[1:6]==0)
  min<-npm
  max<-min+nweib*(nweib+1)/2+nweib*(npm-nweib)
  Vweib<-matrix(0,npm,npm)
  
  fu<-output[1:min]
  if(nweib>0){
    
    val<-c(output[(min+1):(max)],rep(0,(npm-nweib)*(npm-nweib+1)/2))
    Vweib[lower.tri(Vweib,diag=TRUE)] <- val
    Vweib[1:nweib,1:nweib]<-4*matrix(b[which(fix[1:6]==0)],ncol=1)%*%b[which(fix[1:6]==0)]*Vweib[1:nweib,1:nweib]+diag(output[1:nweib])*2
    
    #try put no CV in BETA and weibull
    Vweib[(nweib+1):npm,1:nweib]<-2*matrix(rep(b[which(fix[1:6]==0)],npm-nweib),ncol=nweib,byrow=T)*
    Vweib[(nweib+1):npm,1:nweib]
    #Vweib[(nweib+1):npm,1:nweib]<-0
    fu[1:nweib]<-2*b[which(fix[1:6]==0)]*fu[1:nweib]
    
  }
  
  min<-max
  V01<- matrix(0,nva01,nva01)
  V01[lower.tri(V01,diag=TRUE)] <- output[(min+1):(min+nva01*(nva01+1)/2)]
  
  
  min<-min+(nva01*(nva01+1)/2)
  if(nva01>0&nva02>0){
    V0102<- matrix(data=output[(min+1):(min+nva02*nva01)],
                   nrow=nva02,ncol=nva01)
  }
  
  min<-min+nva02*nva01
  
  if(nva01>0&nva12>0){
    V0112<- matrix(data=output[(min+1):(min+nva12*nva01)],
                   nrow=nva12,ncol=nva01)
  }
  
  
  min<-min+nva12*nva01
  V02<- matrix(0,nva02,nva02)
  V02[lower.tri(V02,diag=TRUE)] <- output[(min+1):(min+nva02*(nva02+1)/2)]
  
  
  min<-min+(nva02*(nva02+1)/2)
  
  if(nva02>0&nva12>0){
    V0212<- matrix(data=output[(min+1):(min+nva12*nva02)],
                   nrow=nva12,ncol=nva02)
  }
  
  
  min<-min+nva12*nva02
  V12<- matrix(0,nva12,nva12)
  V12[lower.tri(V12,diag=TRUE)] <- output[(min+1):length(output)]
  
  
  
  if(nva01>0){
    Vweib[(nweib+1):(nweib+nva01),(nweib+1):(nweib+nva01)]<-V01
    if(nva02>0){
      Vweib[(nweib+nva01+1):(nweib+nva01+nva02),(nweib+1):(nweib+nva01)]<-V0102
    }
    if(nva12>0){
      Vweib[(nweib+nva01+nva02+1):npm,(nweib+1):(nweib+nva01)]<-V0112
    }
  }
  if(nva02>0){
    Vweib[(nva01+nweib+1):(nva01+nva02+nweib),(nva01+nweib+1):(nva01+nva02+nweib)]<-V02
    if(nva12>0){
      Vweib[(nva01+nva02+nweib+1):npm,(nva01+nweib+1):(nva01+nva02+nweib)]<-V0212
    }
  }
  
  if(nva12>0){
    Vweib[(nva01+nva02+nweib+1):npm,(nva01+nva02+1+nweib):npm]<-V12
  }
  

  # hessian is - second derivatives 
  #V<-V+t(V)
  #diag(V)<-diag(V)/2
  # hessian is - second derivatives 
  #browser()
  return(list(v=-t(Vweib)[upper.tri(Vweib,diag=T)],
         fu=fu))
  
}
