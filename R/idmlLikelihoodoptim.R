### Code:
##' @title Penalised log likelihood with M-splines baseline risk
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
##' @param lambda matrix of penalty lambda with column one for transiton from 
##' 0 --> 1, column two from 0 --> 2 and column 3 for 1 --> 2
##' @param alpha alpha on all transitions 
##' @param penalty which penalty to consider
##' @param penalty.factor which variable should be penalised
#' @useDynLib SmoothHazardoptim9
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
#'
idmlLikelihoodoptim<-function(b,npm,npar,bfix,fix,zi01,zi02,zi12,ctime,no,nz01,nz02,nz12,ve01,ve02,ve12,
                              dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                              t0,t1,t2,t3,troncature,gausspoint,lambda,alpha,penalty.factor,penalty){
  
  #browser()
  size_spline<-nz01+nz02+nz12+6
  start<-sum(fix[1:(size_spline)]==1)
  if(start==size_spline){
    b0<-b[1:npm]-b[(1+npm):(2*npm)]
  }else{
    b0<-c(b[1:(size_spline-start)],b[(size_spline-start+1):(npm)]-b[(npm+1):(length(b))])
  }
  
  
  res<-idmlLikelihood(b0,npm,npar,bfix,fix,zi01,zi02,zi12,ctime,no,nz01,nz02,nz12,ve01,ve02,ve12,
    dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
    t0,t1,t2,t3,troncature,gausspoint)
  
  ball<-rep(NA,npar)
  ball[which(fix==0)]<-b0
  ball[which(fix==1)]<-bfix
  
  if(nva01>0){
    b01<-ball[(size_spline+1):(size_spline+nva01)][penalty.factor[1:nva01]==1]
    
  }else{b01<-0}
  
  if(nva02>0){
    b02<-ball[(size_spline+1+nva01):(size_spline+nva01+nva02)][penalty.factor[(nva01+1):(nva01+nva02)]==1]
    
  }else{b02<-0}
  
  if(nva12>0){
    b12<-ball[(size_spline+1+nva01+nva02):npar][penalty.factor[(nva01+nva02+1):(nva01+nva02+nva12)]==1]
    
  }else{b12<-0}
  # lpen = l-pen
  res<-res-lambda[,1]*(1-alpha)*sum(b01*b01)-lambda[,1]*alpha*sum(abs(b01))
  res<-res-lambda[,2]*(1-alpha)*sum(b02*b02)-lambda[,2]*alpha*sum(abs(b02))
  res<-res-lambda[,3]*(1-alpha)*sum(b12*b12)-lambda[,3]*alpha*sum(abs(b12))
  
  
  return(as.double(res))
}

groptim<-function(b,npm,npar,bfix,fix,zi01,zi02,zi12,ctime,no,nz01,nz02,nz12,ve01,ve02,ve12,
                  dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                  t0,t1,t2,t3,troncature,gausspoint,lambda,alpha,penalty.factor,penalty){
  
  size_spline<-nz01+nz02+nz12+6
  start<-sum(fix[1:size_spline]==1)
  svar<-b[1:(size_spline-start)]
  # test<-grad(x=b,
  #            f=idmlLikelihoodoptim,
  #            npm=sum(fix==0),
  #            npar=length(fix),
  #            bfix=bfix,
  #            fix=fix,
  #            ctime=ctime,
  #            no=no,
  #            ve01=ve01,
  #            ve02=ve02,
  #            ve12=ve12,
  #            dimnva01=dimnva01,
  #            dimnva02=dimnva02,
  #            dimnva12=dimnva12,
  #            nva01=nva01,
  #            nva02=nva02,
  #            nva12=nva12,
  #            t0=t0,
  #            t1=t1,
  #            t2=t2,
  #            t3=t3,
  #            troncature=troncature,
  #            lambda=lambda,
  #            alpha=alpha,
  #            penalty.factor=penalty.factor,
  #            penalty=penalty,
  #            gausspoint=gausspoint,weib=weib)
  # testbis<-deriva.gradient(b=b,
  #                       nproc=1,
  #                       funcpa=idmlLikelihoodweiboptim,
  #                       npm=sum(fix==0),
  #                       npar=length(fix),
  #                       bfix=bfix,
  #                       fix=fix,
  #                       ctime=ctime,
  #                       no=no,
  #                       ve01=ve01,
  #                       ve02=ve02,
  #                       ve12=ve12,
  #                       dimnva01=dimnva01,
  #                       dimnva02=dimnva02,
  #                       dimnva12=dimnva12,
  #                       nva01=nva01,
  #                       nva02=nva02,
  #                       nva12=nva12,
  #                       t0=t0,
  #                       t1=t1,
  #                       t2=t2,
  #                       t3=t3,
  #                       troncature=troncature,
  #                       lambda=lambda,
  #                       alpha=alpha,
  #                       penalty.factor=penalty.factor,
  #                       penalty=penalty,
  #                       gausspoint=gausspoint,weib=weib)
  if(start==size_spline){
    ball<-b[1:npm]-b[(1+npm):(2*npm)]
  }else{
    ball<-b[(size_spline-start+1):(npm)]-b[(npm+1):(length(b))]
  }
  
  npm_all<-length(ball)
  grbeta<-rep(0,(npm_all*(npm_all+1)/2)+npm_all)
  
  
  fixbeta<-fix
  fixbeta[1:size_spline]<-1
  bb<-rep(NA,npar)
  bb[which(fix==1)]<-bfix
  bb[which(fixbeta==1 & fix==0)]<-svar
  bb<-na.omit(bb)
  # return first and second derivatives of the loglik
  grbeta<-.Fortran("firstderivaspline",
                   ## input
                   as.double(ball),
                   as.integer(npm_all),
                   as.integer(npar),
                   as.double(bb),
                   as.integer(fixbeta),
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
                   likelihood_deriv=as.double(grbeta),
                   PACKAGE="SmoothHazardoptim9")$likelihood_deriv[1:npm_all]
  

  bb<-rep(NA,npar)
  bb[fix==0]<-c(svar,ball)
  bb[fix==1]<-bfix
  
  if(nva01>0){
    b01<-bb[(size_spline+1):(size_spline+nva01)][penalty.factor[1:nva01]==1]
    
    grbeta[1:(nva01)][penalty.factor[1:nva01]==1]<-grbeta[1:(nva01)][penalty.factor[1:nva01]==1]-2*lambda[,1]*(1-alpha)*b01-sign(b01)*lambda[,1]*(alpha)
    
  }
  
  if(nva02>0){
    b02<-bb[(size_spline+1+nva01):(size_spline+nva01+nva02)][penalty.factor[(nva01+1):(nva01+nva02)]==1]
    
    grbeta[(nva01+1):(nva01+nva02)][penalty.factor[(nva01+1):(nva01+nva02)]==1]<-grbeta[(nva01+1):(nva01+nva02)][penalty.factor[(nva01+1):(nva01+nva02)]==1]-2*lambda[,2]*(1-alpha)*b02-sign(b02)*lambda[,2]*(alpha)
    
    
  }
  
  if(nva12>0){
    b12<-bb[(size_spline+1+nva01+nva02):npar][penalty.factor[(nva01+nva02+1):(nva01+nva02+nva12)]==1]
    
    grbeta[(nva01+nva02+1):(nva01+nva02+nva12)][penalty.factor[(nva01+nva02+1):(nva01+nva02+nva12)]==1]<-grbeta[(nva01+nva02+1):(nva01+nva02+nva12)][penalty.factor[(nva01+nva02+1):(nva01+nva02+nva12)]==1]-2*lambda[,3]*(1-alpha)*b12-sign(b12)*lambda[,3]*(alpha)
    
  }
  
  
  fixbeta<-fix
  fixbeta[(size_spline+1):length(fixbeta)]<-1
  bb<-rep(NA,npar)
  bb[which(fix==1)]<-bfix
  bb[which(fixbeta==1 & fix==0)]<-ball
  bb<-na.omit(bb)
  
  grs<-deriva.gradient(b=svar,
                       nproc=1,
                       funcpa=idmlLikelihood,
                       npm=length(svar),
                       npar=npar,
                       bfix=bb,
                       fix=fixbeta,
                       zi01=knots01,
                       zi02=knots02,
                       zi12=knots12,
                       ctime=ctime,
                       no=no,
                       nz01=nknots01,
                       nz02=nknots02,
                       nz12=nknots12,
                       ve01=ve01,
                       ve02=ve02,
                       ve12=ve12,
                       dimnva01=dimnva01,
                       dimnva02=dimnva02,
                       dimnva12=dimnva12,
                       nva01=nva01,
                       nva02=nva02,
                       nva12=nva12,
                       t0=t0,
                       t1=t1,
                       t2=t2,
                       t3=t3,
                       troncature=troncature,
                       gausspoint=gausspoint)
  
  
  #browser()
  sol<-c(-grs$v,grbeta,-grbeta)
  #<-c(-grs$v,-grbeta,grbeta)
  return(as.double(sol))
}