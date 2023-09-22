idmlLikelihoodweibpena<-function(b,np0,npar0,bfix0,fix0,c0,no0,ve010,ve020,ve120,
                         dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                         t00,t10,t20,t30,troncature0,lambda,alpha,penalty.factor,penalty){
  res<-0
  b0<-b
  #browser()
  res<-.Fortran("idmlikelihoodweib",
           ## input
           as.double(b),
           as.integer(np0),
           as.integer(npar0),
           as.double(bfix0),
           as.integer(fix0),
           as.integer(c0),
           as.integer(no0),
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
           likelihood_res=as.double(res),
           PACKAGE="SmoothHazardoptim9")$likelihood_res
  #browser()
  b<-rep(NA,npar0)
  b[fix0==0]<-b0
  b[fix0==1]<-bfix0
  if(nva01>0){
  b01<-b[(6+1):(6+nva01)][penalty.factor[1:nva01]==1]
  }else{b01<-0}
  #b01<-b01[fix0[(6+1):(6+nva01)]==0]
  if(nva02>0){
  b02<-b[(6+1+nva01):(6+nva01+nva02)][penalty.factor[(nva01+1):(nva01+nva02)]==1]
  }else{b02<-0}
  #b02<-b02[fix0[(6+1+nva01):(6+nva01+nva02)]==0]
  if(nva12>0){
  b12<-b[(6+1+nva01+nva02):npar0][penalty.factor[(nva01+nva02+1):(nva01+nva02+nva12)]==1]
  }else{b12<-0}
  #b12<-b12[fix0[(6+1+nva01+nva02):npar0]==0]
  if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
    res<-res+lambda[,1]*alpha*sum(abs(b01))+lambda[,1]*(1-alpha)*sum(b01*b01)
    res<-res+lambda[,2]*alpha*sum(abs(b02))+lambda[,2]*(1-alpha)*sum(b02*b02)
    res<-res+lambda[,3]*alpha*sum(abs(b12))+lambda[,3]*(1-alpha)*sum(b12*b12)
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

    
    res<-res+sum(p01)+sum(p02)+sum(p12)
    
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
   
    
    res<-res+sum(p01)+sum(p02)+sum(p12)
  }
  return(as.double(res))
}

