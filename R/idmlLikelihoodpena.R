idmlLikelihoodpena<-function(b0,np0,npar0,bfix0,fix0,zi010,zi020,zi120,c0,no0,nz010,nz020,nz120,ve010,ve020,ve120,
                         dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                         t00,t10,t20,t30,troncature0,gausspoint0,lambda,alpha,penalty.factor,penalty){
  res<-0
  #browser()
  res<-.Fortran("idmlikelihood",
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
           as.integer(nz010),
           as.integer(nz120),
           as.integer(nz020),
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
  #browser()
  b<-rep(NA,npar0)
  b[fix0==0]<-b0
  b[fix0==1]<-bfix0
  size_spline<-nz010+nz020+nz120+6
  if(nva01>0){
  b01<-b[(size_spline+1):(size_spline+nva01)][penalty.factor[1:nva01]==1]
  }else{b01<-0}
  #b01<-b01[fix0[(size_spline+1):(size_spline+nva01)]==0]
  if(nva02>0){
  b02<-b[(size_spline+1+nva01):(size_spline+nva01+nva02)][penalty.factor[(nva01+1):(nva01+nva02)]==1]
  }else{b02<-0}
  #b02<-b02[fix0[(size_spline+1+nva01):(size_spline+nva01+nva02)]==0]
  if(nva12>0){
  b12<-b[(size_spline+1+nva01+nva02):npar0][penalty.factor[(nva01+nva02+1):(nva01+nva02+nva12)]==1]
  }else{b12<-0}
  #b12<-b12[fix0[(size_spline+1+nva01+nva02):npar0]==0]
  if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
    res<-res+lambda[1]*alpha*sum(abs(b01))+lambda[1]*(1-alpha)*sum(b01*b01)
    res<-res+lambda[2]*alpha*sum(abs(b02))+lambda[2]*(1-alpha)*sum(b02*b02)
    res<-res+lambda[3]*alpha*sum(abs(b12))+lambda[3]*(1-alpha)*sum(b12*b12)
  }
  
  
  if(penalty=="mcp"){
    p01<-unlist(sapply(1:length(b01),FUN=function(x){
      if(b01[x]<=alpha*lambda[1]){
        return(lambda[1]*abs(b01[x])-((b01[x]*b01[x])/2*alpha))
      }else{return(alpha*lambda[1]*lambda[1]/2)}
    }))
    
    p02<-unlist(sapply(1:length(b02),FUN=function(x){
      if(b02[x]<=alpha*lambda[2]){
        return(lambda[2]*abs(b02[x])-((b02[x]*b02[x])/2*alpha))
      }else{return(alpha*lambda[2]*lambda[2]/2)}
    }))
    p12<-unlist(sapply(1:length(b12),FUN=function(x){
      if(b12[x]<=alpha*lambda[3]){
        return(lambda[3]*abs(b12[x])-((b12[x]*b12[x])/2*alpha))
      }else{return(alpha*lambda[3]*lambda[3]/2)}
    }))
    res<-res+sum(p01)+sum(p02)+sum(p12)
    
  }
  
  if(penalty=="scad"){
    
    p01<-unlist(sapply(1:length(b01),FUN=function(x){
      if(b01[x]<=lambda[1]){
        return(lambda[1]*abs(b01[x]))
      }else{
        if(abs(b01[x])<lambda[1]*alpha){
          return((2*alpha*lambda[1]*abs(b01[x])-b01[x]^2-lambda[1]^2)/(2*(alpha-1)))
        }else{
          return((lambda[1]^2)*(alpha+1)/2)
        }
      }
    }))
    
    p02<-unlist(sapply(1:length(b02),FUN=function(x){
      if(b02[x]<=lambda[2]){
        return(lambda[2]*abs(b02[x]))
      }else{
        if(abs(b02[x])<lambda[2]*alpha){
          return((2*alpha*lambda[2]*abs(b02[x])-b02[x]^2-lambda[2]^2)/(2*(alpha-1)))
        }else{
          return((lambda[2]^2)*(alpha+1)/2)
        }
      }
    }))
    
    p12<-unlist(sapply(1:length(b12),FUN=function(x){
      if(b12[x]<=lambda[3]){
        return(lambda[3]*abs(b12[x]))
      }else{
        if(abs(b12[x])<lambda[3]*alpha){
          return((2*alpha*lambda[3]*abs(b12[x])-b12[x]^2-lambda[3]^2)/(2*(alpha-1)))
        }else{
          return((lambda[3]^2)*(alpha+1)/2)
        }
      }
    }))
    res<-res+sum(p01)+sum(p02)+sum(p12)
    
  }
  return(as.double(res))
}


