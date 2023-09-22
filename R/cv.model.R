
cv.model<-function(beta,
                      nva01,
                      nva02,
                      nva12,
                      
                   fix,
                   penalty.factor,
                   penalty,
                   v,
                   fu,
                   lambda,
                   alpha){
  

  BETA<-beta[fix==0]
  NEW.BETA.all<-beta
  penalty.factor<-penalty.factor[fix==0]
  
  num<-sapply(c(1:dim(v)[1]),FUN=function(x){
    fu[x]-sum(v[x,-x]*BETA[-x])+sum(BETA*v[x,])
  })

  
  sign<-ifelse(num<0,-1,
               ifelse(num>0,1,0))
  denum<-diag(v)
  num<-abs(num)
  if(nva01>0){
    num01<-num[1:nva01]
    denum01<-denum[1:nva01]
    sign01<-sign[1:nva01]
  }else{num01<-denum01<-sign01<-NULL}
  
  if(nva02>0){
    num02<-num[(nva01+1):(nva01+nva02)]
    denum02<-denum[(nva01+1):(nva01+nva02)]
    sign02<-sign[(nva01+1):(nva01+nva02)]
    
  }else{num02<-denum02<-sign02<-NULL}
  
  if(nva12>0){
    num12<-num[(nva01+nva02+1):length(num)]
    denum12<-denum[(nva01+nva02+1):length(denum)]
    sign12<-sign[(nva01+nva02+1):length(num)]
  }else{num12<-denum12<-sign12<-NULL}
  
  NEWBETA<-rep(NA,length(num))
  NEWBETA[penalty.factor==0]<-as.double(sign[penalty.factor==0]*num[penalty.factor==0]/denum[penalty.factor==0])
  
  idbeta<-NULL
  if(penalty%in%c("lasso","ridge","elasticnet")){
    # 0 -> 1
    idbeta<-which(num01>lambda[,1]*alpha)
    NEWBETA[idbeta]<-as.double(sign01[idbeta]*(num01[idbeta]-lambda[,1]*alpha)/(denum01[idbeta]+2*lambda[,1]*(1-alpha)))
    idbeta<-which(num01<=lambda[,1]*alpha)
    NEWBETA[idbeta]<-0
    
    # 0 ->2
    idbeta<-which(num02>lambda[,2]*alpha)
    NEWBETA[idbeta+nva01]<-as.double(sign02[idbeta]*(num02[idbeta]-lambda[,2]*alpha)/(denum02[idbeta]+2*lambda[,2]*(1-alpha)))
    idbeta<-which(num02<=lambda[,2]*alpha)
    NEWBETA[idbeta+nva01]<-0
    
    # 1 ->2
    idbeta<-which(num12>lambda[,3]*alpha)
    NEWBETA[idbeta+nva01+nva02]<-as.double(sign12[idbeta]*(num12[idbeta]-lambda[,3]*alpha)/(denum12[idbeta]+2*lambda[,3]*(1-alpha)))
    idbeta<-which(num12<=lambda[,3]*alpha)
    NEWBETA[idbeta+nva01+nva02]<-0
    }
  
  if(penalty%in%c("corrected.elasticnet")){
    
    # 0 ->1
    idbeta<-which(num01>lambda[,1]*alpha)
    NEWBETA[idbeta]<-as.double(sign01[idbeta]*(1+lambda[,1]*(1-alpha))*(num01[idbeta]-lambda[,1]*alpha)/(denum01[idbeta]+2*lambda[,1]*(1-alpha)))
    idbeta<-which(num01<=lambda[,1]*alpha)
    NEWBETA[idbeta]<-0
    
    # 0 ->2
    idbeta<-which(num02>lambda[,2]*alpha)+
    NEWBETA[idbeta+nva01]<-as.double(sign02[idbeta]*(1+lambda[,2]*(1-alpha))*(num02[idbeta]-lambda[,2]*alpha)/(denum02[idbeta]+2*lambda[,2]*(1-alpha)))
    idbeta<-which(num02<=lambda[,2]*alpha)
    NEWBETA[idbeta+nva01]<-0
    
    # 1 ->2
    idbeta<-which(num12>lambda[,3]*alpha)+
    NEWBETA[idbeta+nva01+nva02]<-as.double(sign12[idbeta]*(1+lambda[,3]*(1-alpha))*(num12[idbeta]-lambda[,3]*alpha)/(denum12[idbeta]+2*lambda[,3]*(1-alpha)))
    idbeta<-which(num12<=lambda[,3]*alpha)
    NEWBETA[idbeta+nva01+nva02]<-0
  }
  
  
  if(penalty=="mcp"){
    
    # 0 -> 1
    idbeta<-which(num01>alpha*lambda[,1]*denum01)
    NEWBETA[idbeta]<-as.double(sign01[idbeta]*num01[idbeta]/denum01[idbeta])
    idbeta<-which((num01<=alpha*lambda[,1]*denum01) & num01>lambda[,1])
    NEWBETA[idbeta]<-as.double(sign01[idbeta]*(num01[idbeta]-lambda[,1])/(denum01[idbeta]-(1/alpha)))
    idbeta<-which((num01<=alpha*lambda[,1]*denum01) & num01<=lambda[,1])
    NEWBETA[idbeta]<-0
    
    # 0 -> 2
    idbeta<-which(num02>alpha*lambda[,2]*denum02)
    NEWBETA[idbeta+nva01]<-as.double(sign02[idbeta]*num02[idbeta]/denum02[idbeta])
    idbeta<-which((num02<=alpha*lambda[,2]*denum02) & num02>lambda[,2])
    NEWBETA[idbeta+nva01]<-as.double(sign02[idbeta]*(num02[idbeta]-lambda[,2])/(denum02[idbeta]-(1/alpha)))
    idbeta<-which((num02<=alpha*lambda[,2]*denum02) & num02<=lambda[,2])
    NEWBETA[idbeta+nva01]<-0
    
    
    # 1 -> 2
    idbeta<-which(num12>alpha*lambda[,3]*denum12)
    NEWBETA[idbeta+nva01+nva02]<-as.double(sign12[idbeta]*num12[idbeta]/denum12[idbeta])
    idbeta<-which((num12<=alpha*lambda[,3]*denum12) & num12>lambda[,3])
    NEWBETA[idbeta+nva01+nva02]<-as.double(sign12[idbeta]*(num12[idbeta]-lambda[,3])/(denum12[idbeta]-(1/alpha)))
    idbeta<-which((num12<=alpha*lambda[,3]*denum12) & num12<=lambda[,3])
    NEWBETA[idbeta+nva01+nva02]<-0
    
    }
  
  
  if(penalty=="scad"){
    
    # 0 ->1 
    idbeta<-which(num01>alpha*lambda[,1]*denum01)
    NEWBETA[idbeta]<-as.double(sign01[idbeta]*num01[idbeta]/denum01[idbeta])
    idbeta<-which((num01>alpha*lambda[,1]*denum01) & (num01<=lambda[,1]*(1+denum01)) & (num01>lambda[,1]))
    NEWBETA[idbeta]<-as.double(sign01[idbeta]*(num01[idbeta]-lambda[,1])/denum01[idbeta])
    idbeta<-which((num01>alpha*lambda[,1]*denum01) & (num01<=lambda[,1]*(1+denum01)) & (num01<=lambda[,1]))
    NEWBETA[idbeta]<-0
    idbeta<-which((num01>alpha*lambda[,1]*denum01) & (num01>lambda[,1]*(1+denum01)) & (num01>(lambda[,1]*alpha/(1-alpha))))
    NEWBETA[idbeta]<-as.double(sign01[idbeta]*(num01[idbeta]-lambda[,1]*alpha/(1-alpha))/(denum01[idbeta]-(1/(alpha-1))))
    idbeta<-which((num01>alpha*lambda[,1]*denum01) & (num01>lambda[,1]*(1+denum01)) & (num01<=(lambda[,1]*alpha/(1-alpha))))
    NEWBETA[idbeta]<-0
    
    # 0 ->2
    idbeta<-which(num02>alpha*lambda[,2]*denum02)
    NEWBETA[idbeta+nva01]<-as.double(sign02[idbeta]*num02[idbeta]/denum02[idbeta])
    idbeta<-which((num02>alpha*lambda[,2]*denum02) & (num02<=lambda[,2]*(1+denum02)) & (num02>lambda[,2]))
    NEWBETA[idbeta+nva01]<-as.double(sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta])
    idbeta<-which((num02>alpha*lambda[,2]*denum02) & (num02<=lambda[,2]*(1+denum02)) & (num02<=lambda[,2]))
    NEWBETA[idbeta+nva01]<-0
    idbeta<-which((num02>alpha*lambda[,2]*denum02) & (num02>lambda[,2]*(1+denum02)) & (num02>(lambda[,2]*alpha/(1-alpha))))
    NEWBETA[idbeta+nva01]<-as.double(sign02[idbeta]*(num02[idbeta]-lambda[,2]*alpha/(1-alpha))/(denum02[idbeta]-(1/(alpha-1))))
    idbeta<-which((num02>alpha*lambda[,2]*denum02) & (num02>lambda[,2]*(1+denum02)) & (num02<=(lambda[,2]*alpha/(1-alpha))))
    NEWBETA[idbeta+nva01]<-0
    
    # 1 ->2
    idbeta<-which(num12>alpha*lambda[,3]*denum12)
    NEWBETA[idbeta+nva01+nva02]<-as.double(sign12[idbeta]*num12[idbeta]/denum12[idbeta])
    idbeta<-which((num12>alpha*lambda[,3]*denum12) & (num12<=lambda[,3]*(1+denum12)) & (num12>lambda[,3]))
    NEWBETA[idbeta+nva01+nva02]<-as.double(sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta])
    idbeta<-which((num12>alpha*lambda[,3]*denum12) & (num12<=lambda[,3]*(1+denum12)) & (num12<=lambda[,3]))
    NEWBETA[idbeta+nva01+nva02]<-0
    idbeta<-which((num12>alpha*lambda[,3]*denum12) & (num12>lambda[,3]*(1+denum12)) & (num12>(lambda[,3]*alpha/(1-alpha))))
    NEWBETA[idbeta+nva01+nva02]<-as.double(sign12[idbeta]*(num12[idbeta]-lambda[,3]*alpha/(1-alpha))/(denum12[idbeta]-(1/(alpha-1))))
    idbeta<-which((num12>alpha*lambda[,3]*denum12) & (num12>lambda[,3]*(1+denum12)) & (num12<=(lambda[,3]*alpha/(1-alpha))))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    
    }

  
  NEW.BETA.all[fix==0]<-NEWBETA
  
  return(list(b=as.double(NEW.BETA.all)))
}

  