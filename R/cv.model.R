
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
  NEWBETA<-sapply(c(1:dim(v)[1]),FUN=function(x){

    if(penalty.factor[x]==0){
      return(as.double(sign[x]*num[x]/denum[x]))
    }else{
      if(x<=nva01){
        if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
          if(num[x]>lambda[1]*alpha){
            if(penalty%in%c("lasso","ridge","elasticnet")){
            return(as.double(sign[x]*(num[x]-lambda[1]*alpha)/(denum[x]+2*lambda[1]*(1-alpha))))}
            if(penalty=="corrected.elasticnet"){
              return(as.double(sign[x]*(1+lambda[1]*(1-alpha))*(num[x]-lambda[1]*alpha)/(denum[x]+2*lambda[1]*(1-alpha))))}
            
            
          }else{return(0)}
        }
        
        if(penalty=="mcp"){
          if(num[x]>alpha*lambda[1]*denum[x]){
            return(as.double(sign[x]*num[x]/denum[x]))
          }else{
            if(num[x]>lambda[1]){
              return(as.double(sign[x]*(num[x]-lambda[1])/(denum[x]-(1/alpha))))
            }else{return(0)}
          }
        }
        
        if(penalty=="scad"){
          if(num[x]>alpha*lambda[1]*denum[x]){
            return(as.double(sign[x]*num[x]/denum[x]))
          }else{
          if(num[x]<=lambda[1]*(1+denum[x])){
            if(num[x]>lambda[1]){
              return(as.double(sign[x]*(num[x]-lambda[1])/denum[x]))
            }else{return(0)}
          }else{
            if(num[x]>(lambda[1]*alpha/(1-alpha))){
              return(as.double(sign[x]*(num[x]-lambda[1]*alpha/(1-alpha))/(denum[x]-(1/(alpha-1)))))
            }else{return(0)}
          }
          }
        }
      }else{
        if(x>nva01 & x<=(nva01+nva02)){
          if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
            if(num[x]>lambda[2]*alpha){
              if(penalty%in%c("lasso","ridge","elasticnet")){
              return(as.double(sign[x]*(num[x]-lambda[2]*alpha)/(denum[x]+2*lambda[2]*(1-alpha))))}
              if(penalty=="corrected.elasticnet"){
                return(as.double(sign[x]*(1+lambda[2]*(1-alpha))*(num[x]-lambda[2]*alpha)/(denum[x]+2*lambda[2]*(1-alpha))))}
            }else{return(0)}
          }
          
          if(penalty=="mcp"){
            if(num[x]>alpha*lambda[2]*denum[x]){
              return(as.double(sign[x]*num[x]/denum[x]))
            }else{
              if(num[x]>lambda[2]){
                return(as.double(sign[x]*(num[x]-lambda[2])/(denum[x]-(1/alpha))))
              }else{return(0)}
            }
          }
          
          if(penalty=="scad"){
            if(num[x]>alpha*lambda[2]*denum[x]){
              return(as.double(sign[x]*num[x]/denum[x]))
            }else{
              if(num[x]<=lambda[2]*(1+denum[x])){
                if(num[x]>lambda[2]){
                  return(as.double(sign[x]*(num[x]-lambda[2])/denum[x]))
                }else{return(0)}
              }else{
                if(num[x]>(lambda[2]*alpha/(1-alpha))){
                  return(as.double(sign[x]*(num[x]-lambda[2]*alpha/(1-alpha))/(denum[x]-(1/(alpha-1)))))
                }else{return(0)}
              }
            }
          }
        }else{
          
          if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
            if(num[x]>lambda[3]*alpha){
              if(penalty%in%c("lasso","ridge","elasticnet")){
              return(as.double(sign[x]*(num[x]-lambda[3]*alpha)/(denum[x]+2*lambda[3]*(1-alpha))))}
              if(penalty=="corrected.elasticnet"){
                return(as.double(sign[x]*(1+lambda[3]*(1-alpha))*(num[x]-lambda[3]*alpha)/(denum[x]+2*lambda[3]*(1-alpha))))}
            }else{return(0)}
          }
          
          
          if(penalty=="mcp"){
            if(num[x]>alpha*lambda[3]*denum[x]){
              return(as.double(sign[x]*num[x]/denum[x]))
            }else{
              if(num[x]>lambda[3]){
                return(as.double(sign[x]*(num[x]-lambda[3])/(denum[x]-(1/alpha))))
              }else{return(0)}
            }
          }
          
          if(penalty=="scad"){
            if(num[x]>alpha*lambda[3]*denum[x]){
              return(as.double(sign[x]*num[x]/denum[x]))
            }else{
              if(num[x]<=lambda[3]*(1+denum[x])){
                if(num[x]>lambda[3]){
                  return(as.double(sign[x]*(num[x]-lambda[3])/denum[x]))
                }else{return(0)}
              }else{
                if(num[x]>(lambda[3]*alpha/(1-alpha))){
                  return(as.double(sign[x]*(num[x]-lambda[3]*alpha/(1-alpha))/(denum[x]-(1/(alpha-1)))))
                }else{return(0)}
              }
            }
          }
        }
      }
    }
  }
  )
  
  NEW.BETA.all[fix==0]<-NEWBETA
  
  return(list(b=as.double(NEW.BETA.all)))
}

  