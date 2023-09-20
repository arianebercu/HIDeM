# 
# 
# # UPDATE BK
# cv.model<-function(beta,
#                    fix,
#                    penalty.factor,
#                    v,
#                    fu,
#                    lambda,
#                    ind.lambda,
#                    alpha,
#                    size_spline){
#   penalty.factor<-c(rep(0,size_spline),penalty.factor)
#   BETA<-previous.beta<-beta[fix==0]
#   NEW.BETA.all<-rep(NA,length(fix))
#   k<-1
#   for(i in 1:length(fix)){
#     #browser()
#     num<-fu[k]-sum(v[k,-k]*BETA[-k])+sum(previous.beta*v[k,])
#     sign<-ifelse(num<0,-1,
#                  ifelse(num>0,1,0))
#     denum<-v[k,k]
#     num<-abs(num)
#     if( i<=size_spline & fix[i]==0){ # unpenalised splines 
#       NEWBETA<-beta[i]
#       k<-k+1
#     }else{
#       if(i<=size_spline & fix[i]==1){
#         NEWBETA<-beta[i]
#         }else{# penalised other parameters
#           if(fix[i]==1 | penalty.factor[i]==0){
#             
#             if(fix[i]==1){
#               NEWBETA<-beta[i]
#             }else{NEWBETA<-sign*num/denum}
#             
#             k<-k+(1-fix[i])
#           }else{
#             if(ind.lambda[i]=="01"){
#               NEWBETA<-ifelse(num>lambda[1]*alpha,sign*(num-lambda[1])/(denum-2*lambda[1]*(1-alpha)),0)
#             }else{
#               if(ind.lambda[i]=="02"){
#                 NEWBETA<-ifelse(num>lambda[2]*alpha,sign*(num-lambda[2])/(denum-2*lambda[2]*(1-alpha)),0)
#               }else{
#                 NEWBETA<-ifelse(num>lambda[3]*alpha,sign*(num-lambda[3])/(denum-2*lambda[3]*(1-alpha)),0)
#               }
#             }
#             k<-k+1
#           }
#         }
#     }
#     NEW.BETA.all[i]<-NEWBETA
#   }
#     
#   
#   return(list(b=as.double(NEW.BETA.all)))
# }
#   
# 




# UPDATE BK
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

  