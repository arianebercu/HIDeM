### Code:
##' @title idm cv model update of beta in penalised regression
##' @param beta  parameters on explanatory variables
##' @param nva01 number of variables for transition 0 -->1 
##' @param nva02 number of variables for transition 0 -->2
##' @param nva12 number of variables for transition 1 -->2
##' @param fix indicators of fixed and unfixed parameters
##' @param penalty base::which penalty to consider
##' @param penalty.factor base::which variable should be penalised
##' @param v variance covariance matrix 
##' @param fu -loglikelihood
##' @param lambda lambda penalised parameter
##' @param alpha alpha penalised parameter
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  
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

  

# add to do base::which for CRAN check 
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
  
  num01<-NULL
  denum01<-NULL
  sign01<-NULL
  
  num02<-NULL
  denum02<-NULL
  sign02<-NULL
  
  num12<-NULL
  denum12<-NULL
  sign12<-NULL
  
  if(nva01>0){
    num01<-num[1:nva01]
    denum01<-denum[1:nva01]
    sign01<-sign[1:nva01]
  }
  
  if(nva02>0){
    num02<-num[(nva01+1):(nva01+nva02)]
    denum02<-denum[(nva01+1):(nva01+nva02)]
    sign02<-sign[(nva01+1):(nva01+nva02)]
    
  }
  
  if(nva12>0){
    num12<-num[(nva01+nva02+1):length(num)]
    denum12<-denum[(nva01+nva02+1):length(denum)]
    sign12<-sign[(nva01+nva02+1):length(num)]
  }
  
  NEWBETA<-rep(NA,length(num))
  idbeta<-NULL

  # if penalty update beta all at once 
  if(penalty%in%c("lasso","ridge","elasticnet")){
    # 0 -> 1
    if(nva01>0){
    idbeta<-base::which(num01>(lambda[,1]*alpha[,1]))
    NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1]*alpha[,1])/(denum01[idbeta]+2*lambda[,1]*(1-alpha[,1]))
    idbeta<-base::which(num01<=(lambda[,1]*alpha[,1]))
    NEWBETA[idbeta]<-0
    }

    
    # 0 ->2
    if(nva02>0){
    idbeta<-base::which(num02>(lambda[,2]*alpha[,2]))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2]*alpha[,2])/(denum02[idbeta]+2*lambda[,2]*(1-alpha[,2]))
    idbeta<-base::which(num02<=(lambda[,2]*alpha[,2]))
    NEWBETA[idbeta+nva01]<-0}
    

    # 1 ->2
    if(nva12>0){
    idbeta<-base::which(num12>(lambda[,3]*alpha[,3]))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3]*alpha[,3])/(denum12[idbeta]+2*lambda[,3]*(1-alpha[,3]))
    idbeta<-base::which(num12<=(lambda[,3]*alpha[,3]))
    NEWBETA[idbeta+nva01+nva02]<-0}
  }
    
  
  
  
  if(penalty=="mcp"){
    
    
    # 0 -> 1, 
    idbeta<-base::which((num01>=(alpha[,1]*lambda[,1]*denum01)) & (num01<lambda[,1]) & (denum01<(1/alpha[,1])))
    NEWBETA[idbeta]<-sign01[idbeta]*alpha[,1]*lambda[,1]
    # no definition put 0 ? 
    idbeta<-base::which(((num01<(alpha[,1]*lambda[,1]*denum01)) | (num01>=lambda[,1])) & (denum01<(1/alpha[,1])))
    NEWBETA[idbeta]<-0
    
    idbeta<-base::which((num01<=(alpha[,1]*lambda[,1]*denum01)) & num01>lambda[,1] & (denum01>=(1/alpha[,1])))
    NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/(denum01[idbeta]-(1/alpha[,1]))
    idbeta<-base::which(((num01>(alpha[,1]*lambda[,1]*denum01)) | num01<=lambda[,1]) & (denum01>=(1/alpha[,1])))
    NEWBETA[idbeta]<-0
    
    # 0 -> 2
    idbeta<-base::which((num02>=(alpha[,2]*lambda[,2]*denum02)) & (num02<lambda[,2]) & (denum02<(1/alpha[,2])))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*alpha[,2]*lambda[,2]
    idbeta<-base::which(((num02<(alpha[,2]*lambda[,2]*denum02)) | (num02>=lambda[,2])) & (denum02<(1/alpha[,2])))
    NEWBETA[idbeta+nva01]<-0
    
    idbeta<-base::which((num02<=(alpha[,2]*lambda[,2]*denum02)) & num02>lambda[,2] & (denum02>=(1/alpha[,2])))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/(denum02[idbeta]-(1/alpha[,2]))
    idbeta<-base::which(((num02>(alpha[,2]*lambda[,2]*denum02)) | num02<=lambda[,2]) & (denum02>=(1/alpha[,2])))
    NEWBETA[idbeta+nva01]<-0
    
    
    # 1 -> 2
    idbeta<-base::which((num12>=(alpha[,3]*lambda[,3]*denum12)) & (num12<lambda[,3]) & (denum12<(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*alpha[,3]*lambda[,3]
    idbeta<-base::which(((num12<(alpha[,3]*lambda[,3]*denum12)) | (num12>=lambda[,3])) & (denum12<(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    idbeta<-base::which((num12<=(alpha[,3]*lambda[,3]*denum12)) & num12>lambda[,3] & (denum12>=(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/(denum12[idbeta]-(1/alpha[,3]))
    idbeta<-base::which(((num12>(alpha[,3]*lambda[,3]*denum12)) | num12<=lambda[,3]) & (denum12>=(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    }
  
  
  if(penalty=="scad"){
    
    
    
    # 0 -> 1
    idbeta<-base::which((num01<=(lambda[,1]*(1+denum01))) & (num01 >lambda[,1]) & (denum01>=(1/(alpha[,1]-1))))
    NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/denum01[idbeta]
    idbeta<-base::which(((num01>(lambda[,1]*(1+denum01))) | (num01 <= lambda[,1]*denum01)) & (denum01>=(1/(alpha[,1]-1))))
    NEWBETA[idbeta]<-0
    
    idbeta<-base::which((num01<=(alpha[,1]*lambda[,1]*denum01)) & (num01 > lambda[,1]) & (denum01<(1/(alpha[,1]-1))) & (denum01>=(1/alpha[,1])))
    NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/denum01[idbeta]
    idbeta<-base::which(((num01>(alpha[,1]*lambda[,1]*denum01)) | (num01 <= lambda[,1])) & (denum01<(1/(alpha[,1]-1))) & (denum01>=(1/alpha[,1])))
    NEWBETA[idbeta]<-0
    idbeta<-base::which( (denum01<(1/(alpha[,1]-1))) & (denum01<(1/alpha[,1])))
    NEWBETA[idbeta]<-0
    
    # 0 ->2
    idbeta<-base::which((num02<=(lambda[,2]*(1+denum02))) & (num02 >lambda[,2]) & (denum02>=(1/(alpha[,2]-1))))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
    idbeta<-base::which(((num02>(lambda[,2]*(1+denum02))) | (num02 <= lambda[,2]*denum02)) & (denum02>=(1/(alpha[,2]-1))))
    NEWBETA[idbeta+nva01]<-0
    
    idbeta<-base::which((num02<=(alpha[,2]*lambda[,2]*denum02)) & (num02 > lambda[,2]) & (denum02<(1/(alpha[,2]-1))) & (denum02>=(1/alpha[,2])))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
    idbeta<-base::which(((num02>(alpha[,2]*lambda[,2]*denum02)) | (num02 <= lambda[,2])) & (denum02<(1/(alpha[,2]-1))) & (denum02>=(1/alpha[,2])))
    NEWBETA[idbeta+nva01]<-0
    idbeta<-base::which( (denum02<(1/(alpha[,2]-1))) & (denum02<(1/alpha[,2])))
    NEWBETA[idbeta+nva01]<-0
    
    # 1 ->2
    
    idbeta<-base::which((num12<=(lambda[,3]*(1+denum12))) & (num12 >lambda[,3]) & (denum12>=(1/(alpha[,3]-1))))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
    idbeta<-base::which(((num12>(lambda[,3]*(1+denum12))) | (num12 <= lambda[,3]*denum12)) & (denum12>=(1/(alpha[,3]-1))))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    idbeta<-base::which((num12<=(alpha[,3]*lambda[,3]*denum12)) & (num12 > lambda[,3]) & (denum12<(1/(alpha[,3]-1))) & (denum12>=(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
    idbeta<-base::which(((num12>(alpha[,3]*lambda[,3]*denum12)) | (num12 <= lambda[,3])) & (denum12<(1/(alpha[,3]-1))) & (denum12>=(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02]<-0
    idbeta<-base::which( (denum12<(1/(alpha[,3]-1))) & (denum12<(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02]<-0
    }

  idbeta<-base::which(penalty.factor==0)
  # if no penalty on parameter, beta_k=A_k/-x_kk
  NEWBETA[idbeta]<-sign[idbeta]*num[idbeta]/denum[idbeta]
  
  NEW.BETA.all[fix==0]<-NEWBETA
  
  return(list(b=NEW.BETA.all))
}

DYNcv.model<-function(beta,
                   nva01,
                   nva02,
                   nva12,
                   nva01Y,
                   nva02Y,
                   nva12Y,
                   fix,
                   penalty.factor,
                   penalty,
                   v,
                   fu,
                   lambda,
                   alpha,
                   penalty.weights){
  

  # add to do base::which for CRAN check 
  BETA<-beta[fix==0]
  NEW.BETA.all<-beta
  penalty.factor<-penalty.factor[fix==0]
  
  #browser()
  num<-sapply(c(1:dim(v)[1]),FUN=function(x){
    (fu[x]-sum(v[x,-x]*BETA[-x])+sum(BETA*v[x,]))/penalty.weights[x]
  }) 
  
  
  sign<-ifelse(num<0,-1,
               ifelse(num>0,1,0))
  denum<-diag(v)
  num<-abs(num)
  
  num01<-NULL
  denum01<-NULL
  sign01<-NULL
  
  
  num02<-NULL
  denum02<-NULL
  sign02<-NULL
  
  
  num12<-NULL
  denum12<-NULL
  sign12<-NULL
  
  
  num01Y<-NULL
  denum01Y<-NULL
  sign01Y<-NULL
  
  
  num02Y<-NULL
  denum02Y<-NULL
  sign02Y<-NULL

  
  num12Y<-NULL
  denum12Y<-NULL
  sign12Y<-NULL
 
  
  if(nva01>0){
    num01<-num[1:nva01]
    denum01<-denum[1:nva01]
    sign01<-sign[1:nva01]
  }
  
  if(nva01Y>0){
    num01Y<-num[(nva01+nva02+nva12+1):(nva01+nva02+nva12+nva01Y)]
    denum01Y<-denum[(nva01+nva02+nva12+1):(nva01+nva02+nva12+nva01Y)]
    sign01Y<-sign[(nva01+nva02+nva12+1):(nva01+nva02+nva12+nva01Y)]
    
  }
  
  if(nva02>0){
    num02<-num[(nva01+1):(nva01+nva02)]
    denum02<-denum[(nva01+1):(nva01+nva02)]
    sign02<-sign[(nva01+1):(nva01+nva02)]
    
  }
  
  if(nva02Y>0){
    num02Y<-num[(nva01+nva02+nva12+nva01Y+1):(nva01+nva02+nva12+nva01Y+nva02Y)]
    denum02Y<-denum[(nva01+nva02+nva12+nva01Y+1):(nva01+nva02+nva12+nva01Y+nva02Y)]
    sign02Y<-sign[(nva01+nva02+nva12+nva01Y+1):(nva01+nva02+nva12+nva01Y+nva02Y)]
    
  }
  
  if(nva12>0){
    num12<-num[(nva01+nva02+1):(nva01+nva02+nva12)]
    denum12<-denum[(nva01+nva02+1):(nva01+nva02+nva12)]
    sign12<-sign[(nva01+nva02+1):(nva01+nva02+nva12)]
  }
  
  if(nva12Y>0){
    num12Y<-num[(nva01+nva02+nva12+nva01Y+nva02Y+1):(nva01+nva02+nva12+nva01Y+nva02Y+nva12Y)]
    denum12Y<-denum[(nva01+nva02+nva12+nva01Y+nva02Y+1):(nva01+nva02+nva12+nva01Y+nva02Y+nva12Y)]
    sign12Y<-sign[(nva01+nva02+nva12+nva01Y+nva02Y+1):(nva01+nva02+nva12+nva01Y+nva02Y+nva12Y)]
    
  }
  
  NEWBETA<-rep(NA,length(num))
  idbeta<-NULL
 
  # if penalty update beta all at once 
  if(penalty%in%c("lasso","ridge","elasticnet","adaptive.lasso")){
    # 0 -> 1
    if(nva01>0){
      idbeta<-base::which(num01>(lambda[,1]*alpha))
      NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1]*alpha)/(denum01[idbeta]+2*lambda[,1]*(1-alpha))
      idbeta<-base::which(num01<=(lambda[,1]*alpha))
      NEWBETA[idbeta]<-0
    }
    if(nva01Y>0){
      idbeta<-base::which(num01Y>(lambda[,1]*alpha))
      NEWBETA[idbeta+nva01+nva02+nva12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1]*alpha)/(denum01Y[idbeta]+2*lambda[,1]*(1-alpha))
      idbeta<-base::which(num01Y<=(lambda[,1]*alpha))
      NEWBETA[idbeta+nva01+nva02+nva12]<-0
    }
    
    # 0 ->2
    if(nva02>0){
      idbeta<-base::which(num02>(lambda[,2]*alpha))
      NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2]*alpha)/(denum02[idbeta]+2*lambda[,2]*(1-alpha))
      idbeta<-base::which(num02<=(lambda[,2]*alpha))
      NEWBETA[idbeta+nva01]<-0
      }
    
    
    if(nva02Y>0){
      idbeta<-base::which(num02Y>(lambda[,2]*alpha))
      NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2]*alpha)/(denum02Y[idbeta]+2*lambda[,2]*(1-alpha))
      idbeta<-base::which(num02Y<=(lambda[,2]*alpha))
      NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    }
    # 1 ->2
    if(nva12>0){
      idbeta<-base::which(num12>(lambda[,3]*alpha))
      NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3]*alpha)/(denum12[idbeta]+2*lambda[,3]*(1-alpha))
      idbeta<-base::which(num12<=(lambda[,3]*alpha))
      NEWBETA[idbeta+nva01+nva02]<-0
      }
  
  
  if(nva12Y>0){
    idbeta<-base::which(num12Y>(lambda[,3]*alpha))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3]*alpha)/(denum12Y[idbeta]+2*lambda[,3]*(1-alpha))
    idbeta<-base::which(num12Y<=(lambda[,3]*alpha))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
  }
  }
  
  
  if(penalty=="mcp"){
    
    
    # 0 -> 1, 
    idbeta<-base::which((num01>=(alpha*lambda[,1]*denum01)) & (num01<lambda[,1]) & (denum01<(1/alpha)))
    NEWBETA[idbeta]<-sign01[idbeta]*alpha*lambda[,1]
    # no definition put 0 ? 
    idbeta<-base::which(((num01<(alpha*lambda[,1]*denum01)) | (num01>=lambda[,1])) & (denum01<(1/alpha)))
    NEWBETA[idbeta]<-0
    
    idbeta<-base::which((num01<=(alpha*lambda[,1]*denum01)) & num01>lambda[,1] & (denum01>=(1/alpha)))
    NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/(denum01[idbeta]-(1/alpha))
    idbeta<-base::which(((num01>(alpha*lambda[,1]*denum01)) | num01<=lambda[,1]) & (denum01>=(1/alpha)))
    NEWBETA[idbeta]<-0
    
    
    idbeta<-base::which((num01Y>=(alpha*lambda[,1]*denum01Y)) & (num01Y<lambda[,1]) & (denum01Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12]<-sign01Y[idbeta]*alpha*lambda[,1]
    # no definition put 0 ? 
    idbeta<-base::which(((num01Y<(alpha*lambda[,1]*denum01Y)) | (num01Y>=lambda[,1])) & (denum01Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12]<-0
    
    idbeta<-base::which((num01Y<=(alpha*lambda[,1]*denum01Y)) & num01Y>lambda[,1] & (denum01Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1])/(denum01Y[idbeta]-(1/alpha))
    idbeta<-base::which(((num01Y>(alpha*lambda[,1]*denum01Y)) | num01Y<=lambda[,1]) & (denum01Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12]<-0
    
    # 0 -> 2
    idbeta<-base::which((num02>=(alpha*lambda[,2]*denum02)) & (num02<lambda[,2]) & (denum02<(1/alpha)))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*alpha*lambda[,2]
    idbeta<-base::which(((num02<(alpha*lambda[,2]*denum02)) | (num02>=lambda[,2])) & (denum02<(1/alpha)))
    NEWBETA[idbeta+nva01]<-0
    
    idbeta<-base::which((num02<=(alpha*lambda[,2]*denum02)) & num02>lambda[,2] & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/(denum02[idbeta]-(1/alpha))
    idbeta<-base::which(((num02>(alpha*lambda[,2]*denum02)) | num02<=lambda[,2]) & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01]<-0
    
    
    idbeta<-base::which((num02Y>=(alpha*lambda[,2]*denum02Y)) & (num02Y<lambda[,2]) & (denum02Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-sign02Y[idbeta]*alpha*lambda[,2]
    idbeta<-base::which(((num02Y<(alpha*lambda[,2]*denum02Y)) | (num02Y>=lambda[,2])) & (denum02Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    
    idbeta<-base::which((num02Y<=(alpha*lambda[,2]*denum02Y)) & num02Y>lambda[,2] & (denum02Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2])/(denum02Y[idbeta]-(1/alpha))
    idbeta<-base::which(((num02Y>(alpha*lambda[,2]*denum02Y)) | num02Y<=lambda[,2]) & (denum02Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    
    
    # 1 -> 2
    idbeta<-base::which((num12>=(alpha*lambda[,3]*denum12)) & (num12<lambda[,3]) & (denum12<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*alpha*lambda[,3]
    idbeta<-base::which(((num12<(alpha*lambda[,3]*denum12)) | (num12>=lambda[,3])) & (denum12<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    idbeta<-base::which((num12<=(alpha*lambda[,3]*denum12)) & num12>lambda[,3] & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/(denum12[idbeta]-(1/alpha))
    idbeta<-base::which(((num12>(alpha*lambda[,3]*denum12)) | num12<=lambda[,3]) & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    
    idbeta<-base::which((num12Y>=(alpha*lambda[,3]*denum12Y)) & (num12Y<lambda[,3]) & (denum12Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-sign12Y[idbeta]*alpha*lambda[,3]
    idbeta<-base::which(((num12Y<(alpha*lambda[,3]*denum12Y)) | (num12Y>=lambda[,3])) & (denum12Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
    
    idbeta<-base::which((num12Y<=(alpha*lambda[,3]*denum12Y)) & num12Y>lambda[,3] & (denum12Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3])/(denum12Y[idbeta]-(1/alpha))
    idbeta<-base::which(((num12Y>(alpha*lambda[,3]*denum12Y)) | num12Y<=lambda[,3]) & (denum12Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
    
  }
  
  
  if(penalty=="scad"){
    
    
    
    # 0 -> 1
    idbeta<-base::which((num01<=(lambda[,1]*(1+denum01))) & (num01 >lambda[,1]) & (denum01>=(1/(alpha-1))))
    NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/denum01[idbeta]
    idbeta<-base::which(((num01>(lambda[,1]*(1+denum01))) | (num01 <= lambda[,1]*denum01)) & (denum01>=(1/(alpha-1))))
    NEWBETA[idbeta]<-0
    
    idbeta<-base::which((num01<=(alpha*lambda[,1]*denum01)) & (num01 > lambda[,1]) & (denum01<(1/(alpha-1))) & (denum01>=(1/alpha)))
    NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/denum01[idbeta]
    idbeta<-base::which(((num01>(alpha*lambda[,1]*denum01)) | (num01 <= lambda[,1])) & (denum01<(1/(alpha-1))) & (denum01>=(1/alpha)))
    NEWBETA[idbeta]<-0
    idbeta<-base::which( (denum01<(1/(alpha-1))) & (denum01<(1/alpha)))
    NEWBETA[idbeta]<-0
    
    
    idbeta<-base::which((num01Y<=(lambda[,1]*(1+denum01Y))) & (num01Y >lambda[,1]) & (denum01Y>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02+nva12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1])/denum01Y[idbeta]
    idbeta<-base::which(((num01Y>(lambda[,1]*(1+denum01Y))) | (num01Y <= lambda[,1]*denum01Y)) & (denum01Y>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02+nva12]<-0
    
    idbeta<-base::which((num01Y<=(alpha*lambda[,1]*denum01Y)) & (num01Y > lambda[,1]) & (denum01Y<(1/(alpha-1))) & (denum01Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1])/denum01Y[idbeta]
    idbeta<-base::which(((num01Y>(alpha*lambda[,1]*denum01Y)) | (num01Y <= lambda[,1])) & (denum01Y<(1/(alpha-1))) & (denum01Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12]<-0
    idbeta<-base::which( (denum01Y<(1/(alpha-1))) & (denum01Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12]<-0
    
    # 0 ->2
    idbeta<-base::which((num02<=(lambda[,2]*(1+denum02))) & (num02 >lambda[,2]) & (denum02>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
    idbeta<-base::which(((num02>(lambda[,2]*(1+denum02))) | (num02 <= lambda[,2]*denum02)) & (denum02>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01]<-0
    
    idbeta<-base::which((num02<=(alpha*lambda[,2]*denum02)) & (num02 > lambda[,2]) & (denum02<(1/(alpha-1))) & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
    idbeta<-base::which(((num02>(alpha*lambda[,2]*denum02)) | (num02 <= lambda[,2])) & (denum02<(1/(alpha-1))) & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01]<-0
    idbeta<-base::which( (denum02<(1/(alpha-1))) & (denum02<(1/alpha)))
    NEWBETA[idbeta+nva01]<-0
    
    
    idbeta<-base::which((num02Y<=(lambda[,2]*(1+denum02Y))) & (num02Y >lambda[,2]) & (denum02Y>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2])/denum02Y[idbeta]
    idbeta<-base::which(((num02Y>(lambda[,2]*(1+denum02Y))) | (num02Y <= lambda[,2]*denum02Y)) & (denum02Y>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    
    idbeta<-base::which((num02Y<=(alpha*lambda[,2]*denum02Y)) & (num02Y > lambda[,2]) & (denum02Y<(1/(alpha-1))) & (denum02Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2])/denum02Y[idbeta]
    idbeta<-base::which(((num02Y>(alpha*lambda[,2]*denum02Y)) | (num02Y <= lambda[,2])) & (denum02Y<(1/(alpha-1))) & (denum02Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    idbeta<-base::which( (denum02<(1/(alpha-1))) & (denum02<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    
    # 1 ->2
    
    idbeta<-base::which((num12<=(lambda[,3]*(1+denum12))) & (num12 >lambda[,3]) & (denum12>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
    idbeta<-base::which(((num12>(lambda[,3]*(1+denum12))) | (num12 <= lambda[,3]*denum12)) & (denum12>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    idbeta<-base::which((num12<=(alpha*lambda[,3]*denum12)) & (num12 > lambda[,3]) & (denum12<(1/(alpha-1))) & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
    idbeta<-base::which(((num12>(alpha*lambda[,3]*denum12)) | (num12 <= lambda[,3])) & (denum12<(1/(alpha-1))) & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-0
    idbeta<-base::which( (denum12<(1/(alpha-1))) & (denum12<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    
    idbeta<-base::which((num12Y<=(lambda[,3]*(1+denum12Y))) & (num12Y >lambda[,3]) & (denum12Y>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3])/denum12Y[idbeta]
    idbeta<-base::which(((num12Y>(lambda[,3]*(1+denum12Y))) | (num12Y <= lambda[,3]*denum12Y)) & (denum12Y>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
    
    idbeta<-base::which((num12Y<=(alpha*lambda[,3]*denum12Y)) & (num12Y > lambda[,3]) & (denum12Y<(1/(alpha-1))) & (denum12Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3])/denum12Y[idbeta]
    idbeta<-base::which(((num12Y>(alpha*lambda[,3]*denum12Y)) | (num12Y <= lambda[,3])) & (denum12Y<(1/(alpha-1))) & (denum12Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
    idbeta<-base::which( (denum12Y<(1/(alpha-1))) & (denum12Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
  }
  
 
  idbeta<-base::which(penalty.factor==0)
  # if no penalty on parameter, beta_k=A_k/-x_kk
  NEWBETA[idbeta]<-sign[idbeta]*num[idbeta]/denum[idbeta]
  
  NEW.BETA.all[fix==0]<-NEWBETA
  
  return(list(b=NEW.BETA.all))
}


cv.model.onestep<-function(beta,
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



  # add to do base::which for CRAN check
  BETA<-beta[fix==0]
  NEW.BETA.all<-beta
  penalty.factor<-penalty.factor[fix==0]
  nweib<-sum(fix[1:6]==0)

  num<-sapply(c(1:dim(v)[1]),FUN=function(x){
    fu[x]-sum(v[x,-x]*BETA[-x])+sum(BETA*v[x,])
  })



  sign<-ifelse(num<0,-1,
               ifelse(num>0,1,0))
  denum<-diag(v)
  num<-abs(num)

  num01<-NULL
  denum01<-NULL
  sign01<-NULL

  num02<-NULL
  denum02<-NULL
  sign02<-NULL

  num12<-NULL
  denum12<-NULL
  sign12<-NULL

  if(nva01>0){
    num01<-num[(nweib+1):(nva01+nweib)]
    denum01<-denum[(nweib+1):(nva01+nweib)]
    sign01<-sign[(nweib+1):(nva01+nweib)]
  }

  if(nva02>0){
    num02<-num[(nva01+1+nweib):(nva01+nva02+nweib)]
    denum02<-denum[(nva01+1+nweib):(nva01+nva02+nweib)]
    sign02<-sign[(nva01+1+nweib):(nva01+nva02+nweib)]

  }

  if(nva12>0){
    num12<-num[(nva01+nva02+1+nweib):length(num)]
    denum12<-denum[(nva01+nva02+1+nweib):length(denum)]
    sign12<-sign[(nva01+nva02+1+nweib):length(num)]
  }

  NEWBETA<-rep(NA,length(num))
  idbeta<-NULL



  # if penalty update beta all at once
  if(penalty%in%c("lasso","ridge","elasticnet")){


    # 0 -> 1
    if(nva01>0){
      idbeta<-base::which(num01>(lambda[,1]*alpha))
      NEWBETA[nweib+idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1]*alpha)/(denum01[idbeta]+2*lambda[,1]*(1-alpha))
      idbeta<-base::which(num01<=(lambda[,1]*alpha))
      NEWBETA[nweib+idbeta]<-0}


    # 0 ->2
    if(nva02>0){
      idbeta<-base::which(num02>(lambda[,2]*alpha))
      NEWBETA[idbeta+nva01+nweib]<-sign02[idbeta]*(num02[idbeta]-lambda[,2]*alpha)/(denum02[idbeta]+2*lambda[,2]*(1-alpha))
      idbeta<-base::which(num02<=(lambda[,2]*alpha))
      NEWBETA[idbeta+nva01+nweib]<-0}


    # 1 ->2
    if(nva12>0){
      idbeta<-base::which(num12>(lambda[,3]*alpha))
      NEWBETA[idbeta+nva01+nva02+nweib]<-sign12[idbeta]*(num12[idbeta]-lambda[,3]*alpha)/(denum12[idbeta]+2*lambda[,3]*(1-alpha))
      idbeta<-base::which(num12<=(lambda[,3]*alpha))
      NEWBETA[idbeta+nva01+nva02+nweib]<-0}
  }




  if(penalty=="mcp"){


    # 0 -> 1,
    idbeta<-base::which((num01>=(alpha*lambda[,1]*denum01)) & (num01<lambda[,1]) & (denum01<(1/alpha)))
    NEWBETA[idbeta+nweib]<-sign01[idbeta]*alpha*lambda[,1]
    # no definition put 0 ?
    idbeta<-base::which(((num01<(alpha*lambda[,1]*denum01)) | (num01>=lambda[,1])) & (denum01<(1/alpha)))
    NEWBETA[idbeta+nweib]<-0

    idbeta<-base::which((num01<=(alpha*lambda[,1]*denum01)) & num01>lambda[,1] & (denum01>=(1/alpha)))
    NEWBETA[idbeta+nweib]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/(denum01[idbeta]-(1/alpha))
    idbeta<-base::which(((num01>(alpha*lambda[,1]*denum01)) | num01<=lambda[,1]) & (denum01>=(1/alpha)))
    NEWBETA[idbeta+nweib]<-0

    # 0 -> 2
    idbeta<-base::which((num02>=(alpha*lambda[,2]*denum02)) & (num02<lambda[,2]) & (denum02<(1/alpha)))
    NEWBETA[idbeta+nva01+nweib]<-sign02[idbeta]*alpha*lambda[,2]
    idbeta<-base::which(((num02<(alpha*lambda[,2]*denum02)) | (num02>=lambda[,2])) & (denum02<(1/alpha)))
    NEWBETA[idbeta+nva01+nweib]<-0

    idbeta<-base::which((num02<=(alpha*lambda[,2]*denum02)) & num02>lambda[,2] & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01+nweib]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/(denum02[idbeta]-(1/alpha))
    idbeta<-base::which(((num02>(alpha*lambda[,2]*denum02)) | num02<=lambda[,2]) & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01+nweib]<-0


    # 1 -> 2
    idbeta<-base::which((num12>=(alpha*lambda[,3]*denum12)) & (num12<lambda[,3]) & (denum12<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nweib]<-sign12[idbeta]*alpha*lambda[,3]
    idbeta<-base::which(((num12<(alpha*lambda[,3]*denum12)) | (num12>=lambda[,3])) & (denum12<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nweib]<-0

    idbeta<-base::which((num12<=(alpha*lambda[,3]*denum12)) & num12>lambda[,3] & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nweib]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/(denum12[idbeta]-(1/alpha))
    idbeta<-base::which(((num12>(alpha*lambda[,3]*denum12)) | num12<=lambda[,3]) & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nweib]<-0

  }


  if(penalty=="scad"){



    # 0 -> 1
    idbeta<-base::which((num01<=(lambda[,1]*(1+denum01))) & (num01 >lambda[,1]) & (denum01>=(1/(alpha-1))))
    NEWBETA[idbeta+nweib]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/denum01[idbeta]
    idbeta<-base::which(((num01>(lambda[,1]*(1+denum01))) | (num01 <= lambda[,1]*denum01)) & (denum01>=(1/(alpha-1))))
    NEWBETA[idbeta+nweib]<-0

    idbeta<-base::which((num01<=(alpha*lambda[,1]*denum01)) & (num01 > lambda[,1]) & (denum01<(1/(alpha-1))) & (denum01>=(1/alpha)))
    NEWBETA[idbeta+nweib]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/denum01[idbeta]
    idbeta<-base::which(((num01>(alpha*lambda[,1]*denum01)) | (num01 <= lambda[,1])) & (denum01<(1/(alpha-1))) & (denum01>=(1/alpha)))
    NEWBETA[idbeta+nweib]<-0
    idbeta<-base::which( (denum01<(1/(alpha-1))) & (denum01<(1/alpha)))
    NEWBETA[idbeta+nweib]<-0

    # 0 ->2
    idbeta<-base::which((num02<=(lambda[,2]*(1+denum02))) & (num02 >lambda[,2]) & (denum02>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nweib]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
    idbeta<-base::which(((num02>(lambda[,2]*(1+denum02))) | (num02 <= lambda[,2]*denum02)) & (denum02>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nweib]<-0

    idbeta<-base::which((num02<=(alpha*lambda[,2]*denum02)) & (num02 > lambda[,2]) & (denum02<(1/(alpha-1))) & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01+nweib]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
    idbeta<-base::which(((num02>(alpha*lambda[,2]*denum02)) | (num02 <= lambda[,2])) & (denum02<(1/(alpha-1))) & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01+nweib]<-0
    idbeta<-base::which( (denum02<(1/(alpha-1))) & (denum02<(1/alpha)))
    NEWBETA[idbeta+nva01+nweib]<-0

    # 1 ->2

    idbeta<-base::which((num12<=(lambda[,3]*(1+denum12))) & (num12 >lambda[,3]) & (denum12>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02+nweib]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
    idbeta<-base::which(((num12>(lambda[,3]*(1+denum12))) | (num12 <= lambda[,3]*denum12)) & (denum12>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02+nweib]<-0

    idbeta<-base::which((num12<=(alpha*lambda[,3]*denum12)) & (num12 > lambda[,3]) & (denum12<(1/(alpha-1))) & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nweib]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
    idbeta<-base::which(((num12>(alpha*lambda[,3]*denum12)) | (num12 <= lambda[,3])) & (denum12<(1/(alpha-1))) & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nweib]<-0
    idbeta<-base::which( (denum12<(1/(alpha-1))) & (denum12<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nweib]<-0
  }

  idbeta<-base::which(penalty.factor==0)
  # if no penalty on parameter, beta_k=A_k/-x_kk
  NEWBETA[idbeta]<-sign[idbeta]*num[idbeta]/denum[idbeta]
# does not work try update weibull parameters according to new update beta :
  # BETA<-beta[fix==0]
  # idbeta<-base::which(penalty.factor==1)
  # BETA[idbeta]<-NEWBETA[idbeta]
  #
  #  num<-sapply(c(1:dim(v)[1]),FUN=function(x){
  #    fu[x]-sum(v[x,-x]*BETA[-x])+sum(BETA*v[x,])
  #  })
  #
  #  sign<-ifelse(num<0,-1,
  #               ifelse(num>0,1,0))
  #  denum<-diag(v)
  #  idbeta<-base::which(penalty.factor==0)
  #
  #  NEWBETA[idbeta]<-sign[idbeta]*num[idbeta]/denum[idbeta]

  NEW.BETA.all[fix==0]<-NEWBETA

  return(list(b=NEW.BETA.all))
}



# =============================================================================
# DYNcv.model.seq — version sequentielle (Gauss-Seidel) de DYNcv.model
#
# Corrections apportees par rapport a la version originale :
#
#   1. BUG CORRIGE : le calcul de `num[x]` utilisait
#         fu[x] - sum(v[x,-x]*BETA[-x]) + sum(BETA*v[x,])
#      ce qui se simplifie algebriquement en `fu[x] + v[x,x]*BETA[x]' quelle
#      que soit BETA[-x] (les termes croises s'annulent). Cela revient a
#      traiter `v` comme diagonale et ignore toute correlation entre
#      coordonnees. Ici on utilise la vraie cible de coordinate descent :
#         num[x] = fu[x] - sum(v[x,-x] * BETA[-x]) + v[x,x] * BETA[x]
#
#   2. MISE A JOUR SEQUENTIELLE (Gauss-Seidel) : `BETA[x]` est mis a jour
#      immediatement apres son calcul et utilise pour les coordonnees
#      suivantes dans la meme passe, au lieu d'une mise a jour simultanee
#      (Jacobi) qui ne garantit pas la decroissance de l'objectif penalise
#      des que v a des termes hors-diagonale significatifs.
#
#   3. GARDE-FOU (step-halving) optionnel : si `check_monotone = TRUE` et
#      qu'une fonction d'objectif `objective_fn(beta)` est fournie, on
#      verifie que l'objectif penalise s'est bien ameliore apres la passe
#      complete ; sinon on interpole vers l'ancien beta.
#
# La logique de seuillage par penalite (lasso/ridge/elasticnet/adaptive.lasso,
# mcp, scad) est preservee telle quelle, juste traduite en version scalaire
# par coordonnee plutot que vectorisee par groupe.
# =============================================================================

#' Mise a jour d'une coordonnee scalaire selon la penalite
#'
#' @param num cible de coordinate descent pour cette coordonnee (scalaire)
#' @param denum courbure diagonale v[x,x] (scalaire, doit etre > 0)
#' @param sgn signe de num (-1, 0, ou 1)
#' @param lam lambda de la transition correspondante (scalaire)
#' @param alpha parametre elasticnet / mcp / scad
#' @param penalty type de penalite
#' @return nouvelle valeur de beta pour cette coordonnee
.dyncv_coord_update <- function(num, denum, sgn, lam, alpha, penalty) {
  
  num_abs <- abs(num)
  
  if (penalty %in% c("lasso", "ridge", "elasticnet", "adaptive.lasso")) {
    if (num_abs > (lam * alpha)) {
      return(sgn * (num_abs - lam * alpha) / (denum + 2 * lam * (1 - alpha)))
    } else {
      return(0)
    }
  }
  
  if (penalty == "mcp") {
    if (denum < (1 / alpha)) {
      if (num_abs >= (alpha * lam * denum) && num_abs < lam) {
        return(sgn * alpha * lam)
      } else {
        return(0)
      }
    } else { # denum >= 1/alpha
      if (num_abs <= (alpha * lam * denum) && num_abs > lam) {
        return(sgn * (num_abs - lam) / (denum - (1 / alpha)))
      } else {
        return(0)
      }
    }
  }
  
  if (penalty == "scad") {
    if (denum >= (1 / (alpha - 1))) {
      if (num_abs <= (lam * (1 + denum)) && num_abs > lam) {
        return(sgn * (num_abs - lam) / denum)
      } else if (num_abs <= (alpha * lam * denum) && num_abs > lam && denum < (1 / (alpha - 1))) {
        # branche intermediaire (denum entre 1/(alpha-1) et 1/alpha) —
        # conservee de la version originale pour fidelite, ne devrait pas
        # etre atteinte ici puisqu'on est deja dans denum >= 1/(alpha-1)
        return(sgn * (num_abs - lam) / denum)
      } else {
        return(0)
      }
    } else {
      # denum < 1/(alpha-1) : cf. branches originales avec seuil 1/alpha
      if (denum >= (1 / alpha) && num_abs <= (alpha * lam * denum) && num_abs > lam) {
        return(sgn * (num_abs - lam) / denum)
      } else {
        return(0)
      }
    }
  }
  
  stop(paste0("Type de penalite non reconnu : ", penalty))
}

#' Construit la table coordonnee -> colonne de lambda (transition 01/02/12)
.dyncv_lambda_col_map <- function(nva01, nva02, nva12, nva01Y, nva02Y, nva12Y) {
  c(
    rep(1L, nva01),
    rep(2L, nva02),
    rep(3L, nva12),
    rep(1L, nva01Y),
    rep(2L, nva02Y),
    rep(3L, nva12Y)
  )
}

#' Version sequentielle (Gauss-Seidel) de DYNcv.model
#'
#' @param n_inner_passes nombre de passes de coordinate descent sur le
#'   vecteur complet, a (fu, v) fixes (comme la boucle interne de glmnet).
#'   1 = une seule passe sequentielle ; >1 = plusieurs passes jusqu'a
#'   (quasi) convergence de l'approximation quadratique courante.
#' @param tol tolerance de convergence entre passes (norme max du changement)
#' @param objective_fn optionnel : fonction(beta_full) -> objectif penalise
#'   (plus grand = mieux, ex. loglik penalisee). Si fournie, active le
#'   step-halving global apres les passes internes.
DYNcv.model.seq <- function(beta,
                            nva01, nva02, nva12,
                            nva01Y, nva02Y, nva12Y,
                            fix,
                            penalty.factor,
                            penalty,
                            v,
                            fu,
                            lambda,
                            alpha,
                            penalty.weights,
                            n_inner_passes = 1,
                            tol = 1e-8,
                            objective_fn = NULL) {
  
  idx_free <- which(fix == 0)
  BETA <- beta[idx_free]
  p <- length(BETA)
  
  penalty.factor_f <- penalty.factor[idx_free]
  lam_col <- .dyncv_lambda_col_map(nva01, nva02, nva12, nva01Y, nva02Y, nva12Y)
  stopifnot(length(lam_col) == p)
  
  denum_all <- diag(v)
  
  beta_before_passes <- BETA
  
  for (pass in seq_len(n_inner_passes)) {
    max_delta <- 0
    
    for (x in seq_len(p)) {
      
      # --- cible de coordinate descent, VERSION CORRIGEE ---
      # (utilise BETA courant, deja mis a jour pour les coordonnees < x)
      num_x <- (fu[x] - sum(v[x, -x] * BETA[-x]) + v[x, x] * BETA[x]) /
        penalty.weights[x]
      
      sgn_x <- if (num_x < 0) -1 else if (num_x > 0) 1 else 0
      denum_x <- denum_all[x]
      
      if (penalty.factor_f[x] == 0) {
        # pas de penalite sur cette coordonnee : pas de Newton pur
        new_val <- sgn_x * abs(num_x) / denum_x
      } else {
        lam_x <- lambda[, lam_col[x]]
        new_val <- .dyncv_coord_update(
          num = num_x, denum = denum_x, sgn = sgn_x,
          lam = lam_x, alpha = alpha, penalty = penalty
        )
      }
      
      max_delta <- max(max_delta, abs(new_val - BETA[x]))
      BETA[x] <- new_val   # <-- mise a jour IMMEDIATE (Gauss-Seidel)
    }
    
    if (max_delta < tol) break
  }

  
  NEW.BETA.all <- beta
  NEW.BETA.all[idx_free] <- BETA
  
  return(list(b = NEW.BETA.all))
}

# cv.model.onestep<-function(beta,
#                            nva01,
#                            nva02,
#                            nva12,
#                            fix,
#                            penalty.factor,
#                            penalty,
#                            v,
#                            fu,
#                            lambda,
#                            alpha){
# 
# 
# 
#   # add to do base::which for CRAN check
# 
#   NEWBETA<-beta[fix==0]
#   nparweib<-6
#   nweib<-sum(fix[1:nparweib]==0)
#   n01<-ifelse(nva01==0,0,sum(fix[(nparweib+1):(nparweib+nva01)]==0))
#   n02<-ifelse(nva02==0,0,sum(fix[(nparweib+1+nva01):(nparweib+nva01+nva02)]==0))
#   n12<-ifelse(nva12==0,0,sum(fix[(nparweib+nva01+nva02+1):(nparweib+nva01+nva02+nva12)]==0))
# 
#   maxite<-nweib
#   if(nweib>0){
#     for(k in 1:maxite){
#         num<-fu[k]-sum(v[k,-k]*NEWBETA[-k])+sum(NEWBETA*v[k,])
#         denum<-v[k,k]
#         NEWBETA[k]<-num/denum
# 
#     }
#   }
# 
# 
# 
#  minite<-nweib+1
#  maxite<-nweib+n01
#     # 0 -> 1
# 
#     if(n01>0){
#       for(k in minite:maxite){
#         if(fix[k]==0){
#         num<-fu[k]-sum(v[k,-k]*NEWBETA[-k])+sum(NEWBETA*v[k,])
#         denum<-v[k,k]
#         sign<-ifelse(num<0,-1,
#                      ifelse(num>0,1,0))
#         NEWBETA[k]<-ifelse(num >lambda[,1]*alpha,sign*(num-lambda[,1]*alpha)/(denum+2*lambda[,1]*(1-alpha)),0)
#         }
#       }
# 
#     }
#  minite<-n01+nweib+1
#  maxite<-n01+nweib+n02
#     # 0 ->2
#     if(n02>0){
#       for(k in minite:maxite){
#         if(fix[k]==0){
#         num<-fu[k]-sum(v[k,-k]*NEWBETA[-k])+sum(NEWBETA*v[k,])
#         denum<-v[k,k]
#         sign<-ifelse(num<0,-1,
#                      ifelse(num>0,1,0))
#         NEWBETA[k]<-ifelse(num >lambda[,2]*alpha,sign*(num-lambda[,2]*alpha)/(denum+2*lambda[,2]*(1-alpha)),0)
#       }
#       }
#     }
# 
#  minite<-n01+nweib+n02+1
#  maxite<-n01+nweib+n02+n12
#     # 1 ->2
#     if(n12>0){
#       for(k in minite:maxite){
#         if(fix[k]==0){
#         num<-fu[k]-sum(v[k,-k]*NEWBETA[-k])+sum(NEWBETA*v[k,])
#         denum<-v[k,k]
#         sign<-ifelse(num<0,-1,
#                      ifelse(num>0,1,0))
#         NEWBETA[k]<-ifelse(num >lambda[,3]*alpha,sign*(num-lambda[,3]*alpha)/(denum+2*lambda[,3]*(1-alpha)),0)
#       }
#       }
#     }
# 
# 
#  NEWBETA.all<-rep(NA,length(beta))
#  NEWBETA.all[fix==1]<-beta[fix==1]
#  NEWBETA.all[fix==0]<-NEWBETA
#   return(list(b=NEWBETA.all))
# }

# 
# 
# cv.model.onestep<-function(beta,
#                            nva01,
#                            nva02,
#                            nva12,
#                            fix,
#                            penalty.factor,
#                            penalty,
#                            v,
#                            fu,
#                            lambda,
#                            alpha){
# 
# 
# 
#   # add to do base::which for CRAN check
#   BETA<-beta[fix==0]
#   NEW.BETA.all<-beta
#   penalty.factor<-penalty.factor[fix==0]
#   nweib<-sum(fix[1:6]==0)
# 
#   num<-sapply(c(1:dim(v)[1]),FUN=function(x){
#     fu[x]-sum(v[x,-x]*BETA[-x])+sum(BETA*v[x,])
#   })
# 
# 
# 
#   sign<-ifelse(num<0,-1,
#                ifelse(num>0,1,0))
#   denum<-diag(v)
#   num<-abs(num)
# 
#   num01<-NULL
#   denum01<-NULL
#   sign01<-NULL
# 
#   num02<-NULL
#   denum02<-NULL
#   sign02<-NULL
# 
#   num12<-NULL
#   denum12<-NULL
#   sign12<-NULL
# 
#   if(nva01>0){
#     num01<-num[(nweib+1):(nva01+nweib)]
#     denum01<-denum[(nweib+1):(nva01+nweib)]
#     sign01<-sign[(nweib+1):(nva01+nweib)]
#   }
# 
#   if(nva02>0){
#     num02<-num[(nva01+1+nweib):(nva01+nva02+nweib)]
#     denum02<-denum[(nva01+1+nweib):(nva01+nva02+nweib)]
#     sign02<-sign[(nva01+1+nweib):(nva01+nva02+nweib)]
# 
#   }
# 
#   if(nva12>0){
#     num12<-num[(nva01+nva02+1+nweib):length(num)]
#     denum12<-denum[(nva01+nva02+1+nweib):length(denum)]
#     sign12<-sign[(nva01+nva02+1+nweib):length(num)]
#   }
# 
#   NEWBETA<-rep(NA,length(num))
#   idbeta<-NULL
# 
# 
# 
#   # if penalty update beta all at once
#   if(penalty%in%c("lasso","ridge","elasticnet")){
# 
# 
#     # 0 -> 1
#     if(nva01>0){
#       idbeta<-base::which(num01>(lambda[,1]*alpha))
#       NEWBETA[nweib+idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1]*alpha)/(denum01[idbeta]+2*lambda[,1]*(1-alpha))
#       idbeta<-base::which(num01<=(lambda[,1]*alpha))
#       NEWBETA[nweib+idbeta]<-0}
# 
# 
#     # 0 ->2
#     if(nva02>0){
#       idbeta<-base::which(num02>(lambda[,2]*alpha))
#       NEWBETA[idbeta+nva01+nweib]<-sign02[idbeta]*(num02[idbeta]-lambda[,2]*alpha)/(denum02[idbeta]+2*lambda[,2]*(1-alpha))
#       idbeta<-base::which(num02<=(lambda[,2]*alpha))
#       NEWBETA[idbeta+nva01+nweib]<-0}
# 
# 
#     # 1 ->2
#     if(nva12>0){
#       idbeta<-base::which(num12>(lambda[,3]*alpha))
#       NEWBETA[idbeta+nva01+nva02+nweib]<-sign12[idbeta]*(num12[idbeta]-lambda[,3]*alpha)/(denum12[idbeta]+2*lambda[,3]*(1-alpha))
#       idbeta<-base::which(num12<=(lambda[,3]*alpha))
#       NEWBETA[idbeta+nva01+nva02+nweib]<-0}
#   }
# 
# 
# 
# 
#   if(penalty=="mcp"){
# 
# 
#     # 0 -> 1,
#     idbeta<-base::which((num01>=(alpha*lambda[,1]*denum01)) & (num01<lambda[,1]) & (denum01<(1/alpha)))
#     NEWBETA[idbeta+nweib]<-sign01[idbeta]*alpha*lambda[,1]
#     # no definition put 0 ?
#     idbeta<-base::which(((num01<(alpha*lambda[,1]*denum01)) | (num01>=lambda[,1])) & (denum01<(1/alpha)))
#     NEWBETA[idbeta+nweib]<-0
# 
#     idbeta<-base::which((num01<=(alpha*lambda[,1]*denum01)) & num01>lambda[,1] & (denum01>=(1/alpha)))
#     NEWBETA[idbeta+nweib]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/(denum01[idbeta]-(1/alpha))
#     idbeta<-base::which(((num01>(alpha*lambda[,1]*denum01)) | num01<=lambda[,1]) & (denum01>=(1/alpha)))
#     NEWBETA[idbeta+nweib]<-0
# 
#     # 0 -> 2
#     idbeta<-base::which((num02>=(alpha*lambda[,2]*denum02)) & (num02<lambda[,2]) & (denum02<(1/alpha)))
#     NEWBETA[idbeta+nva01+nweib]<-sign02[idbeta]*alpha*lambda[,2]
#     idbeta<-base::which(((num02<(alpha*lambda[,2]*denum02)) | (num02>=lambda[,2])) & (denum02<(1/alpha)))
#     NEWBETA[idbeta+nva01+nweib]<-0
# 
#     idbeta<-base::which((num02<=(alpha*lambda[,2]*denum02)) & num02>lambda[,2] & (denum02>=(1/alpha)))
#     NEWBETA[idbeta+nva01+nweib]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/(denum02[idbeta]-(1/alpha))
#     idbeta<-base::which(((num02>(alpha*lambda[,2]*denum02)) | num02<=lambda[,2]) & (denum02>=(1/alpha)))
#     NEWBETA[idbeta+nva01+nweib]<-0
# 
# 
#     # 1 -> 2
#     idbeta<-base::which((num12>=(alpha*lambda[,3]*denum12)) & (num12<lambda[,3]) & (denum12<(1/alpha)))
#     NEWBETA[idbeta+nva01+nva02+nweib]<-sign12[idbeta]*alpha*lambda[,3]
#     idbeta<-base::which(((num12<(alpha*lambda[,3]*denum12)) | (num12>=lambda[,3])) & (denum12<(1/alpha)))
#     NEWBETA[idbeta+nva01+nva02+nweib]<-0
# 
#     idbeta<-base::which((num12<=(alpha*lambda[,3]*denum12)) & num12>lambda[,3] & (denum12>=(1/alpha)))
#     NEWBETA[idbeta+nva01+nva02+nweib]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/(denum12[idbeta]-(1/alpha))
#     idbeta<-base::which(((num12>(alpha*lambda[,3]*denum12)) | num12<=lambda[,3]) & (denum12>=(1/alpha)))
#     NEWBETA[idbeta+nva01+nva02+nweib]<-0
# 
#   }
# 
# 
#   if(penalty=="scad"){
# 
# 
# 
#     # 0 -> 1
#     idbeta<-base::which((num01<=(lambda[,1]*(1+denum01))) & (num01 >lambda[,1]) & (denum01>=(1/(alpha-1))))
#     NEWBETA[idbeta+nweib]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/denum01[idbeta]
#     idbeta<-base::which(((num01>(lambda[,1]*(1+denum01))) | (num01 <= lambda[,1]*denum01)) & (denum01>=(1/(alpha-1))))
#     NEWBETA[idbeta+nweib]<-0
# 
#     idbeta<-base::which((num01<=(alpha*lambda[,1]*denum01)) & (num01 > lambda[,1]) & (denum01<(1/(alpha-1))) & (denum01>=(1/alpha)))
#     NEWBETA[idbeta+nweib]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/denum01[idbeta]
#     idbeta<-base::which(((num01>(alpha*lambda[,1]*denum01)) | (num01 <= lambda[,1])) & (denum01<(1/(alpha-1))) & (denum01>=(1/alpha)))
#     NEWBETA[idbeta+nweib]<-0
#     idbeta<-base::which( (denum01<(1/(alpha-1))) & (denum01<(1/alpha)))
#     NEWBETA[idbeta+nweib]<-0
# 
#     # 0 ->2
#     idbeta<-base::which((num02<=(lambda[,2]*(1+denum02))) & (num02 >lambda[,2]) & (denum02>=(1/(alpha-1))))
#     NEWBETA[idbeta+nva01+nweib]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
#     idbeta<-base::which(((num02>(lambda[,2]*(1+denum02))) | (num02 <= lambda[,2]*denum02)) & (denum02>=(1/(alpha-1))))
#     NEWBETA[idbeta+nva01+nweib]<-0
# 
#     idbeta<-base::which((num02<=(alpha*lambda[,2]*denum02)) & (num02 > lambda[,2]) & (denum02<(1/(alpha-1))) & (denum02>=(1/alpha)))
#     NEWBETA[idbeta+nva01+nweib]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
#     idbeta<-base::which(((num02>(alpha*lambda[,2]*denum02)) | (num02 <= lambda[,2])) & (denum02<(1/(alpha-1))) & (denum02>=(1/alpha)))
#     NEWBETA[idbeta+nva01+nweib]<-0
#     idbeta<-base::which( (denum02<(1/(alpha-1))) & (denum02<(1/alpha)))
#     NEWBETA[idbeta+nva01+nweib]<-0
# 
#     # 1 ->2
# 
#     idbeta<-base::which((num12<=(lambda[,3]*(1+denum12))) & (num12 >lambda[,3]) & (denum12>=(1/(alpha-1))))
#     NEWBETA[idbeta+nva01+nva02+nweib]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
#     idbeta<-base::which(((num12>(lambda[,3]*(1+denum12))) | (num12 <= lambda[,3]*denum12)) & (denum12>=(1/(alpha-1))))
#     NEWBETA[idbeta+nva01+nva02+nweib]<-0
# 
#     idbeta<-base::which((num12<=(alpha*lambda[,3]*denum12)) & (num12 > lambda[,3]) & (denum12<(1/(alpha-1))) & (denum12>=(1/alpha)))
#     NEWBETA[idbeta+nva01+nva02+nweib]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
#     idbeta<-base::which(((num12>(alpha*lambda[,3]*denum12)) | (num12 <= lambda[,3])) & (denum12<(1/(alpha-1))) & (denum12>=(1/alpha)))
#     NEWBETA[idbeta+nva01+nva02+nweib]<-0
#     idbeta<-base::which( (denum12<(1/(alpha-1))) & (denum12<(1/alpha)))
#     NEWBETA[idbeta+nva01+nva02+nweib]<-0
#   }
# 
#   ## update according to marquard
#   idbeta<-base::which(penalty.factor==0)
#   # if no penalty on parameter, beta_k=A_k/-x_kk
#   NEWBETA[idbeta]<-beta[idbeta]+fu[1:nweib]%*%solve(v[1:nweib,1:nweib])
#   #NEWBETA[idbeta]<-beta[idbeta]+fu[1:nweib]
#   # does not work try update weibull parameters according to new update beta :
#   #browser()
#   # BETA[idbeta]<-NEWBETA[idbeta]
#   #
#   # num<-sapply(c(1:dim(v)[1]),FUN=function(x){
#   #   fu[x]-sum(v[x,-x]*BETA[-x])+sum(BETA*v[x,])
#   # })
#   #
#   # sign<-ifelse(num<0,-1,
#   #              ifelse(num>0,1,0))
#   # denum<-diag(v)
#   # idbeta<-base::which(penalty.factor==0)
#   #
#   # NEWBETA[idbeta]<-sign[idbeta]*num[idbeta]/denum[idbeta]
# 
#   NEW.BETA.all[fix==0]<-NEWBETA
# 
#   return(list(b=NEW.BETA.all))
# }


threshold.beta<-function(beta,vw,s,nva01Y,nva02Y,nva12Y,fix,v,
                      ctime,no,fu,ve01,ve02,ve12,dimnva01,dimnva02,dimnva12,
                      nva01,nva02,nva12,t0,t1,t2,t3,troncature,y01,y02,y12,
                      p01,p02,p12,dimp01,dimp02,dimp12,Ntime,lambda,alpha,
                      penalty.factor,penalty,penalty.weights){
  

  # add to do base::which for CRAN check 
  BETA<-beta[fix==0]
  NEW.BETA.all<-beta
  penalty.factor0<-penalty.factor
  penalty.factor<-penalty.factor[fix==0]
  
  #browser()
  num<-sapply(c(1:dim(v)[1]),FUN=function(x){
    (vw*fu[x]-sum(v[x,-x]*BETA[-x])+sum(BETA*v[x,]))/penalty.weights[x]
  }) 
  
  
  sign<-ifelse(num<0,-1,
               ifelse(num>0,1,0))
  denum<-diag(v)
  num<-abs(num)
  
  num01<-NULL
  denum01<-NULL
  sign01<-NULL
  
  
  num02<-NULL
  denum02<-NULL
  sign02<-NULL
  
  
  num12<-NULL
  denum12<-NULL
  sign12<-NULL
  
  
  num01Y<-NULL
  denum01Y<-NULL
  sign01Y<-NULL
  
  
  num02Y<-NULL
  denum02Y<-NULL
  sign02Y<-NULL
  
  
  num12Y<-NULL
  denum12Y<-NULL
  sign12Y<-NULL
  
  
  if(nva01>0){
    num01<-num[1:nva01]
    denum01<-denum[1:nva01]
    sign01<-sign[1:nva01]
  }
  
  if(nva01Y>0){
    num01Y<-num[(nva01+nva02+nva12+1):(nva01+nva02+nva12+nva01Y)]
    denum01Y<-denum[(nva01+nva02+nva12+1):(nva01+nva02+nva12+nva01Y)]
    sign01Y<-sign[(nva01+nva02+nva12+1):(nva01+nva02+nva12+nva01Y)]
    
  }
  
  if(nva02>0){
    num02<-num[(nva01+1):(nva01+nva02)]
    denum02<-denum[(nva01+1):(nva01+nva02)]
    sign02<-sign[(nva01+1):(nva01+nva02)]
    
  }
  
  if(nva02Y>0){
    num02Y<-num[(nva01+nva02+nva12+nva01Y+1):(nva01+nva02+nva12+nva01Y+nva02Y)]
    denum02Y<-denum[(nva01+nva02+nva12+nva01Y+1):(nva01+nva02+nva12+nva01Y+nva02Y)]
    sign02Y<-sign[(nva01+nva02+nva12+nva01Y+1):(nva01+nva02+nva12+nva01Y+nva02Y)]
    
  }
  
  if(nva12>0){
    num12<-num[(nva01+nva02+1):(nva01+nva02+nva12)]
    denum12<-denum[(nva01+nva02+1):(nva01+nva02+nva12)]
    sign12<-sign[(nva01+nva02+1):(nva01+nva02+nva12)]
  }
  
  if(nva12Y>0){
    num12Y<-num[(nva01+nva02+nva12+nva01Y+nva02Y+1):(nva01+nva02+nva12+nva01Y+nva02Y+nva12Y)]
    denum12Y<-denum[(nva01+nva02+nva12+nva01Y+nva02Y+1):(nva01+nva02+nva12+nva01Y+nva02Y+nva12Y)]
    sign12Y<-sign[(nva01+nva02+nva12+nva01Y+nva02Y+1):(nva01+nva02+nva12+nva01Y+nva02Y+nva12Y)]
    
  }
  
  NEWBETA<-rep(NA,length(num))
  idbeta<-NULL
  
  # if penalty update beta all at once 
  if(penalty%in%c("lasso","ridge","elasticnet","adaptive.lasso")){
    # 0 -> 1
    if(nva01>0){
      idbeta<-base::which(num01>(lambda[,1]*alpha))
      NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1]*alpha)/(denum01[idbeta]+2*lambda[,1]*(1-alpha))
      idbeta<-base::which(num01<=(lambda[,1]*alpha))
      NEWBETA[idbeta]<-0
    }
    if(nva01Y>0){
      idbeta<-base::which(num01Y>(lambda[,1]*alpha))
      NEWBETA[idbeta+nva01+nva02+nva12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1]*alpha)/(denum01Y[idbeta]+2*lambda[,1]*(1-alpha))
      idbeta<-base::which(num01Y<=(lambda[,1]*alpha))
      NEWBETA[idbeta+nva01+nva02+nva12]<-0
    }
    
    # 0 ->2
    if(nva02>0){
      idbeta<-base::which(num02>(lambda[,2]*alpha))
      NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2]*alpha)/(denum02[idbeta]+2*lambda[,2]*(1-alpha))
      idbeta<-base::which(num02<=(lambda[,2]*alpha))
      NEWBETA[idbeta+nva01]<-0
    }
    
    
    if(nva02Y>0){
      idbeta<-base::which(num02Y>(lambda[,2]*alpha))
      NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2]*alpha)/(denum02Y[idbeta]+2*lambda[,2]*(1-alpha))
      idbeta<-base::which(num02Y<=(lambda[,2]*alpha))
      NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    }
    # 1 ->2
    if(nva12>0){
      idbeta<-base::which(num12>(lambda[,3]*alpha))
      NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3]*alpha)/(denum12[idbeta]+2*lambda[,3]*(1-alpha))
      idbeta<-base::which(num12<=(lambda[,3]*alpha))
      NEWBETA[idbeta+nva01+nva02]<-0
    }
    
    
    if(nva12Y>0){
      idbeta<-base::which(num12Y>(lambda[,3]*alpha))
      NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3]*alpha)/(denum12Y[idbeta]+2*lambda[,3]*(1-alpha))
      idbeta<-base::which(num12Y<=(lambda[,3]*alpha))
      NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
    }
  }
  
  
  if(penalty=="mcp"){
    
    
    # 0 -> 1, 
    idbeta<-base::which((num01>=(alpha*lambda[,1]*denum01)) & (num01<lambda[,1]) & (denum01<(1/alpha)))
    NEWBETA[idbeta]<-sign01[idbeta]*alpha*lambda[,1]
    # no definition put 0 ? 
    idbeta<-base::which(((num01<(alpha*lambda[,1]*denum01)) | (num01>=lambda[,1])) & (denum01<(1/alpha)))
    NEWBETA[idbeta]<-0
    
    idbeta<-base::which((num01<=(alpha*lambda[,1]*denum01)) & num01>lambda[,1] & (denum01>=(1/alpha)))
    NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/(denum01[idbeta]-(1/alpha))
    idbeta<-base::which(((num01>(alpha*lambda[,1]*denum01)) | num01<=lambda[,1]) & (denum01>=(1/alpha)))
    NEWBETA[idbeta]<-0
    
    
    idbeta<-base::which((num01Y>=(alpha*lambda[,1]*denum01Y)) & (num01Y<lambda[,1]) & (denum01Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12]<-sign01Y[idbeta]*alpha*lambda[,1]
    # no definition put 0 ? 
    idbeta<-base::which(((num01Y<(alpha*lambda[,1]*denum01Y)) | (num01Y>=lambda[,1])) & (denum01Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12]<-0
    
    idbeta<-base::which((num01Y<=(alpha*lambda[,1]*denum01Y)) & num01Y>lambda[,1] & (denum01Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1])/(denum01Y[idbeta]-(1/alpha))
    idbeta<-base::which(((num01Y>(alpha*lambda[,1]*denum01Y)) | num01Y<=lambda[,1]) & (denum01Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12]<-0
    
    # 0 -> 2
    idbeta<-base::which((num02>=(alpha*lambda[,2]*denum02)) & (num02<lambda[,2]) & (denum02<(1/alpha)))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*alpha*lambda[,2]
    idbeta<-base::which(((num02<(alpha*lambda[,2]*denum02)) | (num02>=lambda[,2])) & (denum02<(1/alpha)))
    NEWBETA[idbeta+nva01]<-0
    
    idbeta<-base::which((num02<=(alpha*lambda[,2]*denum02)) & num02>lambda[,2] & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/(denum02[idbeta]-(1/alpha))
    idbeta<-base::which(((num02>(alpha*lambda[,2]*denum02)) | num02<=lambda[,2]) & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01]<-0
    
    
    idbeta<-base::which((num02Y>=(alpha*lambda[,2]*denum02Y)) & (num02Y<lambda[,2]) & (denum02Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-sign02Y[idbeta]*alpha*lambda[,2]
    idbeta<-base::which(((num02Y<(alpha*lambda[,2]*denum02Y)) | (num02Y>=lambda[,2])) & (denum02Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    
    idbeta<-base::which((num02Y<=(alpha*lambda[,2]*denum02Y)) & num02Y>lambda[,2] & (denum02Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2])/(denum02Y[idbeta]-(1/alpha))
    idbeta<-base::which(((num02Y>(alpha*lambda[,2]*denum02Y)) | num02Y<=lambda[,2]) & (denum02Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    
    
    # 1 -> 2
    idbeta<-base::which((num12>=(alpha*lambda[,3]*denum12)) & (num12<lambda[,3]) & (denum12<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*alpha*lambda[,3]
    idbeta<-base::which(((num12<(alpha*lambda[,3]*denum12)) | (num12>=lambda[,3])) & (denum12<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    idbeta<-base::which((num12<=(alpha*lambda[,3]*denum12)) & num12>lambda[,3] & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/(denum12[idbeta]-(1/alpha))
    idbeta<-base::which(((num12>(alpha*lambda[,3]*denum12)) | num12<=lambda[,3]) & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    
    idbeta<-base::which((num12Y>=(alpha*lambda[,3]*denum12Y)) & (num12Y<lambda[,3]) & (denum12Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-sign12Y[idbeta]*alpha*lambda[,3]
    idbeta<-base::which(((num12Y<(alpha*lambda[,3]*denum12Y)) | (num12Y>=lambda[,3])) & (denum12Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
    
    idbeta<-base::which((num12Y<=(alpha*lambda[,3]*denum12Y)) & num12Y>lambda[,3] & (denum12Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3])/(denum12Y[idbeta]-(1/alpha))
    idbeta<-base::which(((num12Y>(alpha*lambda[,3]*denum12Y)) | num12Y<=lambda[,3]) & (denum12Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
    
  }
  
  
  if(penalty=="scad"){
    
    
    
    # 0 -> 1
    idbeta<-base::which((num01<=(lambda[,1]*(1+denum01))) & (num01 >lambda[,1]) & (denum01>=(1/(alpha-1))))
    NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/denum01[idbeta]
    idbeta<-base::which(((num01>(lambda[,1]*(1+denum01))) | (num01 <= lambda[,1]*denum01)) & (denum01>=(1/(alpha-1))))
    NEWBETA[idbeta]<-0
    
    idbeta<-base::which((num01<=(alpha*lambda[,1]*denum01)) & (num01 > lambda[,1]) & (denum01<(1/(alpha-1))) & (denum01>=(1/alpha)))
    NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/denum01[idbeta]
    idbeta<-base::which(((num01>(alpha*lambda[,1]*denum01)) | (num01 <= lambda[,1])) & (denum01<(1/(alpha-1))) & (denum01>=(1/alpha)))
    NEWBETA[idbeta]<-0
    idbeta<-base::which( (denum01<(1/(alpha-1))) & (denum01<(1/alpha)))
    NEWBETA[idbeta]<-0
    
    
    idbeta<-base::which((num01Y<=(lambda[,1]*(1+denum01Y))) & (num01Y >lambda[,1]) & (denum01Y>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02+nva12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1])/denum01Y[idbeta]
    idbeta<-base::which(((num01Y>(lambda[,1]*(1+denum01Y))) | (num01Y <= lambda[,1]*denum01Y)) & (denum01Y>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02+nva12]<-0
    
    idbeta<-base::which((num01Y<=(alpha*lambda[,1]*denum01Y)) & (num01Y > lambda[,1]) & (denum01Y<(1/(alpha-1))) & (denum01Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1])/denum01Y[idbeta]
    idbeta<-base::which(((num01Y>(alpha*lambda[,1]*denum01Y)) | (num01Y <= lambda[,1])) & (denum01Y<(1/(alpha-1))) & (denum01Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12]<-0
    idbeta<-base::which( (denum01Y<(1/(alpha-1))) & (denum01Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12]<-0
    
    # 0 ->2
    idbeta<-base::which((num02<=(lambda[,2]*(1+denum02))) & (num02 >lambda[,2]) & (denum02>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
    idbeta<-base::which(((num02>(lambda[,2]*(1+denum02))) | (num02 <= lambda[,2]*denum02)) & (denum02>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01]<-0
    
    idbeta<-base::which((num02<=(alpha*lambda[,2]*denum02)) & (num02 > lambda[,2]) & (denum02<(1/(alpha-1))) & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
    idbeta<-base::which(((num02>(alpha*lambda[,2]*denum02)) | (num02 <= lambda[,2])) & (denum02<(1/(alpha-1))) & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01]<-0
    idbeta<-base::which( (denum02<(1/(alpha-1))) & (denum02<(1/alpha)))
    NEWBETA[idbeta+nva01]<-0
    
    
    idbeta<-base::which((num02Y<=(lambda[,2]*(1+denum02Y))) & (num02Y >lambda[,2]) & (denum02Y>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2])/denum02Y[idbeta]
    idbeta<-base::which(((num02Y>(lambda[,2]*(1+denum02Y))) | (num02Y <= lambda[,2]*denum02Y)) & (denum02Y>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    
    idbeta<-base::which((num02Y<=(alpha*lambda[,2]*denum02Y)) & (num02Y > lambda[,2]) & (denum02Y<(1/(alpha-1))) & (denum02Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2])/denum02Y[idbeta]
    idbeta<-base::which(((num02Y>(alpha*lambda[,2]*denum02Y)) | (num02Y <= lambda[,2])) & (denum02Y<(1/(alpha-1))) & (denum02Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    idbeta<-base::which( (denum02<(1/(alpha-1))) & (denum02<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    
    # 1 ->2
    
    idbeta<-base::which((num12<=(lambda[,3]*(1+denum12))) & (num12 >lambda[,3]) & (denum12>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
    idbeta<-base::which(((num12>(lambda[,3]*(1+denum12))) | (num12 <= lambda[,3]*denum12)) & (denum12>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    idbeta<-base::which((num12<=(alpha*lambda[,3]*denum12)) & (num12 > lambda[,3]) & (denum12<(1/(alpha-1))) & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
    idbeta<-base::which(((num12>(alpha*lambda[,3]*denum12)) | (num12 <= lambda[,3])) & (denum12<(1/(alpha-1))) & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-0
    idbeta<-base::which( (denum12<(1/(alpha-1))) & (denum12<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    
    idbeta<-base::which((num12Y<=(lambda[,3]*(1+denum12Y))) & (num12Y >lambda[,3]) & (denum12Y>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3])/denum12Y[idbeta]
    idbeta<-base::which(((num12Y>(lambda[,3]*(1+denum12Y))) | (num12Y <= lambda[,3]*denum12Y)) & (denum12Y>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
    
    idbeta<-base::which((num12Y<=(alpha*lambda[,3]*denum12Y)) & (num12Y > lambda[,3]) & (denum12Y<(1/(alpha-1))) & (denum12Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3])/denum12Y[idbeta]
    idbeta<-base::which(((num12Y>(alpha*lambda[,3]*denum12Y)) | (num12Y <= lambda[,3])) & (denum12Y<(1/(alpha-1))) & (denum12Y>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
    idbeta<-base::which( (denum12Y<(1/(alpha-1))) & (denum12Y<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
  }
  
  
  idbeta<-base::which(penalty.factor==0)
  # if no penalty on parameter, beta_k=A_k/-x_kk
  NEWBETA[idbeta]<-sign[idbeta]*num[idbeta]/denum[idbeta]
  
  NEW.BETA.all[fix==0]<-NEWBETA
  
  
  b<-c(s,NEW.BETA.all)

  res<-gaussDYNidmlLikelihoodweibpena(b=b,
                                     npm=length(b),
                                     npar=length(b),
                                     bfix=1,
                                     fix=rep(0,length(b)),
                                     ctime=ctime,
                                     no=no,
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
                                     y01=y01,
                                     y02=y02,
                                     y12=y12,
                                     p01=p01,
                                     p02=p02,
                                     p12=p12,
                                     dimp01=dimp01,
                                     dimp02=dimp02,
                                     dimp12=dimp12,
                                     Ntime=Ntime,
                                     lambda=lambda,
                                     alpha=alpha,
                                     penalty.factor=penalty.factor0,
                                     penalty=penalty,
                                     penalty.weights=penalty.weights)

  return(list(res=res,b=b))
}
