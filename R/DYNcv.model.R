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
  penalty.weights<-penalty.factor[fix==0]
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
      idbeta<-base::which(num01>(lambda[,1]*alpha[,1]))
      NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1]*alpha[,1])/(denum01[idbeta]+2*lambda[,1]*(1-alpha[,1]))
      idbeta<-base::which(num01<=(lambda[,1]*alpha[,1]))
      NEWBETA[idbeta]<-0
    }
    if(nva01Y>0){
      idbeta<-base::which(num01Y>(lambda[,1]*alpha[,1]))
      NEWBETA[idbeta+nva01+nva02+nva12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1]*alpha[,1])/(denum01Y[idbeta]+2*lambda[,1]*(1-alpha[,1]))
      idbeta<-base::which(num01Y<=(lambda[,1]*alpha[,1]))
      NEWBETA[idbeta+nva01+nva02+nva12]<-0
    }
    
    # 0 ->2
    if(nva02>0){
      idbeta<-base::which(num02>(lambda[,2]*alpha[,2]))
      NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2]*alpha[,2])/(denum02[idbeta]+2*lambda[,2]*(1-alpha[,2]))
      idbeta<-base::which(num02<=(lambda[,2]*alpha[,2]))
      NEWBETA[idbeta+nva01]<-0
    }
    
    
    if(nva02Y>0){
      idbeta<-base::which(num02Y>(lambda[,2]*alpha[,2]))
      NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2]*alpha[,2])/(denum02Y[idbeta]+2*lambda[,2]*(1-alpha[,2]))
      idbeta<-base::which(num02Y<=(lambda[,2]*alpha[,2]))
      NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    }
    # 1 ->2
    if(nva12>0){
      idbeta<-base::which(num12>(lambda[,3]*alpha[,3]))
      NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3]*alpha[,3])/(denum12[idbeta]+2*lambda[,3]*(1-alpha[,3]))
      idbeta<-base::which(num12<=(lambda[,3]*alpha[,3]))
      NEWBETA[idbeta+nva01+nva02]<-0
    }
    
    
    if(nva12Y>0){
      idbeta<-base::which(num12Y>(lambda[,3]*alpha[,3]))
      NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3]*alpha[,3])/(denum12Y[idbeta]+2*lambda[,3]*(1-alpha[,3]))
      idbeta<-base::which(num12Y<=(lambda[,3]*alpha[,3]))
      NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
    }
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
    
    
    idbeta<-base::which((num01Y>=(alpha[,1]*lambda[,1]*denum01Y)) & (num01Y<lambda[,1]) & (denum01Y<(1/alpha[,1])))
    NEWBETA[idbeta+nva01+nva02+nva12]<-sign01Y[idbeta]*alpha[,1]*lambda[,1]
    # no definition put 0 ? 
    idbeta<-base::which(((num01Y<(alpha[,1]*lambda[,1]*denum01Y)) | (num01Y>=lambda[,1])) & (denum01Y<(1/alpha[,1])))
    NEWBETA[idbeta+nva01+nva02+nva12]<-0
    
    idbeta<-base::which((num01Y<=(alpha[,1]*lambda[,1]*denum01Y)) & num01Y>lambda[,1] & (denum01Y>=(1/alpha[,1])))
    NEWBETA[idbeta+nva01+nva02+nva12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1])/(denum01Y[idbeta]-(1/alpha[,1]))
    idbeta<-base::which(((num01Y>(alpha[,1]*lambda[,1]*denum01Y)) | num01Y<=lambda[,1]) & (denum01Y>=(1/alpha[,1])))
    NEWBETA[idbeta+nva01+nva02+nva12]<-0
    
    # 0 -> 2
    idbeta<-base::which((num02>=(alpha[,2]*lambda[,2]*denum02)) & (num02<lambda[,2]) & (denum02<(1/alpha[,2])))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*alpha[,2]*lambda[,2]
    idbeta<-base::which(((num02<(alpha[,2]*lambda[,2]*denum02)) | (num02>=lambda[,2])) & (denum02<(1/alpha[,2])))
    NEWBETA[idbeta+nva01]<-0
    
    idbeta<-base::which((num02<=(alpha[,2]*lambda[,2]*denum02)) & num02>lambda[,2] & (denum02>=(1/alpha[,2])))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/(denum02[idbeta]-(1/alpha[,2]))
    idbeta<-base::which(((num02>(alpha[,2]*lambda[,2]*denum02)) | num02<=lambda[,2]) & (denum02>=(1/alpha[,2])))
    NEWBETA[idbeta+nva01]<-0
    
    
    idbeta<-base::which((num02Y>=(alpha[,2]*lambda[,2]*denum02Y)) & (num02Y<lambda[,2]) & (denum02Y<(1/alpha[,2])))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-sign02Y[idbeta]*alpha[,2]*lambda[,2]
    idbeta<-base::which(((num02Y<(alpha[,2]*lambda[,2]*denum02Y)) | (num02Y>=lambda[,2])) & (denum02Y<(1/alpha[,2])))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    
    idbeta<-base::which((num02Y<=(alpha[,2]*lambda[,2]*denum02Y)) & num02Y>lambda[,2] & (denum02Y>=(1/alpha[,2])))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2])/(denum02Y[idbeta]-(1/alpha[,2]))
    idbeta<-base::which(((num02Y>(alpha[,2]*lambda[,2]*denum02Y)) | num02Y<=lambda[,2]) & (denum02Y>=(1/alpha[,2])))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    
    
    # 1 -> 2
    idbeta<-base::which((num12>=(alpha[,3]*lambda[,3]*denum12)) & (num12<lambda[,3]) & (denum12<(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*alpha[,3]*lambda[,3]
    idbeta<-base::which(((num12<(alpha[,3]*lambda[,3]*denum12)) | (num12>=lambda[,3])) & (denum12<(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    idbeta<-base::which((num12<=(alpha[,3]*lambda[,3]*denum12)) & num12>lambda[,3] & (denum12>=(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/(denum12[idbeta]-(1/alpha[,3]))
    idbeta<-base::which(((num12>(alpha[,3]*lambda[,3]*denum12)) | num12<=lambda[,3]) & (denum12>=(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    
    idbeta<-base::which((num12Y>=(alpha[,3]*lambda[,3]*denum12Y)) & (num12Y<lambda[,3]) & (denum12Y<(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-sign12Y[idbeta]*alpha[,3]*lambda[,3]
    idbeta<-base::which(((num12Y<(alpha[,3]*lambda[,3]*denum12Y)) | (num12Y>=lambda[,3])) & (denum12Y<(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
    
    idbeta<-base::which((num12Y<=(alpha[,3]*lambda[,3]*denum12Y)) & num12Y>lambda[,3] & (denum12Y>=(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3])/(denum12Y[idbeta]-(1/alpha[,3]))
    idbeta<-base::which(((num12Y>(alpha[,3]*lambda[,3]*denum12Y)) | num12Y<=lambda[,3]) & (denum12Y>=(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
    
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
    
    
    idbeta<-base::which((num01Y<=(lambda[,1]*(1+denum01Y))) & (num01Y >lambda[,1]) & (denum01Y>=(1/(alpha[,1]-1))))
    NEWBETA[idbeta+nva01+nva02+nva12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1])/denum01Y[idbeta]
    idbeta<-base::which(((num01Y>(lambda[,1]*(1+denum01Y))) | (num01Y <= lambda[,1]*denum01Y)) & (denum01Y>=(1/(alpha[,1]-1))))
    NEWBETA[idbeta+nva01+nva02+nva12]<-0
    
    idbeta<-base::which((num01Y<=(alpha[,1]*lambda[,1]*denum01Y)) & (num01Y > lambda[,1]) & (denum01Y<(1/(alpha[,1]-1))) & (denum01Y>=(1/alpha[,1])))
    NEWBETA[idbeta+nva01+nva02+nva12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1])/denum01Y[idbeta]
    idbeta<-base::which(((num01Y>(alpha[,1]*lambda[,1]*denum01Y)) | (num01Y <= lambda[,1])) & (denum01Y<(1/(alpha[,1]-1))) & (denum01Y>=(1/alpha[,1])))
    NEWBETA[idbeta+nva01+nva02+nva12]<-0
    idbeta<-base::which( (denum01Y<(1/(alpha[,1]-1))) & (denum01Y<(1/alpha[,1])))
    NEWBETA[idbeta+nva01+nva02+nva12]<-0
    
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
    
    
    idbeta<-base::which((num02Y<=(lambda[,2]*(1+denum02Y))) & (num02Y >lambda[,2]) & (denum02Y>=(1/(alpha[,2]-1))))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2])/denum02Y[idbeta]
    idbeta<-base::which(((num02Y>(lambda[,2]*(1+denum02Y))) | (num02Y <= lambda[,2]*denum02Y)) & (denum02Y>=(1/(alpha[,2]-1))))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    
    idbeta<-base::which((num02Y<=(alpha[,2]*lambda[,2]*denum02Y)) & (num02Y > lambda[,2]) & (denum02Y<(1/(alpha[,2]-1))) & (denum02Y>=(1/alpha[,2])))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2])/denum02Y[idbeta]
    idbeta<-base::which(((num02Y>(alpha[,2]*lambda[,2]*denum02Y)) | (num02Y <= lambda[,2])) & (denum02Y<(1/(alpha[,2]-1))) & (denum02Y>=(1/alpha[,2])))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    idbeta<-base::which( (denum02<(1/(alpha[,2]-1))) & (denum02<(1/alpha[,2])))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y]<-0
    
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
    
    
    idbeta<-base::which((num12Y<=(lambda[,3]*(1+denum12Y))) & (num12Y >lambda[,3]) & (denum12Y>=(1/(alpha[,3]-1))))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3])/denum12Y[idbeta]
    idbeta<-base::which(((num12Y>(lambda[,3]*(1+denum12Y))) | (num12Y <= lambda[,3]*denum12Y)) & (denum12Y>=(1/(alpha[,3]-1))))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
    
    idbeta<-base::which((num12Y<=(alpha[,3]*lambda[,3]*denum12Y)) & (num12Y > lambda[,3]) & (denum12Y<(1/(alpha[,3]-1))) & (denum12Y>=(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3])/denum12Y[idbeta]
    idbeta<-base::which(((num12Y>(alpha[,3]*lambda[,3]*denum12Y)) | (num12Y <= lambda[,3])) & (denum12Y<(1/(alpha[,3]-1))) & (denum12Y>=(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
    idbeta<-base::which( (denum12Y<(1/(alpha[,3]-1))) & (denum12Y<(1/alpha[,3])))
    NEWBETA[idbeta+nva01+nva02+nva12+nva01Y+nva02Y]<-0
  }
  
  
  idbeta<-base::which(penalty.factor==0)
  # if no penalty on parameter, beta_k=A_k/-x_kk
  NEWBETA[idbeta]<-sign[idbeta]*num[idbeta]/denum[idbeta]
  
  NEW.BETA.all[fix==0]<-NEWBETA
  
  return(list(b=NEW.BETA.all))
}


threshold.beta_weib<-function(beta,vw,s,nva01Y,nva02Y,nva12Y,npm01,npm02,npm12,fix,v,
                              ctime,no,fu,ve01,ve02,ve12,dimnva01,dimnva02,dimnva12,
                              nva01,nva02,nva12,t0,t1,t2,t3,troncature,y01,y02,y12,
                              p01,p02,p12,dimp01,dimp02,dimp12,Ntime,lambda,alpha,
                              penalty.factor,penalty,penalty.weights){
  
  
  # add to do base::which for CRAN check 
  BETA<-beta[fix==0]
  NEW.BETA.all<-beta
  penalty.factor0<-penalty.factor
  penalty.factor<-penalty.factor[fix==0]
  penalty.weights0<-penalty.weights
  penalty.weights<-penalty.weights[fix==0]
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
  
  
  if(npm01>0){
    num01<-num[1:npm01]
    denum01<-denum[1:npm01]
    sign01<-sign[1:npm01]
  }
  
  if(nva01Y>0){
    num01Y<-num[(npm01+npm02+npm12+1):(npm01+npm02+npm12+nva01Y)]
    denum01Y<-denum[(npm01+npm02+npm12+1):(npm01+npm02+npm12+nva01Y)]
    sign01Y<-sign[(npm01+npm02+npm12+1):(npm01+npm02+npm12+nva01Y)]
    
  }
  
  if(npm02>0){
    num02<-num[(npm01+1):(npm01+npm02)]
    denum02<-denum[(npm01+1):(npm01+npm02)]
    sign02<-sign[(npm01+1):(npm01+npm02)]
    
  }
  
  if(nva02Y>0){
    num02Y<-num[(npm01+npm02+npm12+nva01Y+1):(npm01+npm02+npm12+nva01Y+nva02Y)]
    denum02Y<-denum[(npm01+npm02+npm12+nva01Y+1):(npm01+npm02+npm12+nva01Y+nva02Y)]
    sign02Y<-sign[(npm01+npm02+npm12+nva01Y+1):(npm01+npm02+npm12+nva01Y+nva02Y)]
    
  }
  
  if(npm12>0){
    num12<-num[(npm01+npm02+1):(npm01+npm02+npm12)]
    denum12<-denum[(npm01+npm02+1):(npm01+npm02+npm12)]
    sign12<-sign[(npm01+npm02+1):(npm01+npm02+npm12)]
  }
  
  if(nva12Y>0){
    num12Y<-num[(npm01+npm02+npm12+nva01Y+nva02Y+1):(npm01+npm02+npm12+nva01Y+nva02Y+nva12Y)]
    denum12Y<-denum[(npm01+npm02+npm12+nva01Y+nva02Y+1):(npm01+npm02+npm12+nva01Y+nva02Y+nva12Y)]
    sign12Y<-sign[(npm01+npm02+npm12+nva01Y+nva02Y+1):(npm01+npm02+npm12+nva01Y+nva02Y+nva12Y)]
    
  }
  
  NEWBETA<-rep(NA,length(num))
  idbeta<-NULL
  
  # if penalty update beta all at once 
  if(penalty%in%c("lasso","ridge","elasticnet","adaptive.lasso")){
    # 0 -> 1
    if(npm01>0){
      idbeta<-base::which(num01>(lambda[,1]*alpha[,1]))
      NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1]*alpha[,1])/(denum01[idbeta]+2*lambda[,1]*(1-alpha[,1]))
      idbeta<-base::which(num01<=(lambda[,1]*alpha[,1]))
      NEWBETA[idbeta]<-0
    }
    if(nva01Y>0){
      idbeta<-base::which(num01Y>(lambda[,1]*alpha[,1]))
      NEWBETA[idbeta+npm01+npm02+npm12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1]*alpha[,1])/(denum01Y[idbeta]+2*lambda[,1]*(1-alpha[,1]))
      idbeta<-base::which(num01Y<=(lambda[,1]*alpha[,1]))
      NEWBETA[idbeta+npm01+npm02+npm12]<-0
    }
    
    # 0 ->2
    if(npm02>0){
      idbeta<-base::which(num02>(lambda[,2]*alpha[,2]))
      NEWBETA[idbeta+npm01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2]*alpha[,2])/(denum02[idbeta]+2*lambda[,2]*(1-alpha[,2]))
      idbeta<-base::which(num02<=(lambda[,2]*alpha[,2]))
      NEWBETA[idbeta+npm01]<-0
    }
    
    
    if(nva02Y>0){
      idbeta<-base::which(num02Y>(lambda[,2]*alpha[,2]))
      NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2]*alpha[,2])/(denum02Y[idbeta]+2*lambda[,2]*(1-alpha[,2]))
      idbeta<-base::which(num02Y<=(lambda[,2]*alpha[,2]))
      NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-0
    }
    # 1 ->2
    if(nva12>0){
      idbeta<-base::which(num12>(lambda[,3]*alpha[,3]))
      NEWBETA[idbeta+npm01+npm02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3]*alpha[,3])/(denum12[idbeta]+2*lambda[,3]*(1-alpha[,3]))
      idbeta<-base::which(num12<=(lambda[,3]*alpha[,3]))
      NEWBETA[idbeta+npm01+npm02]<-0
    }
    
    
    if(nva12Y>0){
      idbeta<-base::which(num12Y>(lambda[,3]*alpha[,3]))
      NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3]*alpha[,3])/(denum12Y[idbeta]+2*lambda[,3]*(1-alpha[,3]))
      idbeta<-base::which(num12Y<=(lambda[,3]*alpha[,3]))
      NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-0
    }
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
    
    
    idbeta<-base::which((num01Y>=(alpha[,1]*lambda[,1]*denum01Y)) & (num01Y<lambda[,1]) & (denum01Y<(1/alpha[,1])))
    NEWBETA[idbeta+npm01+npm02+npm12]<-sign01Y[idbeta]*alpha[,1]*lambda[,1]
    # no definition put 0 ? 
    idbeta<-base::which(((num01Y<(alpha[,1]*lambda[,1]*denum01Y)) | (num01Y>=lambda[,1])) & (denum01Y<(1/alpha[,1])))
    NEWBETA[idbeta+npm01+npm02+npm12]<-0
    
    idbeta<-base::which((num01Y<=(alpha[,1]*lambda[,1]*denum01Y)) & num01Y>lambda[,1] & (denum01Y>=(1/alpha[,1])))
    NEWBETA[idbeta+npm01+npm02+npm12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1])/(denum01Y[idbeta]-(1/alpha[,1]))
    idbeta<-base::which(((num01Y>(alpha[,1]*lambda[,1]*denum01Y)) | num01Y<=lambda[,1]) & (denum01Y>=(1/alpha[,1])))
    NEWBETA[idbeta+npm01+npm02+npm12]<-0
    
    # 0 -> 2
    idbeta<-base::which((num02>=(alpha[,2]*lambda[,2]*denum02)) & (num02<lambda[,2]) & (denum02<(1/alpha[,2])))
    NEWBETA[idbeta+npm01]<-sign02[idbeta]*alpha[,2]*lambda[,2]
    idbeta<-base::which(((num02<(alpha[,2]*lambda[,2]*denum02)) | (num02>=lambda[,2])) & (denum02<(1/alpha[,2])))
    NEWBETA[idbeta+npm01]<-0
    
    idbeta<-base::which((num02<=(alpha[,2]*lambda[,2]*denum02)) & num02>lambda[,2] & (denum02>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/(denum02[idbeta]-(1/alpha[,2]))
    idbeta<-base::which(((num02>(alpha[,2]*lambda[,2]*denum02)) | num02<=lambda[,2]) & (denum02>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01]<-0
    
    
    idbeta<-base::which((num02Y>=(alpha[,2]*lambda[,2]*denum02Y)) & (num02Y<lambda[,2]) & (denum02Y<(1/alpha[,2])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-sign02Y[idbeta]*alpha[,2]*lambda[,2]
    idbeta<-base::which(((num02Y<(alpha[,2]*lambda[,2]*denum02Y)) | (num02Y>=lambda[,2])) & (denum02Y<(1/alpha[,2])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-0
    
    idbeta<-base::which((num02Y<=(alpha[,2]*lambda[,2]*denum02Y)) & num02Y>lambda[,2] & (denum02Y>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2])/(denum02Y[idbeta]-(1/alpha[,2]))
    idbeta<-base::which(((num02Y>(alpha[,2]*lambda[,2]*denum02Y)) | num02Y<=lambda[,2]) & (denum02Y>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-0
    
    
    # 1 -> 2
    idbeta<-base::which((num12>=(alpha[,3]*lambda[,3]*denum12)) & (num12<lambda[,3]) & (denum12<(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02]<-sign12[idbeta]*alpha[,3]*lambda[,3]
    idbeta<-base::which(((num12<(alpha[,3]*lambda[,3]*denum12)) | (num12>=lambda[,3])) & (denum12<(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02]<-0
    
    idbeta<-base::which((num12<=(alpha[,3]*lambda[,3]*denum12)) & num12>lambda[,3] & (denum12>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/(denum12[idbeta]-(1/alpha[,3]))
    idbeta<-base::which(((num12>(alpha[,3]*lambda[,3]*denum12)) | num12<=lambda[,3]) & (denum12>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02]<-0
    
    
    idbeta<-base::which((num12Y>=(alpha[,3]*lambda[,3]*denum12Y)) & (num12Y<lambda[,3]) & (denum12Y<(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-sign12Y[idbeta]*alpha[,3]*lambda[,3]
    idbeta<-base::which(((num12Y<(alpha[,3]*lambda[,3]*denum12Y)) | (num12Y>=lambda[,3])) & (denum12Y<(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-0
    
    idbeta<-base::which((num12Y<=(alpha[,3]*lambda[,3]*denum12Y)) & num12Y>lambda[,3] & (denum12Y>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3])/(denum12Y[idbeta]-(1/alpha[,3]))
    idbeta<-base::which(((num12Y>(alpha[,3]*lambda[,3]*denum12Y)) | num12Y<=lambda[,3]) & (denum12Y>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-0
    
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
    
    
    idbeta<-base::which((num01Y<=(lambda[,1]*(1+denum01Y))) & (num01Y >lambda[,1]) & (denum01Y>=(1/(alpha[,1]-1))))
    NEWBETA[idbeta+npm01+npm02+npm12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1])/denum01Y[idbeta]
    idbeta<-base::which(((num01Y>(lambda[,1]*(1+denum01Y))) | (num01Y <= lambda[,1]*denum01Y)) & (denum01Y>=(1/(alpha[,1]-1))))
    NEWBETA[idbeta+npm01+npm02+npm12]<-0
    
    idbeta<-base::which((num01Y<=(alpha[,1]*lambda[,1]*denum01Y)) & (num01Y > lambda[,1]) & (denum01Y<(1/(alpha[,1]-1))) & (denum01Y>=(1/alpha[,1])))
    NEWBETA[idbeta+npm01+npm02+npm12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1])/denum01Y[idbeta]
    idbeta<-base::which(((num01Y>(alpha[,1]*lambda[,1]*denum01Y)) | (num01Y <= lambda[,1])) & (denum01Y<(1/(alpha[,1]-1))) & (denum01Y>=(1/alpha[,1])))
    NEWBETA[idbeta+npm01+npm02+npm12]<-0
    idbeta<-base::which( (denum01Y<(1/(alpha[,1]-1))) & (denum01Y<(1/alpha[,1])))
    NEWBETA[idbeta+npm01+npm02+npm12]<-0
    
    # 0 ->2
    idbeta<-base::which((num02<=(lambda[,2]*(1+denum02))) & (num02 >lambda[,2]) & (denum02>=(1/(alpha[,2]-1))))
    NEWBETA[idbeta+npm01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
    idbeta<-base::which(((num02>(lambda[,2]*(1+denum02))) | (num02 <= lambda[,2]*denum02)) & (denum02>=(1/(alpha[,2]-1))))
    NEWBETA[idbeta+npm01]<-0
    
    idbeta<-base::which((num02<=(alpha[,2]*lambda[,2]*denum02)) & (num02 > lambda[,2]) & (denum02<(1/(alpha[,2]-1))) & (denum02>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
    idbeta<-base::which(((num02>(alpha[,2]*lambda[,2]*denum02)) | (num02 <= lambda[,2])) & (denum02<(1/(alpha[,2]-1))) & (denum02>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01]<-0
    idbeta<-base::which( (denum02<(1/(alpha[,2]-1))) & (denum02<(1/alpha[,2])))
    NEWBETA[idbeta+npm01]<-0
    
    
    idbeta<-base::which((num02Y<=(lambda[,2]*(1+denum02Y))) & (num02Y >lambda[,2]) & (denum02Y>=(1/(alpha[,2]-1))))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2])/denum02Y[idbeta]
    idbeta<-base::which(((num02Y>(lambda[,2]*(1+denum02Y))) | (num02Y <= lambda[,2]*denum02Y)) & (denum02Y>=(1/(alpha[,2]-1))))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-0
    
    idbeta<-base::which((num02Y<=(alpha[,2]*lambda[,2]*denum02Y)) & (num02Y > lambda[,2]) & (denum02Y<(1/(alpha[,2]-1))) & (denum02Y>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2])/denum02Y[idbeta]
    idbeta<-base::which(((num02Y>(alpha[,2]*lambda[,2]*denum02Y)) | (num02Y <= lambda[,2])) & (denum02Y<(1/(alpha[,2]-1))) & (denum02Y>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-0
    idbeta<-base::which( (denum02<(1/(alpha[,2]-1))) & (denum02<(1/alpha[,2])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-0
    
    # 1 ->2
    
    idbeta<-base::which((num12<=(lambda[,3]*(1+denum12))) & (num12 >lambda[,3]) & (denum12>=(1/(alpha[,3]-1))))
    NEWBETA[idbeta+npm01+npm02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
    idbeta<-base::which(((num12>(lambda[,3]*(1+denum12))) | (num12 <= lambda[,3]*denum12)) & (denum12>=(1/(alpha[,3]-1))))
    NEWBETA[idbeta+npm01+npm02]<-0
    
    idbeta<-base::which((num12<=(alpha[,3]*lambda[,3]*denum12)) & (num12 > lambda[,3]) & (denum12<(1/(alpha[,3]-1))) & (denum12>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
    idbeta<-base::which(((num12>(alpha[,3]*lambda[,3]*denum12)) | (num12 <= lambda[,3])) & (denum12<(1/(alpha[,3]-1))) & (denum12>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02]<-0
    idbeta<-base::which( (denum12<(1/(alpha[,3]-1))) & (denum12<(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02]<-0
    
    
    idbeta<-base::which((num12Y<=(lambda[,3]*(1+denum12Y))) & (num12Y >lambda[,3]) & (denum12Y>=(1/(alpha[,3]-1))))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3])/denum12Y[idbeta]
    idbeta<-base::which(((num12Y>(lambda[,3]*(1+denum12Y))) | (num12Y <= lambda[,3]*denum12Y)) & (denum12Y>=(1/(alpha[,3]-1))))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-0
    
    idbeta<-base::which((num12Y<=(alpha[,3]*lambda[,3]*denum12Y)) & (num12Y > lambda[,3]) & (denum12Y<(1/(alpha[,3]-1))) & (denum12Y>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3])/denum12Y[idbeta]
    idbeta<-base::which(((num12Y>(alpha[,3]*lambda[,3]*denum12Y)) | (num12Y <= lambda[,3])) & (denum12Y<(1/(alpha[,3]-1))) & (denum12Y>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-0
    idbeta<-base::which( (denum12Y<(1/(alpha[,3]-1))) & (denum12Y<(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-0
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
                                      penalty.weights=penalty.weights0)
  
  return(list(res=res,b=b))
}


threshold.beta<-function(beta,vw,s,nva01Y,nva02Y,nva12Y,npm01,npm02,npm12,fix,v,zi01,zi02,zi12,
                         ctime,no,nz01,nz02,nz12,
                         fu,ve01,ve02,ve12,dimnva01,dimnva02,dimnva12,
                         nva01,nva02,nva12,t0,t1,t2,t3,troncature,y01,y02,y12,
                         p01,p02,p12,dimp01,dimp02,dimp12,Ntime,lambda,alpha,
                         penalty.factor,penalty,penalty.weights){
  
  
  # add to do base::which for CRAN check 
  BETA<-beta[fix==0]
  NEW.BETA.all<-beta
  penalty.factor0<-penalty.factor
  penalty.factor<-penalty.factor[fix==0]
  penalty.weights0<-penalty.weights
  penalty.weights<-penalty.weights[fix==0]
  
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
  
  
  if(npm01>0){
    num01<-num[1:npm01]
    denum01<-denum[1:npm01]
    sign01<-sign[1:npm01]
  }
  
  if(nva01Y>0){
    num01Y<-num[(npm01+npm02+npm12+1):(npm01+npm02+npm12+nva01Y)]
    denum01Y<-denum[(npm01+npm02+npm12+1):(npm01+npm02+npm12+nva01Y)]
    sign01Y<-sign[(npm01+npm02+npm12+1):(npm01+npm02+npm12+nva01Y)]
    
  }
  
  if(npm02>0){
    num02<-num[(npm01+1):(npm01+npm02)]
    denum02<-denum[(npm01+1):(npm01+npm02)]
    sign02<-sign[(npm01+1):(npm01+npm02)]
    
  }
  
  if(nva02Y>0){
    num02Y<-num[(npm01+npm02+npm12+nva01Y+1):(npm01+npm02+npm12+nva01Y+nva02Y)]
    denum02Y<-denum[(npm01+npm02+npm12+nva01Y+1):(npm01+npm02+npm12+nva01Y+nva02Y)]
    sign02Y<-sign[(npm01+npm02+npm12+nva01Y+1):(npm01+npm02+npm12+nva01Y+nva02Y)]
    
  }
  
  if(npm12>0){
    num12<-num[(npm01+npm02+1):(npm01+npm02+npm12)]
    denum12<-denum[(npm01+npm02+1):(npm01+npm02+npm12)]
    sign12<-sign[(npm01+npm02+1):(npm01+npm02+npm12)]
  }
  
  if(nva12Y>0){
    num12Y<-num[(npm01+npm02+npm12+nva01Y+nva02Y+1):(npm01+npm02+npm12+nva01Y+nva02Y+nva12Y)]
    denum12Y<-denum[(npm01+npm02+npm12+nva01Y+nva02Y+1):(npm01+npm02+npm12+nva01Y+nva02Y+nva12Y)]
    sign12Y<-sign[(npm01+npm02+npm12+nva01Y+nva02Y+1):(npm01+npm02+npm12+nva01Y+nva02Y+nva12Y)]
    
  }
  
  NEWBETA<-rep(NA,length(num))
  idbeta<-NULL
  
  # if penalty update beta all at once 
  if(penalty%in%c("lasso","ridge","elasticnet","adaptive.lasso")){
    # 0 -> 1
    if(npm01>0){
      idbeta<-base::which(num01>(lambda[,1]*alpha[,1]))
      NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1]*alpha[,1])/(denum01[idbeta]+2*lambda[,1]*(1-alpha[,1]))
      idbeta<-base::which(num01<=(lambda[,1]*alpha[,1]))
      NEWBETA[idbeta]<-0
    }
    if(nva01Y>0){
      idbeta<-base::which(num01Y>(lambda[,1]*alpha[,1]))
      NEWBETA[idbeta+npm01+npm02+npm12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1]*alpha[,1])/(denum01Y[idbeta]+2*lambda[,1]*(1-alpha[,1]))
      idbeta<-base::which(num01Y<=(lambda[,1]*alpha[,1]))
      NEWBETA[idbeta+npm01+npm02+npm12]<-0
    }
    
    # 0 ->2
    if(npm02>0){
      idbeta<-base::which(num02>(lambda[,2]*alpha[,2]))
      NEWBETA[idbeta+npm01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2]*alpha[,2])/(denum02[idbeta]+2*lambda[,2]*(1-alpha[,2]))
      idbeta<-base::which(num02<=(lambda[,2]*alpha[,2]))
      NEWBETA[idbeta+npm01]<-0
    }
    
    
    if(nva02Y>0){
      idbeta<-base::which(num02Y>(lambda[,2]*alpha[,2]))
      NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2]*alpha[,2])/(denum02Y[idbeta]+2*lambda[,2]*(1-alpha[,2]))
      idbeta<-base::which(num02Y<=(lambda[,2]*alpha[,2]))
      NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-0
    }
    # 1 ->2
    if(npm12>0){
      idbeta<-base::which(num12>(lambda[,3]*alpha[,3]))
      NEWBETA[idbeta+npm01+npm02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3]*alpha[,3])/(denum12[idbeta]+2*lambda[,3]*(1-alpha[,3]))
      idbeta<-base::which(num12<=(lambda[,3]*alpha[,3]))
      NEWBETA[idbeta+npm01+npm02]<-0
    }
    
    
    if(nva12Y>0){
      idbeta<-base::which(num12Y>(lambda[,3]*alpha[,3]))
      NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3]*alpha[,3])/(denum12Y[idbeta]+2*lambda[,3]*(1-alpha[,3]))
      idbeta<-base::which(num12Y<=(lambda[,3]*alpha[,3]))
      NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-0
    }
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
    
    
    idbeta<-base::which((num01Y>=(alpha[,1]*lambda[,1]*denum01Y)) & (num01Y<lambda[,1]) & (denum01Y<(1/alpha[,1])))
    NEWBETA[idbeta+npm01+npm02+npm12]<-sign01Y[idbeta]*alpha[,1]*lambda[,1]
    # no definition put 0 ? 
    idbeta<-base::which(((num01Y<(alpha[,1]*lambda[,1]*denum01Y)) | (num01Y>=lambda[,1])) & (denum01Y<(1/alpha[,1])))
    NEWBETA[idbeta+npm01+npm02+npm12]<-0
    
    idbeta<-base::which((num01Y<=(alpha[,1]*lambda[,1]*denum01Y)) & num01Y>lambda[,1] & (denum01Y>=(1/alpha[,1])))
    NEWBETA[idbeta+npm01+npm02+npm12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1])/(denum01Y[idbeta]-(1/alpha[,1]))
    idbeta<-base::which(((num01Y>(alpha[,1]*lambda[,1]*denum01Y)) | num01Y<=lambda[,1]) & (denum01Y>=(1/alpha[,1])))
    NEWBETA[idbeta+npm01+npm02+npm12]<-0
    
    # 0 -> 2
    idbeta<-base::which((num02>=(alpha[,2]*lambda[,2]*denum02)) & (num02<lambda[,2]) & (denum02<(1/alpha[,2])))
    NEWBETA[idbeta+npm01]<-sign02[idbeta]*alpha[,2]*lambda[,2]
    idbeta<-base::which(((num02<(alpha[,2]*lambda[,2]*denum02)) | (num02>=lambda[,2])) & (denum02<(1/alpha[,2])))
    NEWBETA[idbeta+npm01]<-0
    
    idbeta<-base::which((num02<=(alpha[,2]*lambda[,2]*denum02)) & num02>lambda[,2] & (denum02>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/(denum02[idbeta]-(1/alpha[,2]))
    idbeta<-base::which(((num02>(alpha[,2]*lambda[,2]*denum02)) | num02<=lambda[,2]) & (denum02>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01]<-0
    
    
    idbeta<-base::which((num02Y>=(alpha[,2]*lambda[,2]*denum02Y)) & (num02Y<lambda[,2]) & (denum02Y<(1/alpha[,2])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-sign02Y[idbeta]*alpha[,2]*lambda[,2]
    idbeta<-base::which(((num02Y<(alpha[,2]*lambda[,2]*denum02Y)) | (num02Y>=lambda[,2])) & (denum02Y<(1/alpha[,2])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-0
    
    idbeta<-base::which((num02Y<=(alpha[,2]*lambda[,2]*denum02Y)) & num02Y>lambda[,2] & (denum02Y>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2])/(denum02Y[idbeta]-(1/alpha[,2]))
    idbeta<-base::which(((num02Y>(alpha[,2]*lambda[,2]*denum02Y)) | num02Y<=lambda[,2]) & (denum02Y>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-0
    
    
    # 1 -> 2
    idbeta<-base::which((num12>=(alpha[,3]*lambda[,3]*denum12)) & (num12<lambda[,3]) & (denum12<(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02]<-sign12[idbeta]*alpha[,3]*lambda[,3]
    idbeta<-base::which(((num12<(alpha[,3]*lambda[,3]*denum12)) | (num12>=lambda[,3])) & (denum12<(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02]<-0
    
    idbeta<-base::which((num12<=(alpha[,3]*lambda[,3]*denum12)) & num12>lambda[,3] & (denum12>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/(denum12[idbeta]-(1/alpha[,3]))
    idbeta<-base::which(((num12>(alpha[,3]*lambda[,3]*denum12)) | num12<=lambda[,3]) & (denum12>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02]<-0
    
    
    idbeta<-base::which((num12Y>=(alpha[,3]*lambda[,3]*denum12Y)) & (num12Y<lambda[,3]) & (denum12Y<(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-sign12Y[idbeta]*alpha[,3]*lambda[,3]
    idbeta<-base::which(((num12Y<(alpha[,3]*lambda[,3]*denum12Y)) | (num12Y>=lambda[,3])) & (denum12Y<(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-0
    
    idbeta<-base::which((num12Y<=(alpha[,3]*lambda[,3]*denum12Y)) & num12Y>lambda[,3] & (denum12Y>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3])/(denum12Y[idbeta]-(1/alpha[,3]))
    idbeta<-base::which(((num12Y>(alpha[,3]*lambda[,3]*denum12Y)) | num12Y<=lambda[,3]) & (denum12Y>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-0
    
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
    
    
    idbeta<-base::which((num01Y<=(lambda[,1]*(1+denum01Y))) & (num01Y >lambda[,1]) & (denum01Y>=(1/(alpha[,1]-1))))
    NEWBETA[idbeta+npm01+npm02+npm12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1])/denum01Y[idbeta]
    idbeta<-base::which(((num01Y>(lambda[,1]*(1+denum01Y))) | (num01Y <= lambda[,1]*denum01Y)) & (denum01Y>=(1/(alpha[,1]-1))))
    NEWBETA[idbeta+npm01+npm02+npm12]<-0
    
    idbeta<-base::which((num01Y<=(alpha[,1]*lambda[,1]*denum01Y)) & (num01Y > lambda[,1]) & (denum01Y<(1/(alpha[,1]-1))) & (denum01Y>=(1/alpha[,1])))
    NEWBETA[idbeta+npm01+npm02+npm12]<-sign01Y[idbeta]*(num01Y[idbeta]-lambda[,1])/denum01Y[idbeta]
    idbeta<-base::which(((num01Y>(alpha[,1]*lambda[,1]*denum01Y)) | (num01Y <= lambda[,1])) & (denum01Y<(1/(alpha[,1]-1))) & (denum01Y>=(1/alpha[,1])))
    NEWBETA[idbeta+npm01+npm02+npm12]<-0
    idbeta<-base::which( (denum01Y<(1/(alpha[,1]-1))) & (denum01Y<(1/alpha[,1])))
    NEWBETA[idbeta+npm01+npm02+npm12]<-0
    
    # 0 ->2
    idbeta<-base::which((num02<=(lambda[,2]*(1+denum02))) & (num02 >lambda[,2]) & (denum02>=(1/(alpha[,2]-1))))
    NEWBETA[idbeta+npm01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
    idbeta<-base::which(((num02>(lambda[,2]*(1+denum02))) | (num02 <= lambda[,2]*denum02)) & (denum02>=(1/(alpha[,2]-1))))
    NEWBETA[idbeta+npm01]<-0
    
    idbeta<-base::which((num02<=(alpha[,2]*lambda[,2]*denum02)) & (num02 > lambda[,2]) & (denum02<(1/(alpha[,2]-1))) & (denum02>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
    idbeta<-base::which(((num02>(alpha[,2]*lambda[,2]*denum02)) | (num02 <= lambda[,2])) & (denum02<(1/(alpha[,2]-1))) & (denum02>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01]<-0
    idbeta<-base::which( (denum02<(1/(alpha[,2]-1))) & (denum02<(1/alpha[,2])))
    NEWBETA[idbeta+npm01]<-0
    
    
    idbeta<-base::which((num02Y<=(lambda[,2]*(1+denum02Y))) & (num02Y >lambda[,2]) & (denum02Y>=(1/(alpha[,2]-1))))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2])/denum02Y[idbeta]
    idbeta<-base::which(((num02Y>(lambda[,2]*(1+denum02Y))) | (num02Y <= lambda[,2]*denum02Y)) & (denum02Y>=(1/(alpha[,2]-1))))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-0
    
    idbeta<-base::which((num02Y<=(alpha[,2]*lambda[,2]*denum02Y)) & (num02Y > lambda[,2]) & (denum02Y<(1/(alpha[,2]-1))) & (denum02Y>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-sign02Y[idbeta]*(num02Y[idbeta]-lambda[,2])/denum02Y[idbeta]
    idbeta<-base::which(((num02Y>(alpha[,2]*lambda[,2]*denum02Y)) | (num02Y <= lambda[,2])) & (denum02Y<(1/(alpha[,2]-1))) & (denum02Y>=(1/alpha[,2])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-0
    idbeta<-base::which( (denum02<(1/(alpha[,2]-1))) & (denum02<(1/alpha[,2])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y]<-0
    
    # 1 ->2
    
    idbeta<-base::which((num12<=(lambda[,3]*(1+denum12))) & (num12 >lambda[,3]) & (denum12>=(1/(alpha[,3]-1))))
    NEWBETA[idbeta+npm01+npm02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
    idbeta<-base::which(((num12>(lambda[,3]*(1+denum12))) | (num12 <= lambda[,3]*denum12)) & (denum12>=(1/(alpha[,3]-1))))
    NEWBETA[idbeta+npm01+npm02]<-0
    
    idbeta<-base::which((num12<=(alpha[,3]*lambda[,3]*denum12)) & (num12 > lambda[,3]) & (denum12<(1/(alpha[,3]-1))) & (denum12>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
    idbeta<-base::which(((num12>(alpha[,3]*lambda[,3]*denum12)) | (num12 <= lambda[,3])) & (denum12<(1/(alpha[,3]-1))) & (denum12>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02]<-0
    idbeta<-base::which( (denum12<(1/(alpha[,3]-1))) & (denum12<(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02]<-0
    
    
    idbeta<-base::which((num12Y<=(lambda[,3]*(1+denum12Y))) & (num12Y >lambda[,3]) & (denum12Y>=(1/(alpha[,3]-1))))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3])/denum12Y[idbeta]
    idbeta<-base::which(((num12Y>(lambda[,3]*(1+denum12Y))) | (num12Y <= lambda[,3]*denum12Y)) & (denum12Y>=(1/(alpha[,3]-1))))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-0
    
    idbeta<-base::which((num12Y<=(alpha[,3]*lambda[,3]*denum12Y)) & (num12Y > lambda[,3]) & (denum12Y<(1/(alpha[,3]-1))) & (denum12Y>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-sign12Y[idbeta]*(num12Y[idbeta]-lambda[,3])/denum12Y[idbeta]
    idbeta<-base::which(((num12Y>(alpha[,3]*lambda[,3]*denum12Y)) | (num12Y <= lambda[,3])) & (denum12Y<(1/(alpha[,3]-1))) & (denum12Y>=(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-0
    idbeta<-base::which( (denum12Y<(1/(alpha[,3]-1))) & (denum12Y<(1/alpha[,3])))
    NEWBETA[idbeta+npm01+npm02+npm12+nva01Y+nva02Y]<-0
  }
  
  
  idbeta<-base::which(penalty.factor==0)
  # if no penalty on parameter, beta_k=A_k/-x_kk
  NEWBETA[idbeta]<-sign[idbeta]*num[idbeta]/denum[idbeta]
  
  NEW.BETA.all[fix==0]<-NEWBETA
  
  
  b<-c(s,NEW.BETA.all)
  
  res<-gaussDYNidmlLikelihoodpena(b=b,
                                  npm=length(b),
                                  npar=length(b),
                                  bfix=1,
                                  fix=rep(0,length(b)),
                                  zi01=zi01,
                                  zi02=zi02,
                                  zi12=zi12,
                                  ctime=ctime,
                                  no=no,
                                  nz01=nz01,
                                  nz02=nz02,
                                  nz12=nz12,
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
                                  penalty.weights=penalty.weights0)
  
  return(list(res=res,b=b))
}
