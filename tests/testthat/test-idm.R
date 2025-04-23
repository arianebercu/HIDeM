

test_that("HIDeM idm", {
  
data("mydata", package = "HIDeM") 

mydata.scale<-scale(mydata[,paste0("X",seq(1:50))])
mydata[,paste0("X",seq(1:50))]<-mydata.scale

k<-1
penalty<-ifelse(k==1,"lasso","elasticnet")

f02<-as.formula(paste0("Hist(time = observed.lifetime, event = seen.exit) ~",paste0(paste0("X",seq(1:50)),collapse="+")))
f01<-as.formula(paste0("Hist(time = list(L,R), event = seen.ill) ~",paste0(paste0("X",seq(1:50)),collapse="+")))
f12<-as.formula(paste0("~",paste0(paste0("X",seq(1:50)),collapse="+")))


model <-HIDeM::  idm(formula02 = f02,
                         formula01 = f01,
                         formula12 = f12,data=mydata,
                         penalty = penalty,eps=c(7,4,2),
                         method="Weib",scale.X=F,
                         lambda01 = 40,
                         lambda02=45,
                         lambda12=30,nproc=1,
                         alpha=k,maxiter=100,maxiter.pena = 10)

expect_true(exists("model"))


})



test_that("HIDeM grid search", {
  data("mydata", package = "HIDeM") 
  
  mydata.scale<-scale(mydata[,paste0("X",seq(1:50))])
  mydata[,paste0("X",seq(1:50))]<-mydata.scale
  
  ############################## 01 ##############################################
  
  k<-1
  f02<-as.formula("Hist(time = observed.lifetime, event = seen.exit) ~1")
  f01<-as.formula(paste0("Hist(time = list(L,R), event = seen.ill) ~",paste0(paste0("X",seq(1:50)),collapse="+")))
  f12<-as.formula("~1")
  
  #setwd("C:/Users/ab17/Desktop/CURTA/SMOOTH_HAZARD/OPTIM_B_DIFFLAMBDA_derivaana_all/Weibull_simu/Scenario1")
  
  nlambda<-20
  lambdatest<-seq(0.01,100,length.out=nlambda)
  penalty<-ifelse(k==1,"lasso","elasticnet")
  
  lambda01<-lambda02<-lambda12<-rep(NA,3)
  
  # fixed base risk for 02 - 12
  
  #load(paste0("solution_compare_model_",rep,"_mlaana_splines.RData"))
  
  m01 <-HIDeM::  idm(formula02 = f02,
                     formula01 = f01,
                     formula12 = f12,data=mydata,
                     nproc=1, penalty=penalty,eps=c(7,4,2),
                     method="Weib",scale.X=F,
                     lambda01 = lambdatest/k,
                     lambda02=0.0001,
                     lambda12=0.0001,
                     alpha=k,maxiter=100,maxiter.pena = 10)
  
  
  
  
  BIC01<-sort(m01$BIC[m01$converged==1])[1:3]
  BIC01<-na.omit(unique(BIC01))
  
  id<-which(m01$BIC%in%BIC01)
  if(length(id)>3){
    lambda<-m01$lambda[1,id]
    lambda01[1:3]<-sort(lambda)[1:3]
  }else{
    if(length(id)==0){
      BIC01<-sort(m01$BIC)[1:3]
      BIC01<-na.omit(unique(BIC01))
      lambda01<-m01$lambda[1,which(m01$BIC%in%BIC01)[1:3]]
      lambda01<-na.omit(lambda01)
    }else{
      lambda01[1:3]<-m01$lambda[1,id]
      lambda01<-na.omit(lambda01)}
  }
  
  
  ################################ 02 ############################################
  
  f02<-as.formula(paste0("Hist(time = observed.lifetime, event = seen.exit) ~",paste0(paste0("X",seq(1:50)),collapse="+")))
  f01<-as.formula("Hist(time = list(L,R), event = seen.ill) ~ 1")
  f12<-as.formula("~1")
  #setwd("C:/Users/ab17/Desktop/CURTA/SMOOTH_HAZARD/OPTIM_B_DIFFLAMBDA_derivaana_all/Weibull_simu/Scenario1")
  
  m02 <-HIDeM::  idm(formula02 = f02,
                     formula01 = f01,
                     formula12 = f12,data=mydata,
                     nproc=1, penalty=penalty,eps=c(7,4,2),
                     method="Weib",scale.X=F,
                     lambda01 = 0.0001,
                     lambda02= lambdatest/k,
                     lambda12=0.0001,
                     alpha=k,maxiter=100,maxiter.pena = 10)
  
  
  
  BIC02<-sort(m02$BIC[m02$converged==1])[1:3]
  BIC02<-na.omit(unique(BIC02))
  
  id<-which(m02$BIC%in%BIC02)
  if(length(id)>3){
    lambda<-m02$lambda[2,id]
    lambda02[1:3]<-sort(lambda)[1:3]
  }else{
    if(length(id)==0){
      BIC02<-sort(m02$BIC)[1:3]
      BIC02<-na.omit(unique(BIC02))
      lambda02<-m02$lambda[2,which(m02$BIC%in%BIC02)[1:3]]
      lambda02<-na.omit(lambda02)
    }else{
      lambda02[1:3]<-m02$lambda[2,id]
      lambda02<-na.omit(lambda02)}
  }
  
  
  
  ######################################## 12 ###################################
  
  f02<-as.formula("Hist(time = observed.lifetime, event = seen.exit) ~1")
  f01<-as.formula("Hist(time = list(L,R), event = seen.ill) ~ 1")
  f12<-as.formula(paste0("~",paste0(paste0("X",seq(1:50)),collapse="+")))
  
  
  m12 <-HIDeM::  idm(formula02 = f02,
                     formula01 = f01,
                     formula12 = f12,data=mydata,
                     nproc=1, penalty=penalty,eps=c(7,4,2),
                     method="Weib",scale.X=F,
                     lambda01 = 0.0001,
                     lambda12= lambdatest/k,
                     lambda02=0.0001,
                     alpha=k,maxiter=100,maxiter.pena = 10)
  
  
  BIC12<-sort(m12$BIC[m12$converged==1])[1:3]
  BIC12<-na.omit(unique(BIC12))
  
  
  
  id<-which(m12$BIC%in%BIC12)
  if(length(id)>3){
    lambda<-m12$lambda[3,id]
    lambda12[1:3]<-sort(lambda)[1:3]
  }else{
    if(length(id)==0){
      BIC12<-sort(m12$BIC)[1:3]
      BIC12<-na.omit(unique(BIC12))
      lambda12<-m12$lambda[3,which(m12$BIC%in%BIC12)[1:3]]
      lambda12<-na.omit(lambda12)
    }else{
      lambda12[1:3]<-m12$lambda[3,id]
      lambda12<-na.omit(lambda12)}
  }
  
  
  #setwd("C:/Users/ab17/Desktop/CURTA/SMOOTH_HAZARD/OPTIM_B_DIFFLAMBDA_derivaana")
  
  expect_true(exists("m12") & exists("m02") & exists("m01"))
  
  
})

