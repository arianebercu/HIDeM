# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(SmoothHazardoptim9)

test_check("SmoothHazardoptim9")

test_that("Illness-death model with weibull baseline risk", {
  library(lava)
  library(prodlim)
  set.seed(17)
  d <- simulateIDM(n=1000)$data

  fitweib <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2,
               formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
               formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d)
  fitweib
  summary(fitweib)
  int<-intensity(times=fitweib$time[,1],
            theta = fitweib$modelPar[1:2]^2)
  int
  predict(fitweib,s=10,t=15)
  plot(fitweib)
  
})

test_that("Illness-death model with M-spline baseline risk", {
  library(lava)
  library(prodlim)
  set.seed(17)
  d <- simulateIDM(n=1000)$data

  fitsplines <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2,
                 formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
                 formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d,
                 method="splines")
  fitsplines
  summary(fitsplines)
  int<-intensity(times=fitsplines$time[,1],
                 theta = fitsplines$theta01^2,
                 knots = fitsplines$knots01,
                 number.knots = fitsplines$nknots01,
                 method = "splines")
  int
  predict(fitsplines,s=10,t=15)
  plot(fitsplines)
  
})

test_that("Penalised illness-death model with weibull baseline risk", {
  library(lava)
  library(prodlim)
  set.seed(17)
  d <- simulateIDM(n=1000,
                   beta01=c(1,1,0,0.5,0.5,rep(0,5)),
                   beta02=c(1,0,0,0,0.5,rep(0,5)),
                   beta12=c(1,0,0,0,0.5,rep(0,5)))$data
  

  
  fitpenweib <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                    formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                    formula12=~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                    data=d,penalty="lasso",lambda01 = c(10,20),lambda02 = 10, lambda12 = 10)
  fitpenweib
  summary(fitpenweib)
  int<-intensity(times=fitpenweib$time[,1],
                 theta = fitpenweib$modelPar[1:2]^2)
  int
  plot(fitpenweib,lambda=c(10,10,10))
  predict(fitpenweib,s=10,t=15,lambda = "BIC")
  
  
  fitpenweibunique <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                          formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                          formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                          data=d,penalty="lasso",lambda01 = 10,lambda02 = 10, lambda12 = 10)
  fitpenweibunique
  summary(fitpenweibunique)
  int<-intensity(times=fitpenweibunique$time[,1],
                 theta = fitpenweibunique$modelPar[1:2]^2)
  int
  predict(fitpenweibunique,s=10,t=15)
  plot(fitpenweibunique)
  
  
})

#here 
test_that("Penalised illness-death model with M-splines baseline risk", {
  library(lava)
  library(prodlim)
  set.seed(17)
  d <- simulateIDM(n=1000,
                   beta01=c(1,1,0,0.5,0.5,rep(0,5)),
                   beta02=c(1,0,0,0,0.5,rep(0,5)),
                   beta12=c(1,0,0,0,0.5,rep(0,5)))$data
  
  
  
  fitpenspline <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                    formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                    formula12=~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,method="splines",
                    data=d,penalty="lasso",lambda01 = c(10,20),lambda02 = 10, lambda12 = 10)
  fitpenspline
  summary(fitpenspline)
  int<-intensity(times=fitpenspline$time[,1],
                 theta = fitpenspline$theta01[,1]^2,method="splines",
                 knots=fitpenspline$knots01,number.knots=fitpenspline$nknots01)
  int
  plot(fitpenspline,lambda=c(10,10,10))
  pp<-predict(fitpenspline,s=10,t=15,lambda = "BIC")
  pp
  
  
  fitpensplineunique <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                          formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                          formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                          data=d,penalty="lasso",lambda01 = 10,lambda02 = 10, lambda12 = 10,method="splines")
  fitpensplineunique
  summary(fitpensplineunique)
  int<-intensity(times=fitpensplineunique$time[,1],
                 theta = fitpensplineunique$theta01[,1]^2,method="splines",
                 knots=fitpensplineunique$knots01,number.knots=fitpensplineunique$nknots01)
  int
  plot(fitpensplineunique)
  pp<-predict(fitpensplineunique,s=10,t=15,lambda = rep(10,3))
  pp
  
  
})