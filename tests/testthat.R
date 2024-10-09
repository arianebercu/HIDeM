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
  summary(fitweib)
  int<-intensity(times=fitweib$time[,1],
            theta = fitweib$modelPar[1:2]^2)
  int
  predict(fitweib,s=10,t=15)
  
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
  summary(fitsplines)
  int<-intensity(times=fitsplines$time[,1],
                 theta = fitsplines$theta01^2,
                 knots = fitsplines$knots01,
                 number.knots = fitsplines$nknots01,
                 method = "splines")
  int
  predict(fitsplines,s=10,t=15)
  
})
#here
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
                 formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                 data=d,penalty="lasso",lambda01 = 10,lambda02 = 10, lambda12 = 10)
  summary(fitpenweib)
  int<-intensity(times=fitpenweib$time[,1],
                 theta = fitpenweib$modelPar[1:2]^2)
  int
  predict(fitpenweib,s=10,t=15)
  
  
  fitpenweib <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                    formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                    formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                    data=d,penalty="lasso",lambda01 = c(10,100),lambda02 = 10, lambda12 = 10)
  summary(fitpenweib)
  int<-intensity(times=fitpenweib$time[,1],
                 theta = fitpenweib$modelPar[1:2]^2)
  int
  plot(fitpenweib,lambda="BIC")
  predict(fitpenweib,s=10,t=15,lambda = c(10,10,100))
  
})


test_that("Illness-death model with M-spline baseline risk", {
  library(lava)
  library(prodlim)
  set.seed(17)
  d <- simulateIDM(n=1000)$data
  # right censored data
  fitsplines <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2,
                    formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
                    formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d,
                    method="splines")
  summary(fitsplines)
  int<-intensity(times=fitsplines$time[,1],
                 theta = fitsplines$theta01^2,
                 knots = fitsplines$knots01,
                 number.knots = fitsplines$nknots01,
                 method = "splines")
  int
  predict(fitsplines,s=10,t=15)
  
})