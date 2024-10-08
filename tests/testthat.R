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
  # right censored data
  fitweib <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2,
               formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
               formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d)
  summary(fitweib)
  int<-intensity(times=fitweib$time[,1],
            theta = fitweib$modelPar[1:2]^2)
  predict(fitweib,s=10,t=15)
  
})