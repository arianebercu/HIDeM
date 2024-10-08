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
  d <- simulateIDM(n=1000)
  # right censored data
  fitRC <- idm(formula01=Hist(time=observed.illtime,event=seen.ill)~X1+X2,
               formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
               formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d,
               conf.int=FALSE)
  fitRC
  
})