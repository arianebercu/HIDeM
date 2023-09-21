### test-idm.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Oct 22 2015 (13:57) 
## Version: 
## last-updated: Feb 25 2016 (08:43) 
##           By: Thomas Alexander Gerds
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(SmoothHazardoptim9)
#library(testthat)
#library(epiDisplay)
#library(doParallel)

data(Paq1000)
head(Paq1000)

set.seed(17)
d <- simulateIDM(100)$data
# right censored data
fitRC <- idm(formula01=Hist(time=observed.illtime,event=seen.ill)~X1+X2,
             formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
             formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,data=d
             )
fitRC
