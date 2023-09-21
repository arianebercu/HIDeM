### test-idm.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Oct 22 2015 (13:57) 
## Version: 
## last-updated: Feb 25 2016 (08:43) 
##           By: Thomas Alexander Gerds
##     Update #: 9
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(SmoothHazardoptim9)
library(testthat)

testthat::test_that("idm weibull paquid data with covariates",{
    data(Paq1000)
    set.seed(17)
    paq <- Paq1000[sample(1:NROW(Paq1000),size=200),]
    fit.weib <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
                    formula01=Hist(time=list(l,r),event=dementia)~certif,
                    data=paq)
    fit.weib2 <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
                    formula01=Hist(time=list(l,r),event=dementia)~certif,
                    formula12 = ~ 1,
                    data=paq)
    pred <- predict(fit.weib,70,t=80,newdata=data.frame(certif=1), conf.int = TRUE, nsim=4)
    pred2 <- predict(fit.weib2,70,t=80,newdata=data.frame(certif=1), conf.int = TRUE, nsim=4,lifeExpect=TRUE)
    expect_output(pred)
    expect_output(pred2)
}) 


#----------------------------------------------------------------------
### test-idm.R ends here
