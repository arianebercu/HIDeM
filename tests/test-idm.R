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
library(SmoothHazard)
library(testthat)
library(epiDisplay)
library(doParallel)

data(Paq1000)
head(Paq1000)

load("~/Test_smooth_hazard/fit.idm.quant4.type.equi.RData")
fit.idm.nopena.quant.4.type.equi$cv
fit.idm.nopena.quant.4.type.equi$loglik
fit.idm.nopena.quant.4.type.equi$call

Paq1000$l_delai<-Paq1000$l-Paq1000$e
Paq1000$r_delai<-Paq1000$r-Paq1000$e
Paq1000$t_delai<-Paq1000$t-Paq1000$e

subPaq1000<-Paq1000[!(Paq1000$dementia==0 & Paq1000$death==1),]

subPaq1000[395,]
Paq1000[1000,]

Paq1000[Paq1000$death==1,]

################################################################################
#################### NO TWO STEP ###############################################
################################################################################


###################### 8 #######################################################

fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       knots="quantile",eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",type.quantile=3,twostep.spline = F,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.type.quantile3.3.no.step.RData")




fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",twostep.spline = F,maxiter=50)
save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.equi.3.no.step.RData")









###################### 6 #######################################################


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       knots="quantile",eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",type.quantile=3,twostep.spline = F,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.type.quantile3.3.no.step.RData")




fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",twostep.spline = F,maxiter=50)
save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.equi.3.no.step.RData")




###################### 4 #######################################################

Rprof()
fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       knots="quantile",eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",type.quantile=3,twostep.spline = F,maxiter=50)
Rprof(NULL)
summaryRprof()

fit.idm.paquid$runtime
fit.idm.paquid$runtime.beforemla
fit.idm.paquid$runtime.beforemla*100/fit.idm.paquid$runtime
fit.idm.paquid$runtime.aftermla
fit.idm.paquid$runtime.aftermla*100/fit.idm.paquid$runtime
fit.idm.paquid$runtime.intensity
fit.idm.paquid$runtime.intensity*100/fit.idm.paquid$runtime

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.type.quantile3.3.no.step.RData")




fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",twostep.spline = F,maxiter=50)
save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.equi.3.no.step.RData")





####################### FORTRAN ################################################
####################### ATTENTION INTENSITY IC ARE WRONGED #####################

########################## 8 ##################################################
fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       knots="quantiles",eps=c(3,3,2),method="sparse.splines",
                       type.quantile=3,twostep.spline = F,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.type.quantile3.3.no.step.FORTRAN.RData")



fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       eps=c(3,3,2),method="sparse.splines",
                       twostep.spline = F,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.type.equi.3.no.step.FORTRAN.RData")


########################## 6 ##################################################


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       knots="quantiles",eps=c(3,3,2),method="sparse.splines",
                       type.quantile=3,twostep.spline = F,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.type.quantile3.3.no.step.FORTRAN.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       eps=c(3,3,2),method="sparse.splines",
                       twostep.spline = F,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.type.equi.3.no.step.FORTRAN.RData")



########################## 4 ##################################################
fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       knots="quantiles",eps=c(3,3,2),method="sparse.splines",
                       type.quantile=3,twostep.spline = F,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.type.quantile3.3.no.step.FORTRAN.RData")



fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       eps=c(3,3,2),method="sparse.splines",
                       twostep.spline = F,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.type.equi.3.no.step.FORTRAN.RData")


################################################################################
####################  TWO STEP #################################################
################################################################################


###################### 8 #######################################################

fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       knots="quantile",eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",type.quantile=3,step.sequential = T,maxiter=50,option.sequential=list(step=5,
                                                                                                       cutoff.spline = 0.01,
                                                                                                       min=20,
                                                                                                       value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.type.quantile3.3.step.5.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       knots="quantile",eps=c(3,3,2),type.quantile=3,step.sequential = T,maxiter=50,option.sequential=list(step=5,
                                                                                                                   cutoff.spline = 0.01,
                                                                                                                   min=20,
                                                                                                                   value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.type.quantile3.3.step.5.nproc1.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       knots="quantile",eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",type.quantile=3,step.sequential = T,maxiter=50,option.sequential=list(step=10,
                                                                                                                   cutoff.spline = 0.01,
                                                                                                                   min=20,
                                                                                                                   value.zero=F))



save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.type.quantile3.3.step.10.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       knots="quantile",eps=c(3,3,2),type.quantile=3,step.sequential = T,maxiter=50,option.sequential=list(step=10,
                                                                                                                           cutoff.spline = 0.01,
                                                                                                                           min=20,
                                                                                                                           value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.type.quantile3.3.step.10.nproc1.RData")



fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",step.sequential = T,maxiter=50,option.sequential=list(step=5,
                                                                                                   cutoff.spline = 0.01,
                                                                                                   min=20,
                                                                                                   value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.equi.3.step.5.RData")

fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       eps=c(3,3,2),step.sequential = T,maxiter=50,option.sequential=list(step=5,
                                                                                                   cutoff.spline = 0.01,
                                                                                                   min=20,
                                                                                                   value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.equi.3.step.5.nproc1.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",step.sequential = T,maxiter=50,option.sequential=list(step=10,
                                                                                                   cutoff.spline = 0.01,
                                                                                                   min=20,
                                                                                                   value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.equi.3.step.10.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       eps=c(3,3,2),step.sequential = T,maxiter=50,option.sequential=list(step=10,
                                                                                                   cutoff.spline = 0.01,
                                                                                                   min=20,
                                                                                                   value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.equi.3.step.10.nproc1.RData")





###################### 6 #######################################################

fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       knots="quantile",eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",type.quantile=3,step.sequential = T,maxiter=50,option.sequential=list(step=5,
                                                                                                                   cutoff.spline = 0.01,
                                                                                                                   min=20,
                                                                                                                   value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.type.quantile3.3.step.5.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       knots="quantile",eps=c(3,3,2),type.quantile=3,step.sequential = T,maxiter=50,option.sequential=list(step=5,
                                                                                                                           cutoff.spline = 0.01,
                                                                                                                           min=20,
                                                                                                                           value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.type.quantile3.3.step.5.nproc1.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       knots="quantile",eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",type.quantile=3,step.sequential = T,maxiter=50,option.sequential=list(step=10,
                                                                                                                   cutoff.spline = 0.01,
                                                                                                                   min=20,
                                                                                                                   value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.type.quantile3.3.step.10.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       knots="quantile",eps=c(3,3,2),type.quantile=3,step.sequential = T,maxiter=50,option.sequential=list(step=10,
                                                                                                                           cutoff.spline = 0.01,
                                                                                                                           min=20,
                                                                                                                           value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.type.quantile3.3.step.10.nproc1.RData")



fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",step.sequential = T,maxiter=50,option.sequential=list(step=5,
                                                                                                   cutoff.spline = 0.01,
                                                                                                   min=20,
                                                                                                   value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.equi.3.step.5.RData")

fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       eps=c(3,3,2),step.sequential = T,maxiter=50,option.sequential=list(step=5,
                                                                                          cutoff.spline = 0.01,
                                                                                          min=20,
                                                                                          value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.equi.3.step.5.nproc1.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",step.sequential = T,maxiter=50,option.sequential=list(step=10,
                                                                                                   cutoff.spline = 0.01,
                                                                                                   min=20,
                                                                                                   value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.equi.3.step.10.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       eps=c(3,3,2),step.sequential = T,maxiter=50,option.sequential=list(step=10,
                                                                                          cutoff.spline = 0.01,
                                                                                          min=20,
                                                                                          value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.equi.3.step.10.nproc1.RData")




###################### 4 #######################################################

fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       knots="quantile",eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",type.quantile=3,step.sequential = T,maxiter=50,option.sequential=list(step=5,
                                                                                                                   cutoff.spline = 0.01,
                                                                                                                   min=20,
                                                                                                                   value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.type.quantile3.3.step.5.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       knots="quantile",eps=c(3,3,2),type.quantile=3,step.sequential = T,maxiter=50,option.sequential=list(step=5,
                                                                                                                           cutoff.spline = 0.01,
                                                                                                                           min=20,
                                                                                                                           value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.type.quantile3.3.step.5.nproc1.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       knots="quantile",eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",type.quantile=3,step.sequential = T,maxiter=50,option.sequential=list(step=10,
                                                                                                                   cutoff.spline = 0.01,
                                                                                                                   min=20,
                                                                                                                   value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.type.quantile3.3.step.10.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       knots="quantile",eps=c(3,3,2),type.quantile=3,step.sequential = T,maxiter=50,option.sequential=list(step=10,
                                                                                                                           cutoff.spline = 0.01,
                                                                                                                           min=20,
                                                                                                                           value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.type.quantile3.3.step.10.nproc1.RData")



fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",step.sequential = T,maxiter=50,option.sequential=list(step=5,
                                                                                                   cutoff.spline = 0.01,
                                                                                                   min=20,
                                                                                                   value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.equi.3.step.5.RData")

fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       eps=c(3,3,2),step.sequential = T,maxiter=50,option.sequential=list(step=5,
                                                                                          cutoff.spline = 0.01,
                                                                                          min=20,
                                                                                          value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.equi.3.step.5.nproc1.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",step.sequential = T,maxiter=50,option.sequential=list(step=10,
                                                                                                   cutoff.spline = 0.01,
                                                                                                   min=20,
                                                                                                   value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.equi.3.step.10.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       eps=c(3,3,2),step.sequential = T,maxiter=50,option.sequential=list(step=10,
                                                                                          cutoff.spline = 0.01,
                                                                                          min=20,
                                                                                          value.zero=F))

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.equi.3.step.10.nproc1.RData")









####################### FORTRAN ################################################
####################### ATTENTION INTENSITY IC ARE WRONGED #####################

########################## 8 ##################################################
fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       knots="quantiles",eps=c(3,3,2),method="sparse.splines",
                       type.quantile=3,twostep.spline = T,step.spline = 5,cutoff.spline = 0.01,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.type.quantile3.3.step.5.FORTRAN.RData")



fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       knots="quantiles",eps=c(3,3,2),method="sparse.splines",
                       type.quantile=3,twostep.spline = T,step.spline = 10,cutoff.spline = 0.01,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.type.quantile3.3.step.10.FORTRAN.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       eps=c(3,3,2),method="sparse.splines",
                       twostep.spline = T,step.spline = 5,cutoff.spline = 0.01,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.type.equi.3.step.5.FORTRAN.RData")



fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       eps=c(3,3,2),method="sparse.splines",
                       twostep.spline = T,step.spline = 10,cutoff.spline = 0.01,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.type.equi.3.step.10.FORTRAN.RData")




########################## 6 ##################################################
fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       knots="quantiles",eps=c(3,3,2),method="sparse.splines",
                       type.quantile=3,twostep.spline = T,step.spline = 5,cutoff.spline = 0.01,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.type.quantile3.3.step.5.FORTRAN.RData")



fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       knots="quantiles",eps=c(3,3,2),method="sparse.splines",
                       type.quantile=3,twostep.spline = T,step.spline = 10,cutoff.spline = 0.01,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.type.quantile3.3.step.10.FORTRAN.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       eps=c(3,3,2),method="sparse.splines",
                       twostep.spline = T,step.spline = 5,cutoff.spline = 0.01,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.type.equi.3.step.5.FORTRAN.RData")



fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(6,3),
                       eps=c(3,3,2),method="sparse.splines",
                       twostep.spline = T,step.spline = 10,cutoff.spline = 0.01,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant6.type.equi.3.step.10.FORTRAN.RData")


########################## 4 ##################################################
fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       knots="quantiles",eps=c(3,3,2),method="sparse.splines",
                       type.quantile=3,twostep.spline = T,step.spline = 5,cutoff.spline = 0.01,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.type.quantile3.3.step.5.FORTRAN.RData")



fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       knots="quantiles",eps=c(3,3,2),method="sparse.splines",
                       type.quantile=3,twostep.spline = T,step.spline = 10,cutoff.spline = 0.01,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.type.quantile3.3.step.10.FORTRAN.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       eps=c(3,3,2),method="sparse.splines",
                       twostep.spline = T,step.spline = 5,cutoff.spline = 0.01,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.type.equi.3.step.5.FORTRAN.RData")



fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(4,3),
                       eps=c(3,3,2),method="sparse.splines",
                       twostep.spline = T,step.spline = 10,cutoff.spline = 0.01,maxiter=50)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant4.type.equi.3.step.10.FORTRAN.RData")



#################################################################################
######################## ANALYSER NITER.SPLINE ##################################
#################################################################################

fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       knots="quantile",eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",type.quantile=3,twostep.spline = T,step.spline = 5,cutoff.spline = 0.01,maxiter=20)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.type.quantile3.3.step.5.maxit20.RData")


fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       knots="quantile",eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",type.quantile=3,twostep.spline = T,step.spline = 1,cutoff.spline = 0.01,maxiter=20)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.type.quantile3.3.step.1.maxit20.RData")

fit.idm.paquid <-  idm(formula02 = Hist(time = t_delai, event = death) ~ certif,
                       formula01 = Hist(time = list(l_delai,r_delai), event = dementia) ~ certif,
                       formula12 = ~ certif, data = Paq1000,n.knots=rep(8,3),
                       knots="quantile",eps=c(3,3,2),nproc=10,
                       clustertype = "PSOCK",type.quantile=3,twostep.spline = T,step.spline = 20,cutoff.spline = 0.01,maxiter=20)

save(fit.idm.paquid ,file="C:/Users/ab17/Documents/Test_smooth_hazard/allR/fit.idm.quant8.type.quantile3.3.step.20.maxit20.RData")


