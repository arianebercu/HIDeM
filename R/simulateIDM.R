##' Function to generate an illness-death model for simulation.
##'
##' Based on the functionality of the lava PACKAGE the function generates
##' a latent variable model (latent illtime, waittime and lifetime)
##' and censoring mechanism (censtime, inspection1,inspection2,...,inspectionK).
##' 
##' The function \code{\link{sim.idmModel}} then simulates
##' right censored lifetimes and interval censored illness times.
##'
##' @title Generate illness-death model objects
##' @aliases idmModel
##' @param scale.illtime Weilbull scale for latent illness time
##' @param shape.illtime Weilbull shape for latent illness time
##' @param scale.lifetime Weilbull scale for latent life time
##' @param shape.lifetime Weilbull shape for latent life time
##' @param scale.waittime Weilbull scale for latent life time
##' @param shape.waittime Weilbull shape for latent life time
##' @param scale.censtime Weilbull scale for censoring time
##' @param shape.censtime Weilbull shape for censoring time
##' @param n.inspections Number of inspection times
##' @param schedule Mean of the waiting time between adjacent
##' inspections.
##' @param punctuality Standard deviation of waiting time between
##' inspections.
##' @examples
##' \dontrun{
##' library(lava)
##' library(prodlim)
##' # generate illness-death model based on exponentially
##' # distributed times
##' m <- idmModel(scale.illtime=1/70,
##'               shape.illtime=1.8,
##'               scale.lifetime=1/50,
##'               shape.lifetime=0.7,
##'               scale.waittime=1/30,
##'               shape.waittime=0.7)
##' round(sim(m,6),1)
##' 
##' # Estimate the parameters of the Weibull models
##' # based on the uncensored exact event times
##' # and the uncensored illstatus.
##' set.seed(18)
##' d <- sim(m,100,latent=FALSE)
##' d$uncensored.status <- 1
##' f <- idm(formula01=Hist(time=illtime,event=illstatus)~1,
##'          formula02=Hist(time=lifetime,event=uncensored.status)~1,
##'          data=d,
##'          conf.int=FALSE)
##' print(f)
##' 
##' # Change the rate of the 0->2 and 0->1 transitions
##' # also the rate of the 1->2 transition
##' # and also lower the censoring rate
##' m <- idmModel(scale.lifetime=1/2000,
##'               scale.waittime=1/30,
##'               scale.illtime=1/1000,
##'               scale.censtime=1/1000)
##' set.seed(18)
##' d <- sim(m,50,latent=TRUE)
##' d$uncensored.status <- 1
##' 
##' f <- idm(formula01=Hist(time=observed.illtime,event=illstatus)~1,
##'          formula02=Hist(time=observed.lifetime,event=uncensored.status)~1,
##'          data=d,
##'          conf.int=FALSE)
##' print(f)
##' 
##' # Estimate based on the right censored observations
##' fc <- idm(formula01=Hist(time=illtime,event=seen.ill)~1,
##'           formula02=Hist(time=observed.lifetime,event=seen.exit)~1,
##'           data=d,
##'           conf.int=FALSE)
##' print(fc)
##' 
##' # Estimate based on interval censored and right censored observations
##' fi <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~1,
##'           formula02=Hist(time=observed.lifetime,event=seen.exit)~1,
##'           data=d,
##'           conf.int=FALSE)
##' print(fi)
##' 
##' # Estimation of covariate effects:
##' # X1, X2, X3
##' m <- idmModel(shape.waittime=2,
##'               scale.lifetime=1/2000,
##'               scale.waittime=1/300,
##'               scale.illtime=1/10000,
##'               scale.censtime=1/10000)
##' distribution(m,"X1") <- binomial.lvm(p=0.3)
##' distribution(m,"X2") <- normal.lvm(mean=120,sd=20)
##' distribution(m,"X3") <- normal.lvm(mean=50,sd=20)
##' regression(m,to="latent.illtime",from="X1") <- 1.7
##' regression(m,to="latent.illtime",from="X2") <- 0.07
##' regression(m,to="latent.illtime",from="X3") <- -0.1
##' regression(m,to="latent.waittime",from="X1") <- 1.8
##' regression(m,to="latent.lifetime",from="X1") <- 0.7
##' set.seed(28)
##' d <- sim(m,100,latent=TRUE)
##' head(d)
##' table(ill=d$seen.ill,death=d$seen.exit)
##' 
##' # Estimation based on uncensored data
##' d$uncensored.status <- 1
##' # uncensored data
##' F1 <- idm(formula01=Hist(time=illtime,event=illstatus)~X1+X2+X3,
##'           formula02=Hist(time=lifetime,event=uncensored.status)~X1+X2+X3,
##'           data=d,conf.int=FALSE)
##' print(F1)
##' 
##' # Estimation based on right censored data
##' F2 <- idm(formula01=Hist(time=illtime,event=seen.ill)~X1+X2+X3,
##'           formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2+X3,
##'           data=d,conf.int=FALSE)
##' print(F2)
##' 
##' # Estimation based on interval censored and right censored data
##' F3 <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2+X3,
##'           formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2+X3,
##'           data=d,conf.int=FALSE)
##' print(F3)
##' cbind(uncensored=F1$coef,right.censored=F2$coef,interval.censored=F3$coef)
##' }
##' @return A latent variable model object \code{lvm}
##' @author Thomas Alexander Gerds
##' @useDynLib SmoothHazardoptim9
##' @export
idmModel <- function(scale.illtime=1/100,
                     shape.illtime=1,
                     scale.lifetime=1/100,
                     shape.lifetime=1,
                     scale.waittime=1/100,
                     shape.waittime=1,
                     scale.censtime=1/100,
                     shape.censtime=1,
                     n.inspections=5,
                     schedule=10,
                     punctuality=5){

    ## illness-death-model
    ##
    ## model waiting time in state 0
    ## based on latent "ill" and "life" times
    ## separate model for the waiting time in state "ill"
    idm <- lava::lvm()
    lava::latent(idm) <- c("latent.illtime","latent.lifetime","latent.waittime","censtime")
    lava::distribution(idm,"latent.illtime") <- lava::coxWeibull.lvm(scale=scale.illtime,
                                                                     shape=shape.illtime)
    lava::distribution(idm,"latent.lifetime") <- lava::coxWeibull.lvm(scale=scale.lifetime,
                                                                      shape=shape.lifetime)
    ## idm <- lava::eventTime(idm,illtime~min(latent.illtime=1,latent.lifetime=0),"ill")
    ## idm <- lava::eventTime(idm,lifetime~min(illtime=1,latent.lifetime=2),"event")
    lava::distribution(idm,"latent.waittime") <- lava::coxWeibull.lvm(scale=scale.waittime,
                                                                      shape=shape.waittime)
    ## transform(m,lifetime~latent.lifetime+latent.illtime+ill+latent.waittime) <- function(x){
    ## if (x$ill==1)
    ## x$latent.illtime+x$latent.waittime
    ## else
    ## x$latent.lifetime
    ## }
    ## ===============================================
    ## right censoring time
    ## ===============================================
    lava::distribution(idm,"censtime") <- lava::coxWeibull.lvm(shape=shape.censtime,
                                                               scale=scale.censtime)
    if (n.inspections>0)
        for (k in 1:n.inspections){
            lava::distribution(idm,paste("inspection",k,sep="")) <- lava::normal.lvm(mean=schedule,sd=punctuality)
        }
    ## }
    class(idm) <- c("idmModel","lvm")
    ## class(idm) <- c("lvm")
    idm
}
##' Function to simulate illness-death model data
##'
##' Based on the functionality of the lava PACKAGE 
##' @title Simulate illness-death model data
##' @param x An \code{idmModel} object as obtained with
##' \code{idmModel}
##' @param n Number of observations
##' @param illness.known.at.death Affects the value of variable
##' seen.ill
##' @param compliance Probability of missing an inspection time.
##' @param latent if TRUE keep the latent event times
##' @param keep.inspectiontimes if \code{TRUE} keep the inspection
##' times.
##' @param ... Extra arguments given to \code{sim}
##' @return A data set with interval censored observations from an illness-death model
##' @examples
##' example(idmModel)
##' help(idmModel)
##' @author Thomas Alexander Gerds
##' @importFrom lava sim
##' @useDynLib SmoothHazardoptim9
##' @export
sim.idmModel <- function(x,
                         n,
                         plot,
                         pen,
                         illness.known.at.death=TRUE,
                         compliance=1,
                         latent=FALSE,
                         keep.inspectiontimes=FALSE,
                         ...){
    # simulate latent data
    #class(x) <- "lvm"
    #dat <- lava::sim(x,n=n,...)
    # construct illtime and true illness status

  if(pen==T){
    dat<-x
    T01<-dat$latent.illtime
    T02<-dat$latent.lifetime
    dat$illtime <- dat$latent.illtime
    dat$illstatus <- 1*(dat$illtime<=dat$latent.lifetime)
    dat$illtime[dat$illtime>dat$latent.lifetime] <- 0
    # construct lifetime
    # for ill subjects as the sum of the time to illness (illtime) and
    # the time spent in the illness state (waittime)
    dat$lifetime <- dat$latent.lifetime
    dat$lifetime[dat$illstatus==1] <- dat$illtime[dat$illstatus==1]+dat$latent.waittime[dat$illstatus==1]
    # interval censored illtime
    ipos <- grep("inspection[0-9]+",names(dat))
  }else{
    class(x) <- "lvm"
    dat <- lava::sim(x,n=n,...)
    # construct illtime and true illness status
    dat$illtime <- dat$latent.illtime
    dat$illstatus <- 1*(dat$illtime<=dat$latent.lifetime)
    dat$illtime[dat$illtime>dat$latent.lifetime] <- 0
    # construct lifetime
    # for ill subjects as the sum of the time to illness (illtime) and
    # the time spent in the illness state (waittime)
    dat$lifetime <- dat$latent.lifetime
    dat$lifetime[dat$illstatus==1] <- dat$illtime[dat$illstatus==1]+dat$latent.waittime[dat$illstatus==1]
    # interval censored illtime
    ipos <- grep("inspection[0-9]+",names(dat))
  }
    if (length(ipos)>0) {
        # compute inspection times
        # make sure all inspection times are in the future
        # of the previous inspection time
        iframe <- dat[,ipos]
        dat <- dat[,-ipos]
        iframe <- do.call("rbind",lapply(1:n,function(i){
            cumsum(pmax(unlist(iframe[i,]),0))
        }))
        
        interval <- do.call("rbind",lapply(1:n,function(i){
            ## remove duplicates
            itimes <- unique(iframe[i,])
            ## remove inspections where compliance is
            ## sampled as zero
            if (compliance<1 & compliance>0){
                comp <- rbinom(length(itimes),1,compliance)
                itimes <- itimes[comp==1]
            }
            ## remove inspection times that are 
            ## larger than the individual lifetime
            itimes <- itimes[itimes<dat$lifetime[i]]
            ## and those larger than the right censoring time
            itimes <- itimes[itimes<dat$censtime[i]]
            ## if all inspection times are censored
            ## set a single one at 0
            if (length(itimes)==0) itimes <- 0
            ## mark the last inspection time 
            last.inspection <- itimes[length(itimes)]
            ## find the interval where illness happens
            if (dat$illstatus[i]){
                ## subject was ill
                if (dat$illtime[i] > last.inspection){
                    ## no illness observed at last inspection
                    if (dat$censtime[i]<dat$lifetime[i]){
                        ## right censored: no illness observed
                        c(last.inspection,dat$censtime[i],0)
                    }else{
                        ## user option decides if illness is recorded at death
                        c(last.inspection,dat$lifetime[i],illness.known.at.death)
                    }
                }else{ ## illtime is smaller or equal to last inspection time
                    if (length(itimes)==1){
                        c(0,itimes,1)
                    } else{
                        hit <- prodlim::sindex(eval.times=dat$illtime[i],
                                               jump.times=itimes,
                                               strict=TRUE)
                        c(c(0,itimes)[c(1+hit,2+hit)],1)
                    }
                }
            } else {
                ## subject was never ill
                if (dat$censtime[i]<dat$lifetime[i]){
                    ## right censored: no illness observed
                    ## until last inspection
                    c(last.inspection,dat$censtime[i],0)
                } else{
                    ## no illness observed until death
                    if(illness.known.at.death==1){
                        c(dat$lifetime[i],dat$lifetime[i],0)
                    }else{
                        ## no illness observed
                        ## until last inspection
                        ## provide [L=last insp.; R=lifetime]
                        c(last.inspection,dat$lifetime[i],0)
                    }
                }
            }
        }))
        colnames(interval) <- c("L","R","seen.ill")
        dat <- cbind(dat,interval)
        if (latent==FALSE)
            dat <- dat[,-grep("latent\\.",names(dat))]
        if (keep.inspectiontimes) dat <- cbind(dat,iframe)
    }
    dat$seen.exit <- 1*(dat$lifetime<dat$censtime)
    dat$observed.lifetime <- pmin(dat$lifetime,dat$censtime)
    dat$observed.illtime <- pmin(dat$illtime,dat$censtime)
    dat$observed.illtime[dat$illstatus==0] <- -9
    dat$lifetime <- x$latent.lifetime
    dat$illtime[dat$illstatus==0] <- -9
    if(pen==T){
    dat$T01<-T01
    dat$T02<-T02}
    return(list(data=dat,
                plot=plot))
}


##' Simulate data from an illness-death model with interval censored event times and covariates 
##'
##' Simulate data from an illness-death model with interval censored event times
##' and covariates for the purpose of illustrating the help pages of the SmoothHazard package.
##' See the body of the function for details, i.e., evaluate simulateIDM
##' @seealso idmModel sim.idmModel
##' @title Sample illness-death model data
##' @examples
##' # simulateIDM
##' simulateIDM(100)
##' @useDynLib SmoothHazardoptim9
#' @export
#' @param n number of observations
simulateIDM <- function(n=100){
    m <- idmModel(shape.waittime=2,
                  scale.lifetime=1/2000,
                  scale.waittime=1/3000,
                  scale.illtime=1/7000,
                  scale.censtime=1/1000,
                  n.inspections=5,schedule=10,punctuality=5)
    lava::distribution(m,"X1") <- lava::binomial.lvm(p=0.3)
    lava::distribution(m,"X2") <- lava::normal.lvm(mean=120,sd=20)
    lava::distribution(m,"X3") <- lava::normal.lvm(mean=50,sd=20)
    lava::regression(m,to="latent.illtime",from="X1") <- 0.7
    lava::regression(m,to="latent.illtime",from="X2") <- 0.07
    lava::regression(m,to="latent.illtime",from="X3") <- -0.1
    lava::regression(m,to="latent.waittime",from="X1") <- 0.8
    lava::regression(m,to="latent.lifetime",from="X1") <- 0.7
    
    sim.idmModel(m,n,plot=NULL,pen = F)
}



##' Simulate data from an illness-death model with interval censored event times and penalized covariates 
##'
##' Simulate data from an illness-death model with interval censored event times
##' and covariates for the purpose of illustrating the help pages of the SmoothHazard package.
##' See the body of the function for details, i.e., evaluate simulateIDM
##' @seealso idmModel sim.idmModel simulateIDM
##' @title Sample illness-death model data
##' @examples
##' # simulatepenIDM
##' simulatepenIDM(100)
##' @useDynLib SmoothHazardoptim9
#' @export
#' @param n number of observations

simulatepenIDM <- function(n=100,seed,scale.illtime,shape.illtime,
                           scale.lifetime,shape.lifetime,
                           scale.waittime,shape.waittime,
                           #scale.censtime,shape.censtime,
                           prob.censoring,administrative.censoring,
                           n.inspections,
                           schedule,punctuality,nvar,mean,sd,cov,
                           x01,x02,x12,beta01,beta02,beta12){

  # Set the seed for reproducibility
  set.seed(seed)
  
  # Define the number of latent processes
  # "latent.illtime","latent.lifetime","latent.waittime","censtime"
  num_latent_processes <- 4
  
  # Define the number of exogenous variables
  # nvar and inspections
  num_exogenous_vars <- nvar+n.inspections
  
  # Generate random data for exogenous variables X1 to XN
  #exogenous_data <- as.data.frame(matrix(rnorm(nvar * n, mean = mean, sd = sd), ncol = nvar))
  exogenous_data<-mvrnorm(n = n, mu = mean, Sigma = cov)
  # Set the column names of the exogenous data frame
  #colnames(exogenous_data) <- paste0("X", 1:nvar)
  
  # Create a covariance matrix for the exogenous variables
  #cov(exogenous_data) <- cov
  #exogenous_data<-cbind(exogenous_data,as.data.frame(matrix(rnorm(nvar * n, mean = schedule, sd = punctuality), ncol = n.inspections)))
  #cov<-matrix(0,n.inspections,n.inspections)
  #diag(cov)<-punctuality
  #exogenous_data<-cbind(exogenous_data,mvrnorm(n = n, mu = rep(schedule,n.inspections), Sigma = cov))

  V<-matrix(NA,nrow=n,ncol=n.inspections)
  N_max<-n
  still<-c(1:N_max)
  C<-rep(NA,n)
  
   for(j in 1:n.inspections){
    if(j==1){
      V[,j]<-0
    }else{
      #V[,j]<-runif(n=N[i],min=(j-1)*step[i],max=(step[i]*(j-1)+var.step[i]))}
      V[,j]<-stats::runif(n=n,min=(j-1)*schedule,max=(schedule*(j-1)+punctuality))}
     
     if(j>1){
       id.censoring<-stats::rbinom(N_max,1,prob.censoring)
       C[still]<-ifelse(id.censoring==1,V[still,j-1],NA)
       N_max<-N_max-sum(id.censoring)
       still<-still[id.censoring==0]
       
     }
   }
  # C at NA put last visit 
  C[is.na(C)]<-V[is.na(C),dim(V)[2]]
  
  exogenous_data<-as.data.frame(cbind(exogenous_data,V))
  colnames(exogenous_data) <- c(paste0("X", 1:nvar),paste0("inspection",1:n.inspections))
  #cov(exogenous_data)<-matrix(0,num_exogenous_vars ,num_exogenous_vars )
  #cov(exogenous_data)[1:nvar,1:nvar] <- cov
  #diag(cov(exogenous_data)[(nvar+1):num_exogenous_vars,(nvar+1):num_exogenous_vars]) <- punctuality
  # Define the structural model for latent processes following a Weibull distribution
  U<-stats::runif(n=n)
  X01<-as.matrix(exogenous_data[,colnames(exogenous_data)%in%x01])
  X02<-as.matrix(exogenous_data[,colnames(exogenous_data)%in%x02])
  X12<-as.matrix(exogenous_data[,colnames(exogenous_data)%in%x12])
  
  latent.illtime<-((-log(1-U)*exp(-X01%*%beta01))^(1/shape.illtime))/scale.illtime
  latent.lifetime<-((-log(1-U)*exp(-X02%*%beta02))^(1/shape.lifetime))/scale.lifetime
  latent.waittime<-((-log(1-U)*exp(-X12%*%beta12))^(1/shape.waittime))/scale.waittime
  #censtime<-((-log(1-U))^(shape.censtime))/shape.censtime
  
  censtime<-C
  censtime<-pmin(censtime,administrative.censoring)
  
  fit <- data.frame(latent.illtime=latent.illtime,
                    latent.lifetime=latent.lifetime,
                    latent.waittime=latent.waittime,
                    censtime=censtime)
  fit<-cbind(fit,exogenous_data)
  
  time<-c(1:max(V))
  S01<-exp(-(scale.illtime*time)^shape.illtime)
  S02<-exp(-(scale.lifetime*time)^shape.lifetime)
  S12<-exp(-(scale.waittime*time)^shape.waittime)
  #cens<-exp(-(scale.censtime*time)^shape.censtime)
  
  data.weibull<-base::data.frame(survie=c(S01,S02,S12),
                           type=c(rep("01",length(S01)),
                                  rep("02",length(S02)),
                                  rep("12",length(S12))),
                           time=rep(time,3))
  
  p2<-ggplot2::ggplot(data=data.weibull,aes(y=survie,x=time,color=type))+ggplot2::geom_point()+ggplot2::geom_line()+
    ggplot2::theme_classic()+ggplot2::ylab("Survival")
  
  
  sim.idmModel(x=fit,n=n,plot=p2,pen=T,illness.known.at.death=F)
  
}


