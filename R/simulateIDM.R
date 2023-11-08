
##' Function to simulate illness-death model data
##'
##' Based on the functionality of the lava PACKAGE 
##' @title Simulate illness-death model data
##' @param x An \code{idmModel} object as obtained with
##' \code{idmModel}
##' @param n Number of observations
##' seen.ill
##' @param compliance Probability of missing an inspection time.
##' @param latent if TRUE keep the latent event times
##' @param keep.inspectiontimes if \code{TRUE} keep the inspection
##' times.
##' @param pen if the simulated data are for a penalised version of illness-death model
##' @param plot plot of base survival for all transition
##' @param ... Extra arguments given to \code{sim}
##' @return A data set with interval censored observations from an illness-death model
##' example(idmModel)
##' help(idmModel)
##' @author Thomas Alexander Gerds
##' @importFrom lava sim
##' @useDynLib SmoothHazardoptim9
sim.idmModel <- function(x,
                         n,
                         plot,
                         pen,
                         latent=FALSE,
                         keep.inspectiontimes=FALSE,
                         ...){
    # simulate latent data
    #class(x) <- "lvm"
    #dat <- lava::sim(x,n=n,...)
    # construct illtime and true illness status

    dat<-x
    T01<-dat$latent.illtime
    T02<-dat$latent.lifetime
    T12<-dat$latent.waittime
    dat$illtime <- dat$latent.illtime
    dat$illstatus <- 1*((dat$illtime<dat$latent.lifetime)& (dat$illtime<dat$censtime))
    #dat$illtime[dat$illtime>dat$latent.lifetime] <- 0
    # construct lifetime
    # for ill subjects as the sum of the time to illness (illtime) and
    # the time spent in the illness state (waittime)
    dat$lifetime <- dat$latent.lifetime
    dat$lifetime[dat$illstatus==1]<-dat$latent.waittime[dat$illstatus==1]
    id.nodem.death<-rep(0,n)
    
    

    # interval censored illtime
    ipos <- grep("inspection[0-9]+",names(dat))

    browser()
    if (length(ipos)>0) {
        # compute inspection times
        # make sure all inspection times are in the future
        # of the previous inspection time

        iframe <- dat[,ipos]
        dat <- dat[,-ipos]
        
        interval <- do.call("rbind",lapply(1:n,function(i){
          
            ## remove duplicates
            itimes <- unique(iframe[i,])
            
            ## remove inspection times that are 
            ## larger than the individual lifetime
            itimes <- itimes[itimes<dat$lifetime[i]]
            ## and those larger than the right censoring time
            itimes <- itimes[itimes<=dat$censtime[i]]
            ## if all inspection times are censored
            ## set a single one at 0
            #if (length(itimes)==0) {
            #  itimes <- 0}
            
            ## mark the last inspection time 
            #last.inspection <- itimes[length(itimes)]
            ## find the interval where illness happens
            
            if (dat$illstatus[i]){
              
              idL<-which(itimes<dat$illtime[i])
              idR<-which(itimes>=dat$illtime[i])
      
              # only one visit at inclusion
              if(length(idR)==0){
                R<-L<-itimes[max(idL)]
                if(dat$lifetime[i]>dat$administrative.censoring[i]){
                  c(L,R,dat$administrative.censoring[i],0,0)
                }else{
                  c(L,R,dat$lifetime[i],0,1)
                }
              }else{
                L<-itimes[max(idL)]
                R<-itimes[min(idR)]
                  ## subject was ill
                if (dat$lifetime[i]>dat$administrative.censoring[i]){
                  ## illness observed but not death
                     c(L,R,dat$administrative.censoring[i],1,0)
                 }else{
                   ## death observed all times
                   browser()
                   idVis<-which(itimes>dat$illtime[i])
                   # no visit between the two times : illness not observed
                   if(length(idVis)==0){
                     c(itimes[length(itimes)],itimes[length(itimes)],dat$lifetime[i],-1,1)
                      }else{ ## we observe illness 
                        c(L,R,dat$lifetime[i],1,1)
                      }
                 }
              }
              }else{
                
                # check administrative censoring for death
                if(dat$lifetime[i]<dat$administrative.censoring[i]){
                  
                  c(itimes[length(itimes)],itimes[length(itimes)],dat$lifetime[i],0,1)
                }else{
                  
                  c(itimes[length(itimes)],itimes[length(itimes)],dat$administrative.censoring[i],0,0)}
              }
          }))
        colnames(interval) <- c("L","R","observed.lifetime","seen.ill","seen.exit")
        # count illness not observed due to death 
        dat <- cbind(dat,interval)
        if (latent==FALSE)
            dat <- dat[,-grep("latent\\.",names(dat))]
        if (keep.inspectiontimes) dat <- cbind(dat,iframe)
    }
    browser()
    id.nodem.death[which(dat$seen.ill==-1)]<-1
    dat$seen.ill[dat$seen.ill==-1]<-0

    dat$T01<-T01
    dat$T02<-T02
    dat$T12<-T12
    dat$id.nodem.death<-id.nodem.death
    return(list(data=dat,
                plot=plot))
}





##' Simulate data from an illness-death model with interval censored event times and penalized covariates 
##'
##' Simulate data from an illness-death model with interval censored event times
##' and covariates for the purpose of illustrating the help pages of the SmoothHazard package.
##' See the body of the function for details, i.e., evaluate simulateIDM
##' @seealso idmModel sim.idmModel simulateIDM
##' @title Sample illness-death model data for penalised regression
##' @param scale.illtime Weilbull scale for latent illness time
##' @param shape.illtime Weilbull shape for latent illness time
##' @param scale.lifetime Weilbull scale for latent life time
##' @param shape.lifetime Weilbull shape for latent life time
##' @param scale.waittime Weilbull scale for latent life time
##' @param shape.waittime Weilbull shape for latent life time
##' @param seed specify the seed 
##' @param prob.censoring probability of censoring at each visit 
##' @param administrative.censoring specify time of administrative censoring
##' @param n.inspections Number of inspection times
##' @param schedule Mean of the waiting time between adjacent
##' inspections.
##' @param punctuality Standard deviation of waiting time between
##' inspections.
##' @param nvar number of variables
##' @param mean mean of each explanatory variables
##' @param sd standard-error of each explanatory variables
##' @param cov covariance matrix of explanatory variables
##' @param x01 names of variables on transition 0 --> 1
##' @param x02 names of variables on transition 0 --> 2
##' @param x12 names of variables on transition 1 --> 2
##' @param beta01 value of beta on transition 0 --> 1
##' @param beta02 value of beta on transition 0 --> 2
##' @param beta12 value of beta on transition 1 --> 2
##' @param n number of observations
#' @importFrom ggplot2 ggplot geom_line geom_point theme_classic ylab aes_string aes facet_grid
#' @useDynLib SmoothHazardoptim9
#' @export


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
  exogenous_data<-MASS::mvrnorm(n = n, mu = mean, Sigma = cov)
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
  U01<-stats::runif(n=n)
  U02<-stats::runif(n=n)
  U12<-stats::runif(n=n)
  # center and reduce var : 
  #exogenous_data<-scale(exogenous_data)
  
  X01<-as.matrix(exogenous_data[,colnames(exogenous_data)%in%x01])
  X02<-as.matrix(exogenous_data[,colnames(exogenous_data)%in%x02])
  X12<-as.matrix(exogenous_data[,colnames(exogenous_data)%in%x12])
  
  latent.illtime<-((-log(U01)*exp(-X01%*%beta01))^(1/shape.illtime))/scale.illtime
  latent.lifetime<-((-log(U02)*exp(-X02%*%beta02))^(1/shape.lifetime))/scale.lifetime
  latent.waittime<-rep(0,n)
  censtime<-C
  administrative.censoring<-rep(administrative.censoring,n)
  if(any(censtime>administrative.censoring)){stop("All visit should be performed before the administrative censoring")}
  
  cumulative.intensity<-(scale.waittime*latent.illtime)^shape.waittime
  e <- exp(X12%*%beta12)
  cumulative.intensity<-cumulative.intensity*e
  S12 <- exp(-cumulative.intensity)
  illstatus <-1*((latent.illtime<latent.lifetime)&(latent.illtime<censtime))
  latent.waittime[illstatus==1]<-(((-log(U12*S12)*exp(-X12%*%beta12))^(1/shape.waittime))/scale.waittime)[illstatus==1]
  
  
  fit <- data.frame(latent.illtime=latent.illtime,
                    latent.lifetime=latent.lifetime,
                    latent.waittime=latent.waittime,
                    administrative.censoring=administrative.censoring,
                    censtime=censtime)
  fit<-cbind(fit,exogenous_data)
  
  time<-c(1:unique(administrative.censoring))
  
  cumulative.intensity<-(scale.illtime*time)^shape.illtime
  e <- exp(X01%*%beta01)
  cumulative.intensity<-cumulative.intensity%*%t(e)
  survival01 <- exp(-cumulative.intensity)
  
  cumulative.intensity<-(scale.lifetime*time)^shape.lifetime
  e <- exp(X02%*%beta02)
  cumulative.intensity<-cumulative.intensity%*%t(e)
  survival02 <- exp(-cumulative.intensity)
  
  cumulative.intensity<-(scale.waittime*time)^shape.waittime
  e <- exp(X12%*%beta12)
  cumulative.intensity<-cumulative.intensity%*%t(e)
  survival12 <- exp(-cumulative.intensity)
  
  surv01<-data.frame(time=rep(time,n),
             surv=as.vector(survival01),
             id=sort(rep(c(1:(n)),length(time))))
  
  surv01$id<-as.factor(surv01$id)
  p01<-ggplot(surv01[surv01$id%in%c(1:3),],aes(x=time,y=surv,color=id))+
    geom_point()+
    geom_line() +facet_grid(~id)
  
  surv02<-data.frame(time=rep(time,n),
                     surv=as.vector(survival02),
                     id=sort(rep(c(1:(n)),length(time))))
  
  surv02$id<-as.factor(surv02$id)
  p02<-ggplot(surv02[surv02$id%in%c(1:3),],aes(x=time,y=surv,color=id))+
    geom_point()+
    geom_line() +facet_grid(~id)
  
  surv12<-data.frame(time=rep(time,n),
                     surv=as.vector(survival12),
                     id=sort(rep(c(1:(n)),length(time))))
  surv12$id<-as.factor(surv12$id)
  p12<-ggplot(surv12[surv12$id%in%c(1:3),],aes(x=time,y=surv,color=id))+
    geom_point()+
    geom_line() +facet_grid(~id)
  
  
  S01<-exp(-(scale.illtime*time)^shape.illtime)
  S02<-exp(-(scale.lifetime*time)^shape.lifetime)
  S12<-exp(-(scale.waittime*time)^shape.waittime)
  #cens<-exp(-(scale.censtime*time)^shape.censtime)
  
  data.weibull<-matrix(nrow=length(time)*3,ncol=3)
  
  data.weibull[,1]<-c(S01,S02,S12)
  data.weibull[,2]<-c(rep("01",length(S01)),
                       rep("02",length(S02)),
                       rep("12",length(S12)))
  data.weibull[,3]<-rep(time,3)
  colnames(data.weibull)<-c("survie","type","time")
  data.weibull<-as.data.frame(data.weibull)
  data.weibull$survie<-as.numeric(data.weibull$survie)
  data.weibull$time<-as.numeric(data.weibull$time)
  p2<-ggplot2::ggplot(data=data.weibull,aes_string(y="survie",x="time",color="type"))+geom_point()+geom_line()+
    theme_classic()+ylab("Survival")
  
  
  sim.idmModel(x=fit,n=n,plot=list(p2,surv01,p01,surv02,p02,surv12,p12),pen=T)
  
}


