### Code:
##' @title Sample illness-death model data 
##' @description
##'  Simulate data from an illness-death model with interval censored event times
##' and covariates for the purpose of illustrating the help pages of the HIDeM package.
##' See the body of the function for details, i.e., evaluate simulateIDM
##' @param x An \code{idmModel} object as obtained with
##' \code{idmModel}
##' @param n Number of observations
##' seen.ill
##' @param latent if TRUE keep the latent event times
##' @param keep.inspectiontimes if \code{TRUE} keep the inspection
##' times.
##' @param plot plot of base survival for all transition
##' @param ... Extra arguments given to \code{sim}
##' @return A data set with interval censored observations from an illness-death model
##' example(idmModel)
##' help(idmModel)
##' @examples
##'  \dontrun{
#' library(lava)
#' library(prodlim)
#' set.seed(17)
#' d <- simulateIDM(n=1000,
#'                  beta01=c(1,1,0,0.5,0.5,rep(0,5)),
#'                  beta02=c(1,0,0,0,0.5,rep(0,5)),
#'                  beta12=c(1,0,0,0,0.5,rep(0,5)))$data
#' }
##' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
##' @importFrom lava sim
##' @useDynLib HIDeM




sim.idmModel <- function(x,
                         n,
                         plot,
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
  dat$illstatus <- 1*((dat$illtime<dat$latent.lifetime)& (dat$illtime<dat$administrative.censoring))
  
  #dat$illtime[dat$illtime>dat$latent.lifetime] <- 0
  # construct lifetime
  # for ill subjects as the sum of the time to illness (illtime) and
  # the time spent in the illness state (waittime)
  dat$lifetime <- dat$latent.lifetime
  dat$lifetime[dat$illstatus==1]<-dat$latent.waittime[dat$illstatus==1]
  id.nodem.death<-rep(0,n)
  
  #id<-which(T01<T02 & T01<18 & T12<18 & T01>dat$censtime)
  # interval censored illtime
  ipos <- grep("inspection[0-9]+",names(dat))
  
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
        
        if(dat$lifetime[i]>=dat$administrative.censoring[i]){
          if(dat$latent.illtime[i]<dat$censtime[i]){
            idL<-which(itimes<dat$illtime[i])
            idR<-which(itimes>=dat$illtime[i])
            L<-itimes[max(idL)]
            R<-itimes[min(idR)]
            c(L,R,dat$administrative.censoring[i],1,0)
          }else{
            c(dat$censtime[i],dat$censtime[i],dat$administrative.censoring[i],0,0)
          }
        }else{
          if(dat$latent.illtime[i]<dat$censtime[i]){
            idL<-which(itimes<dat$illtime[i])
            idR<-which(itimes>=dat$illtime[i])
            
            #dat$lifetime[i]<R is equivalent to length(idR)==0
            if(length(idR)==0){
              L<-itimes[max(idL)]
              c(L,L,dat$lifetime[i],-1,1)
            }else{
              R<-itimes[min(idR)]
              L<-itimes[max(idL)]
              c(L,R,dat$lifetime[i],1,1)
            }
          }else{
            c(dat$censtime[i],dat$censtime[i],dat$lifetime[i],0,1)
          }
        }
        
        
      }else{
        
        # check administrative censoring for death
        if(dat$lifetime[i]<=dat$administrative.censoring[i]){
          
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
  id.nodem.death[which(dat$seen.ill==-1)]<-1
  dat$seen.ill[dat$seen.ill==-1]<-0
  
  dat$T01<-T01
  dat$T02<-T02
  dat$T12<-T12
  dat$id.nodem.death<-id.nodem.death
  return(list(data=dat,
              plot=plot))
}

simdep.idmModel <- function(x,
                            n,
                            latent=FALSE,
                            keep.inspectiontimes=FALSE,
                            ...){
  # simulate latent data
  #class(x) <- "lvm"
  #dat <- lava::sim(x,n=n,...)

  dat<-x[x$num.visit==1,]
  T01<-dat$latent.illtime
  T02<-dat$latent.lifetime
  T12<-dat$latent.waittime
  dat$illtime <- dat$latent.illtime
  dat$illstatus <- 1*((dat$illtime<dat$latent.lifetime)& (dat$illtime<dat$administrative.censoring))
  #dat$illtime[dat$illtime>dat$latent.lifetime] <- 0
  # construct lifetime
  # for ill subjects as the sum of the time to illness (illtime) and
  # the time spent in the illness state (waittime)
  dat$lifetime <- dat$latent.lifetime
  dat$lifetime[dat$illstatus==1]<-dat$latent.waittime[dat$illstatus==1]
  id.nodem.death<-rep(0,n)

  #id<-which(T01<T02 & T01<18 & T12<18 & T01>dat$censtime)
  # interval censored illtime
  ipos <- grep("inspection[0-9]+",names(dat))

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
      if (length(itimes)==0) {
        itimes <- 0}
      
      ## mark the last inspection time 
      #last.inspection <- itimes[length(itimes)]
      ## find the interval where illness happens
      
      if (dat$illstatus[i]){
        
        if(dat$lifetime[i]>=dat$administrative.censoring[i]){
          if(dat$latent.illtime[i]<dat$censtime[i]){
            idL<-which(itimes<dat$illtime[i])
            idR<-which(itimes>=dat$illtime[i])
            L<-itimes[max(idL)]
            R<-itimes[min(idR)]
            c(L,R,dat$administrative.censoring[i],1,0)
          }else{
            c(dat$censtime[i],dat$censtime[i],dat$administrative.censoring[i],0,0)
          }
        }else{
          if(dat$latent.illtime[i]<dat$censtime[i]){
            idL<-which(itimes<dat$illtime[i])
            idR<-which(itimes>=dat$illtime[i])
            
            #dat$lifetime[i]<R is equivalent to length(idR)==0
            if(length(idR)==0){
              L<-itimes[max(idL)]
              c(L,L,dat$lifetime[i],-1,1)
            }else{
              R<-itimes[min(idR)]
              L<-itimes[max(idL)]
              c(L,R,dat$lifetime[i],1,1)
            }
          }else{
            c(dat$censtime[i],dat$censtime[i],dat$lifetime[i],0,1)
          }
        }
        
        
      }else{
        
        # check administrative censoring for death
        if(dat$lifetime[i]<=dat$administrative.censoring[i]){
          
          c(itimes[length(itimes)],itimes[length(itimes)],dat$lifetime[i],0,
            ifelse(dat$lifetime[i]==0,0,1))
        }else{
          
          itimes<-itimes[itimes<=dat$administrative.censoring[i]]
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
  

  id.nodem.death[which(dat$seen.ill==-1)]<-1
  dat$seen.ill[dat$seen.ill==-1]<-0
  
  dat$T01<-T01
  dat$T02<-T02
  dat$T12<-T12
  dat$id.nodem.death<-id.nodem.death

  
  #delete measure of longitudinal marker after Tdeath or censoring

  col<-colnames(x)[1:(which(colnames(x)=="latent.illtime")-1)]
  dat<-merge(dat[,!colnames(dat)%in%c("num.visit","visit",col[4:length(col)])],x[,c(1:(which(colnames(x)=="latent.illtime")-1))],by=c("ID"))
  
  # Y censoring after diagnosis visit, death or censoring 
  dat<-do.call(rbind,apply(dat,MARGIN = 1,FUN=function(x){
    if(x[names(x)=="seen.ill"]==1){
      
      if(x[names(x)=="R"]<x[names(x)=="visit"] |x[names(x)=="censtime"]<x[names(x)=="visit"]){
        return(NULL)
      }else{
        return(x)
      }
      
    }else{
      if(x[names(x)=="seen.exit"]==1){
      if(x[names(x)=="observed.lifetime"]<x[names(x)=="visit"] |x[names(x)=="censtime"]<x[names(x)=="visit"]){
        return(NULL)
      }else{
        return(x)
      }
    }else{
      if(x[names(x)=="censtime"]<x[names(x)=="visit"]){
        return(NULL)
      }else{
      return(x)}
    }
    }
  }))

  # censoring Y after first event or censoring 
  # dat<-do.call(rbind,apply(dat,MARGIN = 1,FUN=function(x){
  #   if(x[names(x)=="seen.ill"]==1){
  #     
  #     if(x[names(x)=="R"]<x[names(x)=="visit"] ){
  #       return(NULL)
  #     }else{
  #       return(x)
  #     }
  #   }else{
  #   if(x[names(x)=="seen.exit"]==1){
  #     if(x[names(x)=="observed.lifetime"]<x[names(x)=="visit"] |x[names(x)=="censtime"]<x[names(x)=="visit"]){
  #       return(NULL)
  #     }else{
  #       return(x)
  #     }
  #   }else{
  #     if(x[names(x)=="censtime"]<x[names(x)=="visit"]){
  #       return(NULL)
  #     }else{
  #       return(x)}
  #   }
  #   }
  # }))
  dat<-as.data.frame(dat)
  return(dat)
}

##' @title Sample illness-death model data 
##' @description
##'  Simulate data from an illness-death model with interval censored event times
##' and covariates for the purpose of illustrating the help pages of the HIDeM package.
##' See the body of the function for details, i.e., evaluate simulateIDM
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
##' @param cov covariance matrix of explanatory variables
##' @param x01 names of variables on transition 0 --> 1
##' @param x02 names of variables on transition 0 --> 2
##' @param x12 names of variables on transition 1 --> 2
##' @param beta01 value of beta on transition 0 --> 1
##' @param beta02 value of beta on transition 0 --> 2
##' @param beta12 value of beta on transition 1 --> 2
##' @param n number of observations
#' @importFrom ggplot2 ggplot geom_line geom_point theme_classic ylab aes_string aes facet_grid
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
#' @export


simulateIDM <- function(n=100,
                        seed=1,
                        scale.illtime=2.5,
                        shape.illtime=8/100,
                        scale.lifetime=2.5,
                        shape.lifetime=8/100,
                        scale.waittime=2.5,
                        shape.waittime=8/100,
                        prob.censoring=0.05,
                        administrative.censoring=18,
                        n.inspections=8,
                        schedule=2.5,
                        punctuality=0.5,
                        nvar=10,
                        mean=rep(0,10),
                        cov=matrix(c(1,0,0,0,0,0,0,0,0,0,
                                     0,1,0,0,0,0,0,0,0,0,
                                     0,0,1,0,0,0,0,0,0,0,
                                     0,0,0,1,0,0,0,0,0,0,
                                     0,0,0,0,1,0,0,0,0,0,
                                     0,0,0,0,0,1,0,0,0,0,
                                     0,0,0,0,0,0,1,0,0,0,
                                     0,0,0,0,0,0,0,1,0,0,
                                     0,0,0,0,0,0,0,0,1,0,
                                     0,0,0,0,0,0,0,0,0,1),nrow=10,ncol=10),
                        x01=paste0("X",1:10),
                        x02=paste0("X",1:10),
                        x12=paste0("X",1:10),
                        beta01=rep(0.5,10),
                        beta02=rep(0.5,10),
                        beta12=rep(0.5,10)){
  
  ##############################################################################
  ####################### check entry parameters ###############################
  ##############################################################################
  
  if(!inherits(seed,c("numeric","integer")) | length(seed)!=1|seed<=0)stop("The seed has to be a numeric or integer higher than 0.")
  if(!inherits(n,c("numeric","integer")) | round(n)!=n |length(n)!=1|n<=0)stop("The number of subject has to be an integer higher than 0.")
  if(!inherits(prob.censoring,c("numeric","integer")) | prob.censoring > 1 |prob.censoring<0 |length(prob.censoring)!=1)stop("The prob.censoring has to a numeric or integer between 0 and 1.")
  if(!inherits(administrative.censoring,c("numeric","integer")) | administrative.censoring <=0 |length(administrative.censoring)!=1)stop("The administrative.censoring has to be a numeric or integer higher than 0.")
  if(!inherits(n.inspections,c("numeric","integer")) | round(n.inspections)!=n.inspections |length(n.inspections)!=1 | n.inspections <= 0)stop("The n.inspections has to be an integer higher than 0.")
  if(!inherits(schedule,c("numeric","integer")) |length(schedule)!=1 | schedule <= 0)stop("The schedule has to be an integer higher than 0.")
  if(!inherits(punctuality,c("numeric","integer")) | length(punctuality)!=1 | punctuality <= 0)stop("The punctuality has to be an integer or numeric higher than 0.")
  if(!inherits(nvar,c("numeric","integer")) | round(nvar)!=nvar |length(nvar)!=1|nvar<=0)stop("The number of variables has to be an integer higher than 0.")
  if(!inherits(mean,c("numeric","integer")) | length(mean)!=nvar)stop("The mean value of variables has to be an integer or numeric of length nvar.")
  if(!inherits(as.vector(cov),c("numeric","integer")) | dim(cov)[1]!=nvar | dim(cov)[2]!=nvar)stop("The covariance matrix of variables has to contain integer or numeric of dimension nvar*nvar.")
  if(!inherits(x01,c("character")) | any(!x01%in%paste0("X",1:nvar)) | length(x01)>nvar | length(x01)<=0)stop(paste0("The x01 has to contain either X1, X2, ..., until X",nvar))
  if(!inherits(x02,c("character")) | any(!x02%in%paste0("X",1:nvar)) | length(x02)>nvar | length(x02)<=0)stop(paste0("The x02 has to contain either X1, X2, ..., until X",nvar))
  if(!inherits(x12,c("character")) | any(!x12%in%paste0("X",1:nvar)) | length(x12)>nvar | length(x12)<=0)stop(paste0("The x12 has to contain either X1, X2, ..., until X",nvar))
  if(!inherits(beta01,c("numeric","integer")) |  length(beta01)>nvar | length(beta01)<=0)stop("The effets on transition 0 -> 1 beta01 has to be an integer or numeric of length nvar.")
  if(!inherits(beta02,c("numeric","integer")) |  length(beta02)>nvar | length(beta02)<=0)stop("The effets on transition 0 -> 2 beta02 has to be an integer or numeric of length nvar.")
  if(!inherits(beta12,c("numeric","integer")) |  length(beta12)>nvar | length(beta12)<=0)stop("The effets on transition 1 -> 2 beta12 has to be an integer or numeric of length nvar.")
  
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
  # change 23/09/2024 
  #illstatus <-1*((latent.illtime<latent.lifetime)&(latent.illtime<censtime))
  illstatus <-1*((latent.illtime<latent.lifetime)&(latent.illtime<administrative.censoring))
  latent.waittime[illstatus==1]<-(((-log(U12*S12)*exp(-X12%*%beta12))^(1/shape.waittime))/scale.waittime)[illstatus==1]
  latent.waittime[illstatus!=1]<-latent.lifetime[illstatus!=1]
  
  #for(i in 1:n){
  #   if(illstatus[i]==1){
  #     cv<-F
  #     while(cv==F){
  #       U12<-stats::runif(n=1)
  #       latent.waittime[i]<-(((-log(U12)*exp(-X12%*%beta12))^(1/shape.waittime))/scale.waittime)[i]
  #       if(latent.waittime[i]>latent.illtime[i]){cv<-T}
  #     }
  #   }else{
  #     latent.waittime[i]<-latent.lifetime[i]
  #   }
  # }
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
  
  surv<-id<-NULL
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
  
  
  sim.idmModel(x=fit,n=n,plot=list(p2,surv01,p01,surv02,p02,surv12,p12))
  
}

##' @title Sample illness-death model data 
##' @description
##'  Simulate data from an illness-death model with interval censored event times
##' and covariates for the purpose of illustrating the help pages of the HIDeM package.
##' See the body of the function for details, i.e., evaluate simulateIDM
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
##' @param cov covariance matrix of explanatory variables
##' @param x01 names of variables on transition 0 --> 1
##' @param x02 names of variables on transition 0 --> 2
##' @param x12 names of variables on transition 1 --> 2
##' @param beta01 value of beta on transition 0 --> 1
##' @param beta02 value of beta on transition 0 --> 2
##' @param beta12 value of beta on transition 1 --> 2
##' @param n number of observations
##' @param B correlation matrix between random effects of all longitudinal markers
##' @param beta0 Intercept regression parameter for each longitudinal marker
##' @param beta1 Slope regression parameter for each longitudinal marker
##' @param alpha_y_01 Regression parameters of each current value of longitudinal marker on transition 0 to 1
##' @param alpha_y_02 Regression parameters of each current value of longitudinal marker on transition 0 to 2
##' @param alpha_y_12 Regression parameters of each current value of longitudinal marker on transition 1 to 2
##' @param alpha_slope_01 Regression parameters of each slope of longitudinal marker on transition 0 to 1
##' @param alpha_slope_02 Regression parameters of each slope of longitudinal marker on transition 0 to 2
##' @param alpha_slope_12 Regression parameters of each slope of longitudinal marker on transition 1 to 2
#' @importFrom ggplot2 ggplot geom_line geom_point theme_classic ylab aes_string aes facet_grid
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
#' @export


simulateDYNIDM <- function(n=100,
                           seed=1,
                           scale.illtime=2.5,
                           shape.illtime=8/100,
                           scale.lifetime=2.5,
                           shape.lifetime=8/100,
                           scale.waittime=2.5,
                           shape.waittime=8/100,
                           prob.censoring=0.05,
                           administrative.censoring=18,
                           n.inspections=8,
                           schedule=2.5,
                           punctuality=0.5,
                           nvar=10,
                           scale.Y=T,
                           mean=rep(0,10),
                           cov=matrix(c(1,0,0,0,0,0,0,0,0,0,
                                        0,1,0,0,0,0,0,0,0,0,
                                        0,0,1,0,0,0,0,0,0,0,
                                        0,0,0,1,0,0,0,0,0,0,
                                        0,0,0,0,1,0,0,0,0,0,
                                        0,0,0,0,0,1,0,0,0,0,
                                        0,0,0,0,0,0,1,0,0,0,
                                        0,0,0,0,0,0,0,1,0,0,
                                        0,0,0,0,0,0,0,0,1,0,
                                        0,0,0,0,0,0,0,0,0,1),nrow=10,ncol=10),
                           x01=paste0("X",1:10),
                           x02=paste0("X",1:10),
                           x12=paste0("X",1:10),
                           beta01=rep(0.5,10),
                           beta02=rep(0.5,10),
                           beta12=rep(0.5,10),
                           B=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3),
                           sigma=1,
                           beta0=0.5, beta1=0.5,
                           alpha_y_01=0.5, alpha_slope_01=0,
                           alpha_y_02=0.5, alpha_slope_02=0,
                           alpha_y_12=0.5, alpha_slope_12=0){
  
  ##############################################################################
  ####################### check entry parameters ###############################
  ##############################################################################

  if(!inherits(seed,c("numeric","integer")) | length(seed)!=1|seed<=0)stop("The seed has to be a numeric or integer higher than 0.")
  if(!inherits(n,c("numeric","integer")) | round(n)!=n |length(n)!=1|n<=0)stop("The number of subject has to be an integer higher than 0.")
  if(!inherits(prob.censoring,c("numeric","integer")) | prob.censoring > 1 |prob.censoring<0 |length(prob.censoring)!=1)stop("The prob.censoring has to a numeric or integer between 0 and 1.")
  if(!inherits(administrative.censoring,c("numeric","integer")) | administrative.censoring <=0 |length(administrative.censoring)!=1)stop("The administrative.censoring has to be a numeric or integer higher than 0.")
  if(!inherits(n.inspections,c("numeric","integer")) | round(n.inspections)!=n.inspections |length(n.inspections)!=1 | n.inspections <= 0)stop("The n.inspections has to be an integer higher than 0.")
  if(!inherits(schedule,c("numeric","integer")) |length(schedule)!=1 | schedule <= 0)stop("The schedule has to be an integer higher than 0.")
  if(!inherits(punctuality,c("numeric","integer")) | length(punctuality)!=1 | punctuality <= 0)stop("The punctuality has to be an integer or numeric higher than 0.")
  
  if(!inherits(mean,c("numeric","integer")) | length(mean)!=nvar)stop("The mean value of variables has to be an integer or numeric of length nvar.")
  if(!inherits(as.vector(cov),c("numeric","integer")) | dim(cov)[1]!=nvar | dim(cov)[2]!=nvar)stop("The covariance matrix of variables has to contain integer or numeric of dimension nvar*nvar.")
  if(nvar>0){
    if(!inherits(x01,c("character")) | any(!x01%in%paste0("X",1:nvar)) | length(x01)>nvar | length(x01)<=0)stop(paste0("The x01 has to contain either X1, X2, ..., until X",nvar))
    if(!inherits(x02,c("character")) | any(!x02%in%paste0("X",1:nvar)) | length(x02)>nvar | length(x02)<=0)stop(paste0("The x02 has to contain either X1, X2, ..., until X",nvar))
    if(!inherits(x12,c("character")) | any(!x12%in%paste0("X",1:nvar)) | length(x12)>nvar | length(x12)<=0)stop(paste0("The x12 has to contain either X1, X2, ..., until X",nvar))
    if(!inherits(beta01,c("numeric","integer")) |  length(beta01)>nvar | length(beta01)<=0)stop("The effets on transition 0 -> 1 beta01 has to be an integer or numeric of length nvar.")
    if(!inherits(beta02,c("numeric","integer")) |  length(beta02)>nvar | length(beta02)<=0)stop("The effets on transition 0 -> 2 beta02 has to be an integer or numeric of length nvar.")
    if(!inherits(beta12,c("numeric","integer")) |  length(beta12)>nvar | length(beta12)<=0)stop("The effets on transition 1 -> 2 beta12 has to be an integer or numeric of length nvar.")}else{
      beta12<-beta01<-beta02<-0
    }

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
  if(nvar>0){
    exogenous_data<-MASS::mvrnorm(n = n, mu = mean, Sigma = cov)
  }else{exogenous_data<-NULL}
  
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
  tseq <- seq(schedule, by =schedule,length.out=n.inspections-1)
  tunif <- runif(n*length(tseq), -punctuality, punctuality)
  
  for(j in 1:n.inspections){
    if(j==1){
      V[,j]<-0
    }else{
      V[,j]<-tseq[j-1]+tunif[(n*(j-2)+1):(n*(j-1))]
    }
  }
  exogenous_data<-as.data.frame(cbind(exogenous_data,V))
  
  #cov(exogenous_data)<-matrix(0,num_exogenous_vars ,num_exogenous_vars )
  #cov(exogenous_data)[1:nvar,1:nvar] <- cov
  #diag(cov(exogenous_data)[(nvar+1):num_exogenous_vars,(nvar+1):num_exogenous_vars]) <- punctuality
  # Define the structural model for latent processes following a Weibull distribution
  U01<-stats::runif(n=n)
  U02<-stats::runif(n=n)
  U12<-stats::runif(n=n)
  # center and reduce var : 
  #exogenous_data<-scale(exogenous_data)
  if(nvar>0){
    colnames(exogenous_data) <- c(paste0("X", 1:nvar),paste0("inspection",1:n.inspections))
    X01<-as.matrix(exogenous_data[,colnames(exogenous_data)%in%x01])
    X02<-as.matrix(exogenous_data[,colnames(exogenous_data)%in%x02])
    X12<-as.matrix(exogenous_data[,colnames(exogenous_data)%in%x12])
  }else{
    colnames(exogenous_data) <- paste0("inspection",1:n.inspections)
    X01<-X02<-X12<-as.matrix(rep(0,n))
  }
  
  data_long <- c()
  mu_sigma<- 0 #random error are normally distributed around 0 
  B0<-B1<-c()
  

  
  random.effects <- mvtnorm::rmvnorm(n, sigma = B)
  brandom <- data.frame(ID = 1:n, random.effects)
  
  sigma_mat<-matrix(0,ncol=length(beta0),nrow=length(beta0))
  diag(sigma_mat)<-sigma
  error_all<-mvtnorm::rmvnorm(n.inspections*n,sigma =sigma_mat)
  

  for( i in 1:n){
    
    #On tire les effets aléatoires
    b0 <- random.effects[i,seq(1,dim(B)[2],by=2)]
    b1 <- random.effects[i,seq(2,dim(B)[2],by=2)]
    
    B0<-rbind(B0,b0)
    B1<-rbind(B1,b1)
    
    visit<-as.numeric(exogenous_data[i,(dim(exogenous_data)[2]-n.inspections+1):dim(exogenous_data)[2]])
    data_long_i <- c()
    data_long_i <- as.data.frame(cbind(rep(i, n.inspections),
                                       c(1:n.inspections),
                                       visit))
    colnames(data_long_i) <- c("ID", "num.visit",  "visit")
   
    error1<-error_all[((i-1)*n.inspections+1):(i*n.inspections),]
    
    y <- matrix(rep(beta0+b0,n.inspections),ncol=dim(B)[2]/2,nrow=n.inspections,byrow = T)+(visit)%*%t(beta1+b1) + error1
    ty <- matrix(rep(beta0+b0,n.inspections),ncol=dim(B)[2]/2,nrow=n.inspections,byrow = T)+(visit)%*%t(beta1+b1) 
    b0y<-matrix(rep(b0,n.inspections),ncol=dim(B)[2]/2,nrow=n.inspections,byrow = T)
    b1y<-matrix(rep(b1,n.inspections),ncol=dim(B)[2]/2,nrow=n.inspections,byrow = T)
    
    mY<-rep(0,ny)
    sdY<-rep(1,ny)
    if(scale.Y==T){
      
      y0 <- matrix(rep(beta0+b0,n.inspections),ncol=dim(B)[2]/2,nrow=n.inspections,byrow = T) + error1
      mY<-apply(y0,MARGIN=2,FUN = mean)
      sdY<-apply(y0,MARGIN=2,FUN = sd)
      y<-(y-mY)/sdY
    }
    


    colnames(y)<-paste0("Y",c(1:dim(y)[2]))
    colnames(ty)<-paste0("TRUE_Y",c(1:dim(y)[2]))
    colnames(b0y)<-paste0("B0i_Y",c(1:dim(y)[2]))
    colnames(b1y)<-paste0("B1i_Y",c(1:dim(y)[2]))
    data_long<-rbind(data_long,cbind(data_long_i,y,ty,b0y,b1y))
    
  }
  


  data_long$administrative.censoring<-administrative.censoring
  data_long$latent.illtime<-NA
  data_long$latent.lifetime<-NA
  data_long$latent.waittime<-NA
  
  ## generate dropout as Viviane 
  Tdrop <- rep(administrative.censoring,n)
  if(prob.censoring>0)
  {
    idrop <- rbinom(n = n, size = 1, prob = prob.censoring)
    Tdrop[which(idrop == 1)] <- runif(n = length(which(idrop == 1)),
                                      min = schedule,
                                      max = schedule* (n.inspections-1)+punctuality)
  }
 
  newdata <- data.frame(ID = 1:n, Tdrop = Tdrop)
  data_long<-merge(data_long,newdata,by="ID")
  C <- data.frame(
    ID = names(tapply(
      seq_len(nrow(data_long)),
      data_long$ID,
      function(idx) {
        max(data_long$visit[idx][data_long$visit[idx] <= data_long$Tdrop[idx]], na.rm = TRUE)
      }
    )),
    censtime = tapply(
      seq_len(nrow(data_long)),
      data_long$ID,
      function(idx) {
        max(data_long$visit[idx][data_long$visit[idx] <= data_long$Tdrop[idx]], na.rm = TRUE)
      }
    )
  )
  data_long<-merge(data_long,C,by="ID")
  data_long$administrative
  
  computeCL <- function(fixed,random,time,ny,alpha_y,mY,sdY)
  {
    ntime<-length(time)
    di <- matrix(c(rep(1,ntime),time),ncol=2)
    fixed<-matrix(as.vector(fixed),ncol=2,nrow=ny)
    random<-matrix(as.vector(random),ncol=2,nrow=ny)
    cl<-matrix(NA,ncol=ntime,nrow=ny)
    matmY<-matrix(rep(mY,ntime),ncol=ntime,nrow=ny)
    matsdY<-matrix(rep(sdY,ntime),ncol=ntime,nrow=ny)
    
    
    for(k in 1:ny){
      fixedk<-fixed[k,]
      randomk<-random[k,]
      cl[k,] <- di %*% fixedk + di%*%randomk
    }

    cl<-(cl-matmY)/matsdY
    cl<-colSums(alpha_y%*%cl)
  }

  for( i in 1:n){
    
  
    S_01 <- function(tstar,i){
     
      CLy<-computeCL(fixed=c(beta0,beta1),random=c(B0[i,],B1[i,]),ny=ny,time=tstar,alpha_y=alpha_y_01,mY=mY,sdY=sdY)
      
      return(shape.illtime*(scale.illtime^shape.illtime)*((tstar)^(shape.illtime-1))*exp(CLy))
      
      
    }
    
    S_02 <- function(tstar,i){
      
      
      CLy<-computeCL(fixed=c(beta0,beta1),random=c(B0[i,],B1[i,]),ny=ny,time=tstar,alpha_y=alpha_y_02,mY=mY,sdY=sdY)
      return(shape.lifetime*(scale.lifetime^shape.lifetime)*((tstar)^(shape.lifetime-1))*exp(CLy))
      
      
    }
    
    S_12 <- function(tstar,i){
      
      
      CLy<-computeCL(fixed=c(beta0,beta1),random=c(B0[i,],B1[i,]),ny=ny,time=tstar,alpha_y=alpha_y_12,mY=mY,sdY=sdY)
      return(shape.waittime*(scale.waittime^shape.waittime)*((tstar)^(shape.waittime-1))*exp(CLy))
      
      
    }
    
    
    
    
    S_inv <- function(time, f, val, di){
      integrand <- function(t) f(t, di)
      res <- integrate(integrand, lower = 0, upper = time)$value
      res + log(val)
    }
    

    minInv01 <- S_inv(time = 0,f=S_01, val = U01[i], di = i)
    maxInv01 <- S_inv(time = administrative.censoring,f=S_01, val = U01[i], di = i)
    minInv02 <- S_inv(time = 0,f=S_02, val = U02[i], di = i)
    maxInv02 <- S_inv(time = administrative.censoring,f=S_02, val = U02[i], di = i)
    minInv12 <- S_inv(time = 0,f=S_12, val = U12[i], di = i)
    
    if(minInv01>0 | minInv02>0 | minInv12>0 ){
    
      data_long$latent.illtime[data_long$ID==i]<-0
      data_long$latent.lifetime[data_long$ID==i]<-0
      data_long$latent.waittime[data_long$ID==i]<-0
      next
    }
    
    ## if S(Tdrop) > u, ie maxInv < 0, then subject is censored
    if(maxInv01<0 & maxInv02<0){
      
      data_long$latent.illtime[data_long$ID==i]<-100000000
      data_long$latent.lifetime[data_long$ID==i]<-100000000
      data_long$latent.waittime[data_long$ID==i]<-100000000
      
      
    }
    # minIn > 0 then event at inclusion 
    
    T_01 <- try(
      uniroot(
        function(time) S_inv(time, f = S_01, val = U01[i], di = i),
        interval = c(0, administrative.censoring)
      )$root,
      silent = TRUE
    )
    
    T_02 <- try(
      uniroot(
        function(time) S_inv(time, f = S_02, val = U02[i], di = i),
        interval = c(0, administrative.censoring)
      )$root,
      silent = TRUE
    )
    
    if(inherits(T_01, "try-error")){
      T_01 <- 1000000
    }
    if(inherits(T_02, "try-error")){
      T_02 <- 1000000
    }
    
    data_long$latent.illtime[data_long$ID==i]<-T_01
    data_long$latent.lifetime[data_long$ID==i]<-T_02
    data_long$latent.waittime[data_long$ID==i]<-T_02
    illstatus <-1*((T_01<T_02)&(T_01<administrative.censoring))
   
    if(illstatus==1){

      # attention needed integrate otherwise did not correct for survival 
      # 23/04/2026
 
      u12_corrige <- U12[i]*exp(-integrate(S_12,lower=0,upper=T_01,i=i)$value)
      
      minInv12 <- S_inv(time = T_01,f=S_12, val = u12_corrige, di = i)
      maxInv12 <- S_inv(time = administrative.censoring,f=S_12, val = u12_corrige, di = i)

      if(minInv12>0){
      
        data_long$latent.illtime[data_long$ID==i]<-0
        data_long$latent.lifetime[data_long$ID==i]<-0
        data_long$latent.waittime[data_long$ID==i]<-0
        next} 
      
      if(maxInv12<0){ # censored before death 
        data_long$latent.waittime[data_long$ID==i]<-100000000
        next} # had the event prior to the time 
      
      T_12 <- try(
        uniroot(
          function(time) S_inv(time, f = S_12, val = u12_corrige, di = i),
          interval = c(0, administrative.censoring)
        )$root,
        silent = TRUE
      )
      
      if(inherits(T_12, "try-error")){
        T_12 <- 100000000
      }
      
      data_long$latent.waittime[data_long$ID==i]<-T_12
    }else{
      

      # check if infinity in integrand over C --> T_02
      # if it is the case do not consider this case 
      S_surv <- Vectorize(function(tstar, i,ci) {
        lower<-ci
        
        res12 <- integrate(S_12, lower = lower, upper = tstar, i = i)$value
        res01 <- integrate(S_01, lower = lower, upper = tstar, i = i)$value
        res02 <- integrate(S_02, lower = lower, upper = tstar, i = i)$value
        
        res <- exp(-res01 - res02 + res12) * S_01(tstar, i)
        return(res)
      })
      
      ci<-data_long$visit[data_long$ID==i & data_long$visit<=T_02]
      ci<-ifelse(length(ci)==0,0,max(ci))
      ci<-min(ci,administrative.censoring,T_02)
      check<-try(integrate(S_surv,lower=ci,upper=min(T_02,administrative.censoring),i=i,ci=ci)$value
                 ,
                 silent = TRUE
      )
      #check T_12 identifiable à T_02
      
      if(inherits(check, "try-error")){
        check <- 1e+8
      }
      
      if(check>1e+6 | check < 0){
   
        data_long$latent.illtime[data_long$ID==i]<-0
        data_long$latent.lifetime[data_long$ID==i]<-0
        data_long$latent.waittime[data_long$ID==i]<-0
        next} # had the event prior to the time 
      
    }
   
    
  }
  
  exogenous_data$ID<-c(1:n)
  data_long<-merge(x=data_long,exogenous_data,by="ID",all.x=T)

  simdep.idmModel(x=data_long,n=n)
  
}



