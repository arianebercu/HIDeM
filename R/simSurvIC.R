##' Function to generate a latent variable model for interval censored survival times.
##'
##' Based on the functionality of the lava PACKAGE the function generates
##' a latent variable model with a latent time
##' and a censoring mechanism (censtime, inspection1,inspection2,...,inspectionK).
##' 
##' The function \code{\link{sim.survIC}} then simulates
##' interval censored times.
##'
##' @title Generate survival model objects
##' @aliases survModelIC
##' @param scale.time Weilbull scale for latent time
##' @param shape.time Weilbull shape for latent time
##' @param n.inspections Number of inspection times
##' @param schedule Mean of the waiting time between adjacent
##' inspections.
##' @param punctuality Standard deviation of waiting time between
##' inspections.
##' @examples
##' \dontrun{
##' library(lava)
##' library(prodlim)
##' # generate survival model based on exponentially
##' # distributed times
##' m <- survIC(scale.time=1/50, shape.time=0.7)
##' round(sim(m,6),1)
##' 
##' # Estimate the parameters of the Weibull models
##' # based on the uncensored exact event times
##' # and the uncensored illstatus.
##' set.seed(18)
##' d <- sim(m,100,latent=FALSE)
##' d$uncensored.status <- 1
##' f <- shr(Hist(time=list(L,R),event=uncensored.status)~1,
##'          data=d,
##'          conf.int=FALSE)
##' print(f)
##' }
##' @return A latent variable model object \code{lvm}
##' @author Thomas Alexander Gerds
##' @useDynLib SmoothHazardoptim9
##' @export
survIC <- function(scale.time=1/100,
                        shape.time=1,
                        n.inspections=5,
                        schedule=10,
                        punctuality=5){

    ## survival-model
    m <- lava::lvm()
    lava::latent(m) <- c("latent.time")
    lava::distribution(m,"latent.time") <- lava::coxWeibull.lvm(scale=scale.time,shape=shape.time)
    if (n.inspections>0)
        for (k in 1:n.inspections){
            lava::distribution(m,paste("inspection",k,sep="")) <- lava::normal.lvm(mean=schedule,sd=punctuality)
        }
    else{
        stop("Number of inspection times must be an integer greater than 0.")
    }
    class(m) <- c("survIC","lvm")
    m
}
##' Function to simulate interval censored survival data
##'
##' Based on the functionality of the lava PACKAGE 
##' @title Simulate interval censored survival data
##' @param x An \code{survIC} object as obtained with
##' \code{survIC}
##' @param n Number of observations
##' @param compliance Probability of missing an inspection time.
##' @param latent if TRUE keep the latent event times
##' @param keep.inspectiontimes if \code{TRUE} keep the inspection
##' times.
##' @param ... Extra arguments given to \code{sim}
##' @return A data set with interval censored observations 
##' @examples
##' library(lava)
##' example(survIC)
##' help(survIC)
##' ol <- survIC()
##' dat.ol <- sim(ol,10)
##' @author Thomas Alexander Gerds
##' @importFrom lava sim
##' @export
sim.survIC <- function(x,n,compliance=1,latent=TRUE,keep.inspectiontimes=FALSE,...){
    # simulate latent data
    class(x) <- "lvm"
    dat <- lava::sim(x,n=n,...)
    ipos <- grep("inspection[0-9]+",names(dat))
    if (length(ipos)>0) {
        # compute inspection times
        # make sure all inspection times are in the future
        # of the previous inspection time
        iframe <- dat[,ipos]
        dat <- dat[,-ipos,drop=FALSE]
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
                                               ## mark the last inspection time 
                                               last.inspection <- itimes[length(itimes)]
                                               ## find the interval where event happens
                                               if (dat$latent.time[i] > last.inspection){
                                                   ## no event observed until last inspection
                                                   c(last.inspection,Inf,0)
                                               }else{ ## time is smaller or equal to last inspection time
                                                    if (length(itimes)==1){
                                                        c(0,itimes,1)
                                                    } else{
                                                          hit <- prodlim::sindex(eval.times=dat$latent.time[i],
                                                                                 jump.times=itimes,
                                                                                 strict=TRUE)
                                                          c(c(0,itimes)[c(1+hit,2+hit)],1)
                                                      }
                                                }
                                           }))
        colnames(interval) <- c("L","R","status")
        dat <- cbind(dat,interval)
        if (latent==FALSE)
            dat <- dat[,-grep("latent\\.",names(dat))]
        if (keep.inspectiontimes) dat <- cbind(dat,iframe)
        dat
    }
}

#----------------------------------------------------------------------
### simSurvIC.R ends here
