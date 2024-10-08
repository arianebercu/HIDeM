### Code:
# 0 : health state
# 1 : illness state
# 2 : death state
#' Predictions for an illness-death model using either a penalized likelihood
#' approach or a Weibull parametrization.
#' 
#' Predict transition probabilities and cumulative probabilities from an object
#' of class \code{idmSplines} with confidence intervals are calculated.
#' 
#' @param object an \code{idm} class objects returned by a call to the
#'     \code{\link{idm}} function
#' @param s time point at which prediction is made.
#' @param t time horizon for prediction.
#' @param newdata A data frame with covariate values for prediction.
#' @param nsim number of simulations for the confidence intervals
#'     calculations.  The default is 200.
#' @param seed Seed passed to \code{set.seed} for Monte Carlo
#'     simulation of confidence intervals.
#' @param conf.int Level of confidence, i.e., a value between 0 and 1,
#'     the default is \code{0.95}.  The default is also used when
#'     \code{conf.int=TRUE}.  To avoid computation of confidence
#'     intervals, set \code{conf.int} to FALSE or NULL.
#' @param lifeExpect Logical. If \code{TRUE} compute life
#'     expectancies, i.e., \code{t=Inf}.
#' @param maxtime The upper limit of integration for calculations of life expectancies from Weibull parametrizations.
#' @param lambda If penalised illness-death model define if prediction need to be performed 
#' looking at the minimum of BIC or GCV or a vector of three values for each transition 0 -> 1, 0 -> 2 and 1 -> 2.
#' @param ... other parameters.
#' @return a list containing the following predictions with pointwise
#'     confidence intervals: \item{p00}{the transition probability
#'     \eqn{p_{00}}.}  \item{p01}{the transition probability
#'     \eqn{p_{01}}.} \item{p11}{the transition probability
#'     \eqn{p_{11}}.} \item{p12}{the transition probability
#'     \eqn{p_{12}}.} \item{p02_0}{the probability of direct
#'     transition from state 0 to state 2.} \item{p02_1}{the
#'     probability of transition from state 0 to state 2 via state 1.}
#'     \item{p02}{transition probability \eqn{p_{02}}. Note that
#'     \code{p02}=\code{p_02_0}+\code{p02_1}.} \item{F01}{the lifetime
#'     risk of disease. \code{F01}=\code{p01}+\code{p02_1}.}
#'     \item{F0.}{the probability of exit from state
#'     0. \code{F0.}=\code{p02_0}+\code{p01}+\code{p02_1}.}
#' @seealso \code{\link{idm}}
#' @keywords methods
#' @examples
#'
#' \dontrun{
#' set.seed(100)
#' d=simulateIDM(n = 100)$data
#' fit <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2+X3,
#'                formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2+X3,
#'                data=d,conf.int=FALSE)
#' predict(fit,s=0,t=80,conf.int=FALSE,lifeExpect=FALSE)
#' predict(fit,s=0,t=80,nsim=4,conf.int=TRUE,lifeExpect=FALSE)
#' predict(fit,s=0,t=80,nsim=4,conf.int=FALSE,lifeExpect=TRUE)
#' 
#' data(Paq1000)
#' library(prodlim)
#' fit.paq <- idm(formula02=Hist(time=t,event=death,entry=e)~certif,
#' formula01=Hist(time=list(l,r),event=dementia)~certif,data=Paq1000)
#' 
#' predict(fit.paq,s=70,t=80,newdata=data.frame(certif=1))
#' predict(fit.paq,s=70,lifeExpect=TRUE,newdata=data.frame(certif=1))
#' 
#' fit.splines <-  idm(formula02=Hist(time=t,event=death,entry=e)~certif,
#' 		formula01=Hist(time=list(l,r),event=dementia)~certif,
#'                 formula12=~1,
#'                 method="Splines",
#' 		data=Paq1000)
#' 
#' predict(fit.splines,s=70,t=80,newdata=data.frame(certif=1))
#' predict(fit.splines,s=70,t=80,lifeExpect=TRUE,newdata=data.frame(certif=1),nsim=20)
#' 
#' 
#' }
##' @importFrom pracma
#'@useDynLib SmoothHazardoptim9
#' @export
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
predict.idm <- function(object,s,
                        t,newdata,nsim=200,seed=21,conf.int=.95,lifeExpect=FALSE,maxtime,
                        lambda="BIC",...) {
    # check if model has weibull or splines baseline risk 
  if(!is.null(object$modelPar)){
    object$method<-"weib"
    }else{object$method<-"splines"
  }

  if((length(s)>1)|(length(t)>1)) {stop("s and t must be a single value")}
  if((object$penalty!="none") & nsim<=2){stop("With penalised model prediction can only be performed by Monte-Carlo")}
  if(object$penalty!="none"){ nsim<-1}
    if (lifeExpect==TRUE) {
      if (!missing(maxtime) && is.numeric(maxtime)) {
        maxtime <- min(maxtime,object$maxtime)
      } else {
        maxtime <- object$maxtime
      }
      if(nsim>2){ # need value for non montecarole approxi
        t <- Inf}else{
            t<-maxtime
          }
        }
        
    
    if (any(s>t)) {stop("You must respect the condition 's<t' to calculate p(s,t)")}
    do.conf.int <- !is.null(conf.int) && !is.na(conf.int) && !conf.int==FALSE
    if (is.logical(conf.int)) conf.int <- .95
    if (do.conf.int == TRUE){
        stopifnot(0<conf.int && conf.int<1)
        #if (nsim < 2) stop("Need at least two simulations to construct confidence limits.")
    }
    if (missing(t) && lifeExpect==FALSE) stop("Argument t is missing.")
    if ((lifeExpect==TRUE) & (nsim>2)) t <- Inf
    if (missing(s)) stop("Argument s is missing.")
    ## if (missing(t) || is.infinite(t)) lifeExpect <- TRUE
    # if covariates: cov=c(cov1,cov2,cov3,...)
    nvar01 <- object$NC[1]
    nvar02 <- object$NC[2]
    nvar12 <- object$NC[3]

    # keep previous levels 
    object$levels$class[sapply(object$levels$class, is.null)] <- NA
    if(length(object$levels$values)>0){
      object$levels$values[sapply(object$levels$values, is.null)] <- NA}
    xlevels<-unlist(object$levels$class)
    if(any(xlevels%in%c("factor","character"))){
      xnames<-c(all.vars(object$terms$Formula01)[all.vars(object$terms$Formula01)%in%labels(terms(object$terms$Formula01))],
                      all.vars(object$terms$Formula02)[all.vars(object$terms$Formula02)%in%labels(terms(object$terms$Formula02))],
                      all.vars(object$terms$Formula12)[all.vars(object$terms$Formula12)%in%labels(terms(object$terms$Formula12))])
      xnamesfactor<-xnames[which(xlevels%in%c("factor","character"))]
      if(length(xnamesfactor)>0){
        id<-which(xlevels%in%c("factor","character"))
        m<-1
        for(k in xnamesfactor){
          newdata[,k] <- factor(newdata[,k], levels=object$levels$values[[id[m]]])
          m<-m+1 }
      }
      
      
    }
    
    #################### prediction if model not from penalty ##################
    if(object$penalty=="none"){
    #update dataset from the formula 
    if (!missing(newdata)){
        if (NROW(newdata)>1) stop("Argument newdata has more than one row\n.Currently this function works only for one covariate constallation at a time.")
        if (length(object$Xnames01)>0){
          if(length(grep(":",names(object$coef[1:object$NC[1]])))>0){
          #Z01 <- model.matrix(object$terms$Formula01,data=newdata)[, -1, drop = FALSE]
          Z01 <- model.matrix(update.formula(formula(object$terms$Formula01),NULL~.),data=newdata)[, -1, drop = FALSE]
          }else{
            Z01 <-as.matrix(model.frame(formula=update.formula(formula(object$terms$Formula01),NULL~.),data=newdata))}
          }else{
            Z01 <- 0}
        if (length(object$Xnames02)>0){
          
          if(length(grep(":",names(object$coef[(1+object$NC[1]):(object$NC[2]+object$NC[1])])))>0){
            #Z02 <- model.matrix(object$terms$Formula02,data=newdata)[, -1, drop = FALSE]
            Z02 <- model.matrix(update.formula(formula(object$terms$Formula02),NULL~.),data=newdata)[, -1, drop = FALSE]
          }else{
            Z02 <-as.matrix(model.frame(formula=update.formula(formula(object$terms$Formula02),NULL~.),data=newdata))}
        }else{
            Z02 <- 0}
        if (length(object$Xnames12)>0){
          if(length(grep(":",names(object$coef[(1+object$NC[1]+object$NC[2]):(object$NC[3]+object$NC[2]+object$NC[1])])))>0){
            #Z12 <- model.matrix(object$terms$Formula12,data=newdata)[, -1, drop = FALSE]
            Z12 <- model.matrix(update.formula(formula(object$terms$Formula12),NULL~.),data=newdata)[, -1, drop = FALSE]
          }else{
            Z12 <-as.matrix(model.frame(formula=update.formula(formula(object$terms$Formula12),NULL~.),data=newdata))}
          }else{
            Z12 <- 0}
    }else{
        vars <- unique(c(object$Xnames01,object$Xnames02,object$Xnames12))
        newdata <- data.frame(matrix(0,ncol=pmax(1,length(vars))))
        names(newdata) <- vars
        Z01 <- matrix(rep(0,length(object$Xnames01)),nrow=1)
        Z02 <- matrix(rep(0,length(object$Xnames02)),nrow=1)
        Z12 <- matrix(rep(0,length(object$Xnames12)),nrow=1)
    }
    
    if(nvar01 > 0){
        beta01 <- object$coef[1:nvar01,]
        names(beta01) <- paste0("beta01.",names(beta01))
        bZ01 <- sum(Z01 * beta01)
    }else{
        bZ01 <- 0
        beta01 <- NULL
    }
    if (nvar02 != 0) {
        beta02 <- object$coef[(nvar01+1):(nvar01+nvar02),]
        names(beta02) <- paste0("beta02.",names(beta02))
        bZ02 <- sum(Z02 * beta02)
    }else{
        beta02 <- NULL
        bZ02 <- 0
    }
    if (nvar12 != 0) {
        beta12 <- object$coef[(nvar01+nvar02+1):(nvar01+nvar02+nvar12),]
        names(beta12) <- paste0("beta12.",names(beta12))
        bZ12 <- sum(Z12 * beta12)
    }else{
        beta12 <- NULL
        bZ12 <- 0
    }

    uppercumulative.intensity<-lowercumulative.intensity<-cumulative.intensity<-NULL
    upperintensity<-lowerintensity<-intensity<-NULL
    ############### splines ####################################################
    if (object$method=="splines"){

        nknots01 <- object$nknots01
        nknots02 <- object$nknots02
        nknots12 <- object$nknots12
        knots01 <- object$knots01
        knots02 <- object$knots02
        knots12 <- object$knots12
        the01 <- object$theta01
        
        if(s<min(knots01,knots02,knots12))stop(paste0("The argument s must be at least equal to :",min(knots01,knots02,knots12)))
        if(t<min(knots01,knots02,knots12))stop(paste0("The argument t must be at least equal to :",min(knots01,knots02,knots12)))
        
        names(the01) <- paste0("the01.",1:length(the01))
        the02 <- object$theta02
        names(the02) <- paste0("the02.",1:length(the02))
        the12 <- object$theta12
        names(the12) <- paste0("the12.",1:length(the12))
        if (do.conf.int == TRUE & nsim>2){
            ### conf.int prediction by Monte-Carlo
            Vmean <- c(the01,the02,the12,beta01,beta02,beta12) # vector of estimates
            #Vmean<-Vmean[object$fix==0]   # vector of estimates not fixed
            set.seed(seed)
            X <- mvtnorm::rmvnorm(nsim,Vmean,object$V)
            colnames(X) <- names(Vmean)
            # 1 set of simulated parameters for each element of the list
            Xtheta01=X[,names(the01)]
            Xtheta02=X[,names(the02)]
            Xtheta12=X[,names(the12)]
            if (!is.null(beta01))
                linPred01=X[,names(beta01),drop=FALSE] %*% t(Z01)
            else
                linPred01=matrix(0,nrow=nsim,ncol=1)
            if (!is.null(beta02))
                linPred02=X[,names(beta02),drop=FALSE] %*% t(Z02)
            else
                linPred02=matrix(0,nrow=nsim,ncol=1)
            if (!is.null(beta12))
                linPred12=X[,names(beta12),drop=FALSE] %*% t(Z12)
            else
                linPred12=matrix(0,nrow=nsim,ncol=1)
           
            if (lifeExpect==TRUE){
                simResults <- do.call("rbind",
                                      lapply(1:nsim,function(i){
                                          lifexpect0.idmPl(s,
                                                           knots01,
                                                           nknots01,
                                                           Xtheta01[i,],
                                                           knots12,
                                                           nknots12,
                                                           Xtheta12[i,],
                                                           knots02,
                                                           nknots02,
                                                           Xtheta02[i,],
                                                           linPred01[i,],
                                                           linPred12[i,],
                                                           linPred02[i,])
                                      }))
            }else{
                simResults <- do.call("rbind",
                                      lapply(1:nsim,function(i){
                                          Predict0.idmPl(s,
                                                         t,
                                                         knots01,
                                                         nknots01,
                                                         Xtheta01[i,],
                                                         knots12,
                                                         nknots12,
                                                         Xtheta12[i,],
                                                         knots02,
                                                         nknots02,
                                                         Xtheta02[i,],
                                                         linPred01[i,],
                                                         linPred12[i,],
                                                         linPred02[i,])
                                      }))
            }
            q.lower <- (1-conf.int)/2
            q.upper <- 1-q.lower
            ci <- apply(simResults,2,function(x)quantile(unlist(x),c(q.lower,q.upper)))
        }
        #browser()
        # want to calculate variability based on variance-covariance matrix
        if (do.conf.int == TRUE & nsim<=2){

          knots.unique<-unique(object$knots01)
          knots.bound<-knots.unique[c(1,length(knots.unique))]
          knots.int<-knots.unique[-c(1,length(knots.unique))]
          msplines01<-splinesMI(x=t,knots=sort(knots.int),Boundary.knots=knots.bound)$Mspline
          isplines01<-(splinesMI(x=t,knots=sort(knots.int),Boundary.knots=knots.bound)$Ispline-splinesMI(x=s,knots=sort(knots.int),Boundary.knots=knots.bound)$Ispline)
          
          knots.unique<-unique(object$knots02)
          knots.bound<-knots.unique[c(1,length(knots.unique))]
          knots.int<-knots.unique[-c(1,length(knots.unique))]
          msplines02<-splinesMI(x=t,knots=sort(knots.int),Boundary.knots=knots.bound)$Mspline
          isplines02<-(splinesMI(x=t,knots=sort(knots.int),Boundary.knots=knots.bound)$Ispline-splinesMI(x=s,knots=sort(knots.int),Boundary.knots=knots.bound)$Ispline)
          
          knots.unique<-unique(object$knots12)
          knots.bound<-knots.unique[c(1,length(knots.unique))]
          knots.int<-knots.unique[-c(1,length(knots.unique))]
          msplines12<-splinesMI(x=t,knots=sort(knots.int),Boundary.knots=knots.bound)$Mspline
          isplines12<-(splinesMI(x=t,knots=sort(knots.int),Boundary.knots=knots.bound)$Ispline-splinesMI(x=s,knots=sort(knots.int),Boundary.knots=knots.bound)$Ispline)
          
          
          theta.square01<-the01^2
          intensity01<-msplines01%*%theta.square01
          cumulative.intensity01<-isplines01%*%theta.square01
          
          theta.square02<-the02^2
          intensity02<-msplines02%*%theta.square02
          cumulative.intensity02<-isplines02%*%theta.square02
          
          theta.square12<-the12^2
          intensity12<-msplines12%*%theta.square12
          cumulative.intensity12<-isplines12%*%theta.square12
          
          if (!is.null(beta01))
            linPred01<-beta01 %*% t(Z01)
          else
            linPred01<-0
          if (!is.null(beta02))
            linPred02<-beta02 %*% t(Z02)
          else
            linPred02<-0
          if (!is.null(beta12))
            linPred12<-beta12 %*% t(Z12)
          else
            linPred12<-0
          
          
          e01 <- exp(linPred01)
          intensity01<-intensity01*e01
          cumulative.intensity01<-cumulative.intensity01*e01
          
          e02 <- exp(linPred02)
          intensity02<-intensity02*e02
          cumulative.intensity02<-cumulative.intensity02*e02
          
          e12 <- exp(linPred12)
          intensity12<-intensity12*e12
          cumulative.intensity12<-cumulative.intensity12*e12
          
          intensity<-c(intensity01,intensity02,intensity12)
          cumulative.intensity<-c(cumulative.intensity01,cumulative.intensity02,cumulative.intensity12)
          
          lowerintensity<-rep(NA,3)
          upperintensity<-rep(NA,3)
          lowercumulative.intensity<-rep(NA,3)
          uppercumulative.intensity<-rep(NA,3)
          
           #V<-object$V[object$fix==0,object$fix==0]
           Nspline<-object$nknots01+object$nknots02+object$nknots12+6+1
           V01<-object$V[c(1:(object$nknots01+2),Nspline:(Nspline+object$NC[1]-1)),c(1:(object$nknots01+2),Nspline:(Nspline+object$NC[1]-1))]
           Nspline<-object$nknots01+object$nknots02+object$nknots12+6+object$NC[1]+1
           V02<-object$V[c((object$nknots01+3):(object$nknots01+object$nknots02+4),Nspline:(Nspline+object$NC[2]-1)),c((object$nknots01+3):(object$nknots01+object$nknots02+4),Nspline:(Nspline+object$NC[2]-1))]
           Nspline<-object$nknots01+object$nknots02+object$nknots12+6+object$NC[1]+object$NC[2]+1
           V12<-object$V[c((object$nknots01+object$nknots02+5):(object$nknots01+object$nknots02+object$nknots12+6),Nspline:(Nspline+object$NC[3]-1)),c((object$nknots01+object$nknots02+5):(object$nknots01+object$nknots02+object$nknots12+6),Nspline:(Nspline+object$NC[3]-1))]
           
           
            deriv01<-c(2*the01*msplines01*as.numeric(e01),as.numeric(intensity01)*Z01)
            se01<-sqrt(deriv01%*%V01%*%deriv01)
            
            deriv02<-c(2*the02*msplines02*as.numeric(e02),as.numeric(intensity02)*Z02)
            se02<-sqrt(deriv02%*%V02%*%deriv02)
            
            deriv12<-c(2*the12*msplines12*as.numeric(e12),as.numeric(intensity12)*Z12)
            se12<-sqrt(deriv12%*%V12%*%deriv12)
            
            lowerintensity[1]<-intensity01+qnorm((1-conf.int)/2)*se01
            lowerintensity[2]<-intensity02+qnorm((1-conf.int)/2)*se02
            lowerintensity[3]<-intensity12+qnorm((1-conf.int)/2)*se12
            upperintensity[1]<-intensity01-qnorm((1-conf.int)/2)*se01
            upperintensity[2]<-intensity02-qnorm((1-conf.int)/2)*se02
            upperintensity[3]<-intensity12-qnorm((1-conf.int)/2)*se12
            
            deriv01<-c(2*the01*isplines01*as.numeric(e01),as.numeric(cumulative.intensity01)*Z01)
            se01<-sqrt(deriv01%*%V01%*%deriv01)
            
            deriv02<-c(2*the02*isplines02*as.numeric(e02),as.numeric(cumulative.intensity02)*Z02)
            se02<-sqrt(deriv02%*%V02%*%deriv02)
            
            deriv12<-c(2*the12*isplines12*as.numeric(e12),as.numeric(cumulative.intensity12)*Z12)
            se12<-sqrt(deriv12%*%V12%*%deriv12)
            
            lowercumulative.intensity[1]<-cumulative.intensity01+qnorm((1-conf.int)/2)*se01
            lowercumulative.intensity[2]<-cumulative.intensity02+qnorm((1-conf.int)/2)*se02
            lowercumulative.intensity[3]<-cumulative.intensity12+qnorm((1-conf.int)/2)*se12
            uppercumulative.intensity[1]<-cumulative.intensity01-qnorm((1-conf.int)/2)*se01
            uppercumulative.intensity[2]<-cumulative.intensity02-qnorm((1-conf.int)/2)*se02
            uppercumulative.intensity[3]<-cumulative.intensity12-qnorm((1-conf.int)/2)*se12
            
              }
         
        
        if (lifeExpect==TRUE){
            transprob <- unlist(lifexpect0.idmPl(s,
                                                 knots01,
                                                 nknots01,
                                                 the01,
                                                 knots12,
                                                 nknots12,
                                                 the12,
                                                 knots02,
                                                 nknots02,
                                                 the02,
                                                 bZ01,
                                                 bZ12,
                                                 bZ02))
        }else{
            transprob <- unlist(Predict0.idmPl(s,t,knots01,nknots01,the01,knots12,nknots12,the12,knots02,nknots02,the02,bZ01,bZ12,bZ02))
        }
    }else {
      ############################ weibull #####################################
        a01 <- object$modelPar[1]
        b01 <- object$modelPar[2]
        a02 <- object$modelPar[3]
        b02 <- object$modelPar[4]
        a12 <- object$modelPar[5]
        b12 <- object$modelPar[6]
        if (do.conf.int==TRUE & nsim>2) {
            ## conf.int prediction by Monte-Carlo
            ## vector of parameter estimates
            Vmean <- c(a01,b01,a02,b02,a12,b12,beta01,beta02,beta12)
            names(Vmean) <- c("a01","b01","a02","b02","a12","b12",names(beta01),names(beta02),names(beta12))
            set.seed(seed)
            X <- mvtnorm::rmvnorm(nsim,Vmean,object$V)
            colnames(X) <- names(Vmean)
            # set of simulated parameters for each element of the list
            Xa01=X[,"a01"]^2
            Xb01=1/(X[,"b01"]^2)
            Xa02=X[,"a02"]^2
            Xb02=1/(X[,"b02"]^2)
            Xa12=X[,"a12"]^2
            Xb12=1/(X[,"b12"]^2)
            if (!is.null(beta01))
                linPred01=X[,names(beta01),drop=FALSE] %*% t(Z01)
            else
                linPred01=matrix(0,nrow=nsim,ncol=1)
            if (!is.null(beta02))
                linPred02=X[,names(beta02),drop=FALSE] %*% t(Z02)
            else
                linPred02=matrix(0,nrow=nsim,ncol=1)
            if (!is.null(beta12))
                linPred12=X[,names(beta12),drop=FALSE] %*% t(Z12)
            else
                linPred12=matrix(0,nrow=nsim,ncol=1)
            
            if (lifeExpect==TRUE){
                simResults <- do.call("rbind",
                                      lapply(1:nsim,function(i){
                                          lifexpect0.idmWeib(s,
                                                             a01=Xa01[[i]],
                                                             b01=Xb01[[i]],
                                                             a02=Xa02[[i]],
                                                             b02=Xb02[[i]],
                                                             a12=Xa12[[i]],
                                                             b12=Xb12[[i]],
                                                             bZ01=linPred01[[i]],
                                                             bZ02=linPred02[[i]],
                                                             bZ12=linPred12[[i]],max=maxtime)
                                      }))
            }else{
                simResults <- do.call("rbind",
                                      lapply(1:nsim,function(i){
                                          Predict0.idmWeib(s,
                                                           t,
                                                           a01=Xa01[[i]],
                                                           b01=Xb01[[i]],
                                                           a02=Xa02[[i]],
                                                           b02=Xb02[[i]],
                                                           a12=Xa12[[i]],
                                                           b12=Xb12[[i]],
                                                           bZ01=linPred01[[i]],
                                                           bZ02=linPred02[[i]],
                                                           bZ12=linPred12[[i]])
                                      }))
            }
            q.lower <- (1-conf.int)/2
            q.upper <- 1-q.lower
            ci <- apply(simResults,2,function(x)quantile(unlist(x),c(q.lower,q.upper)))
        }
        # want to calculate variability based on variance-covariance matrix
        if (do.conf.int==TRUE & nsim<=2) {
          # need to test
          
          
          modelPar01<-object$modelPar[1:2]^2
          intensity01<-modelPar01[1]*(modelPar01[2]^modelPar01[1])*t^(modelPar01[1]-1)
          cumulative.intensity01<-(modelPar01[2]*t)^modelPar01[1]
          
          modelPar02<-object$modelPar[3:4]^2
          intensity02<-modelPar02[1]*(modelPar02[2]^modelPar02[1])*t^(modelPar02[1]-1)
          cumulative.intensity02<-(modelPar02[2]*t)^modelPar02[1]
          
          modelPar12<-object$modelPar[5:6]^2
          intensity12<-modelPar12[1]*(modelPar12[2]^modelPar12[1])*t^(modelPar12[1]-1)
          cumulative.intensity12<-(modelPar12[2]*t)^modelPar12[1]
          
          if (!is.null(beta01))
            linPred01<-beta01 %*% t(Z01)
          else
            linPred01<-0
          if (!is.null(beta02))
            linPred02<-beta02 %*% t(Z02)
          else
            linPred02<-0
          if (!is.null(beta12))
            linPred12<-beta12 %*% t(Z12)
          else
            linPred12<-0
          
          
          e01 <- exp(linPred01)
          intensity01<-intensity01*e01
          cumulative.intensity01<-cumulative.intensity01*e01
          
          e02 <- exp(linPred02)
          intensity02<-intensity02*e02
          cumulative.intensity02<-cumulative.intensity02*e02
          
          e12 <- exp(linPred12)
          intensity12<-intensity12*e12
          cumulative.intensity12<-cumulative.intensity12*e12
          
          intensity<-c(intensity01,intensity02,intensity12)
          cumulative.intensity<-c(cumulative.intensity01,cumulative.intensity02,cumulative.intensity12)
          lowerintensity<-rep(NA,3)
          upperintensity<-rep(NA,3)
          lowercumulative.intensity<-rep(NA,3)
          uppercumulative.intensity<-rep(NA,3)
          
          
          V01<-object$V[c(1:2,7:(7+object$NC[1]-1)),c(1:2,7:(7+object$NC[1]-1))]
          V02<-object$V[c(3:4,(6+object$NC[1]+1):(6+object$NC[1]+object$NC[2])),c(3:4,(6+object$NC[1]+1):(6+object$NC[1]+object$NC[2]))]
          V12<-object$V[c(5:6,(6+object$NC[1]+object$NC[2]+1):(6+object$NC[1]+object$NC[2]+object$NC[3])),c(5:6,(6+object$NC[1]+object$NC[2]+1):(6+object$NC[1]+object$NC[2]+object$NC[3]))]

          
          deriv01<-c((t^(modelPar01[1]-1))*(modelPar01[2]^modelPar01[1])*(modelPar01[1]*log(t)+modelPar01[1]*log(modelPar01[2])+1),
                     (modelPar01[1]^2)*(t^(modelPar01[1]-1))*(modelPar01[2]^(modelPar01[1]-1)),
                     as.numeric(intensity01)*Z01)
          deriv01[c(1,2)]<-deriv01[c(1,2)]*(2*sqrt(modelPar01))
          se01<-sqrt(deriv01%*%V01%*%deriv01)
          
          deriv02<-c((t^(modelPar02[1]-1))*(modelPar02[2]^modelPar02[1])*(modelPar02[1]*log(t)+modelPar02[1]*log(modelPar02[2])+1),
                     (modelPar02[1]^2)*(t^(modelPar02[1]-1))*(modelPar02[2]^(modelPar02[1]-1)),
                     as.numeric(intensity02)*Z02)
          deriv02[c(1,2)]<-deriv02[c(1,2)]*(2*sqrt(modelPar02))
          se02<-sqrt(deriv02%*%V02%*%deriv02)
          
          deriv12<-c((t^(modelPar12[1]-1))*(modelPar12[2]^modelPar12[1])*(modelPar12[1]*log(t)+modelPar12[1]*log(modelPar12[2])+1),
                     (modelPar12[1]^2)*(t^(modelPar12[1]-1))*(modelPar12[2]^(modelPar12[1]-1)),
                     as.numeric(intensity12)*Z12)
          deriv12[c(1,2)]<-deriv12[c(1,2)]*(2*sqrt(modelPar12))
          se12<-sqrt(deriv12%*%V12%*%deriv12)
          
          lowerintensity[1]<-intensity01+qnorm((1-conf.int)/2)*se01
          lowerintensity[2]<-intensity02+qnorm((1-conf.int)/2)*se02
          lowerintensity[3]<-intensity12+qnorm((1-conf.int)/2)*se12
          upperintensity[1]<-intensity01-qnorm((1-conf.int)/2)*se01
          upperintensity[2]<-intensity02-qnorm((1-conf.int)/2)*se02
          upperintensity[3]<-intensity12-qnorm((1-conf.int)/2)*se12
          
          deriv01<-c(((t*modelPar01[2])^modelPar01[1])*log(t*modelPar01[2]),
                     (modelPar01[1]*((t*modelPar01[2])^modelPar01[1]))/modelPar01[2],
                     as.numeric(cumulative.intensity01)*Z01)
          deriv01[c(1,2)]<-deriv01[c(1,2)]*(2*sqrt(modelPar01))
          se01<-sqrt(deriv01%*%V01%*%deriv01)
          
          deriv02<-c(((t*modelPar02[2])^modelPar02[1])*log(t*modelPar02[2]),
                     (modelPar02[1]*((t*modelPar02[2])^modelPar02[1]))/modelPar02[2],
                     as.numeric(cumulative.intensity02)*Z02)
          deriv02[c(1,2)]<-deriv02[c(1,2)]*(2*sqrt(modelPar02))
          se02<-sqrt(deriv02%*%V02%*%deriv02)
          
          deriv12<-c(((t*modelPar12[2])^modelPar12[1])*log(t*modelPar12[2]),
                     (modelPar12[1]*((t*modelPar12[2])^modelPar12[1]))/modelPar12[2],
                     as.numeric(cumulative.intensity12)*Z12)
          deriv12[c(1,2)]<-deriv12[c(1,2)]*(2*sqrt(modelPar12))
          se12<-sqrt(deriv12%*%V12%*%deriv12)
          
          lowercumulative.intensity[1]<-cumulative.intensity01+qnorm((1-conf.int)/2)*se01
          lowercumulative.intensity[2]<-cumulative.intensity02+qnorm((1-conf.int)/2)*se02
          lowercumulative.intensity[3]<-cumulative.intensity12+qnorm((1-conf.int)/2)*se12
          uppercumulative.intensity[1]<-cumulative.intensity01-qnorm((1-conf.int)/2)*se01
          uppercumulative.intensity[2]<-cumulative.intensity02-qnorm((1-conf.int)/2)*se02
          uppercumulative.intensity[3]<-cumulative.intensity12-qnorm((1-conf.int)/2)*se12
          
          
          
        }
        if (lifeExpect==TRUE){
            transprob <- unlist(lifexpect0.idmWeib(s,
                                                   a01^2,
                                                   1/b01^2,
                                                   a02^2,
                                                   1/b02^2,
                                                   a12^2,
                                                   1/b12^2,
                                                   bZ01,
                                                   bZ02,
                                                   bZ12,max=maxtime))
        }else{
            transprob <- unlist(Predict0.idmWeib(s,
                                                 t,
                                                 a01^2,
                                                 1/b01^2,
                                                 a02^2,
                                                 1/b02^2,
                                                 a12^2,
                                                 1/b12^2,
                                                 bZ01,
                                                 bZ02,
                                                 bZ12))
        }
    }

    if (do.conf.int==TRUE & nsim>2){
        transprob <- data.frame(cbind(transprob,t(ci)))
        names(transprob) <- c("Estimate",paste("Lower",round(100*conf.int),sep="."),paste("Upper",round(100*conf.int),sep="."))
        transprob <- cbind("Parameter"=rownames(transprob),transprob)
        rownames(transprob) <- NULL
    }else{
        transprob <- data.frame(cbind(transprob))
        names(transprob) <- c("Estimate")
        transprob <- cbind("Parameter"=rownames(transprob),transprob)
        rownames(transprob) <- NULL
    }
  
    }else{
      #################### predict with penalty ################################
    # lambda need to be either a vector of three values (01,02,12) or BIC or GCV

      if(is.null(lambda)){lambda<-"BIC"}
      if(length(lambda)==1 & (!lambda%in%c("GCV","BIC"))){stop("Lambda need to be either a vector of three values (01,02 and 12) or BIC or GCV")}
      if(!length(lambda)%in%c(1,3)){stop("Lambda need to be either a vector of three values (01,02 and 12) or BIC or GCV")}

      if(length(lambda)==3){
        if(any(sapply(object$lambda,FUN=function(x){sum(x==lambda)})==3)){stop("Lambda need to be either a vector of three values (01,02 and 12) from object$lambda")}
        id<-which(any(sapply(object$lambda,FUN=function(x){sum(x==lambda)})==3))[1]
      }
      if(length(lambda)==1){
        if(lambda=="BIC"){id<-which.min(object$BIC)[1]
        }else{id<-which.min(object$GCV)[1]}
      }
      #update dataset from the formula 
      if (!missing(newdata)){
        if (NROW(newdata)>1) stop("Argument newdata has more than one row\n.Currently this function works only for one covariate constallation at a time.")
        if (length(object$Xnames01)>0){
          if(length(grep(":",names(object$coef[1:object$NC[1]])))>0){
            #Z01 <- model.matrix(object$terms$Formula01,data=newdata)[, -1, drop = FALSE]
            Z01 <- model.matrix(update.formula(formula(object$terms$Formula01),NULL~.),data=newdata)[, -1, drop = FALSE]
          }else{
            Z01 <-as.matrix(model.frame(formula=update.formula(formula(object$terms$Formula01),NULL~.),data=newdata))}
        }else{
          Z01 <- 0}
        if (length(object$Xnames02)>0){
          
          if(length(grep(":",names(object$coef[(1+object$NC[1]):(object$NC[2]+object$NC[1])])))>0){
            #Z02 <- model.matrix(object$terms$Formula02,data=newdata)[, -1, drop = FALSE]
            Z02 <- model.matrix(update.formula(formula(object$terms$Formula02),NULL~.),data=newdata)[, -1, drop = FALSE]
          }else{
            Z02 <-as.matrix(model.frame(formula=update.formula(formula(object$terms$Formula02),NULL~.),data=newdata))}
        }else{
          Z02 <- 0}
        if (length(object$Xnames12)>0){
          if(length(grep(":",names(object$coef[(1+object$NC[1]+object$NC[2]):(object$NC[3]+object$NC[2]+object$NC[1])])))>0){
            #Z12 <- model.matrix(object$terms$Formula12,data=newdata)[, -1, drop = FALSE]
            Z12 <- model.matrix(update.formula(formula(object$terms$Formula12),NULL~.),data=newdata)[, -1, drop = FALSE]
          }else{
            Z12 <-as.matrix(model.frame(formula=update.formula(formula(object$terms$Formula12),NULL~.),data=newdata))}
        }else{
          Z12 <- 0}
      }else{
        vars <- unique(c(object$Xnames01,object$Xnames02,object$Xnames12))
        newdata <- data.frame(matrix(0,ncol=pmax(1,length(vars))))
        names(newdata) <- vars
        Z01 <- matrix(rep(0,length(object$Xnames01)),nrow=1)
        Z02 <- matrix(rep(0,length(object$Xnames02)),nrow=1)
        Z12 <- matrix(rep(0,length(object$Xnames12)),nrow=1)
      }
      
      if(nvar01 > 0){
        beta01 <- object$coef[1:nvar01,id]
        names(beta01) <- paste0("beta01.",names(beta01))
        bZ01 <- sum(Z01 * beta01)
      }else{
        bZ01 <- 0
        beta01 <- NULL
      }
      if (nvar02 != 0) {
        beta02 <- object$coef[(nvar01+1):(nvar01+nvar02),id]
        names(beta02) <- paste0("beta02.",names(beta02))
        bZ02 <- sum(Z02 * beta02)
      }else{
        beta02 <- NULL
        bZ02 <- 0
      }
      if (nvar12 != 0) {
        beta12 <- object$coef[(nvar01+nvar02+1):(nvar01+nvar02+nvar12),id]
        names(beta12) <- paste0("beta12.",names(beta12))
        bZ12 <- sum(Z12 * beta12)
      }else{
        beta12 <- NULL
        bZ12 <- 0
      }

      uppercumulative.intensity<-lowercumulative.intensity<-cumulative.intensity<-NULL
      upperintensity<-lowerintensity<-intensity<-NULL
      ############### splines ####################################################
      if (object$method=="splines"){
        
        nknots01 <- object$nknots01
        nknots02 <- object$nknots02
        nknots12 <- object$nknots12
        knots01 <- object$knots01
        knots02 <- object$knots02
        knots12 <- object$knots12
        the01 <- object$theta01[,id]
        
        if(s<min(knots01,knots02,knots12))stop(paste0("The argument s must be at least equal to :",min(knots01,knots02,knots12)))
        if(t<min(knots01,knots02,knots12))stop(paste0("The argument t must be at least equal to :",min(knots01,knots02,knots12)))
        
        names(the01) <- paste0("the01.",1:length(the01))
        the02 <- object$theta02[,id]
        names(the02) <- paste0("the02.",1:length(the02))
        the12 <- object$theta12[,id]
        names(the12) <- paste0("the12.",1:length(the12))
        
        #browser()
        
        knots.unique<-unique(object$knots01)
        knots.bound<-knots.unique[c(1,length(knots.unique))]
        knots.int<-knots.unique[-c(1,length(knots.unique))]
        msplines01<-splinesMI(x=t,knots=sort(knots.int),Boundary.knots=knots.bound)$Mspline
        isplines01<-(splinesMI(x=t,knots=sort(knots.int),Boundary.knots=knots.bound)$Ispline-splinesMI(x=s,knots=sort(knots.int),Boundary.knots=knots.bound)$Ispline)
        
        knots.unique<-unique(object$knots02)
        knots.bound<-knots.unique[c(1,length(knots.unique))]
        knots.int<-knots.unique[-c(1,length(knots.unique))]
        msplines02<-splinesMI(x=t,knots=sort(knots.int),Boundary.knots=knots.bound)$Mspline
        isplines02<-(splinesMI(x=t,knots=sort(knots.int),Boundary.knots=knots.bound)$Ispline-splinesMI(x=s,knots=sort(knots.int),Boundary.knots=knots.bound)$Ispline)
        
        knots.unique<-unique(object$knots12)
        knots.bound<-knots.unique[c(1,length(knots.unique))]
        knots.int<-knots.unique[-c(1,length(knots.unique))]
        msplines12<-splinesMI(x=t,knots=sort(knots.int),Boundary.knots=knots.bound)$Mspline
        isplines12<-(splinesMI(x=t,knots=sort(knots.int),Boundary.knots=knots.bound)$Ispline-splinesMI(x=s,knots=sort(knots.int),Boundary.knots=knots.bound)$Ispline)
        
          
        theta.square01<-the01^2
        intensity01<-msplines01%*%theta.square01
        cumulative.intensity01<-isplines01%*%theta.square01
        
        theta.square02<-the02^2
        intensity02<-msplines02%*%theta.square02
        cumulative.intensity02<-isplines02%*%theta.square02
        
        theta.square12<-the12^2
        intensity12<-msplines12%*%theta.square12
        cumulative.intensity12<-isplines12%*%theta.square12
        
        if (!is.null(beta01))
          linPred01<-beta01 %*% t(Z01)
        else
          linPred01<-0
        if (!is.null(beta02))
          linPred02<-beta02 %*% t(Z02)
        else
          linPred02<-0
        if (!is.null(beta12))
          linPred12<-beta12 %*% t(Z12)
        else
          linPred12<-0
          
          
        e01 <- exp(linPred01)
        intensity01<-intensity01*e01
        cumulative.intensity01<-cumulative.intensity01*e01
        
        e02 <- exp(linPred02)
        intensity02<-intensity02*e02
        cumulative.intensity02<-cumulative.intensity02*e02
        
        e12 <- exp(linPred12)
        intensity12<-intensity12*e12
        cumulative.intensity12<-cumulative.intensity12*e12
        
        intensity<-c(intensity01,intensity02,intensity12)
        cumulative.intensity<-c(cumulative.intensity01,cumulative.intensity02,cumulative.intensity12)

        if (lifeExpect==TRUE){
          transprob <- unlist(lifexpect0.idmPl(s,
                                               knots01,
                                               nknots01,
                                               the01,
                                               knots12,
                                               nknots12,
                                               the12,
                                               knots02,
                                               nknots02,
                                               the02,
                                               bZ01,
                                               bZ12,
                                               bZ02))
        }else{
          transprob <- unlist(Predict0.idmPl(s,t,knots01,nknots01,the01,knots12,nknots12,the12,knots02,nknots02,the02,bZ01,bZ12,bZ02))
        }
      }else {
        ############################ weibull #####################################
        a01 <- object$modelPar[1,id]
        b01 <- object$modelPar[2,id]
        a02 <- object$modelPar[3,id]
        b02 <- object$modelPar[4,id]
        a12 <- object$modelPar[5,id]
        b12 <- object$modelPar[6,id]
        # want to calculate variability based on variance-covariance matrix
        modelPar01<-object$modelPar[1:2,id]^2
        intensity01<-modelPar01[1]*(modelPar01[2]^modelPar01[1])*t^(modelPar01[1]-1)
        cumulative.intensity01<-(modelPar01[2]*t)^modelPar01[1]
        
        modelPar02<-object$modelPar[3:4,id]^2
        intensity02<-modelPar02[1]*(modelPar02[2]^modelPar02[1])*t^(modelPar02[1]-1)
        cumulative.intensity02<-(modelPar02[2]*t)^modelPar02[1]
        
        modelPar12<-object$modelPar[5:6]^2
        intensity12<-modelPar12[1]*(modelPar12[2]^modelPar12[1])*t^(modelPar12[1]-1)
        cumulative.intensity12<-(modelPar12[2]*t)^modelPar12[1]
        
        if (!is.null(beta01))
          linPred01<-beta01 %*% t(Z01)
        else
          linPred01<-0
        if (!is.null(beta02))
          linPred02<-beta02 %*% t(Z02)
        else
          linPred02<-0
        if (!is.null(beta12))
          linPred12<-beta12 %*% t(Z12)
        else
          linPred12<-0
        
        
        e01 <- exp(linPred01)
        intensity01<-intensity01*e01
        cumulative.intensity01<-cumulative.intensity01*e01
        
        e02 <- exp(linPred02)
        intensity02<-intensity02*e02
        cumulative.intensity02<-cumulative.intensity02*e02
        
        e12 <- exp(linPred12)
        intensity12<-intensity12*e12
        cumulative.intensity12<-cumulative.intensity12*e12
        
        intensity<-c(intensity01,intensity02,intensity12)
        cumulative.intensity<-c(cumulative.intensity01,cumulative.intensity02,cumulative.intensity12)
    
        
        if (lifeExpect==TRUE){
          transprob <- unlist(lifexpect0.idmWeib(s,
                                                 a01^2,
                                                 1/b01^2,
                                                 a02^2,
                                                 1/b02^2,
                                                 a12^2,
                                                 1/b12^2,
                                                 bZ01,
                                                 bZ02,
                                                 bZ12,max=maxtime))
        }else{
          transprob <- unlist(Predict0.idmWeib(s,
                                               t,
                                               a01^2,
                                               1/b01^2,
                                               a02^2,
                                               1/b02^2,
                                               a12^2,
                                               1/b12^2,
                                               bZ01,
                                               bZ02,
                                               bZ12))
        }
      }
      
     
      transprob <- data.frame(cbind(transprob))
      names(transprob) <- c("Estimate")
      transprob <- cbind("Parameter"=rownames(transprob),transprob)
      rownames(transprob) <- NULL
      uppercumulative.intensity<-lowercumulative.intensity<-NULL
      upperintensity<-lowerintensity<-NULL
    }
    
    out <- list(transprob=transprob)
    out <- c(out,list(newdata=newdata))
    out <- c(out,list(s=s,t=t,conf.int=ifelse(do.conf.int,conf.int,FALSE)))
    out<- c(out,list(uppercumulative.intensity= uppercumulative.intensity,
                     lowercumulative.intensity= lowercumulative.intensity,
                     cumulative.intensity=cumulative.intensity,
                     upperintensity= upperintensity,
                     lowerintensity=lowerintensity,
                     intensity=intensity))
    class(out) <- "predict.idm"
    out
}

## prediction indicators with splines 
Predict0.idmPl <- function(s,t,knots01,nknots01,the01,knots12,nknots12,the12,knots02,nknots02,the02,bZ01=0,bZ12=0,bZ02=0) {
    if (s>(min(knots01[nknots01+6],knots02[nknots02+6],knots12[nknots12+6]))) {stop("argument s is off")}    
    if (any(t>knots12[nknots12+6])) {stop("argument t is off")}
    if (any(s<knots01[1])) {stop("argument s is off")}
    p11 <- S.pl(s,t,knots12,nknots12,the12,bZ12)
    p12 <- 1-p11
    p00 <- S.pl(s,t,knots01,nknots01,the01,bZ01)*S.pl(s,t,knots02,nknots02,the02,bZ02)
    
    if(t==Inf){
      p02_0 <- sapply(t,function(t) {integrate(f=function(x){
        S.pl(s,x,knots01,nknots01,the01,bZ01)*S.pl(s,x,knots02,nknots02,the02,bZ02)*intensity(times=x,knots=knots02,number.knots=nknots02,theta=the02,linear.predictor=bZ02,method="splines")$intensity},lower=s,upper=t)$value})
      p01 <- sapply(t,function(t) {integrate(f=function(x){S.pl(s,x,knots01,nknots01,the01,bZ01)*S.pl(s,x,knots02,nknots02,the02,bZ02)*intensity(times=x,knots=knots01,number.knots=nknots01,theta=the01,linear.predictor=bZ01,method="splines")$intensity*S.pl(x,t,knots12,nknots12,the12,bZ12)},lower=s,upper=t)$value})
      p02_1 <- 1-p00-p02_0-p01
      p02 <- p02_0+p02_1
      RM<- integrate(f = function(x) {S.pl(s, x, knots01, nknots01, the01, bZ01) * S.pl(s, x, knots02, nknots02, the02, bZ02) }, lower = s, upper = t)$value
      F01<- integrate(f = function(x) {S.pl(s, x, knots01, nknots01, the01, bZ01) * S.pl(s, x, knots02, nknots02, the02, bZ02) * intensity(times = x, knots = knots01, number.knots = nknots01, theta = the01,linear.predictor = bZ01,method="splines")$intensity }, lower = s, upper = t)$value

    }else{
      p02_0 <- sapply(t,function(t) {gauss_kronrod(f=function(x){
        S.pl(s,x,knots01,nknots01,the01,bZ01)*S.pl(s,x,knots02,nknots02,the02,bZ02)*intensity(times=x,knots=knots02,number.knots=nknots02,theta=the02,linear.predictor=bZ02,method="splines")$intensity},a=s,b=t)$value})
      p01 <- sapply(t,function(t) {gauss_kronrod(f=function(x){S.pl(s,x,knots01,nknots01,the01,bZ01)*S.pl(s,x,knots02,nknots02,the02,bZ02)*intensity(times=x,knots=knots01,number.knots=nknots01,theta=the01,linear.predictor=bZ01,method="splines")$intensity*S.pl(x,t,knots12,nknots12,the12,bZ12)},a=s,b=t)$value})
      p02_1 <- 1-p00-p02_0-p01
      p02 <- p02_0+p02_1
      RM<- gauss_kronrod(f = function(x) {S.pl(s, x, knots01, nknots01, the01, bZ01) * S.pl(s, x, knots02, nknots02, the02, bZ02) }, a = s, b = t)$value
      F01<- gauss_kronrod(f = function(x) {S.pl(s, x, knots01, nknots01, the01, bZ01) * S.pl(s, x, knots02, nknots02, the02, bZ02) * intensity(times = x, knots = knots01, number.knots = nknots01, theta = the01,linear.predictor = bZ01,method="splines")$intensity }, a = s, b = t)$value

    }
    list(p00=p00,p01=p01,p11=p11,p12=p12,p02_0=p02_0,p02_1=p02_1,p02=p02,F01=F01,F0.=p02_0+p01+p02_1, RM=RM)
}

## prediction indicators with weibull 
Predict0.idmWeib <- function(s,t,a01,b01,a02,b02,a12,b12,bZ01=0,bZ02=0,bZ12=0) {

  if(t==Inf){
    p11 = S.weib(s,t,a12,b12,bZ12)
    p12 = 1-p11
    p00 = S.weib(s,t,a01,b01,bZ01)*S.weib(s,t,a02,b02,bZ02)
    p02_0 = sapply(t,function(t) {integrate(f=function(x){S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)*iweibull(x,a02,b02,bZ02)},lower=s,upper=t)$value })
    p01 = sapply(t,function(t) {integrate(f=function(x){S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)*iweibull(x,a01,b01,bZ01)*S.weib(x,t,a12,b12,bZ12)},lower=s,upper=t)$value})
    p02_1 = 1-p00-p02_0-p01
    p02 = p02_0+p02_1
    RM= integrate(f = function(x) {S.weib(s, x, a01, b01, bZ01) * S.weib(s, x, a02, b02, bZ02)}, lower = s, upper = t)$value
    F01= integrate(f = function(x) {S.weib(s, x, a01, b01, bZ01) * S.weib(s, x, a02, b02, bZ02)*iweibull(x,a01,b01,bZ01)}, lower = s, upper = t)$value
  }else{
    p11 = S.weib(s,t,a12,b12,bZ12)
    p12 = 1-p11
    p00 = S.weib(s,t,a01,b01,bZ01)*S.weib(s,t,a02,b02,bZ02)
    p02_0 = sapply(t,function(t) {gauss_kronrod(f=function(x){S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)*iweibull(x,a02,b02,bZ02)},a=s,b=t)$value })
    p01 = sapply(t,function(t) {gauss_kronrod(f=function(x){S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)*iweibull(x,a01,b01,bZ01)*S.weib(x,t,a12,b12,bZ12)},a=s,b=t)$value})
    p02_1 = 1-p00-p02_0-p01
    p02 = p02_0+p02_1
    RM= gauss_kronrod(f = function(x) {S.weib(s, x, a01, b01, bZ01) * S.weib(s, x, a02, b02, bZ02)},a=s,b=t)$value
    F01=  gauss_kronrod(f = function(x) {S.weib(s, x, a01, b01, bZ01) * S.weib(s, x, a02, b02, bZ02)*iweibull(x,a01,b01,bZ01)}, a = s, b = t)$value
    
  }
  list(p00=p00,p01=p01,p11=p11,p12=p12,p02_0=p02_0,p02_1=p02_1,p02=p02,F01=F01,F0.=p02_0+p01+p02_1, RM=RM)
}

# a = shape parameter
# b = scale parameter

iweibull <- function(x,a,b,bZ=0) {
    res = (a/b) * (x/b)**(a-1) * exp(bZ)
    return(res) 
}

# S(s,t) = S(t)/S(s)
# if S(s)=0, S(s,t)=0
# survival weibull
S.weib <- function(s,t,a,b,bZ=0) {	
    res <- 0
    St <- (1-pweibull(t,shape=a,scale=b))^(exp(bZ))
    Ss <- (1-pweibull(s,shape=a,scale=b))^(exp(bZ))
    if (length(s)==1){
        if (Ss==0){res <- 0}
        else{res <- St/Ss}
    }else{
         idx0 <- which(Ss==0)
         idx <- which(Ss!=0)
         res[idx0] <- 0 
         res[idx] <- St/Ss[idx]
     }
    return(res)
}


# Ok for new version
A <- function(s,t,zi,nknots,the,bZ=0) {
    res=rep(0,length(t))
    TF = (t>=zi[length(zi)])
    ind = which(TF)
    if (sum(TF)!=0) {res[ind]=intensity(
      times=(zi[nknots+6]-10^-5),
      knots=zi,
      number.knots=nknots,
      method="splines",
      theta=the,
      linear.predictor=bZ)$cumul.intensity-intensity(times=s,
                                                     knots=zi,
                                                     number.knots=nknots,
                                                     method="splines",
                                                     theta=the,
                                                     linear.predictor=bZ)$cumul.intensity}
    TF = (t<zi[length(zi)])
    ind = which(TF)
    if (sum(TF)!=0) {res[ind]=intensity(times=t[ind],
                                        knots=zi,
                                        number.knots=nknots,
                                        method="splines",
                                        theta=the,
                                        linear.predictor=bZ)$cumul.intensity-intensity(times=s,
                                                                                       knots=zi,
                                                                                       number.knots=nknots,
                                                                                       method="splines",
                                                                                       theta=the,
                                                                                       linear.predictor=bZ)$cumul.intensity}
    return(res)
}

### Survival function with two time s, t for splines
# S(s,t) = S(t)/S(s)
#        = exp(-A(s,t))

S.pl <- function(s,t,zi,nknots,the,bZ=0) {
    if (length(t)>=length(s)){
        res=rep(0,length(t))
        TF = (t>zi[length(zi)])
        ind = which(TF)
        if (sum(TF)!=0) {res[ind]=0}
        TF = (t<=zi[length(zi)])
        ind = which(TF)
        if (sum(TF)!=0) {res[ind]=intensity(times=t[ind],
                                            knots=zi,
                                            number.knots=nknots,
                                            method="splines",
                                            theta=the,
                                            linear.predictor=bZ)$survival/intensity(times=s,
                                                                                    knots=zi,
                                                                                    number.knots=nknots,
                                                                                    method="splines",
                                                                                    theta=the,
                                                                                    linear.predictor=bZ)$survival}
    }else{		
         res=rep(0,length(s))
         if (t>zi[length(zi)]) {res=0}
         else {res=intensity(times=t,
                             knots=zi,
                             number.knots=nknots,
                             method="splines",
                             theta=the,
                             linear.predictor=bZ)$survival/intensity(times=s,
                                                                     knots=zi,
                                                                     number.knots=nknots,
                                                                     method="splines",
                                                                     theta=the,
                                                                     linear.predictor=bZ)$survival}
     }
    return(res)
}


# Life expectency with weibull 
lifexpect0.idmWeib <- function(s,a01,b01,a02,b02,a12,b12,bZ01=0,bZ02=0,bZ12=0,max) {

  if(max==Inf){
    ET12 = integrate(
        f=function(x) {
            S.weib(s,x,a12,b12,bZ12)
        },s,max)
    ET0dot = integrate(f=function(x) {
        S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)
    },s,max)
    ET01 = integrate(f=function(x){
        sapply(x,function(x){
            integrate(f=function(y){
                S.weib(s,y,a01,b01,bZ01)*S.weib(s,y,a02,b02,bZ02)*iweibull(y,a01,b01,bZ01)*S.weib(y,x,a12,b12,bZ12)},
                lower=s,
                upper=x)$value})},s,max)
    LTR=integrate(f=function(x){
      S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)*iweibull(x,a01,b01,bZ01)},s,max)
  }else{
    ET12 = gauss_kronrod(
      f=function(x) {
        S.weib(s,x,a12,b12,bZ12)
      },a=s,b=max)
    ET0dot = gauss_kronrod(f=function(x) {
      S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)
    },a=s,b=max)
    ET01 = gauss_kronrod(f=function(x){
      sapply(x,function(x){
        gauss_kronrod(f=function(y){
          S.weib(s,y,a01,b01,bZ01)*S.weib(s,y,a02,b02,bZ02)*iweibull(y,a01,b01,bZ01)*S.weib(y,x,a12,b12,bZ12)},
          a=s,
          b=x)$value})},a=s+0.0001,b=max) # as a < b 
    LTR=gauss_kronrod(f=function(x){
      S.weib(s,x,a01,b01,bZ01)*S.weib(s,x,a02,b02,bZ02)*iweibull(x,a01,b01,bZ01)},a=s,b=max)
  }
    list(LE.00=ET0dot$value,
         LE.0.=ET01$value+ET0dot$value,
         LE.01=ET01$value,
         LE.11=ET12$value,
         LTR=LTR$value)

}
# Life expectency with splines
lifexpect0.idmPl <- function(s,knots01,nknots01,the01,knots12,nknots12,the12,knots02,nknots02,the02,bZ01=0,bZ12=0,bZ02=0) {
  if(any(c(knots12[nknots12+6],knots02[nknots02+6],knots01[nknots01+6])==Inf)){
    ET12 = integrate(f=function(x) {
      S.pl(s,x,knots12,nknots12,the12,bZ12)},s,knots12[nknots12+6])
    ET0dot = integrate(f=function(x) {
      S.pl(s,x,knots01,nknots01,the01,bZ01)*S.pl(s,x,knots02,nknots02,the02,bZ02)  },s,knots02[nknots02+6])
    ET01 = integrate(f=function(x) {
      sapply(x,function(x) {integrate(f=function(y){
        (S.pl(s,y,knots01,nknots01,the01,bZ01)
         *S.pl(s,y,knots02,nknots02,the02,bZ02)*
           intensity(times=y,
                     knots=knots01,
                     number.knots=nknots01,
                     theta=the01,
                     linear.predictor=bZ01,
                     method="splines")$intensity
         *S.pl(y,x,knots12,nknots12,the12,bZ12))},
        lower=s,upper=x)$value})},s,knots01[nknots01+6])
    LTR=integrate(f = function(x) {S.pl(s, x, knots01, nknots01, the01, bZ01) * S.pl(s, x, knots02, nknots02, the02, bZ02) * intensity(times = x,
                                                                                                                                       knots = knots01,
                                                                                                                                       number.knots = nknots01,
                                                                                                                                       theta = the01,
                                                                                                                                       linear.predictor = bZ01,
                                                                                                                                       method="splines")$intensity }, lower = s, upper = knots01[nknots01+6])$value
  }else{
    ET12 = gauss_kronrod(f=function(x) {
      S.pl(s,x,knots12,nknots12,the12,bZ12)},a=s,b=knots12[nknots12+6])
    ET0dot = gauss_kronrod(f=function(x) {
      S.pl(s,x,knots01,nknots01,the01,bZ01)*S.pl(s,x,knots02,nknots02,the02,bZ02)  },a=s,b=knots02[nknots02+6])
    ET01 = gauss_kronrod(f=function(x) {
      sapply(x,function(x) {gauss_kronrod(f=function(y){
        (S.pl(s,y,knots01,nknots01,the01,bZ01)
         *S.pl(s,y,knots02,nknots02,the02,bZ02)*
           intensity(times=y,
                     knots=knots01,
                     number.knots=nknots01,
                     theta=the01,
                     linear.predictor=bZ01,
                     method="splines")$intensity
         *S.pl(y,x,knots12,nknots12,the12,bZ12))},
        a=s,b=x)$value})},a=s+0.0001,b=knots01[nknots01+6])# as a < b 
    LTR=gauss_kronrod(f = function(x) {S.pl(s, x, knots01, nknots01, the01, bZ01) * S.pl(s, x, knots02, nknots02, the02, bZ02) * intensity(times = x,
                                                                                                                                       knots = knots01,
                                                                                                                                       number.knots = nknots01,
                                                                                                                                       theta = the01,
                                                                                                                                       linear.predictor = bZ01,
                                                                                                                                       method="splines")$intensity }, a = s, b = knots01[nknots01+6])$value
    
  }
    list(LE.00=ET0dot$value,
       LE.0.=ET01$value+ET0dot$value,
       LE.01=ET01$value,
       LE.11=ET12$value,
       LTR=LTR)
}

# 
#  data(Paq1000)
# library(prodlim)
#  library(SmoothHazardoptim9)
#  fit.paq <- SmoothHazardoptim9::idm(formula02=Hist(time=t,event=death,entry=e)~certif:gender,
#                 formula01=Hist(time=list(l,r),event=dementia)~certif:gender,data=Paq1000)
# 
#  p1<-predict(fit.paq,
#              s=70,t=80,lifeExpect=F,conf.int=T,
#              newdata=data.frame(certif=1,gender=1),nsim=2)
#  p1<-predict(fit.paq,
#              s=70,t=80,lifeExpect=T,conf.int=T,
#              newdata=data.frame(certif=1,gender=0))
#  p1bis<-predict(fit.paq,s=70,t=80,lifeExpect=F,conf.int=T,newdata=data.frame(certif=1,gender=1),
#                 nsim = 1)
#  p1bis<-predict(fit.paq,s=70,t=80,lifeExpect=T,conf.int=T,newdata=data.frame(certif=1,gender=1),
#                 nsim = 1)
# 
#  library(SmoothHazardoptim9)
#  fit.paq2 <- SmoothHazardoptim9::idm(formula02=Hist(time=t,event=death,entry=e)~certif:gender,
#                 formula01=Hist(time=list(l,r),event=dementia)~certif,data=Paq1000,method="splines")
# 
#  p1<-predict(fit.paq2,s=70,t=80,lifeExpect=F,newdata=data.frame(certif=1,gender=1))
#  p1bis<-predict(fit.paq2,s=70,t=80,lifeExpect=F,newdata=data.frame(certif=1),nsim=1)

#  p1<-predict(fit.paq2,s=70,t=80,lifeExpect=T,newdata=data.frame(certif=1))
#  p1bis<-predict(fit.paq2,s=70,t=80,lifeExpect=T,newdata=data.frame(certif=1),nsim=1)
 # pbr with print 
# 
# p3<- SmoothHazardoptim9::predict.idm(fit.splines,s=70,t=80,newdata=data.frame(certif=1))
# predict(fit.splines,s=70,t=80,lifeExpect=TRUE,newdata=data.frame(certif=1),nsim=20)
# 
 
 # pspline<-SmoothHazardoptim9::  idm(formula02 = f02,
 #                           formula01 = f01,
 #                           formula12 = f12,data=simu,
 #                           nproc=1,penalty = penalty,eps=c(7,2,1),
 #                           method="splines",scale.X=F,
 #                           lambda01 = lambda01,
 #                           lambda02=lambda02,
 #                           lambda12=lambda12,
 #                           alpha=k,maxiter=100,maxiter.pena = 10)
 # 
 # for(k in 1:dim(simu)[1]){
 # print(k)
 # p1<-predict(fit.idm.C,
 #             s=2,t=10,lifeExpect=T,conf.int=T,
 #             newdata=simu[k,])
 # p1bis<-predict(fit.idm.C,
 #             s=2,t=10,lifeExpect=F,conf.int=T,
 #             newdata=simu[k,])
 # p1<-predict( pspline,
 #             s=2,t=10,lifeExpect=T,conf.int=T,
 #             newdata=simu[k,])
 # p1bis<-predict( pspline,
 #                s=2,t=10,lifeExpect=F,conf.int=T,
 #                newdata=simu[k,])
 # }
