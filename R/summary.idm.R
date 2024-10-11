### Code:
#' @title Summary of a fitted illness-death model
#' @description
#' Summarize the event history data of an illness-death regression model
#' and show regression coefficients for transition intensities
#' @param object a \code{idmSplines} object, i.e., the result of a call to the
#' \code{\link{idm}} function with \code{intensities}="Splines".
#' @param conf.int  The level of confidence for the hazard ratios. The default is \code{0.95}.
#' @param digits number of digits to print.
#' @param pvalDigits number of digits to print for p-values.
#' @param eps convergence criterion used for p-values.
#' @param \dots other unusued arguments.
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
#' @seealso \code{\link{idm}}, \code{\link{print.idm}},
#' \code{\link{plot.idm}} 
#' @keywords methods
#' @examples
#'
#' \dontrun{
#' library(lava)
#' library(prodlim)
#' set.seed(17)
#' d <- simulateIDM(n=1000,
#'                  beta01=c(1,1,0,0.5,0.5,rep(0,5)),
#'                  beta02=c(1,0,0,0,0.5,rep(0,5)),
#'                  beta12=c(1,0,0,0,0.5,rep(0,5)))$data
#'                  
#' fitpenweib <- idm(formula01=Hist(time=list(L,R),
#' event=seen.ill)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
#'                   formula02=Hist(time=observed.lifetime,
#'                   event=seen.exit)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
#'                   formula12=~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
#'                   data=d,penalty="lasso",lambda01 = c(10,20),
#'                   lambda02 = 10, lambda12 = 10)
#' summary(fitpenweib)
#' }
#' @useDynLib SmoothHazardoptim9
#' @export
summary.idm <- function(object,conf.int=.95,digits=4,pvalDigits=4,eps=.0001, ...){
    if (!inherits(object,"idm")) stop("Object must be of class 'idm'")
  
    if(!is.null(object$modelPar)){
      object$method<-"Weib"
    }else{object$method<-"splines"
    }
  
  if(object$penalty=="none"){
        cat("Method:",switch(object$method,
                             "splines"="M-splines based on likelihood",
                             "Weib"="Weibull parametrization"),"\n")
        cat("\n")
        cat("number of subjects: ", object$N,"\n")
        cat("number of events '0-->1': ", object$events1,"\n")
        cat("number of events '0-->2 or 0-->1-->2': ", object$events2,"\n")
        cat("number of covariates: ", object$NC,"\n")
        if(length(object$na.action))cat("observation deleted due to missing: ",length(object$na.action),"\n")
        if((sum(object$NC)>0)&&(object$converged==1)){
            n_base<-ifelse(object$method=="Weib",7,object$nknots01+object$nknots02+object$nknots12+6+1)
            se<-object$V[n_base:dim(object$V)[1],n_base:dim(object$V)[2]]
            se<-sqrt(diag(se))
            wald <- (object$coef/se)**2
            z <- abs(qnorm((1 + conf.int)/2))
            out <- data.frame("Hazard ratio"=format(round(exp(object$coef),digits)),
                              "Standard error"=format(round(se,digits)),
                              "CI"=paste("[",format(round(exp(object$coef - z * se),2)),";",format(round(exp(object$coef + z *se),2)),"]",sep=""),
                              "P-value"=format.pval(1 - pchisq(wald, 1),digits=pvalDigits,eps=eps))
            names(out)[3] <- paste("CI",round(100*conf.int),sep=".")
            Xnames <- NULL
            if(object$NC[1]>0) Xnames <- c(Xnames,paste(object$Xnames01,"_01",sep=""))
            if(object$NC[2]>0) Xnames <- c(Xnames,paste(object$Xnames02,"_02",sep=""))
            if(object$NC[3]>0) Xnames <- c(Xnames,paste(object$Xnames12,"_12",sep=""))
            rownames(out) <- Xnames
            print(out,row.names=T)
        }

  }else{
    BICmin<-min(object$BIC[object$converged==1])
    if(length(BICmin)==0){BICmin<-"No convergence"}else{id<-which(object$BIC==BICmin)}
    
    cat("Method:",switch(object$method,
                         "splines"="M-splines based on likelihood",
                         "Weib"="Weibull parametrization"),"\n")
    cat("\n")
    cat("number of subjects: ", object$N,"\n")
    cat("number of events '0-->1': ", object$events1,"\n")
    cat("number of events '0-->2 or 0-->1-->2': ", object$events2,"\n")
    cat("number of covariates: ", object$NC,"\n")
    cat("type of penalty: ", object$penalty,"\n")
    cat("range of lambda value '0-->1': ", min(object$lambda[1,]),";",max(object$lambda[1,]),"\n")
    cat("range of lambda value '0-->2': ", min(object$lambda[2,]),";",max(object$lambda[2,]),"\n")
    cat("range of lambda value '1-->2': ", min(object$lambda[3,]),";",max(object$lambda[3,]),"\n")
    if(length(object$na.action))cat("observation deleted due to missing: ",length(object$na.action),"\n \n")
    
    if(BICmin!="No convergence"){
      cat("Minimum BIC with convergence : BIC=", min(object$BIC[object$converged==1])," \n with lambda (",object$lambda[1,id],";",
          object$lambda[2,id],";",object$lambda[3,id],") and alpha =",object$alpha[id],"\n")
      
      out01 <- data.frame("Variable"=object$Xnames01,
                        "Hazard_ratio"=format(round(exp(object$coef[1:object$NC[1],id]),digits)),
                        "Selected"=ifelse(object$coef[1:object$NC[1],id]==0,"No","Yes"))
      out02 <- data.frame("Variable"=object$Xnames02,
                        "Hazard_ratio"=format(round(exp(object$coef[(object$NC[1]+1):(object$NC[1]+object$NC[2]),id]),digits)),
                        "Selected"=ifelse(object$coef[(object$NC[1]+1):(object$NC[1]+object$NC[2]),id]==0,"No","Yes"))
      out12 <- data.frame("Variable"=object$Xnames12,
                        "Hazard_ratio"=format(round(exp(object$coef[(object$NC[1]+object$NC[2]+1):(sum(object$NC)),id]),digits)),
                        "Selected"=ifelse(object$coef[(object$NC[1]+object$NC[2]+1):(sum(object$NC)),id]==0,"No","Yes"))
      cat("For transition '0-->1': \n \n")
      print(out01,row.names=F)
      cat(" \n \n For transition '0-->2': \n \n")
      print(out02,row.names=F)
      cat("\n \n For transition '1-->2': \n \n")
      print(out12,row.names=F)
    }else{cat("No convergence was obtained for all lambda values, consider increasing maxiter or/and changing values of lambda \n")
    }
  }
}
