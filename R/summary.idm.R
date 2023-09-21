#' Summary of a fitted illness-death model
#' 
#' Summarize the event history data of an illness-death regression model
#' and show regression coefficients for transition intensities
#' 
#' @param object a \code{idmSplines} object, i.e., the result of a call to the
#' \code{\link{idm}} function with \code{intensities}="Splines".
#' @param conf.int  The level of confidence for the hazard ratios. The default is \code{0.95}.
#' @param digits number of digits to print.
#' @param pvalDigits number of digits to print for p-values.
#' @param eps convergence criterion used for p-values.
#' @param \dots other unusued arguments.
#' @author R: Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr> Fortran:
#' Pierre Joly <Pierre.Joly@@isped.u-bordeaux2.fr>
#' @seealso \code{\link{idm}}, \code{\link{print.idm}},
#' \code{\link{plot.idm}} 
#' @keywords methods
#' @examples
#' 
#' \dontrun{
#' library(prodlim)
#' data(Paq1000)
#' fit.splines <-  idm(formula02=Hist(time=t,event=death,entry=e)~certif,
#' 		formula01=Hist(time=list(l,r),event=dementia)~certif,
#'                 formula12=~1,
#'                 method="splines",
#' 		data=Paq1000)
#' summary(fit.splines) 
#' }
#' @useDynLib SmoothHazardoptim9
#' @export
summary.idm <- function(object,conf.int=.95,digits=4,pvalDigits=4,eps=.0001, ...){
    if (!inherits(object,"idm")) stop("Object must be of class 'idm'")
  
    if(!is.null(object$modelPar)){
      object$method<-"weib"
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
            n_base<-ifelse(object$method=="weib",7,object$nknots01+object$nknots02+object$nknots12+6+1)
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
    
    cat("Method:",switch(object$method,
                         "splines"="M-splines based on likelihood",
                         "Weib"="Weibull parametrization"),"\n")
    cat("\n")
    cat("number of subjects: ", object$N,"\n")
    cat("number of events '0-->1': ", object$events1,"\n")
    cat("number of events '0-->2 or 0-->1-->2': ", object$events2,"\n")
    cat("number of covariates: ", object$NC,"\n")
    if(length(object$na.action))cat("observation deleted due to missing: ",length(object$na.action),"\n")
    
  }
}
