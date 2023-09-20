#' Print method for \code{idm} objects
#' 
#' Print a summary of a fitted illness-death model 
#' 
#' 
#' @param x Class \code{idm} object, i.e. the result of a call to the
#' \code{\link{idm}} function with \code{intensities}="Weib".
#' @param conf.int The level of confidence for the hazard ratios. The default is \code{0.95}.
#' @param digits Number of digits to print.
#' @param pvalDigits Number of digits to print for p-values.
#' @param eps Passed to \code{format.pval}.
#' @param \dots Not used.
#' @author Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr>, Thomas A. Gerds <tag@@biostat.ku.dk> 
#' @seealso \code{\link{summary.idm}}, \code{\link{plot.idm}}
#' @keywords methods
#' @examples
#' 
#' \dontrun{
#' data(Paq1000)
#' library(prodlim)
#' fit.splines <-  idm(formula02=Hist(time=t,event=death,entry=t0)~certif,
#' 		formula01=Hist(time=list(l,r),event=dementia)~certif,
#'                 formula12=~1,
#'                 method="Splines",
#' 		data=Paq1000)
#' print(fit.splines)
#' 
#' }
#' @useDynLib SmoothHazardoptim9
#' @export
print.idm <- function(x,conf.int=.95,digits=4,pvalDigits=4,eps=0.0001,...){
    # {{{  call
    cl <- x$call
    if(is.null(x$BIC)){
      n_model<-1
    }else{
      n_model<-length(x$BIC)
    }
    
    if(!is.null(x$modelPar)){
     method<-"weib"}else{method<-"splines"}
    
    
    cat("Call:\n")
    dput(cl)
    cat("\n")
    # }}}
    # {{{ number of subjects etc
    cat("Illness-death regression model using",
        ifelse(method=="splines"," M-spline approximations","Weibull parametrization"),
        "\nto estimate the baseline transition intensities.\n")
    cat("\n")
    cat("number of subjects: ", x$N,"\n")
    cat("number of events '0 -> 1': ", x$events1,"\n")
    cat("number of events '0 -> 2' or '0 -> 1 -> 2': ", x$events2,"\n")
    cat("number of covariates: ", x$NC,"\n")
    # }}}
    # {{{ convergence
    # FIXME: what is the difference between maximum number of iterations reached and
    #        model did not converge?
    
    if( sum(x$converged!=1)>0 ){
      if(length(x$converged)==1){
        warning("The model did not converge.","\n")
        switch(as.character(x$converged[1]),
               "2"={ warning("Maximum number of iterations reached.",call.=FALSE)},
               "3"={ warning("Fisher information matrix non-positive definite.",call.=FALSE)})
      }else{
        warning(paste0("Among the ",n_model," models, ",sum(x$converged!=1)," did not converged"),"\n")
        warning(paste0("Maximum number of iterations reached for ",sum(x$converged==2)," models"))}
    }else{
        if(length(x$converged)==1 ){
            warning("The model did converge.","\n")
            cat("Log-likelihood : ",x$loglik[1], "\n")
            cat("----\nModel converged.\n")
            cat("number of iterations: ", x$niter,"\n")
            cat("convergence criteria: parameters=", signif(x$cv$cb,2), "\n")
            cat("                    : likelihood=", signif(x$cv$ca,2), "\n") 
            cat("                    : second derivatives=", signif(x$cv$rdm,2), "\n")
  
        }else{
          warning(paste0("All ",n_model," models did converge."),"\n")
          cat("Log-likelihood : ",x$loglik[1:n_model], "\n")
          cat("----\nModel converged.\n")
          cat("number of iterations: ", x$niter,"\n")
          if(x$penalty=="none"){
            
          cat("convergence criteria: parameters=", signif(x$cv$cb,2), "\n")
          cat("                    : likelihood=", signif(x$cv$ca,2), "\n") 
          cat("                    : second derivatives=", signif(x$cv$rdm,2), "\n")
          }else{
            ca.beta<-apply(m01$cv$ca.beta,MARGIN=2,FUN=function(x){x<-na.omit(x) 
            return(x[length(x)])})
            ca.spline<-apply(m01$cv$ca.spline,MARGIN=2,FUN=function(x){x<-na.omit(x) 
            return(x[length(x)])})
            cb<-apply(m01$cv$cb,MARGIN=2,FUN=function(x){x<-na.omit(x) 
            return(abs(x[length(x)]-x[length(x)-1]))})
            cat("convergence criteria: parameters beta=", ca.beta, "\n")
            cat("                    : parameters base risk=", ca.spline, "\n") 
            cat("                    : likelihood=", cb, "\n")
          }
        }
          
          
    # }}}
    # {{{ Spline: baseline parameters
    if (method=="splines"){
        n_spline<-x$nknots01+x$nknots02+x$nknots12+6
        splinepars <- data.frame("transition01"=x$nknots01,
                                 "transition02"=x$nknots02,
                                 "transition12"=x$nknots12)
        rownames(splinepars) <- c("knots")
        cat("\n")
        cat("Splines parameters:\n")
        print(splinepars,row.names=TRUE)
    
        
        # }}}
        # {{{ Weibull: baseline parameters
    }
        # }}}
        # {{{ Weibull: baseline parameter
    if (method=="weib"){
      n_spline<-6
        cat("Parameters of the Weibull distributions: 'S(t) = exp(-(b*t)^a)'\n")
      if(n_model==1){
        wpars <- matrix(x$modelPar,nrow=2)
        dimnames(wpars) <- list(c("shape (a)","scale (b)"),
                                c("transition 0 -> 1",
                                  "transition 0 -> 2",
                                  "transition 1 -> 2"))
        prmatrix(wpars)
      }else{
        rownames(x$modelPar)<-rep(c("shape (a)","scale (b)"),3)
        cat("transition 0 ->1 : \n")
        print(x$modelPar[1:2,],row.names=TRUE)
        cat("\n")
        cat("transition 0 ->2 : \n")
        print(x$modelPar[3:4,],row.names=TRUE)
        cat("\n")
        cat("transition 1 ->2 : \n")
        print(x$modelPar[5:6,],row.names=TRUE)
        cat("\n")
        
      }
    }
    # }}}
    # {{{  Regression coefficients
    if(sum(x$NC)>0){
      
      if(n_model==1){
        
        se<-sqrt(diag(x$V[(n_spline+1):(dim(x$V)[1]),(n_spline+1):(dim(x$V)[1])]))
        
      }else{
        se<-matrix(NA,ncol=n_model,nrow=dim(x$V)[1])
        for(k in 1:n_model){
          se[,k]<-sqrt(diag(x$V[1:(dim(x$V)[1]),((k-1)*dim(x$V)[1]+1):(k*dim(x$V)[1])]))
        }
      }
        
        wald <- (x$coef/se)**2
        z <- abs(qnorm((1 + conf.int)/2))
        coefmat <- data.frame("coef"=format(round(x$coef,digits)),
                              "SE coef"=format(round(se,digits)),
                              "exp(coef)"=format(round(x$HR,digits)),
                              "CI"=matrix(paste("[",
                                         format(round(exp(x$coef - z * se),2)),
                                         ";",
                                         format(round(exp(x$coef + z * se),2)),
                                         "]",
                                         sep=""),ncol=n_model,nrow=dim(x$coef)),
                              ## "Wald"=format(wald,digits),
                              "p-value"=matrix(format.pval(1 - pchisq(wald, 1),digits=pvalDigits,eps=eps),
                                               ,ncol=n_model,nrow=dim(x$coef)),
                              check.names=FALSE)
        coefmat <- cbind(Factor=ifelse(n_model==1,names(x$coef),
                                       rownames(x$coef)),coefmat)
        coeflist <- split(coefmat,rep(c("transition 0 -> 1","transition 0 -> 2","transition 1 -> 2"),x$NC))
        cat("\n\nRegression coefficients:\n\n")
        print(coeflist)
      }
       
    }
}

