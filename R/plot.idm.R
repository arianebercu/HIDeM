### Code:
#' @title Plot method for an illness-death model 
#' @description
#' Plot estimated baseline transition intensities from an object of class
#' \code{idm} optionally with confidence limits.
#' @param x a \code{idmWeib} class object (output from calling
#' \code{idm} with the (default) option \code{intensities}="Weib".
#' @param conf.int If TRUE show confidence limits
#' @param citype Type of confidence limits, can be "shadow" or "bars"
#' @param add If TRUE add to existing plot
#' @param axes If TRUE axes are drawn
#' @param col Color of the lines
#' @param lwd Width of the lines
#' @param lty Type of the lines
#' @param xlim Limits for x-axis
#' @param ylim Limits for y-axis
#' @param xlab Label for x-axis
#' @param ylab Label for y-axis
#' @param cex character expansion factor relative to current par("cex"). Used for text, and provides the default for pt.cex. 
#' @param y.intersp vertical (y) distances (in lines of text shared above/below each legend entry). A vector with one element for each row of the legend can be used.
#' @param lambda give the lambda for which you want to plot the penalised illness-death model
#' @param legend If TRUE a legend is drawn, which can be further controlled via \code{\link{SmartControl}}.
#' @param transition Choose one of the transition intensities: \code{c("01","02","12")}.
#' @param ... Passed to \code{\link{SmartControl}}
#' @return Print a plot of the baseline transition intensities of an
#' illness-death model estimated using a Weibull approach.
#' @seealso
#' \code{\link{print.idm}},\code{\link{summary.idm}},\code{\link{idm}},
#' @seealso \code{\link{idm}}
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
#' fitsplines <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2,
#' formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
#' formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
#' data=d,method="splines")
#' plot(fitsplines)
#' 
#' fitpenspline <- idm(formula01=Hist(time=list(L,R),
#' event=seen.ill)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
#'                     formula02=Hist(time=observed.lifetime,
#'                     event=seen.exit)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
#'                     formula12=~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
#'                     method="splines",
#'                     data=d,penalty="lasso",
#'                     lambda01 = c(10,20),
#'                     lambda02 = 10, lambda12 = 10)
#' plot(fitpenspline,lambda=c(10,10,10))
#' }
#'  
#'@useDynLib HIDeM
#' @export
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
plot.idm <- function(x,
                     lambda=NULL,
                     conf.int=FALSE,
                     citype="shadow",
                     add=FALSE,
                     axes=TRUE,
                     col,
                     lwd,
                     lty,
                     xlim,
                     ylim,
                     xlab,
                     ylab,
                     cex,
                     y.intersp,
                     legend=TRUE,
                     transition=c("01","02","12"),
                     ...){ 
  
  # {{{ collecting the (X, Y)-values of the lines
  
  
  if(any(x$converged!=1)){
    conf.int<-F
  }
  
  
  if(conf.int==T){
    conf<-0.95
  }else{conf<-F}
  
  if(x$penalty=="none"){
    
    if(x$converged!=1){
      stop("The model did not converged, plot can not be displayed")
    }
    if(!is.null(x$modelPar)){
      
      x$method<-"weib"
      
      X<-X01 <- X02 <- X12 <- seq(x$mintime,x$maxtime,length.out=100)
      
      
      intensity01<-intensity(times=X01,theta=x$modelPar[1:2]^2,
                             fix=x$fix[1:2],
                             conf.int=conf,
                             converged=x$converged,
                             V=x$V[1:2,1:2])
      x$intensity01<-intensity01$intensity
      x$lowerIntensity01<-intensity01$lowerintensity
      x$upperIntensity01<-intensity01$upperintensity
      
      intensity02<-intensity(times=X02,theta=x$modelPar[3:4]^2,
                             fix=x$fix[3:4],
                             conf.int=conf,
                             converged=x$converged,
                             V=x$V[3:4,3:4])
      x$intensity02<-intensity02$intensity
      x$lowerIntensity02<-intensity02$lowerintensity
      x$upperIntensity02<-intensity02$upperintensity
      
      intensity12<-intensity(times=X12,theta=x$modelPar[5:6]^2,
                             fix=x$fix[5:6],
                             conf.int=conf,
                             converged=x$converged,
                             V=x$V[5:6,5:6])
      x$intensity12<-intensity12$intensity
      x$lowerIntensity12<-intensity12$lowerintensity
      x$upperIntensity12<-intensity12$upperintensity
      
      
      
    }else{
      
      x$method<-"splines"
      
      X01 <-  seq(min(x$knots01),max(x$knots01),length.out=100)
      X02 <- seq(min(x$knots02),max(x$knots02),length.out=100)
      X12 <- seq(min(x$knots12),max(x$knots12),length.out=100)
      X<-unique(c(X01,X02,X12))
      
      intensity01<-intensity(times=X01,theta=x$theta01,
                             knots=x$knots01,
                             number.knots=x$nknots01,
                             fix=x$fix[1:(x$nknots01+2)],
                             conf.int=conf,
                             converged=x$converged,
                             V=x$V[1:(x$nknots01+2),1:(x$nknots01+2)],method="splines")
      
      x$intensity01<-intensity01$intensity
      x$lowerIntensity01<-intensity01$lowerintensity
      x$upperIntensity01<-intensity01$upperintensity
      
      
      intensity02<-intensity(times=X02,theta=x$theta02,
                             knots=x$knots02,
                             number.knots=x$nknots02,
                             fix=x$fix[(x$nknots01+3):(x$nknots01+x$nknots02+4)],
                             conf.int=conf,
                             converged=x$converged,
                             V=x$V[(x$nknots01+3):(x$nknots01+x$nknots02+4),(x$nknots01+3):(x$nknots01+x$nknots02+4)],
                             method="splines")
      
      x$intensity02<-intensity02$intensity
      x$lowerIntensity02<-intensity02$lowerintensity
      x$upperIntensity02<-intensity02$upperintensity
      
      intensity12<-intensity(times=X12,theta=x$theta12,
                             knots=x$knots12,
                             fix=x$fix[(x$nknots01+x$nknots02+5):(x$nknots01+x$nknots02+x$nknots12+6)],
                             conf.int=conf,
                             converged=x$converged,
                             V=x$V[(x$nknots01+x$nknots02+5):(x$nknots01+x$nknots02+x$nknots12+6),(x$nknots01+x$nknots02+5):(x$nknots01+x$nknots02+x$nknots12+6)],
                             number.knots=x$nknots12,method="splines")
      
      x$intensity12<-intensity12$intensity
      x$lowerIntensity12<-intensity12$lowerintensity
      x$upperIntensity12<-intensity12$upperintensity
      
      
      
    }
  }else{
    
    if(sum(x$converged==1)==0){
      stop("None of the model converged, no plot can be displayed")
    }
    
    
    if(is.null(lambda)){lambda<-"BIC"}
    if(length(lambda)==1){if(!lambda%in%c("GCV","BIC")){stop("Lambda need to be either a vector of three values (01,02 and 12) or BIC or GCV")}}
    if(!length(lambda)%in%c(1,3)){stop("Lambda need to be either a vector of three values (01,02 and 12) or BIC or GCV")}
    
    if(length(lambda)==3){
      if(!any(apply(x$lambda,FUN=function(x){sum(x==lambda)},MARGIN=2)==3)){stop("Lambda need to be either a vector of three values (01,02 and 12) from x$lambda")}
      id<-which(apply(x$lambda,FUN=function(x){sum(x==lambda)},MARGIN=2)==3)[1]
    }
    if(length(lambda)==1){
      if(lambda=="BIC"){
        BIC<-min(x$BIC[x$converged==1])
        id<-which(x$BIC==BIC)[1]
      }else{
        GCV<-min(x$GCV[x$converged==1])
        id<-which(x$GCV==GCV)[1]}
    }
    
    if(!is.null(x$modelPar)){
      
      x$method<-"weib"
      
      X<-X01 <- X02 <- X12 <- seq(x$mintime,x$maxtime,length.out=100)
      
      
      
      intensity01<-intensity(times=X01,theta=x$modelPar[1:2,id]^2,
                             fix=x$fix[1:2],
                             conf.int=F,
                             converged=x$converged[id])
      #V=x$V[(6*(id-1)+1):(6*(id-1)+2),(6*(id-1)+1):(6*(id-1)+2)])
      x$intensity01<-intensity01$intensity
      #x$lowerIntensity01<-intensity01$lowerintensity
      #x$upperIntensity01<-intensity01$upperintensity
      
      intensity02<-intensity(times=X02,theta=x$modelPar[3:4,id]^2,
                             fix=x$fix[3:4],
                             conf.int=F,
                             converged=x$converged[id])
      #V=x$V[(6*(id-1)+3):(6*(id-1)+4),(6*(id-1)+3):(6*(id-1)+4)])
      x$intensity02<-intensity02$intensity
      #x$lowerIntensity02<-intensity02$lowerintensity
      #x$upperIntensity02<-intensity02$upperintensity
      
      intensity12<-intensity(times=X12,theta=x$modelPar[5:6,id]^2,
                             fix=x$fix[5:6],
                             conf.int=F,
                             converged=x$converged[id])
      #V=x$V[(6*(id-1)+5):(6*(id-1)+6),(6*(id-1)+5):(6*(id-1)+6)])
      x$intensity12<-intensity12$intensity
      #x$lowerIntensity12<-intensity12$lowerintensity
      #x$upperIntensity12<-intensity12$upperintensity
      
      
    }else{
      
      x$method<-"splines"
      
      
      X01 <-  seq(min(x$knots01),max(x$knots01),length.out=100)
      X02 <- seq(min(x$knots02),max(x$knots02),length.out=100)
      X12 <- seq(min(x$knots12),max(x$knots12),length.out=100)
      X<-unique(c(X01,X02,X12))
      
      n_spline<-x$nknots01+2+x$nknots02+2+x$nknots12+2
      
      
      intensity01<-intensity(times=X01,theta=x$theta01[,id],
                             knots=x$knots01,
                             number.knots=x$nknots01,
                             fix=x$fix[1:(x$nknots01+2)],
                             conf.int=F,
                             converged=x$converged[id],
                             #V=x$V[(n_spline*(id-1)+1):(n_spline*(id-1)+x$nknots01+2),(n_spline*(id-1)+1):(n_spline*(id-1)+x$nknots01+2)],
                             method="splines")
      
      x$intensity01<-intensity01$intensity
      #x$lowerIntensity01<-intensity01$lowerintensity
      #x$upperIntensity01<-intensity01$upperintensity
      
      
      intensity02<-intensity(times=X02,theta=x$theta02[,id],
                             knots=x$knots02,
                             number.knots=x$nknots02,
                             fix=x$fix[(x$nknots01+3):(x$nknots01+x$nknots02+4)],
                             conf.int=F,
                             converged=x$converged[id],
                             #V=x$V[(n_spline*(id-1)+x$nknots01+3):(n_spline*(id-1)+x$nknots01+x$nknots02+4),(n_spline*(id-1)+x$nknots01+3):(n_spline*(id-1)+x$nknots01+x$nknots02+4)],
                             method="splines")
      
      x$intensity02<-intensity02$intensity
      #x$lowerIntensity02<-intensity02$lowerintensity
      #x$upperIntensity02<-intensity02$upperintensity
      
      intensity12<-intensity(times=X12,theta=x$theta12[,id],
                             knots=x$knots12,
                             fix=x$fix[(x$nknots01+x$nknots02+5):(n_spline)],
                             conf.int=F,
                             converged=x$converged[id],
                             #V=x$V[(n_spline*(id-1)+x$nknots01+x$nknots02+5):(n_spline*id),(n_spline*(id-1)+x$nknots01+x$nknots02+5):(n_spline*id)],
                             number.knots=x$nknots12,method="splines")
      
      x$intensity12<-intensity12$intensity
      #x$lowerIntensity12<-intensity12$lowerintensity
      #x$upperIntensity12<-intensity12$upperintensity
      
      
      
    }
    
  }
  
  
  
  if(sum(c("01","02","12") %in% transition)!=length(transition)){
    stop("'transition' must be a subset of 'c('01','02','12')'")}
  
  Y <- list("01"=x$intensity01,"02"=x$intensity02,"12"=x$intensity12)
  Y <- Y[transition]
  nlines <- length(Y)
  # }}}
  # {{{ setting default arguments for plot, axes, legend, confint 
  if (missing(ylab)) ylab <- "Transition intensity"
  if (missing(xlab)) xlab <- "Time"
  if (missing(xlim)) xlim <- c(min(X), max(X))
  if (missing(ylim)) ylim <- c(0, 1)
  if (missing(lwd)) lwd <- rep(3,nlines)
  if (missing(col)) col <- 1:nlines
  if (missing(lty)) lty <- rep(1, nlines)
  if (missing(cex)) cex<-1.5
  if (missing(y.intersp)) y.intersp<-0.5
  if (length(lwd) < nlines) lwd <- rep(lwd, nlines)
  if (length(lty) < nlines) lty <- rep(lty, nlines)
  if (length(col) < nlines) col <- rep(col, nlines)
  axis1.DefaultArgs <- list()
  axis2.DefaultArgs <- list(at=seq(0,1,.25),side=2,percent=TRUE)
  lines.DefaultArgs <- list(type="l")
  plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,xlab = xlab,ylab = ylab)
  legend.DefaultArgs <- list(legend=paste("Transition",names(Y)),lwd=lwd,col=col,lty=lty,cex=cex,bty="n",y.intersp=y.intersp,x="topleft")
  confint.DefaultArgs <- list(x=x,citype="shadow",times=X,density=55,col=col[1:nlines],lwd=rep(2,nlines),lty=rep(3,nlines))
  # }}}
  control <- prodlim::SmartControl(call=  list(...),
                                   keys=c("plot","lines","legend","confint","axis1","axis2"),
                                   ignore=c("x","transition","add","col","lty","lwd","ylim","xlim","xlab","ylab","legend","conf.int","axes"),
                                   defaults=list("plot"=plot.DefaultArgs,"lines"=lines.DefaultArgs,"legend"=legend.DefaultArgs,"confint"=confint.DefaultArgs,"axis1"=axis1.DefaultArgs,"axis2"=axis2.DefaultArgs),
                                   forced=list("plot"=list(axes=FALSE),"axis1"=list(side=1)),
                                   ignore.case=TRUE,
                                   replaceDefaults=FALSE,
                                   verbose=TRUE)
  # {{{  plot and backGround
  if (!add) {
    do.call("plot",control$plot)
  }
  # }}}
  # {{{  axes
  if (!add) {
    if (axes){
      do.call("axis",control$axis1)
      if (control$axis2$percent & is.null(control$axis2$labels))
        control$axis2$labels <- paste(100*control$axis2$at,"%")
      do.call("axis",control$axis2[-match("percent",names(control$axis2),nomatch=0)])
    }
  }
  # }}}
  # {{{confidence intervals
  
  nix <- lapply(1:nlines,function(i){
    ci.lower <- x[[paste("lowerIntensity",names(Y)[i],sep="")]]
    ci.upper <- x[[paste("upperIntensity",names(Y)[i],sep="")]]
    time <- switch(i, "1"= X01,"2"=X02,"3"=X12)
    switch(citype,
           "bars"={
             segments(x0=time,x1=time,y0=ci.lower,y1=ci.upper,lwd=lwd[i],col=col[i],lty=lty[i],...)
           },
           "shadow"={
             ccrgb=as.list(col2rgb(col[i],alpha=TRUE))
             names(ccrgb) <- c("red","green","blue","alpha")
             ccrgb$alpha=control$confint$density
             cc=do.call("rgb",c(ccrgb,list(max=255)))
             ttt <- time
             nt <- length(ttt)
             ttt <- c(ttt,ttt)
             uuu <- c(0,ci.upper[-nt],ci.upper)
             lll <- c(0,ci.lower[-nt],ci.lower)
             neworder <- order(ttt)
             uuu <- uuu[neworder]
             lll <- lll[neworder]
             ttt <- sort(ttt)
             polygon(x=c(ttt,rev(ttt)),y=c(lll,rev(uuu)),col=cc,border=NA)
           },{
             lines(x=time,ci.lower,type="l",lwd=lwd[i],col=col[i],lty=3,...)
             lines(x=time,ci.upper,type="l",lwd=lwd[i],col=col[i],lty=3,...)
           })
  })
  # }}}
  # {{{  adding the lines
  lines.type <- control$lines$type
  nix <- lapply(1:nlines, function(s) {
    time <- switch(s, "1"= X01,"2"=X02,"3"=X12)
    lines(x = time,y = Y[[s]],type = lines.type,col = col[s],lty = lty[s],lwd = lwd[s])
  })
  # {{{  legend
  if(legend==TRUE && !add && !is.null(names(Y))){
    if (is.null(control$legend$title)){
      if (x$method=="splines")
        control$legend$title <- "M-spline intensity model"
      else
        control$legend$title <- "Weibull model"
    }
    save.xpd <- par()$xpd
    par(xpd=TRUE)
    do.call("legend",control$legend)
    par(xpd=save.xpd)
  }
  # }}}
  invisible(x)
}
# }}}
