% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.idm.R
\name{plot.idm}
\alias{plot.idm}
\title{Plot method for an illness-death model}
\usage{
\method{plot}{idm}(
  x,
  lambda = NULL,
  conf.int = FALSE,
  citype = "shadow",
  add = FALSE,
  axes = TRUE,
  col,
  lwd,
  lty,
  xlim,
  ylim,
  xlab,
  ylab,
  cex,
  y.intersp,
  legend = TRUE,
  transition = c("01", "02", "12"),
  ...
)
}
\arguments{
\item{x}{a \code{idmWeib} class object (output from calling
\code{idm} with the (default) option \code{intensities}="Weib".}

\item{lambda}{give the lambda for which you want to plot the penalised illness-death model}

\item{conf.int}{If TRUE show confidence limits}

\item{citype}{Type of confidence limits, can be "shadow" or "bars"}

\item{add}{If TRUE add to existing plot}

\item{axes}{If TRUE axes are drawn}

\item{col}{Color of the lines}

\item{lwd}{Width of the lines}

\item{lty}{Type of the lines}

\item{xlim}{Limits for x-axis}

\item{ylim}{Limits for y-axis}

\item{xlab}{Label for x-axis}

\item{ylab}{Label for y-axis}

\item{cex}{character expansion factor relative to current par("cex"). Used for text, and provides the default for pt.cex.}

\item{y.intersp}{vertical (y) distances (in lines of text shared above/below each legend entry). A vector with one element for each row of the legend can be used.}

\item{legend}{If TRUE a legend is drawn, which can be further controlled via \code{\link{SmartControl}}.}

\item{transition}{Choose one of the transition intensities: \code{c("01","02","12")}.}

\item{...}{Passed to \code{\link{SmartControl}}}
}
\value{
Print a plot of the baseline transition intensities of an
illness-death model estimated using a Weibull approach.
}
\description{
Plot estimated baseline transition intensities from an object of class
\code{idm} optionally with confidence limits.
}
\examples{

\dontrun{
library(lava)
library(prodlim)
set.seed(17)
d <- simulateIDM(n=1000,
                 beta01=c(1,1,0,0.5,0.5,rep(0,5)),
                 beta02=c(1,0,0,0,0.5,rep(0,5)),
                 beta12=c(1,0,0,0,0.5,rep(0,5)))$data
                 
fitsplines <- idm(formula01=Hist(time=list(L,R),event=seen.ill)~X1+X2,
formula02=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
formula12=Hist(time=observed.lifetime,event=seen.exit)~X1+X2,
data=d,method="splines")
plot(fitsplines)

fitpenspline <- idm(formula01=Hist(time=list(L,R),
event=seen.ill)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                    formula02=Hist(time=observed.lifetime,
                    event=seen.exit)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                    formula12=~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
                    method="splines",
                    data=d,penalty="lasso",
                    lambda01 = c(10,20),
                    lambda02 = 10, lambda12 = 10)
plot(fitpenspline,lambda=c(10,10,10))
}
 
}
\seealso{
\code{\link{print.idm}},\code{\link{summary.idm}},\code{\link{idm}},

\code{\link{idm}}
}
\author{
R: Ariane Bercu <ariane.bercu@u-bordeaux.fr>
}
\keyword{methods}
