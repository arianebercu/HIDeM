% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gridsearch.R
\name{gridsearch}
\alias{gridsearch}
\title{Grid-search for lambda penalty parameters in illness-death model}
\usage{
gridsearch(
  gridmethod = "BIC",
  sizegrid = c(3, 3, 3),
  formula01,
  formula02,
  formula12,
  data,
  method = "Weib",
  scale.X = T,
  maxiter = 100,
  maxiter.pena = 10,
  eps = c(5, 5, 3),
  n.knots = NULL,
  knots = "equidistant",
  type.quantile = 1,
  na.action = na.fail,
  B = NULL,
  posfix = NULL,
  gausspoint = 10,
  lambda01 = NULL,
  lambda02 = NULL,
  lambda12 = NULL,
  nlambda01 = 50,
  nlambda02 = 50,
  nlambda12 = 50,
  penalty = "lasso",
  penalty.factor = NULL,
  alpha = ifelse(penalty == "scad", 3.7, ifelse(penalty == "mcp", 3, ifelse(penalty
    \%in\% c("elasticnet", "corrected.elasticnet"), 0.5, 1))),
  nproc = 1,
  clustertype = "FORK"
)
}
\arguments{
\item{gridmethod}{On which indicator the grid should be based, either BIC or GCV}

\item{sizegrid}{The number of values that we want for each lambda, its a vector 
of three elements, \code{0 --> 1}, \code{0 --> 2} and \code{1 --> 2}. The size of each
element should not exceed respectively nlambda01, nlambda02 and nlambda12}

\item{formula01}{A formula specifying a regression model for the
\code{0 --> 1} transition from the initial state to the transient
state of the illness-death model.  The right hand side of the
formula specifies the covariate terms, and the left hand side must
be an event history object as returned by the function \code{Hist}.}

\item{formula02}{A formula specifying a regression model for the
\code{0 --> 2} transition from the initial state to the absorbing
state. The left hand side must be equal to the left hand side of
\code{formula01}. If missing it is set to \code{formula01}.}

\item{formula12}{A formula specifying a regression model for the
\code{1 --> 2} transition from the transient state to the absorbing
state.  operator is not required. If missing it is set to
\code{formula01}.}

\item{data}{A data frame in which to interpret the variables of
\code{formula01}, \code{formula02} and \code{formula12}.}

\item{method}{type of estimation method: "splines" for a likelihood approach with approximation of the transition
intensities by M-splines, "Weib" for a parametric approach with a
Weibull distribution on the transition intensities. Default is
"Weib".}

\item{scale.X}{do you want to center and reduce your explanatory variables}

\item{maxiter}{Maximum number of iterations. The default is 200.}

\item{maxiter.pena}{Maximum number of iterations for penalised coefficients}

\item{eps}{A vector of 3 integers >0 used to define the power of
three convergence criteria: 1. for the regression parameters,
2. for the likelihood, 3. for the second derivatives. The default
is \code{c(5,5,3)} which is translated into convergence if the
respective values change less then \eqn{10^{-5}} (for regression
parameters and likelihood) and \eqn{10^{-3}} for the second
derivatives between two iterations.}

\item{n.knots}{For \code{method="splines"} only, a vector of length
3 specifing the number of knots, one for each transition, for the
M-splines estimate of the baseline intensities in the order \code{0
--> 1}, \code{0 --> 2}, \code{1 --> 2}. The default is c(7,7,7). When \code{knots}
are specified as a list this argument is ignored.
The algorithm needs least 5 knots and at most 20 knots.}

\item{knots}{Argument only active for the penalized likelihood approach \code{method="Splines"}.
There are three ways to control the placement of the knots between the smallest and the largest
of all time points:
\describe{
 \item{\code{knots="equidistant"}}{Knots are placed with same distance on the time scale.}
 \item{\code{knots="quantile"}}{Knots are placed such that the number of observations is roughly the same between knots.}
\item{knots=list()}{List of 1 or 2 or three vectors. The list elements are the actual placements
(timepoints) of the knots for the M-spline. The list may contain
one vector of placements for each transition in the order \code{0 --> 1}, \code{0 --> 2}, \code{1 --> 2}.
If only vector is specifified the knots are used for all transitions. If only 2 vectors are specifified, the
knots for the \code{0 --> 1} transition are also used for the \code{1 --> 2} transition.}
}
The algorithm needs at least 3 knots in spline and allows no more than 20 knots.}

\item{type.quantile}{Argument only active for the likelihood approach \code{method="splines"}.
There are three ways to control the placement of the knots  according to the time considered between states :=
\describe{
 \item{\code{type.quantile=1}}{Time for \code{0 --> 1} is the imputed to the middle of the interval left and right for demence . Time for \code{0 --> 2}
 and \code{1 --> 2} is the same t, time of news. }
 \item{\code{type.quantile=2}}{Time for \code{0 --> 1} is the imputed to the middle of the interval left and right. Time for \code{0 --> 2}
 and \code{1 --> 2} is the same t, time of news. }
\item{\code{type.quantile=3}}{Time for \code{0 --> 1} is the imputed to the middle of the interval left and right. Time for \code{0 --> 2}
 is time of death for non demented sujects only. Time for \code{1 --> 2} is time of death for suject diagnose with dementia. }
 \item{\code{type.quantile=4}}{Time for \code{0 --> 1} is left and right. Time for \code{0 --> 2}
 is time of death for non demented sujects only. Time for \code{1 --> 2} is time of death for suject diagnose with dementia. }
}}

\item{na.action}{how NAs are treated. The default is first, any
na.action attribute of data, second a na.action setting of options,
and third 'na.fail' if that is unset. The 'factory-fresh' default
is na.omit. Another possible value is NULL.}

\item{B}{vector of size the number of parameters, in the following order, first the parameters of splines \code{0 --> 1}, \code{0 --> 2}, \code{1 --> 2},
second the parameters of explanatory variables in order  \code{0 --> 1}, \code{0 --> 2}, \code{1 --> 2}.
This argument is only used for models with M-splines.}

\item{posfix}{index of fixed parameters}

\item{gausspoint}{gauss quadrature points in the approximation of integrals}

\item{lambda01}{Lambda on transition 0 --> 1}

\item{lambda02}{Lambda on transition 0 --> 2}

\item{lambda12}{Lambda on transition 1 --> 2}

\item{nlambda01}{number of Lambda on transition 0 --> 1}

\item{nlambda02}{number of Lambda on transition 0 --> 2}

\item{nlambda12}{number of Lambda on transition 1 --> 2}

\item{penalty}{which penalty to consider}

\item{penalty.factor}{which variable should be penalised}

\item{alpha}{alpha on all transitions}

\item{nproc}{number of cluster}

\item{clustertype}{in which cluster to work}
}
\value{
\item{m01}{Model estimated on 0 --> 1} \item{m02}{Model estimated on 0 --> 2} 
\item{m12}{Model estimated on 1 --> 2} \item{lambda01}{vector of lambda penalty parameters
on transition 0 --> 1 minimising the BIC or GCV in model m01} \item{lambda02}{vector of lambda penalty parameters
on transition 0 --> 2 minimising the BIC or GCV in model m02} \item{lambda12}{vector of lambda penalty parameters
on transition 1 --> 2 minimising the BIC or GCV in model m12}
\item{alpha}{the penalty parameter alpha} 
\item{sizegrid}{the size of lambda penalty parameters for each transition 0 -->1, 0 -->2 and 1 -->2}
\item{gridmethod}{On which indicator the grid should be based, either BIC or GCV}
}
\description{
Grid-search for lambda penalty parameters in illness-death model
}
\examples{
{
\dontrun{
library(lava)
library(prodlim)
set.seed(17)
d <- simulateIDM(n=1000,beta01=c(1,1,0,0.5,0.5,rep(0,5)),
beta02=c(1,0,0,0,0.5,rep(0,5)),beta12=c(1,0,0,0,0.5,rep(0,5)))$data
fitgrid<- gridsearch(formula01=Hist(time=list(L,R),
event=seen.ill)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
formula02=Hist(time=observed.lifetime,
event=seen.exit)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
formula12=~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,
 nlambda01=20, nlambda02=20, nlambda12=20,
data=d,penalty="lasso")
}}
}
\references{
D. Marquardt (1963). An algorithm for least-squares estimation
of nonlinear parameters.  \emph{SIAM Journal of Applied Mathematics},
431-441.
}
\seealso{
\code{\link{print.idm}}\code{\link{summary.idm}} \code{\link{predict.idm}}
}
\author{
R: Ariane Bercu <ariane.bercu@u-bordeaux.fr>
}
