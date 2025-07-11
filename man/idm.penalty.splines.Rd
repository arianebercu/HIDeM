% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/idm.penalty.splines.R
\name{idm.penalty.splines}
\alias{idm.penalty.splines}
\title{Penalised Illness-death model algorithm with M-splines baseline risk}
\usage{
idm.penalty.splines(
  b,
  fix0,
  size_V,
  size_spline,
  clustertype,
  epsa,
  epsb,
  epsd,
  eps.eigen,
  nproc,
  maxiter,
  maxiter.pena,
  knots01,
  knots02,
  knots12,
  ctime,
  N,
  nknots01,
  nknots02,
  nknots12,
  ve01,
  ve02,
  ve12,
  dimnva01,
  dimnva02,
  dimnva12,
  nvat01,
  nvat02,
  nvat12,
  t0,
  t1,
  t2,
  t3,
  troncature,
  gausspoint,
  nlambda01,
  lambda01,
  nlambda02,
  lambda02,
  nlambda12,
  lambda12,
  alpha,
  penalty.factor,
  penalty,
  analytics,
  partialH
)
}
\arguments{
\item{b}{parameters not fixed}

\item{fix0}{indicators of fixed and unfixed parameters}

\item{size_V}{number of parameters}

\item{size_spline}{number of splines parameters}

\item{clustertype}{in which cluster to work}

\item{epsa}{control convergence parameter for beta}

\item{epsb}{control convergence parameter for loglik}

\item{epsd}{control convergence for distance to minimum rdm}

\item{eps.eigen}{the power of convergence for eigen values of covariance matrix only}

\item{nproc}{number of cluster}

\item{maxiter}{Maximum number of iterations. The default is 200.}

\item{maxiter.pena}{Maximum number of iterations for penalised coefficients}

\item{knots01}{knots of transition 0 --> 1}

\item{knots02}{knots of transition 0 --> 2}

\item{knots12}{knots of transition 1 --> 2}

\item{ctime}{classification of subject according to their observations}

\item{N}{number of subjects}

\item{nknots01}{number of knots for transition 0 -->1}

\item{nknots02}{number of knots for transition 0 -->2}

\item{nknots12}{number of knots for transition 1 -->2}

\item{ve01}{variables for transition 0 -->1}

\item{ve02}{variables for transition 0 -->2}

\item{ve12}{variables for transition 1 -->2}

\item{dimnva01}{number of variables for transition 0 -->1}

\item{dimnva02}{number of variables for transition 0 -->2}

\item{dimnva12}{number of variables for transition 1 -->2}

\item{nvat01}{number of variables for transition 0 -->1}

\item{nvat02}{number of variables for transition 0 -->2}

\item{nvat12}{number of variables for transition 1 -->2}

\item{t0}{time entry}

\item{t1}{time L}

\item{t2}{time R}

\item{t3}{time of event/out}

\item{troncature}{indicator if troncature or not}

\item{gausspoint}{number of points in gauss quadrature}

\item{nlambda01}{number of Lambda on transition 0 --> 1}

\item{lambda01}{Lambda on transition 0 --> 1}

\item{nlambda02}{number of Lambda on transition 0 --> 2}

\item{lambda02}{Lambda on transition 0 --> 2}

\item{nlambda12}{number of Lambda on transition 1 --> 2}

\item{lambda12}{Lambda on transition 1 --> 2}

\item{alpha}{alpha on all transitions}

\item{penalty.factor}{which variable should be penalised}

\item{penalty}{which penalty to consider}

\item{partialH}{which hessian is computed for Newton-Raphson path of penalised parameters}

\item{methodCV}{methodolgy to optimise penalised parameters}
}
\description{
Penalised Illness-death model algorithm with M-splines baseline risk
}
\author{
R: Ariane Bercu <ariane.bercu@u-bordeaux.fr>
}
