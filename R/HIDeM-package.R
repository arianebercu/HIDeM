#' HIDeM: High-dimensional Illness-death Models
#'
#' This package provides tools for fitting regularized illness-death models 
#' with interval-censored data.
#' 
#' \strong{Main Functions:}
#' \itemize{
#'   \item{\code{idm()}: Fit the illness-death model with variable selection.}
#'   \item{\code{gridsearch()}: Fit transition specific regularized illness-death model to obtain grid values of penalty parameters in order to run complete regularized illness-death model -see paper for detailed approach.}
#'   \item{\code{boostrap()}: Perform boostrap on individual to evaluate variable selection stability of regularized illness-death model.}
#'   \item{\code{simulateIDM()}: Simulate data from an illness-death model with interval censored event times and covariates.}
#'   \item{\code{intensity()}: Estimate transition intensities}
#'   \item{\code{predict()}: Predict life-time expectancy, survival and transition intensities.}
#'   \item{\code{plot}: Plot estimated baseline transition intensities.}
#'   \item{\code{summary()}: Summarize model results.}
#' }
#' @name HIDeM
#' @keywords package
"_PACKAGE"