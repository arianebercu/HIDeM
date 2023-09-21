

#' Data set for survival models: right-censored and interval-censored data.
#' 
#' A simulated data frame for survival models composed of right-censored and
#' interval-censored data.
#' 
#' 
#' @name testdata
#' @docType data
#' @format A data frame with 936 observations on the following 4 variables.
#' \describe{ \item{l}{for diseased subjects: left endpoint of
#' censoring interval; for non-diseased subjects: right censoring time}
#' \item{r}{for diseased subjects: right endpoint of censoring
#' interval; for non-diseased subjects: right censoring time for the disease
#' event} \item{id}{disease status} \item{cov}{covariate} }
#' @keywords datasets
#' @examples
#' 
#' data(testdata)
#' head(testdata)
#' 
NULL