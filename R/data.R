#' Paq1000 dataset
#'
#' Jeu de données simulé pour HIDeM
#'
#' Un jeu de données contenant 2000 observations simulées avec 50 covariables, utilisé pour ajuster un modèle illness-death avec des intensité de transition de base suivant une loi Weibull.
#'
#' @format Un data frame avec 1000 lignes et 53 colonnes :
#' \describe{
#'   \item{observed.lifetime}{Time of death if seen.exit=1 otherwise last time with information on the vital status}
#'   \item{L}{Time of last visit illness-free if seen.ill=1 otherwise time of last visit}
#'   \item{R}{Diagnosis visit time if seen.ill=1 otherwise time of last visit}
#'   \item{seen.ill}{illness indicator, 0=censored, 1=ill}
#'   \item{seen.exit}{death indicator, 0=censored, 1=dead}
#'   \item{X1 à X50}{Covariables simulées}
#' }
#' @docType data
#' @source Simulations internes
#' @usage data(Paq1000)
#' @keywords datasets
"Paq1000"

#' m01example dataset
#'
#' Description of the m01example dataset.
#' @docType data
#' @format model on transition specific 0 to 1 using data set simulated
#' @usage data(m01example)
#' @keywords datasets
"m01example"

#' m02example dataset
#'
#' Description of the m02example dataset.
#' @docType data
#' @format model on transition specific 0 to 2 using data set simulated
#' @usage data(m02example)
#' @keywords datasets
"m02example"

#' m12example dataset
#'
#' Description of the m12example dataset.
#' @docType data
#' @format model on transition specific 1 to 2 using data set simulated
#' @usage data(m12example)
#' @keywords datasets
"m12example"

#' m01Paq1000 dataset
#'
#' Description of the m01Paq1000 dataset.
#' @docType data
#' @format model on transition specific 0 to 1 using Paquid data set
#' @usage data(m01Paq1000)
#' @keywords datasets
"m01Paq1000"

#' m02Paq1000 dataset
#'
#' Description of the m02Paq1000 dataset.
#' @docType data
#' @format model on transition specific 0 to 2 using Paquid data set
#' @usage data(m02Paq1000)
#' @keywords datasets
"m02Paq1000"

#' m12Paq1000 dataset
#'
#' Description of the m12Paq1000 dataset.
#' @docType data
#' @format model on transition specific 1 to 2 using Paquid data set
#' @usage data(m12Paq1000)
#' @keywords datasets
"m12Paq1000"

#' bootPaq1000 dataset
#'
#' Description of the bootPaq1000 dataset.
#' @docType data
#' @format boostrap model on Paquid data set
#' @usage data(bootPaq1000)
#' @keywords datasets
"bootPaq1000"

#' modelexample dataset
#'
#' Description of the modelexample dataset.
#' @docType data
#' @format Final model on simulated data set
#' @usage data(modelexample)
#' @keywords datasets
"modelexample"

#' modelPaq1000 dataset
#'
#' Description of the modelPaq1000 dataset.
#' @docType data
#' @format Final model on Paquid data set
#' @usage data(modelPaq1000)
#' @keywords datasets
"modelPaq1000"

#' mvoidexample dataset
#'
#' Description of the mvoidexample dataset.
#' @docType data
#' @format void model on simulated data set
#' @usage data(mvoidexample)
#' @keywords datasets
"mvoidexample"

#' mvoidPaq1000 dataset
#'
#' Description of the mvoidPaq1000 dataset.
#' @docType data
#' @format void model on Paquid data set
#' @usage data(mvoidPaq1000)
#' @keywords datasets
"mvoidPaq1000"

#' remodelexample dataset
#'
#' Description of the remodelexample dataset.
#' @docType data
#' @format Re estimated final model on simulated data set
#' @usage data(remodelexample)
#' @keywords datasets
"remodelexample"

#' remodelPaq1000 dataset
#'
#' Description of the remodelPaq1000 dataset.
#' @docType data
#' @format Re estimated final model on Paquid data set
#' @usage data(remodelPaq1000)
#' @keywords datasets
"remodelPaq1000"
