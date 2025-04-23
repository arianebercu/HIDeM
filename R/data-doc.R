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
#' @source Simulations internes
"mydata"
