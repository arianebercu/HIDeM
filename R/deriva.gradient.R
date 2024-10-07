#' Numerical derivatives
#'
#' The function computes the first derivates and the information score matrix.
#' Central finite-differences and forward finite-differences are used for the first
#' and second derivatives respectively.
#'
#' @param nproc number of processors for parallel computing
#' @param b value of parameters to be optimized over
#' @param funcpa function to be minimized (or maximized), with argument the vector
#' of parameters over which minimization isto take place.
#' It should return a scalar result.
#' @param .packages character vector of packages that funcpa depends on
#' @param \dots other arguments of the funcpa function
#'
#' @return \item{v}{vector containing the upper part of the information score
#' matrix and the first derivatives} \item{rl}{the value of the funcpa function
#' at point b}
#' @author Viviane Philipps, Boris Hejblum, Cecile Proust-Lima, Daniel Commenges
#' @references Donald W. Marquardt An algorithm for Least-Squares Estimation of Nonlinear Parameters. Journal of the Society for Industrial and Applied Mathematics, Vol. 11, No. 2. (Jun, 1963), pp. 431-441.
#'
#' @examples
#' b <- 0.1
#' f <- function(b){return((2*b[1]**2+3*b[1]))}
#' d <- deriva(b=b,funcpa=f)
#'
#'@useDynLib SmoothHazardoptim9
#' @export
#'


deriva.gradient <- function(nproc=1,b,funcpa,.packages=NULL,...){
  
  m <- length(b)
  bh2 <- bh <- rep(0,m)
  v <- rep(0,(m*(m+3)/2))
  fcith <- fcith2 <- rep(0,m)
  ## function 
  rl <- funcpa(b,...)
  
  if(nproc>1)
  {
    ### remplacer les 2 boucles par une seule
    grid <- cbind(c(rep(1:m,1:m),1:m),c(unlist(sapply(1:m,function(k) seq(1,k))),rep(0,m)))
    mm <- nrow(grid)
    h <- sapply(b,function(x){max(1E-7,(1E-4*abs(x)))})
    
    
    k<-NULL # for cran check 
    ## derivees premieres:
    ll <- foreach::foreach(k=(m*(m+1)/2)+1:m,
                  .combine=cbind,
                  .packages=.packages) %dopar%
      {
        i <- grid[k,1]
        
        bp <- b
        bp[i] <- b[i]+h[i]
        av <- funcpa(bp,...)
        
        bm <- b
        bm[i] <- b[i]-h[i]
        ar <- funcpa(bm,...)
        
        d <- (av-ar)/(2*h[i])
        
        c(av,d)
      }
    v <- ll[2,] 

  }
  else
  {
    ## gradient null
    v<-rep(NA,m)
    for(i in 1:m){
      bh <- bh2 <- b
      th <- max(1E-7,(1E-4*abs(b[i])))
      bh[i] <- bh[i] + th
      bh2[i] <- bh2[i] - th
      
      fcith[i] <- funcpa(bh,...)
      fcith2[i] <- funcpa(bh2,...)
      thn<-max(1E-7,(1E-4*abs(b[i])))
      v[i]<--(fcith[i]-fcith2[i])/(2*thn)
    }
  }
    
   
  result <- list(v=v,rl=rl)
  return(result)
}
