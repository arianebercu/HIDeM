### Code:
##' @title Define the time points at which Y need to be estimated
##' @param lower.intdouble lower bound of integral for double integral estimation
##' @param upper.intdouble upper bound of integral for double integral estimation
##' @param end.time upper bound of integral from 0 to end of follow-up (death or censoring)
##' @param truncation TRUE if left-truncation, otherwise FALSE
##' @param entry.time time of entry in the study (values 0 if no left truncation, otherwise >0)
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  

gauss_kronrod_points <- function(lower.intdouble, 
                                 upper.intdouble,
                                 end.time,
                                 truncated,
                                 entry.time) {
  
  # Gauss-Kronrod 15-point nodes and weights for interval [-1, 1]
  gk15_nodes <- c(
    0.0000000000000000,
    -0.2077849550078985,  0.2077849550078985,
    -0.4058451513773972,  0.4058451513773972,
    -0.5860872354676911,  0.5860872354676911,
    -0.7415311855993945,  0.7415311855993945,
    -0.8648644233597691,  0.8648644233597691,
    -0.9491079123427585,  0.9491079123427585,
    -0.9914553711208126,  0.9914553711208126
  )
  # Gauss-Kronrod 15-point nodes and weights for interval [-1, 1]
  # and external point 1
  
  gk15_nodes_ext <- c(
    0.0000000000000000,
    -0.2077849550078985,  0.2077849550078985,
    -0.4058451513773972,  0.4058451513773972,
    -0.5860872354676911,  0.5860872354676911,
    -0.7415311855993945,  0.7415311855993945,
    -0.8648644233597691,  0.8648644233597691,
    -0.9491079123427585,  0.9491079123427585,
    -0.9914553711208126,  0.9914553711208126,
    1
  )
  # Transform from [-1,1] to [a,b]
  x<-0.5 * ((upper.intdouble - lower.intdouble) * gk15_nodes + (upper.intdouble + lower.intdouble)) # timepoint necessary to estimate the outer integral of int(a,b)int(0,t)f(u)dug(t)dt
  x<-0.5 * (matrix(x,ncol=1)%*% gk15_nodes_ext + x) # timepoint necessary to estimate the inner integral of int(a,b)int(0,t)f(u)dug(t)dt + needed added one as we have g(t) 
  
  # add integral points from 0 to end of study and need end.time value also (death or censoring): 
  # attention to take row by row t(x)
  x<-c(t(x),0.5 * (end.time * gk15_nodes_ext + end.time)) 
  if(truncated==T){
    x<-c(x,0.5 * (entry.time * gk15_nodes + entry.time) )
  }
  return(as.vector(x))
}

  

  


