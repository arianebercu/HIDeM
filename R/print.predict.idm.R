
### Code:
#' @param object an \code{idm} class objects returned by a call to the
#'     \code{\link{idm}} function
#' @useDynLib SmoothHazardoptim9
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
#' @export

print.predict.idm <- function(x,digits=3,...){
    lifeExpect <- (is.infinite(x$t))|("LE.00"%in%x$transprob$Parameter)
    cat("Predictions of an irreversible illness-death model with states (0,1,2).\n\n")
    cat("For covariate values:\n\n")
    print(x$newdata,row.names=FALSE)
    cat("\n")
    fmt <- paste0("%1.", digits[[1]], "f")
    px <- x$transprob
    for (j in 2:NCOL(px))
        px[,j] <- sprintf(px[,j],fmt=fmt)
    rownames(px) <- NULL
    if (lifeExpect==TRUE){
        cat("Remaining life expected sojourn times (starting at time ",x$s,"):\n\n",sep="")
        print(cbind("State at time s"=c("0","0","0","1"),"Expected years in states 0,1"=c("In state 0","Total","In state 1","Total"),px[px$Parameter %in% c("LE.00","LE.0.","LE.01","LE.11"),]),row.names=FALSE)
        cat("Life time risk at time ",x$s,":\n",sep="")
        print(cbind(px[px$Parameter %in% c("LTR"),]),row.names=FALSE)
    }else{
        cat("For a subject in state '0' at time ",x$s,",\npredicted state occupation probability at time ",x$t,":\n\n",sep="")
        print(cbind("State"=c(0,1,2,0),px[px$Parameter %in% c("p00","p01","p02","RM"),]),row.names=FALSE)
        cat("\nThe probability p02 can be further decomposed into\ndirect and indirect transition probabilities:\n\n")
        print(cbind("Path"=c("direct","via state 1","total"),px[px$Parameter %in% c("p02_0","p02_1","p02"),]),row.names=FALSE)
        cat("\nFor a subject in state '0' at time ",x$s,",\npredicted probability of exit from state 0 until time ",x$t,":\n\n",sep="")
        print(cbind("Path"=c("via state 1","any"),px[px$Parameter %in% c("F01","F0."),]),row.names=FALSE)
        cat("\nFor a subject in state '1' at time ",x$s,",\npredicted state occupation probability at time ",x$t,":\n\n",sep="")
        print(cbind("State"=c(1,2),px[px$Parameter %in% c("p11","p12"),]),row.names=FALSE)
    }
    invisible(px)
}



#----------------------------------------------------------------------
### print.predict.idm.R ends here
