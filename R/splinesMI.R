#' @param x vector in which splines will be evaluated
#' @param knots internal nodes values
#' @param Boundary.knots external nodes values
#' @keywords methods
#' @examples
#' @useDynLib SmoothHazardoptim9
#' @export
splinesMI<-function(x,knots,Boundary.knots){
  if(any(Boundary.knots%in%knots)){stop("All interval nodes must be different from external nodes")}
  knots<-c(rep(knots.bound[1],4),knots,rep(knots.bound[2],4))
  number.knots<-length(knots)-2
  Mspline<-matrix(0,nrow=length(x),ncol=number.knots-2)
  Ispline<-matrix(0,nrow=length(x),ncol=number.knots-2)
  k<-1
  for (i in 5:(number.knots+3)) {
    TF = ( (knots[i-1]<=x) & (x<knots[i]) )
    if (sum(TF) != 0) { 
      ind = which(TF) 
      mm3=rep(0,length(ind))
      mm2=rep(0,length(ind))
      mm1=rep(0,length(ind))
      mm=rep(0,length(ind))
      im3=rep(0,length(ind))
      im2=rep(0,length(ind))
      im1=rep(0,length(ind))
      im=rep(0,length(ind))
      j = i-1
      if (j>4) { 
        som = sum(theta[1:(j-4)])
      }
      ht = x[ind]-knots[j] #
      htm = x[ind]-knots[j-1] #
      h2t = x[ind]-knots[j+2] #
      ht2 = knots[j+1]-x[ind] #
      ht3 = knots[j+3]-x[ind] #
      hht = x[ind]-knots[j-2] #
      h = knots[j+1]-knots[j]
      hh = knots[j+1]-knots[j-1]
      h2 = knots[j+2]-knots[j]
      h3 = knots[j+3]-knots[j]
      h4 = knots[j+4]-knots[j]
      h3m = knots[j+3]-knots[j-1]
      h2n = knots[j+2]-knots[j-1]
      hn= knots[j+1]-knots[j-2]
      hh3 = knots[j+1]-knots[j-3]
      hh2 = knots[j+2]-knots[j-2]
      mm3[ind] = ((4*ht2*ht2*ht2)/(h*hh*hn*hh3))
      mm2[ind] = ((4*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4*h2t*htm*ht2)/(hh2*h2n*hh*h))+((4*h2t*h2t*ht)/(hh2*h2*h*h2n))
      mm1[ind] = (4*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4*htm*ht*h2t)/(h3m*h2*h*h2n))+((4*ht3*ht*ht)/(h3m*h3*h2*h))
      mm[ind] = 4*(ht*ht*ht)/(h4*h3*h2*h)
      im3[ind] = (0.25*(x[ind]-knots[j-3])*mm3[ind])+(0.25*hh2*mm2[ind])+(0.25*h3m*mm1[ind])+(0.25*h4*mm[ind])
      im2[ind] = (0.25*hht*mm2[ind])+(h3m*mm1[ind]*0.25)+(h4*mm[ind]*0.25)
      im1[ind] = (htm*mm1[ind]*0.25)+(h4*mm[ind]*0.25)
      im[ind] = ht*mm[ind]*0.25
      
    } # fin if (sum(TF) != 0)
  } # fin for
  Mspline[,1]<-mm3
  Mspline[,2]<-mm2
  Mspline[,3]<-mm1
  Mspline[,4]<-mm
  Ispline[,1]<-im3
  Ispline[,2]<-im2
  Ispline[,3]<-im1
  Ispline[,4]<-im
  return(list(Mspline=Mspline,Ispline=Ispline))
}
