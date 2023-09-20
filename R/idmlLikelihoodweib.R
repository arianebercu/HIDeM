#' @useDynLib SmoothHazardoptim8
#' @export


idmlLikelihoodweib<-function(b0,np0,npar0,bfix0,fix0,c0,no0,ve010,ve020,ve120,
                         dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                         t00,t10,t20,t30,troncature0){
  res<-0

  .Fortran("idmlikelihoodweib",
           ## input
           as.double(b0),
           as.integer(np0),
           as.integer(npar0),
           as.double(bfix0),
           as.integer(fix0),
           as.integer(c0),
           as.integer(no0),
           as.double(ve010),
           as.double(ve120),
           as.double(ve020),
           as.integer(dimnva01),
           as.integer(dimnva12),
           as.integer(dimnva02),
           as.integer(nva01),
           as.integer(nva12),
           as.integer(nva02),
           as.double(t00),
           as.double(t10),
           as.double(t20),
           as.double(t30),
           as.integer(troncature0),
           likelihood_res=as.double(res),
           PACKAGE="SmoothHazardoptim9")$likelihood_res
}


