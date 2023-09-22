derivaspline<-function(b,np0,npar0,bfix0,fix0,zi01,zi02,zi12,c0,no0,nz01,nz02,nz12,ve010,ve020,ve120,
                       dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                       t0,t1,t2,t3,troncature){
  res<-rep(0,(np0*(np0+1)/2)+np0)
  .Fortran("derivaspline",
           ## input
           as.double(b),
           as.integer(np0),
           as.integer(npar0),
           as.double(bfix0),
           as.integer(fix0),
           as.double(zi01),
           as.double(zi12),
           as.double(zi02),
           as.integer(c0),
           as.integer(no0),
           as.integer(nz01),
           as.integer(nz12),
           as.integer(nz02),
           as.double(ve010),
           as.double(ve120),
           as.double(ve020),
           as.integer(dimnva01),
           as.integer(dimnva12),
           as.integer(dimnva02),
           as.integer(nva01),
           as.integer(nva12),
           as.integer(nva02),
           as.double(t0),
           as.double(t1),
           as.double(t2),
           as.double(t3),
           as.integer(troncature),
           likelihood_deriv=as.double(res),
           PACKAGE="SmoothHazardoptim9")$likelihood_deriv
}


