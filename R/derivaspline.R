derivaspline<-function(b0,np0,npar0,bfix0,fix0,zi010,zi020,zi120,c0,no0,nz010,nz020,nz120,ve010,ve020,ve120,
                       dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                       t0,t1,t2,t3,troncature){
  res<-rep(0,(np0*(np0+1)/2)+np0)
  .Fortran("derivaspline",
           ## input
           as.double(b0),
           as.integer(np0),
           as.integer(npar0),
           as.double(bfix0),
           as.integer(fix0),
           as.double(zi010),
           as.double(zi120),
           as.double(zi020),
           as.integer(c0),
           as.integer(no0),
           as.integer(nz010),
           as.integer(nz120),
           as.integer(nz020),
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


