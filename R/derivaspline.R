derivaspline<-function(b,np,npar,bfix,fix,zi01,zi02,zi12,ctime,no,nz01,nz02,nz12,ve01,ve02,ve12,
                       dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                       t0,t1,t2,t3,troncature){
  res<-rep(0,(np*(np+1)/2)+np)
  .Fortran("derivaspline",
           ## input
           as.double(b),
           as.integer(np),
           as.integer(npar),
           as.double(bfix),
           as.integer(fix),
           as.double(zi01),
           as.double(zi12),
           as.double(zi02),
           as.integer(ctime),
           as.integer(no),
           as.integer(nz01),
           as.integer(nz12),
           as.integer(nz02),
           as.double(ve01),
           as.double(ve12),
           as.double(ve02),
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


