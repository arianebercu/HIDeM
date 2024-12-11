### Code:
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
#' @references Donald W. Marquardt An algorithm for Least-Squares Estimation of Nonlinear Parameters. Journal of the Society for Industrial and Applied Mathematics, Vol. 11, No. 2. (Jun, 1963), pp. 431-441.
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
    #h <- sapply(b,function(x){.Machine$double.eps^(1/4)})
  
    
    k<-NULL # for cran check 
    ## derivees premieres:
    ll <- foreach::foreach(k=(m*(m+1)/2)+1:m,
                  .combine=cbind,
                  .packages=.packages) %dopar%
      {
        i <- grid[k,1]
        
        bp <- b
        #bp[i] <- sqrt((b[i] + h[i])^2 -2*b[i]*h[i]-h[i])
        bp[i] <- b[i] + h[i]
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
    ## gradient null attention : as put to square in definition
    ## we want for b^2+th not (b+th)^2
    v<-rep(NA,m)
    for(i in 1:m){
      bh <- bh2 <- b
      th <- max(1E-7,(1E-4*abs(b[i])))
      #th<-.Machine$double.eps^(1/4)
      #bh[i] <-sqrt((bh[i] + th)^2 -2*bh[i]*th-th)
      bh[i] <- bh[i] + th
      #bh2[i] <-sqrt((bh2[i] -th)^2 +2*bh2[i]*th+th)
      bh2[i] <-bh2[i]-th
      
      fcith[i] <- funcpa(bh,...)
      fcith2[i] <- funcpa(bh2,...)
      thn<-max(1E-7,(1E-4*abs(b[i])))
      #thn<-.Machine$double.eps^(1/4)
      v[i]<--(fcith[i]-fcith2[i])/(2*thn)
    }
  }
    
   
  result <- list(v=v,rl=rl)
  return(result)
}


deriva.hessianweib <- function(nproc=1,b,fix,funcpa,.packages=NULL,...){
     nvarmodelPar<-sum(fix[c(1:6)]==0)
      m <- length(b)
      bh <- rep(0,m)
      fcith <- rep(0,m)
      ## function 
      rl <- funcpa(b,fix=fix,...)
      ## gradient null
      for(i in 1:m){
        bh <- b
        th <- max(1E-7,(1E-4*abs(b[i])))
        #th <- .Machine$double.eps^(1/4)
        #if(i<=nvarmodelPar){
          #bh[i] <- sqrt((bh[i] + th)^2 -2*bh[i]*th-th)
          #fcith[i] <- funcpa(bh,fix=fix,...)
        #}else{
        bh[i] <- bh[i] + th
        fcith[i] <- funcpa(bh,fix=fix,...)
        #}
      }
      
      # mm<-0
      # for(k in 1:nvarmodelPar){
      #   mm<-mm+(m-k+1)
      # }
      # VmodelPar<-rep(NA,mm)
      # k<-1
      
      VmodelPar<-matrix(0,m,m)
      k<-1
      for(i in 1:m){
        bhi <- b
        thi <- max(1E-7,(1E-4*abs(b[i])))
        bhi[i] <-bhi[i]+thi
        maxj<-i
        if(i>nvarmodelPar){maxj<-nvarmodelPar}
        for(j in 1:maxj){
          
          bhj <- bhi
          thj <- max(1E-7,(1E-4*abs(b[j])))
          #thi<-thj<-.Machine$double.eps^(1/4)
          th <- thi * thj
          #bh[i] <-sqrt((bh[i] + thi)^2 -2*bh[i]*thi-thi)
          
          bhj[j] <- bhj[j]+thj
          temp <-funcpa(bhj,fix=fix,...)
          #VmodelPar[k] <- -(temp-(fcith[j])-(fcith[i])+rl)/th
          VmodelPar[j,i] <- -(temp-(fcith[j])-(fcith[i])+rl)/th
          k<-k+1
          
        }
      }
      
      return(VmodelPar)
}

deriva <- function(nproc=1,b,funcpa,.packages=NULL,...){
  
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
    
    
    ## derivees premieres:
    ll <- foreach(k=(m*(m+1)/2)+1:m,
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
    
    fcith <- ll[1,]
    v1 <- ll[2,] 
    
    ## derivees secondes:
    v2 <- foreach(k=1:(m*(m+1)/2),
                  .combine=c,
                  .packages=.packages) %dopar%
      {
        i <- grid[k,1]
        j <- grid[k,2]
        bij <- b
        bij[i] <- bij[i]+h[i]
        bij[j] <- bij[j]+h[j]
        
        res <- -(funcpa(bij,...)-fcith[i]-fcith[j]+rl)/(h[i]*h[j])
        res
      }
    
    v <- c(v2,v1)
    
  }
  else
  {
    ## gradient null
    for(i in 1:m){
      bh <- bh2 <- b
      th <- max(1E-7,(1E-4*abs(b[i])))
      bh[i] <- bh[i] + th
      bh2[i] <- bh2[i] - th
      
      fcith[i] <- funcpa(bh,...)
      fcith2[i] <- funcpa(bh2,...)
    }
    
    k <- 0
    m1 <- m*(m+1)/2
    l <- m1
    for(i in 1:m){
      l <- l+1
      bh <- b
      thn <- - max(1E-7,(1E-4*abs(b[i])))
      v[l] <- -(fcith[i]-fcith2[i])/(2*thn)
      for(j in 1:i){
        bh <- b
        k <- k+1
        thi <- max(1E-7,(1E-4*abs(b[i])))
        thj <- max(1E-7,(1E-4*abs(b[j])))
        th <- thi * thj
        bh[i] <- bh[i]+thi
        bh[j] <- bh[j]+thj
        temp <-funcpa(bh,...)
        v[k] <- -(temp-(fcith[j])-(fcith[i])+rl)/th
      }
    }
  }
  
  result <- list(v=v,rl=rl)
  return(result)
}



grmlaweib<-function(b,npm,npar,bfix,fix,ctime,no,ve01,ve02,ve12,
                      dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                      t0,t1,t2,t3,troncature,lambda,alpha,penalty.factor,penalty,gausspoint,weib){
  
  #browser()
  
  # test<-deriva(b=b,funcpa=idmlLikelihoodweib,npm=length(b),
  #              npar=size_V,
  #              bfix=bfix,
  #              fix=fix0,
  #              ctime=ctime,
  #              no=N,
  #              ve01=ve01,
  #              ve02=ve02,
  #              ve12=ve12,
  #              dimnva01=dimnva01,
  #              dimnva02=dimnva02,
  #              dimnva12=dimnva12,
  #              nva01=nvat01,
  #              nva02=nvat02,
  #              nva12=nvat12,
  #              t0=t0,
  #              t1=t1,
  #              t2=t2,
  #              t3=t3,
  #              troncature=troncature,
  #              gausspoint=gausspoint,weib=weib)
  # 
  # 
  # testbis<-grad(x0=b,f=idmlLikelihoodweib,npm=length(b),
  #              npar=size_V,
  #              bfix=bfix,
  #              fix=fix0,
  #              ctime=ctime,
  #              no=N,
  #              ve01=ve01,
  #              ve02=ve02,
  #              ve12=ve12,
  #              dimnva01=dimnva01,
  #              dimnva02=dimnva02,
  #              dimnva12=dimnva12,
  #              nva01=nvat01,
  #              nva02=nvat02,
  #              nva12=nvat12,
  #              t0=t0,
  #              t1=t1,
  #              t2=t2,
  #              t3=t3,
  #              troncature=troncature,
  #              gausspoint=gausspoint,weib=weib)
  start<-sum(fix[1:6]==1)
  if(sum(fix[1:6])==1){
    svar<-NULL
  }else{svar<-b[1:(6-start)]}
  ball<-b[(6-start+1):(npm)]
  npm_all<-length(ball)
  grbeta<-rep(0,npm_all)
  
  fixbeta<-fix
  fixbeta[1:6]<-1
  bb<-rep(NA,npar)
  bb[which(fix==1)]<-bfix
  bb[which(fixbeta==1 & fix==0)]<-svar
  bb<-na.omit(bb)
  # b_positive<-pmax(ball,0)
  # return first and second derivatives of the loglik
  grbeta<-.Fortran("firstderivaweib",
                   ## input
                   as.double(ball),
                   as.integer(npm_all),
                   as.integer(npar),
                   as.double(bb),
                   as.integer(fixbeta),
                   as.integer(ctime),
                   as.integer(no),
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
                   as.integer(weib),
                   likelihood_deriv=as.double(grbeta),
                   PACKAGE="SmoothHazardoptim9")$likelihood_deriv
  
  if(any(grbeta==Inf)| any(grbeta==-Inf) | any(is.na(grbeta)) | any(is.nan(grbeta))){
    cat("Problem of computation on the analytical gradient, perform finite differencies. \n")
    sol<-deriva.gradient(b=b,
                         nproc=1,
                         funcpa=idmlLikelihoodweib,
                         npm=length(b),
                         npar=npar,
                         bfix=bfix,
                         fix=fix,
                         ctime=ctime,
                         no=no,
                         ve01=ve01,
                         ve02=ve02,
                         ve12=ve12,
                         dimnva01=dimnva01,
                         dimnva02=dimnva02,
                         dimnva12=dimnva12,
                         nva01=nva01,
                         nva02=nva02,
                         nva12=nva12,
                         t0=t0,
                         t1=t1,
                         t2=t2,
                         t3=t3,
                         troncature=troncature,
                         gausspoint=gausspoint,weib=weib)$v
    return(-sol)
  }else{
  bb<-rep(NA,npar)
  bb[fix==0]<-c(svar,ball)
  bb[fix==1]<-bfix
  
  fixbeta<-fix
  fixbeta[7:length(fixbeta)]<-1
  bb<-rep(NA,npar)
  bb[which(fix==1)]<-bfix
  bb[which(fixbeta==1 & fix==0)]<-ball
  bb<-na.omit(bb)
  
  grs<-deriva.gradient(b=svar,
                       nproc=1,
                       funcpa=idmlLikelihoodweib,
                       npm=length(svar),
                       npar=npar,
                       bfix=bb,
                       fix=fixbeta,
                       ctime=ctime,
                       no=no,
                       ve01=ve01,
                       ve02=ve02,
                       ve12=ve12,
                       dimnva01=dimnva01,
                       dimnva02=dimnva02,
                       dimnva12=dimnva12,
                       nva01=nva01,
                       nva02=nva02,
                       nva12=nva12,
                       t0=t0,
                       t1=t1,
                       t2=t2,
                       t3=t3,
                       troncature=troncature,
                       gausspoint=gausspoint,weib=weib)
  
  
  sol<-c(-grs$v,grbeta)
  if(any(sol==Inf)| any(sol==-Inf) | any(is.na(sol)) | any(is.nan(sol))){
    stop(paste0("Problem of computation on the gradient. Verify your function specification...\n. Infinite value with finite parameters : b=",round(b,4),"\n"))
  }
  return(as.double(sol))
  }
}





hessianmlaweib<-function(b,npm,npar,bfix,fix,ctime,no,ve01,ve02,ve12,
                    dimnva01,dimnva02,dimnva12,nva01,nva02,nva12,
                    t0,t1,t2,t3,troncature,gausspoint,weib){
  

  start<-sum(fix[1:6]==1)
  
  if(sum(fix[1:6])==1){
    svar<-NULL
  }else{svar<-b[1:(6-start)]}
  
  ball<-b[(6-start+1):(npm)]
  npm_all<-length(ball)
  output<-rep(0,(npm_all*(npm_all+1)/2)+npm_all)
  
  fixbeta<-fix
  fixbeta[1:6]<-1
  
  bb<-rep(NA,npar)
  bb[which(fix==1)]<-bfix
  bb[which(fixbeta==1 & fix==0)]<-svar
  bb<-na.omit(bb)
  
  output<-.Fortran("derivaweib",
                   ## input
                   as.double(ball),
                   as.integer(npm_all),
                   as.integer(npar),
                   as.double(bb),
                   as.integer(fixbeta),
                   as.integer(ctime),
                   as.integer(no),
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
                   as.integer(weib),
                   likelihood_deriv=as.double(output),
                   PACKAGE="SmoothHazardoptim9")$likelihood_deriv
  
  if(any(output[(npm_all+1):length(output)]==Inf)| any(output[(npm_all+1):length(output)]==-Inf) | any(is.na(output[(npm_all+1):length(output)])) | any(is.nan(output[(npm_all+1):length(output)]))){
    cat("Problem of computation on the analytical hessian, thus compute finite differences.\n")
    Vall<-deriva(b=b,f=idmlLikelihoodweib,npm=length(b),
                             npar=npar,
                             bfix=bfix,
                             fix=fix,
                             ctime=ctime,
                             no=no,
                             ve01=ve01,
                             ve02=ve02,
                             ve12=ve12,
                             dimnva01=dimnva01,
                             dimnva02=dimnva02,
                             dimnva12=dimnva12,
                             nva01=nva01,
                             nva02=nva02,
                             nva12=nva12,
                             t0=t0,
                             t1=t1,
                             t2=t2,
                             t3=t3,
                             troncature=troncature,
                             gausspoint=gausspoint,weib=weib)
    
    if(any(Vall$v==Inf)| any(Vall$v==-Inf) | any(is.na(Vall$v)) | any(is.nan(Vall$v))){
      stop(paste0("Problem of computation on the hessian with finite differences. Verify your function specification...\n.
                  Infinite value with finite parameters : b=",round(b,4),"\n"))
    }
    return(Vall$v[1:(length(b)*(length(b)+1)/2)])
  }else{
  min<-npm_all
  V01<- matrix(0,nva01,nva01)
  V01[lower.tri(V01,diag=TRUE)] <- output[(min+1):(min+nva01*(nva01+1)/2)]
  
  
  min<-min+(nva01*(nva01+1)/2)
  if(nva01>0&nva02>0){
    V0102<- matrix(data=output[(min+1):(min+nva02*nva01)],
                   nrow=nva02,ncol=nva01)
  }
  
  min<-min+nva02*nva01
  
  if(nva01>0&nva12>0){
    V0112<- matrix(data=output[(min+1):(min+nva12*nva01)],
                   nrow=nva12,ncol=nva01)
  }
  
  
  min<-min+nva12*nva01
  V02<- matrix(0,nva02,nva02)
  V02[lower.tri(V02,diag=TRUE)] <- output[(min+1):(min+nva02*(nva02+1)/2)]
  
  
  min<-min+(nva02*(nva02+1)/2)
  
  if(nva02>0&nva12>0){
    V0212<- matrix(data=output[(min+1):(min+nva12*nva02)],
                   nrow=nva12,ncol=nva02)
    }
  
  
  min<-min+nva12*nva02
  V12<- matrix(0,nva12,nva12)
  V12[lower.tri(V12,diag=TRUE)] <- output[(min+1):length(output)]
  
  
  V<- matrix(0,npm_all,npm_all)
  if(nva01>0){
    V[1:nva01,1:nva01]<-V01
    if(nva02>0){V[(nva01+1):(nva01+nva02),1:nva01]<-V0102}
    if(nva12>0){V[(nva01+nva02+1):npm_all,1:nva01]<-V0112}
  }
  if(nva02>0){
    V[(nva01+1):(nva01+nva02),(nva01+1):(nva01+nva02)]<-V02
    if(nva12>0){V[(nva01+nva02+1):npm_all,(nva01+1):(nva01+nva02)]<-V0212}
  }
  
  if(nva12>0){
    V[(nva01+nva02+1):npm_all,(nva01+nva02+1):(npm_all)]<-V12
    }
  
  # hessian is - second derivatives 
  #V<-V+t(V)
  #diag(V)<-diag(V)/2
  # hessian is - second derivatives 
  V<--V
  #browser()
  if(start==6){
    return(V)
    }else{
      Vall<-deriva.hessianweib(b=b,f=idmlLikelihoodweib,npm=length(b),
              npar=npar,
              bfix=bfix,
              fix=fix,
              ctime=ctime,
              no=no,
              ve01=ve01,
              ve02=ve02,
              ve12=ve12,
              dimnva01=dimnva01,
              dimnva02=dimnva02,
              dimnva12=dimnva12,
              nva01=nva01,
              nva02=nva02,
              nva12=nva12,
              t0=t0,
              t1=t1,
              t2=t2,
              t3=t3,
              troncature=troncature,
              gausspoint=gausspoint,weib=weib)
      Vall[c((length(svar)+1):dim(Vall)[1]),c((length(svar)+1):dim(Vall)[1])]<-t(V)
    
      if(any(Vall==Inf)| any(Vall==-Inf) | any(is.na(Vall)) | any(is.nan(Vall))){
        stop(paste0("Problem of computation on the hessian for weibull parameters. Verify your function specification...\n.
        Infinite value with finite parameters : b=",round(b,4),"\n"))
     
      }
      
      
      
  return(Vall[upper.tri(Vall,diag=T)])
     
    }
  }

 
}

