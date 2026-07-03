### Code:
#' @title Covariance for an illness-death model 
#' @description
#' Plot estimated baseline transition intensities from an object of class
#' \code{idm} optionally with confidence limits.
#' @param x a \code{idmWeib} class object (output from calling
#' \code{idm} with the (default) option \code{intensities}="Weib".
#' @param formula01 A formula specifying a regression model for the
#' \code{0 --> 1} transition from the initial state to the transient
#' state of the illness-death model.  The right hand side of the
#' formula specifies the covariate terms, and the left hand side must
#' be an event history object as returned by the function \code{Hist}.
#' @param formula02 A formula specifying a regression model for the
#' \code{0 --> 2} transition from the initial state to the absorbing
#' state. The left hand side must be equal to the left hand side of
#' \code{formula01}. If missing it is set to \code{formula01}.
#' @param formula12 A formula specifying a regression model for the
#' \code{1 --> 2} transition from the transient state to the absorbing
#' state.  operator is not required. If missing it is set to
#' \code{formula01}.
#' @param data A data frame in which to interpret the variables of
#' \code{formula01}, \code{formula02} and \code{formula12}.
#' @return Return the covariance matrix and the first derivatives of the parameters accounting for the selection process if needed 
#' @seealso
#' \code{\link{print.idm}},\code{\link{summary.idm}},\code{\link{idm}},
#' @seealso \code{\link{idm}}
#' @keywords methods
#'
#'@useDynLib HIDeM
#' @export
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr> 
covIDM <- function(object,
                    data,
                    formula01,
                    formula02,
                    formula12,
                   lambda="BIC",
                   only.beta=F){ 
  
  call <- match.call()
  # }}}
  # {{{ evaluate formula in data
  if(missing(object)) stop("Need an object HIDeM.")
  if(missing(data)) stop("Need a data frame.")
  if(sum(is.na(data))>0)stop("Need a data frame with no missing data.")
  
  if(!inherits(data,"data.frame"))stop("Argument 'data' must be a data.frame")
  
  if(missing(formula01))stop("Argument formula01 is missing.")
  if(missing(formula02))stop("Argument formula02 is missing.")
  if(!inherits(formula01,"formula"))stop("The argument formula01 must be a class 'formula'.")
  if(!inherits(formula02,"formula"))stop("The argument formula02 must be a class 'formula'.")
  if(missing(formula12)) formula12 <- formula02
  
  
  if(!is.null(object$modelPar)){
    object$method<-"weib"
  }else{object$method<-"splines"
  }
  
  #############################################################################
  ########################### get baseline ####################################
  #############################################################################
  
  if(object$penalty!="none"){
    if(is.null(lambda)){lambda<-"BIC"}
    if(length(lambda)==1){if(!lambda%in%c("GCV","BIC")){stop("Lambda need to be either a vector of three values (01,02 and 12) or BIC or GCV")}}
    if(!length(lambda)%in%c(1,3)){stop("Lambda need to be either a vector of three values (01,02 and 12) or BIC or GCV")}
    
    if(length(lambda)==3){
      if(!any(apply(object$lambda,FUN=function(x){sum(x==lambda)},MARGIN=2)==3)){stop("Lambda need to be either a vector of three values (01,02 and 12) from object$lambda")}
      id<-which(apply(object$lambda,FUN=function(x){sum(x==lambda)},MARGIN=2)==3)[1]
    }
    if(length(lambda)==1){
      if(lambda=="BIC"){
        BIC<-min(object$BIC[object$converged==1])
        id<-which(object$BIC==BIC)[1]
      }else{
        GCV<-min(object$GCV[object$converged==1])
        id<-which(object$GCV==GCV)[1]}
    }
  }
  # even if parameters are fixed they are given back
  if (object$method=="splines"){
    
    nknots01 <- object$nknots01
    nknots02 <- object$nknots02
    nknots12 <- object$nknots12
    knots01 <- object$knots01
    knots02 <- object$knots02
    knots12 <- object$knots12
    if(object$penalty=="none"){
    basepar <- c(object$theta01,
                 object$theta02,
                 object$theta12)
    }else{
      basepar<-c(object$theta01[,id],
                 object$theta02[,id],
                 object$theta12[,id])
    }
    size_spline<-nknots01+nknots02+nknots12 + 6
  }else {
    ############################ weibull #####################################
    if(object$penalty=="none"){
    basepar <- object$modelPar
    }else{
      basepar <- object$modelPar[,id]
    }
    size_spline<-6
  }
  
  sm_arg <- object$call$semiMarkov
  is_semiMarkov <- isTRUE(sm_arg) || identical(sm_arg, quote(T)) || identical(sm_arg, TRUE)
  semiMark <- as.numeric(is_semiMarkov)
  
  ############################################################################
  ############################### get database defined by formulas ###########
  ############################################################################
  
  ################################################################################
  ########################## distinguish 2 cases #################################
  ################### - penalized parameters #####################################
  ################### - reestimated parameters ###################################
  ################################################################################
  
  # if reestimated parameters, we need to know the formula before the reestimation, 
  # if null we would consider that both are the same and we do not add zero coef where 
  # variable not selected #########################################################

 
  if(object$penalty!="none"){

  if(!formulas_equal_terms(object$terms$Formula01,formula01))stop("Formula 01 do not match the object formula of HIDeM object")
  if(!formulas_equal_terms(object$terms$Formula02,formula02))stop("Formula 02 do not match the object formula of HIDeM object")
  if(!formulas_equal_terms(object$terms$Formula12,formula12))stop("Formula 12 do not match the object formula of HIDeM object")

  }else{
    
    m <- match.call()
    m01 <- m02 <- m12 <- m[match(c("","data","na.action"),names(m),nomatch=0)]
    
    #get number of variables on each transition based on formula terms 
    m01$formula <- object$terms$Formula01
    m02$formula <- object$terms$Formula02
    m12$formula <- object$terms$Formula12

    m01[[1]] <- m02[[1]] <- m12[[1]] <- as.name("model.frame")
    m01 <- eval(m01,parent.frame())
    m02 <- eval(m02,parent.frame())
    m12 <- eval(m12,parent.frame())
    
    
    x01 <- model.matrix(object$terms$Formula01,data=m01)[, -1, drop = FALSE]
    p01 <- NCOL(x01)
    
    x02 <- model.matrix(object$terms$Formula02,data=m02)[, -1, drop = FALSE]
    p02 <- NCOL(x02)
    
    x12 <- model.matrix(object$terms$Formula12,data=m12)[, -1, drop = FALSE]
    p12 <- NCOL(x12)
    
    
  }
    
  m <- match.call()
  m01 <- m02 <- m12 <- m[match(c("","data","na.action"),names(m),nomatch=0)]
  
  
    m01$formula <- formula01
    m02$formula <- formula02
    m12$formula <- formula12
    
    
  
  m01[[1]] <- m02[[1]] <- m12[[1]] <- as.name("model.frame")
    if(anyNA(data)){
      variables=unique(c(all.vars(formula01),all.vars(formula02),all.vars(formula12)))
      data=data[,variables]
      data=na.omit(data)
      m01[[2]] <- m02[[2]] <- m12[[2]] <- data
    }
    
    m01 <- eval(m01,parent.frame())
    m02 <- eval(m02,parent.frame())
    m12 <- eval(m12,parent.frame())
    
    
    responseTrans <- stats::model.response(m01)
    responseAbs <- stats::model.response(m02)
    
    #################################################################################
    #####################   extract covariates   ####################################
    #################################################################################
    
    x01 <- model.matrix(formula01,data=m01)[, -1, drop = FALSE]
    NC01 <- NCOL(x01)
    
    if (NC01>0)
      Xnames01 <- colnames(x01)
    else
      Xnames01 <- NULL
    ## formula02
    x02 <- model.matrix(formula02,data=m02)[, -1, drop = FALSE]
    NC02 <- NCOL(x02)
    
    
    if (NC02>0)
      Xnames02 <- colnames(x02)
    else
      Xnames02 <- NULL

    
    sm_arg <- object$call$semiMarkov
    is_semiMarkov <- isTRUE(sm_arg) || identical(sm_arg, quote(T)) || identical(sm_arg, TRUE)
    semiMark <- as.numeric(is_semiMarkov)
    
    x12 <- model.matrix(formula12,data=m12)[, -1, drop = FALSE]
    NC12 <- NCOL(x12)
    
    
    if (NC12>0)
      Xnames12 <- colnames(x12) 
    else
        Xnames12 <- NULL
    
  
  #################################################################################
  ####################  prepare censored event times  #############################
  #################################################################################
  
  isIntervalCensored <- attr(responseTrans,"cens.type")=="intervalCensored"
  truncated <- nchar(attr(responseAbs,"entry.type"))>1
  abstime <- as.double(responseAbs[,"time"])
  ## It may happen that the illness time is observed exactly, in which case
  ## the status is 1, thus we need two criteria to declare illness status:
  ## 1. exact observations with illness status ==1
  ## 2. interval censored with any illness status. FIXME: check the corresponding likelihood
  
  idm <- responseTrans[,"status"]==(as.integer(isIntervalCensored)+1)
  if (isIntervalCensored)
    idm[(responseTrans[,"status"]==1 & (responseTrans[,"L"]==responseTrans[,"R"]))] <- 1
  ## exit status
  idd <- responseAbs[,"status"]==1
  
  
  N <- length(abstime)
  if (truncated==0){
    entrytime <- as.double(NULL)
  }else{
    entrytime <- as.double(responseAbs[,"entry"])
  }
  if (isIntervalCensored){
    Ltime <- as.double(responseTrans[,"L",drop=TRUE])
    Rtime <- as.double(responseTrans[,"R",drop=TRUE])
    ## if (any(Rtime<abstime & idm ==0))
    ## warning(paste("For ",
    ## sum(Rtime<abstime & idm ==0),
    ## " cases where the ill status is not observed\n and the last inspection time (R) is smaller than the right censored time (T)\n the time R is set to T."))
  }else{# exactly observed transition times
    Ltime <- as.double(responseTrans[,"time",drop=TRUE])
    Rtime <- as.double(responseTrans[,"time",drop=TRUE])
    Ltime[idm==0] <- abstime[idm==0]
    Rtime[idm==0] <- abstime[idm==0]
  }

  
  
  if (attr(responseAbs,"cens.type")=="intervalCensored") stop("No method available when the transtion to the absorbing state is interval censored.")
  if (isIntervalCensored && any(Rtime<Ltime)) stop("Misspecified transitition times:\nSome left interval limits are greater than the corresponding right limits.")
  
  #################################################################################
  ########################### Define number of parameters per transitions #########
  #################################################################################
  
  size1 <- NC01 + NC02 + NC12
  size_V <- size1 + size_spline
  
  noVar<-c(ifelse(as.integer(NC01)>0,0,1),
           ifelse(as.integer(NC02)>0,0,1),
           ifelse(as.integer(NC12)>0,0,1))
  
  
  nvat01 <- ifelse(noVar[1]==1,0,NC01)
  nvat02 <- ifelse(noVar[2]==1,0,NC02)
  nvat12 <- ifelse(noVar[3]==1,0,NC12)
  
  dimnva01<-ifelse(nvat01==0,1,nvat01)
  dimnva02<-ifelse(nvat02==0,1,nvat02)
  dimnva12<-ifelse(nvat12==0,1,nvat12)
  
  
  NC<-c(NC01,NC02,NC12)
  
  if(noVar[1]==1){ve01<-as.double(rep(0,N))}else{ve01<-as.double(x01)}
  if(noVar[2]==1){ve02<-as.double(rep(0,N))}else{ve02<-as.double(x02)}
  if(noVar[3]==1){ve12<-as.double(rep(0,N))}else{ve12<-as.double(x12)}
  
  #################################################################################
  ####################### define profile of subjects ##############################
  #################################################################################
  
  troncature<-ifelse(truncated==T,1,0)
  
  if(truncated==1){
    t0<-entrytime
  }else{t0<-rep(0,N)}
  
  t1<-Ltime
  t2<-Rtime
  t3<-abstime
  t4<-rep(NA,N)
  ctime<-rep(NA,N)
  
  
  ctime<-ifelse(idm==0 & idd==0 & t1==t3,1,NA)
  ctime<-ifelse(idm==1 & idd==0 & t1<t2,2,ctime)
  ctime<-ifelse(idm==1 & idd==0 & t1==t2,3,ctime)
  ctime<-ifelse(idm==1 & idd==1 & t1<t2,4,ctime)
  ctime<-ifelse(idm==1 & idd==1 & t1==t2,5,ctime)
  ctime<-ifelse(idm==0 & idd==0 & t1<t3,6,ctime)
  ctime<-ifelse(idm==0 & idd==1,7,ctime)
  
  if(sum(is.na(ctime))>0){stop("For subject with no event, time for event 01 cannot be equal to time for 12 and 02")}
  
  t2<-ifelse(ctime==1 | ctime==3 | ctime==5,t1,
             ifelse(ctime==2 | ctime==4,t2,
                    ifelse(ctime==6 | ctime==7, t3,NA)
             )
  )
  
  t3<-ifelse(ctime==1, t1,
             ifelse(ctime==2 | ctime==3 | ctime==4 | ctime==5 | ctime==6 | ctime==7,t3,NA))
  
  t4<-ifelse(ctime==1 | ctime==2 | ctime==3 | ctime==4 | ctime==5, t1,
             ifelse( ctime==6 | ctime==7,t3,NA))
  
  
  ######################## initialize parameters ###############################
  
  if(object$penalty!="none"){
    if(only.beta==F){
    b<-c(basepar,object$coef[,id])
    }else{
      bfix<-basepar
      fix0<-rep(0,length(c(basepar,object$coef[,id])))
      fix0[1:length(basepar)]<-1
      b<-object$coef[,id]
    }
  }else{
    
    b01<-rep(0,nvat01)
    b02<-rep(0,nvat02)
    b12<-rep(0,nvat12)
    
    if(nvat01 >0){
    ######################### for 01 #########################################
    # Extract term labels from each formula object

    terms1 <- attr(terms(object$terms$Formula01), "term.labels")
    terms2 <- attr(terms(formula01), "term.labels")
    
    
    # Common variables
    common <- intersect(terms1, terms2)
    
    # Variables only in object$terms$Formula01
    only_in_1 <- setdiff(terms1, terms2)
    
    # Variables only in formula01
    only_in_2 <- setdiff(terms2, terms1)
    
    if(length(only_in_1)>0){stop("Need to have in formula01 all the variables used for the model")}
    
    if(length(only_in_2)>0){
      idx <- which(!(terms2 %in% only_in_2))
      b01[idx]<-object$coef[1:p01]
    }else{
      b01<-object$coef[1:p01]
      
    }
      
    
    }else{
      b01<-NULL
    }
  
    
    ######################### for 02 #########################################
    # Extract term labels from each formula object
    
    if(nvat02>0){
      
    terms1 <- attr(terms(object$terms$Formula02), "term.labels")
    terms2 <- attr(terms(formula02), "term.labels")
    
    
    # Common variables
    common <- intersect(terms1, terms2)
    
    # Variables only in object$terms$Formula02
    only_in_1 <- setdiff(terms1, terms2)
    
    # Variables only in formula02
    only_in_2 <- setdiff(terms2, terms1)
    
    if(length(only_in_1)>0){stop("Need to have in formula02 all the variables used for the model")}
    
    if(length(only_in_2)>0){
        idx <- which(!(terms2 %in% only_in_2))
        b02[idx]<-object$coef[(p01+1):(p01+p02)]
    }else{
        b02<-object$coef[(p01+1):(p01+p02)]
      
    }
        
    
    }else{
      b02<-NULL
    }
    
    ######################### for 12 #########################################
    # Extract term labels from each formula object
    
    if(nvat12>0){
      
        
    terms1 <- attr(terms(object$terms$Formula12), "term.labels")
    terms2 <- attr(terms(formula12), "term.labels")
    
    
    # Common variables
    common <- intersect(terms1, terms2)
    
    # Variables only in object$terms$Formula12
    only_in_1 <- setdiff(terms1, terms2)
    
    # Variables only in formula12
    only_in_2 <- setdiff(terms2, terms1)
    
    if(length(only_in_1)>0){stop("Need to have in formula12 all the variables used for the model")}
    
    if(length(only_in_2)>0){
      idx <- which(!(terms2 %in% only_in_2))
      b12[idx]<-object$coef[(p01+p02+1):(p01+p02+p12)]
    }else{
      b12<-object$coef[(p01+p02+1):(p01+p02+p12)]
      
    }
      
    
    }else{
      b12<-NULL
    }
    
    if(only.beta==F){
    b<-c(basepar,b01,b02,b12)
    }else{
      bfix<-basepar
      fix0<-rep(0,length(c(basepar,b01,b02,b12)))
      fix0[1:length(basepar)]<-1
      b<-c(b01,b02,b12)
    }
  }


# do finite differencies 
  if((only.beta==F)|(semiMark==T)){
    if(object$method=="weib"){
     output<-deriva(h=1e-8,b=b,
                           funcpa=idmlLikelihoodweib,
                           npm=length(b),
                           npar=size_V,
                           bfix=1,
                           fix=rep(0,size_V),
                           ctime=ctime,
                           no=N,
                           ve01=ve01,
                           ve02=ve02,
                           ve12=ve12,
                           dimnva01=dimnva01,
                           dimnva02=dimnva02,
                           dimnva12=dimnva12,
                           nva01=nvat01,
                           nva02=nvat02,
                           nva12=nvat12,
                           t0=t0,
                           t1=t1,
                           t2=t2,
                           t3=t3,
                           troncature=troncature,
                           gausspoint=15)
    }else{
      
     output<-deriva(h=.Machine$double.eps,
                 b=b,
          funcpa=idmlLikelihood,
          npm=length(b),
          npar=size_V,
          bfix=1,
          fix=rep(0,size_V),
          zi01=knots01,
          zi02=knots02,
          zi12=knots12,
          ctime=ctime,
          no=N,
          nz01=nknots01,
          nz02=nknots02,
          nz12=nknots12,
          ve01=ve01,
          ve02=ve02,
          ve12=ve12,
          dimnva01=dimnva01,
          dimnva02=dimnva02,
          dimnva12=dimnva12,
          nva01=nvat01,
          nva02=nvat02,
          nva12=nvat12,
          t0=t0,
          t1=t1,
          t2=t2,
          t3=t3,
          troncature=troncature,
          gausspoint=15)
      
      
    }
    
    fu <- output$v[((size_V*(size_V+1)/2)+1):(length(output$v))]
    
    V <- matrix(0, size_V, size_V)
    V[upper.tri(V, diag = TRUE)] <- output$v[1:(size_V*(size_V+1)/2)]
    V <- V + t(V)
    diag(V) <- diag(V)/2
    
    if(any(eigen(V)$values<0)){
      natureV<-"Fischer"
    }else{
      V <- solve(V)
      natureV<-"Cov"
    }
   
    
  }else{
    if((semiMark==F) & (only.beta==F)){
      
      if(object$method=="weib"){
        output<-deriva(h=.Machine$double.eps,b=b,
                    funcpa=idmlLikelihoodweib,
                    npm=length(b),
                    npar=size_V,
                    bfix=1,
                    fix=rep(0,size_V),
                    ctime=ctime,
                    no=N,
                    ve01=ve01,
                    ve02=ve02,
                    ve12=ve12,
                    dimnva01=dimnva01,
                    dimnva02=dimnva02,
                    dimnva12=dimnva12,
                    nva01=nvat01,
                    nva02=nvat02,
                    nva12=nvat12,
                    t0=t0,
                    t1=t1,
                    t2=t2,
                    t3=t3,
                    troncature=troncature,
                    gausspoint=15)
      }else{
        
        output<-deriva(h=.Machine$double.eps,
                    b=b,
                    npm=length(b),
                    funcpa=idmlLikelihood,
                    npar=size_V,
                    bfix=1,
                    fix=rep(0,size_V),
                    zi01=knots01,
                    zi02=knots02,
                    zi12=knots12,
                    ctime=ctime,
                    no=N,
                    nz01=nknots01,
                    nz02=nknots02,
                    nz12=nknots12,
                    ve01=ve01,
                    ve02=ve02,
                    ve12=ve12,
                    dimnva01=dimnva01,
                    dimnva02=dimnva02,
                    dimnva12=dimnva12,
                    nva01=nvat01,
                    nva02=nvat02,
                    nva12=nvat12,
                    t0=t0,
                    t1=t1,
                    t2=t2,
                    t3=t3,
                    troncature=troncature,
                    gausspoint=15)
        
        
      }
      
      fu <- output$v[((size_V*(size_V+1)/2)+1):(length(output$v))]
      
      V <- matrix(0, size_V, size_V)
      V[upper.tri(V, diag = TRUE)] <- output$v[1:(size_V*(size_V+1)/2)]
      V <- V + t(V)
      diag(V) <- diag(V)/2
      if(any(eigen(V)$values<0)){
        natureV<-"Fischer"
      }else{
        V <- solve(V)
        natureV<-"Cov"
      }
      
    }else{
      
      if(semiMark==T){
        
        if(object$method=="weib"){
          output<-derivaweibsemiMarkov(b=b,
                      npm=length(b),
                      npar=size_V,
                      bfix=bfix,
                      fix=fix0,
                      ctime=ctime,
                      no=N,
                      ve01=ve01,
                      ve02=ve02,
                      ve12=ve12,
                      dimnva01=dimnva01,
                      dimnva02=dimnva02,
                      dimnva12=dimnva12,
                      nva01=nvat01,
                      nva02=nvat02,
                      nva12=nvat12,
                      t0=t0,
                      t1=t1,
                      t2=t2,
                      t3=t3,
                      troncature=troncature)
        }else{
          
          output<-derivasplinesemiMarkov(b=b,
                      npm=length(b),
                      npar=size_V,
                      bfix=bfix,
                      fix=fix0,
                      zi01=knots01,
                      zi02=knots02,
                      zi12=knots12,
                      ctime=ctime,
                      no=N,
                      nz01=nknots01,
                      nz02=nknots02,
                      nz12=nknots12,
                      ve01=ve01,
                      ve02=ve02,
                      ve12=ve12,
                      dimnva01=dimnva01,
                      dimnva02=dimnva02,
                      dimnva12=dimnva12,
                      nva01=nvat01,
                      nva02=nvat02,
                      nva12=nvat12,
                      t0=t0,
                      t1=t1,
                      t2=t2,
                      t3=t3,
                      troncature=troncature)
          
          
        }
        
        fu <- output[1:length(b)]
        V<- matrix(0,length(b),length(b))
        V[lower.tri(V,diag=TRUE)] <- output[(length(b)+1):length(output)]
        V<-V+t(V)
        diag(V)<-diag(V)/2
        V<--V
        if(any(eigen(V)$values<0)){
          natureV<-"Fischer"
        }else{
          V <- solve(V)
          natureV<-"Cov"
        }
      }else{
        
        if(object$method=="weib"){
          output<-derivaweib(b=b,
                                    npm=length(b),
                                    npar=size_V,
                                    bfix=bfix,
                                    fix=fix0,
                                    ctime=ctime,
                                    no=N,
                                    ve01=ve01,
                                    ve02=ve02,
                                    ve12=ve12,
                                    dimnva01=dimnva01,
                                    dimnva02=dimnva02,
                                    dimnva12=dimnva12,
                                    nva01=nvat01,
                                    nva02=nvat02,
                                    nva12=nvat12,
                                    t0=t0,
                                    t1=t1,
                                    t2=t2,
                                    t3=t3,
                                    troncature=troncature)
        }else{
          
          output<-derivaspline(b=b,
                                      npm=length(b),
                                      npar=size_V,
                                      bfix=bfix,
                                      fix=fix0,
                                      zi01=knots01,
                                      zi02=knots02,
                                      zi12=knots12,
                                      ctime=ctime,
                                      no=N,
                                      nz01=nknots01,
                                      nz02=nknots02,
                                      nz12=nknots12,
                                      ve01=ve01,
                                      ve02=ve02,
                                      ve12=ve12,
                                      dimnva01=dimnva01,
                                      dimnva02=dimnva02,
                                      dimnva12=dimnva12,
                                      nva01=nvat01,
                                      nva02=nvat02,
                                      nva12=nvat12,
                                      t0=t0,
                                      t1=t1,
                                      t2=t2,
                                      t3=t3,
                                      troncature=troncature)
          
          
          
          
        }
        fu <- output[1:length(b)]
        V<- matrix(0,length(b),length(b))
        V[lower.tri(V,diag=TRUE)] <- output[(length(b)+1):length(output)]
        V<-V+t(V)
        diag(V)<-diag(V)/2
        V<--V
        if(any(eigen(V)$values<0)){
          natureV<-"Fischer"
        }else{
          V <- solve(V)
          natureV<-"Cov"
        }
      }
    }
    
  }
    res<-list(fu=fu,
              V=V,
              natureV=natureV)
  return(res)
  
}


formulas_equal_terms <- function(t1, t2) {
  identical(attr(t1, "term.labels"), attr(t2, "term.labels"))
}

