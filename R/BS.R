### Code:
##' @title Calculate predictions for time-depend covariates using INLA
##' @param pred A vector containing the prediction for the probability of having the illness between s and horizon time
##' @param objectSurvival A idm object from HIDeM package containing the illness-death model estimation with time-fixed covariates 
##' @param newdata The newdata on which we want to calculate the BS, note : it must be the same as the one calculated for the predictions and in the same order as the predictions 
##' @param s entry time 
##' @param horizon horizon time 
##' @param envir working environment 
#' @importFrom foreach "%do%"
#' @importFrom foreach "%dopar%"
#' @importFrom Deriv "Deriv"
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  
#' @export


BS<-function(pred,
             objectSurvival,
             formula01,
             formula02,
             formula12,
             id,
             s,
             horizon,
             newdata,
             info.death=T,
             info.after.horizon=F,
             envir=parent.frame(),
             k=get("k",envir=envir),
             ncores=NULL)
{
  
  
  call <- match.call()
  ptm <- proc.time()
  
  ids_new<-NULL
  
  if(missing(pred)){stop("Need to specify a vector of predictions in pred")}
  if(!(class(pred)%in%"numeric")){stop("Need to specify a vector of numeric predictions in pred")}
  if(missing(objectSurvival)){stop("Need to specify objectSurvival as a idm or transitions predictions given by times elements of DYNpredIDM object")}
  if(!inherits(objectSurvival,c("list","idm"))){stop("Need to specify objectSurvival as a idm or transitions predictions given by times elements of DYNpredIDM object")}
  if(inherits(objectSurvival,"list")){
    
    if(missing(formula01)){stop("Need to formula01 left side if objectSurvival is a list")}
    if(missing(formula02)){stop("Need to formula02 left side if objectSurvival is a list")}
    if(missing(formula12)){stop("Need to formula12 left side if objectSurvival is a list")}
    if(missing(id)){stop("Need to id variable of newdata")}
    
    
    if(length(objectSurvival)!=15){stop("objectSurvival must be a list of 15 elements")}

    expected <- c("start", "end", "id", "p00", "p01", "p02_0","p02_1","p11","p12","prediction_p00","prediction_p01", "prediction_p02_0","prediction_p02_1","prediction_p11","prediction_p12")
    
    if (!identical(names(objectSurvival), expected)) {
      stop(
        "Invalid names in list.\n",
        "Expected: ", paste(expected, collapse = ", "), "\n",
        "Found:    ", paste(names(objectSurvival), collapse = ", ")
      )
    }
    if(!inherits(id,"character")|| length(id) != 1){stop("id must be a character")}
    if (!id %in% colnames(newdata)) stop("id '", id, "' is not a column of newdata")
    forms <- list(formula01 = formula01, formula02 = formula02)
    for (nm in names(forms)) {
      f <- forms[[nm]]
      
      if (!inherits(f, "formula")) stop(nm, " must be a formula")
      if (length(f) != 3) stop(nm, " must have a left-hand side (response)")
      
      lhs_vars <- all.vars(f[[2]])             # variables on the left side only
      miss <- setdiff(lhs_vars, colnames(newdata))
      if (length(miss) > 0)
        stop("Variable(s) on the left side of ", nm, " not found in newdata: ",
             paste(miss, collapse = ", "))
    }
    for(nm in c("p00","p01","p02_0","p02_1","p11","p12"))
      if(!is.logical(objectSurvival[[nm]])) stop("objectSurvival$", nm, " must be a logical")
    ids_model <- objectSurvival[["id"]]
    
    if (is.null(ids_model)) {
      stop("objectSurvival does not contain an 'id' element")
    }
    
    ids_new <- newdata[[id]]
    miss <- setdiff(ids_model, ids_new)
    
    if (length(miss) > 0) {
      stop("All id in objectSurvival must be in newdata. Missing: ",
           paste(head(miss, 10), collapse = ", "),
           if (length(miss) > 10) ", ..." else "")
    }
    lens <- lengths(objectSurvival)
    if (length(unique(lens)) != 1) {
      stop(
        "All elements must have the same length.\n",
        "Lengths found: ", paste(names(objectSurvival)[1:15], lens, sep = " = ", collapse = ", ")
      )
    }
    
    
    if(!inherits(objectSurvival$start,c("numeric","integer")))stop("objectSurvival$start need to be an integer or numeric")
    if(any(objectSurvival$start < 0))stop("objectSurvival$start need to be numeric superior or equal to 0 with s < horizon")
    if(!inherits(objectSurvival$end,c("numeric","integer")))stop("objectSurvival$end need to be an integer or numeric")
    if(any(objectSurvival$end < 0))stop("objectSurvival$end need to be numeric superior or equal to 0 with s < horizon")
    
  }
  
  
  
  if(missing(newdata)){stop("Need to provide newdata as a data.frame")}
  if(!inherits(newdata,"data.frame")){stop("Need to provide newdata as a data.frame")}
  
  if(sum(is.na(newdata))>0)stop("Need a new data frame with no missing data.")
  if(missing(s)|missing(horizon))stop("Need to specify s and horizon")
  
  if(!inherits(s,c("numeric","integer")))stop("s need to be an integer or numeric")
  if(length(s)!=1)stop("Length of s need to be 1")
  if((s <= 0))stop("s need to be numeric superior to 0 with s < horizon")
  
  
  if(!inherits(horizon,c("numeric","integer")))stop("horizon need to be an integer or numeric")
  if(length(horizon)!=1)stop("Length of horizon need to be 1")
  if(is.null(s))stop("landmark time need to be provided")
  if((s <= 0) | (horizon <= 0) | (s >= horizon))stop("s and horizon need to be numeric superior to 0 with s < horizon")
  
  newdata<-unique(newdata)
  if(dim(newdata)[1]!=length(pred)){stop("The number of lines in newdata needs to be the size of pred")}
  N<-dim(newdata)[1]
  
  
  #################################################################################
  ####################  prepare censored event times  #############################
  #################################################################################
  
  if(inherits(objectSurvival,"idm")){
  m01 <- model.frame(objectSurvival$terms$Formula01, data = newdata)
  m02 <- model.frame(objectSurvival$terms$Formula02, data = newdata)
  m12 <- model.frame(objectSurvival$terms$Formula12, data = newdata)
  }else{
    m01 <- model.frame(formula01, data = newdata)
    m02 <- model.frame(formula02, data = newdata)
    m12 <- model.frame(formula12, data = newdata)
  }
  
  responseTrans <- stats::model.response(m01)
  responseAbs <- stats::model.response(m02)
  
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
  
  
  #N <- length(abstime)
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
  ## find time boundaries
  if (length(entrytime)>0){
    alltimes <- sort(unique(c(Ltime, Rtime,entrytime,abstime)))
    amax <- max(alltimes)
    amin <- min(alltimes)
  }
  else{
    alltimes <- sort(unique(c(Ltime, Rtime,abstime)))
    amax <- max(alltimes)
    amin <- 0
  }
  
  t0<-rep(s,N)
  t1<-Ltime
  t2<-Rtime
  t3<-abstime
  t4<-rep(horizon,N)
  ctime<-rep(NA,N)
  
  
  if (isIntervalCensored){
    
    ctime<-ifelse(t1>horizon,0,NA)
    ctime<-ifelse((t1>s) & (t2<=horizon) & (idm==1),1,ctime)
    ctime<-ifelse((t1<=horizon) & (t1>s) & (t3>horizon) & (t2>horizon) & (idm==1),2,ctime)
    ctime<-ifelse((t1<=horizon) & (t1>s) & (t3>horizon) & (idm==0),2,ctime)
    ctime<-ifelse((idm==1) & (t1<=s) & (t2>s) & (t2<=horizon),3,ctime)
    ctime<-ifelse((idm==0) & (t1>s) & (t1<=horizon) & (idd==1) & (t3<=horizon) & (t3>s),4,ctime)
    ctime<-ifelse((idd==1) & (t1<=s) & (idm==0) & (t3<=horizon) & (t3>s),5,ctime)
    ctime<-ifelse((t1<=s) & (t3>horizon) & (idm==0),6,ctime)
    ctime<-ifelse((t1<=s) & (t3>horizon) & (t2>horizon) & (idm==1),6,ctime)
    
    if(sum(is.na(ctime))>0){
      id<-which(is.na(ctime))
      id<-paste0(id,collapse=";")
      stop(paste0("For subject(s) : ",id," error of classification ctime, please check if the subjects are at risk between s and horizon"))
    }
    
    if(any(t3 <= s)) stop("Subjects dead or censored at or before the landmark time s must be excluded: ",
                          paste(which(t3 <= s), collapse=";"))
  }else{ #t1=t2
    
    stop("No implementation of BS when time to event are not interval-censored")
    
  }
  
  if (is.null(ncores)){
    ncores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
    if (is.na(ncores)) ncores <- max(1L, parallel::detectCores() - 1L)
  }
  ncores <- max(1L, min(as.integer(ncores), N))
  use_cl <- ncores > 1L
  
  ## ---- memoised predict(): each distinct (s,t) pair is computed once -------------
  ## ncores = 1 : cache kept in this function (lcache).
  ## ncores > 1 : one cache per worker, in the worker's global environment (.BScache),
  ##              so that it survives from the w1 loop to the w0 loop (parLapply gives
  ##              the same subjects to the same worker each time).
  ## Helps for the subject-independent calls (0,horizon), (0,s), (s,horizon) and lets
  ## w0 reuse what w1 already computed for the same subject.
  lcache <- new.env(hash = TRUE)
  
  ## P00(0,a) without calling predict() when a == 0
  P00_0 <- function(a, x){
    if(a <= 0) return(1)
    pt <- PR(0, a, x)$transprob
    pt[pt$Parameter == "p00", 2]
  }
  
  PR <- function(s_, t_, x){
    cache <- if (use_cl) get(".BScache", envir = .GlobalEnv) else lcache
    key <- sprintf("%d|%.15g|%.15g", x, s_, t_)
    v <- cache[[key]]
    if (is.null(v)){
      v <- predict(objectSurvival, s = s_, t = t_, nsim = 1,newdata = newdata[x, , drop = FALSE],conf.int = F)
      assign(key, v, envir = cache)
    }
    v
  }
  browser()
  ####################### calculate weights over subjects ######################

  if(inherits(objectSurvival,"idm")){
    
    w0w1func<-function(x){
      
      tryCatch({
        if(ctime[x]==0){return(list(w1=0,w0=1))}
        if(ctime[x]==1){return(list(w1=1,w0=0))}
        if(ctime[x]%in%c(2,6)){
          
          if(info.death==T){
            
            if(info.after.horizon==F){
              
              if(ctime[x]==2){ #check 9:00 09/10 no change for ctime = 2 !!
                predtiming<-PR(t1[x],horizon,x)
                p01<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                p00_1<-P00_0(t1[x],x)
                p01<-p00_1*p01
                
                predtiming<-PR(0,horizon,x)
                p00<-predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
                
                
                return(list(w1=(p01/(p00+p01)),
                            w0=(1-(p01/(p00+p01)))))#must be [0;1]
                
              }else{#check 9:00 09/10 changes for ctime = 6 !!
                
                predtiming<-PR(s,horizon,x)
                p01_1<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                predtiming<-PR(0,s,x)
                p01_1<-p01_1*predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
                
                
                predtiming<-PR(t1[x],s,x)
                p01<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                p00_1<-P00_0(t1[x],x)
                p01<-p00_1*p01
                predtiming<-PR(s,horizon,x)
                p01<-p01*predtiming$transprob[predtiming$transprob$Parameter=="p11",2]
                
                predtiming<-PR(0,horizon,x)
                p00<-predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
                
                
                
                return(list(w1=(p01_1/(p01_1+p00+p01)),
                            w0=(p00/(p01_1+p00+p01))))#must be : p01_1 < p01 and [0;1]
              }
            }else{
              
              if(idm[x]==1){
                
                
                
                if(ctime[x]==2){ # check 08:57 09/10 !!
                  
                  predtiming<-PR(t1[x],horizon,x)
                  p01<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                  p00_1<-P00_0(t1[x],x)
                  p01<-p00_1*p01
                  predtiming<-PR(horizon,t2[x],x)
                  p01<-p01*predtiming$transprob[predtiming$transprob$Parameter=="p11",2]
                  
                  
                  predtiming<-PR(horizon,t2[x],x)
                  p01_1<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                  predtiming<-PR(0,horizon,x)
                  p01_1<-p01_1*predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
                  
                  return(list(w1=(p01/(p01+p01_1)),
                              w0=(p01_1/(p01+p01_1)))) #must be [0;1]_
                  
                }else{ # check 09:00 09/10 !!
                  
                  predtiming<-PR(t1[x],s,x)
                  p01<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                  p00_1<-P00_0(t1[x],x)
                  p01<-p00_1*p01
                  predtiming<-PR(s,t2[x],x)
                  p01<-p01*predtiming$transprob[predtiming$transprob$Parameter=="p11",2]
                  
                  
                  predtiming<-PR(s,horizon,x)
                  p01_1<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                  predtiming<-PR(0,s,x)
                  p01_1<-p01_1*predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
                  predtiming<-PR(horizon,t2[x],x)
                  p01_1<-p01_1*predtiming$transprob[predtiming$transprob$Parameter=="p11",2]
                  
                  
                  p01_2<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                  predtiming<-PR(0,horizon,x)
                  p01_2<-p01_2*predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
                  
                  
                  
                  return(list(w1=(p01_1/(p01_2+p01_1+p01)),
                              w0=(p01_2/(p01_2+p01_1+p01))))#must be [0;1]
                }
                
              }else{
                if(idd[x]==1){
                  
                  if(ctime[x]==2){ #check 08:55 09/10 !!
                    
                    predtiming<-PR(t1[x],horizon,x)
                    p01<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                    p00_1<-P00_0(t1[x],x)
                    p01<-p00_1*p01
                    predtiming<-PR(horizon,t3[x],x)
                    p01<-p01*predtiming$transprob[predtiming$transprob$Parameter=="p11",2]*predtiming$intensity[3]
                    
                    p01_1<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]*predtiming$intensity[3]
                    predtiming<-PR(0,horizon,x)
                    p01_1<-p01_1*predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
                    predtiming<-PR(0,t3[x],x)
                   p02<-predtiming$transprob[predtiming$transprob$Parameter=="p00",2]*predtiming$intensity[2]
                                                      
                                                      
                   return(list(w1=(p01/(p01+p01_1+p02)),
                               w0=((p02+p01_1)/(p01+p01_1+p02))))
                                                      
                  }else{ #check 08:50 09/10 !!
                    
                    predtiming<-PR(s,horizon,x)
                    p01<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                    predtiming<-PR(0,s,x)
                    p01<-p01*predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
                    predtiming<-PR(horizon,t3[x],x)
                    p01<-p01*predtiming$transprob[predtiming$transprob$Parameter=="p11",2]*predtiming$intensity[3]
                    
                    p01_2<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]*predtiming$intensity[3]
                    predtiming<-PR(0,horizon,x)
                    p01_2<-p01_2*predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
                    
                    predtiming<-PR(t1[x],s,x)
                    p01_1<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                    p00_1<-P00_0(t1[x],x)
                    p01_1<-p00_1*p01_1
                    predtiming<-PR(s,t3[x],x)
                    p01_1<-p01_1*predtiming$transprob[predtiming$transprob$Parameter=="p11",2]*predtiming$intensity[3]
                    
                    predtiming<-PR(0,t3[x],x)
                    p02<-predtiming$transprob[predtiming$transprob$Parameter=="p00",2]*predtiming$intensity[2]
                    
                    return(list(w1=(p01/(p01+p01_1+p01_2+p02)),
                                w0=((p02+p01_2)/(p01+p01_1+p01_2+p02))))
                    
                  }
                  
                }else{
                  
                  if(ctime[x]==2){ #check 08:55 09/10
                    
                    predtiming<-PR(t1[x],horizon,x)
                    p01<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                    p00_1<-P00_0(t1[x],x)
                    p01<-p00_1*p01
                    predtiming<-PR(horizon,t3[x],x)
                    p01<-p01*predtiming$transprob[predtiming$transprob$Parameter=="p11",2]
                    
                    p01_1<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                    predtiming<-PR(0,horizon,x)
                    p01_1<-p01_1*predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
                                                      
                    predtiming<-PR(0,t3[x],x)
                    p02<-predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
                                                      
                                                      
                   return(list(w1=(p01/(p01+p01_1+p02)),
                               w0=((p02+p01_1)/(p01+p01_1+p02))))
                                                      
                  }else{ #check 08:50 09/10
                    
                    predtiming<-PR(s,horizon,x)
                    p01<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                    predtiming<-PR(0,s,x)
                    p01<-p01*predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
                    predtiming<-PR(horizon,t3[x],x)
                    p01<-p01*predtiming$transprob[predtiming$transprob$Parameter=="p11",2]
                    
                    p01_2<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                    predtiming<-PR(0,horizon,x)
                    p01_2<-p01_2*predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
                    
                    predtiming<-PR(t1[x],s,x)
                    p01_1<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
                    p00_1<-P00_0(t1[x],x)
                    p01_1<-p00_1*p01_1
                    predtiming<-PR(s,t3[x],x)
                    p01_1<-p01_1*predtiming$transprob[predtiming$transprob$Parameter=="p11",2]
                    
                    predtiming<-PR(0,t3[x],x)
                    p02<-predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
                    
                    return(list(w1=(p01/(p01+p01_1+p01_2+p02)),
                                w0=((p02+p01_2)/(p01+p01_1+p01_2+p02))))
                    
                  }
                }
              }
            }
            
          }else{
            
            predtiming<-PR(t1[x],horizon,x)
            F01<-predtiming$transprob[predtiming$transprob$Parameter=="F01",2]
            
            if(ctime[x]==2){
              return(list(w1=F01,
                          w0=(1-F01))) #must be [0;1]
            }else{
              
              predtiming<-PR(t1[x],s,x)
              F01_1<-predtiming$transprob[predtiming$transprob$Parameter=="F01",2]
              return(list(w1=(F01-F01_1),
                          w0=1-F01))#must be : p01_1 < p01 and [0;1]
            }
          }
        }
        
        if(ctime[x]%in%c(3)){
          
          
          if(info.death==T){ # check 09:08 09/10 !!
            
            predtiming<-PR(s,t2[x],x)
            p01_1<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
            predtiming<-PR(0,s,x)
            p01_1<-p01_1*predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
            
            predtiming<-PR(t1[x],s,x)
            p01_2<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
            p00_2<-P00_0(t1[x],x)
            p01_2<-p00_2*p01_2
            predtiming<-PR(s,t2[x],x)
            p01_2<-p01_2*predtiming$transprob[predtiming$transprob$Parameter=="p11",2]
            
            
            return(list(w1=(p01_1/(p01_2+p01_1)),
                        w0=0))
          }else{
            
            predtiming<-PR(t1[x],t2[x],x)
            F01<-predtiming$transprob[predtiming$transprob$Parameter=="F01",2]
            
            predtiming<-PR(t1[x],s,x)
            F01_1<-predtiming$transprob[predtiming$transprob$Parameter=="F01",2]
            
            return(list(w1=((F01-F01_1)/F01),
                        w0=0))
          }
          
        }
        
        if(ctime[x]%in%c(4,5)){
          
          if((ctime[x]==4) & (t3[x] <= t1[x])) return(list(w1=0,w0=1)) 
          
          if(info.death==T){
            
            if(ctime[x]==4){
              # check 09:08 09/10 no change for ctime = 4 !! 
              predtiming<-PR(t1[x],t3[x],x)
              p01<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
              p01<-p01*predtiming$intensity[3]
              p00<-P00_0(t1[x],x)
              p01<-p01*p00
              
              predtiming<-PR(0,t3[x],x)
              p02<-predtiming$transprob[predtiming$transprob$Parameter=="p00",2]*predtiming$intensity[2]
              
              
              return(list(w1=(p01/(p02+p01)),
                          w0=(1-(p01/(p02+p01)))))
            }else{ #check 09:13 09/10 changes made for ctime=5 !!
              
              predtiming<-PR(t1[x],s,x)
              p01<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
              p00<-P00_0(t1[x],x)
              p01<-p01*p00
              predtiming<-PR(s,t3[x],x)
              p01<-p01*predtiming$transprob[predtiming$transprob$Parameter=="p11",2]*predtiming$intensity[3]
              
              
              p01_1<-predtiming$transprob[predtiming$transprob$Parameter=="p01",2]
              p01_1<-p01_1*predtiming$intensity[3]
              predtiming<-PR(0,s,x)
              p01_1<-p01_1*predtiming$transprob[predtiming$transprob$Parameter=="p00",2]
              
              
              predtiming<-PR(0,t3[x],x)
              p02<-predtiming$transprob[predtiming$transprob$Parameter=="p00",2]*predtiming$intensity[2]
              
              return(list(w1=(p01_1/(p01_1+p01+p02)),
                          w0=(p02/(p01_1+p01+p02))))#must be : p01_1 < p01 and [0;1]
            }
            
          }else{
            
            predtiming<-PR(t1[x],t3[x],x)
            F01<-predtiming$transprob[predtiming$transprob$Parameter=="F01",2]
            
            if(ctime[x]==4){
              return(list(w1=F01,
                          w0=(1-F01)))
            }else{
              
              predtiming<-PR(t1[x],s,x)
              F01_1<-predtiming$transprob[predtiming$transprob$Parameter=="F01",2]
              return(list(w1=(F01-F01_1),
                          w0=(1-F01)))
            }
          }
          
        }
        
      }, error = function(e){
        message(sprintf("x=%d failed: %s", x, conditionMessage(e)))
        return(NA)
      })
      
    }

  
  }else{
    
    w0w1func<-function(x){
      
      tryCatch({
        if(ctime[x]==0){return(list(w1=0,w0=1))}
        if(ctime[x]==1){return(list(w1=1,w0=0))}
        if(ctime[x]%in%c(2,6)){
          
          if(info.after.horizon==F){
            
            if(ctime[x]==2){ #check 09:39 09/10 !!
              
            index_x<-which(objectSurvival$start==t1[x] & objectSurvival$end==horizon & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
            p01<-unique(objectSurvival$prediction_p01[index_x])[1]
            
            index_x<-which(objectSurvival$start==0 & objectSurvival$end==horizon & objectSurvival$id==ids_new[x] & objectSurvival$p00==T)
            p00<-unique(objectSurvival$prediction_p00[index_x])[1]
            
              return(list(w1=(p01/(p01+p00)),
                          w0=(1-(p01/(p01+p00))))) #must be [0;1]
              
            }else{
              
              index_x<-which(objectSurvival$start==s & objectSurvival$end==horizon & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
              p01_1<-unique(objectSurvival$prediction_p01[index_x])[1]
              
              
              index_x<-which(objectSurvival$start==t1[x] & objectSurvival$end==s & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
              p01<-unique(objectSurvival$prediction_p01[index_x])[1]
              index_x<-which(objectSurvival$start==s & objectSurvival$end==horizon & objectSurvival$id==ids_new[x] & objectSurvival$p11==T)
              p01<-p01*unique(objectSurvival$prediction_p11[index_x])[1]
              
              
              index_x<-which(objectSurvival$start==0 & objectSurvival$end==horizon & objectSurvival$id==ids_new[x] & objectSurvival$p00==T)
              p00<-unique(objectSurvival$prediction_p00[index_x])[1]
              
              return(list(w1=(p01_1/(p01+p00+p01_1)),
                          w0=(p00/(p01+p00+p01_1))))#must be : p01_1 < p01 and [0;1]
            }
            
            
            
          }else{
            
            if(idm[x]==1){ 
              
              if(ctime[x]==2){#check 09:39 09/10 !!
                
              index_x<-which(objectSurvival$start==horizon & objectSurvival$end==t2[x] & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
              p01<-unique(objectSurvival$prediction_p01[index_x])[1]
              
              
              index_x<-which(objectSurvival$start==t1[x] & objectSurvival$end==horizon & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
              p01_1<-unique(objectSurvival$prediction_p01[index_x])[1]
              index_x<-which(objectSurvival$start==horizon & objectSurvival$end==t2[x] & objectSurvival$id==ids_new[x] & objectSurvival$p11==T)
              p01_1<-p01_1*unique(objectSurvival$prediction_p11[index_x])[1]
                
                return(list(w1=(p01_1/(p01_1+p01)),
                            w0=(p01/(p01_1+p01))))
              }else{
                #check 09:39 09/10 !!
                
                index_x<-which(objectSurvival$start==horizon & objectSurvival$end==t2[x] & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
                p01<-unique(objectSurvival$prediction_p01[index_x])[1]
                
                
                index_x<-which(objectSurvival$start==t1[x] & objectSurvival$end==s & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
                p01_1<-unique(objectSurvival$prediction_p01[index_x])[1]
                index_x<-which(objectSurvival$start==s & objectSurvival$end==t2[x] & objectSurvival$id==ids_new[x] & objectSurvival$p11==T)
                p01_1<-p01_1*unique(objectSurvival$prediction_p11[index_x])[1]
                
                
                index_x<-which(objectSurvival$start==s & objectSurvival$end==horizon & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
                p01_2<-unique(objectSurvival$prediction_p01[index_x])[1]
                index_x<-which(objectSurvival$start==horizon & objectSurvival$end==t2[x] & objectSurvival$id==ids_new[x] & objectSurvival$p11==T)
                p01_2<-p01_2*unique(objectSurvival$prediction_p11[index_x])[1]
                
                return(list(w1=(p01_2/(p01_1+p01+p01_2)),
                            w0=(p01/(p01_1+p01+p01_2))))
              }
              
            }else{
              
              if(idd[x]==1){
                
                if(ctime[x]==2){ #check 09:54 09/10 !!
                
                  index_x<-which(objectSurvival$start==horizon & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p02_1==T)
                  p02_1<-unique(objectSurvival$prediction_p02_1[index_x])[1]
                  
                  index_x<-which(objectSurvival$start==0 & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p02_0==T)
                  p02_0<-unique(objectSurvival$prediction_p02_0[index_x])[1]
                  
               index_x<-which(objectSurvival$start==t1[x] & objectSurvival$end==horizon & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
                  p01<-unique(objectSurvival$prediction_p01[index_x])[1]
               index_x<-which(objectSurvival$start==horizon & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p12==T)
              p01<-p01*unique(objectSurvival$prediction_p12[index_x])[1]
                  
                  return(list(w1=(p01/(p01+p02_1+p02_0)),
                              w0=((p02_0+p02_1)/(p01+p02_1+p02_0))))
                  
                  
                }else{ #check 09:54 09/10 !!
                  
                  index_x<-which(objectSurvival$start==horizon & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p02_1==T)
                  p02_1<-unique(objectSurvival$prediction_p02_1[index_x])[1]
                  
                  index_x<-which(objectSurvival$start==0 & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p02_0==T)
                  p02_0<-unique(objectSurvival$prediction_p02_0[index_x])[1]
                  
                  index_x<-which(objectSurvival$start==s & objectSurvival$end==horizon & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
                  p01<-unique(objectSurvival$prediction_p01[index_x])[1]
                  index_x<-which(objectSurvival$start==horizon & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p12==T)
                  p01<-p01*unique(objectSurvival$prediction_p12[index_x])[1]
                  
                  
                  index_x<-which(objectSurvival$start==t1[x] & objectSurvival$end==s & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
                  p01_1<-unique(objectSurvival$prediction_p01[index_x])[1]
                  index_x<-which(objectSurvival$start==s & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p12==T)
                  p01_1<-p01_1*unique(objectSurvival$prediction_p12[index_x])[1]
                  
                  return(list(w1=(p01/(p01_1+p01+p02_0+p02_1)),
                              w0=((p02_0+p02_1)/(p01_1+p01+p02_0+p02_1))))
                }
                
                
              }else{
                
                if(ctime[x]==2){ #check 09:54 09/10 !!
                  
                  index_x<-which(objectSurvival$start==horizon & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
                  p01_1<-unique(objectSurvival$prediction_p01[index_x])[1]
                  
                  index_x<-which(objectSurvival$start==0 & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p00==T)
                  p00<-unique(objectSurvival$prediction_p00[index_x])[1]
                  
                  index_x<-which(objectSurvival$start==t1[x] & objectSurvival$end==horizon & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
                  p01<-unique(objectSurvival$prediction_p01[index_x])[1]
                  index_x<-which(objectSurvival$start==horizon & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p11==T)
                  p01<-p01*unique(objectSurvival$prediction_p11[index_x])[1]
                  
                  return(list(w1=(p01/(p01+p01_1+p00)),
                              w0=((p00+p01_1)/(p01+p01_1+p00))))
                  
                  
                }else{ #check 09:54 09/10 !!
                  
                  index_x<-which(objectSurvival$start==horizon & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
                  p01_2<-unique(objectSurvival$prediction_p01[index_x])[1]
                  
                  index_x<-which(objectSurvival$start==0 & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p00==T)
                  p00<-unique(objectSurvival$prediction_p00[index_x])[1]
                  
                  index_x<-which(objectSurvival$start==s & objectSurvival$end==horizon & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
                  p01<-unique(objectSurvival$prediction_p01[index_x])[1]
                  index_x<-which(objectSurvival$start==horizon & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p11==T)
                  p01<-p01*unique(objectSurvival$prediction_p11[index_x])[1]
                  
                  
                  index_x<-which(objectSurvival$start==t1[x] & objectSurvival$end==s & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
                  p01_1<-unique(objectSurvival$prediction_p01[index_x])[1]
                  index_x<-which(objectSurvival$start==s & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p11==T)
                  p01_1<-p01_1*unique(objectSurvival$prediction_p11[index_x])[1]
                  
                  return(list(w1=(p01/(p01_1+p01+p00+p01_2)),
                              w0=((p01_2+p00)/(p01_1+p01+p00+p01_2))))
                }
                
                
                
              }
            }
            
            
          }
        }
        
        if(ctime[x]%in%c(3)){ # check 10:06 09/10 !!
          
          
          index_x<-which(objectSurvival$start==s & objectSurvival$end==t2[x] & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
          p01_1<-unique(objectSurvival$prediction_p01[index_x])[1]
          
          
          index_x<-which(objectSurvival$start==t1[x] & objectSurvival$end==s & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
          p01_2<-unique(objectSurvival$prediction_p01[index_x])[1]
          index_x<-which(objectSurvival$start==s & objectSurvival$end==t2[x] & objectSurvival$id==ids_new[x] & objectSurvival$p11==T)
          p01_2<-p01_2*unique(objectSurvival$prediction_p11[index_x])[1]
          
          return(list(w1=(p01_1/(p01_2+p01_1)),
                      w0=0))
          
          
        }
        
        if(ctime[x]%in%c(4,5)){
          
          
          if((ctime[x]==4) & (t3[x] <= t1[x])) return(list(w1=0,w0=1)) 
          
          if(ctime[x]==4){ # check 09/10 at 10:07 !!
          index_x<-which(objectSurvival$start==t1[x] & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p02_1==T)
          p01<-unique(objectSurvival$prediction_p02_1[index_x])[1]
          
          index_x<-which(objectSurvival$start==0 & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p02_0==T)
          p02<-unique(objectSurvival$prediction_p02_0[index_x])[1]
          
          
          
            return(list(w1=(p01/(p02+p01)),
                        w0=(1-(p01/(p02+p01)))))
          }else{
            
            
            index_x<-which(objectSurvival$start==s & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p02_1==T)
            p02_1<-unique(objectSurvival$prediction_p02_1[index_x])[1]
            
            index_x<-which(objectSurvival$start==0 & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p02_0==T)
            p02<-unique(objectSurvival$prediction_p02_0[index_x])[1]
            
            
            index_x<-which(objectSurvival$start==t1[x] & objectSurvival$end==s & objectSurvival$id==ids_new[x] & objectSurvival$p01==T)
            p01_1<-unique(objectSurvival$prediction_p01[index_x])[1]
            index_x<-which(objectSurvival$start==s & objectSurvival$end==t3[x] & objectSurvival$id==ids_new[x] & objectSurvival$p12==T)
            p01_1<-p01_1*unique(objectSurvival$prediction_p12[index_x])[1]
            
            return(list(w1=(p02_1/(p02_1+p02+p01_1)),
                        w0=(p02/(p02_1+p02+p01_1))))#must be : p01_1 < p01 and [0;1]
          }
          
          
          
        }
        
      }, error = function(e){
        message(sprintf("x=%d failed: %s", x, conditionMessage(e)))
        return(NA)
      })
      
    }

  }
  

  if(ncores==1){
  if (isIntervalCensored){
    
    res <- lapply(seq_len(N), w0w1func)   # or parLapply
    w1 <- vapply(res, function(r) if (is.list(r)) r$w1 else NA_real_, numeric(1))
    w0 <- vapply(res, function(r) if (is.list(r)) r$w0 else NA_real_, numeric(1))

    
  }else{
    
    stop("No implementation of BS when time to event are not interval-censored")
    
    
  }
    
  }else{
    
    cl <- parallel::makeCluster(ncores, outfile="")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(HIDeM))
    parallel::clusterEvalQ(cl, .BScache <- new.env(hash = TRUE))
    parallel::clusterExport(cl,
                            c("objectSurvival","ctime","t1","t2","t3","s","horizon","info.death","info.after.horizon","newdata","idm","idd","ids_new"),
                            envir = environment())
    
    if (isIntervalCensored){
      
      res<-parallel::parLapply(cl, seq_len(N),w0w1func)
      w1 <- vapply(res, function(r) if (is.list(r)) r$w1 else NA_real_, numeric(1))
      w0 <- vapply(res, function(r) if (is.list(r)) r$w0 else NA_real_, numeric(1))
      
      
    }else{
      
      stop("No implementation of BS when time to event are not interval-censored")
      
      
    }
  }
  
  
  if(any(is.na(w1))|any(is.nan(w1))|any(is.na(w0))|any(is.nan(w0))){
    BS<-NULL
  }else{
    
    if(any(w1>1)){w1[which(w1>1)]<-1}
    if(any(w0>1)){w0[which(w0>1)]<-1}
    if(any(w0<0)){w0[which(w0<0)]<-0}
    BS<-sum(w1*(1-pred)^2+w0*(0-pred)^2)/N
  }
  end <- proc.time()
  res<-list(BS=BS,
            w1=w1,
            w0=w0,
            pred=pred,
            s=s,
            ctime=ctime,
            idd=idd,
            idm=idm,
            t1=t1,
            t2=t2,
            t3=t3,
            horizon=horizon,
            time=end-ptm)
  return(res)
  
}


