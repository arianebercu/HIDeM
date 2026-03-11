### Code:
##' @title Calculate predictions for time-depend covariates using INLA
##' @param objectY  A predYidm object from HIDeM package containing the univariate joint model with competing risk from step 1
##' @param objectSurvival A DYNidm object from HIDeM package containing the illness-death model estimation with time-dependent covariates 
##' @param newdata The newdata on which we want to do survival predictions 
##' @param s entry time 
##' @param horizon horizon time 
##' @param envir working environment 
#' @importFrom foreach "%do%"
#' @importFrom foreach "%dopar%"
#' @importFrom Deriv "Deriv"
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  
#' @export


DYNpredIDM<-function(objectY,
                          objectSurvival,
                          newdata,s,horizon,
                          envir=parent.frame()){

  call <- match.call()
  ptm <- proc.time()

  if(missing(objectY)){stop("Need to specify objectY as a predYidm object")}
  if(!inherits(objectY,"predYidm")){stop("Need to specify objectY as a predYidm object")}
  if(missing(objectSurvival)){stop("Need to specify objectSurvival as a regDYNidm or DYNidm object")}
  if(!inherits(objectSurvival,"DYNidm")){stop("Need to specify objectSurvival as a DYNidm object")}

  
  if(missing(newdata)){stop("Need to provide newData as a data.frame")}
  if(!inherits(newdata,"data.frame")){stop("Need to provide newData as a data.frame")}
  
  if(sum(is.na(newdata))>0)stop("Need a new data frame with no missing data.")
  if(missing(s)|missing(horizon))stop("Need to specify s and horizon")
  if(!inherits(s,c("numeric","integer")))stop("s need to be an integer or numeric")
  if(!inherits(horizon,c("numeric","integer")))stop("horizon need to be an integer or numeric")
  if(length(s)!=1)stop("Length of s need to be 1")
  if(length(horizon)!=1)stop("Length of horizon need to be 1")
  if((s < 0) | (horizon < 0) | (s >= horizon))stop("s and horizon need to be numeric superior or equal to 0 with s < horizon")
  
  timeVar<-objectY$timeVar
  id<-objectY$id
  
  variables<-c(id,
               timeVar,
               objectSurvival$Xnames01,
               objectSurvival$Xnames02,
               objectSurvival$Xnames12,
               objectSurvival$timedepXnames01,
               objectSurvival$timedepXnames02,
               objectSurvival$timedepXnames12)
  variables<-unique(variables)
  
  if(any((variables%in%colnames(newdata))==F)){stop(paste0("Need a new data frame with no missing data on the variables : ",paste0(variables,collapse = ", ")))}
  
  if(objectY$method!="INLA"){stop("Prediction not available with JMBayes 2 method")}
  
  ############################## supress info before s #########################
  
  newdata<-newdata[newdata[,colnames(newdata)%in%timeVar]<=s,]
  newdata<-newdata[,colnames(newdata)%in%variables]
  newdata<-na.omit(newdata)
  if(dim(newdata)[1]==0){stop("No follow-up after s, should proceed with smaller values for s")}
  N<-length(unique(newdata[,colnames(newdata)%in%objectY$id]))
  
  
  ######################## prepare for survival ################################
  
  ######################### association check RE ###############################
  assoc<-objectSurvival$assoc
  verifRE<-0
  for(k in 1:length(assoc)){
      vv<-na.omit(unlist(assoc[[k]]))
      if(sum(vv%in%"RE")==length(vv)){verifRE<-verifRE+1}
    }
    
  ####################### time fixed var #######################################
  x01<-newdata[,colnames(newdata)%in%c(objectSurvival$Xnames01,id)]
  x01<-unique(newdata)
  x01<-x01[,colnames(x01)%in%c(objectSurvival$Xnames01)]
  NC01<-NCOL(x01)
  
  x02<-newdata[,colnames(newdata)%in%c(objectSurvival$Xnames02,id)]
  x02<-unique(newdata)
  x02<-x02[,colnames(x02)%in%c(objectSurvival$Xnames02)]
  NC02<-NCOL(x02)
  
  size1 <- NC01+NC02
  
  noVar<-c(ifelse(as.integer(NC01)>0,0,1),
           ifelse(as.integer(NC02)>0,0,1))
  
  nvat01 <- ifelse(noVar[1]==1,0,NC01)
  nvat02 <- ifelse(noVar[2]==1,0,NC02)
  nvat12<-length(unique(objectSurvival$Xnames12))
  
  dimnva01<-ifelse(nvat01==0,1,nvat01)
  dimnva02<-ifelse(nvat02==0,1,nvat02)
  
  NC<-c(NC01,NC02)
  
  if(noVar[1]==1){ve01<-as.double(rep(0,N))}else{ve01<-as.double(x01)}
  if(noVar[2]==1){ve02<-as.double(rep(0,N))}else{ve02<-as.double(x02)}
  
  t0<-rep(s,N)
  t1<-rep(horizon,N)
  
  ########################## time dep var ######################################
  outcome01<-objectSurvival$linktimedepXnames01
  outcome02<-objectSurvival$linktimedepXnames02
  
  p01<-length(outcome01)
  dimp01<-ifelse(length(outcome01)>0,length(outcome01),1)
  p02<-length(outcome02)
  dimp02<-ifelse(length(outcome02)>0,length(outcome02),1)
  if(p01==0 & p02==0){stop("No time dependent variable please refer to HIDeM::intensity for calculation of intensities with only time fixed covariates")}
  
  if(objectSurvival$method=="splines"){
    
   nknots01 <- length(unique(objectSurvival$knots01))
   nknots02<- length(unique(objectSurvival$knots02))
   amin<-min(newdata[,colnames(newdata)%in%timeVar])
   amax<-max(newdata[,colnames(newdata)%in%timeVar])
   
   if(knots01[1]> amin) stop(paste("Transition 0->1: Smallest time point should not appear before the starting knots at :",knots01[1]))
   if (knots01[length(knots01)] < amax) stop(paste("Transition 0->1: Highest time point should not appear after the last knots at :",knots01[length(knots01)]))
   
   if(knots02[1]> amin) stop(paste("Transition 0->2: Smallest time point should not appear before the starting knots at :",knots02[1]))
   if (knots02[length(knots02)] < amax) stop(paste("Transition 0->2: Highest time point should not appear after the last knots at :",knots02[length(knots02)]))
   
      
    size_spline<-nknots01+nknots02+4
    index_size1<-c(1:(nknots01+nknots02+4))
    start<-nknots01+nknots02+nknots12+6
    if(nvat01>0|nvat02>0){
    index_size1<-c(index_size1,((start+1):(nvat01+nvat02+start)))
    }
    start<-nvat01+nvat02+nvat12+start

  }else{ 
    index_size1<-c(1:4)
    
    if(nvat01>0|nvat02>0){
    index_size1<-c(index_size1,(7:(nvat01+nvat02+6)))
    }
    start<-nvat01+nvat02+nvat12+6
    size_spline<-4 }
  

  size1<-size1+size_spline
  
  ### index to keep with binit
  # jump x12 and end after p02
  index<-c(index_size1,(start+1):(start+p01+p02))
  
  size_V<-length(index)
  ############################################################################
  #################### defines initiate values ###############################
  ############################################################################
  
  # if length of regDYNidm or DYNidm superior to 1 need to summaries over replicates 
  if(length(objectSurvival$DYNidm)!=1){
    istop<-lapply(objectSurvival$DYNidm,FUN = function(x){
      if(x$istop%in%c(1,3)){return(T)}else{return(F)}
    })
    istop<-do.call(c,istop)
    if(sum(istop==F)==length(istop)){stop("All the survival models did not converged")}
    binit<-prepareData(object=objectSurvival$DYNidm,istop=istop,index=index)
  }else{
    if(!objectSurvival$DYNidm[[1]]$istop%in%c(1,3)){stop("The survival model did not converged")}
    
    binit<-objectSurvival$DYNidm[[1]]$b[index]
  }

  ############################ keep index of variables selected ################
  
  varY<-unique(c(objectSurvival$timedepXnames01,
           objectSurvival$timedepXnames02,
           objectSurvival$timedepXnames12))
  Yindex<-lapply(objectY$formLong,FUN=function(x){as.character(x[[2]])})
  Yindex<-do.call(c,Yindex)
  index<-which(Yindex%in%varY)
  
  ############################ perform prediction of Y #########################
  
  dataY<-DYNINLAidmpredY(object=objectY$modelY,
                         newdata=newdata,
                         s=s,
                         horizon=horizon,
                         scale.X=objectSurvival$scale.X,
                         assoc=objectY$assoc,
                         assocSurv=objectSurvival$assoc,
                         id=id,
                         timeVar=timeVar,
                         formLong=objectY$formLong,
                         formSurv=objectY$formSurv,
                         basRisk=objectY$basRisk,
                         index=index,
                         family=objectY$family,
                         envir=envir)
  
  
  ########################## check prediction ##################################
  NtimePoints<-255
  for( m in unique(c(outcome01,outcome02))){
    subdata<-dataY[dataY$Outcome==m,]
    x<-table(subdata[,colnames(subdata)%in%id])
    if(any(x!=NtimePoints)){stop("Prediction of marker ",m," could not be perform for each quadrature points")}
    
  }
  
  dataY$Outcome<-as.character(dataY$Outcome)
  # attention if NtimePoints equidistant with INLA then NtimePoints takes 
  # need ID to be numeric -- then 
  dataY[,colnames(dataY)%in%id]<-as.numeric(dataY[,colnames(dataY)%in%id])
  # to keep tracks of time order for each individual 
  dataY$order<-as.numeric(ave(dataY[,colnames(dataY)%in%id], cumsum(c(TRUE, diff(dataY[,colnames(dataY)%in%id]) != 0)), FUN = seq_along))
  
  
  if(length(outcome01)>=1){
    y01<-dataY[dataY$Outcome%in%outcome01,]
    # order  by individual and timeline 
    y01<-y01[order(y01[,colnames(y01)%in%id],y01$order),4]
    
    
  }else{
    y01<-rep(0,N*NtimePoints)
  }
  
  if(length(outcome02)>=1){
    y02<-dataY[dataY$Outcome%in%outcome02,]
    # order  by individual and timeline 
    y02<-y02[order(y02[,colnames(y02)%in%id],y02$order),4]
    
  }else{
    y02<-rep(0,N*NtimePoints)
  }
  
 
  if(objectSurvival$method=="splines"){
    res<-rep(0,N)
    out<- tryCatch({  .Fortran("citimedep",
                             ## input
                             as.double(binit),
                             as.integer(sizeV),
                             as.double(knots01),
                             as.double(knots02),
                             as.integer(N),
                             as.integer(nknots01),
                             as.integer(nknots02),
                             as.double(ve01),
                             as.double(ve02),
                             as.double(y01),
                             as.double(y02),
                             as.integer(p01),
                             as.integer(p02),
                             as.integer(dimp01),
                             as.integer(dimp02),
                             as.integer(NtimePoints),
                             as.integer(dimnva01),
                             as.integer(dimnva12),
                             as.integer(nvat01),
                             as.integer(nvat02),
                             as.double(t0),
                             as.double(t1),
                             likelihood_res=as.double(res),
                             PACKAGE="HIDeM")$likelihood_res
  }, error = function(e) {
    # Return NULL on error to skip this patient
    NULL
  })
  }else{
    res<-rep(0,N)
    out<- tryCatch({   .Fortran("ciweibtimedep",
                                ## input
                                as.double(binit),
                                as.integer(size_V),
                                as.integer(N),
                                as.double(ve01),
                                as.double(ve02),
                                as.double(y01),
                                as.double(y02),
                                as.integer(p01),
                                as.integer(p02),
                                as.integer(dimp01),
                                as.integer(dimp02),
                                as.integer(NtimePoints),
                                as.integer(dimnva01),
                                as.integer(dimnva02),
                                as.integer(nvat01),
                                as.integer(nvat02),
                                as.double(t0),
                                as.double(t1),
                                likelihood_res=as.double(res),
                                PACKAGE="HIDeM")$likelihood_res
    }, error = function(e) {
      # Return NULL on error to skip this patient
      NULL
    })
  }
  
  return(out)
 
}


prepareData<-function(object,istop,index){
  
  Nrep<-sum(istop==T)
  b<-lapply(object,FUN=function(x){
    if(x$istop%in%c(1,3)){return(x$b[index])}else{return(rep(NA,length(x$b[index])))}})
  b<-do.call(rbind,b)
  b<-na.omit(b)
  b<-colSums(b)/Nrep
  return(b)
}
