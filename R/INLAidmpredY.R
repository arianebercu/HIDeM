### Code:
##' @title Calculate predictions for time-depend covariates using INLA
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  


INLAidmpredY<-function(timeVar,truncated,formLong,dataSurv,dataLongi,id,
                  Nsample,t0,t1,t2,t3,
                  ctime,modelY,seed,BLUP){
  
  # define timePoints of prediction : 

  set.seed(seed)
  
  idsubjects<-unique(dataSurv[,colnames(dataSurv)%in%id])
  

  
  timePointsdata<-do.call(rbind, lapply(idsubjects,FUN=function(x){
    index<-which(dataSurv[,colnames(dataSurv)%in%id]==x)
    max.int<-ifelse(ctime[index]%in%c(2:5),t2[index],t3[index])
    min.int<-t1[index]
  
    
    timePoints<-gauss_kronrod_points(lower.intdouble=min.int,
                                     upper.intdouble=max.int,
                                     end.time=t3[index],
                                     truncated=truncated,
                                     entry.time=t0[index])
    return(data.frame(index=x,timePoints=timePoints))}))
  
  
  

  
  ## augmentation of the data 
  colnames(timePointsdata)<-c(id,timeVar)
  dataLongi_augmented<-merge(timePointsdata,dataLongi,by=c(id,timeVar),all.x=T,all.y=T)
  rownames(dataLongi_augmented)<-NULL
  dataLongi_augmented<-dataLongi_augmented[order(dataLongi_augmented[,colnames(dataLongi_augmented)%in%id],
                                                 dataLongi_augmented[,colnames(dataLongi_augmented)%in%timeVar]),]
  
  Yall<-list()
  length(Yall)<-length(formLong)
  dataSurv<-dataSurv[order(dataSurv[,colnames(dataSurv)%in%id]),]
  
  

    for(indice in 1:length(formLong)){


      INLAmodel<-modelY$modelY[[indice]]
      
      if(is.null(INLAmodel)){stop("The inla model for your marker could not be run, see above warnings.")}
      
      
      if(BLUP==F){
        
        #samples seed=seed cannot do parallel estimation on it 
        SMP <- inla.posterior.sample(Nsample, INLAmodel,seed=seed)
        linPred <- sapply(SMP, function(x) x$latent) 

        # keep only indice we want: 
        # Collapse each row into a string
        key1 <- do.call(paste, c(dataLongi_augmented[,colnames(dataLongi_augmented)%in%c(id,timeVar)], sep = "\r"))
        key2 <- do.call(paste, c(timePointsdata, sep = "\r"))
        
        # Find indices of rows from data1 that are in data2
        # while keeping order of key2
        indices <- match(key2, key1)
        PredYx<-linPred[indices,]
        
        #add BLUP first column
        PredYx<-cbind(INLAmodel$summary.linear.predictor$mean[indices],PredYx)
        #add informations 
        Outcome<-all.vars(terms(formLong[[indice]]))[1]
        PredYx<-cbind(timePointsdata,Outcome,PredYx)
        colnames(PredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
        
        
        
        
        
      }else{
        
        # keep only indice we want: 
        # Collapse each row into a string
        key1 <- do.call(paste, c(dataLongi_augmented[,colnames(dataLongi_augmented)%in%c(id,timeVar)], sep = "\r"))
        key2 <- do.call(paste, c(timePointsdata, sep = "\r"))
        
        # Find indices of rows from data1 that are in data2
        # while keeping order of key2
        indices <- match(key2, key1)
        #add BLUP first column
        #add informations 
        Outcome<-all.vars(terms(formLong[[indice]]))[1]
        PredYx<-cbind(timePointsdata,Outcome,INLAmodel$summary.linear.predictor$mean[indices])
        colnames(PredYx)[4]<-paste0("Sample_1")
        
      }

      Yall[[indice]]<-PredYx
      
  
      }

  Yall<-do.call(rbind,Yall)
  
  return(Yall[,colnames(Yall)%in%c(id,timeVar,"Outcome",paste0("Sample_",c(1:Nsample)))])
  
}


make_dXINLA <- function(formula, X, timevar, data,use_splines = FALSE, ...) {
  
  
  # Identify terms
  terms_labels <- attr(terms(formula), "term.labels")
  
  # Initialize derivative matrix
  dX <- X * 0
  
  for (lab in terms_labels) {
    if (lab == timevar) {
      # simple linear time
      dX[, lab] <- 1
      
    } else if (grepl("^I\\(", lab)) {
      # e.g. I(time^2)
      expr <- parse(text = gsub("I\\((.*)\\)", "\\1", lab))[[1]]
      dexpr <- D(expr, timevar)
      dX[, lab] <- eval(dexpr,envir = data)
      
    } else if (use_splines && grepl("bs\\(|ns\\(", lab)) {
      
      # Extract knots etc. from the original call if needed
      stop("Spline derivative not supported so far")
      
    } else {
      # other covariates (not functions of time) → derivative is 0
      dX[, lab] <- 0
    }
  }
  
  # Intercept derivative is always 0
  if ("(Intercept)" %in% colnames(X)) {
    dX[, "(Intercept)"] <- 0
  }
  
  dX
}
