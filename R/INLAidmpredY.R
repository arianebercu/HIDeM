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
      
      # data structure
      ct <- INLAmodel$misc$configs$contents
      
      if(is.null(INLAmodel)){stop("The inla model for your marker could not be run, see above warnings.")}
      
      
      if(BLUP==F){
        
        #samples seed=seed cannot do parallel estimation on it 
        

        
        SMP <- inla.posterior.sample(Nsample, INLAmodel,seed=seed)
        
        make_dXINLA(formula=formLong[[indice]], X=attr(SMP, ".contents")$A, timevar=timevar, data=dataLongi_augmented,ct=ct,id=id) 
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

make_XINLA <- function(formula, X, timevar, data,use_splines = FALSE,ct,id, ...) {
  
  
  terms_labels <- attr(terms(formula), "term.labels")
  terms_fixed <- terms_labels[grepl(id, terms_labels)==F]
  
  terms_RE <- terms_labels[grepl(id, terms_labels)==T][1]
  terms_RE <- gsub("\\|.*", "", terms_RE)
  terms_RE <- as.formula(paste("~", terms_RE))
  terms_RE <- attr(terms(terms_RE), "term.labels")
  
  eval(dexpr, envir = data)
  browser()
  # Identify terms
  
  # Initialize derivative matrix
  X <- NULL
  
  for (lab in terms_fixed) {
    if (lab == timevar) {
      # simple linear time
      Xlab <- matrix(data[,colnames(data)%in%timevar],nrow=dim(data)[1],ncol=1)
      colnames(Xlab)<-paste0(lab,"_L1")
      
    } else if (grepl("\\(", lab) && grepl(timevar, lab)) {
      #INLA takes only function no : I(time^2)
      expr <- parse(text = lab)[[1]]
      Xlab <- matrix(eval(expr, envir = data),ncol=1,nrow=dim(data)[1])
      colnames(dXlab)<-paste0(gsub("[())]","",lab),"_L1")
      
    } else if (use_splines && grepl("bs\\(|ns\\(", lab)) {
      
      # Extract knots etc. from the original call if needed
      stop("Spline derivative not supported so far")
      
    } 
    X<-cbind(X,Xlab)
  }
  
  XRE <- NULL
  k<-1
  for (lab in terms_RE) {
    if (lab == timevar) {
      # simple linear time
      Xlab <- matrix(data[,colnames(data)%in%timevar],nrow=dim(data)[1],ncol=1)
      nameX<-paste0("ID",lab,"_L1")
      
    } else if (grepl("\\(", lab) && grepl(timevar, lab)) {
      #INLA takes only function no : I(time^2)
      
      expr <- parse(text = lab)[[1]]
      Xlab <- matrix(eval(expr, envir = data),ncol=1,nrow=dim(data)[1])
      nameX<-paste0("ID",gsub("[())]","",lab),"_L1")
      
    } else if (use_splines && grepl("bs\\(|ns\\(", lab)) {
      
      # Extract knots etc. from the original call if needed
      stop("Spline derivative not supported so far")
      
    } 
    XRE<-cbind(XRE,Xlab)
    colnames(XRE)[k]<-nameX
    k<-k+1
  }
  
  return(list(X=X,XRE=XRE))
}

make_dXINLA <- function(formula, X, timevar, data,use_splines = FALSE,ct,id, ...) {
  
  
  terms_labels <- attr(terms(formula), "term.labels")
  terms_fixed <- terms_labels[grepl(id, terms_labels)==F]
  
  terms_RE <- terms_labels[grepl(id, terms_labels)==T][1]
  terms_RE <- gsub("\\|.*", "", terms_RE)
  terms_RE <- as.formula(paste("~", terms_RE))
  terms_RE <- attr(terms(terms_RE), "term.labels")
  
  browser()
  # Identify terms
  
  # Initialize derivative matrix
  dX <- X * 0
  
  for (lab in terms_fixed) {
    if (lab == timevar) {
      # simple linear time
      dXlab <- matrix(1,nrow=dim(data)[1],ncol=1)
      colnames(dXlab)<-paste0(lab,"_L1")
      
    } else if (grepl("\\(", lab) && grepl(timevar, lab)) {
      #INLA takes only function no : I(time^2)
      expr <- parse(text = lab)[[1]]
      dexpr <- Deriv::Deriv(expr, timevar)
      dXlab <- matrix(eval(dexpr, envir = data),ncol=1,nrow=dim(data)[1])
      colnames(dXlab)<-paste0(gsub("[())]","",lab),"_L1")
      
    } else if (use_splines && grepl("bs\\(|ns\\(", lab)) {
      
      # Extract knots etc. from the original call if needed
      stop("Spline derivative not supported so far")
      
    } 
    dX<-cbind(dX,dXlab)
  }
  
  dXRE <- list()
  length(dXRE)<-length(terms_RE)
  k<-1
  idd<-unique(data[,colnames(data)%in% id])
  for (lab in terms_RE) {
    if (lab == timevar) {
      # simple linear time
      dXlab <- do.call(cbind,lapply(idd,FUN=function(x){
        vec<-rep(0,dim(data)[1])
        vec[which(data[,colnames(data)%in%id]==x)]<-1
        return(vec)
      }))
      nameX<-paste0("ID",lab,"_L1")
      
    } else if (grepl("\\(", lab) && grepl(timevar, lab)) {
      #INLA takes only function no : I(time^2)
      
      expr <- parse(text = lab)[[1]]
      dexpr <- Deriv::Deriv(expr, timevar)
      xx<-eval(dexpr, envir = data)
      dXlab <- do.call(cbind,lapply(idd,FUN=function(x){
        vec<-rep(0,dim(data)[1])
        vec[which(data[,colnames(data)%in%id]==x)]<-xx[which(data[,colnames(data)%in%id]==x)]
        return(vec)
      }))
      nameX<-paste0("ID",gsub("[())]","",lab),"_L1")
      
    } else if (use_splines && grepl("bs\\(|ns\\(", lab)) {
      
      # Extract knots etc. from the original call if needed
      stop("Spline derivative not supported so far")
      
    } 
    dXRE[[k]]<-dXlab
    names(dXRE)[k]<-nameX
    k<-k+1
  }
  
  return(list(dX=dX,dXRE=dXRE))
}
