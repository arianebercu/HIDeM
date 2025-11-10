### Code:
##' @title Calculate predictions for time-depend covariates using INLA
#' @importFrom foreach "%do%"
#' @importFrom foreach "%dopar%"
#' @importFrom Deriv "Deriv"
#' @importFrom INLA "inla.posterior.sample"
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  


INLAidmpredY<-function(timeVar,truncated,formLong,dataSurv,dataLongi,id,
                  Nsample,t0,t1,t2,t3,assoc,
                  ctime,modelY,seed,BLUP,nproc,clustertype){

  
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
  
  NtimePoints<-ifelse(truncated==F,256,271)
  
  if(nproc==1){
    
    
    for(indice in 1:length(formLong)){


      browser()
      INLAmodel<-modelY$modelY[[indice]]
      
      # data structure
      ct <- INLAmodel$misc$configs$contents
      
      if(is.null(INLAmodel)){stop("The inla model for your marker could not be run, see above warnings.")}
      choiceY<-na.omit(unlist(assoc[[indice]]))
      
      
      if(BLUP==F){
        
        
        #samples seed=seed cannot do parallel estimation on it 
        SMP <- INLA::inla.posterior.sample(Nsample, INLAmodel,seed=seed)
        
        res<-NULL
        key1 <- do.call(paste, c(dataLongi_augmented[,colnames(dataLongi_augmented)%in%c(id,timeVar)], sep = "\r"))
        key2 <- do.call(paste, c(timePointsdata, sep = "\r"))
        # keep only indice we want: 
        # Collapse each row into a string
        
        # Find indices of rows from data1 that are in data2
        # while keeping order of key2
        indices <- match(key2, key1)
        
        if("value"%in% choiceY){
        
        Y<-do.call(cbind,
                     lapply(c(1:Nsample),FUN=function(x){make_XINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=SMP[[x]])}))
        Y<-Y[indices,]
        Outcome<-all.vars(terms(formLong[[indice]]))[1]
        PredYx<-cbind(timePointsdata,Outcome=Outcome,Y)
        colnames(PredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
        res<-rbind(res,PredYx)
        }
        
        if("slope"%in% choiceY){
        dY<-do.call(cbind,
                    lapply(c(1:Nsample),FUN=function(x){make_dXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=SMP[[x]])}))
        dY<-dY[indices,]
        slopeOutcome<-paste0("slope_",Outcome)
        slopePredYx<-cbind(timePointsdata,Outcome=slopeOutcome,dY)
        colnames(slopePredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
        res<-rbind(res,slopePredYx)
        }
        
        if("RE"%in% choiceY){
          
        REY<-do.call(cbind,
                    lapply(c(1:Nsample),FUN=function(x){make_REXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=SMP[[x]])}))
       REY<-REY[indices,]
       REOutcome<-paste0("RE_",Outcome)
       REPredYx<-cbind(timePointsdata,Outcome=REOutcome,REY)
       colnames(REPredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
       res<-rbind(res,REPredYx)
        }
        
        Yall[[indice]]<- res
        
        
        
      }else{
        
        Nsample<-1
        
        res<-NULL
        key1 <- do.call(paste, c(dataLongi_augmented[,colnames(dataLongi_augmented)%in%c(id,timeVar)], sep = "\r"))
        key2 <- do.call(paste, c(timePointsdata, sep = "\r"))
        # keep only indice we want: 
        # Collapse each row into a string
        
        # Find indices of rows from data1 that are in data2
        # while keeping order of key2
        indices <- match(key2, key1)
        
        if("value"%in% choiceY){
          Y<-as.matrix(make_XINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=INLAmodel))
          Y<-Y[indices,]
          Outcome<-all.vars(terms(formLong[[indice]]))[1]
          PredYx<-cbind(timePointsdata,Outcome=Outcome,Y)
          colnames(PredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
          res<-rbind(res,PredYx)
        }
        
        if("slope"%in% choiceY){
        dY<-as.matrix(make_dXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=INLAmodel))
        dY<-dY[indices,]
        slopeOutcome<-paste0("slope_",Outcome)
        slopePredYx<-cbind(timePointsdata,Outcome=slopeOutcome,dY)
        colnames(slopePredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
        res<-rbind(res,slopePredYx)
        }
        
        if("RE"%in% choiceY){
        REY<-as.matrix(make_REXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=INLAmodel))
        REY<-REY[indices,]
        REOutcome<-paste0("RE_",Outcome)
        REPredYx<-cbind(timePointsdata,Outcome=REOutcome,REY)
        colnames(REPredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
        res<-rbind(res,REPredYx)
        }
        
        Yall[[indice]]<- res
        
        
        
      }

      
  
      }

  Yall<-do.call(rbind,Yall)
  }else{
    
  
      if(is.null(clustertype)){
        clustpar <- parallel::makeCluster(nproc)#, outfile="")
      }
      else{
        clustpar <- parallel::makeCluster(nproc, type=clustertype)#, outfile="")
      }
      
      doParallel::registerDoParallel(clustpar)
      
      
      Yall<-foreach::foreach(indice=4:length(formLong),
                       .combine = 'list',
                       packages=c("INLA","Deriv","HIDeM"))%dopar%{
                         
                         
                         INLAmodel<-modelY$modelY[[indice]]
                         
                         # data structure
                         ct <- INLAmodel$misc$configs$contents
                         
                         choiceY<-na.omit(unlist(assoc[[indice]]))
                         
                         if(is.null(INLAmodel)){stop("The inla model for your marker could not be run, see above warnings.")}
                         
                         
                         if(BLUP==F){
                           
                           
                           #samples seed=seed cannot do parallel estimation on it 
                           SMP <- INLA::inla.posterior.sample(Nsample, INLAmodel,seed=seed)
                           
                           res<-NULL
                           key1 <- do.call(paste, c(dataLongi_augmented[,colnames(dataLongi_augmented)%in%c(id,timeVar)], sep = "\r"))
                           key2 <- do.call(paste, c(timePointsdata, sep = "\r"))
                           # keep only indice we want: 
                           # Collapse each row into a string
                           
                           # Find indices of rows from data1 that are in data2
                           # while keeping order of key2
                           indices <- match(key2, key1)
                           
                           if("value"%in% choiceY){
                             
                             Y<-do.call(cbind,
                                        lapply(c(1:Nsample),FUN=function(x){make_XINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=SMP[[x]])}))
                             Y<-Y[indices,]
                             Outcome<-all.vars(terms(formLong[[indice]]))[1]
                             PredYx<-cbind(timePointsdata,Outcome=Outcome,Y)
                             colnames(PredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
                             res<-rbind(res,PredYx)
                           }
                           
                           if("slope"%in% choiceY){
                             dY<-do.call(cbind,
                                         lapply(c(1:Nsample),FUN=function(x){make_dXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=SMP[[x]])}))
                             dY<-dY[indices,]
                             slopeOutcome<-paste0("slope_",Outcome)
                             slopePredYx<-cbind(timePointsdata,Outcome=slopeOutcome,dY)
                             colnames(slopePredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
                             res<-rbind(res,slopePredYx)
                           }
                           
                           if("RE"%in% choiceY){
                             
                             REY<-do.call(cbind,
                                          lapply(c(1:Nsample),FUN=function(x){make_REXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=SMP[[x]])}))
                             REY<-REY[indices,]
                             REOutcome<-paste0("RE_",Outcome)
                             REPredYx<-cbind(timePointsdata,Outcome=REOutcome,REY)
                             colnames(REPredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
                             res<-rbind(res,REPredYx)
                           }
                           
                           return(res)
                           
                           
                           
                         }else{
                           
                           Nsample<-1
                           
                           res<-NULL
                           key1 <- do.call(paste, c(dataLongi_augmented[,colnames(dataLongi_augmented)%in%c(id,timeVar)], sep = "\r"))
                           key2 <- do.call(paste, c(timePointsdata, sep = "\r"))
                           # keep only indice we want: 
                           # Collapse each row into a string
                           
                           # Find indices of rows from data1 that are in data2
                           # while keeping order of key2
                           indices <- match(key2, key1)
                           
                           if("value"%in% choiceY){
                             Y<-as.matrix(make_XINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=INLAmodel))
                             Y<-Y[indices,]
                             Outcome<-all.vars(terms(formLong[[indice]]))[1]
                             PredYx<-cbind(timePointsdata,Outcome=Outcome,Y)
                             colnames(PredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
                             res<-rbind(res,PredYx)
                           }
                           
                           if("slope"%in% choiceY){
                             dY<-as.matrix(make_dXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=INLAmodel))
                             dY<-dY[indices,]
                             slopeOutcome<-paste0("slope_",Outcome)
                             slopePredYx<-cbind(timePointsdata,Outcome=slopeOutcome,dY)
                             colnames(slopePredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
                             res<-rbind(res,slopePredYx)
                           }
                           
                           if("RE"%in% choiceY){
                             REY<-as.matrix(make_REXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=INLAmodel))
                             REY<-REY[indices,]
                             REOutcome<-paste0("RE_",Outcome)
                             REPredYx<-cbind(timePointsdata,Outcome=REOutcome,REY)
                             colnames(REPredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
                             res<-rbind(res,REPredYx)
                           }
                           
                           return(res)
                           
                           
                           
                         }
                           
                         
                         
                       }
      
      
      parallel::stopCluster(clustpar)
  }
  
  return(Yall[,colnames(Yall)%in%c(id,timeVar,"Outcome",paste0("Sample_",c(1:Nsample)))])
  

  
}

make_XINLA <- function(formula, timeVar, data,use_splines = FALSE,ct,id,SMP, ...) {
  
  
  terms_labels <- attr(terms(formula), "term.labels")
  terms_fixed <- terms_labels[grepl(id, terms_labels)==F]
  
  terms_RE <- terms_labels[grepl(id, terms_labels)==T][1]
  terms_RE <- gsub("\\|.*", "", terms_RE)
  terms_RE <- as.formula(paste("~", terms_RE))
  terms_RE <- attr(terms(terms_RE), "term.labels")

  
  # Identify terms
 if(paste0(id,"Intercept_L1")%in%ct$tag){
   terms_RE<-c("Intercept",terms_RE)
 }
  if("Intercept_L1"%in%ct$tag){
    terms_fixed<-c("Intercept",terms_fixed)
    
  }
  
  X<-NULL
  B<-matrix(0,nrow=dim(data)[1],ncol=length(terms_fixed))
  k<-1
  for (lab in terms_fixed) {
    
    if(lab=="Intercept"){
      
      Xlab <- matrix(1,nrow=dim(data)[1],ncol=1)
      colnames(Xlab)<-paste0(lab,"_L1")
      
      start<-ct$start[which(ct$tag==paste0(lab,"_L1"))]
      B[,k]<-matrix(SMP$latent[start],nrow=dim(data)[1],ncol=1)
      
    } else if (lab == timeVar) {
      # simple linear time
      expr <- parse(text = lab)[[1]]
      
      Xlab<-matrix(eval(expr, envir = data),ncol=1,nrow=dim(data)[1])
      colnames(Xlab)<-paste0(lab,"_L1")
      start<-ct$start[which(ct$tag==paste0(lab,"_L1"))]
      B[,k]<-matrix(SMP$latent[start],nrow=dim(data)[1],ncol=1)
      
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      #INLA takes only function no : I(time^2)
      expr <- parse(text = lab)[[1]]
      
      Xlab<-matrix(eval(expr, envir = data),ncol=1,nrow=dim(data)[1])
      colnames(Xlab)<-paste0(gsub("[())]","",lab),"_L1")
      start<-ct$start[which(ct$tag==paste0(gsub("[())]","",lab),"_L1"))]
      B[,k]<-matrix(SMP$latent[start],nrow=dim(data)[1],ncol=1)
      
    } else if (use_splines && grepl("bs\\(|ns\\(", lab)) {
      
      # Extract knots etc. from the original call if needed
      stop("Spline not supported so far")
      
    } 
    X<-cbind(X,Xlab)
    k<-k+1
  }
  
  
  idd<-unique(data[,colnames(data)%in% id])
  X_RE<-NULL
  B_RE<-matrix(0,nrow=dim(data)[1],ncol=length(terms_RE))
  n_id<-length(unique(data[,colnames(data)%in%id]))
  j<-0
  k<-1
  for (lab in terms_RE) {
    
    if(lab=="Intercept"){
      
      Xlab <- matrix(1,nrow=dim(data)[1],ncol=1)
      colnames(Xlab)<-paste0(id,lab,"_L1")
      
      start<-ct$start[which(ct$tag==paste0(id,lab,"_L1"))]
      start<-start+j*n_id
      end<-start+n_id-1
      B_RE[,k]<-matrix(do.call(c,lapply(c(start:end),FUN=function(x){
        nn<-x-start+1
        nn<-sum(data[,colnames(data)%in%id]%in%idd[nn])
        return(rep(SMP$latent[x],nn))
      })),nrow=dim(data)[1],ncol=1)
      
    } else if (lab == timeVar) {
      # simple linear time
      expr <- parse(text = lab)[[1]]
      Xlab<-matrix(eval(expr, envir = data),ncol=1,nrow=dim(data)[1])
      colnames(Xlab)<-paste0(id,lab,"_L1")
      
      start<-ct$start[which(ct$tag==paste0(id,lab,"_L1"))]
      start<-start+j*n_id
      end<-start+n_id-1
      B_RE[,k]<-matrix(do.call(c,lapply(c(start:end),FUN=function(x){
        nn<-x-start+1
        nn<-sum(data[,colnames(data)%in%id]%in%idd[nn])
        return(rep(SMP$latent[x],nn))
      })),nrow=dim(data)[1],ncol=1)
      
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      #INLA takes only function no : I(time^2)
      
      expr <- parse(text = lab)[[1]]
      Xlab<-matrix(eval(expr, envir = data),ncol=1,nrow=dim(data)[1])
      colnames(Xlab)<-paste0(id,gsub("[())]","",lab),"_L1")
      
      start<-ct$start[which(ct$tag==paste0(id,gsub("[())]","",lab),"_L1"))]
      start<-start+j*n_id
      end<-start+n_id-1
      B_RE[,k]<-matrix(do.call(c,lapply(c(start:end),FUN=function(x){
        nn<-x-start+1
        nn<-sum(data[,colnames(data)%in%id]%in%idd[nn])
        return(rep(SMP$latent[x],nn))
      })),nrow=dim(data)[1],ncol=1)
      
      
    } else if (use_splines && grepl("bs\\(|ns\\(", lab)) {
      
      # Extract knots etc. from the original call if needed
      stop("Spline not supported so far")
      
    } 
    X_RE<-cbind(X_RE,Xlab)
    k<-k+1
    j<-j+1
  }
  
  Y<-X_RE*B_RE + X*B
  Y<-rowSums(Y)
  
  return(Y)
}

make_REXINLA <- function(formula, timeVar, data,use_splines = FALSE,ct,id,SMP, ...) {
  
  
  terms_labels <- attr(terms(formula), "term.labels")
  
  terms_RE <- terms_labels[grepl(id, terms_labels)==T][1]
  terms_RE <- gsub("\\|.*", "", terms_RE)
  terms_RE <- as.formula(paste("~", terms_RE))
  terms_RE <- attr(terms(terms_RE), "term.labels")
  
  
  # Identify terms
  if(paste0(id,"Intercept_L1")%in%ct$tag){
    terms_RE<-c("Intercept",terms_RE)
  }
  
  idd<-unique(data[,colnames(data)%in% id])
  B_RE<-matrix(0,nrow=dim(data)[1],ncol=length(terms_RE))
  n_id<-length(unique(data[,colnames(data)%in%id]))
  j<-0
  k<-1
  names_RE<-rep(NA,length(terms_RE))

  for (lab in terms_RE) {
    
    if(lab=="Intercept"){
      
      start<-ct$start[which(ct$tag==paste0(id,lab,"_L1"))]
      start<-start+j*n_id
      end<-start+n_id-1
      B_RE[,k]<-matrix(do.call(c,lapply(c(start:end),FUN=function(x){
        nn<-x-start+1
        nn<-sum(data[,colnames(data)%in%id]%in%idd[nn])
        return(rep(SMP$latent[x],nn))
      })),nrow=dim(data)[1],ncol=1)
      
      names_RE[k]<-lab
    } else if (lab == timeVar) {
      # simple linear time
      start<-ct$start[which(ct$tag==paste0(id,lab,"_L1"))]
      start<-start+j*n_id
      end<-start+n_id-1
      B_RE[,k]<-matrix(do.call(c,lapply(c(start:end),FUN=function(x){
        nn<-x-start+1
        nn<-sum(data[,colnames(data)%in%id]%in%idd[nn])
        return(rep(SMP$latent[x],nn))
      })),nrow=dim(data)[1],ncol=1)
      
      names_RE[k]<-lab
      
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      #INLA takes only function no : I(time^2)
      
      start<-ct$start[which(ct$tag==paste0(id,gsub("[())]","",lab),"_L1"))]
      start<-start+j*n_id
      end<-start+n_id-1
      B_RE[,k]<-matrix(do.call(c,lapply(c(start:end),FUN=function(x){
        nn<-x-start+1
        nn<-sum(data[,colnames(data)%in%id]%in%idd[nn])
        return(rep(SMP$latent[x],nn))
      })),nrow=dim(data)[1],ncol=1)
      
      names_RE[k]<-lab
      
    } else if (use_splines && grepl("bs\\(|ns\\(", lab)) {
      
      # Extract knots etc. from the original call if needed
      stop("Spline not supported so far")
      
    } 

    k<-k+1
    j<-j+1
  }
  
  colnames(B_RE)<-names_RE

  return(Y)
}

make_dXINLA <- function(formula,timeVar, data,use_splines = FALSE,ct,id,SMP, ...) {
  
  
  terms_labels <- attr(terms(formula), "term.labels")
  terms_fixed <- terms_labels[grepl(id, terms_labels)==F]
  
  terms_RE <- terms_labels[grepl(id, terms_labels)==T][1]
  terms_RE <- gsub("\\|.*", "", terms_RE)
  terms_RE <- as.formula(paste("~", terms_RE))
  terms_RE <- attr(terms(terms_RE), "term.labels")
  

  # Identify terms
  X<-NULL
  B<-matrix(0,nrow=dim(data)[1],ncol=length(terms_fixed))

  k<-1
  for (lab in terms_fixed) {
    if (lab == timeVar) {
      # simple linear time
      Xlab <- matrix(1,nrow=dim(data)[1],ncol=1)
      colnames(Xlab)<-paste0(lab,"_L1")
      
      start<-ct$start[which(ct$tag==paste0(lab,"_L1"))]
      B[,k]<-matrix(SMP$latent[start],nrow=dim(data)[1],ncol=1)
      
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      #INLA takes only function no : I(time^2)
      expr <- parse(text = lab)[[1]]
      dexpr <- Deriv::Deriv(expr, timeVar)
      
      Xlab<-matrix(eval(dexpr, envir = data),ncol=1,nrow=dim(data)[1])
      colnames(Xlab)<-paste0(gsub("[())]","",lab),"_L1")
      start<-ct$start[which(ct$tag==paste0(gsub("[())]","",lab),"_L1"))]
      B[,k]<-matrix(SMP$latent[start],nrow=dim(data)[1],ncol=1)
      
    } else if (use_splines && grepl("bs\\(|ns\\(", lab)) {
      
      # Extract knots etc. from the original call if needed
      stop("Spline not supported so far")
      
    } 
    X<-cbind(X,Xlab)
    k<-k+1
  }
  

  idd<-unique(data[,colnames(data)%in% id])
  X_RE<-NULL
  B_RE<-matrix(0,nrow=dim(data)[1],ncol=length(terms_RE))
  n_id<-length(unique(data[,colnames(data)%in%id]))
  j<-ifelse(paste0(id,"Intercept_L1")%in%ct$tag,1,0)

  k<-1
  for (lab in terms_RE) {
    
    if (lab == timeVar) {
      # simple linear time
      Xlab <- matrix(1,nrow=dim(data)[1],ncol=1)
      colnames(Xlab)<-paste0(id,lab,"_L1")
      
      start<-ct$start[which(ct$tag==paste0(id,lab,"_L1"))]
      start<-start+j*n_id
      end<-start+n_id-1
      B_RE[,k]<-matrix(do.call(c,lapply(c(start:end),FUN=function(x){
        nn<-x-start+1
        nn<-sum(data[,colnames(data)%in%id]%in%idd[nn])
        return(rep(SMP$latent[x],nn))
      })),nrow=dim(data)[1],ncol=1)
      
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      #INLA takes only function no : I(time^2)
      
      expr <- parse(text = lab)[[1]]
      dexpr <- Deriv::Deriv(expr, timeVar)
      Xlab<-matrix(eval(dexpr, envir = data),ncol=1,nrow=dim(data)[1])
      colnames(Xlab)<-paste0(id,gsub("[())]","",lab),"_L1")
      
      start<-ct$start[which(ct$tag==paste0(id,gsub("[())]","",lab),"_L1"))]
      start<-start+j*n_id
      end<-start+n_id-1
      B_RE[,k]<-matrix(do.call(c,lapply(c(start:end),FUN=function(x){
        nn<-x-start+1
        nn<-sum(data[,colnames(data)%in%id]%in%idd[nn])
        return(rep(SMP$latent[x],nn))
      })),nrow=dim(data)[1],ncol=1)
      
      
    } else if (use_splines && grepl("bs\\(|ns\\(", lab)) {
      
      # Extract knots etc. from the original call if needed
      stop("Spline not supported so far")
      
    } 
    X_RE<-cbind(X_RE,Xlab)
    k<-k+1
    j<-j+1
  }
  
  dY<-X_RE*B_RE + X*B
  dY<-rowSums(dY)
  
  return(dY)
}

make_XINLA_BLUP <- function(formula, timeVar, data,use_splines = FALSE,ct,id,SMP, ...) {
  
  
  terms_labels <- attr(terms(formula), "term.labels")
  terms_fixed <- terms_labels[grepl(id, terms_labels)==F]
  
  terms_RE <- terms_labels[grepl(id, terms_labels)==T][1]
  terms_RE <- gsub("\\|.*", "", terms_RE)
  terms_RE <- as.formula(paste("~", terms_RE))
  terms_RE <- attr(terms(terms_RE), "term.labels")
  
  
  # Identify terms
  if(paste0(id,"Intercept_L1")%in%ct$tag){
    terms_RE<-c("Intercept",terms_RE)
  }
  if("Intercept_L1"%in%ct$tag){
    terms_fixed<-c("Intercept",terms_fixed)
    
  }
  
  X<-NULL
  B<-matrix(0,nrow=dim(data)[1],ncol=length(terms_fixed))
  k<-1
  for (lab in terms_fixed) {
    
    if(lab=="Intercept"){
      
      Xlab <- matrix(1,nrow=dim(data)[1],ncol=1)
      colnames(Xlab)<-paste0(lab,"_L1")
      
      start<-which(rownames(SMP$summary.fixed)==paste0(lab,"_L1"))
      B[,k]<-matrix(SMP$summary.fixed[start,"mode"],nrow=dim(data)[1],ncol=1)
      
      
    } else if (lab == timeVar) {
      # simple linear time
      expr <- parse(text = lab)[[1]]
      
      Xlab<-matrix(eval(expr, envir = data),ncol=1,nrow=dim(data)[1])
      colnames(Xlab)<-paste0(lab,"_L1")
      
      start<-which(rownames(SMP$summary.fixed)==paste0(lab,"_L1"))
      B[,k]<-matrix(SMP$summary.fixed[start,"mode"],nrow=dim(data)[1],ncol=1)
      
      
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      #INLA takes only function no : I(time^2)
      expr <- parse(text = lab)[[1]]
      
      Xlab<-matrix(eval(expr, envir = data),ncol=1,nrow=dim(data)[1])
      colnames(Xlab)<-paste0(gsub("[())]","",lab),"_L1")
      
      start<-which(rownames(SMP$summary.fixed)==paste0(gsub("[())]","",lab),"_L1"))
      B[,k]<-matrix(SMP$summary.fixed[start,"mode"],nrow=dim(data)[1],ncol=1)
      
      
    } else if (use_splines && grepl("bs\\(|ns\\(", lab)) {
      
      # Extract knots etc. from the original call if needed
      stop("Spline not supported so far")
      
    } 
    X<-cbind(X,Xlab)
    k<-k+1
  }
  
  
  idd<-unique(data[,colnames(data)%in% id])
  X_RE<-NULL
  B_RE<-matrix(0,nrow=dim(data)[1],ncol=length(terms_RE))
  n_id<-length(unique(data[,colnames(data)%in%id]))
  j<-0
  k<-1
  for (lab in terms_RE) {
    
    if(lab=="Intercept"){
      
      Xlab <- matrix(1,nrow=dim(data)[1],ncol=1)
      colnames(Xlab)<-paste0(id,lab,"_L1")
      
      start<-which(names(SMP$summary.random)==paste0(id,lab,"_L1"))
      
      B_RE[,k]<-matrix(do.call(c,lapply(c(1:n_id),FUN=function(x){
        nn<-sum(data[,colnames(data)%in%id]%in%idd[x])
        return(rep(SMP$summary.random[[start]][x,"mode"],nn))
      })),nrow=dim(data)[1],ncol=1)
      
      
    } else if (lab == timeVar) {
      # simple linear time
      expr <- parse(text = lab)[[1]]
      Xlab<-matrix(eval(expr, envir = data),ncol=1,nrow=dim(data)[1])
      colnames(Xlab)<-paste0(id,lab,"_L1")
      
      start<-which(names(SMP$summary.random)==paste0(id,lab,"_L1"))
      
      B_RE[,k]<-matrix(do.call(c,lapply(c(1:n_id),FUN=function(x){
        nn<-sum(data[,colnames(data)%in%id]%in%idd[x])
        return(rep(SMP$summary.random[[start]][x,"mode"],nn))
      })),nrow=dim(data)[1],ncol=1)
      
      
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      #INLA takes only function no : I(time^2)
      
      expr <- parse(text = lab)[[1]]
      Xlab<-matrix(eval(expr, envir = data),ncol=1,nrow=dim(data)[1])
      colnames(Xlab)<-paste0(id,gsub("[())]","",lab),"_L1")
      

      start<-which(names(SMP$summary.random)==paste0(id,gsub("[())]","",lab),"_L1"))
      
      B_RE[,k]<-matrix(do.call(c,lapply(c(1:n_id),FUN=function(x){
        nn<-sum(data[,colnames(data)%in%id]%in%idd[x])
        return(rep(SMP$summary.random[[start]][x,"mode"],nn))
      })),nrow=dim(data)[1],ncol=1)
      
      
      
    } else if (use_splines && grepl("bs\\(|ns\\(", lab)) {
      
      # Extract knots etc. from the original call if needed
      stop("Spline not supported so far")
      
    } 
    X_RE<-cbind(X_RE,Xlab)
    k<-k+1
    j<-j+1
  }
  
  Y<-X_RE*B_RE + X*B
  Y<-rowSums(Y)
  
  return(Y)
}

make_REXINLA_BLUP <- function(formula, timeVar, data,use_splines = FALSE,ct,id,SMP, ...) {
  
  
  terms_labels <- attr(terms(formula), "term.labels")
  
  terms_RE <- terms_labels[grepl(id, terms_labels)==T][1]
  terms_RE <- gsub("\\|.*", "", terms_RE)
  terms_RE <- as.formula(paste("~", terms_RE))
  terms_RE <- attr(terms(terms_RE), "term.labels")
  
  
  # Identify terms
  if(paste0(id,"Intercept_L1")%in%ct$tag){
    terms_RE<-c("Intercept",terms_RE)
  }

  
  idd<-unique(data[,colnames(data)%in% id])
  B_RE<-matrix(0,nrow=dim(data)[1],ncol=length(terms_RE))
  n_id<-length(unique(data[,colnames(data)%in%id]))
  j<-0
  k<-1

  names_RE<-rep(NA,length(terms_RE))
  for (lab in terms_RE) {
    browser()
    if(lab=="Intercept"){
      
      start<-which(names(SMP$summary.random)==paste0(id,lab,"_L1"))
      
      B_RE[,k]<-matrix(do.call(c,lapply(c(1:n_id),FUN=function(x){
        nn<-sum(data[,colnames(data)%in%id]%in%idd[x])
        return(rep(SMP$summary.random[[start]][x,"mode"],nn))
      })),nrow=dim(data)[1],ncol=1)
      
      names_RE[k]<-lab
    } else if (lab == timeVar) {
      # simple linear time

      start<-which(names(SMP$summary.random)==paste0(id,lab,"_L1"))
      
      B_RE[,k]<-matrix(do.call(c,lapply(c(1:n_id),FUN=function(x){
        nn<-sum(data[,colnames(data)%in%id]%in%idd[x])
        return(rep(SMP$summary.random[[start]][x,"mode"],nn))
      })),nrow=dim(data)[1],ncol=1)
      
      names_RE[k]<-lab
      
      
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      #INLA takes only function no : I(time^2)
      
      
      start<-which(names(SMP$summary.random)==paste0(id,gsub("[())]","",lab),"_L1"))
      
      B_RE[,k]<-matrix(do.call(c,lapply(c(1:n_id),FUN=function(x){
        nn<-sum(data[,colnames(data)%in%id]%in%idd[x])
        return(rep(SMP$summary.random[[start]][x,"mode"],nn))
      })),nrow=dim(data)[1],ncol=1)
      names_RE[k]<-lab
      
      
      
    } else if (use_splines && grepl("bs\\(|ns\\(", lab)) {
      
      # Extract knots etc. from the original call if needed
      stop("Spline not supported so far")
      
    } 
    k<-k+1
    j<-j+1
  }
  browser()
  colnames(B_RE)<-names_RE
  
  return(Y)
}

make_dXINLA_BLUP <- function(formula,timeVar, data,use_splines = FALSE,ct,id,SMP, ...) {
  
  
  terms_labels <- attr(terms(formula), "term.labels")
  terms_fixed <- terms_labels[grepl(id, terms_labels)==F]
  
  terms_RE <- terms_labels[grepl(id, terms_labels)==T][1]
  terms_RE <- gsub("\\|.*", "", terms_RE)
  terms_RE <- as.formula(paste("~", terms_RE))
  terms_RE <- attr(terms(terms_RE), "term.labels")
  
  
  # Identify terms
  X<-NULL
  B<-matrix(0,nrow=dim(data)[1],ncol=length(terms_fixed))
  
  k<-1
  for (lab in terms_fixed) {
    if (lab == timeVar) {
      # simple linear time
      Xlab <- matrix(1,nrow=dim(data)[1],ncol=1)
      colnames(Xlab)<-paste0(lab,"_L1")
      
      start<-which(rownames(SMP$summary.fixed)==paste0(lab,"_L1"))
      B[,k]<-matrix(SMP$summary.fixed[start,"mode"],nrow=dim(data)[1],ncol=1)
      
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      #INLA takes only function no : I(time^2)
      expr <- parse(text = lab)[[1]]
      dexpr <- Deriv::Deriv(expr, timeVar)
      
      Xlab<-matrix(eval(dexpr, envir = data),ncol=1,nrow=dim(data)[1])
      colnames(Xlab)<-paste0(gsub("[())]","",lab),"_L1")
      start<-which(rownames(SMP$summary.fixed)==paste0(gsub("[())]","",lab),"_L1"))
      B[,k]<-matrix(SMP$summary.fixed[start,"mode"],nrow=dim(data)[1],ncol=1)
      
    } else if (use_splines && grepl("bs\\(|ns\\(", lab)) {
      
      # Extract knots etc. from the original call if needed
      stop("Spline not supported so far")
      
    } 
    X<-cbind(X,Xlab)
    k<-k+1
  }
  
  
  idd<-unique(data[,colnames(data)%in% id])
  X_RE<-NULL
  B_RE<-matrix(0,nrow=dim(data)[1],ncol=length(terms_RE))
  n_id<-length(unique(data[,colnames(data)%in%id]))
  j<-ifelse(paste0(id,"Intercept_L1")%in%ct$tag,1,0)
  
  k<-1
  for (lab in terms_RE) {
    
    if (lab == timeVar) {
      # simple linear time
      Xlab <- matrix(1,nrow=dim(data)[1],ncol=1)
      colnames(Xlab)<-paste0(id,lab,"_L1")
      
      start<-which(names(SMP$summary.random)==paste0(id,lab,"_L1"))
      
      B_RE[,k]<-matrix(do.call(c,lapply(c(1:n_id),FUN=function(x){
        nn<-sum(data[,colnames(data)%in%id]%in%idd[x])
        return(rep(SMP$summary.random[[start]][x,"mode"],nn))
      })),nrow=dim(data)[1],ncol=1)
      
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      #INLA takes only function no : I(time^2)
      
      expr <- parse(text = lab)[[1]]
      dexpr <- Deriv::Deriv(expr, timeVar)
      Xlab<-matrix(eval(dexpr, envir = data),ncol=1,nrow=dim(data)[1])
      colnames(Xlab)<-paste0(id,gsub("[())]","",lab),"_L1")
      
      start<-which(names(SMP$summary.random)==paste0(id,gsub("[())]","",lab),"_L1"))
      
      B_RE[,k]<-matrix(do.call(c,lapply(c(1:n_id),FUN=function(x){
        nn<-sum(data[,colnames(data)%in%id]%in%idd[x])
        return(rep(SMP$summary.random[[start]][x,"mode"],nn))
      })),nrow=dim(data)[1],ncol=1)
      
      
    } else if (use_splines && grepl("bs\\(|ns\\(", lab)) {
      
      # Extract knots etc. from the original call if needed
      stop("Spline not supported so far")
      
    } 
    X_RE<-cbind(X_RE,Xlab)
    k<-k+1
    j<-j+1
  }
  
  dY<-X_RE*B_RE + X*B
  dY<-rowSums(dY)
  
  return(dY)
}
