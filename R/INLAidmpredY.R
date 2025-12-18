### Code:
##' @title Calculate predictions for time-depend covariates using INLA
#' @importFrom foreach "%do%"
#' @importFrom foreach "%dopar%"
#' @importFrom Deriv "Deriv"
#' @importFrom INLA "inla.posterior.sample"
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  


INLAidmpredY<-function(timeVar,truncated,formLong,dataSurv,dataLongi,id,
                  t0,t1,t2,t3,assoc,
                  ctime,modelY,seed,BLUP,nproc,clustertype,scale.X){

  
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

  #
  if(scale.X==T){
    tcenter<-ifelse(truncated==T,0,min(t0))
    dataCenter<-data.frame(ID=idsubjects,time=tcenter)
    colnames(dataCenter)<-c(id,timeVar)
  }
  #browser()
  # should center by median or mean ? 
  if(nproc==1){
    
    
    for(indice in 1:length(formLong)){

      INLAmodel<-modelY$modelY[[indice]]
      
      # data structure
      ct <- INLAmodel$misc$configs$contents
      
      if(is.null(INLAmodel)){stop("The inla model for your marker could not be run, see above warnings.")}
      choiceY<-na.omit(unlist(assoc[[indice]]))
      
      
      if(BLUP==F){
        
        
        #samples seed=seed cannot do parallel estimation on it 
        SMP <- INLA::inla.posterior.sample(1, INLAmodel,seed=seed)
        
        res<-NULL
        key1 <- do.call(paste, c(dataLongi_augmented[,colnames(dataLongi_augmented)%in%c(id,timeVar)], sep = "\r"))
        key2 <- do.call(paste, c(timePointsdata, sep = "\r"))
        # keep only indice we want: 
        # Collapse each row into a string
        
        # Find indices of rows from data1 that are in data2
        # while keeping order of key2
        indices <- match(key2, key1)
        
        if("value"%in% choiceY){
        
        Y<-as.matrix(make_XINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=SMP[[1]]))
        Y<-Y[indices,]
        Outcome<-all.vars(terms(formLong[[indice]]))[1]
        PredYx<-cbind(timePointsdata,Outcome=Outcome,Y)
        if(scale.X==T){
          Ycenter<-make_XINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=SMP[[1]])
          PredYx$Y<-(PredYx$Y-mean(Ycenter))/sd(Ycenter)
        }
        colnames(PredYx)[4]<-"Sample_1"
        res<-rbind(res,PredYx)
        }
        
        if("slope"%in% choiceY){
        dY<-as.matrix(make_dXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=SMP[[1]]))
        dY<-dY[indices,]
        slopeOutcome<-paste0("slope_",Outcome)
        slopePredYx<-cbind(timePointsdata,Outcome=slopeOutcome,dY)
        if(scale.X==T){
          dYcenter<-make_dXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=SMP[[1]])
          slopePredYx$dY<-(slopePredYx$dY-mean(dYcenter))/sd(dYcenter)
        }
        colnames(slopePredYx)[4]<-"Sample_1"
        res<-rbind(res,slopePredYx)
        }
        
        if("RE"%in% choiceY){
          
        REY<-as.matrix(make_REXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=SMP[[1]]))
       REY<-REY[indices,]
       
       if(scale.X==T){
         REYcenter<-as.matrix(make_REXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=SMP[[1]]))
         REY<-do.call(cbind,lapply(c(1:dim(REY)[2]),FUN=function(x){
           (REY[,x]-mean(REYcenter[,x]))/sd(REYcenter[,x])
         }))
       }
       namesREY<-unlist(lapply(1:dim(REY)[2],FUN=function(x){
         rep(paste0("RE_",colnames(REY)[x],"_",Outcome),dim(REY)[1])
       }))
       dataREY<- do.call(rbind, replicate(dim(REY)[2], timePointsdata, simplify = FALSE))
       REPredYx<-cbind(dataREY,Outcome=namesREY,as.vector(REY))
       
       
       colnames(REPredYx)[4]<-"Sample_1"
       res<-rbind(res,REPredYx)
        }
        
        Yall[[indice]]<- res
        
        
        
      }else{
        
        
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
          
          if(scale.X==T){
            Ycenter<-make_XINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=INLAmodel)
            PredYx$Y<-(PredYx$Y-mean(Ycenter))/sd(Ycenter)
          }
          colnames(PredYx)[4]<-"Sample_1"
          res<-rbind(res,PredYx)
        }
        
        if("slope"%in% choiceY){
        dY<-as.matrix(make_dXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=INLAmodel))
        dY<-dY[indices,]
        slopeOutcome<-paste0("slope_",Outcome)
        slopePredYx<-cbind(timePointsdata,Outcome=slopeOutcome,dY)
        
        if(scale.X==T){
          dYcenter<-make_dXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=INLAmodel)
          slopePredYx$dY<-(slopePredYx$dY-mean(dYcenter))/sd(dYcenter)
        }
        colnames(slopePredYx)[4]<-"Sample_1"
        res<-rbind(res,slopePredYx)
        }
        
        if("RE"%in% choiceY){
        REY<-as.matrix(make_REXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=INLAmodel))
        REY<-REY[indices,]
        
        if(scale.X==T){
          REYcenter<-as.matrix(make_REXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=INLAmodel))
          REY<-do.call(cbind,lapply(c(1:dim(REY)[2]),FUN=function(x){
            (REY[,x]-mean(REYcenter[,x]))/sd(REYcenter[,x])
          }))
        }
        namesREY<-unlist(lapply(1:dim(REY)[2],FUN=function(x){
          rep(paste0("RE_",colnames(REY)[x],"_",Outcome),dim(REY)[1])
        }))
        dataREY<- do.call(rbind, replicate(dim(REY)[2], timePointsdata, simplify = FALSE))
        REPredYx<-cbind(dataREY,Outcome=namesREY,as.vector(REY))
        colnames(REPredYx)[4]<-"Sample_1"
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
      
      
      Yall<-foreach::foreach(indice=1:length(formLong),
                       .combine = 'list',
                       packages=c("INLA","Deriv","HIDeM"))%dopar%{
                         
                         
                         INLAmodel<-modelY$modelY[[indice]]
                         
                         # data structure
                         ct <- INLAmodel$misc$configs$contents
                         
                         choiceY<-na.omit(unlist(assoc[[indice]]))
                         
                         if(is.null(INLAmodel)){stop("The inla model for your marker could not be run, see above warnings.")}
                         
                         
                         if(BLUP==F){
                           
                           
                           #samples seed=seed cannot do parallel estimation on it 
                           SMP <- INLA::inla.posterior.sample(1, INLAmodel,seed=seed)
                           
                           res<-NULL
                           key1 <- do.call(paste, c(dataLongi_augmented[,colnames(dataLongi_augmented)%in%c(id,timeVar)], sep = "\r"))
                           key2 <- do.call(paste, c(timePointsdata, sep = "\r"))
                           # keep only indice we want: 
                           # Collapse each row into a string
                           
                           # Find indices of rows from data1 that are in data2
                           # while keeping order of key2
                           indices <- match(key2, key1)
                           
                           if("value"%in% choiceY){
                             
                             Y<-make_XINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=SMP[[1]])
                             Y<-Y[indices,]
                             Outcome<-all.vars(terms(formLong[[indice]]))[1]
                             PredYx<-cbind(timePointsdata,Outcome=Outcome,Y)
                             
                             if(scale.X==T){
                               Ycenter<-make_XINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=SMP[[1]])
                               PredYx$Y<-(PredYx$Y-mean(Ycenter))/sd(Ycenter)
                             }
                             colnames(PredYx)[4]<-"Sample_1"
                             res<-rbind(res,PredYx)
                           }
                           
                           if("slope"%in% choiceY){
                             dY<-make_dXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=SMP[[1]])
                             dY<-dY[indices,]
                             slopeOutcome<-paste0("slope_",Outcome)
                             slopePredYx<-cbind(timePointsdata,Outcome=slopeOutcome,dY)
                             
                             if(scale.X==T){
                               dYcenter<-make_dXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=SMP[[1]])
                               slopePredYx$dY<-(slopePredYx$dY-mean(dYcenter))/sd(dYcenter)
                             }
                             colnames(slopePredYx)[4]<-"Sample_1"
                             res<-rbind(res,slopePredYx)
                           }
                           
                           if("RE"%in% choiceY){
                             
                             REY<-make_REXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=SMP[[1]])
                             REY<-REY[indices,]
                             
                             if(scale.X==T){
                               REYcenter<-as.matrix(make_REXINLA(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=SMP[[1]]))
                               REY<-do.call(cbind,lapply(c(1:dim(REY)[2]),FUN=function(x){
                                 (REY[,x]-mean(REYcenter[,x]))/sd(REYcenter[,x])
                               }))
                             }
                             namesREY<-unlist(lapply(1:dim(REY)[2],FUN=function(x){
                               rep(paste0("RE_",colnames(REY)[x],"_",Outcome),dim(REY)[1])
                             }))
                             dataREY<- do.call(rbind, replicate(dim(REY)[2], timePointsdata, simplify = FALSE))
                             REPredYx<-cbind(dataREY,Outcome=namesREY,as.vector(REY))
                             colnames(REPredYx)[4]<-"Sample_1"
                             res<-rbind(res,REPredYx)
                           }
                           
                           return(res)
                           
                           
                           
                         }else{
                           
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
                             
                             if(scale.X==T){
                               Ycenter<-make_XINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=INLAmodel)
                               PredYx$Y<-(PredYx$Y-mean(Ycenter))/sd(Ycenter)
                             }
                             colnames(PredYx)[4]<-"Sample_1"
                             res<-rbind(res,PredYx)
                           }
                           
                           if("slope"%in% choiceY){
                             dY<-as.matrix(make_dXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=INLAmodel))
                             dY<-dY[indices,]
                             slopeOutcome<-paste0("slope_",Outcome)
                             slopePredYx<-cbind(timePointsdata,Outcome=slopeOutcome,dY)
                             
                             if(scale.X==T){
                               dYcenter<-make_dXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=INLAmodel)
                               slopePredYx$dY<-(slopePredYx$dY-mean(dYcenter))/sd(dYcenter)
                             }
                             colnames(slopePredYx)[4]<-"Sample_1"
                             res<-rbind(res,slopePredYx)
                           }
                           
                           if("RE"%in% choiceY){
                             REY<-as.matrix(make_REXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataLongi_augmented,ct=ct,id=id,SMP=INLAmodel))
                             REY<-REY[indices,]
                             
                             if(scale.X==T){
                               REYcenter<-as.matrix(make_REXINLA_BLUP(formula=formLong[[indice]], timeVar=timeVar, data=dataCenter,ct=ct,id=id,SMP=INLAmodel))
                               REY<-do.call(cbind,lapply(c(1:dim(REY)[2]),FUN=function(x){
                                 (REY[,x]-mean(REYcenter[,x]))/sd(REYcenter[,x])
                               }))
                             }
                             namesREY<-unlist(lapply(1:dim(REY)[2],FUN=function(x){
                               rep(paste0("RE_",colnames(REY)[x],"_",Outcome),dim(REY)[1])
                             }))
                             dataREY<- do.call(rbind, replicate(dim(REY)[2], timePointsdata, simplify = FALSE))
                             REPredYx<-cbind(dataREY,Outcome=namesREY,as.vector(REY))
                             colnames(REPredYx)[4]<-"Sample_1"
                             res<-rbind(res,REPredYx)
                           }
                           
                           return(res)
                           
                           
                           
                         }
                           
                         
                         
                       }
      
      
      parallel::stopCluster(clustpar)
  }
  
  return(Yall[,colnames(Yall)%in%c(id,timeVar,"Outcome","Sample_1")])
  

  
}


make_XINLA <- function(formula, timeVar, data, use_splines = FALSE, ct, id, SMP, ...) {
  
  n <- nrow(data)
  
  # Extract terms
  terms_labels <- attr(terms(formula), "term.labels")
  terms_fixed <- terms_labels[!grepl(id, terms_labels)]
  
  terms_RE <- terms_labels[grepl(id, terms_labels)][1]
  terms_RE <- gsub("\\|.*", "", terms_RE)
  terms_RE <- attr(terms(as.formula(paste("~", terms_RE))), "term.labels")
  
  # Handle intercepts
  if (paste0(id, "Intercept_L1") %in% ct$tag) terms_RE <- c("Intercept", terms_RE)
  if ("Intercept_L1" %in% ct$tag) terms_fixed <- c("Intercept", terms_fixed)
  
  # Precompute
  ct_start <- setNames(ct$start, ct$tag)
  latent <- SMP$latent
  id_values <- data[[id]]
  id_levels <- unique(id_values)
  n_id <- length(id_levels)
  
  # --- Fixed effects ---
  n_fixed <- length(terms_fixed)
  X <- matrix(NA_real_, nrow = n, ncol = n_fixed)
  B <- matrix(0, nrow = n, ncol = n_fixed)
  colnames(X) <- paste0(terms_fixed, "_L1")
  
  for (k in seq_along(terms_fixed)) {
    lab <- terms_fixed[k]
    tag <- paste0(gsub("[()]", "", lab), "_L1")
    
    if (lab == "Intercept") {
      X[, k] <- 1
    } else {
      expr <- if (grepl("\\(", lab)) str2lang(lab) else as.symbol(lab)
      X[, k] <- eval(expr, data)
    }
    
    if (!is.na(ct_start[tag])) {
      B[, k] <- latent[ct_start[tag]]
    }
  }
  
  # --- Random effects ---
  n_re <- length(terms_RE)
  X_RE <- matrix(NA_real_, nrow = n, ncol = n_re)
  B_RE <- matrix(0, nrow = n, ncol = n_re)
  colnames(X_RE) <- paste0(id, terms_RE, "_L1")
  
  for (k in seq_along(terms_RE)) {
    lab <- terms_RE[k]
    tag <- paste0(id, gsub("[()]", "", lab), "_L1")
    base_start <- ct_start[tag]
    if (is.na(base_start)) next
    
    # Design part
    if (lab == "Intercept") {
      X_RE[, k] <- 1
    } else {
      expr <- if (grepl("\\(", lab)) str2lang(lab) else as.symbol(lab)
      X_RE[, k] <- eval(expr, data)
    }
    
    # <- THIS WAS THE MISSING PART!
    j_offset <- (k - 1) * n_id
    start <- base_start + j_offset
    end <- start + n_id - 1
    
    latents_for_ids <- latent[start:end]
    id_to_latent <- setNames(latents_for_ids, id_levels)
    B_RE[, k] <- unname(id_to_latent[id_values])
  }
  
  # Combine
  Y <- rowSums(X * B + X_RE * B_RE)
  return(Y)
}

make_REXINLA <- function(formula, timeVar, data, use_splines = FALSE, ct, id, SMP, ...) {
  n <- nrow(data)
  
  # --- Extract random-effect terms ---
  terms_labels <- attr(terms(formula), "term.labels")
  terms_RE <- terms_labels[grepl(id, terms_labels)][1]
  terms_RE <- gsub("\\|.*", "", terms_RE)
  terms_RE <- attr(terms(as.formula(paste("~", terms_RE))), "term.labels")
  
  # Add intercept if present in ct
  if (paste0(id, "Intercept_L1") %in% ct$tag) {
    terms_RE <- c("Intercept", terms_RE)
  }
  
  # Precompute
  ct_start <- setNames(ct$start, ct$tag)
  latent <- SMP$latent
  id_values <- data[[id]]
  id_levels <- unique(id_values)
  n_id <- length(id_levels)
  
  # --- Allocate output ---
  n_re <- length(terms_RE)
  B_RE <- matrix(0, nrow = n, ncol = n_re)
  colnames(B_RE) <- terms_RE
  
  # --- Main loop (vectorized) ---
  for (k in seq_along(terms_RE)) {
    lab <- terms_RE[k]
    
    # Match correct ct tag
    base_tag <- if (grepl("\\(", lab)) {
      paste0(id, gsub("[())]", "", lab), "_L1")
    } else {
      paste0(id, lab, "_L1")
    }
    base_start <- ct_start[base_tag]
    if (is.na(base_start)) next
    
    # Account for offset (important!)
    j_offset <- (k - 1) * n_id
    start <- base_start + j_offset
    end <- start + n_id - 1
    
    # Latent values for each ID level
    latents_for_ids <- latent[start:end]
    
    # Map each data row to its corresponding latent
    id_to_latent <- setNames(latents_for_ids, id_levels)
    B_RE[, k] <- unname(id_to_latent[id_values])
  }
  
  return(B_RE)
}

make_dXINLA<- function(formula, timeVar, data, use_splines = FALSE, ct, id, SMP, ...) {
  n <- nrow(data)
  
  # --- Extract fixed and random terms ---
  terms_labels <- attr(terms(formula), "term.labels")
  terms_fixed <- terms_labels[!grepl(id, terms_labels)]
  
  terms_RE <- terms_labels[grepl(id, terms_labels)][1]
  terms_RE <- gsub("\\|.*", "", terms_RE)
  terms_RE <- attr(terms(as.formula(paste("~", terms_RE))), "term.labels")
  
  # --- Precompute ---
  ct_start <- setNames(ct$start, ct$tag)
  latent <- SMP$latent
  id_values <- data[[id]]
  id_levels <- unique(id_values)
  n_id <- length(id_levels)
  
  # --- Fixed effects part ---
  n_fixed <- length(terms_fixed)
  X <- matrix(NA_real_, nrow = n, ncol = n_fixed)
  B <- matrix(0, nrow = n, ncol = n_fixed)
  colnames(X) <- paste0(terms_fixed, "_L1")
  
  for (k in seq_along(terms_fixed)) {
    lab <- terms_fixed[k]
    tag <- paste0(gsub("[()]", "", lab), "_L1")
    
    # Handle derivative wrt time
    if (lab == timeVar) {
      X[, k] <- 1
      if (!is.na(ct_start[tag])) B[, k] <- latent[ct_start[tag]]
      
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      expr <- str2lang(lab)
      dexpr <- Deriv::Deriv(expr, timeVar)
      X[, k] <- eval(dexpr, data)
      colnames(X)[k] <- paste0(gsub("[())]", "", lab), "_L1")
      if (!is.na(ct_start[tag])) B[, k] <- latent[ct_start[tag]]
    }
  }
  
  # --- Random effects part ---
  n_re <- length(terms_RE)
  X_RE <- matrix(NA_real_, nrow = n, ncol = n_re)
  B_RE <- matrix(0, nrow = n, ncol = n_re)
  colnames(X_RE) <- paste0(id, terms_RE, "_L1")
  
  j_offset_base <- ifelse(paste0(id, "Intercept_L1") %in% ct$tag, 1, 0)
  
  for (k in seq_along(terms_RE)) {
    lab <- terms_RE[k]
    tag <- paste0(id, gsub("[()]", "", lab), "_L1")
    base_start <- ct_start[tag]
    if (is.na(base_start)) next
    
    if (lab == timeVar) {
      X_RE[, k] <- 1
      
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      expr <- str2lang(lab)
      dexpr <- Deriv::Deriv(expr, timeVar)
      X_RE[, k] <- eval(dexpr, data)
      colnames(X_RE)[k] <- paste0(id, gsub("[())]", "", lab), "_L1")
    }
    
    # Offset logic identical to original
    j_offset <- (j_offset_base + (k - 1)) * n_id
    start <- base_start + j_offset
    end <- start + n_id - 1
    
    latents_for_ids <- latent[start:end]
    id_to_latent <- setNames(latents_for_ids, id_levels)
    B_RE[, k] <- unname(id_to_latent[id_values])
  }
  
  # --- Combine fixed and random derivatives ---
  dY <- rowSums(X * B + X_RE * B_RE)
  return(dY)
}

make_XINLA_BLUP <- function(formula, timeVar, data, use_splines = FALSE, ct, id, SMP, ...) {
  n <- nrow(data)
  N<-length(unique(data[,colnames(data)%in%id]))
  # --- parse terms ---
  terms_labels <- attr(terms(formula), "term.labels")
  terms_fixed <- terms_labels[!grepl(id, terms_labels)]
  terms_RE <- terms_labels[grepl(id, terms_labels)][1]
  terms_RE <- gsub("\\|.*", "", terms_RE)
  terms_RE <- attr(terms(as.formula(paste("~", terms_RE))), "term.labels")
  
  # intercept handling
  if (paste0(id, "Intercept_L1") %in% ct$tag) {
    terms_RE <- c("Intercept", terms_RE)
  }
  if ("Intercept_L1" %in% ct$tag) {
    terms_fixed <- c("Intercept", terms_fixed)
  }
  
  # --- prepare lookups ---
  # fixed summary modes by rowname
  sf <- SMP$summary.fixed
  fixed_mode <- if (!is.null(sf)) setNames(sf[, "mode"], rownames(sf)) else numeric(0)
  
  # random summary list names
  sr_names <- names(SMP$summary.random)
  sr_list <- SMP$summary.random
  
  id_values <- data[[id]]
  id_levels <- unique(id_values)        # preserves original order (like your idd)
  n_id <- length(id_levels)
  id_index <- match(id_values, id_levels) # vector of positions 1..n_id for each row
  
  # --- Fixed effects (design X and coefficient B) ---
  n_fixed <- length(terms_fixed)
  X <- matrix(NA_real_, nrow = n, ncol = n_fixed)
  B <- matrix(0, nrow = n, ncol = n_fixed)
  colnames(X) <- paste0(terms_fixed, "_L1")
  
  #browser()
  for (k in seq_along(terms_fixed)) {
    lab <- terms_fixed[k]
    tag <- paste0(gsub("[()]", "", lab), "_L1")
    
    if (lab == "Intercept") {
      X[, k] <- 1
    } else if (lab == timeVar) {
      # linear time
      expr <- if (grepl("\\(", lab)) str2lang(lab) else as.name(lab)
      X[, k] <- eval(expr, data)
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      expr <- str2lang(lab)
      X[, k] <- eval(expr, data)
      colnames(X)[k] <- paste0(gsub("[())]", "", lab), "_L1")
    } else {
      # other fixed terms (if any)
      X[, k] <- data[[lab]]
    }
    
    if (tag %in% names(fixed_mode)) {
      B[, k] <- fixed_mode[tag]
    }
  }
  
  # --- Random effects (B_RE) ---
  n_re <- length(terms_RE)
  X_RE <- matrix(NA_real_, nrow = n, ncol = n_re)
  B_RE <- matrix(0, nrow = n, ncol = n_re)
  colnames(X_RE) <- paste0(id, terms_RE, "_L1")
  
  for (k in seq_along(terms_RE)) {
    lab <- terms_RE[k]
    # form summary.random tag name exactly like original
    if (lab == "Intercept") {
      sr_tag <- paste0(id, lab, "_L1")
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      sr_tag <- paste0(id, gsub("[())]", "", lab), "_L1")
    } else {
      sr_tag <- paste0(id, lab, "_L1")
    }
    
    pos <- which(sr_names == sr_tag)
    if (length(pos) == 0) {
      # no random summary for that term -> leave zeros
      # still compute X_RE if needed
      if (lab == "Intercept") X_RE[, k] <- 1
      else if (lab == timeVar) X_RE[, k] <- eval(str2lang(lab), data)
      else if (grepl("\\(", lab) && grepl(timeVar, lab)) X_RE[, k] <- eval(str2lang(lab), data)
      next
    }
    
    # design column
    if (lab == "Intercept") {
      X_RE[, k] <- 1
    } else if (lab == timeVar) {
      X_RE[, k] <- eval(str2lang(lab), data)
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      X_RE[, k] <- eval(str2lang(lab), data)
      colnames(X_RE)[k] <- paste0(id, gsub("[())]", "", lab), "_L1")
    } else {
      X_RE[, k] <- data[[lab]]
    }
  
    # extract modes for each id (rows in the summary.random element)
    modes_vec <- sr_list[[pos]][((k-1)*N+1):(k*N), "mode"]
    # names correspond to id_levels order only if sr was built that way; original code used idd order (unique(data[[id]]))
    # so we set names using id_levels to replicate original replication logic
    names(modes_vec) <- id_levels
    B_RE[, k] <- unname(modes_vec[id_values])
  }
  
  # --- combine and return ---
  Y <- rowSums(X * B + X_RE * B_RE)
  return(Y)
}

make_REXINLA_BLUP <- function(formula, timeVar, data, use_splines = FALSE, ct, id, SMP, ...) {
  # --- parse random-effect terms ---
  n <- nrow(data)
  N<-length(unique(data[,colnames(data)%in%id]))
  terms_labels <- attr(terms(formula), "term.labels")
  terms_RE <- terms_labels[grepl(id, terms_labels)][1]
  terms_RE <- gsub("\\|.*", "", terms_RE)
  terms_RE <- attr(terms(as.formula(paste("~", terms_RE))), "term.labels")
  
  # add intercept if present in ct$tag
  if (paste0(id, "Intercept_L1") %in% ct$tag)
    terms_RE <- c("Intercept", terms_RE)
  
  # --- precompute ID mapping ---
  id_values <- data[[id]]
  id_levels <- unique(id_values)       # preserve original order (like your idd)
  n <- nrow(data)
  n_id <- length(id_levels)
  id_index <- match(id_values, id_levels) # each row gets an index 1..n_id
  
  # --- prepare output matrix ---
  n_re <- length(terms_RE)
  B_RE <- matrix(0, nrow = n, ncol = n_re)
  names_RE <- character(n_re)
  
  # --- get list names once ---
  sr_names <- names(SMP$summary.random)
  sr_list  <- SMP$summary.random
  
  # --- main loop (vectorized per term) ---
  for (k in seq_along(terms_RE)) {
    lab <- terms_RE[k]
    
    # Build tag name for summary.random lookup
    if (lab == "Intercept") {
      sr_tag <- paste0(id, lab, "_L1")
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      sr_tag <- paste0(id, gsub("[())]", "", lab), "_L1")
    } else {
      sr_tag <- paste0(id, lab, "_L1")
    }
    
    pos <- which(sr_names == sr_tag)
    if (length(pos) == 0) next  # skip missing RE terms
    
    # Extract BLUP modes for each ID
    modes_vec <- sr_list[[pos]][((k-1)*N+1):(k*N), "mode"]
    names(modes_vec) <- id_levels  # align with unique(data[[id]]) order
    
    # Vectorized mapping
    B_RE[, k] <- unname(modes_vec[id_values])
    
    names_RE[k] <- lab
  }
  
  colnames(B_RE) <- names_RE
  return(B_RE)
}
make_dXINLA_BLUP <- function(formula, timeVar, data, use_splines = FALSE, ct, id, SMP, ...) {
  n <- nrow(data)
  N<-length(unique(data[,colnames(data)%in%id]))
  # parse terms
  terms_labels <- attr(terms(formula), "term.labels")
  terms_fixed <- terms_labels[!grepl(id, terms_labels)]
  terms_RE <- terms_labels[grepl(id, terms_labels)][1]
  terms_RE <- gsub("\\|.*", "", terms_RE)
  terms_RE <- attr(terms(as.formula(paste("~", terms_RE))), "term.labels")
  
  # precompute lookups
  sf <- SMP$summary.fixed
  fixed_mode <- if (!is.null(sf)) setNames(sf[, "mode"], rownames(sf)) else numeric(0)
  sr_names <- names(SMP$summary.random)
  sr_list  <- SMP$summary.random
  
  id_values <- data[[id]]
  id_levels <- unique(id_values)      # preserve original order (like idd)
  n_id <- length(id_levels)
  
  # --- Fixed effects derivative (X and B) ---
  n_fixed <- length(terms_fixed)
  X <- matrix(NA_real_, nrow = n, ncol = n_fixed)
  B <- matrix(0, nrow = n, ncol = n_fixed)
  colnames(X) <- paste0(terms_fixed, "_L1")
  
  for (k in seq_along(terms_fixed)) {
    lab <- terms_fixed[k]
    tag <- paste0(gsub("[()]", "", lab), "_L1")
    
    if (lab == timeVar) {
      # derivative of linear time is 1
      X[, k] <- 1
      if (tag %in% names(fixed_mode)) B[, k] <- fixed_mode[tag]
      
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      # derivative of function of time
      expr <- str2lang(lab)
      dexpr <- Deriv::Deriv(expr, timeVar)
      X[, k] <- eval(dexpr, envir = data)
      colnames(X)[k] <- paste0(gsub("[())]", "", lab), "_L1")
      if (tag %in% names(fixed_mode)) B[, k] <- fixed_mode[tag]
      
    } else {
      # other fixed (not depending on time) -> derivative 0 (leave X as NA? original didn't handle else)
      # original only treated timeVar or function of time; keep behavior by leaving NA -> treat as 0 contribution
      X[, k] <- 0
      if (tag %in% names(fixed_mode)) B[, k] <- fixed_mode[tag]
    }
  }
  
  # --- Random effects derivative (X_RE and B_RE) ---
  n_re <- length(terms_RE)
  X_RE <- matrix(NA_real_, nrow = n, ncol = n_re)
  B_RE <- matrix(0, nrow = n, ncol = n_re)
  colnames(X_RE) <- paste0(id, terms_RE, "_L1")
  
  for (k in seq_along(terms_RE)) {
    lab <- terms_RE[k]
    
    # build summary.random tag exactly as original
    if (lab == "Intercept") {
      sr_tag <- paste0(id, lab, "_L1")
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      sr_tag <- paste0(id, gsub("[())]", "", lab), "_L1")
    } else {
      sr_tag <- paste0(id, lab, "_L1")
    }
    
    pos <- which(sr_names == sr_tag)
    if (length(pos) == 0) {
      # no summary.random for term -> X_RE maybe computed, B_RE remains zero
      if (lab == "Intercept") X_RE[, k] <- 1
      else if (lab == timeVar) X_RE[, k] <- 1
      else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
        expr <- str2lang(lab)
        dexpr <- Deriv::Deriv(expr, timeVar)
        X_RE[, k] <- eval(dexpr, envir = data)
      } else X_RE[, k] <- 0
      next
    }
    
    # design derivative
    if (lab == "Intercept") {
      X_RE[, k] <- 1
    } else if (lab == timeVar) {
      X_RE[, k] <- 1
    } else if (grepl("\\(", lab) && grepl(timeVar, lab)) {
      expr <- str2lang(lab)
      dexpr <- Deriv::Deriv(expr, timeVar)
      X_RE[, k] <- eval(dexpr, envir = data)
      colnames(X_RE)[k] <- paste0(id, gsub("[())]", "", lab), "_L1")
    } else {
      X_RE[, k] <- 0
    }
    
    # extract BLUP modes for each id (the sr_list[[pos]] rows correspond to ids in original order)
    modes_vec <- sr_list[[pos]][((k-1)*N+1):(k*N), "mode"]
    # align names to id_levels so mapping equals original replication logic
    names(modes_vec) <- id_levels
    B_RE[, k] <- unname(modes_vec[id_values])
  }
  
  # combine
  dY <- rowSums(X * B + X_RE * B_RE)
  return(dY)
}
