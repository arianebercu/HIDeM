### Code:
##' @title Calculate predictions for time-depend covariates using INLA
#' @useDynLib HIDeM
#' @importFrom stats D
#' 
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  


## RE NOT OK NEED TO SEE AGAIN
JMidmpredY<-function(timeVar,
                truncated,
                formLong,
                dataSurv,
                dataLongi,
                id,
                t0,t1,t2,t3,ctime,
                assoc,
                modelY,seed,BLUP,scale.X){
  
  # define timePoints of prediction : 
  
  set.seed(seed = seed)

    
  idsubjects<-unique(dataSurv[,colnames(dataSurv)%in%id])
    
    
  timePointsdata<-do.call(rbind, lapply(idsubjects,FUN=function(x){
    index<-which(dataSurv[seq(1,dim(dataSurv)[1],by=2),colnames(dataSurv)%in%id]==x)
    max.int<-ifelse(ctime[index]%in%c(2:5),t2[index],t3[index])
    min.int<-t1[index]
    
    timePoints<-gauss_kronrod_points(lower.intdouble=min.int,
                                     upper.intdouble=max.int,
                                     end.time=t3[index],
                                     truncated=truncated,
                                     entry.time=t0[index])
    return(data.frame(index=x,timePoints=timePoints
                      ))}))
  
  N<-length(unique(dataLongi[,id]))
  newdataLongi<-timePointsdata
  colnames(newdataLongi)<-c(id,timeVar)
  NtimePoints<-ifelse(truncated==F,256,271)
  
  
  Yall<-list()
  length(Yall)<-length(formLong)
  
  if(scale.X==T){
    tcenter<-ifelse(truncated==T,0,min(t0))
    dataCenter<-data.frame(ID=idsubjects,time=tcenter)
    colnames(dataCenter)<-c(id,timeVar)
  }

  for(indice in 1:length(formLong)){


    JMmodel<-modelY$modelY[[indice]]
  
    terms_FE<-JMmodel$model_info$terms$terms_FE[[1]]
    terms_RE<-JMmodel$model_info$terms$terms_RE[[1]]
      
    if(scale.X==F){
      X<-model.matrix(reformulate(attr(terms(terms_FE), "term.labels")),data=newdataLongi)
      Z <- model.matrix(reformulate(attr(terms(terms_RE), "term.labels")),data=newdataLongi)
      
      
      # derivatives of Y(t) design matrix -- does not handle splines so far
      dX <- make_dX(formula(formLong[[indice]]$call$fixed), X=X, timevar = timeVar,data=newdataLongi)
      bars<-findbars(formLong[[indice]]$call$random)
      dZ <- make_dX(reformulate(deparse(bars[[1]][[2]])), X=Z, timevar = timeVar,data=newdataLongi)
    }else{
      X<-model.matrix(reformulate(attr(terms(terms_FE), "term.labels")),data=newdataLongi)
      Z <- model.matrix(reformulate(attr(terms(terms_RE), "term.labels")),data=newdataLongi)
      
      # derivatives of Y(t) design matrix -- does not handle splines so far
      dX <- make_dX(formula(formLong[[indice]]$call$fixed), X=X, timevar = timeVar,data=newdataLongi)
      bars<-findbars(formLong[[indice]]$call$random)
      dZ <- make_dX(reformulate(deparse(bars[[1]][[2]])), X=Z, timevar = timeVar,data=newdataLongi)
      
      X0<-model.matrix(reformulate(attr(terms(terms_FE), "term.labels")),data=dataCenter)
      X0<-do.call(rbind,rep(list(X0,X0),NtimePoints/2))
      Z0 <- model.matrix(reformulate(attr(terms(terms_RE), "term.labels")),data=dataCenter)
      Z0<-do.call(rbind,rep(list(Z0,Z0),NtimePoints/2))
      dX0 <- make_dX(formula(formLong[[indice]]$call$fixed), X=X0, timevar = timeVar,data=dataCenter)
      dZ0 <- make_dX(reformulate(deparse(bars[[1]][[2]])), X=Z0, timevar = timeVar,data=dataCenter)
      
      X[,colnames(X)%in%colnames(X0)]<-X[,colnames(X)%in%colnames(X0)]-X0
      Z[,colnames(Z)%in%colnames(Z0)]<-Z[,colnames(Z)%in%colnames(Z0)]-Z0
      dX[,colnames(dX)%in%colnames(dX0)]<-dX[,colnames(dX)%in%colnames(dX0)]-dX0
      dZ[,colnames(dZ)%in%colnames(dZ0)]<-dZ[,colnames(dZ)%in%colnames(dZ0)]-dZ0
    }
      
      betas<-do.call(rbind,JMmodel$mcmc$betas1)
      
      d1<-dim(JMmodel$mcmc$b[[1]])[3]
      d2<-JMmodel$control$n_chains
      b_mat <- array(0.0, dim = c(dim(JMmodel$mcmc$b[[1]])[1:2], d1*d2))
      b_mat[, , seq(1, d1)]<-JMmodel$mcmc$b[[1]]
      for( i in 2:JMmodel$control$n_chains){
        b_mat[, , seq(d1*(i-1)+1, d1*i)] <- JMmodel$mcmc$b[[i]]
      }
      
      # K <- length(ind_RE)
      # M <- dim(b_mat)[3L]
      M<-dim(betas)[1]
      
      nn<-rep(NtimePoints,N)
      ends   <- cumsum(nn)
      starts <- ends - nn + 1
      choiceY<-na.omit(unlist(assoc[[indice]]))

      if(BLUP==F){
   
      idNsample<-sample(x=c(1:M),size=1)

      if("value"%in% choiceY & "slope"%in%choiceY & "RE"%in%choiceY){
      Fixed <- X %*% t(betas[idNsample, , drop = FALSE])   # (n x m)
      slopeFixed<-dX %*% t(betas[idNsample, , drop = FALSE])
      
      # Terme aléatoire, on remplit directement un tableau vide
      Random_all <- matrix(0, nrow = nrow(Z), ncol = 1)
      Random<-matrix(0,nrow=nrow(Z),ncol=1*dim(b_mat[1, , 1, drop = FALSE])[2])
      slopeRandom_all <- matrix(0, nrow = nrow(Z), ncol = 1)
      
      for (j in seq_len(N)) {
        
        rows <- starts[j]:ends[j]
        Zj   <- Z[rows, , drop = FALSE]                        # (nn[j] x q)
        Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
        
        Random_all[rows, ] <- Zj %*% Bj[1,,]
        Random[rows,]<-matrix(rep(t(Bj[1,,]),NtimePoints),ncol=dim(Bj)[2]*dim(Bj)[3],
                                                nrow=NtimePoints, byrow = TRUE)
        
        dZj   <- dZ[rows, , drop = FALSE]                
        slopeRandom_all[rows, ] <- dZj %*% Bj[1,,]
      }
     
      nameRandom<-colnames(dZ)
      if("(Intercept)"%in%nameRandom){nameRandom[which(nameRandom=="(Intercept)")]<-"Intercept"}
      colnames(Random)<-paste0("RE_",nameRandom,"_",names(functional_forms)[[indice]])
      
      PredYx <- Fixed + Random_all
      slopePredYx<-slopeFixed + slopeRandom_all
      REPredYx<-Random
      
      Outcome<-names(functional_forms)[[indice]]
      slopeOutcome<-paste0("slope_",names(functional_forms)[[indice]])
      REOutcome<-rep(colnames(Random), each = nrow(Random))
      
      PredYx<-cbind(newdataLongi,Outcome=Outcome,PredYx)
      slopePredYx<-cbind(newdataLongi,Outcome=slopeOutcome,slopePredYx)
      REPredYx<-cbind(do.call(rbind, replicate(dim(Random)[2], newdataLongi, simplify = FALSE)),
                      Outcome=REOutcome,
                      as.vector(REPredYx))
      
      colnames(PredYx)[4]<-paste0("Sample_1")
      colnames(slopePredYx)[4]<-paste0("Sample_1")
      colnames(REPredYx)[4]<-paste0("Sample_1")
      
      Yall[[indice]]<- rbind(PredYx,slopePredYx,REPredYx)
      
      }
      
      if("value"%in% choiceY & "slope"%in%choiceY & !"RE"%in%choiceY){
        Fixed <- X %*% t(betas[idNsample, , drop = FALSE])   # (n x m)
        slopeFixed<-dX %*% t(betas[idNsample, , drop = FALSE])
        
        # Terme aléatoire, on remplit directement un tableau vide
        Random_all <- matrix(0, nrow = nrow(Z), ncol = 1)
        
        slopeRandom_all <- matrix(0, nrow = nrow(Z), ncol = 1)
        
        for (j in seq_len(N)) {
          
          rows <- starts[j]:ends[j]
          Zj   <- Z[rows, , drop = FALSE]                        # (nn[j] x q)
          Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
          
          Random_all[rows, ] <- Zj %*% Bj[1,,]
          
          dZj   <- dZ[rows, , drop = FALSE]                
          slopeRandom_all[rows, ] <- dZj %*% Bj[1,,]
        }
        
        PredYx <- Fixed + Random_all
        slopePredYx<-slopeFixed + slopeRandom_all
        
        Outcome<-names(functional_forms)[[indice]]
        slopeOutcome<-paste0("slope_",names(functional_forms)[[indice]])
        
        PredYx<-cbind(newdataLongi,Outcome=Outcome,PredYx)
        slopePredYx<-cbind(newdataLongi,Outcome=slopeOutcome,slopePredYx)
        
        colnames(PredYx)[4]<-"Sample_1"
        colnames(slopePredYx)[4]<-"Sample_1"
        
        Yall[[indice]]<- rbind(PredYx,slopePredYx)
        
      }
      
      if("value"%in% choiceY & !"slope"%in%choiceY & "RE"%in%choiceY){
        Fixed <- X %*% t(betas[idNsample, , drop = FALSE])   # (n x m)
        
        # Terme aléatoire, on remplit directement un tableau vide
        Random_all <- matrix(0, nrow = nrow(Z), ncol = 1)
        Random<-matrix(0,nrow=nrow(Z),ncol=1*dim(b_mat[1, , 1, drop = FALSE])[2])
        
        for (j in seq_len(N)) {
          
          rows <- starts[j]:ends[j]
          Zj   <- Z[rows, , drop = FALSE]                        # (nn[j] x q)
          Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
          
          Random_all[rows, ] <- Zj %*% Bj[1,,]
          Random[rows,]<-matrix(rep(t(Bj[1,,]),NtimePoints),ncol=dim(Bj)[2]*dim(Bj)[3],
                                nrow=NtimePoints, byrow = TRUE)
          
        }
        nameRandom<-colnames(dZ)
        if("(Intercept)"%in%nameRandom){nameRandom[which(nameRandom=="(Intercept)")]<-"Intercept"}
        colnames(Random)<-paste0("RE_",nameRandom,"_",names(functional_forms)[[indice]])
        
        PredYx <- Fixed + Random_all
        REPredYx<-Random
        
        Outcome<-names(functional_forms)[[indice]]
        REOutcome<-rep(colnames(Random), each = nrow(Random))
        
        PredYx<-cbind(newdataLongi,Outcome=Outcome,PredYx)
        REPredYx<-cbind(do.call(rbind, replicate(dim(Random)[2], newdataLongi, simplify = FALSE)),
                        Outcome=REOutcome,
                        as.vector(REPredYx))
        
        colnames(PredYx)[4]<-"Sample_1"
        colnames(REPredYx)[4]<-"Sample_1"
        
        Yall[[indice]]<- rbind(PredYx,REPredYx)
        
      }
      
      if("value"%in% choiceY & !"slope"%in%choiceY & !"RE"%in%choiceY){
        Fixed <- X %*% t(betas[idNsample, , drop = FALSE])   # (n x m)
        
        # Terme aléatoire, on remplit directement un tableau vide
        Random_all <- matrix(0, nrow = nrow(Z), ncol = 1)
        
        for (j in seq_len(N)) {
          
          rows <- starts[j]:ends[j]
          Zj   <- Z[rows, , drop = FALSE]                        # (nn[j] x q)
          Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
          
          Random_all[rows, ] <- Zj %*% Bj[1,,]
         
        }
        
        PredYx <- Fixed + Random_all
        
        Outcome<-names(functional_forms)[[indice]]
        
        PredYx<-cbind(newdataLongi,Outcome=Outcome,PredYx)
        
        colnames(PredYx)[4]<-"Sample_1"
        
        Yall[[indice]]<- PredYx
        
      }

      if(!"value"%in% choiceY & "slope"%in%choiceY & "RE"%in%choiceY){
        
        slopeFixed<-dX %*% t(betas[idNsample, , drop = FALSE])
        
        # Terme aléatoire, on remplit directement un tableau vide
        Random<-matrix(0,nrow=nrow(Z),ncol=1*dim(b_mat[1, , 1, drop = FALSE])[2])
        
        slopeRandom_all <- matrix(0, nrow = nrow(Z), ncol = 1)
        
        for (j in seq_len(N)) {
          
          rows <- starts[j]:ends[j]
          Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
          
          Random[rows,]<-matrix(rep(t(Bj[1,,]),NtimePoints),ncol=dim(Bj)[2]*dim(Bj)[3],
                                nrow=NtimePoints, byrow = TRUE)
          
          dZj   <- dZ[rows, , drop = FALSE]                
          slopeRandom_all[rows, ] <- dZj %*% Bj[1,,]
        }
        nameRandom<-colnames(dZ)
        if("(Intercept)"%in%nameRandom){nameRandom[which(nameRandom=="(Intercept)")]<-"Intercept"}
        colnames(Random)<-paste0("RE_",nameRandom,"_",names(functional_forms)[[indice]])
        
        slopePredYx<-slopeFixed + slopeRandom_all
        REPredYx<-Random
        
        slopeOutcome<-paste0("slope_",names(functional_forms)[[indice]])
        REOutcome<-rep(colnames(Random), each = nrow(Random))
        
        slopePredYx<-cbind(newdataLongi,Outcome=slopeOutcome,slopePredYx)
        REPredYx<-cbind(do.call(rbind, replicate(dim(Random)[2], newdataLongi, simplify = FALSE)),
                        Outcome=REOutcome,
                        as.vector(REPredYx))
        
        colnames(slopePredYx)[4]<-"Sample_1"
        colnames(REPredYx)[4]<-"Sample_1"
        
        Yall[[indice]]<- rbind(slopePredYx,REPredYx)
        
      }
      
      if(!"value"%in% choiceY & "slope"%in%choiceY & !"RE"%in%choiceY){
        slopeFixed<-dX %*% t(betas[idNsample, , drop = FALSE])
        
        for (j in seq_len(N)) {
          
          rows <- starts[j]:ends[j]
          Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
          
          dZj   <- dZ[rows, , drop = FALSE]                
          slopeRandom_all[rows, ] <- dZj %*% Bj[1,,]
        }
        slopePredYx<-slopeFixed + slopeRandom_all
        slopeOutcome<-paste0("slope_",names(functional_forms)[[indice]])
        slopePredYx<-cbind(newdataLongi,Outcome=slopeOutcome,slopePredYx)
        colnames(slopePredYx)[4]<-"Sample_1"
        
        Yall[[indice]]<- slopePredYx
        
      }
      
      if(!"value"%in% choiceY & !"slope"%in%choiceY & "RE"%in%choiceY){
       
        # Terme aléatoire, on remplit directement un tableau vide
        Random<-matrix(0,nrow=nrow(Z),ncol=1*dim(b_mat[1, , 1, drop = FALSE])[2])
        
        for (j in seq_len(N)) {
          
          rows <- starts[j]:ends[j]
          Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
          
          Random[rows,]<-matrix(rep(t(Bj[1,,]),NtimePoints),ncol=dim(Bj)[2]*dim(Bj)[3],
                                nrow=NtimePoints, byrow = TRUE)
          
        }
        nameRandom<-colnames(dZ)
        if("(Intercept)"%in%nameRandom){nameRandom[which(nameRandom=="(Intercept)")]<-"Intercept"}
        colnames(Random)<-paste0("RE_",nameRandom,"_",names(functional_forms)[[indice]])
        
        REPredYx<-Random
        REOutcome<-rep(colnames(Random), each = nrow(Random))
        REPredYx<-cbind(do.call(rbind, replicate(dim(Random)[2], newdataLongi, simplify = FALSE)),
                        Outcome=REOutcome,
                        as.vector(REPredYx))
        colnames(REPredYx)[4]<-"Sample_1"
        
        Yall[[indice]]<- REPredYx
        
      }
      }else{
        REPredYx<-PredYx<-slopePredYx<-NULL
      
        if("value"%in% choiceY & "slope"%in%choiceY & "RE"%in%choiceY){
        # Terme aléatoire, on remplit directement un tableau vide
        Random_mean <- matrix(0,nrow=nrow(Z),ncol=1)
        slopeRandom_mean <- matrix(0,nrow=nrow(Z),ncol=1)
        Random<-matrix(0,nrow=nrow(Z),ncol=length(JMmodel$statistics$Mean$b[1,]))
        for (j in seq_len(N)) {
          rows <- starts[j]:ends[j]
          Zj   <- Z[rows, , drop = FALSE]
          Random_mean[rows, ] <- Zj %*% JMmodel$statistics$Mean$b[j,]
          Random[rows,]<- JMmodel$statistics$Mean$b[j,]
          
          dZj   <- dZ[rows, , drop = FALSE]
          slopeRandom_mean[rows, ] <- dZj %*% JMmodel$statistics$Mean$b[j,]
        }
    
        nameRandom<-names(JMmodel$statistics$Mean$b[j,])
        if("(Intercept)"%in%nameRandom){nameRandom[which(nameRandom=="(Intercept)")]<-"Intercept"}
        colnames(Random)<-paste0("RE_",nameRandom,"_",names(functional_forms)[[indice]])
        
        PredYmean<-X%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean
        slopePredYmean<-dX%*%as.matrix(JMmodel$statistics$Mean$betas1) + slopeRandom_mean
        REPredYx<- Random
        
        Outcome<-names(functional_forms)[[indice]]
        slopeOutcome<-paste0("slope_",names(functional_forms)[[indice]])
        REOutcome<-rep(colnames(Random), each = nrow(Random))
        
        PredYx<-cbind(newdataLongi,Outcome=Outcome,PredYmean)
        slopePredYx<-cbind(newdataLongi,Outcome=slopeOutcome,slopePredYmean)
        REPredYx<-cbind(do.call(rbind, replicate(dim(Random)[2], newdataLongi, simplify = FALSE)),
                        Outcome=REOutcome,
                        as.vector(REPredYx))
        
        colnames(PredYx)[4]<-"Sample_1"
        colnames(slopePredYx)[4]<-"Sample_1"
        colnames(REPredYx)[4]<-"Sample_1"
        
        Yall[[indice]]<- rbind(PredYx,slopePredYx,REPredYx)
        }
        
        if("value"%in% choiceY & "slope"%in%choiceY & !"RE"%in%choiceY){
          # Terme aléatoire, on remplit directement un tableau vide
          Random_mean <- matrix(0,nrow=nrow(Z),ncol=1)
          slopeRandom_mean <- matrix(0,nrow=nrow(Z),ncol=1)
          
          for (j in seq_len(N)) {
            rows <- starts[j]:ends[j]
            Zj   <- Z[rows, , drop = FALSE]
            Random_mean[rows, ] <- Zj %*% JMmodel$statistics$Mean$b[j,]
            
            dZj   <- dZ[rows, , drop = FALSE]
            slopeRandom_mean[rows, ] <- dZj %*% JMmodel$statistics$Mean$b[j,]
          }
          
          
          PredYmean<-X%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean
          slopePredYmean<-dX%*%as.matrix(JMmodel$statistics$Mean$betas1) + slopeRandom_mean
          
          
          Outcome<-names(functional_forms)[[indice]]
          slopeOutcome<-paste0("slope_",names(functional_forms)[[indice]])
          
          PredYx<-cbind(newdataLongi,Outcome=Outcome,PredYmean)
          slopePredYx<-cbind(newdataLongi,Outcome=slopeOutcome,slopePredYmean)
          
          colnames(PredYx)[4]<-"Sample_1"
          colnames(slopePredYx)[4]<-"Sample_1"
          
          Yall[[indice]]<- rbind(PredYx,slopePredYx)
        }
        
        if("value"%in% choiceY & !"slope"%in%choiceY & "RE"%in%choiceY){
          # Terme aléatoire, on remplit directement un tableau vide
          Random_mean <- matrix(0,nrow=nrow(Z),ncol=1)
          Random<-matrix(0,nrow=nrow(Z),ncol=length(JMmodel$statistics$Mean$b[1,]))
          for (j in seq_len(N)) {
            rows <- starts[j]:ends[j]
            Zj   <- Z[rows, , drop = FALSE]
            Random_mean[rows, ] <- Zj %*% JMmodel$statistics$Mean$b[j,]
            Random[rows,]<-JMmodel$statistics$Mean$b[j,]
            
            
          }
          
          nameRandom<-names(JMmodel$statistics$Mean$b[j,])
          if("(Intercept)"%in%nameRandom){nameRandom[which(nameRandom=="(Intercept)")]<-"Intercept"}
          colnames(Random)<-paste0("RE_",nameRandom,"_",names(functional_forms)[[indice]])
          
          PredYmean<-X%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean
          REPredYx<- Random
          
          Outcome<-names(functional_forms)[[indice]]
          REOutcome<-rep(colnames(Random), each = nrow(Random))
          
          PredYx<-cbind(newdataLongi,Outcome=Outcome,PredYmean)
          REPredYx<-cbind(do.call(rbind, replicate(dim(Random)[2], newdataLongi, simplify = FALSE)),
                          Outcome=REOutcome,
                          as.vector(REPredYx))
          
          colnames(PredYx)[4]<-"Sample_1"
          colnames(REPredYx)[4]<-"Sample_1"
          
          Yall[[indice]]<- rbind(PredYx,REPredYx)
        }
        
        if("value"%in% choiceY & !"slope"%in%choiceY & !"RE"%in%choiceY){
          # Terme aléatoire, on remplit directement un tableau vide
          Random_mean <- matrix(0,nrow=nrow(Z),ncol=1)
         
          for (j in seq_len(N)) {
            rows <- starts[j]:ends[j]
            Zj   <- Z[rows, , drop = FALSE]
            Random_mean[rows, ] <- Zj %*% JMmodel$statistics$Mean$b[j,]
            
          }
          
         
          PredYmean<-X%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean
          
          Outcome<-names(functional_forms)[[indice]]
          
          PredYx<-cbind(newdataLongi,Outcome=Outcome,PredYmean)
          
          colnames(PredYx)[4]<-"Sample_1"
          
          Yall[[indice]]<- PredYx
        }
        
        if(!"value"%in% choiceY & "slope"%in%choiceY & "RE"%in%choiceY){
          # Terme aléatoire, on remplit directement un tableau vide
          Random<-matrix(0,nrow=nrow(Z),ncol=length(JMmodel$statistics$Mean$b[1,]))
          slopeRandom_mean <- matrix(0,nrow=nrow(Z),ncol=1)
          
          for (j in seq_len(N)) {
            rows <- starts[j]:ends[j]
            Zj   <- Z[rows, , drop = FALSE]
            Random[rows,]<- JMmodel$statistics$Mean$b[j,]
            
            dZj   <- dZ[rows, , drop = FALSE]
            slopeRandom_mean[rows, ] <- dZj %*% JMmodel$statistics$Mean$b[j,]
          }
          
          nameRandom<-names(JMmodel$statistics$Mean$b[j,])
          if("(Intercept)"%in%nameRandom){nameRandom[which(nameRandom=="(Intercept)")]<-"Intercept"}
          colnames(Random)<-paste0("RE_",nameRandom,"_",names(functional_forms)[[indice]])
          
          slopePredYmean<-dX%*%as.matrix(JMmodel$statistics$Mean$betas1) + slopeRandom_mean
          REPredYx<- Random
          
          slopeOutcome<-paste0("slope_",names(functional_forms)[[indice]])
          REOutcome<-rep(colnames(Random), each = nrow(Random))
          
          slopePredYx<-cbind(newdataLongi,Outcome=slopeOutcome,slopePredYmean)
          REPredYx<-cbind(do.call(rbind, replicate(dim(Random)[2], newdataLongi, simplify = FALSE)),
                          Outcome=REOutcome,
                          as.vector(REPredYx))
          
          colnames(slopePredYx)[4]<-"Sample_1"
          colnames(REPredYx)[4]<-"Sample_1"
          
          Yall[[indice]]<- rbind(slopePredYx,REPredYx)
        }
        
        if(!"value"%in% choiceY & "slope"%in%choiceY & !"RE"%in%choiceY){
          # Terme aléatoire, on remplit directement un tableau vide
         
          slopeRandom_mean <- matrix(0,nrow=nrow(Z),ncol=1)
          
          for (j in seq_len(N)) {
            rows <- starts[j]:ends[j]
            dZj   <- dZ[rows, , drop = FALSE]
            slopeRandom_mean[rows, ] <- dZj %*% JMmodel$statistics$Mean$b[j,]
          }
          
          slopePredYmean<-dX%*%as.matrix(JMmodel$statistics$Mean$betas1) + slopeRandom_mean
          slopeOutcome<-paste0("slope_",names(functional_forms)[[indice]])
          slopePredYx<-cbind(newdataLongi,Outcome=slopeOutcome,slopePredYmean)
          colnames(slopePredYx)[4]<-"Sample_1"
          
          
          Yall[[indice]]<- slopePredYx
        }
        
        if(!"value"%in% choiceY & !"slope"%in%choiceY & "RE"%in%choiceY){
          # Terme aléatoire, on remplit directement un tableau vide
          Random<-matrix(0,nrow=nrow(Z),ncol=length(JMmodel$statistics$Mean$b[1,]))
          for (j in seq_len(N)) {
            rows <- starts[j]:ends[j]
            Zj   <- Z[rows, , drop = FALSE]
            Random[rows,]<- JMmodel$statistics$Mean$b[j,]
            }
          
          nameRandom<-names(JMmodel$statistics$Mean$b[j,])
          if("(Intercept)"%in%nameRandom){nameRandom[which(nameRandom=="(Intercept)")]<-"Intercept"}
          colnames(Random)<-paste0("RE_",nameRandom,"_",names(functional_forms)[[indice]])
          
          REPredYx<- Random
          REOutcome<-rep(colnames(Random), each = nrow(Random))
          REPredYx<-cbind(do.call(rbind, replicate(dim(Random)[2], newdataLongi, simplify = FALSE)),
                          Outcome=REOutcome,
                          as.vector(REPredYx))
          colnames(REPredYx)[4]<-"Sample_1"
          
          Yall[[indice]]<- REPredYx
        }
        
      }
      
      
     
  }
  
  Yall<-do.call(rbind,Yall)
  return(Yall[,colnames(Yall)%in%c(id,timeVar,"Outcome","Sample_1")])
  
}




make_dX <- function(formula, X, timevar, data,use_splines = FALSE, ...) {


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
      dexpr <- stats::D(expr, timevar)
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

