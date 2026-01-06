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

  #browser()
  for(indice in 1:length(formLong)){


    JMmodel<-modelY$modelY[[indice]]
  
    terms_FE<-JMmodel$model_info$terms$terms_FE[[1]]
    terms_RE<-JMmodel$model_info$terms$terms_RE[[1]]
      
  
      X<-model.matrix(reformulate(attr(terms(terms_FE), "term.labels")),data=newdataLongi)
      Z <- model.matrix(reformulate(attr(terms(terms_RE), "term.labels")),data=newdataLongi)
      
      
      # derivatives of Y(t) design matrix -- does not handle splines so far
      dX <- make_dX(formula(formLong[[indice]]$call$fixed), X=X, timevar = timeVar,data=newdataLongi)
      bars<-findbars(formLong[[indice]]$call$random)
      dZ <- make_dX(reformulate(deparse(bars[[1]][[2]])), X=Z, timevar = timeVar,data=newdataLongi)
   
      if(scale.X==T){
      X0<-model.matrix(reformulate(attr(terms(terms_FE), "term.labels")),data=dataCenter)
      X0<-do.call(rbind,rep(list(X0,X0),NtimePoints/2))
      Z0 <- model.matrix(reformulate(attr(terms(terms_RE), "term.labels")),data=dataCenter)
      Z0<-do.call(rbind,rep(list(Z0,Z0),NtimePoints/2))
      dX0 <- make_dX(formula(formLong[[indice]]$call$fixed), X=X0, timevar = timeVar,data=dataCenter)
      dZ0 <- make_dX(reformulate(deparse(bars[[1]][[2]])), X=Z0, timevar = timeVar,data=dataCenter)
      }else{
        X0<-matrix(0,ncol=ncol(X),nrow=nrow(X))
        dX0<-matrix(0,ncol=ncol(dX),nrow=nrow(dX))
        Z0<-matrix(0,ncol=ncol(Z),nrow=nrow(Z))
        dZ0<-matrix(0,ncol=ncol(dZ),nrow=nrow(dZ))
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
      
      Fixed0 <- X0 %*% t(betas[idNsample, , drop = FALSE])   # (n x m)
      slopeFixed0<-dX0 %*% t(betas[idNsample, , drop = FALSE])
      
      # Terme aléatoire, on remplit directement un tableau vide
      Random_all0<-Random_all <- matrix(0, nrow = nrow(Z), ncol = 1)
      Random<-matrix(0,nrow=nrow(Z),ncol=1*dim(b_mat[1, , 1, drop = FALSE])[2])
      slopeRandom_all0<-slopeRandom_all <- matrix(0, nrow = nrow(Z), ncol = 1)
      
      for (j in seq_len(N)) {
        
        rows <- starts[j]:ends[j]
        Zj   <- Z[rows, , drop = FALSE]                        # (nn[j] x q)
        Z0j   <- Z0[rows, , drop = FALSE]  
        Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
        
        Random_all[rows, ] <- Zj %*% Bj[1,,]
        Random_all0[rows, ] <- Z0j %*% Bj[1,,]
        
        Random[rows,]<-matrix(rep(t(Bj[1,,]),NtimePoints),ncol=dim(Bj)[2]*dim(Bj)[3],
                                                nrow=NtimePoints, byrow = TRUE)
        
        dZj   <- dZ[rows, , drop = FALSE] 
        dZ0j   <- dZ0[rows, , drop = FALSE] 
        slopeRandom_all[rows, ] <- dZj %*% Bj[1,,]
        slopeRandom_all0[rows, ] <- dZ0j %*% Bj[1,,]
      }
     
      nameRandom<-colnames(dZ)
      if("(Intercept)"%in%nameRandom){nameRandom[which(nameRandom=="(Intercept)")]<-"Intercept"}
      colnames(Random)<-paste0("RE_",nameRandom,"_",as.character(formula(formLong[[indice]])[[2]]))
      
      PredYx <- Fixed + Random_all
      PredYx0<-Fixed0 + Random_all0
      
      slopePredYx<-slopeFixed + slopeRandom_all
      slopePredYx0<-slopeFixed0 + slopeRandom_all0
      
      if(scale.X==T){
        PredYx<-(PredYx-mean(PredYx0))/sd(PredYx0)
        slopePredYx<-(slopePredYx-mean(slopePredYx0))/sd(slopePredYx0)
      }
      
      REPredYx<-Random
      
      Outcome<-as.character(formula(formLong[[indice]])[[2]])
      slopeOutcome<-paste0("slope_",as.character(formula(formLong[[indice]])[[2]]))
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
        
        Fixed0 <- X0 %*% t(betas[idNsample, , drop = FALSE])   # (n x m)
        slopeFixed0<-dX0 %*% t(betas[idNsample, , drop = FALSE])
        
        # Terme aléatoire, on remplit directement un tableau vide
        Random_all0<-Random_all <- matrix(0, nrow = nrow(Z), ncol = 1)
        
        slopeRandom_all0<-slopeRandom_all <- matrix(0, nrow = nrow(Z), ncol = 1)
        for (j in seq_len(N)) {
          
          rows <- starts[j]:ends[j]
          Zj   <- Z[rows, , drop = FALSE]                        # (nn[j] x q)
          Z0j   <- Z0[rows, , drop = FALSE]                        # (nn[j] x q)
          Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
          
          Random_all[rows, ] <- Zj %*% Bj[1,,]
          Random_all0[rows, ] <- Z0j %*% Bj[1,,]
          
          dZj   <- dZ[rows, , drop = FALSE]                
          slopeRandom_all[rows, ] <- dZj %*% Bj[1,,]
          dZ0j   <- dZ0[rows, , drop = FALSE]                
          slopeRandom_all0[rows, ] <- dZ0j %*% Bj[1,,]
        }
        
    
        PredYx <- Fixed + Random_all
        PredYx0<-Fixed0 + Random_all0
        
        slopePredYx<-slopeFixed + slopeRandom_all
        slopePredYx0<-slopeFixed0 + slopeRandom_all0
        
        if(scale.X==T){
          PredYx<-(PredYx-mean(PredYx0))/sd(PredYx0)
          slopePredYx<-(slopePredYx-mean(slopePredYx0))/sd(slopePredYx0)
        }
        
        Outcome<-as.character(formula(formLong[[indice]])[[2]])
        slopeOutcome<-paste0("slope_",as.character(formula(formLong[[indice]])[[2]]))
        
        PredYx<-cbind(newdataLongi,Outcome=Outcome,PredYx)
        slopePredYx<-cbind(newdataLongi,Outcome=slopeOutcome,slopePredYx)
        
        colnames(PredYx)[4]<-"Sample_1"
        colnames(slopePredYx)[4]<-"Sample_1"
        
        Yall[[indice]]<- rbind(PredYx,slopePredYx)
        
      }
      
      if("value"%in% choiceY & !"slope"%in%choiceY & "RE"%in%choiceY){
        Fixed <- X %*% t(betas[idNsample, , drop = FALSE])   # (n x m)
        Fixed0 <- X0 %*% t(betas[idNsample, , drop = FALSE])   # (n x m)
         
        # Terme aléatoire, on remplit directement un tableau vide
        Random_all0<-Random_all <- matrix(0, nrow = nrow(Z), ncol = 1)
        Random<-matrix(0,nrow=nrow(Z),ncol=1*dim(b_mat[1, , 1, drop = FALSE])[2])
        
        for (j in seq_len(N)) {
          
          rows <- starts[j]:ends[j]
          Zj   <- Z[rows, , drop = FALSE]                        # (nn[j] x q)
          Z0j   <- Z0[rows, , drop = FALSE]                        # (nn[j] x q)
          Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
          
          Random_all[rows, ] <- Zj %*% Bj[1,,]
          Random_all0[rows, ] <- Z0j %*% Bj[1,,]
          Random[rows,]<-matrix(rep(t(Bj[1,,]),NtimePoints),ncol=dim(Bj)[2]*dim(Bj)[3],
                                nrow=NtimePoints, byrow = TRUE)
          
        }
        nameRandom<-colnames(dZ)
        if("(Intercept)"%in%nameRandom){nameRandom[which(nameRandom=="(Intercept)")]<-"Intercept"}
        colnames(Random)<-paste0("RE_",nameRandom,"_",as.character(formula(formLong[[indice]])[[2]]))
        
        
        PredYx <- Fixed + Random_all
        PredYx0<-Fixed0 + Random_all0
        
        if(scale.X==T){
          PredYx<-(PredYx-mean(PredYx0))/sd(PredYx0)
          }
        REPredYx<-Random
        
        Outcome<-as.character(formula(formLong[[indice]])[[2]])
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
        
        
        Fixed0 <- X0 %*% t(betas[idNsample, , drop = FALSE])   # (n x m)
        
        # Terme aléatoire, on remplit directement un tableau vide
        Random_all0<-Random_all <- matrix(0, nrow = nrow(Z), ncol = 1)
       
        
        for (j in seq_len(N)) {
          
          rows <- starts[j]:ends[j]
          Zj   <- Z[rows, , drop = FALSE]                        # (nn[j] x q)
          Z0j   <- Z0[rows, , drop = FALSE]                        # (nn[j] x q)
          Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
          
          Random_all[rows, ] <- Zj %*% Bj[1,,]
          Random_all0[rows, ] <- Z0j %*% Bj[1,,]
         
        }
        

        
        PredYx <- Fixed + Random_all
        PredYx0<-Fixed0 + Random_all0
        
       if(scale.X==T){
          PredYx<-(PredYx-mean(PredYx0))/sd(PredYx0)
          }
        
        Outcome<-as.character(formula(formLong[[indice]])[[2]])
        
        PredYx<-cbind(newdataLongi,Outcome=Outcome,PredYx)
        
        colnames(PredYx)[4]<-"Sample_1"
        
        Yall[[indice]]<- PredYx
        
      }

      if(!"value"%in% choiceY & "slope"%in%choiceY & "RE"%in%choiceY){
        
        slopeFixed<-dX %*% t(betas[idNsample, , drop = FALSE])
        
        slopeFixed0<-dX0 %*% t(betas[idNsample, , drop = FALSE])
        
        Random<-matrix(0,nrow=nrow(Z),ncol=1*dim(b_mat[1, , 1, drop = FALSE])[2])
        slopeRandom_all0<-slopeRandom_all <- matrix(0, nrow = nrow(Z), ncol = 1)
        
        for (j in seq_len(N)) {
          
          rows <- starts[j]:ends[j]
          Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
          
          Random[rows,]<-matrix(rep(t(Bj[1,,]),NtimePoints),ncol=dim(Bj)[2]*dim(Bj)[3],
                                nrow=NtimePoints, byrow = TRUE)
          
          dZj   <- dZ[rows, , drop = FALSE]                
          slopeRandom_all[rows, ] <- dZj %*% Bj[1,,]
          dZ0j   <- dZ0[rows, , drop = FALSE]                
          slopeRandom_all0[rows, ] <- dZ0j %*% Bj[1,,]
        }
        nameRandom<-colnames(dZ)
        if("(Intercept)"%in%nameRandom){nameRandom[which(nameRandom=="(Intercept)")]<-"Intercept"}
        colnames(Random)<-paste0("RE_",nameRandom,"_",as.character(formula(formLong[[indice]])[[2]]))
        
        
        slopePredYx<-slopeFixed + slopeRandom_all
        slopePredYx0<-slopeFixed0 + slopeRandom_all0
        
        if(scale.X==T){
          slopePredYx<-(slopePredYx-mean(slopePredYx0))/sd(slopePredYx0)
        }
        REPredYx<-Random
        
        slopeOutcome<-paste0("slope_",as.character(formula(formLong[[indice]])[[2]]))
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
        slopeFixed0<-dX0 %*% t(betas[idNsample, , drop = FALSE])
        
        
        for (j in seq_len(N)) {
          
          rows <- starts[j]:ends[j]
          Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
          
          dZj   <- dZ[rows, , drop = FALSE]                
          slopeRandom_all[rows, ] <- dZj %*% Bj[1,,]
          dZ0j   <- dZ0[rows, , drop = FALSE]                
          slopeRandom_all0[rows, ] <- dZ0j %*% Bj[1,,]
        }
        
        slopePredYx<-slopeFixed + slopeRandom_all
        slopePredYx0<-slopeFixed0 + slopeRandom_all0
        
        if(scale.X==T){
         slopePredYx<-(slopePredYx-mean(slopePredYx0))/sd(slopePredYx0)
        }
        slopeOutcome<-paste0("slope_",as.character(formula(formLong[[indice]])[[2]]))
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
        colnames(Random)<-paste0("RE_",nameRandom,"_",as.character(formula(formLong[[indice]])[[2]]))
        
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
          Random_mean0<-Random_mean <- matrix(0,nrow=nrow(Z),ncol=1)
          slopeRandom_mean0<-slopeRandom_mean <- matrix(0,nrow=nrow(Z),ncol=1)
        Random<-matrix(0,nrow=nrow(Z),ncol=length(JMmodel$statistics$Mean$b[1,]))
        for (j in seq_len(N)) {
          rows <- starts[j]:ends[j]
          Zj   <- Z[rows, , drop = FALSE]
          Z0j   <- Z0[rows, , drop = FALSE]
          Random_mean[rows, ] <- Zj %*% JMmodel$statistics$Mean$b[j,]
          Random_mean0[rows, ] <- Z0j %*% JMmodel$statistics$Mean$b[j,]
          Random[rows,]<- JMmodel$statistics$Mean$b[j,]
          
          dZj   <- dZ[rows, , drop = FALSE]
          slopeRandom_mean[rows, ] <- dZj %*% JMmodel$statistics$Mean$b[j,]
          dZ0j   <- dZ0[rows, , drop = FALSE]
          slopeRandom_mean0[rows, ] <- dZ0j %*% JMmodel$statistics$Mean$b[j,]
        }
    
        nameRandom<-names(JMmodel$statistics$Mean$b[j,])
        if("(Intercept)"%in%nameRandom){nameRandom[which(nameRandom=="(Intercept)")]<-"Intercept"}
        colnames(Random)<-paste0("RE_",nameRandom,"_",as.character(formula(formLong[[indice]])[[2]]))
        
        PredYmean<-X%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean
        slopePredYmean<-dX%*%as.matrix(JMmodel$statistics$Mean$betas1) + slopeRandom_mean
        
        PredYmean0<-X0%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean0
        slopePredYmean0<-dX0%*%as.matrix(JMmodel$statistics$Mean$betas1) + slopeRandom_mean0
        
        if(scale.X==T){
          PredYmean<-(PredYmean-mean(PredYmean0))/sd(PredYmean0)
          slopePredYmean<-(slopePredYmean-mean(slopePredYmean0))/sd(slopePredYmean0)
        }
        REPredYx<- Random
        
        Outcome<-as.character(formula(formLong[[indice]])[[2]])
        slopeOutcome<-paste0("slope_",as.character(formula(formLong[[indice]])[[2]]))
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
          Random_mean0<-Random_mean <- matrix(0,nrow=nrow(Z),ncol=1)
          slopeRandom_mean0<-slopeRandom_mean <- matrix(0,nrow=nrow(Z),ncol=1)
          for (j in seq_len(N)) {
            rows <- starts[j]:ends[j]
            Zj   <- Z[rows, , drop = FALSE]
            Random_mean[rows, ] <- Zj %*% JMmodel$statistics$Mean$b[j,]
            Z0j   <- Z0[rows, , drop = FALSE]
            Random_mean0[rows, ] <- Z0j %*% JMmodel$statistics$Mean$b[j,]
            
            dZj   <- dZ[rows, , drop = FALSE]
            slopeRandom_mean[rows, ] <- dZj %*% JMmodel$statistics$Mean$b[j,]
            dZ0j   <- dZ0[rows, , drop = FALSE]
            slopeRandom_mean0[rows, ] <- dZ0j %*% JMmodel$statistics$Mean$b[j,]
          }
          
          
          PredYmean<-X%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean
          slopePredYmean<-dX%*%as.matrix(JMmodel$statistics$Mean$betas1) + slopeRandom_mean
          
          PredYmean0<-X0%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean0
          slopePredYmean0<-dX0%*%as.matrix(JMmodel$statistics$Mean$betas1) + slopeRandom_mean0
          
          if(scale.X==T){
            PredYmean<-(PredYmean-mean(PredYmean0))/sd(PredYmean0)
            slopePredYmean<-(slopePredYmean-mean(slopePredYmean0))/sd(slopePredYmean0)
          }
          
          
          Outcome<-as.character(formula(formLong[[indice]])[[2]])
          slopeOutcome<-paste0("slope_",as.character(formula(formLong[[indice]])[[2]]))
          
          PredYx<-cbind(newdataLongi,Outcome=Outcome,PredYmean)
          slopePredYx<-cbind(newdataLongi,Outcome=slopeOutcome,slopePredYmean)
          
          colnames(PredYx)[4]<-"Sample_1"
          colnames(slopePredYx)[4]<-"Sample_1"
          
          Yall[[indice]]<- rbind(PredYx,slopePredYx)
        }
        
        if("value"%in% choiceY & !"slope"%in%choiceY & "RE"%in%choiceY){
          # Terme aléatoire, on remplit directement un tableau vide
          Random<-matrix(0,nrow=nrow(Z),ncol=length(JMmodel$statistics$Mean$b[1,]))
          Random_mean0<-Random_mean <- matrix(0,nrow=nrow(Z),ncol=1)
          
          for (j in seq_len(N)) {
            rows <- starts[j]:ends[j]
            Zj   <- Z[rows, , drop = FALSE]
            Random_mean[rows, ] <- Zj %*% JMmodel$statistics$Mean$b[j,]
            Random[rows,]<-JMmodel$statistics$Mean$b[j,]
            Z0j   <- Z0[rows, , drop = FALSE]
            Random_mean0[rows, ] <- Z0j %*% JMmodel$statistics$Mean$b[j,]
            
            
            
          }
          
          nameRandom<-names(JMmodel$statistics$Mean$b[j,])
          if("(Intercept)"%in%nameRandom){nameRandom[which(nameRandom=="(Intercept)")]<-"Intercept"}
          colnames(Random)<-paste0("RE_",nameRandom,"_",as.character(formula(formLong[[indice]])[[2]]))
          
          PredYmean<-X%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean
          PredYmean0<-X0%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean0
         
          
          if(scale.X==T){
            PredYmean<-(PredYmean-mean(PredYmean0))/sd(PredYmean0)
          }
          REPredYx<- Random
          
          Outcome<-as.character(formula(formLong[[indice]])[[2]])
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
          Random_mean0<-Random_mean <- matrix(0,nrow=nrow(Z),ncol=1)
          
          for (j in seq_len(N)) {
            rows <- starts[j]:ends[j]
            Zj   <- Z[rows, , drop = FALSE]
            Random_mean[rows, ] <- Zj %*% JMmodel$statistics$Mean$b[j,]
            Z0j   <- Z0[rows, , drop = FALSE]
            Random_mean0[rows, ] <- Z0j %*% JMmodel$statistics$Mean$b[j,]
            
          }
          
         
          PredYmean<-X%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean
          
          PredYmean0<-X0%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean0
          
          if(scale.X==T){
            PredYmean<-(PredYmean-mean(PredYmean0))/sd(PredYmean0)
          }
          
          Outcome<-as.character(formula(formLong[[indice]])[[2]])
          
          PredYx<-cbind(newdataLongi,Outcome=Outcome,PredYmean)
          
          colnames(PredYx)[4]<-"Sample_1"
          
          Yall[[indice]]<- PredYx
        }
        
        if(!"value"%in% choiceY & "slope"%in%choiceY & "RE"%in%choiceY){
          # Terme aléatoire, on remplit directement un tableau vide
          Random<-matrix(0,nrow=nrow(Z),ncol=length(JMmodel$statistics$Mean$b[1,]))
          slopeRandom_mean0<-slopeRandom_mean <- matrix(0,nrow=nrow(Z),ncol=1)
          for (j in seq_len(N)) {
            rows <- starts[j]:ends[j]
            Random[rows,]<- JMmodel$statistics$Mean$b[j,]
            
            dZj   <- dZ[rows, , drop = FALSE]
            slopeRandom_mean[rows, ] <- dZj %*% JMmodel$statistics$Mean$b[j,]
            dZ0j   <- dZ0[rows, , drop = FALSE]
            slopeRandom_mean0[rows, ] <- dZ0j %*% JMmodel$statistics$Mean$b[j,]
          }
          
          nameRandom<-names(JMmodel$statistics$Mean$b[j,])
          if("(Intercept)"%in%nameRandom){nameRandom[which(nameRandom=="(Intercept)")]<-"Intercept"}
          colnames(Random)<-paste0("RE_",nameRandom,"_",as.character(formula(formLong[[indice]])[[2]]))
          
          slopePredYmean<-dX%*%as.matrix(JMmodel$statistics$Mean$betas1) + slopeRandom_mean
          
         
          slopePredYmean0<-dX0%*%as.matrix(JMmodel$statistics$Mean$betas1) + slopeRandom_mean0
          
          if(scale.X==T){
            slopePredYmean<-(slopePredYmean-mean(slopePredYmean0))/sd(slopePredYmean0)
          }
          REPredYx<- Random
          
          slopeOutcome<-paste0("slope_",as.character(formula(formLong[[indice]])[[2]]))
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
         
          slopeRandom_mean0<-slopeRandom_mean <- matrix(0,nrow=nrow(Z),ncol=1)
          for (j in seq_len(N)) {
            rows <- starts[j]:ends[j]
            dZj   <- dZ[rows, , drop = FALSE]
            slopeRandom_mean[rows, ] <- dZj %*% JMmodel$statistics$Mean$b[j,]
            dZ0j   <- dZ0[rows, , drop = FALSE]
            slopeRandom_mean0[rows, ] <- dZ0j %*% JMmodel$statistics$Mean$b[j,]
          }
          
          slopePredYmean<-dX%*%as.matrix(JMmodel$statistics$Mean$betas1) + slopeRandom_mean
          
          slopePredYmean0<-dX0%*%as.matrix(JMmodel$statistics$Mean$betas1) + slopeRandom_mean0
          
          if(scale.X==T){
            slopePredYmean<-(slopePredYmean-mean(slopePredYmean0))/sd(slopePredYmean0)
          }
          
          slopeOutcome<-paste0("slope_",as.character(formula(formLong[[indice]])[[2]]))
          slopePredYx<-cbind(newdataLongi,Outcome=slopeOutcome,slopePredYmean)
          colnames(slopePredYx)[4]<-"Sample_1"
          
          
          Yall[[indice]]<- slopePredYx
        }
        
        if(!"value"%in% choiceY & !"slope"%in%choiceY & "RE"%in%choiceY){
          # Terme aléatoire, on remplit directement un tableau vide
          Random<-matrix(0,nrow=nrow(Z),ncol=length(JMmodel$statistics$Mean$b[1,]))
          for (j in seq_len(N)) {
            rows <- starts[j]:ends[j]
            Random[rows,]<- JMmodel$statistics$Mean$b[j,]
            }
          
          nameRandom<-names(JMmodel$statistics$Mean$b[j,])
          if("(Intercept)"%in%nameRandom){nameRandom[which(nameRandom=="(Intercept)")]<-"Intercept"}
          colnames(Random)<-paste0("RE_",nameRandom,"_",as.character(formula(formLong[[indice]])[[2]]))
          
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

