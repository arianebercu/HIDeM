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
                Nsample,t0,t1,t2,t3,ctime,
                modelY,seed,BLUP){
  
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
      
      if(BLUP==F){
        
      idNsample<-sample(x=c(1:M),size=Nsample)

      Fixed <- X %*% t(betas[idNsample, , drop = FALSE])   # (n x m)
      slopeFixed<-dX %*% t(betas[idNsample, , drop = FALSE])
      
      # Terme aléatoire, on remplit directement un tableau vide
      Random<-Random_all <- matrix(0, nrow = nrow(Z), ncol = Nsample)
      
      slopeRandom_all <- matrix(0, nrow = nrow(Z), ncol = Nsample)
      
      for (j in seq_len(N)) {
        
        rows <- starts[j]:ends[j]
        Zj   <- Z[rows, , drop = FALSE]                        # (nn[j] x q)
        Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
        
        Random_all[rows, ] <- Zj %*% Bj[1,,]
        Random[rows,]<-Bj[1,,]
        
        dZj   <- dZ[rows, , drop = FALSE]                
        slopeRandom_all[rows, ] <- dZj %*% Bj[1,,]
      }
      colnames(Random)<-paste0("RE_",c(1:dim(Random)[2]),"_",names(functional_forms)[[indice]])
      
      PredYx <- Fixed + Random_all
      slopePredYx<-slopeFixed + slopeRandom_all
      REPredYx<-Random
      
       browser()
      Outcome<-names(functional_forms)[[indice]]
      slopeOutcome<-paste0("slope_",names(functional_forms)[[indice]])
      REOutcome<-paste0("RE_",names(functional_forms)[[indice]])
      
      PredYx<-cbind(newdataLongi,Outcome=Outcome,PredYx)
      slopePredYx<-cbind(newdataLongi,Outcome=slopeOutcome,slopePredYx)
      REPredYx<-cbind(newdataLongi,Outcome=REOutcome,REPredYx)
      
      colnames(PredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
      colnames(slopePredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
      colnames(REPredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
      
      Yall[[indice]]<- rbind(PredYx,slopePredYx,REPredYx)

      }else{
        REPredYx<-PredYx<-slopePredYx<-NULL
        
        # Terme aléatoire, on remplit directement un tableau vide
        Random<-Random_mean <- matrix(0,nrow=nrow(Z),ncol=1)
        slopeRandom_mean <- matrix(0,nrow=nrow(Z),ncol=1)
        
        for (j in seq_len(N)) {
          rows <- starts[j]:ends[j]
          Zj   <- Z[rows, , drop = FALSE]
          Random_mean[rows, ] <- Zj %*% JMmodel$statistics$Mean$b[j,]
          Random[rows,]<-diag(1,nrow=dim(Zj)[1])%*% JMmodel$statistics$Mean$b[j,]
          
          dZj   <- dZ[rows, , drop = FALSE]
          slopeRandom_mean[rows, ] <- dZj %*% JMmodel$statistics$Mean$b[j,]
        }
        
        colnames(Random)<-paste0("RE_",c(1:dim(Random)[2]),"_",names(functional_forms)[[indice]])
        
        PredYmean<-X%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean
        slopePredYmean<-dX%*%as.matrix(JMmodel$statistics$Mean$betas1) + slopeRandom_mean
        REPredYmean<- Random
        
        Outcome<-names(functional_forms)[[indice]]
        slopeOutcome<-paste0("slope_",names(functional_forms)[[indice]])
        REOutcome<-paste0("RE_",names(functional_forms)[[indice]])
        
        PredYx<-cbind(newdataLongi,Outcome=Outcome,PredYmean)
        slopePredYx<-cbind(newdataLongi,Outcome=slopeOutcome,slopePredYmean)
        REPredYx<-cbind(newdataLongi,Outcome=REOutcome,REPredYxmean)
        
        colnames(PredYx)[4]<-"Sample_1"
        colnames(slopePredYx)[4]<-"Sample_1"
        colnames(REPredYx)[4]<-"Sample_1"
        
        Yall[[indice]]<- rbind(PredYx,slopePredYx,REPredYx)
        
      }
      
      
     
  }
  
  Yall<-do.call(rbind,Yall)
  return(Yall[,colnames(Yall)%in%c(id,timeVar,"Outcome",paste0("Sample_",c(1:Nsample)))])
  
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

