### Code:
##' @title Calculate predictions for time-depend covariates using INLA
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  


JMidm<-function(timeVar,
                functional_forms, 
                truncated,
                formLong,
                formSurv,
                dataSurv,
                dataLongi,
                id,
                Nsample,
                n_iter,
                n_burnin,
                n_thin,
                n_chain,
                nproc,t0,t1,t2,t3,idm,idd,
                clustertype,
                Ypredmethod,
                NtimesPoints,seed){
  
  # define timePoints of prediction : 
  
  set.seed(seed = seed)
  

  if(Ypredmethod=="gauss"){
    
    idsubjects<-unique(dataSurv[,colnames(dataSurv)%in%id])
    
    
  timePointsdata<-do.call(rbind, lapply(idsubjects,FUN=function(x){
    index<-which(dataSurv[seq(1,dim(dataSurv)[1],by=2),colnames(dataSurv)%in%id]==x)
    max.int<-ifelse(idm[index]==1,t2[index],t3[index])
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
  NtimesPoints<-ifelse(truncated==F,256,271)
  
  }else{
    
    times<-seq(min(t0),max(dataLongi[,colnames(dataLongi)%in%timeVar]),length.out=NtimesPoints)
    
    N<-length(unique(dataLongi[,id]))
    
    newdataLongi<-do.call(rbind,lapply(unique(dataLongi[,id]),FUN=function(x){
      data.frame(id=rep(x,NtimesPoints),
                 time=times)}))
    colnames(newdataLongi)<-c(id,timeVar)
  }
  
  print("Start running joint univarite models")
  
  Yall<-list()
  length(Yall)<-length(formLong)
  

  for(indice in 1:length(formLong)){

    # need to have all elements of joint
    # global variables otherwise error in predict
    print(paste0("For marker: ",names(functional_forms)[[indice]]))

    JMmodel<-JMbayes2::jm(formSurv, list(formLong[[indice]]), time_var = timeVar, 
                             functional_forms = functional_forms[[indice]],
                             n_iter =n_iter, n_burnin = n_burnin, n_thin =n_thin,
                             n_chains=n_chain, data_Surv = dataSurv,
                             cores=nproc,save_random_effects=T)
  
      
      
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
      
      nn<-rep(NtimesPoints,N)
      ends   <- cumsum(nn)
      starts <- ends - nn + 1
      
      if(Nsample>1){
        
      idNsample<-sample(x=c(1:M),size=Nsample-1)

      Fixed <- X %*% t(betas[idNsample, , drop = FALSE])   # (n x m)
      slopeFixed<-dX %*% t(betas[idNsample, , drop = FALSE])
      
      # Terme aléatoire, on remplit directement un tableau vide
      Random_all <- matrix(0, nrow = nrow(Z), ncol = Nsample-1)
      Random_mean <- matrix(0,nrow=nrow(Z),ncol=1)
      
      slopeRandom_all <- matrix(0, nrow = nrow(Z), ncol = Nsample-1)
      slopeRandom_mean <- matrix(0,nrow=nrow(Z),ncol=1)
      
      for (j in seq_len(N)) {
        rows <- starts[j]:ends[j]
        Zj   <- Z[rows, , drop = FALSE]                        # (nn[j] x q)
        Bj   <- b_mat[j, , idNsample, drop = FALSE] # (q x m)
        Random_all[rows, ] <- Zj %*% Bj[1,,]
        Random_mean[rows, ] <- Zj %*% JMmodel$statistics$Mean$b[j,]
        
        dZj   <- dZ[rows, , drop = FALSE]                
        slopeRandom_all[rows, ] <- dZj %*% Bj[1,,]
        slopeRandom_mean[rows, ] <- dZj %*% JMmodel$statistics$Mean$b[j,]
      }
      
      PredYx <- Fixed + Random_all
      slopePredYx<-slopeFixed + slopeRandom_all

      }else{
        PredYx<-slopePredYx<-NULL
        
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
  
      }
      
      PredYmean<-X%*%as.matrix(JMmodel$statistics$Mean$betas1) + Random_mean
      slopePredYmean<-dX%*%as.matrix(JMmodel$statistics$Mean$betas1) + slopeRandom_mean
      
      PredYx<-cbind(PredYmean,PredYx)
      slopePredYx<-cbind(slopePredYmean,slopePredYx)

      Outcome<-names(functional_forms)[[indice]]
      slopeOutcome<-paste0("slope_",names(functional_forms)[[indice]])
      
      PredYx<-cbind(newdataLongi,Outcome,PredYx)
      slopePredYx<-cbind(newdataLongi,slopeOutcome,slopePredYx)
      
      colnames(PredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
      colnames(slopePredYx)[4:(Nsample+3)]<-paste0("Sample_",c(1:Nsample))
      
    Yall[[indice]]<- rbind(PredYx,slopePredYx)
  }
  
  Yall<-do.call(rbind,Yall)
  print("End of running joint univarite models")
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

