### Code:
##' @title Calculate predictions for time-depend covariates using INLA
#' @useDynLib HIDeM
#' @importFrom JMbayes2 jm
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  


JMidm<-function(timeVar,
                functional_forms, 
                truncated,
                formLong,
                formSurv,
                dataSurv,
                dataLongi,
                id,
                n_iter,
                n_burnin,
                n_thin,
                n_chain,
                nproc,t0,t1,t2,t3,idm,idd,
                clustertype,lightmode){
  

  
  print("Start running joint univarite models")
  
  modelY<-list()
  length(modelY)<-length(formLong)
  

  for(indice in 1:length(formLong)){

    # need to have all elements of joint
    # global variables otherwise error in predict
    print(paste0("For marker: ",names(functional_forms)[[indice]]))

    JMmodel<-JMbayes2::jm(formSurv, list(formLong[[indice]]), time_var = timeVar, 
                             functional_forms = functional_forms[[indice]],
                             n_iter =n_iter, n_burnin = n_burnin, n_thin =n_thin,
                             n_chains=n_chain, data_Surv = dataSurv,
                             cores=nproc,save_random_effects=T)
  
    browser()
  #keep only necessary info#
    if(lightmode==T){
    JMmodel<-list(model_info=JMmodel$model_info,
                     mcmc=JMmodel$mcmc,
                     control=JMmodel$control,
                     statistics=JMmodel$statistics)
    }
    modelY[[indice]]<-JMmodel
  }
  
  print("End of running univarite models")
  return(modelY)
}
