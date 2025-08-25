### Code:
##' @title Calculate predictions for time-depend covariates using INLA
#' @useDynLib HIDeM
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  


INLAidm<-function(timeVar,family,basRisk,assoc,
                  truncated,formLong,formSurv,dataSurv,dataLongi,id,
                  NsampleHY,NsampleFE,NsampleRE,nproc,t0,t1,t2,t3,
                  idm,idd,clustertype,Ypredmethod,NtimesPoints){
  
  # define timePoints of prediction : 


  idsubjects<-unique(dataSurv[,colnames(dataSurv)%in%id])
  
  timePointsdata<-do.call(rbind, lapply(idsubjects,FUN=function(x){
    index<-which(dataSurv[,colnames(dataSurv)%in%id]==x)
    max.int<-ifelse(idm[index]==1,t2[index],t3[index])
    min.int<-t1[index]
  
    
    timePoints<-gauss_kronrod_points(lower.intdouble=min.int,
                                     upper.intdouble=max.int,
                                     end.time=t3[index],
                                     truncated=truncated,
                                     entry.time=t0[index])
    return(data.frame(timePoints=timePoints,
                      index=x))}))

    
  print("Start running joint univarite models")

  Yall<-list()
  length(Yall)<-length(formLong)
  
  formSurvinla<<-formSurv
  dataLongiinla<<-dataLongi
  dataSurvinla<<-dataSurv
  idinla<<-id
  timeVarinla<<-timeVar
  

    for(indice in 1:length(formLong)){

      # need to have all elements of joint
      # global variables otherwise error in predict
      formLonginla<<-formLong[[indice]]
      familyinla<<-family[indice]
      basRiskinla<<-basRisk[indice]
      associnla<<-assoc[[indice]]

      INLAmodel<-tryCatch({ INLAjoint::joint(formSurv = formSurvinla,
                                       formLong = formLonginla,
                                       dataLong = dataLongiinla, dataSurv=dataSurvinla, id = idinla, timeVar = timeVarinla,
                                       family = familyinla,
                                       basRisk = basRiskinla, NbasRisk = 15, assoc = associnla,
                                       control=list(int.strategy="eb",lightmode=1))
      }, error = function(e) {
        # Return NULL on error to skip this patient
        NULL
      })
      

      
      if(is.null(INLAmodel)){stop("The inla model for your marker could not be run, see above warnings.")}
      

      # attention cannot do parallel over parallel of predict in inla 
     if(Ypredmethod=="gauss"){
       
      PredYx<-do.call(rbind,lapply(idsubjects,FUN=function(x){
        
        indexLongi<-which(dataLongi[,colnames(dataLongi)%in%id]==x)
        datai<-dataLongi[indexLongi,]
        datai<-merge(x=datai,y=dataSurv,by=id,all.x=T)
        timeinla<<-timePointsdata[timePointsdata$index==x,colnames(timePointsdata)=="timePoints"]

        Pred<-tryCatch({
          
            # error if idm==1 thus : 
            # visit in datai needs to have a timePoints greater than its 
            # maximum value --> add one point and erase it after
            # needs to be in timePoints as well as horizon 
            res<-predict(INLAmodel,
                         newData = datai,
                         NtimesPoints=length(timeinla)+1,
                         timePoints = c(timeinla,max(timeinla)+0.1),
                         startTime=0,
                         horizon = max(timeinla)+0.1,
                         NsampleHY =  NsampleHY,
                         NsampleFE =  NsampleFE,
                         NsampleRE =  NsampleRE,
                         return.samples=TRUE)$PredL
          res<-res[-dim(res)[1],]
          res<-res[match(timeinla,res[,colnames(res)%in%timeVar]),]
          res
        }, error = function(e) {
          # Return NULL on error to skip this patient
          NULL
        })
        return(Pred)
        
      }))
      

     }else{
       
       PredYx<-predict(INLAmodel,
                   newData = dataLongi,
                   NtimesPoints=NtimesPoints,
                   timePoints =seq(min(t0),max(dataLongi[,colnames(dataLongi)%in%timeVar]),length.out=NtimesPoints),
                   startTime=0,
                   horizon = max(t3),
                   NsampleHY =  NsampleHY,
                   NsampleFE =  NsampleFE,
                   NsampleRE =  NsampleRE,
                   return.samples=TRUE)$PredL
       
       
     }
      
      Yall[[indice]]<-PredYx
      
  
      }

  browser()
  Yall<-do.call(rbind,Yall)
  print("End of running joint univarite models")
  Nsample<-NsampleHY*NsampleFE*NsampleRE
  return(Yall[,colnames(Yall)%in%c(id,timeVar,"Outcome",paste0("Sample_",c(1:Nsample)))])
  
}
