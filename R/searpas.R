#' Line search
#' @noRd

searpas <- function(vw,step,b,delta,funcpa,res.out.error,...){
  #cat("dans searpas, b=",b,"\n")
  #cat("              vw=",vw,"\n")
  #cat("              delta=",delta,"\n")
  goto50  <- function(step,vlw2,fi1,fi2,fi3,b,delta,funcpa,...){
    vm <- vlw2-(step*(fi1-fi3))/(2*(fi1-2*fi2+fi3)) 
    fim <- valfpa(vm,b,delta,funcpa,...)
    return(list(vm=vm,fim=fim))
  }
  vlw1 <- log(vw)
  vlw2 <- vlw1+step
  fi1 <- valfpa(vlw1,b,delta,funcpa,...)
  fi2 <- valfpa(vlw2,b,delta,funcpa,...)

  if((sum(!is.finite(fi1)) > 0) || (sum(!is.finite(fi2)) > 0)){
    cat("Probably too much accuracy requested...\n")
    cat("Last step values :\n")
    cat("      b :",res.out.error$old.b,"\n")
    cat("      function value :",res.out.error$old.rl,"\n")
    cat("      Convergence criteria: parameters stability on beta=", res.out.error$old.ca, "\n")
    cat("                          : function stability=", res.out.error$old.cb, "\n") 
    stop("")	
  }
  
  if((fi2 >= fi1)){
    vlw3 <- vlw2
    vlw2 <- vlw1
    fi3 <- fi2
    fi2 <- fi1
    step <- -step
    vlw1 <- vlw2+step
    fi1 <- valfpa(vlw1,b,delta,funcpa,...)
    gt50 <- goto50(step,vlw2,fi1,fi2,fi3,b,delta,funcpa,...)
    vm <- gt50$vm
    fim <- gt50$fim
    if(is.na(fim)) fim <- 10E10
    if(fim <= fi2){
      vw <- exp(vm)
    }else{
      vm <- vlw2
      fim <- fi2
      vw <- exp(vm)
    }
    
  }else{
    vlw <- vlw1
    vlw1 <- vlw2
    vlw2 <- vlw
    fim <- fi1
    fi1 <- fi2
    fi2 <- fim
    
    for(i in 1:40){
      vlw3 <- vlw2
      vlw2 <- vlw1
      fi3 <- fi2
      fi2 <- fi1
      vlw1=vlw2+step
      fi1 <- valfpa(vlw1,b,delta,funcpa,...)
      if(fi1 > fi2){
        gt50 <- goto50(step,vlw2,fi1,fi2,fi3,b,delta,funcpa,...) 
        out <- 1
        break
      }
      if(fi1 == fi2){
        fim <- fi2
        vm <- vlw2
        vw <- exp(vm)
        out <- 1
        break
      }
    }
  }	
  return(list(vw=vw,fim=fim))
}


DYNsearpas_weib <- function(step,b,res.out.error,fistart,s,nva01Y,nva02Y,nva12Y,fix,v,
                            ctime,no,fu,ve01,ve02,ve12,dimnva01,dimnva02,dimnva12,
                            nva01,nva02,nva12,t0,t1,t2,t3,troncature,y01,y02,y12,
                            p01,p02,p12,dimp01,dimp02,dimp12,
                            Ntime,lambda,alpha,penalty.factor,penalty,
                            penalty.weights){
  #cat("dans searpas, b=",b,"\n")
  #cat("              vw=",vw,"\n")
  #cat("              delta=",delta,"\n")
  DYNgoto50  <- function(step,vlw2,fi1,fi2,fi3,b,s,nva01Y,nva02Y,nva12Y,fix,v,
                         ctime,no,fu,ve01,ve02,ve12,dimnva01,dimnva02,dimnva12,
                         nva01,nva02,nva12,t0,t1,t2,t3,troncature,y01,y02,y12,
                         p01,p02,p12,dimp01,dimp02,dimp12,Ntime,lambda,alpha,
                         penalty.factor,penalty,penalty.weights){
    vm <- vlw2-(step*(fi1-fi3))/(2*(fi1-2*fi2+fi3)) 
    if(is.na(vm)|(vm=="NaN")) {vm<--2E9}
    fim <- threshold.beta(beta=b,
                          vw=exp(vm),
                          s=s,
                          nva01Y=nva01Y,
                          nva02Y=nva02Y,
                          nva12Y=nva12Y,
                          fix=fix,
                          v=v,
                          ctime=ctime,
                          no=no,
                          fu=fu,
                          ve01=ve01,
                          ve02=ve02,
                          ve12=ve12,
                          dimnva01=dimnva01,
                          dimnva02=dimnva02,
                          dimnva12=dimnva12,
                          nva01=nva01,
                          nva02=nva02,
                          nva12=nva12,
                          t0=t0,
                          t1=t1,
                          t2=t2,
                          t3=t3,
                          troncature=troncature,
                          y01=y01,
                          y02=y02,
                          y12=y12,
                          p01=p01,
                          p02=p02,
                          p12=p12,
                          dimp01=dimp01,
                          dimp02=dimp02,
                          dimp12=dimp12,
                          Ntime=Ntime,
                          lambda=lambda,
                          alpha=alpha,
                          penalty.factor=penalty.factor,
                          penalty=penalty,
                          penalty.weights=penalty.weights)
    bim<-fim$b
    fim<-(-fim$res)
    return(list(vm=vm,fim=fim,bim=bim))
  }
  

  vlw1 <- 0
  vlw2 <- vlw1+step
  fi1 <- threshold.beta(beta=b,
                        vw=exp(vlw1),
                        s=s,
                        nva01Y=nva01Y,
                        nva02Y=nva02Y,
                        nva12Y=nva12Y,
                        fix=fix,
                        v=v,
                        ctime=ctime,
                        no=no,
                        fu=fu,
                        ve01=ve01,
                        ve02=ve02,
                        ve12=ve12,
                        dimnva01=dimnva01,
                        dimnva02=dimnva02,
                        dimnva12=dimnva12,
                        nva01=nva01,
                        nva02=nva02,
                        nva12=nva12,
                        t0=t0,
                        t1=t1,
                        t2=t2,
                        t3=t3,
                        troncature=troncature,
                        y01=y01,
                        y02=y02,
                        y12=y12,
                        p01=p01,
                        p02=p02,
                        p12=p12,
                        dimp01=dimp01,
                        dimp02=dimp02,
                        dimp12=dimp12,
                        Ntime=Ntime,
                        lambda=lambda,
                        alpha=alpha,
                        penalty.factor=penalty.factor,
                        penalty=penalty,
                        penalty.weights=penalty.weights)
  bi1<-fi1$b
  fi1<-(-fi1$res)
  if((fi1<(-fistart)) & !(fi1%in%c(-1e9,1e9))){
    return(list(vw=exp(vlw1),fim=fi1,bim=bi1))
  }

  fi2 <- threshold.beta(beta=b,
                        vw=exp(vlw2),
                        s=s,
                        nva01Y=nva01Y,
                        nva02Y=nva02Y,
                        nva12Y=nva12Y,
                        fix=fix,
                        v=v,
                        ctime=ctime,
                        no=no,
                        fu=fu,
                        ve01=ve01,
                        ve02=ve02,
                        ve12=ve12,
                        dimnva01=dimnva01,
                        dimnva02=dimnva02,
                        dimnva12=dimnva12,
                        nva01=nva01,
                        nva02=nva02,
                        nva12=nva12,
                        t0=t0,
                        t1=t1,
                        t2=t2,
                        t3=t3,
                        troncature=troncature,
                        y01=y01,
                        y02=y02,
                        y12=y12,
                        p01=p01,
                        p02=p02,
                        p12=p12,
                        dimp01=dimp01,
                        dimp02=dimp02,
                        dimp12=dimp12,
                        Ntime=Ntime,
                        lambda=lambda,
                        alpha=alpha,
                        penalty.factor=penalty.factor,
                        penalty=penalty,
                        penalty.weights=penalty.weights)
  bi2<-fi2$b
  fi2<-(-fi2$res)
  
  if((sum(!is.finite(fi1)) > 0) || (sum(!is.finite(fi2)) > 0)){
    cat("Probably too much accuracy requested...\n")
    cat("Last step values :\n")
    cat("      b :",res.out.error$old.b,"\n")
    cat("      function value :",res.out.error$old.rl,"\n")
    cat("      Convergence criteria: parameters stability on beta=", res.out.error$old.ca, "\n")
    cat("                          : function stability=", res.out.error$old.cb, "\n") 
    stop("")	
  }
  
  if(is.na(fi2) | fi2%in%c(1e9,-1e9)){fi2<-1e9}
  if(is.na(fi1) | fi1%in%c(1e9,-1e9)){fi1<-1e9}
  if((fi2 >= fi1)){
    vlw3 <- vlw2
    vlw2 <- vlw1
    fi3 <- fi2
    bi3<-bi2
    fi2 <- fi1
    bi2<-bi1
    step <- -step
    vlw1 <- vlw2+step
    fi1 <- threshold.beta(beta=b,
                          vw=exp(vlw1),
                          s=s,
                          nva01Y=nva01Y,
                          nva02Y=nva02Y,
                          nva12Y=nva12Y,
                          fix=fix,
                          v=v,
                          ctime=ctime,
                          no=no,
                          fu=fu,
                          ve01=ve01,
                          ve02=ve02,
                          ve12=ve12,
                          dimnva01=dimnva01,
                          dimnva02=dimnva02,
                          dimnva12=dimnva12,
                          nva01=nva01,
                          nva02=nva02,
                          nva12=nva12,
                          t0=t0,
                          t1=t1,
                          t2=t2,
                          t3=t3,
                          troncature=troncature,
                          y01=y01,
                          y02=y02,
                          y12=y12,
                          p01=p01,
                          p02=p02,
                          p12=p12,
                          dimp01=dimp01,
                          dimp02=dimp02,
                          dimp12=dimp12,
                          Ntime=Ntime,
                          lambda=lambda,
                          alpha=alpha,
                          penalty.factor=penalty.factor,
                          penalty=penalty,
                          penalty.weights=penalty.weights)
    bi1<-fi1$b
    fi1<-(-fi1$res)
 
    gt50 <- DYNgoto50(step=step,
                      vlw2=vlw2,
                      fi1=fi1,
                      fi2=fi2,
                      fi3=fi3,
                      b=b,
                      s=s,
                      nva01Y=nva01Y,
                      nva02Y=nva02Y,
                      nva12Y=nva12Y,
                      fix=fix,
                      v=v,
                      ctime=ctime,
                      no=no,
                      fu=fu,
                      ve01=ve01,
                      ve02=ve02,
                      ve12=ve12,
                      dimnva01=dimnva01,
                      dimnva02=dimnva02,
                      dimnva12=dimnva12,
                      nva01=nva01,
                      nva02=nva02,
                      nva12=nva12,
                      t0=t0,
                      t1=t1,
                      t2=t2,
                      t3=t3,
                      troncature=troncature,
                      y01=y01,
                      y02=y02,
                      y12=y12,
                      p01=p01,
                      p02=p02,
                      p12=p12,
                      dimp01=dimp01,
                      dimp02=dimp02,
                      dimp12=dimp12,
                      Ntime=Ntime,
                      lambda=lambda,
                      alpha=alpha,
                      penalty.factor=penalty.factor,
                      penalty=penalty,
                      penalty.weights=penalty.weights)
    vm <- gt50$vm
    fim <- gt50$fim
    bim<-gt50$bim
    if(is.na(fim)) fim <- 1e9
    if(fim <= fi2){
      vw <- exp(vm)
    }else{
      vm <- vlw2
      fim <- fi2
      bim<-bi2
      vw <- exp(vm)
    }
    
  }else{
    vlw <- vlw1
    vlw1 <- vlw2
    vlw2 <- vlw
    fim <- fi1
    bim<-bi1
    fi1 <- fi2
    bi1<-bi2
    fi2 <- fim
    bi2<-bim
    
    for(i in 1:40){
      vlw3 <- vlw2
      vlw2 <- vlw1
      fi3 <- fi2
      bi3<-bi2
      fi2 <- fi1
      bi2<-bi1
      vlw1=vlw2+step
      fi1 <- threshold.beta(beta=b,
                            vw=exp(vlw1),
                            s=s,
                            nva01Y=nva01Y,
                            nva02Y=nva02Y,
                            nva12Y=nva12Y,
                            fix=fix,
                            v=v,
                            ctime=ctime,
                            no=no,
                            fu=fu,
                            ve01=ve01,
                            ve02=ve02,
                            ve12=ve12,
                            dimnva01=dimnva01,
                            dimnva02=dimnva02,
                            dimnva12=dimnva12,
                            nva01=nva01,
                            nva02=nva02,
                            nva12=nva12,
                            t0=t0,
                            t1=t1,
                            t2=t2,
                            t3=t3,
                            troncature=troncature,
                            y01=y01,
                            y02=y02,
                            y12=y12,
                            p01=p01,
                            p02=p02,
                            p12=p12,
                            dimp01=dimp01,
                            dimp02=dimp02,
                            dimp12=dimp12,
                            Ntime=Ntime,
                            lambda=lambda,
                            alpha=alpha,
                            penalty.factor=penalty.factor,
                            penalty=penalty,
                            penalty.weights=penalty.weights)
      bi1<-fi1$b
      fi1<-(-fi1$res)
      if(fi1 > fi2){
        gt50 <- DYNgoto50(step=step,
                          vlw2=vlw2,
                          fi1=fi1,
                          fi2=fi2,
                          fi3=fi3,
                          b=b,
                          s=s,
                          nva01Y=nva01Y,
                          nva02Y=nva02Y,
                          nva12Y=nva12Y,
                          fix=fix,
                          v=v,
                          ctime=ctime,
                          no=no,
                          fu=fu,
                          ve01=ve01,
                          ve02=ve02,
                          ve12=ve12,
                          dimnva01=dimnva01,
                          dimnva02=dimnva02,
                          dimnva12=dimnva12,
                          nva01=nva01,
                          nva02=nva02,
                          nva12=nva12,
                          t0=t0,
                          t1=t1,
                          t2=t2,
                          t3=t3,
                          troncature=troncature,
                          y01=y01,
                          y02=y02,
                          y12=y12,
                          p01=p01,
                          p02=p02,
                          p12=p12,
                          dimp01=dimp01,
                          dimp02=dimp02,
                          dimp12=dimp12,
                          Ntime=Ntime,
                          lambda=lambda,
                          alpha=alpha,
                          penalty.factor=penalty.factor,
                          penalty=penalty,
                          penalty.weights=penalty.weights)
        
        out <- 1
        vw<-exp(gt50$vm)
        bim<-gt50$bim
        fim<-(-gt50$fim)
        break
      }
      if(fi1 == fi2){
        fim <- fi2
        bim<-bi2
        vm <- vlw2
        vw <- exp(vm)
        out <- 1
        break
      }
    }
  }	
  return(list(vw=vw,fim=-fim,bim=bim))
}
