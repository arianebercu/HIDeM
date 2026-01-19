


DYNidm.penalty.weib.nonproc<-function(beta.start,
                                    s.start,
                                    npm01,
                                    npm02,
                                    npm12,
                                    npm01Y,
                                    npm02Y,
                                    npm12Y,
                                    npm,npmweib,
                                    fix0,fix00,fix0.beta,
                                    size_V,epsa,epsb,epsd,eps.eigen,maxiter,maxiter.pena,
                                    ctime,N,
                                    ve01,ve02,ve12,dimnva01,dimnva02,dimnva12,nvat01,nvat02,nvat12,
                                    t0,t1,t2,t3,troncature,nlambda,lambda,
                                    alpha,penalty.factor,penalty,partialH,
                                    Nsample,
                                    NtimePoints,
                                    p01,p02,p12,
                                    dimp01,dimp02,dimp12,
                                    defpositive,warmstart,
                                    y01k,y02k,y12k,min){
  
 
  pbr_compu<-0
  # combine model 
  combine_lambda_mla<-function(x,newx){
    
    
      list(b=cbind(x$b,newx$b),
           V=cbind(x$V,newx$V),
           H=cbind(x$H,newx$H),
           fix=cbind(x$fix,newx$fix),
           lambda=cbind(x$lambda,newx$lambda),
           alpha=c(x$alpha,newx$alpha),
           fn.value=c(x$fn.value,newx$fn.value),
           fn.value.pena=c(x$fn.value.pena,newx$fn.value.pena),
           ni=c(x$ni,newx$ni),
           istop=c(x$istop,newx$istop),
           ca.beta=cbind(x$ca.beta,newx$ca.beta),
           ca.spline=cbind(x$ca.spline,newx$ca.spline),
           ca.validity=cbind(x$ca.validity,newx$ca.validity),
           cb=cbind(x$cb,newx$cb))
    
  }
  
  
  id.lambda<-NULL # for cran check 
  V0<-NA
  
  if(warmstart==F){
    if(partialH==F){
      outputNsample<-foreach::foreach(id.lambda=1:nlambda,
                                      .combine = combine_lambda_mla,
                                      .errorhandling = "remove")%do%{
                                        
                                        # computation pbr 
                                        
                                        pbr_compu<-0
                                        
                                        beta<-beta.start
                                        s<-s.start
                                        
                                        
                                        converged<-F
                                        ite<-0
                                        # if beta not change do not need to recalculate weights 
                                        H<-T
                                        
                                        eval.cv.spline<-rep(NA,maxiter+1)
                                        eval.cv.beta<-rep(NA,maxiter+1)
                                        eval.cv.loglik<-rep(NA,maxiter+1)
                                        eval.loglik<-rep(NA,maxiter+1)
                                        eval.validity<-rep(NA,maxiter+1)
                                        
                                        
                                        
                                        
                                        while(converged==F & ite<=maxiter){
                                          
                                          
                                          b<-c(s,beta)
                                          bfix<-b[fix0==1]
                                          b<-b[fix0==0]
                                          # derivative of loglik
                                          
                                          output<-DYNderivaweib( h=1e-8,
                                                                 b=b,
                                                                 npm=npm,
                                                                 npar=size_V,
                                                                 bfix=bfix,
                                                                 fix=fix0,
                                                                 ctime=ctime,
                                                                 no=N,
                                                                 ve01=ve01,
                                                                 ve02=ve02,
                                                                 ve12=ve12,
                                                                 dimnva01=dimnva01,
                                                                 dimnva02=dimnva02,
                                                                 dimnva12=dimnva12,
                                                                 nva01=nvat01,
                                                                 nva02=nvat02,
                                                                 nva12=nvat12,
                                                                 t0=t0,
                                                                 t1=t1,
                                                                 t2=t2,
                                                                 t3=t3,
                                                                 troncature=troncature,
                                                                 y01=y01k,
                                                                 y02=y02k,
                                                                 y12=y12k,
                                                                 p01=p01,
                                                                 p02=p02,
                                                                 p12=p12,
                                                                 dimp01=dimp01,
                                                                 dimp02=dimp02,
                                                                 dimp12=dimp12,
                                                                 Ntime=NtimePoints)
                                          output<-output$v
                                          
                                          if(ite==0){# loglik penalised
                                            fn.value<-gaussDYNidmlLikelihoodweibpena(b=b,
                                                                                     npm=npm,
                                                                                     npar=size_V,
                                                                                     bfix=bfix,
                                                                                     fix=fix0,
                                                                                     ctime=ctime,
                                                                                     no=N,
                                                                                     ve01=ve01,
                                                                                     ve02=ve02,
                                                                                     ve12=ve12,
                                                                                     dimnva01=dimnva01,
                                                                                     dimnva02=dimnva02,
                                                                                     dimnva12=dimnva12,
                                                                                     nva01=nvat01,
                                                                                     nva02=nvat02,
                                                                                     nva12=nvat12,
                                                                                     t0=t0,
                                                                                     t1=t1,
                                                                                     t2=t2,
                                                                                     t3=t3,
                                                                                     troncature=troncature,
                                                                                     y01=y01k,
                                                                                     y02=y02k,
                                                                                     y12=y12k,
                                                                                     p01=p01,
                                                                                     p02=p02,
                                                                                     p12=p12,
                                                                                     dimp01=dimp01,
                                                                                     dimp02=dimp02,
                                                                                     dimp12=dimp12,
                                                                                     Ntime=NtimePoints,
                                                                                     
                                                                                     lambda=lambda[id.lambda,],
                                                                                     alpha=alpha,
                                                                                     penalty.factor=penalty.factor,
                                                                                     penalty=penalty)
                                          }
                                          
                                          if(any(is.na(output))|any(output==Inf) |any(output==-Inf)){
                                            warning("Computational error for calculation of the hessian : division by 0 or Infinite value")
                                            if(ite==0){
                                              
                                              
                                              fu<-output[(min+1):length(output)]
                                              V<-matrix(0,nrow=npm,ncol=npm)
                                              V[upper.tri(V,diag=T)]<-output[1:min]
                                              V<-V+t(V)
                                              diag(V)<-diag(V)/2
                                              # deriva gives information matrix
                                              tr <- sum(diag(V))/npm
                                              V0<-V}
                                            ite<-ite+1
                                            pbr_compu<-1
                                            break
                                          }
                                          
                                          
                                          fu<-output[(min+1):length(output)]
                                          V<-matrix(0,nrow=npm,ncol=npm)
                                          V[upper.tri(V,diag=T)]<-output[1:min]
                                          V<-V+t(V)
                                          diag(V)<-diag(V)/2
                                          
                                          # deriva gives information matrix
                                          tr <- sum(diag(V))/npm
                                          V0<-V
                                          
                                          eigen.values<-eigen(V,symmetric=T,only.values=T)$values
                                          
                                          if(defpositive==T){
                                            idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                          }else{
                                            idpos<-ifelse(any(diag(V)==0),1,0)
                                          }
                                          
                                          
                                          idpos0<-idpos
                                          
                                          ncount<-da<-ga<-0
                                          
                                          while(idpos != 0){
                                            
                                            if(ncount==0){ 
                                              ga <- 0.01
                                              da <- 1E-2
                                            }else{
                                              if(((ncount <= 3) | (ga >= 1)) ){
                                                da <- da * 5
                                              }else{# if ncount > 10 only update ga 
                                                ga <- ga * 5
                                                # do not put ga at 1 as no countmax otherwise infinite while 
                                                if(ga > 1) ga <- 1
                                              }
                                            }
                                            
                                            ncount <- ncount + 1
                                            
                                            diagV <- diag(V)
                                            # put abs (1-ga) better than 1-ga cause ga can now be >1
                                            diagV<-ifelse(diagV!=0,diagV+da*(abs((1.e0-ga))*abs(diagV)+ga*tr),
                                                          da*ga*tr)
                                            
                                            diag(V)<-diagV
                                            # if we have a convex log-vraisemblance in eta then :
                                            # all eigen  values of the hessienne are >0.
                                            
                                            if(sum(V==Inf)>0|sum(V==-Inf)>0){break}
                                            eigen.values<-eigen(V,symmetric=T,only.values=T)$values
                                            # check if hessienne defined positive
                                            
                                            
                                            if(defpositive==T){
                                              idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                            }else{
                                              idpos<-ifelse(any(diag(V)==0),1,0)
                                            }
                                            
                                            # if(def.positive==T){
                                            #   idpos<-ifelse(any(eigen.values<=0),1,0)
                                            # }else{idpos<-ifelse(any(abs(eigen.values)==0),1,0)}
                                            
                                          }
                                          
                                          
                                          if(idpos!=0){
                                            
                                            warning("Hessian not defined positive")
                                            pbr_compu<-2
                                            ite<-ite+1
                                            break
                                          }
                                          
                                          
                                          # update for beta 
                                          output.cv<-DYNcv.model(beta=beta,
                                                                 nva01=npm01,
                                                                 nva02=npm02,
                                                                 nva12=npm12,
                                                                 nva01Y=npm01Y,
                                                                 nva02Y=npm02Y,
                                                                 nva12Y=npm12Y,
                                                                 fix=fix0[7:size_V],
                                                                 penalty.factor=penalty.factor,
                                                                 penalty=penalty,
                                                                 v=V,
                                                                 fu=fu,
                                                                 lambda=lambda[id.lambda,],
                                                                 alpha=alpha
                                          )
                                          
                                          # verify validity of parameters update 
                                          # and that we are better than previous estimates 
                                          
                                          b<-c(s,output.cv$b)
                                          
                                          betanew<-b[(6+1):size_V]
                                          
                                          # penalised loglik see if inferior to previous
                                          res<-gaussDYNidmlLikelihoodweibpena(b=b,
                                                                              npm=size_V,
                                                                              npar=size_V,
                                                                              bfix=1,
                                                                              fix=rep(0,size_V),
                                                                              ctime=ctime,
                                                                              no=N,
                                                                              ve01=ve01,
                                                                              ve02=ve02,
                                                                              ve12=ve12,
                                                                              dimnva01=dimnva01,
                                                                              dimnva02=dimnva02,
                                                                              dimnva12=dimnva12,
                                                                              nva01=nvat01,
                                                                              nva02=nvat02,
                                                                              nva12=nvat12,
                                                                              t0=t0,
                                                                              t1=t1,
                                                                              t2=t2,
                                                                              t3=t3,
                                                                              troncature=troncature,
                                                                              y01=y01k,
                                                                              y02=y02k,
                                                                              y12=y12k,
                                                                              p01=p01,
                                                                              p02=p02,
                                                                              p12=p12,
                                                                              dimp01=dimp01,
                                                                              dimp02=dimp02,
                                                                              dimp12=dimp12,
                                                                              Ntime=NtimePoints,
                                                                              lambda=lambda[id.lambda,],
                                                                              alpha=alpha,
                                                                              penalty.factor=penalty.factor,
                                                                              penalty=penalty)
                                          
                                          
                                          # we want to maximise the loglik thus : 
                                          # we have issue if res is NA or if not higher than previous one 
                                          # if not better or do not exist need to readjust
                                          # value of beta 
                                          if(res %in%c(-1e9,1e9) | res < fn.value){
                                            
                                            th<-1e-5
                                            step<-log(1.5)
                                            delta<-output.cv$b-c(beta)
                                            
                                            maxt <- max(abs(delta)) 
                                            
                                            if(maxt == 0){
                                              vw <- th
                                            }else{
                                              vw <- th/maxt
                                            }
                                            if(ite>0){
                                              res.out.error <- list("old.b"=round(c(s,beta)),
                                                                    "old.rl"=round(fn.value),
                                                                    "old.ca"=round(eval.cv.beta[ite]),
                                                                    "old.cb"=round(eval.cv.loglik[ite]))
                                            }else{
                                              res.out.error <- list("old.b"=round(c(s,beta)),
                                                                    "old.rl"=round(fn.value),
                                                                    "old.ca"=round(1),
                                                                    "old.cb"=round(1))
                                            }
                                            # from mla package 
                                            sears<-searpas(vw=vw,
                                                           step=step,
                                                           b=beta,
                                                           delta=delta,
                                                           funcpa=gaussDYNidmlLikelihoodweibpena,
                                                           res.out.error=res.out.error,
                                                           npm=npm,
                                                           npar=size_V,
                                                           bfix=s,
                                                           fix=fix0,
                                                           ctime=ctime,
                                                           no=N,
                                                           ve01=ve01,
                                                           ve02=ve02,
                                                           ve12=ve12,
                                                           dimnva01=dimnva01,
                                                           dimnva02=dimnva02,
                                                           dimnva12=dimnva12,
                                                           nva01=nvat01,
                                                           nva02=nvat02,
                                                           nva12=nvat12,
                                                           t0=t0,
                                                           t1=t1,
                                                           t2=t2,
                                                           t3=t3,
                                                           troncature=troncature,
                                                           y01=y01k,
                                                           y02=y02k,
                                                           y12=y12k,
                                                           p01=p01,
                                                           p02=p02,
                                                           p12=p12,
                                                           dimp01=dimp01,
                                                           dimp02=dimp02,
                                                           dimp12=dimp12,
                                                           Ntime=NtimePoints,
                                                           lambda=lambda[id.lambda,],
                                                           alpha=alpha,
                                                           penalty.factor=penalty.factor,
                                                           penalty=penalty)
                                            
                                            
                                            betanew<-beta+delta*sears$vw
                                            betanew<-ifelse(abs(betanew)<=0.0001,0,betanew)
                                            b<-c(s,betanew)
                                            
                                            
                                            res<-gaussDYNidmlLikelihoodweibpena(b=b,
                                                                                npm=size_V,
                                                                                npar=size_V,
                                                                                bfix=1,
                                                                                fix=rep(0,size_V),
                                                                                ctime=ctime,
                                                                                no=N,
                                                                                ve01=ve01,
                                                                                ve02=ve02,
                                                                                ve12=ve12,
                                                                                dimnva01=dimnva01,
                                                                                dimnva02=dimnva02,
                                                                                dimnva12=dimnva12,
                                                                                nva01=nvat01,
                                                                                nva02=nvat02,
                                                                                nva12=nvat12,
                                                                                t0=t0,
                                                                                t1=t1,
                                                                                t2=t2,
                                                                                t3=t3,
                                                                                troncature=troncature,
                                                                                y01=y01k,
                                                                                y02=y02k,
                                                                                y12=y12k,
                                                                                p01=p01,
                                                                                p02=p02,
                                                                                p12=p12,
                                                                                dimp01=dimp01,
                                                                                dimp02=dimp02,
                                                                                dimp12=dimp12,
                                                                                Ntime=NtimePoints,
                                                                                lambda=lambda[id.lambda,],
                                                                                alpha=alpha,
                                                                                penalty.factor=penalty.factor,
                                                                                penalty=penalty)
                                          }
                                          # if not better or do not exist need to readjust
                                          # value of beta 
                                          if(res %in%c(-1e9,1e9) | any(is.infinite(c(s,betanew)))){
                                            
                                            ite<-ite+1
                                            validity<-F
                                            eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                            eval.validity[ite]<-validity
                                            # save output over iterations 
                                            
                                            pbr_compu<-3
                                            break
                                          }else{validity<-T}
                                          
                                          
                                          # betanew already include s
                                          b<-c(s,betanew)
                                          
                                          bfix<-b[fix0.beta==1]
                                          b<-b[fix0.beta==0]
                                          # update for modelPar
                                          if(npmweib!=0){
                                            output.mla<- marqLevAlg::mla(b=b,
                                                                         fn=gaussDYNidmlLikelihoodweib,
                                                                         gr=DYNreggrmlaweibana,
                                                                         epsa=epsa,
                                                                         epsb=epsb,
                                                                         epsd=epsd,
                                                                         maxiter=maxiter.pena,
                                                                         minimize=F,
                                                                         npm= npmweib,
                                                                         npar=size_V,
                                                                         bfix=bfix,
                                                                         fix=fix0.beta,
                                                                         ctime=ctime,
                                                                         no=N,
                                                                         ve01=ve01,
                                                                         ve02=ve02,
                                                                         ve12=ve12,
                                                                         dimnva01=dimnva01,
                                                                         dimnva02=dimnva02,
                                                                         dimnva12=dimnva12,
                                                                         nva01=nvat01,
                                                                         nva02=nvat02,
                                                                         nva12=nvat12,
                                                                         t0=t0,
                                                                         t1=t1,
                                                                         t2=t2,
                                                                         t3=t3,
                                                                         troncature=troncature,
                                                                         y01=y01k,
                                                                         y02=y02k,
                                                                         y12=y12k,
                                                                         p01=p01,
                                                                         p02=p02,
                                                                         p12=p12,
                                                                         dimp01=dimp01,
                                                                         dimp02=dimp02,
                                                                         dimp12=dimp12,
                                                                         Ntime=NtimePoints)
                                            
                                            # look at convergence for each lambda :
                                            # mla output is loglik
                                            # new values for splines:
                                            snew<-s
                                            snew[fix00[1:6]==0]<-output.mla$b
                                            
                                            
                                            
                                            if(nvat01>0){
                                              b01<-betanew[1:nvat01][penalty.factor[1:nvat01]==1]
                                              if(p01>0){
                                                b01<-c(b01,betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1])
                                              }
                                            }else{
                                              if(p01>0){
                                                b01<-betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1]
                                              }else{
                                                b01<-0
                                              }
                                            }
                                            
                                            if(nvat02>0){
                                              b02<-betanew[(nvat01+1):(nvat01+nvat02)][penalty.factor[(nvat01+1):(nvat01+nvat02)]==1]
                                              if(p02>0){
                                                b02<-c(b02,betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1])
                                              }
                                            }else{
                                              if(p02>0){
                                                b02<-betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1]
                                              }else{b02<-0}
                                            }
                                            
                                            if(nvat12>0){
                                              b12<-betanew[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)][penalty.factor[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)]==1]
                                              if(p12>0){
                                                b12<-c(b12,betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1])
                                              }
                                            }else{
                                              if(p12>0){
                                                b12<-betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1]
                                              }else{b12<-0}
                                            }
                                            
                                            # calculate loglik pen 
                                            if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
                                              fn.valuenew<-output.mla$fn.value-lambda[id.lambda,1]*alpha*sum(abs(b01))-lambda[id.lambda,1]*(1-alpha)*sum(b01*b01)
                                              fn.valuenew<-fn.valuenew-lambda[id.lambda,2]*alpha*sum(abs(b02))-lambda[id.lambda,2]*(1-alpha)*sum(b02*b02)
                                              fn.valuenew<-fn.valuenew-lambda[id.lambda,3]*alpha*sum(abs(b12))-lambda[id.lambda,3]*(1-alpha)*sum(b12*b12)
                                            }
                                            
                                            
                                            if(penalty=="mcp"){
                                              
                                              p01<-rep(alpha*lambda[id.lambda,1]*lambda[id.lambda,1]/2,length(b01))
                                              idbeta<-which(b01<=alpha*lambda[id.lambda,1])
                                              p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])-((b01[idbeta]*b01[idbeta])/2*alpha)
                                              
                                              p02<-rep(alpha*lambda[id.lambda,2]*lambda[id.lambda,2]/2,length(b02))
                                              idbeta<-which(b02<=alpha*lambda[id.lambda,2])
                                              p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])-((b02[idbeta]*b02[idbeta])/2*alpha)
                                              
                                              p12<-rep(alpha*lambda[id.lambda,3]*lambda[id.lambda,3]/2,length(b12))
                                              idbeta<-which(b12<=alpha*lambda[id.lambda,3])
                                              p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])-((b12[idbeta]*b12[idbeta])/2*alpha)
                                              
                                              fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                              
                                            }
                                            
                                            if(penalty=="scad"){
                                              
                                              p01<-rep((lambda[id.lambda,1]^2)*(alpha+1)/2,length(b01))
                                              idbeta<-which(b01<=lambda[id.lambda,1])
                                              p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])
                                              idbeta<-which(abs(b01)<lambda[id.lambda,1]*alpha)
                                              p01[idbeta]<-(2*alpha*lambda[id.lambda,1]*abs(b01[idbeta])-b01[idbeta]^2-lambda[id.lambda,1]^2)/(2*(alpha-1))
                                              
                                              p02<-rep((lambda[id.lambda,2]^2)*(alpha+1)/2,length(b02))
                                              idbeta<-which(b02<=lambda[id.lambda,2])
                                              p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])
                                              idbeta<-which(abs(b02)<lambda[id.lambda,2]*alpha)
                                              p02[idbeta]<-(2*alpha*lambda[id.lambda,2]*abs(b02[idbeta])-b02[idbeta]^2-lambda[id.lambda,2]^2)/(2*(alpha-1))
                                              
                                              p12<-rep((lambda[id.lambda,3]^2)*(alpha+1)/2,length(b12))
                                              idbeta<-which(b12<=lambda[id.lambda,3])
                                              p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])
                                              idbeta<-which(abs(b12)<lambda[id.lambda,3]*alpha)
                                              p12[idbeta]<-(2*alpha*lambda[id.lambda,3]*abs(b12[idbeta])-b12[idbeta]^2-lambda[id.lambda,3]^2)/(2*(alpha-1))
                                              
                                              fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                              
                                            }
                                          }else{
                                            snew<-s
                                            fn.valuenew<-res
                                          }
                                          
                                          ite<-ite+1
                                          
                                          #check cv 
                                          eval.cv.spline[ite]<-sum((snew-s)^2)
                                          eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                          eval.cv.loglik[ite]<-abs((fn.valuenew-fn.value)/fn.value)
                                          eval.loglik[ite]<-fn.valuenew
                                          eval.validity[ite]<-validity
                                          
                                          s<-snew
                                          beta<-betanew
                                          fn.value<-fn.valuenew
                                          
                                          # eval.cv beta valid only if validity.param=T
                                          if(((eval.cv.beta[ite] + eval.cv.spline[ite])<epsa) & eval.cv.loglik[ite]<epsb & validity==T){
                                            converged<-T}
                                          
                                          
                                        }
                                        
                                        if(maxiter<=ite & converged==F){
                                          istop<-2
                                        }else{
                                          if(ite<=maxiter & converged==T){
                                            istop<-1 
                                            
                                            # need to recalculate second derivatives 
                                            # if converged 
                                            
                                            b<-c(s,beta)
                                            bfix<-b[fix0==1]
                                            b<-b[fix0==0]
                                            
                                            output<-DYNderivaweib( h=1e-8,
                                                                   b=b,
                                                                   npm=npm,
                                                                   npar=size_V,
                                                                   bfix=bfix,
                                                                   fix=fix0,
                                                                   ctime=ctime,
                                                                   no=N,
                                                                   ve01=ve01,
                                                                   ve02=ve02,
                                                                   ve12=ve12,
                                                                   dimnva01=dimnva01,
                                                                   dimnva02=dimnva02,
                                                                   dimnva12=dimnva12,
                                                                   nva01=nvat01,
                                                                   nva02=nvat02,
                                                                   nva12=nvat12,
                                                                   t0=t0,
                                                                   t1=t1,
                                                                   t2=t2,
                                                                   t3=t3,
                                                                   troncature=troncature,
                                                                   y01=y01k,
                                                                   y02=y02k,
                                                                   y12=y12k,
                                                                   p01=p01,
                                                                   p02=p02,
                                                                   p12=p12,
                                                                   dimp01=dimp01,
                                                                   dimp02=dimp02,
                                                                   dimp12=dimp12,
                                                                   Ntime=NtimePoints)
                                            
                                            output<-output$v
                                            
                                            fu<-output[(min+1):length(output)]
                                            V<-matrix(0,nrow=npm,ncol=npm)
                                            V[upper.tri(V,diag=T)]<-output[1:min]
                                            V<-V+t(V)
                                            diag(V)<-diag(V)/2
                                            
                                            # deriva gives information matrix
                                            V0<-V
                                          }else{
                                            if(pbr_compu==1){istop<-3}
                                            if(pbr_compu==2){istop<-4}
                                            if(pbr_compu==3){istop<-5}
                                          }
                                        }
                                        
                                        # if stop==1 we can give matrix of second derivatives 
                                        
                                        return(list(b=c(s,beta),
                                                    H=V0,
                                                    lambda=as.double(lambda[id.lambda,]),
                                                    alpha=alpha,
                                                    fn.value=ifelse(!exists("output.mla"),NA,output.mla$fn.value), # loglik
                                                    fn.value.pena=fn.value, # penalised loglik
                                                    ni=ite,
                                                    ca.beta=eval.cv.beta,
                                                    ca.spline=eval.cv.spline,
                                                    ca.validity=eval.validity,
                                                    cb=eval.loglik,
                                                    istop=istop,
                                                    combine=id.lambda))
                                      }
    }else{
      outputNsample<-foreach::foreach(id.lambda=1:nlambda,
                                      .combine = combine_lambda_mla,
                                      .errorhandling = "remove")%do%{
                                        
                                        # computation pbr 
                                        
                                        pbr_compu<-0
                                        
                                          beta<-beta.start
                                          s<-s.start
                                       
                                        
                                        converged<-F
                                        ite<-0
                                        # if beta not change do not need to recalculate weights 
                                        H<-T
                                        
                                        eval.cv.spline<-rep(NA,maxiter+1)
                                        eval.cv.beta<-rep(NA,maxiter+1)
                                        eval.cv.loglik<-rep(NA,maxiter+1)
                                        eval.loglik<-rep(NA,maxiter+1)
                                        eval.validity<-rep(NA,maxiter+1)
                                        
                                        
                                        while(converged==F & ite<=maxiter){
                                          
                                          
                                          b<-c(s,beta)
                                          bfix<-b[fix0==1]
                                          b<-b[fix0==0]
                                          # derivative of loglik
                                          
                                          
                                          output<-DYNderivaweibdiag( h=1e-8,
                                                                     b=b,
                                                                     npm=npm,
                                                                     npar=size_V,
                                                                     bfix=bfix,
                                                                     fix=fix0,
                                                                     ctime=ctime,
                                                                     no=N,
                                                                     ve01=ve01,
                                                                     ve02=ve02,
                                                                     ve12=ve12,
                                                                     dimnva01=dimnva01,
                                                                     dimnva02=dimnva02,
                                                                     dimnva12=dimnva12,
                                                                     nva01=nvat01,
                                                                     nva02=nvat02,
                                                                     nva12=nvat12,
                                                                     t0=t0,
                                                                     t1=t1,
                                                                     t2=t2,
                                                                     t3=t3,
                                                                     troncature=troncature,
                                                                     y01=y01k,
                                                                     y02=y02k,
                                                                     y12=y12k,
                                                                     p01=p01,
                                                                     p02=p02,
                                                                     p12=p12,
                                                                     dimp01=dimp01,
                                                                     dimp02=dimp02,
                                                                     dimp12=dimp12,
                                                                     Ntime=NtimePoints)
                                          output<-output$v
                                          
                                          if(ite==0){# loglik penalised
                                            fn.value<-gaussDYNidmlLikelihoodweibpena(b=b,
                                                                                     npm=npm,
                                                                                     npar=size_V,
                                                                                     bfix=bfix,
                                                                                     fix=fix0,
                                                                                     ctime=ctime,
                                                                                     no=N,
                                                                                     ve01=ve01,
                                                                                     ve02=ve02,
                                                                                     ve12=ve12,
                                                                                     dimnva01=dimnva01,
                                                                                     dimnva02=dimnva02,
                                                                                     dimnva12=dimnva12,
                                                                                     nva01=nvat01,
                                                                                     nva02=nvat02,
                                                                                     nva12=nvat12,
                                                                                     t0=t0,
                                                                                     t1=t1,
                                                                                     t2=t2,
                                                                                     t3=t3,
                                                                                     troncature=troncature,
                                                                                     y01=y01k,
                                                                                     y02=y02k,
                                                                                     y12=y12k,
                                                                                     p01=p01,
                                                                                     p02=p02,
                                                                                     p12=p12,
                                                                                     dimp01=dimp01,
                                                                                     dimp02=dimp02,
                                                                                     dimp12=dimp12,
                                                                                     Ntime=NtimePoints,
                                                                                     
                                                                                     lambda=lambda[id.lambda,],
                                                                                     alpha=alpha,
                                                                                     penalty.factor=penalty.factor,
                                                                                     penalty=penalty)
                                          }
                                          
                                          if(any(is.na(output))|any(output==Inf) |any(output==-Inf)){
                                            warning("Computational error for calculation of the hessian : division by 0 or Infinite value")
                                            if(ite==0){
                                              
                                              fu<-output[(min+1):length(output)]
                                              V<-matrix(0,nrow=npm,ncol=npm)
                                              diag(V)<-output[1:min]
                                              # deriva gives information matrix
                                              tr <- sum(diag(V))/npm
                                              V0<-V}
                                            ite<-ite+1
                                            pbr_compu<-1
                                            break
                                          }
                                          
                                          
                                          
                                          fu<-output[(min+1):length(output)]
                                          V<-matrix(0,nrow=npm,ncol=npm)
                                          diag(V)<-output[1:min]
                                          # deriva gives information matrix
                                          tr <- sum(diag(V))/npm
                                          V0<-V
                                          
                                          eigen.values<-diag(V)
                                          
                                          if(defpositive==T){
                                            idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                          }else{
                                            idpos<-ifelse(any(diag(V)==0),1,0)
                                          }
                                          
                                          
                                          idpos0<-idpos
                                          
                                          ncount<-da<-ga<-0
                                          
                                          while(idpos != 0){
                                            
                                            if(ncount==0){ 
                                              ga <- 0.01
                                              da <- 1E-2
                                            }else{
                                              if(((ncount <= 3) | (ga >= 1)) ){
                                                da <- da * 5
                                              }else{# if ncount > 10 only update ga 
                                                ga <- ga * 5
                                                # do not put ga at 1 as no countmax otherwise infinite while 
                                                if(ga > 1) ga <- 1
                                              }
                                            }
                                            
                                            ncount <- ncount + 1
                                            
                                            diagV <- diag(V)
                                            # put abs (1-ga) better than 1-ga cause ga can now be >1
                                            diagV<-ifelse(diagV!=0,diagV+da*(abs((1.e0-ga))*abs(diagV)+ga*tr),
                                                          da*ga*tr)
                                            
                                            diag(V)<-diagV
                                            # if we have a convex log-vraisemblance in eta then :
                                            # all eigen  values of the hessienne are >0.
                                            
                                            if(sum(V==Inf)>0|sum(V==-Inf)>0){break}
                                            eigen.values<-diag(V)
                                            # check if hessienne defined positive
                                            
                                            if(defpositive==T){
                                              idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                            }else{
                                              idpos<-ifelse(any(diag(V)==0),1,0)
                                            }
                                            
                                            # if(def.positive==T){
                                            #   idpos<-ifelse(any(eigen.values<=0),1,0)
                                            # }else{idpos<-ifelse(any(abs(eigen.values)==0),1,0)}
                                            
                                          }
                                          
                                          if(idpos!=0){
                                            
                                            warning("Hessian not defined positive")
                                            pbr_compu<-2
                                            ite<-ite+1
                                            break
                                          }
                                          
                                          
                                          
                                          
                                          
                                          output.cv<-DYNcv.model(beta=beta,
                                                                 nva01=npm01,
                                                                 nva02=npm02,
                                                                 nva12=npm12,
                                                                 nva01Y=npm01Y,
                                                                 nva02Y=npm02Y,
                                                                 nva12Y=npm12Y,
                                                                 fix=fix0[7:size_V],
                                                                 penalty.factor=penalty.factor,
                                                                 penalty=penalty,
                                                                 v=V,
                                                                 fu=fu,
                                                                 lambda=lambda[id.lambda,],
                                                                 alpha=alpha
                                          )
                                          
                                          # verify validity of parameters update 
                                          # and that we are better than previous estimates 
                                          
                                          b<-c(s,output.cv$b)
                                          
                                          betanew<-b[(6+1):size_V]
                                          
                                          # penalised loglik see if inferior to previous
                                          res<-gaussDYNidmlLikelihoodweibpena(b=b,
                                                                              npm=size_V,
                                                                              npar=size_V,
                                                                              bfix=1,
                                                                              fix=rep(0,size_V),
                                                                              ctime=ctime,
                                                                              no=N,
                                                                              ve01=ve01,
                                                                              ve02=ve02,
                                                                              ve12=ve12,
                                                                              dimnva01=dimnva01,
                                                                              dimnva02=dimnva02,
                                                                              dimnva12=dimnva12,
                                                                              nva01=nvat01,
                                                                              nva02=nvat02,
                                                                              nva12=nvat12,
                                                                              t0=t0,
                                                                              t1=t1,
                                                                              t2=t2,
                                                                              t3=t3,
                                                                              troncature=troncature,
                                                                              y01=y01k,
                                                                              y02=y02k,
                                                                              y12=y12k,
                                                                              p01=p01,
                                                                              p02=p02,
                                                                              p12=p12,
                                                                              dimp01=dimp01,
                                                                              dimp02=dimp02,
                                                                              dimp12=dimp12,
                                                                              Ntime=NtimePoints,
                                                                              
                                                                              lambda=lambda[id.lambda,],
                                                                              alpha=alpha,
                                                                              penalty.factor=penalty.factor,
                                                                              penalty=penalty)
                                          
                                          
                                          # we want to maximise the loglik thus : 
                                          # we have issue if res is NA or if not higher than previous one 
                                          # if not better or do not exist need to readjust
                                          # value of beta 
                                          if(res %in%c(-1e9,1e9) | res < fn.value){
                                            
                                            th<-1e-5
                                            step<-log(1.5)
                                            delta<-output.cv$b-c(beta)
                                            
                                            maxt <- max(abs(delta)) 
                                            
                                            if(maxt == 0){
                                              vw <- th
                                            }else{
                                              vw <- th/maxt
                                            }
                                            if(ite>0){
                                              res.out.error <- list("old.b"=round(c(s,beta)),
                                                                    "old.rl"=round(fn.value),
                                                                    "old.ca"=round(eval.cv.beta[ite]),
                                                                    "old.cb"=round(eval.cv.loglik[ite]))
                                            }else{
                                              res.out.error <- list("old.b"=round(c(s,beta)),
                                                                    "old.rl"=round(fn.value),
                                                                    "old.ca"=round(1),
                                                                    "old.cb"=round(1))
                                            }
                                            # from mla package 
                                            sears<-searpas(vw=vw,
                                                           step=step,
                                                           b=beta,
                                                           delta=delta,
                                                           funcpa=gaussDYNidmlLikelihoodweibpena,
                                                           res.out.error=res.out.error,
                                                           npm=npm,
                                                           npar=size_V,
                                                           bfix=s,
                                                           fix=fix0,
                                                           ctime=ctime,
                                                           no=N,
                                                           ve01=ve01,
                                                           ve02=ve02,
                                                           ve12=ve12,
                                                           dimnva01=dimnva01,
                                                           dimnva02=dimnva02,
                                                           dimnva12=dimnva12,
                                                           nva01=nvat01,
                                                           nva02=nvat02,
                                                           nva12=nvat12,
                                                           t0=t0,
                                                           t1=t1,
                                                           t2=t2,
                                                           t3=t3,
                                                           troncature=troncature,
                                                           y01=y01k,
                                                           y02=y02k,
                                                           y12=y12k,
                                                           p01=p01,
                                                           p02=p02,
                                                           p12=p12,
                                                           dimp01=dimp01,
                                                           dimp02=dimp02,
                                                           dimp12=dimp12,
                                                           Ntime=NtimePoints,
                                                           lambda=lambda[id.lambda,],
                                                           alpha=alpha,
                                                           penalty.factor=penalty.factor,
                                                           penalty=penalty)
                                            
                                            betanew<-beta+delta*sears$vw
                                            betanew<-ifelse(abs(betanew)<=0.0001,0,betanew)
                                            b<-c(s,betanew)
                                            
                                            
                                            res<-gaussDYNidmlLikelihoodweibpena(b=b,
                                                                                npm=size_V,
                                                                                npar=size_V,
                                                                                bfix=1,
                                                                                fix=rep(0,size_V),
                                                                                ctime=ctime,
                                                                                no=N,
                                                                                ve01=ve01,
                                                                                ve02=ve02,
                                                                                ve12=ve12,
                                                                                dimnva01=dimnva01,
                                                                                dimnva02=dimnva02,
                                                                                dimnva12=dimnva12,
                                                                                nva01=nvat01,
                                                                                nva02=nvat02,
                                                                                nva12=nvat12,
                                                                                t0=t0,
                                                                                t1=t1,
                                                                                t2=t2,
                                                                                t3=t3,
                                                                                troncature=troncature,
                                                                                y01=y01k,
                                                                                y02=y02k,
                                                                                y12=y12k,
                                                                                p01=p01,
                                                                                p02=p02,
                                                                                p12=p12,
                                                                                dimp01=dimp01,
                                                                                dimp02=dimp02,
                                                                                dimp12=dimp12,
                                                                                Ntime=NtimePoints,
                                                                                
                                                                                lambda=lambda[id.lambda,],
                                                                                alpha=alpha,
                                                                                penalty.factor=penalty.factor,
                                                                                penalty=penalty)
                                          }
                                          # if not better or do not exist need to readjust
                                          # value of beta 
                                          if(res %in%c(-1e9,1e9) | any(is.infinite(c(s,betanew)))){
                                            
                                            ite<-ite+1
                                            validity<-F
                                            eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                            eval.validity[ite]<-validity
                                            # save output over iterations 
                                            
                                            pbr_compu<-3
                                            break
                                          }else{validity<-T}
                                          
                                          
                                          # betanew already include s
                                          b<-c(s,betanew)
                                          
                                          bfix<-b[fix0.beta==1]
                                          b<-b[fix0.beta==0]
                                          # update for modelPar
                                          if(npmweib!=0){
                                            output.mla<- marqLevAlg::mla(b=b,
                                                                         fn=gaussDYNidmlLikelihoodweib,
                                                                         gr=DYNreggrmlaweibana,
                                                                         epsa=epsa,
                                                                         epsb=epsb,
                                                                         epsd=epsd,
                                                                         maxiter=maxiter.pena,
                                                                         minimize=F,
                                                                         npm=npmweib,
                                                                         npar=size_V,
                                                                         bfix=bfix,
                                                                         fix=fix0.beta,
                                                                         ctime=ctime,
                                                                         no=N,
                                                                         ve01=ve01,
                                                                         ve02=ve02,
                                                                         ve12=ve12,
                                                                         dimnva01=dimnva01,
                                                                         dimnva02=dimnva02,
                                                                         dimnva12=dimnva12,
                                                                         nva01=nvat01,
                                                                         nva02=nvat02,
                                                                         nva12=nvat12,
                                                                         t0=t0,
                                                                         t1=t1,
                                                                         t2=t2,
                                                                         t3=t3,
                                                                         troncature=troncature,
                                                                         y01=y01k,
                                                                         y02=y02k,
                                                                         y12=y12k,
                                                                         p01=p01,
                                                                         p02=p02,
                                                                         p12=p12,
                                                                         dimp01=dimp01,
                                                                         dimp02=dimp02,
                                                                         dimp12=dimp12,
                                                                         Ntime=NtimePoints)
                                            
                                            # look at convergence for each lambda :
                                            # mla output is loglik
                                            # new values for splines:
                                            snew<-s
                                            snew[fix00[1:6]==0]<-output.mla$b
                                            if(nvat01>0){
                                              b01<-betanew[1:nvat01][penalty.factor[1:nvat01]==1]
                                              if(p01>0){
                                                b01<-c(b01,betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1])
                                              }
                                            }else{
                                              if(p01>0){
                                                b01<-betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1]
                                              }else{
                                                b01<-0
                                              }
                                            }
                                            
                                            if(nvat02>0){
                                              b02<-betanew[(nvat01+1):(nvat01+nvat02)][penalty.factor[(nvat01+1):(nvat01+nvat02)]==1]
                                              if(p02>0){
                                                b02<-c(b02,betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1])
                                              }
                                            }else{
                                              if(p02>0){
                                                b02<-betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1]
                                              }else{b02<-0}
                                            }
                                            
                                            if(nvat12>0){
                                              b12<-betanew[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)][penalty.factor[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)]==1]
                                              if(p12>0){
                                                b12<-c(b12,betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1])
                                              }
                                            }else{
                                              if(p12>0){
                                                b12<-betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1]
                                              }else{b12<-0}
                                            }
                                            
                                            # calculate loglik pen 
                                            if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
                                              fn.valuenew<-output.mla$fn.value-lambda[id.lambda,1]*alpha*sum(abs(b01))-lambda[id.lambda,1]*(1-alpha)*sum(b01*b01)
                                              fn.valuenew<-fn.valuenew-lambda[id.lambda,2]*alpha*sum(abs(b02))-lambda[id.lambda,2]*(1-alpha)*sum(b02*b02)
                                              fn.valuenew<-fn.valuenew-lambda[id.lambda,3]*alpha*sum(abs(b12))-lambda[id.lambda,3]*(1-alpha)*sum(b12*b12)
                                            }
                                            
                                            
                                            if(penalty=="mcp"){
                                              
                                              p01<-rep(alpha*lambda[id.lambda,1]*lambda[id.lambda,1]/2,length(b01))
                                              idbeta<-which(b01<=alpha*lambda[id.lambda,1])
                                              p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])-((b01[idbeta]*b01[idbeta])/2*alpha)
                                              
                                              p02<-rep(alpha*lambda[id.lambda,2]*lambda[id.lambda,2]/2,length(b02))
                                              idbeta<-which(b02<=alpha*lambda[id.lambda,2])
                                              p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])-((b02[idbeta]*b02[idbeta])/2*alpha)
                                              
                                              p12<-rep(alpha*lambda[id.lambda,3]*lambda[id.lambda,3]/2,length(b12))
                                              idbeta<-which(b12<=alpha*lambda[id.lambda,3])
                                              p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])-((b12[idbeta]*b12[idbeta])/2*alpha)
                                              
                                              fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                              
                                            }
                                            
                                            if(penalty=="scad"){
                                              
                                              p01<-rep((lambda[id.lambda,1]^2)*(alpha+1)/2,length(b01))
                                              idbeta<-which(b01<=lambda[id.lambda,1])
                                              p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])
                                              idbeta<-which(abs(b01)<lambda[id.lambda,1]*alpha)
                                              p01[idbeta]<-(2*alpha*lambda[id.lambda,1]*abs(b01[idbeta])-b01[idbeta]^2-lambda[id.lambda,1]^2)/(2*(alpha-1))
                                              
                                              p02<-rep((lambda[id.lambda,2]^2)*(alpha+1)/2,length(b02))
                                              idbeta<-which(b02<=lambda[id.lambda,2])
                                              p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])
                                              idbeta<-which(abs(b02)<lambda[id.lambda,2]*alpha)
                                              p02[idbeta]<-(2*alpha*lambda[id.lambda,2]*abs(b02[idbeta])-b02[idbeta]^2-lambda[id.lambda,2]^2)/(2*(alpha-1))
                                              
                                              p12<-rep((lambda[id.lambda,3]^2)*(alpha+1)/2,length(b12))
                                              idbeta<-which(b12<=lambda[id.lambda,3])
                                              p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])
                                              idbeta<-which(abs(b12)<lambda[id.lambda,3]*alpha)
                                              p12[idbeta]<-(2*alpha*lambda[id.lambda,3]*abs(b12[idbeta])-b12[idbeta]^2-lambda[id.lambda,3]^2)/(2*(alpha-1))
                                              
                                              fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                              
                                            }
                                          }else{
                                            snew<-s
                                            fn.valuenew<-res
                                          }
                                          
                                          ite<-ite+1
                                          
                                          #check cv 
                                          eval.cv.spline[ite]<-sum((snew-s)^2)
                                          eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                          eval.cv.loglik[ite]<-abs((fn.valuenew-fn.value)/fn.value)
                                          eval.loglik[ite]<-fn.valuenew
                                          eval.validity[ite]<-validity
                                          
                                          s<-snew
                                          beta<-betanew
                                          fn.value<-fn.valuenew
                                          
                                          # eval.cv beta valid only if validity.param=T
                                          if(((eval.cv.beta[ite] + eval.cv.spline[ite])<epsa) & eval.cv.loglik[ite]<epsb & validity==T){
                                            converged<-T}
                                          
                                          
                                        }
                                        
                                        if(maxiter<=ite & converged==F){
                                          istop<-2
                                        }else{
                                          if(ite<=maxiter & converged==T){
                                            istop<-1 
                                            
                                            # need to recalculate second derivatives 
                                            # if converged 
                                            
                                            b<-c(s,beta)
                                            bfix<-b[fix0==1]
                                            b<-b[fix0==0]
                                            
                                            output<-DYNderivaweibdiag( h=1e-8,
                                                                       b=b,
                                                                       npm=npm,
                                                                       npar=size_V,
                                                                       bfix=bfix,
                                                                       fix=fix0,
                                                                       ctime=ctime,
                                                                       no=N,
                                                                       ve01=ve01,
                                                                       ve02=ve02,
                                                                       ve12=ve12,
                                                                       dimnva01=dimnva01,
                                                                       dimnva02=dimnva02,
                                                                       dimnva12=dimnva12,
                                                                       nva01=nvat01,
                                                                       nva02=nvat02,
                                                                       nva12=nvat12,
                                                                       t0=t0,
                                                                       t1=t1,
                                                                       t2=t2,
                                                                       t3=t3,
                                                                       troncature=troncature,
                                                                       y01=y01k,
                                                                       y02=y02k,
                                                                       y12=y12k,
                                                                       p01=p01,
                                                                       p02=p02,
                                                                       p12=p12,
                                                                       dimp01=dimp01,
                                                                       dimp02=dimp02,
                                                                       dimp12=dimp12,
                                                                       Ntime=NtimePoints)
                                            
                                            output<-output$v
                                            fu<-output[(npm+1):length(output)]
                                            V<-matrix(0,nrow=npm,ncol=npm)
                                            diag(V)<-output[1:npm]
                                            
                                            # deriva gives information matrix
                                            V0<-V
                                          }else{
                                            if(pbr_compu==1){istop<-3}
                                            if(pbr_compu==2){istop<-4}
                                            if(pbr_compu==3){istop<-5}
                                          }
                                        }
                                        
                                        # if stop==1 we can give matrix of second derivatives 
                                        
                                        
                                        return(list(b=c(s,beta),
                                                    H=V0,
                                                    lambda=as.double(lambda[id.lambda,]),
                                                    alpha=alpha,
                                                    fn.value=ifelse(!exists("output.mla"),NA,output.mla$fn.value), # loglik
                                                    fn.value.pena=fn.value, # penalised loglik
                                                    ni=ite,
                                                    ca.beta=eval.cv.beta,
                                                    ca.spline=eval.cv.spline,
                                                    ca.validity=eval.validity,
                                                    cb=eval.loglik,
                                                    istop=istop,
                                                    combine=id.lambda))
                                      }
    }
  }else{
    outputNsample<-list()
    length(outputNsample)<-nlambda
     if(partialH==F){
       for(id.lambda in 1:nlambda){
         if(id.lambda>1){
           
           #ids <- which(!vapply(outputNsample, is.null, logical(1)))
           ids <- which(unlist(lapply(outputNsample, FUN=function(x){
             if(is.null(x)){return(F)}else{
               if(x$istop==1){
                 return(T)
               }else{return(F)}
             }})))
           last_id <- if (length(ids) == 0) NA_integer_ else max(ids)
           if(!is.na(last_id)){
             beta.start<-outputNsample[[last_id]]$b[7:size_V]
             s.start<-outputNsample[[last_id]]$b[1:6]}
         } 
         outputNsample[[id.lambda]]<-tryCatch({
                                        
                                        # computation pbr 
                                        
                                        pbr_compu<-0
                                        
                                        beta<-beta.start
                                        s<-s.start
                                        
                                        
                                        converged<-F
                                        ite<-0
                                        # if beta not change do not need to recalculate weights 
                                        H<-T
                                        
                                        eval.cv.spline<-rep(NA,maxiter+1)
                                        eval.cv.beta<-rep(NA,maxiter+1)
                                        eval.cv.loglik<-rep(NA,maxiter+1)
                                        eval.loglik<-rep(NA,maxiter+1)
                                        eval.validity<-rep(NA,maxiter+1)
                                        
                                        
                                        
                                        
                                        while(converged==F & ite<=maxiter){
                                          
                                          
                                          b<-c(s,beta)
                                          bfix<-b[fix0==1]
                                          b<-b[fix0==0]
                                          # derivative of loglik
                                          
                                          output<-DYNderivaweib( h=1e-8,
                                                                 b=b,
                                                                 npm=npm,
                                                                 npar=size_V,
                                                                 bfix=bfix,
                                                                 fix=fix0,
                                                                 ctime=ctime,
                                                                 no=N,
                                                                 ve01=ve01,
                                                                 ve02=ve02,
                                                                 ve12=ve12,
                                                                 dimnva01=dimnva01,
                                                                 dimnva02=dimnva02,
                                                                 dimnva12=dimnva12,
                                                                 nva01=nvat01,
                                                                 nva02=nvat02,
                                                                 nva12=nvat12,
                                                                 t0=t0,
                                                                 t1=t1,
                                                                 t2=t2,
                                                                 t3=t3,
                                                                 troncature=troncature,
                                                                 y01=y01k,
                                                                 y02=y02k,
                                                                 y12=y12k,
                                                                 p01=p01,
                                                                 p02=p02,
                                                                 p12=p12,
                                                                 dimp01=dimp01,
                                                                 dimp02=dimp02,
                                                                 dimp12=dimp12,
                                                                 Ntime=NtimePoints)
                                          output<-output$v
                                          
                                          if(ite==0){# loglik penalised
                                            fn.value<-gaussDYNidmlLikelihoodweibpena(b=b,
                                                                                     npm=npm,
                                                                                     npar=size_V,
                                                                                     bfix=bfix,
                                                                                     fix=fix0,
                                                                                     ctime=ctime,
                                                                                     no=N,
                                                                                     ve01=ve01,
                                                                                     ve02=ve02,
                                                                                     ve12=ve12,
                                                                                     dimnva01=dimnva01,
                                                                                     dimnva02=dimnva02,
                                                                                     dimnva12=dimnva12,
                                                                                     nva01=nvat01,
                                                                                     nva02=nvat02,
                                                                                     nva12=nvat12,
                                                                                     t0=t0,
                                                                                     t1=t1,
                                                                                     t2=t2,
                                                                                     t3=t3,
                                                                                     troncature=troncature,
                                                                                     y01=y01k,
                                                                                     y02=y02k,
                                                                                     y12=y12k,
                                                                                     p01=p01,
                                                                                     p02=p02,
                                                                                     p12=p12,
                                                                                     dimp01=dimp01,
                                                                                     dimp02=dimp02,
                                                                                     dimp12=dimp12,
                                                                                     Ntime=NtimePoints,
                                                                                     
                                                                                     lambda=lambda[id.lambda,],
                                                                                     alpha=alpha,
                                                                                     penalty.factor=penalty.factor,
                                                                                     penalty=penalty)
                                          }
                                          
                                          if(any(is.na(output))|any(output==Inf) |any(output==-Inf)){
                                            warning("Computational error for calculation of the hessian : division by 0 or Infinite value")
                                            if(ite==0){
                                              
                                              
                                              fu<-output[(min+1):length(output)]
                                              V<-matrix(0,nrow=npm,ncol=npm)
                                              V[upper.tri(V,diag=T)]<-output[1:min]
                                              V<-V+t(V)
                                              diag(V)<-diag(V)/2
                                              # deriva gives information matrix
                                              tr <- sum(diag(V))/npm
                                              V0<-V}
                                            ite<-ite+1
                                            pbr_compu<-1
                                            break
                                          }
                                          
                                          
                                          fu<-output[(min+1):length(output)]
                                          V<-matrix(0,nrow=npm,ncol=npm)
                                          V[upper.tri(V,diag=T)]<-output[1:min]
                                          V<-V+t(V)
                                          diag(V)<-diag(V)/2
                                          
                                          # deriva gives information matrix
                                          tr <- sum(diag(V))/npm
                                          V0<-V
                                          
                                          eigen.values<-eigen(V,symmetric=T,only.values=T)$values
                                          
                                          if(defpositive==T){
                                            idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                          }else{
                                            idpos<-ifelse(any(diag(V)==0),1,0)
                                          }
                                          
                                          
                                          idpos0<-idpos
                                          
                                          ncount<-da<-ga<-0
                                          
                                          while(idpos != 0){
                                            
                                            if(ncount==0){ 
                                              ga <- 0.01
                                              da <- 1E-2
                                            }else{
                                              if(((ncount <= 3) | (ga >= 1)) ){
                                                da <- da * 5
                                              }else{# if ncount > 10 only update ga 
                                                ga <- ga * 5
                                                # do not put ga at 1 as no countmax otherwise infinite while 
                                                if(ga > 1) ga <- 1
                                              }
                                            }
                                            
                                            ncount <- ncount + 1
                                            
                                            diagV <- diag(V)
                                            # put abs (1-ga) better than 1-ga cause ga can now be >1
                                            diagV<-ifelse(diagV!=0,diagV+da*(abs((1.e0-ga))*abs(diagV)+ga*tr),
                                                          da*ga*tr)
                                            
                                            diag(V)<-diagV
                                            # if we have a convex log-vraisemblance in eta then :
                                            # all eigen  values of the hessienne are >0.
                                            
                                            if(sum(V==Inf)>0|sum(V==-Inf)>0){break}
                                            eigen.values<-eigen(V,symmetric=T,only.values=T)$values
                                            # check if hessienne defined positive
                                            
                                            
                                            if(defpositive==T){
                                              idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                            }else{
                                              idpos<-ifelse(any(diag(V)==0),1,0)
                                            }
                                            
                                            # if(def.positive==T){
                                            #   idpos<-ifelse(any(eigen.values<=0),1,0)
                                            # }else{idpos<-ifelse(any(abs(eigen.values)==0),1,0)}
                                            
                                          }
                                          
                                          
                                          if(idpos!=0){
                                            
                                            warning("Hessian not defined positive")
                                            pbr_compu<-2
                                            ite<-ite+1
                                            break
                                          }
                                          
                                          
                                          # update for beta 
                                          output.cv<-DYNcv.model(beta=beta,
                                                                 nva01=npm01,
                                                                 nva02=npm02,
                                                                 nva12=npm12,
                                                                 nva01Y=npm01Y,
                                                                 nva02Y=npm02Y,
                                                                 nva12Y=npm12Y,
                                                                 fix=fix0[7:size_V],
                                                                 penalty.factor=penalty.factor,
                                                                 penalty=penalty,
                                                                 v=V,
                                                                 fu=fu,
                                                                 lambda=lambda[id.lambda,],
                                                                 alpha=alpha
                                          )
                                          
                                          # verify validity of parameters update 
                                          # and that we are better than previous estimates 
                                          
                                          b<-c(s,output.cv$b)
                                          
                                          betanew<-b[(6+1):size_V]
                                          
                                          # penalised loglik see if inferior to previous
                                          res<-gaussDYNidmlLikelihoodweibpena(b=b,
                                                                              npm=size_V,
                                                                              npar=size_V,
                                                                              bfix=1,
                                                                              fix=rep(0,size_V),
                                                                              ctime=ctime,
                                                                              no=N,
                                                                              ve01=ve01,
                                                                              ve02=ve02,
                                                                              ve12=ve12,
                                                                              dimnva01=dimnva01,
                                                                              dimnva02=dimnva02,
                                                                              dimnva12=dimnva12,
                                                                              nva01=nvat01,
                                                                              nva02=nvat02,
                                                                              nva12=nvat12,
                                                                              t0=t0,
                                                                              t1=t1,
                                                                              t2=t2,
                                                                              t3=t3,
                                                                              troncature=troncature,
                                                                              y01=y01k,
                                                                              y02=y02k,
                                                                              y12=y12k,
                                                                              p01=p01,
                                                                              p02=p02,
                                                                              p12=p12,
                                                                              dimp01=dimp01,
                                                                              dimp02=dimp02,
                                                                              dimp12=dimp12,
                                                                              Ntime=NtimePoints,
                                                                              lambda=lambda[id.lambda,],
                                                                              alpha=alpha,
                                                                              penalty.factor=penalty.factor,
                                                                              penalty=penalty)
                                          
                                          
                                          # we want to maximise the loglik thus : 
                                          # we have issue if res is NA or if not higher than previous one 
                                          # if not better or do not exist need to readjust
                                          # value of beta 
                                          if(res %in%c(-1e9,1e9) | res < fn.value){
                                            
                                            th<-1e-5
                                            step<-log(1.5)
                                            delta<-output.cv$b-c(beta)
                                            
                                            maxt <- max(abs(delta)) 
                                            
                                            if(maxt == 0){
                                              vw <- th
                                            }else{
                                              vw <- th/maxt
                                            }
                                            if(ite>0){
                                              res.out.error <- list("old.b"=round(c(s,beta)),
                                                                    "old.rl"=round(fn.value),
                                                                    "old.ca"=round(eval.cv.beta[ite]),
                                                                    "old.cb"=round(eval.cv.loglik[ite]))
                                            }else{
                                              res.out.error <- list("old.b"=round(c(s,beta)),
                                                                    "old.rl"=round(fn.value),
                                                                    "old.ca"=round(1),
                                                                    "old.cb"=round(1))
                                            }
                                            # from mla package 
                                            sears<-searpas(vw=vw,
                                                           step=step,
                                                           b=beta,
                                                           delta=delta,
                                                           funcpa=gaussDYNidmlLikelihoodweibpena,
                                                           res.out.error=res.out.error,
                                                           npm=npm,
                                                           npar=size_V,
                                                           bfix=s,
                                                           fix=fix0,
                                                           ctime=ctime,
                                                           no=N,
                                                           ve01=ve01,
                                                           ve02=ve02,
                                                           ve12=ve12,
                                                           dimnva01=dimnva01,
                                                           dimnva02=dimnva02,
                                                           dimnva12=dimnva12,
                                                           nva01=nvat01,
                                                           nva02=nvat02,
                                                           nva12=nvat12,
                                                           t0=t0,
                                                           t1=t1,
                                                           t2=t2,
                                                           t3=t3,
                                                           troncature=troncature,
                                                           y01=y01k,
                                                           y02=y02k,
                                                           y12=y12k,
                                                           p01=p01,
                                                           p02=p02,
                                                           p12=p12,
                                                           dimp01=dimp01,
                                                           dimp02=dimp02,
                                                           dimp12=dimp12,
                                                           Ntime=NtimePoints,
                                                           lambda=lambda[id.lambda,],
                                                           alpha=alpha,
                                                           penalty.factor=penalty.factor,
                                                           penalty=penalty)
                                            
                                            
                                            betanew<-beta+delta*sears$vw
                                            betanew<-ifelse(abs(betanew)<=0.0001,0,betanew)
                                            b<-c(s,betanew)
                                            
                                            
                                            res<-gaussDYNidmlLikelihoodweibpena(b=b,
                                                                                npm=size_V,
                                                                                npar=size_V,
                                                                                bfix=1,
                                                                                fix=rep(0,size_V),
                                                                                ctime=ctime,
                                                                                no=N,
                                                                                ve01=ve01,
                                                                                ve02=ve02,
                                                                                ve12=ve12,
                                                                                dimnva01=dimnva01,
                                                                                dimnva02=dimnva02,
                                                                                dimnva12=dimnva12,
                                                                                nva01=nvat01,
                                                                                nva02=nvat02,
                                                                                nva12=nvat12,
                                                                                t0=t0,
                                                                                t1=t1,
                                                                                t2=t2,
                                                                                t3=t3,
                                                                                troncature=troncature,
                                                                                y01=y01k,
                                                                                y02=y02k,
                                                                                y12=y12k,
                                                                                p01=p01,
                                                                                p02=p02,
                                                                                p12=p12,
                                                                                dimp01=dimp01,
                                                                                dimp02=dimp02,
                                                                                dimp12=dimp12,
                                                                                Ntime=NtimePoints,
                                                                                lambda=lambda[id.lambda,],
                                                                                alpha=alpha,
                                                                                penalty.factor=penalty.factor,
                                                                                penalty=penalty)
                                          }
                                          # if not better or do not exist need to readjust
                                          # value of beta 
                                          if(res %in%c(-1e9,1e9) | any(is.infinite(c(s,betanew)))){
                                            
                                            ite<-ite+1
                                            validity<-F
                                            eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                            eval.validity[ite]<-validity
                                            # save output over iterations 
                                            
                                            pbr_compu<-3
                                            break
                                          }else{validity<-T}
                                          
                                          
                                          # betanew already include s
                                          b<-c(s,betanew)
                                          
                                          bfix<-b[fix0.beta==1]
                                          b<-b[fix0.beta==0]
                                          # update for modelPar
                                          if(npmweib!=0){
                                            output.mla<- marqLevAlg::mla(b=b,
                                                                         fn=gaussDYNidmlLikelihoodweib,
                                                                         gr=DYNreggrmlaweibana,
                                                                         epsa=epsa,
                                                                         epsb=epsb,
                                                                         epsd=epsd,
                                                                         maxiter=maxiter.pena,
                                                                         minimize=F,
                                                                         npm= npmweib,
                                                                         npar=size_V,
                                                                         bfix=bfix,
                                                                         fix=fix0.beta,
                                                                         ctime=ctime,
                                                                         no=N,
                                                                         ve01=ve01,
                                                                         ve02=ve02,
                                                                         ve12=ve12,
                                                                         dimnva01=dimnva01,
                                                                         dimnva02=dimnva02,
                                                                         dimnva12=dimnva12,
                                                                         nva01=nvat01,
                                                                         nva02=nvat02,
                                                                         nva12=nvat12,
                                                                         t0=t0,
                                                                         t1=t1,
                                                                         t2=t2,
                                                                         t3=t3,
                                                                         troncature=troncature,
                                                                         y01=y01k,
                                                                         y02=y02k,
                                                                         y12=y12k,
                                                                         p01=p01,
                                                                         p02=p02,
                                                                         p12=p12,
                                                                         dimp01=dimp01,
                                                                         dimp02=dimp02,
                                                                         dimp12=dimp12,
                                                                         Ntime=NtimePoints)
                                            
                                            # look at convergence for each lambda :
                                            # mla output is loglik
                                            # new values for splines:
                                            snew<-s
                                            snew[fix00[1:6]==0]<-output.mla$b
                                            
                                            
                                            
                                            if(nvat01>0){
                                              b01<-betanew[1:nvat01][penalty.factor[1:nvat01]==1]
                                              if(p01>0){
                                                b01<-c(b01,betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1])
                                              }
                                            }else{
                                              if(p01>0){
                                                b01<-betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1]
                                              }else{
                                                b01<-0
                                              }
                                            }
                                            
                                            if(nvat02>0){
                                              b02<-betanew[(nvat01+1):(nvat01+nvat02)][penalty.factor[(nvat01+1):(nvat01+nvat02)]==1]
                                              if(p02>0){
                                                b02<-c(b02,betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1])
                                              }
                                            }else{
                                              if(p02>0){
                                                b02<-betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1]
                                              }else{b02<-0}
                                            }
                                            
                                            if(nvat12>0){
                                              b12<-betanew[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)][penalty.factor[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)]==1]
                                              if(p12>0){
                                                b12<-c(b12,betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1])
                                              }
                                            }else{
                                              if(p12>0){
                                                b12<-betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1]
                                              }else{b12<-0}
                                            }
                                            
                                            # calculate loglik pen 
                                            if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
                                              fn.valuenew<-output.mla$fn.value-lambda[id.lambda,1]*alpha*sum(abs(b01))-lambda[id.lambda,1]*(1-alpha)*sum(b01*b01)
                                              fn.valuenew<-fn.valuenew-lambda[id.lambda,2]*alpha*sum(abs(b02))-lambda[id.lambda,2]*(1-alpha)*sum(b02*b02)
                                              fn.valuenew<-fn.valuenew-lambda[id.lambda,3]*alpha*sum(abs(b12))-lambda[id.lambda,3]*(1-alpha)*sum(b12*b12)
                                            }
                                            
                                            
                                            if(penalty=="mcp"){
                                              
                                              p01<-rep(alpha*lambda[id.lambda,1]*lambda[id.lambda,1]/2,length(b01))
                                              idbeta<-which(b01<=alpha*lambda[id.lambda,1])
                                              p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])-((b01[idbeta]*b01[idbeta])/2*alpha)
                                              
                                              p02<-rep(alpha*lambda[id.lambda,2]*lambda[id.lambda,2]/2,length(b02))
                                              idbeta<-which(b02<=alpha*lambda[id.lambda,2])
                                              p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])-((b02[idbeta]*b02[idbeta])/2*alpha)
                                              
                                              p12<-rep(alpha*lambda[id.lambda,3]*lambda[id.lambda,3]/2,length(b12))
                                              idbeta<-which(b12<=alpha*lambda[id.lambda,3])
                                              p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])-((b12[idbeta]*b12[idbeta])/2*alpha)
                                              
                                              fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                              
                                            }
                                            
                                            if(penalty=="scad"){
                                              
                                              p01<-rep((lambda[id.lambda,1]^2)*(alpha+1)/2,length(b01))
                                              idbeta<-which(b01<=lambda[id.lambda,1])
                                              p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])
                                              idbeta<-which(abs(b01)<lambda[id.lambda,1]*alpha)
                                              p01[idbeta]<-(2*alpha*lambda[id.lambda,1]*abs(b01[idbeta])-b01[idbeta]^2-lambda[id.lambda,1]^2)/(2*(alpha-1))
                                              
                                              p02<-rep((lambda[id.lambda,2]^2)*(alpha+1)/2,length(b02))
                                              idbeta<-which(b02<=lambda[id.lambda,2])
                                              p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])
                                              idbeta<-which(abs(b02)<lambda[id.lambda,2]*alpha)
                                              p02[idbeta]<-(2*alpha*lambda[id.lambda,2]*abs(b02[idbeta])-b02[idbeta]^2-lambda[id.lambda,2]^2)/(2*(alpha-1))
                                              
                                              p12<-rep((lambda[id.lambda,3]^2)*(alpha+1)/2,length(b12))
                                              idbeta<-which(b12<=lambda[id.lambda,3])
                                              p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])
                                              idbeta<-which(abs(b12)<lambda[id.lambda,3]*alpha)
                                              p12[idbeta]<-(2*alpha*lambda[id.lambda,3]*abs(b12[idbeta])-b12[idbeta]^2-lambda[id.lambda,3]^2)/(2*(alpha-1))
                                              
                                              fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                              
                                            }
                                          }else{
                                            snew<-s
                                            fn.valuenew<-res
                                          }
                                          
                                          ite<-ite+1
                                          
                                          #check cv 
                                          eval.cv.spline[ite]<-sum((snew-s)^2)
                                          eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                          eval.cv.loglik[ite]<-abs((fn.valuenew-fn.value)/fn.value)
                                          eval.loglik[ite]<-fn.valuenew
                                          eval.validity[ite]<-validity
                                          
                                          s<-snew
                                          beta<-betanew
                                          fn.value<-fn.valuenew
                                          
                                          # eval.cv beta valid only if validity.param=T
                                          if(((eval.cv.beta[ite] + eval.cv.spline[ite])<epsa) & eval.cv.loglik[ite]<epsb & validity==T){
                                            converged<-T}
                                          
                                          
                                        }
                                        
                                        if(maxiter<=ite & converged==F){
                                          istop<-2
                                        }else{
                                          if(ite<=maxiter & converged==T){
                                            istop<-1 
                                            
                                            # need to recalculate second derivatives 
                                            # if converged 
                                            
                                            b<-c(s,beta)
                                            bfix<-b[fix0==1]
                                            b<-b[fix0==0]
                                            
                                            output<-DYNderivaweib( h=1e-8,
                                                                   b=b,
                                                                   npm=npm,
                                                                   npar=size_V,
                                                                   bfix=bfix,
                                                                   fix=fix0,
                                                                   ctime=ctime,
                                                                   no=N,
                                                                   ve01=ve01,
                                                                   ve02=ve02,
                                                                   ve12=ve12,
                                                                   dimnva01=dimnva01,
                                                                   dimnva02=dimnva02,
                                                                   dimnva12=dimnva12,
                                                                   nva01=nvat01,
                                                                   nva02=nvat02,
                                                                   nva12=nvat12,
                                                                   t0=t0,
                                                                   t1=t1,
                                                                   t2=t2,
                                                                   t3=t3,
                                                                   troncature=troncature,
                                                                   y01=y01k,
                                                                   y02=y02k,
                                                                   y12=y12k,
                                                                   p01=p01,
                                                                   p02=p02,
                                                                   p12=p12,
                                                                   dimp01=dimp01,
                                                                   dimp02=dimp02,
                                                                   dimp12=dimp12,
                                                                   Ntime=NtimePoints)
                                            
                                            output<-output$v
                                            
                                            fu<-output[(min+1):length(output)]
                                            V<-matrix(0,nrow=npm,ncol=npm)
                                            V[upper.tri(V,diag=T)]<-output[1:min]
                                            V<-V+t(V)
                                            diag(V)<-diag(V)/2
                                            
                                            # deriva gives information matrix
                                            V0<-V
                                          }else{
                                            if(pbr_compu==1){istop<-3}
                                            if(pbr_compu==2){istop<-4}
                                            if(pbr_compu==3){istop<-5}
                                          }
                                        }
                                        
                                        # if stop==1 we can give matrix of second derivatives 
                                        
                                        
                                        list(b=c(s,beta),
                                                    H=V0,
                                                    lambda=as.double(lambda[id.lambda,]),
                                                    alpha=alpha,
                                                    fn.value=ifelse(!exists("output.mla"),NA,output.mla$fn.value), # loglik
                                                    fn.value.pena=fn.value, # penalised loglik
                                                    ni=ite,
                                                    ca.beta=eval.cv.beta,
                                                    ca.spline=eval.cv.spline,
                                                    ca.validity=eval.validity,
                                                    cb=eval.loglik,
                                                    istop=istop,
                                                    combine=id.lambda)
         },error=function(e) NULL)
       }
       outputNsample<-Filter(Negate(is.null), outputNsample)
       outputNsample<-Reduce(combine_lambda_mla,outputNsample)
     }else{
       for(id.lambda in 1:nlambda){
         if(id.lambda>1){
           
           #ids <- which(!vapply(outputNsample, is.null, logical(1)))
           ids <- which(unlist(lapply(outputNsample, FUN=function(x){
             if(is.null(x)){return(F)}else{
               if(x$istop==1){
                 return(T)
               }else{return(F)}
             }})))
           last_id <- if (length(ids) == 0) NA_integer_ else max(ids)
           if(!is.na(last_id)){
             beta.start<-outputNsample[[last_id]]$b[7:size_V]
             s.start<-outputNsample[[last_id]]$b[1:6]}
         } 
         outputNsample[[id.lambda]]<-tryCatch({
                                        
                                        # computation pbr 
                                        
                                        pbr_compu<-0
                                        
                                          beta<-beta.start
                                          s<-s.start
                                       
                                        
                                        converged<-F
                                        ite<-0
                                        # if beta not change do not need to recalculate weights 
                                        H<-T
                                        
                                        eval.cv.spline<-rep(NA,maxiter+1)
                                        eval.cv.beta<-rep(NA,maxiter+1)
                                        eval.cv.loglik<-rep(NA,maxiter+1)
                                        eval.loglik<-rep(NA,maxiter+1)
                                        eval.validity<-rep(NA,maxiter+1)
                                        
                                        
                                        while(converged==F & ite<=maxiter){
                                          
                                          
                                          b<-c(s,beta)
                                          bfix<-b[fix0==1]
                                          b<-b[fix0==0]
                                          # derivative of loglik
                                          
                                          
                                          output<-DYNderivaweibdiag( h=1e-8,
                                                                     b=b,
                                                                     npm=npm,
                                                                     npar=size_V,
                                                                     bfix=bfix,
                                                                     fix=fix0,
                                                                     ctime=ctime,
                                                                     no=N,
                                                                     ve01=ve01,
                                                                     ve02=ve02,
                                                                     ve12=ve12,
                                                                     dimnva01=dimnva01,
                                                                     dimnva02=dimnva02,
                                                                     dimnva12=dimnva12,
                                                                     nva01=nvat01,
                                                                     nva02=nvat02,
                                                                     nva12=nvat12,
                                                                     t0=t0,
                                                                     t1=t1,
                                                                     t2=t2,
                                                                     t3=t3,
                                                                     troncature=troncature,
                                                                     y01=y01k,
                                                                     y02=y02k,
                                                                     y12=y12k,
                                                                     p01=p01,
                                                                     p02=p02,
                                                                     p12=p12,
                                                                     dimp01=dimp01,
                                                                     dimp02=dimp02,
                                                                     dimp12=dimp12,
                                                                     Ntime=NtimePoints)
                                          output<-output$v
                                          
                                          if(ite==0){# loglik penalised
                                            fn.value<-gaussDYNidmlLikelihoodweibpena(b=b,
                                                                                     npm=npm,
                                                                                     npar=size_V,
                                                                                     bfix=bfix,
                                                                                     fix=fix0,
                                                                                     ctime=ctime,
                                                                                     no=N,
                                                                                     ve01=ve01,
                                                                                     ve02=ve02,
                                                                                     ve12=ve12,
                                                                                     dimnva01=dimnva01,
                                                                                     dimnva02=dimnva02,
                                                                                     dimnva12=dimnva12,
                                                                                     nva01=nvat01,
                                                                                     nva02=nvat02,
                                                                                     nva12=nvat12,
                                                                                     t0=t0,
                                                                                     t1=t1,
                                                                                     t2=t2,
                                                                                     t3=t3,
                                                                                     troncature=troncature,
                                                                                     y01=y01k,
                                                                                     y02=y02k,
                                                                                     y12=y12k,
                                                                                     p01=p01,
                                                                                     p02=p02,
                                                                                     p12=p12,
                                                                                     dimp01=dimp01,
                                                                                     dimp02=dimp02,
                                                                                     dimp12=dimp12,
                                                                                     Ntime=NtimePoints,
                                                                                     
                                                                                     lambda=lambda[id.lambda,],
                                                                                     alpha=alpha,
                                                                                     penalty.factor=penalty.factor,
                                                                                     penalty=penalty)
                                          }
                                          
                                          if(any(is.na(output))|any(output==Inf) |any(output==-Inf)){
                                            warning("Computational error for calculation of the hessian : division by 0 or Infinite value")
                                            if(ite==0){
                                              
                                              fu<-output[(min+1):length(output)]
                                              V<-matrix(0,nrow=npm,ncol=npm)
                                              diag(V)<-output[1:min]
                                              # deriva gives information matrix
                                              tr <- sum(diag(V))/npm
                                              V0<-V}
                                            ite<-ite+1
                                            pbr_compu<-1
                                            break
                                          }
                                          
                                          
                                          
                                          fu<-output[(min+1):length(output)]
                                          V<-matrix(0,nrow=npm,ncol=npm)
                                          diag(V)<-output[1:min]
                                          # deriva gives information matrix
                                          tr <- sum(diag(V))/npm
                                          V0<-V
                                          
                                          eigen.values<-diag(V)
                                          
                                          if(defpositive==T){
                                            idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                          }else{
                                            idpos<-ifelse(any(diag(V)==0),1,0)
                                          }
                                          
                                          
                                          idpos0<-idpos
                                          
                                          ncount<-da<-ga<-0
                                          
                                          while(idpos != 0){
                                            
                                            if(ncount==0){ 
                                              ga <- 0.01
                                              da <- 1E-2
                                            }else{
                                              if(((ncount <= 3) | (ga >= 1)) ){
                                                da <- da * 5
                                              }else{# if ncount > 10 only update ga 
                                                ga <- ga * 5
                                                # do not put ga at 1 as no countmax otherwise infinite while 
                                                if(ga > 1) ga <- 1
                                              }
                                            }
                                            
                                            ncount <- ncount + 1
                                            
                                            diagV <- diag(V)
                                            # put abs (1-ga) better than 1-ga cause ga can now be >1
                                            diagV<-ifelse(diagV!=0,diagV+da*(abs((1.e0-ga))*abs(diagV)+ga*tr),
                                                          da*ga*tr)
                                            
                                            diag(V)<-diagV
                                            # if we have a convex log-vraisemblance in eta then :
                                            # all eigen  values of the hessienne are >0.
                                            
                                            if(sum(V==Inf)>0|sum(V==-Inf)>0){break}
                                            eigen.values<-diag(V)
                                            # check if hessienne defined positive
                                            
                                            if(defpositive==T){
                                              idpos<-ifelse(any(eigen.values<=eps.eigen),1,0)
                                            }else{
                                              idpos<-ifelse(any(diag(V)==0),1,0)
                                            }
                                            
                                            # if(def.positive==T){
                                            #   idpos<-ifelse(any(eigen.values<=0),1,0)
                                            # }else{idpos<-ifelse(any(abs(eigen.values)==0),1,0)}
                                            
                                          }
                                          
                                          if(idpos!=0){
                                            
                                            warning("Hessian not defined positive")
                                            pbr_compu<-2
                                            ite<-ite+1
                                            break
                                          }
                                          
                                          
                                          
                                          
                                          
                                          output.cv<-DYNcv.model(beta=beta,
                                                                 nva01=npm01,
                                                                 nva02=npm02,
                                                                 nva12=npm12,
                                                                 nva01Y=npm01Y,
                                                                 nva02Y=npm02Y,
                                                                 nva12Y=npm12Y,
                                                                 fix=fix0[7:size_V],
                                                                 penalty.factor=penalty.factor,
                                                                 penalty=penalty,
                                                                 v=V,
                                                                 fu=fu,
                                                                 lambda=lambda[id.lambda,],
                                                                 alpha=alpha
                                          )
                                          
                                          # verify validity of parameters update 
                                          # and that we are better than previous estimates 
                                          
                                          b<-c(s,output.cv$b)
                                          
                                          betanew<-b[(6+1):size_V]
                                          
                                          # penalised loglik see if inferior to previous
                                          res<-gaussDYNidmlLikelihoodweibpena(b=b,
                                                                              npm=size_V,
                                                                              npar=size_V,
                                                                              bfix=1,
                                                                              fix=rep(0,size_V),
                                                                              ctime=ctime,
                                                                              no=N,
                                                                              ve01=ve01,
                                                                              ve02=ve02,
                                                                              ve12=ve12,
                                                                              dimnva01=dimnva01,
                                                                              dimnva02=dimnva02,
                                                                              dimnva12=dimnva12,
                                                                              nva01=nvat01,
                                                                              nva02=nvat02,
                                                                              nva12=nvat12,
                                                                              t0=t0,
                                                                              t1=t1,
                                                                              t2=t2,
                                                                              t3=t3,
                                                                              troncature=troncature,
                                                                              y01=y01k,
                                                                              y02=y02k,
                                                                              y12=y12k,
                                                                              p01=p01,
                                                                              p02=p02,
                                                                              p12=p12,
                                                                              dimp01=dimp01,
                                                                              dimp02=dimp02,
                                                                              dimp12=dimp12,
                                                                              Ntime=NtimePoints,
                                                                              
                                                                              lambda=lambda[id.lambda,],
                                                                              alpha=alpha,
                                                                              penalty.factor=penalty.factor,
                                                                              penalty=penalty)
                                          
                                          
                                          # we want to maximise the loglik thus : 
                                          # we have issue if res is NA or if not higher than previous one 
                                          # if not better or do not exist need to readjust
                                          # value of beta 
                                          if(res %in%c(-1e9,1e9) | res < fn.value){
                                            
                                            th<-1e-5
                                            step<-log(1.5)
                                            delta<-output.cv$b-c(beta)
                                            
                                            maxt <- max(abs(delta)) 
                                            
                                            if(maxt == 0){
                                              vw <- th
                                            }else{
                                              vw <- th/maxt
                                            }
                                            if(ite>0){
                                              res.out.error <- list("old.b"=round(c(s,beta)),
                                                                    "old.rl"=round(fn.value),
                                                                    "old.ca"=round(eval.cv.beta[ite]),
                                                                    "old.cb"=round(eval.cv.loglik[ite]))
                                            }else{
                                              res.out.error <- list("old.b"=round(c(s,beta)),
                                                                    "old.rl"=round(fn.value),
                                                                    "old.ca"=round(1),
                                                                    "old.cb"=round(1))
                                            }
                                            # from mla package 
                                            sears<-searpas(vw=vw,
                                                           step=step,
                                                           b=beta,
                                                           delta=delta,
                                                           funcpa=gaussDYNidmlLikelihoodweibpena,
                                                           res.out.error=res.out.error,
                                                           npm=npm,
                                                           npar=size_V,
                                                           bfix=s,
                                                           fix=fix0,
                                                           ctime=ctime,
                                                           no=N,
                                                           ve01=ve01,
                                                           ve02=ve02,
                                                           ve12=ve12,
                                                           dimnva01=dimnva01,
                                                           dimnva02=dimnva02,
                                                           dimnva12=dimnva12,
                                                           nva01=nvat01,
                                                           nva02=nvat02,
                                                           nva12=nvat12,
                                                           t0=t0,
                                                           t1=t1,
                                                           t2=t2,
                                                           t3=t3,
                                                           troncature=troncature,
                                                           y01=y01k,
                                                           y02=y02k,
                                                           y12=y12k,
                                                           p01=p01,
                                                           p02=p02,
                                                           p12=p12,
                                                           dimp01=dimp01,
                                                           dimp02=dimp02,
                                                           dimp12=dimp12,
                                                           Ntime=NtimePoints,
                                                           lambda=lambda[id.lambda,],
                                                           alpha=alpha,
                                                           penalty.factor=penalty.factor,
                                                           penalty=penalty)
                                            
                                            betanew<-beta+delta*sears$vw
                                            betanew<-ifelse(abs(betanew)<=0.0001,0,betanew)
                                            b<-c(s,betanew)
                                            
                                            
                                            res<-gaussDYNidmlLikelihoodweibpena(b=b,
                                                                                npm=size_V,
                                                                                npar=size_V,
                                                                                bfix=1,
                                                                                fix=rep(0,size_V),
                                                                                ctime=ctime,
                                                                                no=N,
                                                                                ve01=ve01,
                                                                                ve02=ve02,
                                                                                ve12=ve12,
                                                                                dimnva01=dimnva01,
                                                                                dimnva02=dimnva02,
                                                                                dimnva12=dimnva12,
                                                                                nva01=nvat01,
                                                                                nva02=nvat02,
                                                                                nva12=nvat12,
                                                                                t0=t0,
                                                                                t1=t1,
                                                                                t2=t2,
                                                                                t3=t3,
                                                                                troncature=troncature,
                                                                                y01=y01k,
                                                                                y02=y02k,
                                                                                y12=y12k,
                                                                                p01=p01,
                                                                                p02=p02,
                                                                                p12=p12,
                                                                                dimp01=dimp01,
                                                                                dimp02=dimp02,
                                                                                dimp12=dimp12,
                                                                                Ntime=NtimePoints,
                                                                                
                                                                                lambda=lambda[id.lambda,],
                                                                                alpha=alpha,
                                                                                penalty.factor=penalty.factor,
                                                                                penalty=penalty)
                                          }
                                          # if not better or do not exist need to readjust
                                          # value of beta 
                                          if(res %in%c(-1e9,1e9) | any(is.infinite(c(s,betanew)))){
                                            
                                            ite<-ite+1
                                            validity<-F
                                            eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                            eval.validity[ite]<-validity
                                            # save output over iterations 
                                            
                                            pbr_compu<-3
                                            break
                                          }else{validity<-T}
                                          
                                          
                                          # betanew already include s
                                          b<-c(s,betanew)
                                          
                                          bfix<-b[fix0.beta==1]
                                          b<-b[fix0.beta==0]
                                          # update for modelPar
                                          if(npmweib!=0){
                                            output.mla<- marqLevAlg::mla(b=b,
                                                                         fn=gaussDYNidmlLikelihoodweib,
                                                                         gr=DYNreggrmlaweibana,
                                                                         epsa=epsa,
                                                                         epsb=epsb,
                                                                         epsd=epsd,
                                                                         maxiter=maxiter.pena,
                                                                         minimize=F,
                                                                         npm=npmweib,
                                                                         npar=size_V,
                                                                         bfix=bfix,
                                                                         fix=fix0.beta,
                                                                         ctime=ctime,
                                                                         no=N,
                                                                         ve01=ve01,
                                                                         ve02=ve02,
                                                                         ve12=ve12,
                                                                         dimnva01=dimnva01,
                                                                         dimnva02=dimnva02,
                                                                         dimnva12=dimnva12,
                                                                         nva01=nvat01,
                                                                         nva02=nvat02,
                                                                         nva12=nvat12,
                                                                         t0=t0,
                                                                         t1=t1,
                                                                         t2=t2,
                                                                         t3=t3,
                                                                         troncature=troncature,
                                                                         y01=y01k,
                                                                         y02=y02k,
                                                                         y12=y12k,
                                                                         p01=p01,
                                                                         p02=p02,
                                                                         p12=p12,
                                                                         dimp01=dimp01,
                                                                         dimp02=dimp02,
                                                                         dimp12=dimp12,
                                                                         Ntime=NtimePoints)
                                            
                                            # look at convergence for each lambda :
                                            # mla output is loglik
                                            # new values for splines:
                                            snew<-s
                                            snew[fix00[1:6]==0]<-output.mla$b
                                            if(nvat01>0){
                                              b01<-betanew[1:nvat01][penalty.factor[1:nvat01]==1]
                                              if(p01>0){
                                                b01<-c(b01,betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1])
                                              }
                                            }else{
                                              if(p01>0){
                                                b01<-betanew[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)][penalty.factor[(nvat01+nvat02+nvat12+1):(nvat01+nvat02+nvat12+p01)]==1]
                                              }else{
                                                b01<-0
                                              }
                                            }
                                            
                                            if(nvat02>0){
                                              b02<-betanew[(nvat01+1):(nvat01+nvat02)][penalty.factor[(nvat01+1):(nvat01+nvat02)]==1]
                                              if(p02>0){
                                                b02<-c(b02,betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1])
                                              }
                                            }else{
                                              if(p02>0){
                                                b02<-betanew[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)][penalty.factor[(nvat01+nvat02+nvat12+p01+1):(nvat01+nvat02+nvat12+p01+p02)]==1]
                                              }else{b02<-0}
                                            }
                                            
                                            if(nvat12>0){
                                              b12<-betanew[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)][penalty.factor[(nvat01+nvat02+1):(nvat01+nvat02+nvat12)]==1]
                                              if(p12>0){
                                                b12<-c(b12,betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1])
                                              }
                                            }else{
                                              if(p12>0){
                                                b12<-betanew[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)][penalty.factor[(nvat01+nvat02+nvat12+p01+p02+1):(nvat01+nvat02+nvat12+p01+p02+p12)]==1]
                                              }else{b12<-0}
                                            }
                                            
                                            # calculate loglik pen 
                                            if(penalty%in%c("lasso","ridge","elasticnet","corrected.elasticnet")){
                                              fn.valuenew<-output.mla$fn.value-lambda[id.lambda,1]*alpha*sum(abs(b01))-lambda[id.lambda,1]*(1-alpha)*sum(b01*b01)
                                              fn.valuenew<-fn.valuenew-lambda[id.lambda,2]*alpha*sum(abs(b02))-lambda[id.lambda,2]*(1-alpha)*sum(b02*b02)
                                              fn.valuenew<-fn.valuenew-lambda[id.lambda,3]*alpha*sum(abs(b12))-lambda[id.lambda,3]*(1-alpha)*sum(b12*b12)
                                            }
                                            
                                            
                                            if(penalty=="mcp"){
                                              
                                              p01<-rep(alpha*lambda[id.lambda,1]*lambda[id.lambda,1]/2,length(b01))
                                              idbeta<-which(b01<=alpha*lambda[id.lambda,1])
                                              p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])-((b01[idbeta]*b01[idbeta])/2*alpha)
                                              
                                              p02<-rep(alpha*lambda[id.lambda,2]*lambda[id.lambda,2]/2,length(b02))
                                              idbeta<-which(b02<=alpha*lambda[id.lambda,2])
                                              p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])-((b02[idbeta]*b02[idbeta])/2*alpha)
                                              
                                              p12<-rep(alpha*lambda[id.lambda,3]*lambda[id.lambda,3]/2,length(b12))
                                              idbeta<-which(b12<=alpha*lambda[id.lambda,3])
                                              p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])-((b12[idbeta]*b12[idbeta])/2*alpha)
                                              
                                              fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                              
                                            }
                                            
                                            if(penalty=="scad"){
                                              
                                              p01<-rep((lambda[id.lambda,1]^2)*(alpha+1)/2,length(b01))
                                              idbeta<-which(b01<=lambda[id.lambda,1])
                                              p01[idbeta]<-lambda[id.lambda,1]*abs(b01[idbeta])
                                              idbeta<-which(abs(b01)<lambda[id.lambda,1]*alpha)
                                              p01[idbeta]<-(2*alpha*lambda[id.lambda,1]*abs(b01[idbeta])-b01[idbeta]^2-lambda[id.lambda,1]^2)/(2*(alpha-1))
                                              
                                              p02<-rep((lambda[id.lambda,2]^2)*(alpha+1)/2,length(b02))
                                              idbeta<-which(b02<=lambda[id.lambda,2])
                                              p02[idbeta]<-lambda[id.lambda,2]*abs(b02[idbeta])
                                              idbeta<-which(abs(b02)<lambda[id.lambda,2]*alpha)
                                              p02[idbeta]<-(2*alpha*lambda[id.lambda,2]*abs(b02[idbeta])-b02[idbeta]^2-lambda[id.lambda,2]^2)/(2*(alpha-1))
                                              
                                              p12<-rep((lambda[id.lambda,3]^2)*(alpha+1)/2,length(b12))
                                              idbeta<-which(b12<=lambda[id.lambda,3])
                                              p12[idbeta]<-lambda[id.lambda,3]*abs(b12[idbeta])
                                              idbeta<-which(abs(b12)<lambda[id.lambda,3]*alpha)
                                              p12[idbeta]<-(2*alpha*lambda[id.lambda,3]*abs(b12[idbeta])-b12[idbeta]^2-lambda[id.lambda,3]^2)/(2*(alpha-1))
                                              
                                              fn.valuenew<-output.mla$fn.value-sum(p01)-sum(p02)-sum(p12)
                                              
                                            }
                                          }else{
                                            snew<-s
                                            fn.valuenew<-res
                                          }
                                          
                                          ite<-ite+1
                                          
                                          #check cv 
                                          eval.cv.spline[ite]<-sum((snew-s)^2)
                                          eval.cv.beta[ite]<-sum((betanew-beta)^2)
                                          eval.cv.loglik[ite]<-abs((fn.valuenew-fn.value)/fn.value)
                                          eval.loglik[ite]<-fn.valuenew
                                          eval.validity[ite]<-validity
                                          
                                          s<-snew
                                          beta<-betanew
                                          fn.value<-fn.valuenew
                                          
                                          # eval.cv beta valid only if validity.param=T
                                          if(((eval.cv.beta[ite] + eval.cv.spline[ite])<epsa) & eval.cv.loglik[ite]<epsb & validity==T){
                                            converged<-T}
                                          
                                          
                                        }
                                        
                                        if(maxiter<=ite & converged==F){
                                          istop<-2
                                        }else{
                                          if(ite<=maxiter & converged==T){
                                            istop<-1 
                                            
                                            # need to recalculate second derivatives 
                                            # if converged 
                                            
                                            b<-c(s,beta)
                                            bfix<-b[fix0==1]
                                            b<-b[fix0==0]
                                            
                                            output<-DYNderivaweibdiag( h=1e-8,
                                                                       b=b,
                                                                       npm=npm,
                                                                       npar=size_V,
                                                                       bfix=bfix,
                                                                       fix=fix0,
                                                                       ctime=ctime,
                                                                       no=N,
                                                                       ve01=ve01,
                                                                       ve02=ve02,
                                                                       ve12=ve12,
                                                                       dimnva01=dimnva01,
                                                                       dimnva02=dimnva02,
                                                                       dimnva12=dimnva12,
                                                                       nva01=nvat01,
                                                                       nva02=nvat02,
                                                                       nva12=nvat12,
                                                                       t0=t0,
                                                                       t1=t1,
                                                                       t2=t2,
                                                                       t3=t3,
                                                                       troncature=troncature,
                                                                       y01=y01k,
                                                                       y02=y02k,
                                                                       y12=y12k,
                                                                       p01=p01,
                                                                       p02=p02,
                                                                       p12=p12,
                                                                       dimp01=dimp01,
                                                                       dimp02=dimp02,
                                                                       dimp12=dimp12,
                                                                       Ntime=NtimePoints)
                                            
                                            output<-output$v
                                            fu<-output[(npm+1):length(output)]
                                            V<-matrix(0,nrow=npm,ncol=npm)
                                            diag(V)<-output[1:npm]
                                            
                                            # deriva gives information matrix
                                            V0<-V
                                          }else{
                                            if(pbr_compu==1){istop<-3}
                                            if(pbr_compu==2){istop<-4}
                                            if(pbr_compu==3){istop<-5}
                                          }
                                        }
                                        
                                        # if stop==1 we can give matrix of second derivatives 
                                        
                                        
                                        list(b=c(s,beta),
                                                    H=V0,
                                                    lambda=as.double(lambda[id.lambda,]),
                                                    alpha=alpha,
                                                    fn.value=ifelse(!exists("output.mla"),NA,output.mla$fn.value), # loglik
                                                    fn.value.pena=fn.value, # penalised loglik
                                                    ni=ite,
                                                    ca.beta=eval.cv.beta,
                                                    ca.spline=eval.cv.spline,
                                                    ca.validity=eval.validity,
                                                    cb=eval.loglik,
                                                    istop=istop,
                                                    combine=id.lambda)
                                      },error=function(e) NULL)
     }
    outputNsample<-Filter(Negate(is.null), outputNsample)
    outputNsample<-Reduce(combine_lambda_mla,outputNsample)
  }
    }
  
  return(outputNsample)
  
}