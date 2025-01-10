 



!============================================================================================= 
!========================          idmlLikelihood         ====================================
!========================   with baseline M-splines       ==================================== 
!============================================================================================= 

      subroutine idmlikelihood(b0,np0,npar0,bfix0,fix0,zi010,zi120,zi020,c0,&
      no0,nz010,nz120,nz020,ve010,ve120,ve020,&
        dimnva01,dimnva12,dimnva02,nva01,nva12,nva02,t00,t10,t20,&
        t30,troncature0,gausspoint0,likelihood_res)

	use commun
        implicit none
         
        double precision::res,res1,res2,tronc, &
        vet01,vet12,vet02
	double precision, intent(inout)::likelihood_res
        integer::np0,i,j,l,w,k,npar0,nva01,nva12,nva02,no0, &
	nz010,nz020,nz120,troncature0,gausspoint0, & 
        dimnva01,dimnva02,dimnva12

	double precision,dimension(np0)::b0
        double precision,dimension(npar0)::bh
	double precision,dimension(npar0-np0)::bfix0
	integer,dimension(npar0)::fix0
	double precision,dimension(-2:(nz010+3))::zi010
	double precision,dimension(-2:(nz020+3))::zi020
	double precision,dimension(-2:(nz120+3))::zi120
	double precision,dimension(-2:(nz010-1))::the01
	double precision,dimension(-2:(nz120-1))::the12
	double precision,dimension(-2:(nz020-1))::the02
        double precision,dimension(no0,dimnva01)::ve010
	double precision,dimension(no0,dimnva02)::ve020
	double precision,dimension(no0,dimnva12)::ve120
	
	
        double precision::su01,ri01,su12,ri12,su02,ri02,gl01,gl02,gl12
	double precision,dimension(no0)::t00,t10,t20,t30
	integer,dimension(no0)::c0

	allocate(b(np0),bfix(npar0-np0),fix(npar0))
	b=b0
	bfix=bfix0
	fix=fix0
	allocate(zi01(-2:(nz01+3)),zi12(-2:(nz12+3)),zi02(-2:(nz02+3)))
	zi01=zi010
	zi02=zi020
	zi12=zi120

	
	nz01=nz010
	nz02=nz020
	nz12=nz120
	troncature=troncature0
	gausspoint=gausspoint0


	if(nva01.gt.0) then 
		allocate(ve01(no0,nva01))
	else 
		allocate(ve01(no0,1))
	end if 
	
	if(nva02.gt.0) then 
		allocate(ve02(no0,nva02))
	else 
		allocate(ve02(no0,1))
	end if 

	if(nva12.gt.0) then 
		allocate(ve12(no0,nva12))
	else 
		allocate(ve12(no0,1))
	end if 


	ve01=ve010
	ve02=ve020
	ve12=ve120

	allocate(t0(no0),t1(no0),t2(no0),t3(no0),c(no0))
	c=c0
	t0=t00
	t1=t10
	t2=t20
	t3=t30

         
        ! we need to put bh at its original values if in posfix 


       l=0
       w=0


       do k=1,(np0+sum(fix))
         if(fix(k).eq.0) then
            l=l+1
            bh(k)=b(l)
         end if
         if(fix(k).eq.1) then
            w=w+1
            bh(k)=bfix(w)
         end if
      end do
    
	

         do i=1,nz01+2
            the01(i-3)=(bh(i))*(bh(i))
!       the01(i-3)=dexp(bh(i))
         end do
         do i=1,nz02+2
            j = nz01+2+i
            the02(i-3)=(bh(j))*(bh(j))
!       the12(i-3)=dexp(bh(j))
         end do
         do i=1,nz12+2
            j = nz02+2+nz01+2+i
            the12(i-3)=(bh(j))*(bh(j))
!       the02(i-3)=dexp(bh(j))
         end do
	
	
!---------- calcul de la vraisemblance ------------------



        res = 0.d0
	if(gausspoint.eq.10)then
        do i=1,no0

                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
                                vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
                                vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
                                vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)


                res1 = 0.d0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else 
                                call susp(t0(i),the01,nz01,su01,ri01,zi01,gl01) 
                                call susp(t0(i),the02,nz02,su02,ri02,zi02,gl02) 
                                tronc=(gl01*vet01)+(gl02*vet02)
                        end if
                else
                        tronc = 0.d0
                end if

		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
                        call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
                        call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)

                        res1 = (-gl01*vet01)-(gl02*vet02)
                else
                if(c(i).eq.2)then ! cpi 0-->1
                       
                        	call qgaussPL(t1(i),t2(i),the01,the12,&
                        	the02,res2,vet01,vet12,vet02)
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
!               res1=dLOG(res2)-gl12*vet12 (autre ecriture)
		                      res1=dLOG(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
                    call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                    call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                    call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                    res2=-gl01*vet01+dLOG(ri01*vet01) -&
                    gl12*vet12-gl02*vet02
                    call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                    res1=res2 +gl12*vet12
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                          call qgaussPL(t1(i),t2(i),the01,the12,the02,&
                        	res2,vet01,vet12,vet02)
                          res1=dLOG(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
                         call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                         call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                         res2=-gl01*vet01 + dLOG(ri01*vet01) -&
                         gl12*vet12 + dLOG(ri12*vet12) -gl02*vet02
                         call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                         res1=res2+gl12*vet12

                         else
                            if(c(i).eq.6)then ! vivant ???
                          		call qgaussPL(t1(i),t2(i),the01,the12,the02,&
                              res2,vet01,vet12,vet02)
                              call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                              call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                              call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                              res1=dLOG(res2*(su12**vet12) +&
                            (su01**vet01)*(su02**vet02))  

                            else ! passage 0-->2  

                                call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                                call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                        	    	call qgaussPL(t1(i),t3(i),the01,the12,&
                        		the02,res2,vet01,vet12,vet02)
                                res1=dLOG(res2*ri12*vet12*(su12**vet12) +&
				ri02*vet02*(su02**vet02)*(su01**vet01))
                            endif
                          endif  
                        endif    
                    endif
                endif   
             endif   

                res = res + res1  + tronc
             
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        likelihood_res=-1.d9
                        goto 123
                end if
                
          end do 
    else 
        if(gausspoint.eq.15) then 
          do i=1,no0

                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
                                vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
                                vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
                                vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)


                res1 = 0.d0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else 
                                call susp(t0(i),the01,nz01,su01,ri01,zi01,gl01) 
                                call susp(t0(i),the02,nz02,su02,ri02,zi02,gl02) 
                                tronc=(gl01*vet01)+(gl02*vet02)
                        end if
                else
                        tronc = 0.d0
                end if

		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
                        call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
                        call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)

                        res1 = (-gl01*vet01)-(gl02*vet02)
                else
                if(c(i).eq.2)then ! cpi 0-->1
                       
                        	call qgaussPL15(t1(i),t2(i),the01,the12,&
                        	the02,res2,vet01,vet12,vet02)
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
!               res1=dLOG(res2)-gl12*vet12 (autre ecriture)
		                      res1=dLOG(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
                    call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                    call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                    call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                    res2=-gl01*vet01+dLOG(ri01*vet01) -&
                    gl12*vet12-gl02*vet02
                    call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                    res1=res2 +gl12*vet12
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                          call qgaussPL15(t1(i),t2(i),the01,the12,the02,&
                        	res2,vet01,vet12,vet02)
                          res1=dLOG(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
                         call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                         call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                         res2=-gl01*vet01 + dLOG(ri01*vet01) -&
                         gl12*vet12 + dLOG(ri12*vet12) -gl02*vet02
                         call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                         res1=res2+gl12*vet12

                         else
                            if(c(i).eq.6)then ! vivant ???
                          		call qgaussPL15(t1(i),t2(i),the01,the12,the02,&
                              res2,vet01,vet12,vet02)
                              call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                              call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                              call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                              res1=dLOG(res2*(su12**vet12) +&
                            (su01**vet01)*(su02**vet02))  

                            else ! passage 0-->2  

                                call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                                call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                        	    	call qgaussPL15(t1(i),t3(i),the01,the12,&
                        		the02,res2,vet01,vet12,vet02)
                                res1=dLOG(res2*ri12*vet12*(su12**vet12) +&
				ri02*vet02*(su02**vet02)*(su01**vet01))
                            endif
                          endif  
                        endif    
                    endif
                endif   
             endif   

                res = res + res1  + tronc
             
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        likelihood_res=-1.d9
                        goto 123
                end if
                
          end do 
          
      else 
            if(gausspoint.eq.21) then 
            
                do i=1,no0

                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
                                vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
                                vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
                                vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)


                res1 = 0.d0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else 
                                call susp(t0(i),the01,nz01,su01,ri01,zi01,gl01) 
                                call susp(t0(i),the02,nz02,su02,ri02,zi02,gl02) 
                                tronc=(gl01*vet01)+(gl02*vet02)
                        end if
                else
                        tronc = 0.d0
                end if

		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
                        call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
                        call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)

                        res1 = (-gl01*vet01)-(gl02*vet02)
                else
                if(c(i).eq.2)then ! cpi 0-->1
                       
                        	call qgaussPL21(t1(i),t2(i),the01,the12,&
                        	the02,res2,vet01,vet12,vet02)
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
!               res1=dLOG(res2)-gl12*vet12 (autre ecriture)
		                      res1=dLOG(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
                    call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                    call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                    call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                    res2=-gl01*vet01+dLOG(ri01*vet01) -&
                    gl12*vet12-gl02*vet02
                    call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                    res1=res2 +gl12*vet12
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                          call qgaussPL21(t1(i),t2(i),the01,the12,the02,&
                        	res2,vet01,vet12,vet02)
                          res1=dLOG(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
                         call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                         call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                         res2=-gl01*vet01 + dLOG(ri01*vet01) -&
                         gl12*vet12 + dLOG(ri12*vet12) -gl02*vet02
                         call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                         res1=res2+gl12*vet12

                         else
                            if(c(i).eq.6)then ! vivant ???
                          		call qgaussPL21(t1(i),t2(i),the01,the12,the02,&
                              res2,vet01,vet12,vet02)
                              call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                              call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                              call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                              res1=dLOG(res2*(su12**vet12) +&
                            (su01**vet01)*(su02**vet02))  

                            else ! passage 0-->2  

                                call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                                call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                        	    	call qgaussPL21(t1(i),t3(i),the01,the12,&
                        		the02,res2,vet01,vet12,vet02)
                                res1=dLOG(res2*ri12*vet12*(su12**vet12) +&
				ri02*vet02*(su02**vet02)*(su01**vet01))
                            endif
                          endif  
                        endif    
                    endif
                endif   
             endif   

                res = res + res1  + tronc
             
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        likelihood_res=-1.d9
                        goto 123
                end if
                
          end do 
          
      else 
            if(gausspoint.eq.31) then 
            
                do i=1,no0

                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
                                vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
                                vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
                                vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)


                res1 = 0.d0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else 
                                call susp(t0(i),the01,nz01,su01,ri01,zi01,gl01) 
                                call susp(t0(i),the02,nz02,su02,ri02,zi02,gl02) 
                                tronc=(gl01*vet01)+(gl02*vet02)
                        end if
                else
                        tronc = 0.d0
                end if

		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
                        call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
                        call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)

                        res1 = (-gl01*vet01)-(gl02*vet02)
                else
                if(c(i).eq.2)then ! cpi 0-->1
                       
                        	call qgaussPL31(t1(i),t2(i),the01,the12,&
                        	the02,res2,vet01,vet12,vet02)
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
!               res1=dLOG(res2)-gl12*vet12 (autre ecriture)
		                      res1=dLOG(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
                    call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                    call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                    call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                    res2=-gl01*vet01+dLOG(ri01*vet01) -&
                    gl12*vet12-gl02*vet02
                    call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                    res1=res2 +gl12*vet12
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                          call qgaussPL31(t1(i),t2(i),the01,the12,the02,&
                        	res2,vet01,vet12,vet02)
                          res1=dLOG(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
                         call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                         call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                         res2=-gl01*vet01 + dLOG(ri01*vet01) -&
                         gl12*vet12 + dLOG(ri12*vet12) -gl02*vet02
                         call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                         res1=res2+gl12*vet12

                         else
                            if(c(i).eq.6)then ! vivant ???
                          		call qgaussPL31(t1(i),t2(i),the01,the12,the02,&
                              res2,vet01,vet12,vet02)
                              call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                              call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                              call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                              res1=dLOG(res2*(su12**vet12) +&
                            (su01**vet01)*(su02**vet02))  

                            else ! passage 0-->2  

                                call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                                call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                        	    	call qgaussPL31(t1(i),t3(i),the01,the12,&
                        		the02,res2,vet01,vet12,vet02)
                                res1=dLOG(res2*ri12*vet12*(su12**vet12) +&
				ri02*vet02*(su02**vet02)*(su01**vet01))
                            endif
                          endif  
                        endif    
                    endif
                endif   
             endif   

                res = res + res1  + tronc
             
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        likelihood_res=-1.d9
                        goto 123
                end if
                
          end do 
          
      else 
            if(gausspoint.eq.41) then 
            
            do i=1,no0

                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
                                vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
                                vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
                                vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)


                res1 = 0.d0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else 
                                call susp(t0(i),the01,nz01,su01,ri01,zi01,gl01) 
                                call susp(t0(i),the02,nz02,su02,ri02,zi02,gl02) 
                                tronc=(gl01*vet01)+(gl02*vet02)
                        end if
                else
                        tronc = 0.d0
                end if

		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
                        call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
                        call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)

                        res1 = (-gl01*vet01)-(gl02*vet02)
                else
                if(c(i).eq.2)then ! cpi 0-->1
                       
                        	call qgaussPL41(t1(i),t2(i),the01,the12,&
                        	the02,res2,vet01,vet12,vet02)
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
!               res1=dLOG(res2)-gl12*vet12 (autre ecriture)
		                      res1=dLOG(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
                    call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                    call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                    call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                    res2=-gl01*vet01+dLOG(ri01*vet01) -&
                    gl12*vet12-gl02*vet02
                    call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                    res1=res2 +gl12*vet12
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                          call qgaussPL41(t1(i),t2(i),the01,the12,the02,&
                        	res2,vet01,vet12,vet02)
                          res1=dLOG(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
                         call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                         call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                         res2=-gl01*vet01 + dLOG(ri01*vet01) -&
                         gl12*vet12 + dLOG(ri12*vet12) -gl02*vet02
                         call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                         res1=res2+gl12*vet12

                         else
                            if(c(i).eq.6)then ! vivant ???
                          		call qgaussPL41(t1(i),t2(i),the01,the12,the02,&
                              res2,vet01,vet12,vet02)
                              call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                              call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                              call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                              res1=dLOG(res2*(su12**vet12) +&
                            (su01**vet01)*(su02**vet02))  

                            else ! passage 0-->2  

                                call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                                call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                        	    	call qgaussPL41(t1(i),t3(i),the01,the12,&
                        		the02,res2,vet01,vet12,vet02)
                                res1=dLOG(res2*ri12*vet12*(su12**vet12) +&
				ri02*vet02*(su02**vet02)*(su01**vet01))
                            endif
                          endif  
                        endif    
                    endif
                endif   
             endif   

                res = res + res1  + tronc
             
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        likelihood_res=-1.d9
                        goto 123
                end if
                
          end do 

    else 
          if(gausspoint.eq.51) then 
                  
                  do i=1,no0

                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
                                vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
                                vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
                                vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)


                res1 = 0.d0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else 
                                call susp(t0(i),the01,nz01,su01,ri01,zi01,gl01) 
                                call susp(t0(i),the02,nz02,su02,ri02,zi02,gl02) 
                                tronc=(gl01*vet01)+(gl02*vet02)
                        end if
                else
                        tronc = 0.d0
                end if

		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
                        call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
                        call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)

                        res1 = (-gl01*vet01)-(gl02*vet02)
                else
                if(c(i).eq.2)then ! cpi 0-->1
                       
                        	call qgaussPL51(t1(i),t2(i),the01,the12,&
                        	the02,res2,vet01,vet12,vet02)
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
!               res1=dLOG(res2)-gl12*vet12 (autre ecriture)
		                      res1=dLOG(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
                    call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                    call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                    call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                    res2=-gl01*vet01+dLOG(ri01*vet01) -&
                    gl12*vet12-gl02*vet02
                    call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                    res1=res2 +gl12*vet12
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                          call qgaussPL51(t1(i),t2(i),the01,the12,the02,&
                        	res2,vet01,vet12,vet02)
                          res1=dLOG(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
                         call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                         call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                         res2=-gl01*vet01 + dLOG(ri01*vet01) -&
                         gl12*vet12 + dLOG(ri12*vet12) -gl02*vet02
                         call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                         res1=res2+gl12*vet12

                         else
                            if(c(i).eq.6)then ! vivant ???
                          		call qgaussPL51(t1(i),t2(i),the01,the12,the02,&
                              res2,vet01,vet12,vet02)
                              call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                              call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                              call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                              res1=dLOG(res2*(su12**vet12) +&
                            (su01**vet01)*(su02**vet02))  

                            else ! passage 0-->2  

                                call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                                call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                        	    	call qgaussPL51(t1(i),t3(i),the01,the12,&
                        		the02,res2,vet01,vet12,vet02)
                                res1=dLOG(res2*ri12*vet12*(su12**vet12) +&
				ri02*vet02*(su02**vet02)*(su01**vet01))
                            endif
                          endif  
                        endif    
                    endif
                endif   
             endif   

                res = res + res1  + tronc
             
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        likelihood_res=-1.d9
                        goto 123
                end if
                
          end do 
          
     else 
     
         do i=1,no0

                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
                                vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
                                vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
                                vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)


                res1 = 0.d0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else 
                                call susp(t0(i),the01,nz01,su01,ri01,zi01,gl01) 
                                call susp(t0(i),the02,nz02,su02,ri02,zi02,gl02) 
                                tronc=(gl01*vet01)+(gl02*vet02)
                        end if
                else
                        tronc = 0.d0
                end if

		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
                        call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
                        call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)

                        res1 = (-gl01*vet01)-(gl02*vet02)
                else
                if(c(i).eq.2)then ! cpi 0-->1
                       
                        	call qgaussPL61(t1(i),t2(i),the01,the12,&
                        	the02,res2,vet01,vet12,vet02)
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
!               res1=dLOG(res2)-gl12*vet12 (autre ecriture)
		                      res1=dLOG(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
                    call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                    call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                    call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                    res2=-gl01*vet01+dLOG(ri01*vet01) -&
                    gl12*vet12-gl02*vet02
                    call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                    res1=res2 +gl12*vet12
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                          call qgaussPL61(t1(i),t2(i),the01,the12,the02,&
                        	res2,vet01,vet12,vet02)
                          res1=dLOG(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
                         call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                         call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                         res2=-gl01*vet01 + dLOG(ri01*vet01) -&
                         gl12*vet12 + dLOG(ri12*vet12) -gl02*vet02
                         call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                         res1=res2+gl12*vet12

                         else
                            if(c(i).eq.6)then ! vivant ???
                          		call qgaussPL61(t1(i),t2(i),the01,the12,the02,&
                              res2,vet01,vet12,vet02)
                              call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                              call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                              call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                              res1=dLOG(res2*(su12**vet12) +&
                            (su01**vet01)*(su02**vet02))  

                            else ! passage 0-->2  

                                call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                                call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                        	    	call qgaussPL61(t1(i),t3(i),the01,the12,&
                        		the02,res2,vet01,vet12,vet02)
                                res1=dLOG(res2*ri12*vet12*(su12**vet12) +&
				ri02*vet02*(su02**vet02)*(su01**vet01))
                            endif
                          endif  
                        endif    
                    endif
                endif   
             endif   

                res = res + res1  + tronc
             
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        likelihood_res=-1.d9
                        goto 123
                end if
                
          end do 
            
            
            end if 
            end if 
            end if 
            end if
            end if
            end if

        likelihood_res = res


123     continue 
	 
	deallocate(b,bfix,fix,zi01,zi02,zi12,ve01,ve02,ve12,t0,t1,t2,t3,c)

        end subroutine idmlikelihood



!============================================================================================= 
!========================          idmlLikelihood         ====================================
!========================    with baseline weibull        ==================================== 
!============================================================================================= 


      subroutine idmlikelihoodweib(b0,np0,npar0,bfix0,fix0,c0,&
      no0,ve010,ve120,ve020,dimnva01,dimnva12,dimnva02,nva01,&
      nva12,nva02,t00,t10,t20,t30,troncature0,gausspoint0,likelihood_res)

	use commun
        implicit none
         
        double precision::res,res1,res2,tronc, &
        vet01,vet12,vet02
	double precision, intent(inout)::likelihood_res
        integer::np0,i,j,l,w,k,npar0,nva01,nva12,nva02,no0, &
	troncature0,dimnva01,dimnva02,dimnva12,gausspoint0

	double precision,dimension(np0)::b0
        double precision,dimension(npar0)::bh
	double precision,dimension(npar0-np0)::bfix0
	integer,dimension(npar0)::fix0
	double precision,dimension(2)::the01
	double precision,dimension(2)::the12
	double precision,dimension(2)::the02
        double precision,dimension(no0,dimnva01)::ve010
	double precision,dimension(no0,dimnva02)::ve020
	double precision,dimension(no0,dimnva12)::ve120
	
	
        double precision::su01,ri01,su12,ri12,su02,ri02,gl01,gl02,gl12
	double precision,dimension(no0)::t00,t10,t20,t30
	integer,dimension(no0)::c0

	allocate(b(np0),bfix(npar0-np0),fix(npar0))
	b=b0
	bfix=bfix0
	fix=fix0
	
	troncature=troncature0


	if(nva01.gt.0) then 
		allocate(ve01(no0,nva01))
	else 
		allocate(ve01(no0,1))
	end if 
	
	if(nva02.gt.0) then 
		allocate(ve02(no0,nva02))
	else 
		allocate(ve02(no0,1))
	end if 

	if(nva12.gt.0) then 
		allocate(ve12(no0,nva12))
	else 
		allocate(ve12(no0,1))
	end if 


	ve01=ve010
	ve02=ve020
	ve12=ve120

	allocate(t0(no0),t1(no0),t2(no0),t3(no0),c(no0))
	c=c0
	t0=t00
	t1=t10
	t2=t20
	t3=t30

         
        ! we need to put bh at its original values if in posfix 


       l=0
       w=0
       gausspoint=gausspoint0
 

       do k=1,(np0+sum(fix))
         if(fix(k).eq.0) then
            l=l+1
            bh(k)=b(l)
         end if
         if(fix(k).eq.1) then
            w=w+1
            bh(k)=bfix(w)
         end if
      end do
 
	


         do i=1,2
            the01(i)=(bh(i))*(bh(i))
         end do
         do i=1,2
            j = 2+i
            the02(i)=(bh(j))*(bh(j))
         end do
         do i=1,2
            j = 4+i
            the12(i)=(bh(j))*(bh(j))
         end do



	
!---------- calcul de la vraisemblance ------------------



        res = 0.d0
        if(gausspoint.eq.10) then 
        do i=1,no0

                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
                                vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
                                vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
                                vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)


                res1 = 0.d0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else 
                                call fonct(t0(i),the01,ri01,gl01,su01) 
                                call fonct(t0(i),the02,ri02,gl02,su02) 
                                tronc=(gl01*vet01)+(gl02*vet02)
                        end if
                else
                        tronc = 0.d0
                end if

		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2

			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t3(i),the02,ri02,gl02,su02)
                        res1 = -(gl01*vet01)-(gl02*vet02)
                else
                if(c(i).eq.2)then ! cpi 0-->1
			call fonct(t3(i),the12,ri12,gl12,su12)
                         call  qgaussPLweib(t1(i),t2(i),the01,the02,&
                         the12,res2,vet01,vet02,vet12)
                        res1=dLOG(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t1(i),the02,ri02,gl02,su02)
                        call fonct(t1(i),the12,ri12,gl12,su12)
                        res1 = -(gl01*vet01)-(gl02*vet02)+&
                        dLOG(ri01*vet01)+(gl12*vet12)
                        call fonct(t3(i),the12,ri12,gl12,su12)
                        res1 = res1 -(gl12*vet12)
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			call fonct(t3(i),the12,ri12,gl12,su12)
                        call  qgaussPLweib(t1(i),t2(i),the01,the02,the12,&
                        res2,vet01,vet02,vet12)
                        res1=dLOG(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				call fonct(t1(i),the01,ri01,gl01,su01)
                                call fonct(t1(i),the02,ri02,gl02,su02)
                                call fonct(t1(i),the12,ri12,gl12,su12)
                                res1 = -(gl01*vet01)-(gl02*vet02)+&
                                dLOG(ri01*vet01)+(gl12*vet12)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                res1 = res1 -(gl12*vet12) + dLOG(ri12*vet12)
                         else
                            if(c(i).eq.6)then ! vivant ???
				 call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPLweib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12))+&
                                ((su01**vet01)*(su02**vet02))
                                res1 = dLOG(res1)

                            else ! passage 0-->2  

				call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPLweib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12)*ri12*vet12)+&
                                ((su01**vet01)*(su02**vet02)*ri02*vet02)
                                res1 = dLOG(res1)
                            endif
                         endif                        
                      endif
                   endif   
                endif   
                endif   

                res = res + res1 + tronc

                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        likelihood_res=-1.d9
                        goto 123
                end if
        end do  
     else 
       if(gausspoint.eq.15) then 
       
               do i=1,no0

                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
                                vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
                                vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
                                vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)


                res1 = 0.d0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else 
                                call fonct(t0(i),the01,ri01,gl01,su01) 
                                call fonct(t0(i),the02,ri02,gl02,su02) 
                                tronc=(gl01*vet01)+(gl02*vet02)
                        end if
                else
                        tronc = 0.d0
                end if

		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2

			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t3(i),the02,ri02,gl02,su02)
                        res1 = -(gl01*vet01)-(gl02*vet02)
                else
                if(c(i).eq.2)then ! cpi 0-->1
			call fonct(t3(i),the12,ri12,gl12,su12)
                         call  qgaussPL15weib(t1(i),t2(i),the01,the02,&
                         the12,res2,vet01,vet02,vet12)
                        res1=dLOG(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t1(i),the02,ri02,gl02,su02)
                        call fonct(t1(i),the12,ri12,gl12,su12)
                        res1 = -(gl01*vet01)-(gl02*vet02)+&
                        dLOG(ri01*vet01)+(gl12*vet12)
                        call fonct(t3(i),the12,ri12,gl12,su12)
                        res1 = res1 -(gl12*vet12)
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			call fonct(t3(i),the12,ri12,gl12,su12)
                        call  qgaussPL15weib(t1(i),t2(i),the01,the02,the12,&
                        res2,vet01,vet02,vet12)
                        res1=dLOG(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				call fonct(t1(i),the01,ri01,gl01,su01)
                                call fonct(t1(i),the02,ri02,gl02,su02)
                                call fonct(t1(i),the12,ri12,gl12,su12)
                                res1 = -(gl01*vet01)-(gl02*vet02)+&
                                dLOG(ri01*vet01)+(gl12*vet12)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                res1 = res1 -(gl12*vet12) + dLOG(ri12*vet12)
                         else
                            if(c(i).eq.6)then ! vivant ???
				 call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL15weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12))+&
                                ((su01**vet01)*(su02**vet02))
                                res1 = dLOG(res1)

                            else ! passage 0-->2  

				call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL15weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12)*ri12*vet12)+&
                                ((su01**vet01)*(su02**vet02)*ri02*vet02)
                                res1 = dLOG(res1)
                            endif
                         endif                        
                      endif
                   endif   
                endif   
                endif   

                res = res + res1 + tronc

                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        likelihood_res=-1.d9
                        goto 123
                end if
        end do   
 
 else 
       if(gausspoint.eq.21) then 
       
             do i=1,no0

                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
                                vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
                                vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
                                vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)


                res1 = 0.d0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else 
                                call fonct(t0(i),the01,ri01,gl01,su01) 
                                call fonct(t0(i),the02,ri02,gl02,su02) 
                                tronc=(gl01*vet01)+(gl02*vet02)
                        end if
                else
                        tronc = 0.d0
                end if

		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2

			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t3(i),the02,ri02,gl02,su02)
                        res1 = -(gl01*vet01)-(gl02*vet02)
                else
                if(c(i).eq.2)then ! cpi 0-->1
			call fonct(t3(i),the12,ri12,gl12,su12)
                         call  qgaussPL21weib(t1(i),t2(i),the01,the02,&
                         the12,res2,vet01,vet02,vet12)
                        res1=dLOG(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t1(i),the02,ri02,gl02,su02)
                        call fonct(t1(i),the12,ri12,gl12,su12)
                        res1 = -(gl01*vet01)-(gl02*vet02)+&
                        dLOG(ri01*vet01)+(gl12*vet12)
                        call fonct(t3(i),the12,ri12,gl12,su12)
                        res1 = res1 -(gl12*vet12)
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			call fonct(t3(i),the12,ri12,gl12,su12)
                        call  qgaussPL21weib(t1(i),t2(i),the01,the02,the12,&
                        res2,vet01,vet02,vet12)
                        res1=dLOG(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				call fonct(t1(i),the01,ri01,gl01,su01)
                                call fonct(t1(i),the02,ri02,gl02,su02)
                                call fonct(t1(i),the12,ri12,gl12,su12)
                                res1 = -(gl01*vet01)-(gl02*vet02)+&
                                dLOG(ri01*vet01)+(gl12*vet12)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                res1 = res1 -(gl12*vet12) + dLOG(ri12*vet12)
                         else
                            if(c(i).eq.6)then ! vivant ???
				 call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL21weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12))+&
                                ((su01**vet01)*(su02**vet02))
                                res1 = dLOG(res1)

                            else ! passage 0-->2  

				call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL21weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12)*ri12*vet12)+&
                                ((su01**vet01)*(su02**vet02)*ri02*vet02)
                                res1 = dLOG(res1)
                            endif
                         endif                        
                      endif
                   endif   
                endif   
                endif   

                res = res + res1 + tronc

                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        likelihood_res=-1.d9
                        goto 123
                end if
        end do 
        
    else 
          if(gausspoint.eq.31) then 
          
                  do i=1,no0

                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
                                vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
                                vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
                                vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)


                res1 = 0.d0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else 
                                call fonct(t0(i),the01,ri01,gl01,su01) 
                                call fonct(t0(i),the02,ri02,gl02,su02) 
                                tronc=(gl01*vet01)+(gl02*vet02)
                        end if
                else
                        tronc = 0.d0
                end if

		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2

			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t3(i),the02,ri02,gl02,su02)
                        res1 = -(gl01*vet01)-(gl02*vet02)
                else
                if(c(i).eq.2)then ! cpi 0-->1
			call fonct(t3(i),the12,ri12,gl12,su12)
                         call  qgaussPL31weib(t1(i),t2(i),the01,the02,&
                         the12,res2,vet01,vet02,vet12)
                        res1=dLOG(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t1(i),the02,ri02,gl02,su02)
                        call fonct(t1(i),the12,ri12,gl12,su12)
                        res1 = -(gl01*vet01)-(gl02*vet02)+&
                        dLOG(ri01*vet01)+(gl12*vet12)
                        call fonct(t3(i),the12,ri12,gl12,su12)
                        res1 = res1 -(gl12*vet12)
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			call fonct(t3(i),the12,ri12,gl12,su12)
                        call  qgaussPL31weib(t1(i),t2(i),the01,the02,the12,&
                        res2,vet01,vet02,vet12)
                        res1=dLOG(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				call fonct(t1(i),the01,ri01,gl01,su01)
                                call fonct(t1(i),the02,ri02,gl02,su02)
                                call fonct(t1(i),the12,ri12,gl12,su12)
                                res1 = -(gl01*vet01)-(gl02*vet02)+&
                                dLOG(ri01*vet01)+(gl12*vet12)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                res1 = res1 -(gl12*vet12) + dLOG(ri12*vet12)
                         else
                            if(c(i).eq.6)then ! vivant ???
				 call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL31weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12))+&
                                ((su01**vet01)*(su02**vet02))
                                res1 = dLOG(res1)

                            else ! passage 0-->2  

				call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL31weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12)*ri12*vet12)+&
                                ((su01**vet01)*(su02**vet02)*ri02*vet02)
                                res1 = dLOG(res1)
                            endif
                         endif                        
                      endif
                   endif   
                endif   
                endif   

                res = res + res1 + tronc

                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        likelihood_res=-1.d9
                        goto 123
                end if
        end do 
        
   else 
         if(gausspoint.eq.41) then 
         
                 do i=1,no0

                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
                                vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
                                vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
                                vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)


                res1 = 0.d0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else 
                                call fonct(t0(i),the01,ri01,gl01,su01) 
                                call fonct(t0(i),the02,ri02,gl02,su02) 
                                tronc=(gl01*vet01)+(gl02*vet02)
                        end if
                else
                        tronc = 0.d0
                end if

		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2

			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t3(i),the02,ri02,gl02,su02)
                        res1 = -(gl01*vet01)-(gl02*vet02)
                else
                if(c(i).eq.2)then ! cpi 0-->1
			call fonct(t3(i),the12,ri12,gl12,su12)
                         call  qgaussPL41weib(t1(i),t2(i),the01,the02,&
                         the12,res2,vet01,vet02,vet12)
                        res1=dLOG(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t1(i),the02,ri02,gl02,su02)
                        call fonct(t1(i),the12,ri12,gl12,su12)
                        res1 = -(gl01*vet01)-(gl02*vet02)+&
                        dLOG(ri01*vet01)+(gl12*vet12)
                        call fonct(t3(i),the12,ri12,gl12,su12)
                        res1 = res1 -(gl12*vet12)
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			call fonct(t3(i),the12,ri12,gl12,su12)
                        call  qgaussPL41weib(t1(i),t2(i),the01,the02,the12,&
                        res2,vet01,vet02,vet12)
                        res1=dLOG(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				call fonct(t1(i),the01,ri01,gl01,su01)
                                call fonct(t1(i),the02,ri02,gl02,su02)
                                call fonct(t1(i),the12,ri12,gl12,su12)
                                res1 = -(gl01*vet01)-(gl02*vet02)+&
                                dLOG(ri01*vet01)+(gl12*vet12)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                res1 = res1 -(gl12*vet12) + dLOG(ri12*vet12)
                         else
                            if(c(i).eq.6)then ! vivant ???
				 call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL41weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12))+&
                                ((su01**vet01)*(su02**vet02))
                                res1 = dLOG(res1)

                            else ! passage 0-->2  

				call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL41weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12)*ri12*vet12)+&
                                ((su01**vet01)*(su02**vet02)*ri02*vet02)
                                res1 = dLOG(res1)
                            endif
                         endif                        
                      endif
                   endif   
                endif   
                endif   

                res = res + res1 + tronc

                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        likelihood_res=-1.d9
                        goto 123
                end if
        end do   
  
else 
      if(gausspoint.eq.51) then 
            
            do i=1,no0

                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
                                vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
                                vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
                                vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)


                res1 = 0.d0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else 
                                call fonct(t0(i),the01,ri01,gl01,su01) 
                                call fonct(t0(i),the02,ri02,gl02,su02) 
                                tronc=(gl01*vet01)+(gl02*vet02)
                        end if
                else
                        tronc = 0.d0
                end if

		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2

			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t3(i),the02,ri02,gl02,su02)
                        res1 = -(gl01*vet01)-(gl02*vet02)
                else
                if(c(i).eq.2)then ! cpi 0-->1
			call fonct(t3(i),the12,ri12,gl12,su12)
                         call  qgaussPL51weib(t1(i),t2(i),the01,the02,&
                         the12,res2,vet01,vet02,vet12)
                        res1=dLOG(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t1(i),the02,ri02,gl02,su02)
                        call fonct(t1(i),the12,ri12,gl12,su12)
                        res1 = -(gl01*vet01)-(gl02*vet02)+&
                        dLOG(ri01*vet01)+(gl12*vet12)
                        call fonct(t3(i),the12,ri12,gl12,su12)
                        res1 = res1 -(gl12*vet12)
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			call fonct(t3(i),the12,ri12,gl12,su12)
                        call  qgaussPL51weib(t1(i),t2(i),the01,the02,the12,&
                        res2,vet01,vet02,vet12)
                        res1=dLOG(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				call fonct(t1(i),the01,ri01,gl01,su01)
                                call fonct(t1(i),the02,ri02,gl02,su02)
                                call fonct(t1(i),the12,ri12,gl12,su12)
                                res1 = -(gl01*vet01)-(gl02*vet02)+&
                                dLOG(ri01*vet01)+(gl12*vet12)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                res1 = res1 -(gl12*vet12) + dLOG(ri12*vet12)
                         else
                            if(c(i).eq.6)then ! vivant ???
				 call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL51weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12))+&
                                ((su01**vet01)*(su02**vet02))
                                res1 = dLOG(res1)

                            else ! passage 0-->2  

				call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL51weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12)*ri12*vet12)+&
                                ((su01**vet01)*(su02**vet02)*ri02*vet02)
                                res1 = dLOG(res1)
                            endif
                         endif                        
                      endif
                   endif   
                endif   
                endif   

                res = res + res1 + tronc

                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        likelihood_res=-1.d9
                        goto 123
                end if
        end do   
        
  else 
       
       do i=1,no0

                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
                                vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
                                vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
                                vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)


                res1 = 0.d0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else 
                                call fonct(t0(i),the01,ri01,gl01,su01) 
                                call fonct(t0(i),the02,ri02,gl02,su02) 
                                tronc=(gl01*vet01)+(gl02*vet02)
                        end if
                else
                        tronc = 0.d0
                end if

		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2

			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t3(i),the02,ri02,gl02,su02)
                        res1 = -(gl01*vet01)-(gl02*vet02)
                else
                if(c(i).eq.2)then ! cpi 0-->1
			call fonct(t3(i),the12,ri12,gl12,su12)
                         call  qgaussPL61weib(t1(i),t2(i),the01,the02,&
                         the12,res2,vet01,vet02,vet12)
                        res1=dLOG(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t1(i),the02,ri02,gl02,su02)
                        call fonct(t1(i),the12,ri12,gl12,su12)
                        res1 = -(gl01*vet01)-(gl02*vet02)+&
                        dLOG(ri01*vet01)+(gl12*vet12)
                        call fonct(t3(i),the12,ri12,gl12,su12)
                        res1 = res1 -(gl12*vet12)
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			call fonct(t3(i),the12,ri12,gl12,su12)
                        call  qgaussPL61weib(t1(i),t2(i),the01,the02,the12,&
                        res2,vet01,vet02,vet12)
                        res1=dLOG(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				call fonct(t1(i),the01,ri01,gl01,su01)
                                call fonct(t1(i),the02,ri02,gl02,su02)
                                call fonct(t1(i),the12,ri12,gl12,su12)
                                res1 = -(gl01*vet01)-(gl02*vet02)+&
                                dLOG(ri01*vet01)+(gl12*vet12)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                res1 = res1 -(gl12*vet12) + dLOG(ri12*vet12)
                         else
                            if(c(i).eq.6)then ! vivant ???
				 call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL61weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12))+&
                                ((su01**vet01)*(su02**vet02))
                                res1 = dLOG(res1)

                            else ! passage 0-->2  

				call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL61weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12)*ri12*vet12)+&
                                ((su01**vet01)*(su02**vet02)*ri02*vet02)
                                res1 = dLOG(res1)
                            endif
                         endif                        
                      endif
                   endif   
                endif   
                endif   

                res = res + res1 + tronc

                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        likelihood_res=-1.d9
                        goto 123
                end if
        end do   


        end if 
        end if 
        end if
        end if
        end if 
        end if 
  


        likelihood_res = res


123     continue 
	 
	deallocate(b,bfix,fix,ve01,ve02,ve12,t0,t1,t2,t3,c)

        end subroutine idmlikelihoodweib


!============================================================================================= 
!==========================  SUSP  ===========================================================
!===================== calculate splines distribution ========================================
!============================================================================================= 

        subroutine susp(x,the,n,su,lam,zi,gl)

        implicit none
        
        integer::j,k,n,i
        double precision::x,ht,ht2,h2,som,su,lam,htm,h2t,h3,h2n,hn, &
        im,im1,im2,mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2,h,gl,hh
        double precision,dimension(-2:(n+3))::zi
        double precision,dimension(-2:n-1)::the 

        som = 0.d0
        gl = 0.d0 
        
        do k = 2,n
                if ((x.ge.zi(k-1)).and.(x.lt.zi(k)))then
                        j = k-1
                        if (j.gt.1)then
                                do i=2,j
                                som = som+the(i-4)
                                end do  
                        endif   
                        ht = x-zi(j)
                        htm= x-zi(j-1)
                        h2t= x-zi(j+2)
                        ht2 = zi(j+1)-x
                        ht3 = zi(j+3)-x
                        hht = x-zi(j-2)
                        h = zi(j+1)-zi(j)
                        hh= zi(j+1)-zi(j-1)
                        h2= zi(j+2)-zi(j)
                        h3= zi(j+3)-zi(j)
                        h4= zi(j+4)-zi(j)
                        h3m= zi(j+3)-zi(j-1)
                        h2n=zi(j+2)-zi(j-1)
                        hn= zi(j+1)-zi(j-2)
                        hh3 = zi(j+1)-zi(j-3)
                        hh2 = zi(j+2)-zi(j-2)
                        mm3 = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
                        mm2 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
                        *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
                        mm1 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
                        h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
                        mm  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
                        im3 = (0.25d0*(x-zi(j-3))*mm3)+(0.25d0*hh2*mm2) &
                        +(0.25d0*h3m*mm1)+(0.25d0*h4*mm)
                        im2 = (0.25d0*hht*mm2)+(h3m*mm1*0.25d0)+(h4*mm*0.25d0)
                        im1 = (htm*mm1*0.25d0)+(h4*mm*0.25d0)
                        im  = ht*mm*0.25d0
                        gl = som +(the(j-3)*im3)+(the(j-2)*im2)+(the(j-1)*im1)+(the(j)*im)
                        lam = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1)+(the(j)*mm)
                endif
        end do
   
        if(x.ge.zi(n))then
                som = 0.d0
                do i=1,n+2
                        som = som+the(i-3)
                end do
                gl = som
                lam = 4.d0*the(n-1)/(zi(n)-zi(n-1))
        endif

        su  = dexp(-gl)

        return

        end subroutine susp


!============================================================================================= 
!====================================== FONCT  ===============================================
!======================= calculate weibull distribution  =====================================
!============================================================================================= 

subroutine fonct(x,p,risq,glam,surv)

        implicit none

        double precision,dimension(2)::p
        double precision::x,surv,risq,glam


	
        	surv = dexp(-(p(2)*x)**p(1)) !exponential of - integral over 0 to x of base risk
        	glam = (p(2)*x)**p(1) !integral over 0 to x of base risk
        	risq = p(1)*(p(2)**p(1))*(x**(p(1)-1.d0)) ! base risk 

        if (x.le.0.d0) then
                surv = 1.d0
                glam = 0.d0
                risq = 0.d0
        endif

        return

end subroutine fonct
!============================================================================================= 
!================================   QGAUS for weib : CHEBYCHEV   =============================
!============================================================================================= 

subroutine qgaussPLweib(a,b,the01,the02,the12,res,v01,v02,v12)
        implicit none
         double precision a,b,the01(2),the02(2),the12(2)
         double precision dx,xm,xr,w(5),x(5),res,v01,v02,v12
         double precision xx,f1,su01,ri01,ri12,f2,su12,su02,ri02
         double precision gl01,gl12,gl02
         integer j
         save w,x
         data w/0.2955242247d0,0.2692667193d0,0.2190863625d0,0.1494513491d0,0.0666713443d0/
         data x/0.1488743389d0,0.4333953941d0,0.6794095682d0,0.8650633666d0,0.9739065285d0/

         
            xm = 0.5d0*(b+a)
            xr = 0.5d0*(b-a)
            res = 0.d0
            if(a.eq.b)then
               res = 0.d0
            else
               do 11 j=1,5
                  dx=xr*x(j)
                  xx = xm+dx
                  call fonct(xx,the01,ri01,gl01,su01)
                  call fonct(xx,the02,ri02,gl02,su02)
                  call fonct(xx,the12,ri12,gl12,su12)
                  f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
                  xx = xm-dx
                  call fonct(xx,the01,ri01,gl01,su01)
                  call fonct(xx,the02,ri02,gl02,su02)
                  call fonct(xx,the12,ri12,gl12,su12)
                  f2 = ((su01**v01)*(su02**v02)*ri01*v01)/(su12**v12)
                  res = res + w(j)*(f1+f2)
 11            continue
            endif
            res = res*xr

            
          end subroutine qgaussPLweib

!============================================================================================= 
!================================  QGAUS for splines : CHEBYCHEV   ===========================
!============================================================================================= 

      subroutine qgaussPL(a,b,the01,the12,the02,res,v1,v2,v3)

  	use commun,only:zi01,zi12,zi02,nz01,nz12,nz02

         implicit none
         
         integer::j
         double precision::a,b,dx,xm,xr,res,v1,v2,v3
         double precision,dimension(-2:(nz01-1))::the01
        double precision,dimension(-2:(nz12-1))::the12
        double precision,dimension(-2:(nz02-1))::the02
         double precision,dimension(5)::w,x
         double precision::xx,f1,su01,ri01,ri12,f2,su12,su02,ri02,gl01,gl02,gl12
         save w,x
         data w/0.2955242247d0,0.2692667193d0,0.2190863625d0, &
               0.1494513491d0,0.0666713443d0/
         data x/0.1488743389d0,0.4333953941d0,0.6794095682d0, &
               0.8650633666d0,0.9739065285d0/
    
            xm = 0.5d0*(b+a)
            xr = 0.5d0*(b-a)
            res = 0.d0
            do j=1,5
               dx=xr*x(j)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f1 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f2 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               res = res + w(j)*(f1+f2)
            end do
            res = res*xr

         end subroutine qgaussPL


!=============================================================================================  
!==== QGAUS15 out a 15 point Gauss-Kronrod quadrature rule for weib  =========================
!=============================================================================================  
subroutine qgaussPL15weib(a,b,the01,the02,the12,res,v01,v02,v12)
         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resk,v01,v02,v12,&
         fv1,fv2,d1mach(5),epmach,uflow,the01(2),the12(2),the02(2)
         double precision,dimension(8)::xgk,wgk
	 double precision,dimension(4)::wg
         double precision::xx,f1,su01,ri01,ri12,f2,su12,su02,ri02,fc,gl01,gl02,gl12
         save wgk,xgk
	 dimension fv1(7),fv2(7)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)

	wg(1)=0.129484966168869693270611432679082d0
   	wg(2)=0.279705391489276667901467771423780d0
    	wg(3)=0.381830050505118944950369775488975d0
    	wg(4)=0.417959183673469387755102040816327d0

    	xgk(1)=0.991455371120812639206854697526329d0
    	xgk(2)=0.949107912342758524526189684047851d0
    	xgk(3)=0.864864423359769072789712788640926d0
    	xgk(4)=0.741531185599394439863864773280788d0
    	xgk(5)=0.586087235467691130294144838258730d0
    	xgk(6)=0.405845151377397166906606412076961d0
    	xgk(7)=0.207784955007898467600689403773245d0
    	xgk(8)=0.000000000000000000000000000000000d0

    	wgk(1)=0.022935322010529224963732008058970d0
    	wgk(2)=0.063092092629978553290700663189204d0
    	wgk(3)=0.104790010322250183839876322541518d0
    	wgk(4)=0.140653259715525918745189590510238d0
    	wgk(5)=0.169004726639267902826583426598550d0
    	wgk(6)=0.190350578064785409913256402421014d0
    	wgk(7)=0.204432940075298892414161999234649d0
    	wgk(8)=0.209482141084727828012999174891714d0
     

        xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)
        call fonct(xm,the01,ri01,gl01,su01)
        call fonct(xm,the02,ri02,gl02,su02)
        call fonct(xm,the12,ri12,gl12,su12)
        fc = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0

    	
        resk = fc*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
         
            if(a.eq.b)then
               res = 0.d0
            else
               do j=1,3
               	dx=xr*xgk(jtw)
               	xx = xm+dx
               	call fonct(xx,the01,ri01,gl01,su01)
               	call fonct(xx,the02,ri02,gl02,su02)
               	call fonct(xx,the12,ri12,gl12,su12)
               	f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               	xx = xm-dx
               	call fonct(xx,the01,ri01,gl01,su01)
               	call fonct(xx,the02,ri02,gl02,su02)
               	call fonct(xx,the12,ri12,gl12,su12)
               	f2 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               	fv1(jtw) = f1   ! svgrd valeurs fct f a gche du centre
               	fv2(jtw) = f2   ! svgrd valeurs fct f a drte du centre
	       	
               	resk = resk + wgk(jtw)*(f1+f2)

              end do
	      do j=1,4
	       jtwm1 = j*2-1
               dx=xr*xgk(jtwm1)
               xx = xm+dx
               call fonct(xx,the01,ri01,gl01,su01)
               call fonct(xx,the02,ri02,gl02,su02)
               call fonct(xx,the12,ri12,gl12,su12)
               f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               xx = xm-dx
               call fonct(xx,the01,ri01,gl01,su01)
               call fonct(xx,the02,ri02,gl02,su02)
               call fonct(xx,the12,ri12,gl12,su12)
               f2 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               fv1(jtwm1) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtwm1) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtwm1)*(f1+f2)
              end do
	    
    	res = xr*resk
	endif
    
          end subroutine qgaussPL15weib

!=============================================================================================  
!==== QGAUS15 out a 15 point Gauss-Kronrod quadrature rule for splines   =====================
!=============================================================================================  

      subroutine qgaussPL15(a,b,the01,the12,the02,res,v1,v2,v3)

 
	use commun,only:zi01,zi12,zi02,nz01,nz12,nz02

         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resk,v1,v2,v3,&
         fv1,fv2,d1mach(5),epmach,uflow
         double precision,dimension(-2:(nz01-1))::the01
         double precision,dimension(-2:(nz12-1))::the12
         double precision,dimension(-2:(nz02-1))::the02
         double precision,dimension(8)::xgk,wgk
	 double precision,dimension(4)::wg
         double precision::xx,f1,su01,ri01,ri12,f2,su12,su02,ri02,fc,gl01,gl02,gl12
         save wgk,xgk
	 dimension fv1(7),fv2(7)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)

	wg(1)=0.129484966168869693270611432679082d0
   	wg(2)=0.279705391489276667901467771423780d0
    	wg(3)=0.381830050505118944950369775488975d0
    	wg(4)=0.417959183673469387755102040816327d0

    	xgk(1)=0.991455371120812639206854697526329d0
    	xgk(2)=0.949107912342758524526189684047851d0
    	xgk(3)=0.864864423359769072789712788640926d0
    	xgk(4)=0.741531185599394439863864773280788d0
    	xgk(5)=0.586087235467691130294144838258730d0
    	xgk(6)=0.405845151377397166906606412076961d0
    	xgk(7)=0.207784955007898467600689403773245d0
    	xgk(8)=0.000000000000000000000000000000000d0

    	wgk(1)=0.022935322010529224963732008058970d0
    	wgk(2)=0.063092092629978553290700663189204d0
    	wgk(3)=0.104790010322250183839876322541518d0
    	wgk(4)=0.140653259715525918745189590510238d0
    	wgk(5)=0.169004726639267902826583426598550d0
    	wgk(6)=0.190350578064785409913256402421014d0
    	wgk(7)=0.204432940075298892414161999234649d0
    	wgk(8)=0.209482141084727828012999174891714d0

        xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)
        call susp(xm,the01,nz01,su01,ri01,zi01,gl01)
        call susp(xm,the02,nz02,su02,ri02,zi02,gl02)
        call susp(xm,the12,nz12,su12,ri12,zi12,gl12)
        fc = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0

        resk = fc*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        

         do j=1,3
	       jtw = j*2
               dx=xr*xgk(jtw)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f1 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f2 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               fv1(jtw) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtw) = f2   ! svgrd valeurs fct f a drte du centre
               resk = resk + wgk(jtw)*(f1+f2)
         end do
	 do j=1,4
	       jtwm1 = j*2-1
               dx=xr*xgk(jtwm1)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f1 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f2 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               fv1(jtwm1) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtwm1) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtwm1)*(f1+f2)
         end do
    
    res = xr*resk
  
         end subroutine qgaussPL15

!=============================================================================================  
!=====QGAUS21 out a 21 point Gauss-Kronrod quadrature rule for weib ==========================
!=============================================================================================  

subroutine qgaussPL21weib(a,b,the01,the02,the12,res,v01,v02,v12)
         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resk,v01,v02,v12,&
         fv1,fv2,d1mach(5),epmach,uflow,the01(2),the12(2),the02(2)
         double precision,dimension(11)::xgk,wgk
	 double precision,dimension(5)::wg
         double precision::xx,f1,su01,ri01,ri12,f2,su12,su02,ri02,fc,gl01,gl02,gl12
         save wgk,xgk
	 dimension fv1(10),fv2(10)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)

	wg(1)=0.066713443086881375935688098933317928579d0
   	wg(2)=0.149451349150580593145776339657697332403d0
    	wg(3)=0.219086362515982043995534934228163192459d0
    	wg(4)=0.269266719309996355091226921569469352860d0
        wg(5)=0.295524224714752870173892994651338329421d0

    	xgk(1)=0.995657163025808080735527280689002847921d0
    	xgk(2)=0.973906528517171720077964012084452053428d0
    	xgk(3)=0.930157491355708226001207180059508346225d0
    	xgk(4)=0.865063366688984510732096688423493048528d0
    	xgk(5)=0.780817726586416897063717578345042377163d0
    	xgk(6)=0.679409568299024406234327365114873575769d0
    	xgk(7)=0.562757134668604683339000099272694140843d0
        xgk(8)=0.433395394129247190799265943165784162200d0
        xgk(9)=0.294392862701460198131126603103865566163d0
        xgk(10)=0.148874338981631210884826001129719984618d0
    	xgk(11)=0.000000000000000000000000000000000d0
	
    	wgk(1)=0.0116946388673718742780643960621920483962d0
    	wgk(2)=0.0325581623079647274788189724593897606174d0
    	wgk(3)=0.0547558965743519960313813002445801763737d0
    	wgk(4)=0.0750396748109199527670431409161900093952d0
    	wgk(5)=0.0931254545836976055350654650833663443900d0
    	wgk(6)=0.109387158802297641899210590325804960272d0
    	wgk(7)=0.123491976262065851077958109831074159512d0
    	wgk(8)=0.134709217311473325928054001771706832761d0
	wgk(9)=0.142775938577060080797094273138717060886d0
	wgk(10)=0.147739104901338491374841515972068045524d0
	wgk(11)=0.149445554002916905664936468389821203745d0

        xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)
        call fonct(xm,the01,ri01,gl01,su01)
        call fonct(xm,the02,ri02,gl02,su02)
        call fonct(xm,the12,ri12,gl12,su12)
        fc = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0

    	resk = fc*wgk(11)       ! init res Kronrod   ! fc * 8e poids Kronrod
        
            if(a.eq.b)then
               res = 0.d0
            else
               do j=1,5
               	jtw = j*2
               	dx=xr*xgk(jtw)
               	xx = xm+dx
               	call fonct(xx,the01,ri01,gl01,su01)
               	call fonct(xx,the02,ri02,gl02,su02)
               	call fonct(xx,the12,ri12,gl12,su12)
               	f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               	xx = xm-dx
               	call fonct(xx,the01,ri01,gl01,su01)
               	call fonct(xx,the02,ri02,gl02,su02)
               	call fonct(xx,the12,ri12,gl12,su12)
               	f2 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               	fv1(jtw) = f1   ! svgrd valeurs fct f a gche du centre
               	fv2(jtw) = f2   ! svgrd valeurs fct f a drte du centre
	       	resk = resk + wgk(jtw)*(f1+f2)
                
              end do
	      do j=1,5
	       jtwm1 = j*2-1
               dx=xr*xgk(jtwm1)
               xx = xm+dx
               call fonct(xx,the01,ri01,gl01,su01)
               call fonct(xx,the02,ri02,gl02,su02)
               call fonct(xx,the12,ri12,gl12,su12)
               f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               xx = xm-dx
               call fonct(xx,the01,ri01,gl01,su01)
               call fonct(xx,the02,ri02,gl02,su02)
               call fonct(xx,the12,ri12,gl12,su12)
               f2 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               fv1(jtwm1) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtwm1) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtwm1)*(f1+f2)
              end do
	    

    	res = xr*resk
	endif
    
          end subroutine qgaussPL21weib

!=============================================================================================  
!=====QGAUS21 out a 21 point Gauss-Kronrod quadrature rule for splines =======================
!=============================================================================================  


      subroutine qgaussPL21(a,b,the01,the12,the02,res,v1,v2,v3)

	use commun,only:zi01,zi12,zi02,nz01,nz12,nz02

         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resk,v1,v2,v3,&
         fv1,fv2,d1mach(5),epmach,uflow
         double precision,dimension(-2:(nz01-1))::the01
         double precision,dimension(-2:(nz12-1))::the12
         double precision,dimension(-2:(nz02-1))::the02
         double precision,dimension(11)::xgk,wgk
	 double precision,dimension(5)::wg
         double precision::xx,f1,su01,ri01,ri12,f2,su12,su02,ri02,fc,gl01,gl02,gl12
         save wgk,xgk
	 dimension fv1(10),fv2(10)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)

	wg(1)=0.066713443086881375935688098933317928579d0
   	wg(2)=0.149451349150580593145776339657697332403d0
    	wg(3)=0.219086362515982043995534934228163192459d0
    	wg(4)=0.269266719309996355091226921569469352860d0
        wg(5)=0.295524224714752870173892994651338329421d0

    	xgk(1)=0.995657163025808080735527280689002847921d0
    	xgk(2)=0.973906528517171720077964012084452053428d0
    	xgk(3)=0.930157491355708226001207180059508346225d0
    	xgk(4)=0.865063366688984510732096688423493048528d0
    	xgk(5)=0.780817726586416897063717578345042377163d0
    	xgk(6)=0.679409568299024406234327365114873575769d0
    	xgk(7)=0.562757134668604683339000099272694140843d0
        xgk(8)=0.433395394129247190799265943165784162200d0
        xgk(9)=0.294392862701460198131126603103865566163d0
        xgk(10)=0.148874338981631210884826001129719984618d0
    	xgk(11)=0.000000000000000000000000000000000d0
	
    	wgk(1)=0.0116946388673718742780643960621920483962d0
    	wgk(2)=0.0325581623079647274788189724593897606174d0
    	wgk(3)=0.0547558965743519960313813002445801763737d0
    	wgk(4)=0.0750396748109199527670431409161900093952d0
    	wgk(5)=0.0931254545836976055350654650833663443900d0
    	wgk(6)=0.109387158802297641899210590325804960272d0
    	wgk(7)=0.123491976262065851077958109831074159512d0
    	wgk(8)=0.134709217311473325928054001771706832761d0
	wgk(9)=0.142775938577060080797094273138717060886d0
	wgk(10)=0.147739104901338491374841515972068045524d0
	wgk(11)=0.149445554002916905664936468389821203745d0

        xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)
        call susp(xm,the01,nz01,su01,ri01,zi01,gl01)
        call susp(xm,the02,nz02,su02,ri02,zi02,gl02)
        call susp(xm,the12,nz12,su12,ri12,zi12,gl12)
        fc = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0

    	resk = fc*wgk(11)       ! init res Kronrod   ! fc * 8e poids Kronrod
         

         do j=1,5
	       jtw = j*2
               dx=xr*xgk(jtw)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f1 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f2 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               fv1(jtw) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtw) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtw)*(f1+f2)
               
         end do
	 do j=1,5
	       jtwm1 = j*2-1
               dx=xr*xgk(jtwm1)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f1 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f2 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               fv1(jtwm1) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtwm1) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtwm1)*(f1+f2)
               
         end do

    res = xr*resk
    
         end subroutine qgaussPL21

!=============================================================================================  
!=====QGAUS31 out a 31 point Gauss-Kronrod quadrature rule for weib ==========================
!=============================================================================================  


subroutine qgaussPL31weib(a,b,the01,the02,the12,res,v01,v02,v12)
         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resk,v01,v02,v12,&
         fv1,fv2,d1mach(5),epmach,uflow,the01(2),the12(2),the02(2)
         double precision,dimension(16)::xgk,wgk
	 double precision,dimension(8)::wg
         double precision::xx,f1,su01,ri01,ri12,f2,su12,su02,ri02,fc,gl01,gl02,gl12
         save wgk,xgk
	 dimension fv1(15),fv2(15)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)
                                            
	wg(1)=0.0307532419961172683546283935772044177217d0                                   
   	wg(2)=0.0703660474881081247092674164506673384667d0                                   
    	wg(3)=0.107159220467171935011869546685869303416d0                                     
    	wg(4)=0.139570677926154314447804794511028322521d0                                    
        wg(5)=0.166269205816993933553200860481208811131d0
	wg(6)=0.186161000015562211026800561866422824506d0
        wg(7)=0.198431485327111576456118326443839324819d0
        wg(8)=0.202578241925561272880620199967519314839d0

                                             
    	xgk(1)=0.998002298693397060285172840152271209073d0
    	xgk(2)=0.987992518020485428489565718586612581147d0
    	xgk(3)=0.967739075679139134257347978784337225283d0
    	xgk(4)=0.937273392400705904307758947710209471244d0
    	xgk(5)=0.897264532344081900882509656454495882832d0
    	xgk(6)=0.848206583410427216200648320774216851366d0
    	xgk(7)=0.790418501442465932967649294817947346862d0
        xgk(8)=0.724417731360170047416186054613938009631d0
        xgk(9)=0.650996741297416970533735895313274692547d0
        xgk(10)=0.570972172608538847537226737253910641238d0
	xgk(11)=0.485081863640239680693655740232350612866d0
	xgk(12)=0.394151347077563369897207370981045468363d0
	xgk(13)=0.299180007153168812166780024266388962662d0
	xgk(14)=0.201194093997434522300628303394596207813d0
	xgk(15)=0.101142066918717499027074231447392338787d0
    	xgk(16)=0.000000000000000000000000000000000d0
	
    	wgk(1)=0.00537747987292334898779205143012764981831d0
    	wgk(2)=0.0150079473293161225383747630758072680946d0
    	wgk(3)=0.0254608473267153201868740010196533593973d0
    	wgk(4)=0.0353463607913758462220379484783600481226d0
    	wgk(5)=0.0445897513247648766082272993732796902233d0
    	wgk(6)=0.0534815246909280872653431472394302967716d0
    	wgk(7)=0.0620095678006706402851392309608029321904d0
    	wgk(8)=0.0698541213187282587095200770991474757860d0
	wgk(9)=0.0768496807577203788944327774826590067221d0
	wgk(10)=0.0830805028231330210382892472861037896016d0
	wgk(11)=0.0885644430562117706472754436937743032123d0
	wgk(12)=0.0931265981708253212254868727473457185619d0
	wgk(13)=0.0966427269836236785051799076275893351367d0
	wgk(14)=0.0991735987217919593323931734846031310596d0
	wgk(15)=0.100769845523875595044946662617569721916d0
	wgk(16)=0.101330007014791549017374792767492546771d0

        xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)
        call fonct(xm,the01,ri01,gl01,su01)
        call fonct(xm,the02,ri02,gl02,su02)
        call fonct(xm,the12,ri12,gl12,su12)
        fc = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0

    	resk = fc*wgk(16)       ! init res Kronrod   ! fc * 8e poids Kronrod
         
            if(a.eq.b)then
               res = 0.d0
            else
               do j=1,7
               	jtw = j*2
               	dx=xr*xgk(jtw)
               	xx = xm+dx
               	call fonct(xx,the01,ri01,gl01,su01)
               	call fonct(xx,the02,ri02,gl02,su02)
               	call fonct(xx,the12,ri12,gl12,su12)
               	f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               	xx = xm-dx
               	call fonct(xx,the01,ri01,gl01,su01)
               	call fonct(xx,the02,ri02,gl02,su02)
               	call fonct(xx,the12,ri12,gl12,su12)
               	f2 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               	fv1(jtw) = f1   ! svgrd valeurs fct f a gche du centre
               	fv2(jtw) = f2   ! svgrd valeurs fct f a drte du centre
	       	resk = resk + wgk(jtw)*(f1+f2)
              end do
	      do j=1,8
	       jtwm1 = j*2-1
               dx=xr*xgk(jtwm1)
               xx = xm+dx
               call fonct(xx,the01,ri01,gl01,su01)
               call fonct(xx,the02,ri02,gl02,su02)
               call fonct(xx,the12,ri12,gl12,su12)
               f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               xx = xm-dx
               call fonct(xx,the01,ri01,gl01,su01)
               call fonct(xx,the02,ri02,gl02,su02)
               call fonct(xx,the12,ri12,gl12,su12)
               f2 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               fv1(jtwm1) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtwm1) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtwm1)*(f1+f2)
               
              end do
	    
    	res = resk*xr
	endif
	
          end subroutine qgaussPL31weib

!=============================================================================================  
!=====QGAUS31 out a 31 point Gauss-Kronrod quadrature rule for splines =======================
!=============================================================================================  

      subroutine qgaussPL31(a,b,the01,the12,the02,res,v1,v2,v3)

	use commun,only:zi01,zi12,zi02,nz01,nz12,nz02
 

         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resk,v1,v2,v3,&
         fv1,fv2,d1mach(5),epmach,uflow
         double precision,dimension(-2:(nz01-1))::the01
         double precision,dimension(-2:(nz12-1))::the12
         double precision,dimension(-2:(nz02-1))::the02
         double precision,dimension(16)::xgk,wgk
	 double precision,dimension(8)::wg
         double precision::xx,f1,su01,ri01,ri12,f2,su12,su02,ri02,fc,gl01,gl02,gl12
         save wgk,xgk
	 dimension fv1(15),fv2(15)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)
                                            
	wg(1)=0.0307532419961172683546283935772044177217d0                                   
   	wg(2)=0.0703660474881081247092674164506673384667d0                                   
    	wg(3)=0.107159220467171935011869546685869303416d0                                     
    	wg(4)=0.139570677926154314447804794511028322521d0                                    
        wg(5)=0.166269205816993933553200860481208811131d0
	wg(6)=0.186161000015562211026800561866422824506d0
        wg(7)=0.198431485327111576456118326443839324819d0
        wg(8)=0.202578241925561272880620199967519314839d0

                                             
    	xgk(1)=0.998002298693397060285172840152271209073d0
    	xgk(2)=0.987992518020485428489565718586612581147d0
    	xgk(3)=0.967739075679139134257347978784337225283d0
    	xgk(4)=0.937273392400705904307758947710209471244d0
    	xgk(5)=0.897264532344081900882509656454495882832d0
    	xgk(6)=0.848206583410427216200648320774216851366d0
    	xgk(7)=0.790418501442465932967649294817947346862d0
        xgk(8)=0.724417731360170047416186054613938009631d0
        xgk(9)=0.650996741297416970533735895313274692547d0
        xgk(10)=0.570972172608538847537226737253910641238d0
	xgk(11)=0.485081863640239680693655740232350612866d0
	xgk(12)=0.394151347077563369897207370981045468363d0
	xgk(13)=0.299180007153168812166780024266388962662d0
	xgk(14)=0.201194093997434522300628303394596207813d0
	xgk(15)=0.101142066918717499027074231447392338787d0
    	xgk(16)=0.000000000000000000000000000000000d0
	
    	wgk(1)=0.00537747987292334898779205143012764981831d0
    	wgk(2)=0.0150079473293161225383747630758072680946d0
    	wgk(3)=0.0254608473267153201868740010196533593973d0
    	wgk(4)=0.0353463607913758462220379484783600481226d0
    	wgk(5)=0.0445897513247648766082272993732796902233d0
    	wgk(6)=0.0534815246909280872653431472394302967716d0
    	wgk(7)=0.0620095678006706402851392309608029321904d0
    	wgk(8)=0.0698541213187282587095200770991474757860d0
	wgk(9)=0.0768496807577203788944327774826590067221d0
	wgk(10)=0.0830805028231330210382892472861037896016d0
	wgk(11)=0.0885644430562117706472754436937743032123d0
	wgk(12)=0.0931265981708253212254868727473457185619d0
	wgk(13)=0.0966427269836236785051799076275893351367d0
	wgk(14)=0.0991735987217919593323931734846031310596d0
	wgk(15)=0.100769845523875595044946662617569721916d0
	wgk(16)=0.101330007014791549017374792767492546771d0

        xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)
        call susp(xm,the01,nz01,su01,ri01,zi01,gl01)
        call susp(xm,the02,nz02,su02,ri02,zi02,gl02)
        call susp(xm,the12,nz12,su12,ri12,zi12,gl12)
        fc = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0

    	resk = fc*wgk(16)       ! init res Kronrod   ! fc * 8e poids Kronrod
         

         do j=1,7
	       jtw = j*2
               dx=xr*xgk(jtw)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f1 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f2 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               fv1(jtw) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtw) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtw)*(f1+f2)
               
         end do
	 do j=1,8
	       jtwm1 = j*2-1
               dx=xr*xgk(jtwm1)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f1 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f2 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               fv1(jtwm1) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtwm1) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtwm1)*(f1+f2)
               
         end do
    
    
    res = resk*xr
    
         end subroutine qgaussPL31

!=============================================================================================  
!=====QGAUS41 out a 41 point Gauss-Kronrod quadrature rule for weib ==========================
!=============================================================================================  

subroutine qgaussPL41weib(a,b,the01,the02,the12,res,v01,v02,v12)
         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resk,v01,v02,v12,&
         fv1,fv2,d1mach(5),epmach,uflow,the01(2),the12(2),the02(2)
         double precision,dimension(21)::xgk,wgk
	 double precision,dimension(10)::wg
         double precision::xx,f1,su01,ri01,ri12,f2,su12,su02,ri02,fc,gl01,gl02,gl12
         save wgk,xgk
	 dimension fv1(20),fv2(20)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)
                                            
	wg(1)=0.0176140071391521183118619623518528163621d0                                   
   	wg(2)=0.0406014298003869413310399522749321098791d0                                   
    	wg(3)=0.0626720483341090635695065351870416063516d0                                     
    	wg(4)=0.0832767415767047487247581432220462061002d0                                    
        wg(5)=0.101930119817240435036750135480349876167d0
	wg(6)=0.118194531961518417312377377711382287005d0
        wg(7)=0.131688638449176626898494499748163134916d0
        wg(8)=0.142096109318382051329298325067164933035d0
        wg(9)=0.149172986472603746787828737001969436693d0
        wg(10)=0.152753387130725850698084331955097593492d0


                                             
    	xgk(1)=0.998859031588277663838315576545863010000d0
    	xgk(2)=0.993128599185094924786122388471320278223d0
    	xgk(3)=0.981507877450250259193342994720216944567d0
    	xgk(4)=0.963971927277913791267666131197277221912d0
    	xgk(5)=0.940822633831754753519982722212443380274d0
    	xgk(6)=0.912234428251325905867752441203298113049d0
    	xgk(7)=0.878276811252281976077442995113078466711d0
        xgk(8)=0.839116971822218823394529061701520685330d0
        xgk(9)=0.795041428837551198350638833272787942959d0
        xgk(10)=0.746331906460150792614305070355641590311d0
	xgk(11)=0.693237656334751384805490711845931533386d0
	xgk(12)=0.636053680726515025452836696226285936743d0
	xgk(13)=0.575140446819710315342946036586425132814d0
	xgk(14)=0.510867001950827098004364050955250998425d0
	xgk(15)=0.443593175238725103199992213492640107840d0
	xgk(16)=0.373706088715419560672548177024927237396d0
	xgk(17)=0.301627868114913004320555356858592260615d0
	xgk(18)=0.227785851141645078080496195368574624743d0
	xgk(19)=0.152605465240922675505220241022677527912d0
	xgk(20)=0.765265211334973337546404093988382110048d0
    	xgk(21)=0.000000000000000000000000000000000d0
	
    	wgk(1)=0.00307358371852053150121829324603098748803d0
    	wgk(2)=0.00860026985564294219866178795010234725213d0
    	wgk(3)=0.0146261692569712529837879603088683561639d0
    	wgk(4)=0.0203883734612665235980102314327547051228d0
    	wgk(5)=0.0258821336049511588345050670961531429995d0
    	wgk(6)=0.0312873067770327989585431193238007378878d0
    	wgk(7)=0.0366001697582007980305572407072110084875d0
    	wgk(8)=0.0416688733279736862637883059368947380440d0
	wgk(9)=0.0464348218674976747202318809261075168421d0
	wgk(10)=0.0509445739237286919327076700503449486648d0
	wgk(11)=0.0551951053482859947448323724197773291948d0
	wgk(12)=0.0591114008806395723749672206485942171364d0
	wgk(13)=0.0626532375547811680258701221742549805858d0
	wgk(14)=0.0658345971336184221115635569693979431472d0
	wgk(15)=0.0686486729285216193456234118853678017155d0
	wgk(16)=0.0710544235534440683057903617232101674129d0
	wgk(17)=0.0730306903327866674951894176589131127606d0
	wgk(18)=0.0745828754004991889865814183624875286161d0
	wgk(19)=0.0757044976845566746595427753766165582634d0
	wgk(20)=0.0763778676720807367055028350380610018008d0
	wgk(21)=0.0766007119179996564450499015301017408279d0

        xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)
        call fonct(xm,the01,ri01,gl01,su01)
        call fonct(xm,the02,ri02,gl02,su02)
        call fonct(xm,the12,ri12,gl12,su12)
        fc = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0

        resk = fc*wgk(21)       ! init res Kronrod   ! fc * 8e poids Kronrod
        
            if(a.eq.b)then
               res = 0.d0
            else
               do j=1,10
               	jtw = j*2
               	dx=xr*xgk(jtw)
               	xx = xm+dx
               	call fonct(xx,the01,ri01,gl01,su01)
               	call fonct(xx,the02,ri02,gl02,su02)
               	call fonct(xx,the12,ri12,gl12,su12)
               	f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               	xx = xm-dx
               	call fonct(xx,the01,ri01,gl01,su01)
               	call fonct(xx,the02,ri02,gl02,su02)
               	call fonct(xx,the12,ri12,gl12,su12)
               	f2 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               	fv1(jtw) = f1   ! svgrd valeurs fct f a gche du centre
               	fv2(jtw) = f2   ! svgrd valeurs fct f a drte du centre
	       	resk = resk + wgk(jtw)*(f1+f2)
                
              end do
	      do j=1,10
	       jtwm1 = j*2-1
               dx=xr*xgk(jtwm1)
               xx = xm+dx
               call fonct(xx,the01,ri01,gl01,su01)
               call fonct(xx,the02,ri02,gl02,su02)
               call fonct(xx,the12,ri12,gl12,su12)
               f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               xx = xm-dx
               call fonct(xx,the01,ri01,gl01,su01)
               call fonct(xx,the02,ri02,gl02,su02)
               call fonct(xx,the12,ri12,gl12,su12)
               f2 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               fv1(jtwm1) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtwm1) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtwm1)*(f1+f2)
               
              end do
	    
    	res = resk*xr
	endif
    
          end subroutine qgaussPL41weib

!=============================================================================================  
!=====QGAUS41 out a 41 point Gauss-Kronrod quadrature rule for splines =======================
!=============================================================================================  

      subroutine qgaussPL41(a,b,the01,the12,the02,res,v1,v2,v3)

	use commun,only:zi01,zi12,zi02,nz01,nz12,nz02
         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resk,v1,v2,v3,&
         fv1,fv2,d1mach(5),epmach,uflow
         double precision,dimension(-2:(nz01-1))::the01
         double precision,dimension(-2:(nz12-1))::the12
         double precision,dimension(-2:(nz02-1))::the02
         double precision,dimension(21)::xgk,wgk
	 double precision,dimension(10)::wg
         double precision::xx,f1,su01,ri01,ri12,f2,su12,su02,ri02,fc,gl01,gl02,gl12
         save wgk,xgk
	 dimension fv1(20),fv2(20)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)
                                            
	wg(1)=0.0176140071391521183118619623518528163621d0                                   
   	wg(2)=0.0406014298003869413310399522749321098791d0                                   
    	wg(3)=0.0626720483341090635695065351870416063516d0                                     
    	wg(4)=0.0832767415767047487247581432220462061002d0                                    
        wg(5)=0.101930119817240435036750135480349876167d0
	wg(6)=0.118194531961518417312377377711382287005d0
        wg(7)=0.131688638449176626898494499748163134916d0
        wg(8)=0.142096109318382051329298325067164933035d0
        wg(9)=0.149172986472603746787828737001969436693d0
        wg(10)=0.152753387130725850698084331955097593492d0


                                             
    	xgk(1)=0.998859031588277663838315576545863010000d0
    	xgk(2)=0.993128599185094924786122388471320278223d0
    	xgk(3)=0.981507877450250259193342994720216944567d0
    	xgk(4)=0.963971927277913791267666131197277221912d0
    	xgk(5)=0.940822633831754753519982722212443380274d0
    	xgk(6)=0.912234428251325905867752441203298113049d0
    	xgk(7)=0.878276811252281976077442995113078466711d0
        xgk(8)=0.839116971822218823394529061701520685330d0
        xgk(9)=0.795041428837551198350638833272787942959d0
        xgk(10)=0.746331906460150792614305070355641590311d0
	xgk(11)=0.693237656334751384805490711845931533386d0
	xgk(12)=0.636053680726515025452836696226285936743d0
	xgk(13)=0.575140446819710315342946036586425132814d0
	xgk(14)=0.510867001950827098004364050955250998425d0
	xgk(15)=0.443593175238725103199992213492640107840d0
	xgk(16)=0.373706088715419560672548177024927237396d0
	xgk(17)=0.301627868114913004320555356858592260615d0
	xgk(18)=0.227785851141645078080496195368574624743d0
	xgk(19)=0.152605465240922675505220241022677527912d0
	xgk(20)=0.765265211334973337546404093988382110048d0
    	xgk(21)=0.000000000000000000000000000000000d0
	
    	wgk(1)=0.00307358371852053150121829324603098748803d0
    	wgk(2)=0.00860026985564294219866178795010234725213d0
    	wgk(3)=0.0146261692569712529837879603088683561639d0
    	wgk(4)=0.0203883734612665235980102314327547051228d0
    	wgk(5)=0.0258821336049511588345050670961531429995d0
    	wgk(6)=0.0312873067770327989585431193238007378878d0
    	wgk(7)=0.0366001697582007980305572407072110084875d0
    	wgk(8)=0.0416688733279736862637883059368947380440d0
	wgk(9)=0.0464348218674976747202318809261075168421d0
	wgk(10)=0.0509445739237286919327076700503449486648d0
	wgk(11)=0.0551951053482859947448323724197773291948d0
	wgk(12)=0.0591114008806395723749672206485942171364d0
	wgk(13)=0.0626532375547811680258701221742549805858d0
	wgk(14)=0.0658345971336184221115635569693979431472d0
	wgk(15)=0.0686486729285216193456234118853678017155d0
	wgk(16)=0.0710544235534440683057903617232101674129d0
	wgk(17)=0.0730306903327866674951894176589131127606d0
	wgk(18)=0.0745828754004991889865814183624875286161d0
	wgk(19)=0.0757044976845566746595427753766165582634d0
	wgk(20)=0.0763778676720807367055028350380610018008d0
	wgk(21)=0.0766007119179996564450499015301017408279d0

        xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)
        call susp(xm,the01,nz01,su01,ri01,zi01,gl01)
        call susp(xm,the02,nz02,su02,ri02,zi02,gl02)
        call susp(xm,the12,nz12,su12,ri12,zi12,gl12)
        fc = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0

    	resk = fc*wgk(21)       ! init res Kronrod   ! fc * 8e poids Kronrod
        
         do j=1,10
	       jtw = j*2
               dx=xr*xgk(jtw)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f1 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f2 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               fv1(jtw) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtw) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtw)*(f1+f2)
              
         end do
	 do j=1,10
	       jtwm1 = j*2-1
               dx=xr*xgk(jtwm1)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f1 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f2 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               fv1(jtwm1) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtwm1) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtwm1)*(f1+f2)
               
         end do
     
    res = resk*xr
         end subroutine qgaussPL41

!=============================================================================================  
!=====QGAUS51 out a 51 point Gauss-Kronrod quadrature rule for weib ==========================
!=============================================================================================  


subroutine qgaussPL51weib(a,b,the01,the02,the12,res,v01,v02,v12)
         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resk,v01,v02,v12,&
         fv1,fv2,d1mach(5),epmach,uflow,the01(2),the12(2),the02(2)
         double precision,dimension(26)::xgk,wgk
	 double precision,dimension(13)::wg
         double precision::xx,f1,su01,ri01,ri12,f2,su12,su02,ri02,fc,gl01,gl02,gl12
         save wgk,xgk
	 dimension fv1(25),fv2(25)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)
                                            
	wg(1)=0.0113937985010262879479029641132347736033d0                                   
   	wg(2)=0.0263549866150321372619018152952991449360d0                                   
    	wg(3)=0.0409391567013063126556234877116459536608d0                                     
    	wg(4)=0.0549046959758351919259368915404733241601d0                                    
        wg(5)=0.0680383338123569172071871856567079685547d0
	wg(6)=0.0801407003350010180132349596691113022902d0
        wg(7)=0.0910282619829636498114972207028916533810d0
        wg(8)=0.100535949067050644202206890392685826988d0
        wg(9)=0.108519624474263653116093957050116619340d0
        wg(10)=0.114858259145711648339325545869555808641d0
	wg(11)=0.119455763535784772228178126512901047390d0
	wg(12)=0.122242442990310041688959518945851505835d0
	wg(13)=0.123176053726715451203902873079050142438d0


                                             
    	xgk(1)=0.999262104992609834193457486540340593705d0
    	xgk(2)=0.995556969790498097908784946893901617258d0
    	xgk(3)=0.988035794534077247637331014577406227072d0
    	xgk(4)=0.976663921459517511498315386479594067745d0
    	xgk(5)=0.961614986425842512418130033660167241692d0
    	xgk(6)=0.942974571228974339414011169658470531905d0
    	xgk(7)=0.920747115281701561746346084546330631575d0
        xgk(8)=0.894991997878275368851042006782804954175d0
        xgk(9)=0.865847065293275595448996969588340088203d0
        xgk(10)=0.833442628760834001421021108693569569461d0
	xgk(11)=0.797873797998500059410410904994306569409d0
	xgk(12)=0.759259263037357630577282865204360976388d0
	xgk(13)=0.717766406813084388186654079773297780598d0
	xgk(14)=0.673566368473468364485120633247622175883d0
	xgk(15)=0.626810099010317412788122681624517881020d0
	xgk(16)=0.577662930241222967723689841612654067396d0
	xgk(17)=0.526325284334719182599623778158010178037d0
	xgk(18)=0.473002731445714960522182115009192041332d0
	xgk(19)=0.417885382193037748851814394594572487093d0
	xgk(20)=0.361172305809387837735821730127640667422d0
	xgk(21)=0.303089538931107830167478909980339329200d0
	xgk(22)=0.243866883720988432045190362797451586406d0
	xgk(23)=0.183718939421048892015969888759528415785d0
	xgk(24)=0.122864692610710396387359818808036805532d0
	xgk(25)=0.0615444830056850788865463923667966312817d0
    	xgk(26)=0.000000000000000000000000000000000d0
	
    	wgk(1)=0.00198738389233031592650785188284340988943d0
    	wgk(2)=0.00556193213535671375804023690106552207018d0
    	wgk(3)=0.00947397338617415160720771052365532387165d0
    	wgk(4)=0.0132362291955716748136564058469762380776d0
    	wgk(5)=0.0168478177091282982315166675363363158404d0
    	wgk(6)=0.0204353711458828354565682922359389736788d0
    	wgk(7)=0.0240099456069532162200924891648810813929d0
    	wgk(8)=0.0274753175878517378029484555178110786148d0
	wgk(9)=0.0307923001673874888911090202152285856009d0
	wgk(10)=0.0340021302743293378367487952295512032257d0
	wgk(11)=0.0371162714834155435603306253676198759960d0
	wgk(12)=0.0400838255040323820748392844670756464014d0
	wgk(13)=0.0428728450201700494768957924394951611020d0
	wgk(14)=0.0455029130499217889098705847526603930437d0
	wgk(15)=0.0479825371388367139063922557569147549836d0
	wgk(16)=0.0502776790807156719633252594334400844406d0
	wgk(17)=0.0523628858064074758643667121378727148874d0
	wgk(18)=0.0542511298885454901445433704598756068261d0
	wgk(19)=0.0559508112204123173082406863827473468203d0
	wgk(20)=0.0574371163615678328535826939395064719948d0
	wgk(21)=0.0586896800223942079619741758567877641398d0
	wgk(22)=0.0597203403241740599790992919325618538354d0
	wgk(23)=0.0605394553760458629453602675175654271623d0
	wgk(24)=0.0611285097170530483058590304162927119227d0
	wgk(25)=0.0614711898714253166615441319652641775865d0
	wgk(26)=0.0615808180678329350787598242400645531904d0

        xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)
        call fonct(xm,the01,ri01,gl01,su01)
        call fonct(xm,the02,ri02,gl02,su02)
        call fonct(xm,the12,ri12,gl12,su12)
        fc = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0

    	resk = fc*wgk(26)       ! init res Kronrod   ! fc * 8e poids Kronrod
         
            if(a.eq.b)then
               res = 0.d0
            else
               do j=1,12
               	jtw = j*2
               	dx=xr*xgk(jtw)
               	xx = xm+dx
               	call fonct(xx,the01,ri01,gl01,su01)
               	call fonct(xx,the02,ri02,gl02,su02)
               	call fonct(xx,the12,ri12,gl12,su12)
               	f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               	xx = xm-dx
               	call fonct(xx,the01,ri01,gl01,su01)
               	call fonct(xx,the02,ri02,gl02,su02)
               	call fonct(xx,the12,ri12,gl12,su12)
               	f2 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               	fv1(jtw) = f1   ! svgrd valeurs fct f a gche du centre
               	fv2(jtw) = f2   ! svgrd valeurs fct f a drte du centre
	       	resk = resk + wgk(jtw)*(f1+f2)
                
              end do
	      do j=1,13
	       jtwm1 = j*2-1
               dx=xr*xgk(jtwm1)
               xx = xm+dx
               call fonct(xx,the01,ri01,gl01,su01)
               call fonct(xx,the02,ri02,gl02,su02)
               call fonct(xx,the12,ri12,gl12,su12)
               f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               xx = xm-dx
               call fonct(xx,the01,ri01,gl01,su01)
               call fonct(xx,the02,ri02,gl02,su02)
               call fonct(xx,the12,ri12,gl12,su12)
               f2 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               fv1(jtwm1) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtwm1) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtwm1)*(f1+f2)
               
              end do
    
    	res = resk*xr
	endif
    
          end subroutine qgaussPL51weib

!=============================================================================================  
!=====QGAUS51 out a 51 point Gauss-Kronrod quadrature rule for splines =======================
!=============================================================================================  


      subroutine qgaussPL51(a,b,the01,the12,the02,res,v1,v2,v3)

	use commun,only:zi01,zi12,zi02,nz01,nz12,nz02

         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resk,v1,v2,v3,&
         fv1,fv2,d1mach(5),epmach,uflow
         double precision,dimension(-2:(nz01-1))::the01
         double precision,dimension(-2:(nz12-1))::the12
         double precision,dimension(-2:(nz02-1))::the02
         double precision,dimension(26)::xgk,wgk
	 double precision,dimension(13)::wg
         double precision::xx,f1,su01,ri01,ri12,f2,su12,su02,ri02,fc,gl01,gl02,gl12
         save wgk,xgk
	 dimension fv1(25),fv2(25)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)
                                            
	wg(1)=0.0113937985010262879479029641132347736033d0                                   
   	wg(2)=0.0263549866150321372619018152952991449360d0                                   
    	wg(3)=0.0409391567013063126556234877116459536608d0                                     
    	wg(4)=0.0549046959758351919259368915404733241601d0                                    
        wg(5)=0.0680383338123569172071871856567079685547d0
	wg(6)=0.0801407003350010180132349596691113022902d0
        wg(7)=0.0910282619829636498114972207028916533810d0
        wg(8)=0.100535949067050644202206890392685826988d0
        wg(9)=0.108519624474263653116093957050116619340d0
        wg(10)=0.114858259145711648339325545869555808641d0
	wg(11)=0.119455763535784772228178126512901047390d0
	wg(12)=0.122242442990310041688959518945851505835d0
	wg(13)=0.123176053726715451203902873079050142438d0


                                             
    	xgk(1)=0.999262104992609834193457486540340593705d0
    	xgk(2)=0.995556969790498097908784946893901617258d0
    	xgk(3)=0.988035794534077247637331014577406227072d0
    	xgk(4)=0.976663921459517511498315386479594067745d0
    	xgk(5)=0.961614986425842512418130033660167241692d0
    	xgk(6)=0.942974571228974339414011169658470531905d0
    	xgk(7)=0.920747115281701561746346084546330631575d0
        xgk(8)=0.894991997878275368851042006782804954175d0
        xgk(9)=0.865847065293275595448996969588340088203d0
        xgk(10)=0.833442628760834001421021108693569569461d0
	xgk(11)=0.797873797998500059410410904994306569409d0
	xgk(12)=0.759259263037357630577282865204360976388d0
	xgk(13)=0.717766406813084388186654079773297780598d0
	xgk(14)=0.673566368473468364485120633247622175883d0
	xgk(15)=0.626810099010317412788122681624517881020d0
	xgk(16)=0.577662930241222967723689841612654067396d0
	xgk(17)=0.526325284334719182599623778158010178037d0
	xgk(18)=0.473002731445714960522182115009192041332d0
	xgk(19)=0.417885382193037748851814394594572487093d0
	xgk(20)=0.361172305809387837735821730127640667422d0
	xgk(21)=0.303089538931107830167478909980339329200d0
	xgk(22)=0.243866883720988432045190362797451586406d0
	xgk(23)=0.183718939421048892015969888759528415785d0
	xgk(24)=0.122864692610710396387359818808036805532d0
	xgk(25)=0.0615444830056850788865463923667966312817d0
    	xgk(26)=0.000000000000000000000000000000000d0
	
    	wgk(1)=0.00198738389233031592650785188284340988943d0
    	wgk(2)=0.00556193213535671375804023690106552207018d0
    	wgk(3)=0.00947397338617415160720771052365532387165d0
    	wgk(4)=0.0132362291955716748136564058469762380776d0
    	wgk(5)=0.0168478177091282982315166675363363158404d0
    	wgk(6)=0.0204353711458828354565682922359389736788d0
    	wgk(7)=0.0240099456069532162200924891648810813929d0
    	wgk(8)=0.0274753175878517378029484555178110786148d0
	wgk(9)=0.0307923001673874888911090202152285856009d0
	wgk(10)=0.0340021302743293378367487952295512032257d0
	wgk(11)=0.0371162714834155435603306253676198759960d0
	wgk(12)=0.0400838255040323820748392844670756464014d0
	wgk(13)=0.0428728450201700494768957924394951611020d0
	wgk(14)=0.0455029130499217889098705847526603930437d0
	wgk(15)=0.0479825371388367139063922557569147549836d0
	wgk(16)=0.0502776790807156719633252594334400844406d0
	wgk(17)=0.0523628858064074758643667121378727148874d0
	wgk(18)=0.0542511298885454901445433704598756068261d0
	wgk(19)=0.0559508112204123173082406863827473468203d0
	wgk(20)=0.0574371163615678328535826939395064719948d0
	wgk(21)=0.0586896800223942079619741758567877641398d0
	wgk(22)=0.0597203403241740599790992919325618538354d0
	wgk(23)=0.0605394553760458629453602675175654271623d0
	wgk(24)=0.0611285097170530483058590304162927119227d0
	wgk(25)=0.0614711898714253166615441319652641775865d0
	wgk(26)=0.0615808180678329350787598242400645531904d0

        xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)
        call susp(xm,the01,nz01,su01,ri01,zi01,gl01)
        call susp(xm,the02,nz02,su02,ri02,zi02,gl02)
        call susp(xm,the12,nz12,su12,ri12,zi12,gl12)
        fc = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0

    	resk = fc*wgk(26)       ! init res Kronrod   ! fc * 8e poids Kronrod
        
         do j=1,12
	       jtw = j*2
               dx=xr*xgk(jtw)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f1 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f2 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               fv1(jtw) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtw) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtw)*(f1+f2)
               
         end do
	 do j=1,13
	       jtwm1 = j*2-1
               dx=xr*xgk(jtwm1)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f1 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f2 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               fv1(jtwm1) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtwm1) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtwm1)*(f1+f2)
               
         end do
  
    
    res = resk*xr
      end subroutine qgaussPL51

!=============================================================================================  
!=====QGAUS61 out a 61 point Gauss-Kronrod quadrature rule for weib ==========================
!=============================================================================================  


subroutine qgaussPL61weib(a,b,the01,the02,the12,res,v01,v02,v12)
         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resk,v01,v02,v12,&
         fv1,fv2,d1mach(5),epmach,uflow,the01(2),the12(2),the02(2)
         double precision,dimension(31)::xgk,wgk
	 double precision,dimension(15)::wg
         double precision::xx,f1,su01,ri01,ri12,f2,su12,su02,ri02,fc,gl01,gl02,gl12
         save wgk,xgk
	 dimension fv1(30),fv2(30)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)
                                            
	wg(1)=0.00796819249616660561546588347467362245048d0                                   
   	wg(2)=0.0184664683110909591423021319120472690962d0                                   
    	wg(3)=0.0287847078833233693497191796112920436396d0                                     
    	wg(4)=0.0387991925696270495968019364463476920332d0                                    
        wg(5)=0.0484026728305940529029381404228075178153d0
	wg(6)=0.0574931562176190664817216894020561287971d0
        wg(7)=0.0659742298821804951281285151159623612374d0
        wg(8)=0.0737559747377052062682438500221907341538d0
        wg(9)=0.0807558952294202153546949384605297308759d0
        wg(10)=0.0868997872010829798023875307151257025768d0
	wg(11)=0.0921225222377861287176327070876187671969d0
	wg(12)=0.0963687371746442596394686263518098650964d0
	wg(13)=0.0995934205867952670627802821035694765299d0
	wg(14)=0.101762389748405504596428952168554044633d0
	wg(15)=0.102852652893558840341285636705415043868d0


                                             
    	xgk(1)=0.999484410050490637571325895705810819469d0
    	xgk(2)=0.996893484074649540271630050918695283341d0
    	xgk(3)=0.991630996870404594858628366109485724851d0
    	xgk(4)=0.983668123279747209970032581605662801940d0
    	xgk(5)=0.973116322501126268374693868423706884888d0
    	xgk(6)=0.960021864968307512216871025581797662930d0
    	xgk(7)=0.944374444748559979415831324037439121586d0
        xgk(8)=0.926200047429274325879324277080474004086d0
        xgk(9)=0.905573307699907798546522558925958319569d0
        xgk(10)=0.882560535792052681543116462530225590057d0
	xgk(11)=0.857205233546061098958658510658943856821d0
	xgk(12)=0.829565762382768397442898119732501916439d0
	xgk(13)=0.799727835821839083013668942322683240736d0
	xgk(14)=0.767777432104826194917977340974503131695d0
	xgk(15)=0.733790062453226804726171131369527645669d0
	xgk(16)=0.697850494793315796932292388026640068382d0
	xgk(17)=0.660061064126626961370053668149270753038d0
	xgk(18)=0.620526182989242861140477556431189299207d0
	xgk(19)=0.579345235826361691756024932172540495907d0
	xgk(20)=0.536624148142019899264169793311072794164d0
	xgk(21)=0.492480467861778574993693061207708795644d0
	xgk(22)=0.447033769538089176780609900322854000162d0
	xgk(23)=0.400401254830394392535476211542660633611d0
	xgk(24)=0.352704725530878113471037207089373860654d0
	xgk(25)=0.304073202273625077372677107199256553531d0
	xgk(26)=0.254636926167889846439805129817805107883d0
	xgk(27)=0.204525116682309891438957671002024709524d0
	xgk(28)=0.153869913608583546963794672743255920419d0
	xgk(29)=0.102806937966737030147096751318000592472d0
	xgk(30)=0.0514718425553176958330252131667225737491d0
    	xgk(31)=0.000000000000000000000000000000000d0
	
    	wgk(1)=0.00138901369867700762455159122675969968105d0
    	wgk(2)=0.00389046112709988405126720184451550327852d0
    	wgk(3)=0.00663070391593129217331982636975016813363d0
    	wgk(4)=0.00927327965951776342844114689202436042127d0
    	wgk(5)=0.0118230152534963417422328988532505928963d0
    	wgk(6)=0.0143697295070458048124514324435800101958d0
    	wgk(7)=0.0169208891890532726275722894203220923686d0
    	wgk(8)=0.0194141411939423811734089510501284558514d0
	wgk(9)=0.0218280358216091922971674857383389934015d0
	wgk(10)=0.0241911620780806013656863707252320267604d0
	wgk(11)=0.0265099548823331016106017093350754143665d0
	wgk(12)=0.0287540487650412928439787853543342111447d0
	wgk(13)=0.0309072575623877624728842529430922726353d0
	wgk(14)=0.0329814470574837260318141910168539275106d0
	wgk(15)=0.0349793380280600241374996707314678750972d0
	wgk(16)=0.0368823646518212292239110656171359677370d0
	wgk(17)=0.0386789456247275929503486515322810502509d0
	wgk(18)=0.0403745389515359591119952797524681142161d0
	wgk(19)=0.0419698102151642461471475412859697577901d0
	wgk(20)=0.0434525397013560693168317281170732580746d0
	wgk(21)=0.0448148001331626631923555516167232437574d0
	wgk(22)=0.0460592382710069881162717355593735805947d0
	wgk(23)=0.0471855465692991539452614781810994864829d0
	wgk(24)=0.0481858617570871291407794922983045926058d0
	wgk(25)=0.0490554345550297788875281653672381736059d0
	wgk(26)=0.0497956834270742063578115693799423285392d0
	wgk(27)=0.0504059214027823468408930856535850289022d0
	wgk(28)=0.0508817958987496064922974730498046918534d0
	wgk(29)=0.0512215478492587721706562826049442082511d0
	wgk(30)=0.0514261285374590259338628792157812598296d0
	wgk(31)=0.0514947294294515675583404336470993075327d0

        xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)
        call fonct(xm,the01,ri01,gl01,su01)
        call fonct(xm,the02,ri02,gl02,su02)
        call fonct(xm,the12,ri12,gl12,su12)
        fc = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0

    	resk = fc*wgk(31)       ! init res Kronrod   ! fc * 8e poids Kronrod
        
            if(a.eq.b)then
               res = 0.d0
            else
               do j=1,15
               	jtw = j*2
               	dx=xr*xgk(jtw)
               	xx = xm+dx
               	call fonct(xx,the01,ri01,gl01,su01)
               	call fonct(xx,the02,ri02,gl02,su02)
               	call fonct(xx,the12,ri12,gl12,su12)
               	f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               	xx = xm-dx
               	call fonct(xx,the01,ri01,gl01,su01)
               	call fonct(xx,the02,ri02,gl02,su02)
               	call fonct(xx,the12,ri12,gl12,su12)
               	f2 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               	fv1(jtw) = f1   ! svgrd valeurs fct f a gche du centre
               	fv2(jtw) = f2   ! svgrd valeurs fct f a drte du centre
	       	resk = resk + wgk(jtw)*(f1+f2)
              end do
	      do j=1,15
	       jtwm1 = j*2-1
               dx=xr*xgk(jtwm1)
               xx = xm+dx
               call fonct(xx,the01,ri01,gl01,su01)
               call fonct(xx,the02,ri02,gl02,su02)
               call fonct(xx,the12,ri12,gl12,su12)
               f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               xx = xm-dx
               call fonct(xx,the01,ri01,gl01,su01)
               call fonct(xx,the02,ri02,gl02,su02)
               call fonct(xx,the12,ri12,gl12,su12)
               f2 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
               fv1(jtwm1) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtwm1) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtwm1)*(f1+f2)
              end do

    	res = resk*xr
	endif
    
          end subroutine qgaussPL61weib

!=============================================================================================  
!=====QGAUS61 out a 61 point Gauss-Kronrod quadrature rule for splines =======================
!=============================================================================================  


      subroutine qgaussPL61(a,b,the01,the12,the02,res,v1,v2,v3)

	use commun,only:zi01,zi12,zi02,nz01,nz12,nz02

         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resk,v1,v2,v3,&
         fv1,fv2,d1mach(5),epmach,uflow
         double precision,dimension(-2:(nz01-1))::the01
         double precision,dimension(-2:(nz12-1))::the12
         double precision,dimension(-2:(nz02-1))::the02
         double precision,dimension(31)::xgk,wgk
	 double precision,dimension(15)::wg
         double precision::xx,f1,su01,ri01,ri12,f2,su12,su02,ri02,fc,gl01,gl02,gl12
         save wgk,xgk
	 dimension fv1(30),fv2(30)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)
                                            
	wg(1)=0.00796819249616660561546588347467362245048d0                                   
   	wg(2)=0.0184664683110909591423021319120472690962d0                                   
    	wg(3)=0.0287847078833233693497191796112920436396d0                                     
    	wg(4)=0.0387991925696270495968019364463476920332d0                                    
        wg(5)=0.0484026728305940529029381404228075178153d0
	wg(6)=0.0574931562176190664817216894020561287971d0
        wg(7)=0.0659742298821804951281285151159623612374d0
        wg(8)=0.0737559747377052062682438500221907341538d0
        wg(9)=0.0807558952294202153546949384605297308759d0
        wg(10)=0.0868997872010829798023875307151257025768d0
	wg(11)=0.0921225222377861287176327070876187671969d0
	wg(12)=0.0963687371746442596394686263518098650964d0
	wg(13)=0.0995934205867952670627802821035694765299d0
	wg(14)=0.101762389748405504596428952168554044633d0
	wg(15)=0.102852652893558840341285636705415043868d0


                                             
    	xgk(1)=0.999484410050490637571325895705810819469d0
    	xgk(2)=0.996893484074649540271630050918695283341d0
    	xgk(3)=0.991630996870404594858628366109485724851d0
    	xgk(4)=0.983668123279747209970032581605662801940d0
    	xgk(5)=0.973116322501126268374693868423706884888d0
    	xgk(6)=0.960021864968307512216871025581797662930d0
    	xgk(7)=0.944374444748559979415831324037439121586d0
        xgk(8)=0.926200047429274325879324277080474004086d0
        xgk(9)=0.905573307699907798546522558925958319569d0
        xgk(10)=0.882560535792052681543116462530225590057d0
	xgk(11)=0.857205233546061098958658510658943856821d0
	xgk(12)=0.829565762382768397442898119732501916439d0
	xgk(13)=0.799727835821839083013668942322683240736d0
	xgk(14)=0.767777432104826194917977340974503131695d0
	xgk(15)=0.733790062453226804726171131369527645669d0
	xgk(16)=0.697850494793315796932292388026640068382d0
	xgk(17)=0.660061064126626961370053668149270753038d0
	xgk(18)=0.620526182989242861140477556431189299207d0
	xgk(19)=0.579345235826361691756024932172540495907d0
	xgk(20)=0.536624148142019899264169793311072794164d0
	xgk(21)=0.492480467861778574993693061207708795644d0
	xgk(22)=0.447033769538089176780609900322854000162d0
	xgk(23)=0.400401254830394392535476211542660633611d0
	xgk(24)=0.352704725530878113471037207089373860654d0
	xgk(25)=0.304073202273625077372677107199256553531d0
	xgk(26)=0.254636926167889846439805129817805107883d0
	xgk(27)=0.204525116682309891438957671002024709524d0
	xgk(28)=0.153869913608583546963794672743255920419d0
	xgk(29)=0.102806937966737030147096751318000592472d0
	xgk(30)=0.0514718425553176958330252131667225737491d0
    	xgk(31)=0.000000000000000000000000000000000d0
	
    	wgk(1)=0.00138901369867700762455159122675969968105d0
    	wgk(2)=0.00389046112709988405126720184451550327852d0
    	wgk(3)=0.00663070391593129217331982636975016813363d0
    	wgk(4)=0.00927327965951776342844114689202436042127d0
    	wgk(5)=0.0118230152534963417422328988532505928963d0
    	wgk(6)=0.0143697295070458048124514324435800101958d0
    	wgk(7)=0.0169208891890532726275722894203220923686d0
    	wgk(8)=0.0194141411939423811734089510501284558514d0
	wgk(9)=0.0218280358216091922971674857383389934015d0
	wgk(10)=0.0241911620780806013656863707252320267604d0
	wgk(11)=0.0265099548823331016106017093350754143665d0
	wgk(12)=0.0287540487650412928439787853543342111447d0
	wgk(13)=0.0309072575623877624728842529430922726353d0
	wgk(14)=0.0329814470574837260318141910168539275106d0
	wgk(15)=0.0349793380280600241374996707314678750972d0
	wgk(16)=0.0368823646518212292239110656171359677370d0
	wgk(17)=0.0386789456247275929503486515322810502509d0
	wgk(18)=0.0403745389515359591119952797524681142161d0
	wgk(19)=0.0419698102151642461471475412859697577901d0
	wgk(20)=0.0434525397013560693168317281170732580746d0
	wgk(21)=0.0448148001331626631923555516167232437574d0
	wgk(22)=0.0460592382710069881162717355593735805947d0
	wgk(23)=0.0471855465692991539452614781810994864829d0
	wgk(24)=0.0481858617570871291407794922983045926058d0
	wgk(25)=0.0490554345550297788875281653672381736059d0
	wgk(26)=0.0497956834270742063578115693799423285392d0
	wgk(27)=0.0504059214027823468408930856535850289022d0
	wgk(28)=0.0508817958987496064922974730498046918534d0
	wgk(29)=0.0512215478492587721706562826049442082511d0
	wgk(30)=0.0514261285374590259338628792157812598296d0
	wgk(31)=0.0514947294294515675583404336470993075327d0

        xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)
        call susp(xm,the01,nz01,su01,ri01,zi01,gl01)
        call susp(xm,the02,nz02,su02,ri02,zi02,gl02)
        call susp(xm,the12,nz12,su12,ri12,zi12,gl12)
        fc = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0

    	 resk = fc*wgk(31)       ! init res Kronrod   ! fc * 8e poids Kronrod
         

         do j=1,15
	       jtw = j*2
               dx=xr*xgk(jtw)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f1 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f2 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               fv1(jtw) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtw) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtw)*(f1+f2)
               
         end do
	 do j=1,15
	       jtwm1 = j*2-1
               dx=xr*xgk(jtwm1)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f1 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
               call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
               call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
               f2 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               fv1(jtwm1) = f1   ! svgrd valeurs fct f a gche du centre
               fv2(jtwm1) = f2   ! svgrd valeurs fct f a drte du centre
	       resk = resk + wgk(jtwm1)*(f1+f2)
              
         end do
    
    
    res = resk*xr
    
         end subroutine qgaussPL61
         
!=============================================================================================        
!================================  qgaussweibderiv  ==========================================
!=== for derivatives approximation out a 15 point Gauss-Kronrod quadrature rule for weib =====
!=============================================================================================  


subroutine qgaussweibderiv(a,b,the01,the02,the12,resdenum,&
res01num,res02num,res12num,res0101num,res0102num,res0112num,&
res0202num,res0212num,res1212num,v01,v02,v12)

        implicit none
         double precision a,b,the01(2),the02(2),the12(2)
         double precision dx,xm,xr,reskdenum,&
         resdenum,resk01num,res01num,resk02num,res02num, & 
	resk12num,res12num,resk0101num,res0101num,resk0102num,res0102num, &
	resk0112num,res0112num,resk0202num,res0202num, &
	resk0212num,res0212num,resk1212num,res1212num, &
	v01,v02,v12
         double precision xx,f1denum,f2denum, f101num, f102num, f112num, & 
	f10101num,f10102num,f10112num,f10202num,f10212num,f11212num,&
	f201num, f202num, f212num, f20101num,f20102num,f20112num,&
	f20202num,f20212num,f21212num,&
	su01,ri01,ri12,su12,su02,ri02,fv1denum,fv2denum, &
	fv101num,fv102num,fv112num,fv201num,fv202num,fv212num, & 
	fv10101num,fv20101num, fv10102num,fv20102num, & 
	fv10112num,fv20112num, fv10202num,fv20202num, & 
	fv10212num,fv20212num, fv11212num,fv21212num, &  
	d1mach(5),epmach,uflow,fcdenum,fc01num,fc02num,fc12num, & 
	fc0101num,fc0102num,fc0112num,fc0202num,fc0212num,fc1212num

         double precision gl01,gl12,gl02
         integer::j,jtw,jtwm1
         double precision,dimension(8)::xgk,wgk
	 double precision,dimension(4)::wg
	save wgk,xgk

	 dimension fv1denum(7),fv2denum(7),fv101num(7),fv201num(7), &
	fv102num(7),fv202num(7), fv112num(7),fv212num(7), & 
	fv10101num(7),fv20101num(7),fv10102num(7),fv20102num(7),&
	fv10112num(7),fv20112num(7),fv10202num(7),fv20202num(7),&
	fv10212num(7),fv20212num(7),fv11212num(7),fv21212num(7)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)

	wg(1)=0.129484966168869693270611432679082d0
   	wg(2)=0.279705391489276667901467771423780d0
    	wg(3)=0.381830050505118944950369775488975d0
    	wg(4)=0.417959183673469387755102040816327d0

    	xgk(1)=0.991455371120812639206854697526329d0
    	xgk(2)=0.949107912342758524526189684047851d0
    	xgk(3)=0.864864423359769072789712788640926d0
    	xgk(4)=0.741531185599394439863864773280788d0
    	xgk(5)=0.586087235467691130294144838258730d0
    	xgk(6)=0.405845151377397166906606412076961d0
    	xgk(7)=0.207784955007898467600689403773245d0
    	xgk(8)=0.000000000000000000000000000000000d0

    	wgk(1)=0.022935322010529224963732008058970d0
    	wgk(2)=0.063092092629978553290700663189204d0
    	wgk(3)=0.104790010322250183839876322541518d0
    	wgk(4)=0.140653259715525918745189590510238d0
    	wgk(5)=0.169004726639267902826583426598550d0
    	wgk(6)=0.190350578064785409913256402421014d0
    	wgk(7)=0.204432940075298892414161999234649d0
    	wgk(8)=0.209482141084727828012999174891714d0

	xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)

            resdenum = 0.d0
	    res01num = 0.d0
	    res02num = 0.d0
	    res12num = 0.d0
	
	    res0101num = 0.d0
	    res0102num = 0.d0
	    res0112num = 0.d0
	    res0202num = 0.d0
	    res0212num = 0.d0
	    res1212num = 0.d0
        
    	    	call fonct(xm,the01,ri01,gl01,su01)
    		call fonct(xm,the02,ri02,gl02,su02)
   		call fonct(xm,the12,ri12,gl12,su12)
    		fcdenum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0
    		fc01num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
		fc02num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
		fc12num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

		
		fc0101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
		fc0102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl01*v01/(su12**v12)
		fc0112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl01*v01/(su12**v12)
		fc0202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
		fc0212num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl12*v12/(su12**v12)
		fc1212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

    		reskdenum = fcdenum*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk01num = fc01num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk02num = fc02num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk12num = fc12num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0101num = fc0101num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0102num = fc0102num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0112num = fc0112num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0202num = fc0202num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0212num = fc0212num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk1212num = fc1212num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	




		do j=1,3
	       		jtw = j*2
               		dx=xr*xgk(jtw)
               		xx = xm+dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)
                        f1denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f101num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f10101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f10102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl01*v01/(su12**v12)
			f10112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl01*v01/(su12**v12)
			f10202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f10212num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl12*v12/(su12**v12)
			f11212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

			
               		xx = xm-dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)
                        f2denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f201num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f20101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f20102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl01*v01/(su12**v12)
			f20112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl01*v01/(su12**v12)
			f20202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f20212num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl12*v12/(su12**v12)
			f21212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

			
               		fv1denum(jtw) = f1denum   ! svgrd valeurs fct f a gche du centre
               		fv2denum(jtw) = f2denum   ! svgrd valeurs fct f a drte du centre
			reskdenum = reskdenum + wgk(jtw)*(f1denum+f2denum)
               		

			fv101num(jtw) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtw) = f201num   ! svgrd valeurs fct f a drte du centre
			resk01num = resk01num + wgk(jtw)*(f101num+f201num)


			fv102num(jtw) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtw) = f202num   ! svgrd valeurs fct f a drte du centre
			resk02num = resk02num + wgk(jtw)*(f102num+f202num)
               		

			fv112num(jtw) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtw) = f212num   ! svgrd valeurs fct f a drte du centre
			resk12num = resk12num + wgk(jtw)*(f112num+f212num)

			fv10101num(jtw) = f10101num   ! svgrd valeurs fct f a gche du centre
               		fv20101num(jtw) = f20101num   ! svgrd valeurs fct f a drte du centre
			resk0101num = resk0101num + wgk(jtw)*(f10101num+f20101num)
               		

			fv10102num(jtw) = f10102num   ! svgrd valeurs fct f a gche du centre
               		fv20102num(jtw) = f20102num   ! svgrd valeurs fct f a drte du centre
			resk0102num = resk0102num + wgk(jtw)*(f10102num+f20102num)
               		

			fv10112num(jtw) = f10112num   ! svgrd valeurs fct f a gche du centre
               		fv20112num(jtw) = f20112num   ! svgrd valeurs fct f a drte du centre
               		resk0112num = resk0112num + wgk(jtw)*(f10112num+f20112num)
               		

			fv10202num(jtw) = f10202num   ! svgrd valeurs fct f a gche du centre
               		fv20202num(jtw) = f20202num   ! svgrd valeurs fct f a drte du centre
               		resk0202num = resk0202num + wgk(jtw)*(f10202num+f20202num)
               		
			fv10212num(jtw) = f10212num   ! svgrd valeurs fct f a gche du centre
               		fv20212num(jtw) = f20212num   ! svgrd valeurs fct f a drte du centre
               		resk0212num = resk0212num + wgk(jtw)*(f10212num+f20212num)
               		

			fv11212num(jtw) = f11212num   ! svgrd valeurs fct f a gche du centre
               		fv21212num(jtw) = f21212num   ! svgrd valeurs fct f a drte du centre
               		resk1212num = resk1212num + wgk(jtw)*(f11212num+f21212num)
               		

			
         	end do
	 	do j=1,4
			jtwm1 = j*2-1
               		dx=xr*xgk(jtwm1)
               		xx = xm+dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)
      			f1denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f101num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 


			f10101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f10102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl01*v01/(su12**v12)
			f10112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl01*v01/(su12**v12)
			f10202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f10212num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl12*v12/(su12**v12)
			f11212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)


               		xx = xm-dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)
      			f2denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f201num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f20101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f20102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl01*v01/(su12**v12)
			f20112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl01*v01/(su12**v12)
			f20202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f20212num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl12*v12/(su12**v12)
			f21212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)


               		fv1denum(jtwm1) = f1denum   ! svgrd valeurs fct f a gche du centre
               		fv2denum(jtwm1) = f2denum   ! svgrd valeurs fct f a drte du centre
	       		reskdenum = reskdenum + wgk(jtwm1)*(f1denum+f2denum)
               		

			fv101num(jtwm1) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtwm1) = f201num   ! svgrd valeurs fct f a drte du centre
	       		resk01num = resk01num + wgk(jtwm1)*(f101num+f201num)
               		
			fv102num(jtwm1) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtwm1) = f202num   ! svgrd valeurs fct f a drte du centre
	       		resk02num = resk02num + wgk(jtwm1)*(f102num+f202num)
               		

			fv112num(jtwm1) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtwm1) = f212num   ! svgrd valeurs fct f a drte du centre
	       		resk12num = resk12num + wgk(jtwm1)*(f112num+f212num)
               		

			fv10101num(jtwm1) = f10101num   ! svgrd valeurs fct f a gche du centre
               		fv20101num(jtwm1) = f20101num   ! svgrd valeurs fct f a drte du centre
               		resk0101num = resk0101num + wgk(jtwm1)*(f10101num+f20101num)
               		

			fv10102num(jtwm1) = f10102num   ! svgrd valeurs fct f a gche du centre
               		fv20102num(jtwm1) = f20102num   ! svgrd valeurs fct f a drte du centre
               		resk0102num = resk0102num + wgk(jtwm1)*(f10102num+f20102num)
               		

			fv10112num(jtwm1) = f10112num   ! svgrd valeurs fct f a gche du centre
               		fv20112num(jtwm1) = f20112num   ! svgrd valeurs fct f a drte du centre
               		resk0112num = resk0112num + wgk(jtwm1)*(f10112num+f20112num)
               		
			fv10202num(jtwm1) = f10202num   ! svgrd valeurs fct f a gche du centre
               		fv20202num(jtwm1) = f20202num   ! svgrd valeurs fct f a drte du centre
               		resk0202num = resk0202num + wgk(jtwm1)*(f10202num+f20202num)
               		

			fv10212num(jtwm1) = f10212num   ! svgrd valeurs fct f a gche du centre
               		fv20212num(jtwm1) = f20212num   ! svgrd valeurs fct f a drte du centre
               		resk0212num = resk0212num + wgk(jtwm1)*(f10212num+f20212num)
               		

			fv11212num(jtwm1) = f11212num   ! svgrd valeurs fct f a gche du centre
               		fv21212num(jtwm1) = f21212num   ! svgrd valeurs fct f a drte du centre
               		resk1212num = resk1212num + wgk(jtwm1)*(f11212num+f21212num)
               		



			
         	end do
	
    
    		resdenum = reskdenum*xr
    		res01num = resk01num*xr
    		res02num = resk02num*xr
    		res12num = resk12num*xr
		res0101num = resk0101num*xr
		res0102num = resk0102num*xr
		res0112num = resk0112num*xr
		res0202num = resk0202num*xr
		res0212num = resk0212num*xr
		res1212num = resk1212num*xr
	
         
              
          end subroutine qgaussweibderiv

!=============================================================================================        
!================================  qgaussweibderiv  ==========================================
!=== for derivatives approximation out a 15 point Gauss-Kronrod quadrature rule for weib =====
!=============================================================================================  

			
subroutine qgaussweibbetaderiv(a,b,the01,the02,the12,resdenum,&
resthe01,resthe02,resthe12,resthenum,resthenumsquare,&
resthe0101,resthe0102,resthe0112,&
resthe0202,resthe0212,resthe1212,&
resthe0101square,resthe0102square,resthe0112square,&
resthe0202square,resthe0212square,resthe1212square,&
resthe0101dsquare,resthe0202dsquare,resthe1212dsquare,&
res01num,res02num,res12num,res0101num,res0101numbis, &
res0102num,res0112num,&
res0202num,res0212num,res1212num,v01,v02,v12)

        implicit none
         double precision a,b,the01(2),the02(2),the12(2)
         double precision dx,xm,xr,reskdenum,&
         resdenum,resk01num,res01num,resk02num,res02num, & 
	resk12num,res12num,resk0101num,res0101num,resk0102num,res0102num, &
	resk0112num,res0112num,resk0202num,res0202num, &
	resk0212num,res0212num,resk1212num,res1212num, &
	resk0101numbis,res0101numbis,&
	resthe01,resthe02,resthe12,resthe0101,resthe0102,resthe0112,&
	resthenum,resthenumsquare,resthe0202,resthe0212,resthe1212,&
	reskthe01,reskthe02,reskthe12,reskthe0101,reskthe0102,reskthe0112,&
	reskthe0202,reskthe0212,reskthe1212,&
	reskthenum,reskthenumsquare, resthe0101square,reskthe0101square,&
	resthe0102square,reskthe0102square,resthe0112square,reskthe0112square,&
	resthe0202square,reskthe0202square,resthe0212square,reskthe0212square,&
	resthe1212square,reskthe1212square,resthe0101dsquare,reskthe0101dsquare,&
	resthe0202dsquare,reskthe0202dsquare,&
	resthe1212dsquare,reskthe1212dsquare,&
	v01,v02,v12

         double precision xx,f1denum,f2denum, f101num, f102num, f112num, & 
	f10101num,f10102num,f10112num,f10202num,f10212num,f11212num,&
	f201num, f202num, f212num, f20101num,f20102num,f20112num,&
	f20202num,f20212num,f21212num,&
	su01,ri01,ri12,su12,su02,ri02,fv1denum,fv2denum, &
	fv101num,fv102num,fv112num,fv201num,fv202num,fv212num, & 
	fv10101num,fv20101num, fv10102num,fv20102num, & 
	fv10112num,fv20112num, fv10202num,fv20202num, & 
	fv10212num,fv20212num, fv11212num,fv21212num, &  
	d1mach(5),epmach,uflow,fcdenum,fc01num,fc02num,fc12num
	
	double precision fc0101num,fc0102num,fc0112num,fc0202num, &
	fc0212num,fc1212num, &
	f1the01, f2the01,f1the02, f2the02,f1the12, f2the12, &
	fv1the01, fv2the01,fv1the02, fv2the02,fv1the12, fv2the12, &
	f1thenum,f2thenum,fv1thenum,fv2thenum,&
	f1thenumsquare,f2thenumsquare,fv1thenumsquare,fv2thenumsquare,&
	fcthe01, fcthe02, fcthe12,&
	f1the0101, f2the0101, f1the0102, f2the0102, &
	f1the0112, f2the0112, f1the0202, f2the0202,f1the0212, f2the0212,&
	f1the1212, f2the1212,fv1the0101, fv2the0101,fv1the0102, fv2the0102, &
	fv1the0112, fv2the0112, fv1the0202, fv2the0202, &
	fv1the0212, fv2the0212, fv1the1212, fv2the1212
	
	double precision f1the0101square, f2the0101square, &
	f1the0102square, f2the0102square, &
	f1the0112square, f2the0112square, f1the0202square, f2the0202square,&
	f1the0212square, f2the0212square,&
	f1the1212square, f2the1212square,fv1the0101square, fv2the0101square,&
	fv1the0102square, fv2the0102square, &
	fv1the0112square, fv2the0112square, fv1the0202square,&
	fv2the0202square, fv1the0212square, fv2the0212square, fv1the1212square
	
	double precision fv2the1212square, &
	f1the0101dsquare, f2the0101dsquare, f1the0202dsquare, f2the0202dsquare,&
	f1the1212dsquare, f2the1212dsquare,fv1the0101dsquare, fv2the0101dsquare,&
	fv1the0202dsquare,fv2the0202dsquare, fv1the1212dsquare,&
	fv2the1212dsquare, &
	fcthe0101,fcthe0102,fcthe0112, fcthe0202, fcthe0212, fcthe1212,&
	fcthenum,fcthenumsquare,fcthe0101square,fcthe0102square,&
	fcthe0112square,fcthe0202square,fcthe0212square,&
	fcthe1212square,fcthe0101dsquare,fcthe0202dsquare,fcthe1212dsquare

     double precision gl01,gl12,gl02
     integer::j,jtw,jtwm1
     double precision,dimension(8)::xgk,wgk
	 double precision,dimension(4)::wg
	save wgk,xgk

	 dimension fv1denum(7),fv2denum(7),fv101num(7),fv201num(7), &
	fv102num(7),fv202num(7), fv112num(7),fv212num(7), & 
	fv10101num(7),fv20101num(7),fv10102num(7),fv20102num(7),&
	fv10112num(7),fv20112num(7),fv10202num(7),fv20202num(7),&
	fv10212num(7),fv20212num(7),fv11212num(7),fv21212num(7),&
	fv1the01(7),fv2the01(7),fv1the02(7),fv2the02(7),&
	fv1the12(7),fv2the12(7),fv1thenum(7),fv2thenum(7),&
	fv1thenumsquare(7),fv2thenumsquare(7),&
	fv1the0101(7),fv2the0101(7),&
	fv1the0102(7),fv2the0102(7),fv1the0112(7),fv2the0112(7),&
	fv1the0202(7),fv2the0202(7),fv1the0212(7),fv2the0212(7),&
	fv1the1212(7),fv2the1212(7),fv1the0101square(7),fv2the0101square(7),&
	fv1the0102square(7),fv2the0102square(7),fv1the0112square(7),&
	fv2the0112square(7),fv1the0202square(7),fv2the0202square(7),&
	fv1the0212square(7),fv2the0212square(7),&
	fv1the1212square(7),fv2the1212square(7),fv1the0101dsquare(7),&
	fv2the0101dsquare(7),fv1the0202dsquare(7),&
	fv2the0202dsquare(7),fv1the1212dsquare(7),fv2the1212dsquare(7)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)

		wg(1)=0.129484966168869693270611432679082d0
		wg(2)=0.279705391489276667901467771423780d0
    	wg(3)=0.381830050505118944950369775488975d0
    	wg(4)=0.417959183673469387755102040816327d0

    	xgk(1)=0.991455371120812639206854697526329d0
    	xgk(2)=0.949107912342758524526189684047851d0
    	xgk(3)=0.864864423359769072789712788640926d0
    	xgk(4)=0.741531185599394439863864773280788d0
    	xgk(5)=0.586087235467691130294144838258730d0
    	xgk(6)=0.405845151377397166906606412076961d0
    	xgk(7)=0.207784955007898467600689403773245d0
    	xgk(8)=0.000000000000000000000000000000000d0

    	wgk(1)=0.022935322010529224963732008058970d0
    	wgk(2)=0.063092092629978553290700663189204d0
    	wgk(3)=0.104790010322250183839876322541518d0
    	wgk(4)=0.140653259715525918745189590510238d0
    	wgk(5)=0.169004726639267902826583426598550d0
    	wgk(6)=0.190350578064785409913256402421014d0
    	wgk(7)=0.204432940075298892414161999234649d0
    	wgk(8)=0.209482141084727828012999174891714d0

	xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)

        resdenum = 0.d0
	    res01num = 0.d0
	    res02num = 0.d0
	    res12num = 0.d0
	    resthe01 = 0.d0
	    resthe02 = 0.d0
	    resthe12 = 0.d0
	    resthenum = 0.d0
	    resthenumsquare = 0.d0
	
	    res0101num = 0.d0
		res0101numbis = 0.d0
	    res0102num = 0.d0
	    res0112num = 0.d0
	    res0202num = 0.d0
	    res0212num = 0.d0
	    res1212num = 0.d0

	    resthe0101 = 0.d0
	    resthe0102 = 0.d0
	    resthe0112 = 0.d0
	    resthe0202 = 0.d0
	    resthe0212 = 0.d0
	    resthe1212 = 0.d0

	    resthe0101square = 0.d0
	    resthe0102square = 0.d0
	    resthe0112square = 0.d0
	    resthe0202square = 0.d0
	    resthe0212square = 0.d0
	    resthe1212square = 0.d0

	    resthe0101dsquare = 0.d0
	    resthe0202dsquare = 0.d0
	    resthe1212dsquare = 0.d0
		
        !write(6,*) "the01",the01
		!write(6,*) "the02",the02
		!write(6,*) "the12",the12
		!write(6,*) "xm",xm
		 !write(6,*) "weib",weib

    	call fonct(xm,the01,ri01,gl01,su01)
    	call fonct(xm,the02,ri02,gl02,su02)
   		call fonct(xm,the12,ri12,gl12,su12)
		
		!write(6,*) "su01",su01
		!write(6,*) "su02",su02
		!write(6,*) "su12",su12
		!write(6,*) "ri01",ri01
    	fcdenum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0
    	fc01num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
		fc02num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
		fc12num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

		!write(6,*) "xm",xm
		!write(6,*) "LOG(xm)",LOG(xm)
		fcthe01=(su01**v01)*(su02**v02)*ri01*v01*gl01*v01*LOG(xm)/(su12**v12)
		fcthe02=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*LOG(xm)/(su12**v12)
		fcthe12=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*LOG(xm)/(su12**v12)

		fcthe0101=(su01**v01)*(su02**v02)*ri01*v01*((gl01*v01)**2)*LOG(xm)/(su12**v12)
		fcthe0202=(su01**v01)*(su02**v02)*ri01*v01*((gl02*v02)**2)*LOG(xm)/(su12**v12)
		fcthe1212=(su01**v01)*(su02**v02)*ri01*v01*((gl12*v12)**2)*LOG(xm)/(su12**v12)

		fcthe0101dsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*((gl01*v01)**2)*(LOG(xm)**2)/(su12**v12)
		fcthe0202dsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*((gl02*v02)**2)*(LOG(xm)**2)/(su12**v12)
		fcthe1212dsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*((gl12*v12)**2)*(LOG(xm)**2)/(su12**v12)

		fcthe0101square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(LOG(xm)**2)/(su12**v12)
		fcthe0202square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl02*v02)*(LOG(xm)**2)/(su12**v12)
		fcthe1212square=&
		(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*(LOG(xm)**2)/(su12**v12)

		fcthenum=&
		(su01**v01)*(su02**v02)*ri01*v01*LOG(xm)/(su12**v12)
		fcthenumsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*(LOG(xm)**2)/(su12**v12)

		fcthe0102=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl02*v02)*LOG(xm)/(su12**v12)
		fcthe0112=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl12*v12)*LOG(xm)/(su12**v12)
		fcthe0212=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl02*v02)*gl12*v12*LOG(xm)/(su12**v12)

		fcthe0102square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl02*v02)*(LOG(xm)**2)/(su12**v12)
		fcthe0112square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl12*v12)*(LOG(xm)**2)/(su12**v12)
		fcthe0212square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl02*v02)*gl12*v12*(LOG(xm)**2)/(su12**v12)

		fc0101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
		fc0102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl01*v01/(su12**v12)
		fc0112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl01*v01/(su12**v12)
		fc0202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
		fc0212num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl12*v12/(su12**v12)
		fc1212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)


    		reskdenum = fcdenum*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk01num = fc01num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk02num = fc02num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk12num = fc12num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		
		reskthenum = fcthenum*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthenumsquare = fcthenumsquare*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe01 = fcthe01*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe02 = fcthe02*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe12 = fcthe12*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		
        	resk0101num = fc0101num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0102num = fc0102num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0112num = fc0112num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0202num = fc0202num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0212num = fc0212num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk1212num = fc1212num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	
		reskthe0101 = fcthe0101*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe0202 = fcthe0202*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe1212 = fcthe1212*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe0102 = fcthe0102*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe0112 = fcthe0112*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe0212 = fcthe0212*wgk(8) 
		reskthe0101square = fcthe0101square*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe0202square = fcthe0202square*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe1212square = fcthe1212square*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe0102square = fcthe0102square*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe0112square = fcthe0112square*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe0212square = fcthe0212square*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe0101dsquare = fcthe0101dsquare*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe0202dsquare = fcthe0202dsquare*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		reskthe1212dsquare = fcthe1212dsquare*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
		

		 !write(6,*) "reskdenum",reskdenum


		do j=1,3
	       		jtw = j*2
               		dx=xr*xgk(jtw)
               		xx = xm+dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)
            
			f1denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f101num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f10101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f10102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl01*v01/(su12**v12)
			f10112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl01*v01/(su12**v12)
			f10202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f10212num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl12*v12/(su12**v12)
			f11212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

			!write(6,*) "xx",xx
			!write(6,*) "LOG(xx)",LOG(xx)
			f1the01=(su01**v01)*(su02**v02)*ri01*v01*gl01*v01*LOG(xx)/(su12**v12)
			f1the02=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*LOG(xx)/(su12**v12)
			f1the12=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*LOG(xx)/(su12**v12)

			f1the0101=(su01**v01)*(su02**v02)*ri01*v01*((gl01*v01)**2)*LOG(xx)/(su12**v12)
			f1the0202=(su01**v01)*(su02**v02)*ri01*v01*((gl02*v02)**2)*LOG(xx)/(su12**v12)
			f1the1212=(su01**v01)*(su02**v02)*ri01*v01*((gl12*v12)**2)*LOG(xx)/(su12**v12)

			f1the0101dsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*((gl01*v01)**2)*(LOG(xx)**2)/(su12**v12)
			f1the0202dsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*((gl02*v02)**2)*(LOG(xx)**2)/(su12**v12)
			f1the1212dsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*((gl12*v12)**2)*(LOG(xx)**2)/(su12**v12)

			f1the0101square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(LOG(xx)**2)/(su12**v12)
			f1the0202square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl02*v02)*(LOG(xx)**2)/(su12**v12)
			f1the1212square=&
		(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*(LOG(xx)**2)/(su12**v12)

			f1thenum=&
		(su01**v01)*(su02**v02)*ri01*v01*LOG(xx)/(su12**v12)
		f1thenumsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*(LOG(xx)**2)/(su12**v12)

			f1the0102=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl02*v02)*LOG(xx)/(su12**v12)
			f1the0112=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl12*v12)*LOG(xx)/(su12**v12)
			f1the0212=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl02*v02)*gl12*v12*LOG(xx)/(su12**v12)

			f1the0102square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl02*v02)*(LOG(xx)**2)/(su12**v12)
			f1the0112square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl12*v12)*(LOG(xx)**2)/(su12**v12)
			f1the0212square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl02*v02)*gl12*v12*(LOG(xx)**2)/(su12**v12)


			
               		xx = xm-dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)
                        f2denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f201num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f20101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f20102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl01*v01/(su12**v12)
			f20112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl01*v01/(su12**v12)
			f20202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f20212num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl12*v12/(su12**v12)
			f21212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

			!write(6,*) "xx",xx
			!write(6,*) "LOG(xx)",LOG(xx)
			f2the01=(su01**v01)*(su02**v02)*ri01*v01*gl01*v01*LOG(xx)/(su12**v12)
			f2the02=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*LOG(xx)/(su12**v12)
			f2the12=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*LOG(xx)/(su12**v12)

			f2the0101=(su01**v01)*(su02**v02)*ri01*v01*((gl01*v01)**2)*LOG(xx)/(su12**v12)
			f2the0202=(su01**v01)*(su02**v02)*ri01*v01*((gl02*v02)**2)*LOG(xx)/(su12**v12)
			f2the1212=(su01**v01)*(su02**v02)*ri01*v01*((gl12*v12)**2)*LOG(xx)/(su12**v12)

			f2the0101dsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*((gl01*v01)**2)*(LOG(xx)**2)/(su12**v12)
			f2the0202dsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*((gl02*v02)**2)*(LOG(xx)**2)/(su12**v12)
			f2the1212dsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*((gl12*v12)**2)*(LOG(xx)**2)/(su12**v12)

			f2the0101square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(LOG(xx)**2)/(su12**v12)
			f2the0202square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl02*v02)*(LOG(xx)**2)/(su12**v12)
			f2the1212square=&
		(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*(LOG(xx)**2)/(su12**v12)

			f2thenum=&
		(su01**v01)*(su02**v02)*ri01*v01*LOG(xx)/(su12**v12)
		f2thenumsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*(LOG(xx)**2)/(su12**v12)

			f2the0102=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl02*v02)*LOG(xx)/(su12**v12)
			f2the0112=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl12*v12)*LOG(xx)/(su12**v12)
			f2the0212=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl02*v02)*gl12*v12*LOG(xx)/(su12**v12)

			f2the0102square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl02*v02)*(LOG(xx)**2)/(su12**v12)
			f2the0112square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl12*v12)*(LOG(xx)**2)/(su12**v12)
			f2the0212square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl02*v02)*gl12*v12*(LOG(xx)**2)/(su12**v12)

			
               		fv1denum(jtw) = f1denum   ! svgrd valeurs fct f a gche du centre
               		fv2denum(jtw) = f2denum   ! svgrd valeurs fct f a drte du centre
			reskdenum = reskdenum + wgk(jtw)*(f1denum+f2denum)
			
               		

			fv101num(jtw) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtw) = f201num   ! svgrd valeurs fct f a drte du centre
			resk01num = resk01num + wgk(jtw)*(f101num+f201num)


			fv102num(jtw) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtw) = f202num   ! svgrd valeurs fct f a drte du centre
			resk02num = resk02num + wgk(jtw)*(f102num+f202num)
               		

			fv112num(jtw) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtw) = f212num   ! svgrd valeurs fct f a drte du centre
			resk12num = resk12num + wgk(jtw)*(f112num+f212num)

			fv10101num(jtw) = f10101num   ! svgrd valeurs fct f a gche du centre
               		fv20101num(jtw) = f20101num   ! svgrd valeurs fct f a drte du centre
			resk0101num = resk0101num + wgk(jtw)*(f10101num+f20101num)
               		

			fv10102num(jtw) = f10102num   ! svgrd valeurs fct f a gche du centre
               		fv20102num(jtw) = f20102num   ! svgrd valeurs fct f a drte du centre
			resk0102num = resk0102num + wgk(jtw)*(f10102num+f20102num)
               		

			fv10112num(jtw) = f10112num   ! svgrd valeurs fct f a gche du centre
               		fv20112num(jtw) = f20112num   ! svgrd valeurs fct f a drte du centre
               		resk0112num = resk0112num + wgk(jtw)*(f10112num+f20112num)
               		

			fv10202num(jtw) = f10202num   ! svgrd valeurs fct f a gche du centre
               		fv20202num(jtw) = f20202num   ! svgrd valeurs fct f a drte du centre
               		resk0202num = resk0202num + wgk(jtw)*(f10202num+f20202num)
               		
			fv10212num(jtw) = f10212num   ! svgrd valeurs fct f a gche du centre
               		fv20212num(jtw) = f20212num   ! svgrd valeurs fct f a drte du centre
               		resk0212num = resk0212num + wgk(jtw)*(f10212num+f20212num)
               		

			fv11212num(jtw) = f11212num   ! svgrd valeurs fct f a gche du centre
               		fv21212num(jtw) = f21212num   ! svgrd valeurs fct f a drte du centre
               		resk1212num = resk1212num + wgk(jtw)*(f11212num+f21212num)
					
			fv1thenum(jtw) = f1thenum   ! svgrd valeurs fct f a gche du centre
               		fv2thenum(jtw) = f2thenum   ! svgrd valeurs fct f a drte du centre
			reskthenum = reskthenum + wgk(jtw)*(f1thenum+f2thenum)
               		
				fv1thenumsquare(jtw) = f1thenumsquare   ! svgrd valeurs fct f a gche du centre
               		fv2thenumsquare(jtw) = f2thenumsquare   ! svgrd valeurs fct f a drte du centre
			reskthenumsquare = reskthenumsquare + wgk(jtw)*(f1thenumsquare+f2thenumsquare)
			
				fv1the01(jtw) = f1the01   ! svgrd valeurs fct f a gche du centre
               		fv2the01(jtw) = f2the01  ! svgrd valeurs fct f a drte du centre
			reskthe01 = reskthe01 + wgk(jtw)*(f1the01+f2the01)
			
			fv1the02(jtw) = f1the02   ! svgrd valeurs fct f a gche du centre
               		fv2the02(jtw) = f2the02  ! svgrd valeurs fct f a drte du centre
			reskthe02 = reskthe02 + wgk(jtw)*(f1the02+f2the02)
			
			fv1the12(jtw) = f1the12   ! svgrd valeurs fct f a gche du centre
               		fv2the12(jtw) = f2the12  ! svgrd valeurs fct f a drte du centre
			reskthe12 = reskthe12 + wgk(jtw)*(f1the12+f2the12)
			
			fv1the0101(jtw) = f1the0101   ! svgrd valeurs fct f a gche du centre
               		fv2the0101(jtw) = f2the0101  ! svgrd valeurs fct f a drte du centre
			reskthe0101 = reskthe0101 + wgk(jtw)*(f1the0101+f2the0101)
			
			fv1the0202(jtw) = f1the0202   ! svgrd valeurs fct f a gche du centre
               		fv2the0202(jtw) = f2the0202  ! svgrd valeurs fct f a drte du centre
			reskthe0202 = reskthe0202 + wgk(jtw)*(f1the0202+f2the0202)
			
			fv1the1212(jtw) = f1the1212   ! svgrd valeurs fct f a gche du centre
               		fv2the1212(jtw) = f2the1212  ! svgrd valeurs fct f a drte du centre
			reskthe1212 = reskthe1212 + wgk(jtw)*(f1the1212+f2the1212)
			
			fv1the0102(jtw) = f1the0102   ! svgrd valeurs fct f a gche du centre
               		fv2the0102(jtw) = f2the0102  ! svgrd valeurs fct f a drte du centre
			reskthe0102 = reskthe0102 + wgk(jtw)*(f1the0102+f2the0102)
			
			fv1the0112(jtw) = f1the0112   ! svgrd valeurs fct f a gche du centre
               		fv2the0112(jtw) = f2the0112  ! svgrd valeurs fct f a drte du centre
			reskthe0112 = reskthe0112 + wgk(jtw)*(f1the0112+f2the0112)
			
			fv1the0212(jtw) = f1the0212   ! svgrd valeurs fct f a gche du centre
               		fv2the0212(jtw) = f2the0212  ! svgrd valeurs fct f a drte du centre
			reskthe0212 = reskthe0212 + wgk(jtw)*(f1the0212+f2the0212)
			
			
			fv1the0101square(jtw) = f1the0101square   ! svgrd valeurs fct f a gche du centre
               		fv2the0101square(jtw) = f2the0101square  ! svgrd valeurs fct f a drte du centre
			reskthe0101square = reskthe0101square + wgk(jtw)*(f1the0101square+f2the0101square)
			
			fv1the0202square(jtw) = f1the0202square   ! svgrd valeurs fct f a gche du centre
               		fv2the0202square(jtw) = f2the0202square  ! svgrd valeurs fct f a drte du centre
			reskthe0202square = reskthe0202square + wgk(jtw)*(f1the0202square+f2the0202square)
			
			fv1the1212square(jtw) = f1the1212square   ! svgrd valeurs fct f a gche du centre
               		fv2the1212square(jtw) = f2the1212square  ! svgrd valeurs fct f a drte du centre
			reskthe1212square = reskthe1212square + wgk(jtw)*(f1the1212square+f2the1212square)
			
			fv1the0102square(jtw) = f1the0102square   ! svgrd valeurs fct f a gche du centre
               		fv2the0102square(jtw) = f2the0102square  ! svgrd valeurs fct f a drte du centre
			reskthe0102square = reskthe0102square + wgk(jtw)*(f1the0102square+f2the0102square)
			
			fv1the0112square(jtw) = f1the0112square   ! svgrd valeurs fct f a gche du centre
               		fv2the0112square(jtw) = f2the0112square  ! svgrd valeurs fct f a drte du centre
			reskthe0112square = reskthe0112square + wgk(jtw)*(f1the0112square+f2the0112square)
			
			fv1the0212square(jtw) = f1the0212square   ! svgrd valeurs fct f a gche du centre
               		fv2the0212square(jtw) = f2the0212square  ! svgrd valeurs fct f a drte du centre
			reskthe0212square = reskthe0212square + wgk(jtw)*(f1the0212square+f2the0212square)
             
			fv1the0101dsquare(jtw) = f1the0101dsquare   ! svgrd valeurs fct f a gche du centre
               		fv2the0101dsquare(jtw) = f2the0101dsquare  ! svgrd valeurs fct f a drte du centre
			reskthe0101dsquare = reskthe0101dsquare + wgk(jtw)*(f1the0101dsquare+f2the0101dsquare)
			
			fv1the0202dsquare(jtw) = f1the0202dsquare   ! svgrd valeurs fct f a gche du centre
               		fv2the0202dsquare(jtw) = f2the0202dsquare  ! svgrd valeurs fct f a drte du centre
			reskthe0202dsquare = reskthe0202dsquare + wgk(jtw)*(f1the0202dsquare+f2the0202dsquare)
			
			fv1the1212dsquare(jtw) = f1the1212dsquare   ! svgrd valeurs fct f a gche du centre
               		fv2the1212dsquare(jtw) = f2the1212dsquare  ! svgrd valeurs fct f a drte du centre
			reskthe1212dsquare = reskthe1212dsquare + wgk(jtw)*(f1the1212dsquare+f2the1212dsquare)
			
			 !write(6,*) "reskdenum",reskdenum

         	end do
	 	do j=1,4
			jtwm1 = j*2-1
               		dx=xr*xgk(jtwm1)
               		xx = xm+dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)
      			f1denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f101num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 


			f10101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f10102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl01*v01/(su12**v12)
			f10112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl01*v01/(su12**v12)
			f10202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f10212num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl12*v12/(su12**v12)
			f11212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

			!write(6,*) "xx",xx
			!write(6,*) "LOG(xx)",LOG(xx)
			f1the01=(su01**v01)*(su02**v02)*ri01*v01*gl01*v01*LOG(xx)/(su12**v12)
			f1the02=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*LOG(xx)/(su12**v12)
			f1the12=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*LOG(xx)/(su12**v12)

			f1the0101=(su01**v01)*(su02**v02)*ri01*v01*((gl01*v01)**2)*LOG(xx)/(su12**v12)
			f1the0202=(su01**v01)*(su02**v02)*ri01*v01*((gl02*v02)**2)*LOG(xx)/(su12**v12)
			f1the1212=(su01**v01)*(su02**v02)*ri01*v01*((gl12*v12)**2)*LOG(xx)/(su12**v12)

			f1the0101dsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*((gl01*v01)**2)*(LOG(xx)**2)/(su12**v12)
			f1the0202dsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*((gl02*v02)**2)*(LOG(xx)**2)/(su12**v12)
			f1the1212dsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*((gl12*v12)**2)*(LOG(xx)**2)/(su12**v12)

			f1the0101square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(LOG(xx)**2)/(su12**v12)
			f1the0202square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl02*v02)*(LOG(xx)**2)/(su12**v12)
			f1the1212square=&
		(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*(LOG(xx)**2)/(su12**v12)

			f1thenum=&
		(su01**v01)*(su02**v02)*ri01*v01*LOG(xx)/(su12**v12)
		f1thenumsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*(LOG(xx)**2)/(su12**v12)

			f1the0102=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl02*v02)*LOG(xx)/(su12**v12)
			f1the0112=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl12*v12)*LOG(xx)/(su12**v12)
			f1the0212=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl02*v02)*gl12*v12*LOG(xx)/(su12**v12)

			f1the0102square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl02*v02)*(LOG(xx)**2)/(su12**v12)
			f1the0112square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl12*v12)*(LOG(xx)**2)/(su12**v12)
			f1the0212square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl02*v02)*gl12*v12*(LOG(xx)**2)/(su12**v12)

               		xx = xm-dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)
      			f2denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f201num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f20101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f20102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl01*v01/(su12**v12)
			f20112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl01*v01/(su12**v12)
			f20202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f20212num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl12*v12/(su12**v12)
			f21212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

			!write(6,*) "xx",xx
			!write(6,*) "LOG(xx)",LOG(xx)
			f2the01=(su01**v01)*(su02**v02)*ri01*v01*gl01*v01*LOG(xx)/(su12**v12)
			f2the02=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*LOG(xx)/(su12**v12)
			f2the12=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*LOG(xx)/(su12**v12)

			f2the0101=(su01**v01)*(su02**v02)*ri01*v01*((gl01*v01)**2)*LOG(xx)/(su12**v12)
			f2the0202=(su01**v01)*(su02**v02)*ri01*v01*((gl02*v02)**2)*LOG(xx)/(su12**v12)
			f2the1212=(su01**v01)*(su02**v02)*ri01*v01*((gl12*v12)**2)*LOG(xx)/(su12**v12)

			f2the0101dsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*((gl01*v01)**2)*(LOG(xx)**2)/(su12**v12)
			f2the0202dsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*((gl02*v02)**2)*(LOG(xx)**2)/(su12**v12)
			f2the1212dsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*((gl12*v12)**2)*(LOG(xx)**2)/(su12**v12)

			f2the0101square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(LOG(xx)**2)/(su12**v12)
			f2the0202square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl02*v02)*(LOG(xx)**2)/(su12**v12)
			f2the1212square=&
		(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*(LOG(xx)**2)/(su12**v12)

			f2thenum=&
		(su01**v01)*(su02**v02)*ri01*v01*LOG(xx)/(su12**v12)
		f2thenumsquare=&
		(su01**v01)*(su02**v02)*ri01*v01*(LOG(xx)**2)/(su12**v12)

			f2the0102=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl02*v02)*LOG(xx)/(su12**v12)
			f2the0112=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl12*v12)*LOG(xx)/(su12**v12)
			f2the0212=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl02*v02)*gl12*v12*LOG(xx)/(su12**v12)

			f2the0102square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl02*v02)*(LOG(xx)**2)/(su12**v12)
			f2the0112square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl01*v01)*(gl12*v12)*(LOG(xx)**2)/(su12**v12)
			f2the0212square=&
		(su01**v01)*(su02**v02)*ri01*v01*(gl02*v02)*gl12*v12*(LOG(xx)**2)/(su12**v12)

               		fv1denum(jtwm1) = f1denum   ! svgrd valeurs fct f a gche du centre
               		fv2denum(jtwm1) = f2denum   ! svgrd valeurs fct f a drte du centre
	       		reskdenum = reskdenum + wgk(jtwm1)*(f1denum+f2denum)
               		

			fv101num(jtwm1) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtwm1) = f201num   ! svgrd valeurs fct f a drte du centre
	       		resk01num = resk01num + wgk(jtwm1)*(f101num+f201num)
               		
			fv102num(jtwm1) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtwm1) = f202num   ! svgrd valeurs fct f a drte du centre
	       		resk02num = resk02num + wgk(jtwm1)*(f102num+f202num)
               		

			fv112num(jtwm1) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtwm1) = f212num   ! svgrd valeurs fct f a drte du centre
	       		resk12num = resk12num + wgk(jtwm1)*(f112num+f212num)
               		

			fv10101num(jtwm1) = f10101num   ! svgrd valeurs fct f a gche du centre
               		fv20101num(jtwm1) = f20101num   ! svgrd valeurs fct f a drte du centre
               		resk0101num = resk0101num + wgk(jtwm1)*(f10101num+f20101num)
               		

			fv10102num(jtwm1) = f10102num   ! svgrd valeurs fct f a gche du centre
               		fv20102num(jtwm1) = f20102num   ! svgrd valeurs fct f a drte du centre
               		resk0102num = resk0102num + wgk(jtwm1)*(f10102num+f20102num)
               		

			fv10112num(jtwm1) = f10112num   ! svgrd valeurs fct f a gche du centre
               		fv20112num(jtwm1) = f20112num   ! svgrd valeurs fct f a drte du centre
               		resk0112num = resk0112num + wgk(jtwm1)*(f10112num+f20112num)
               		
			fv10202num(jtwm1) = f10202num   ! svgrd valeurs fct f a gche du centre
               		fv20202num(jtwm1) = f20202num   ! svgrd valeurs fct f a drte du centre
               		resk0202num = resk0202num + wgk(jtwm1)*(f10202num+f20202num)
               		

			fv10212num(jtwm1) = f10212num   ! svgrd valeurs fct f a gche du centre
               		fv20212num(jtwm1) = f20212num   ! svgrd valeurs fct f a drte du centre
               		resk0212num = resk0212num + wgk(jtwm1)*(f10212num+f20212num)
               		

			fv11212num(jtwm1) = f11212num   ! svgrd valeurs fct f a gche du centre
               		fv21212num(jtwm1) = f21212num   ! svgrd valeurs fct f a drte du centre
               		resk1212num = resk1212num + wgk(jtwm1)*(f11212num+f21212num)
               		


			fv1thenum(jtwm1) = f1thenum   ! svgrd valeurs fct f a gche du centre
               		fv2thenum(jtwm1) = f2thenum   ! svgrd valeurs fct f a drte du centre
			reskthenum = reskthenum + wgk(jtwm1)*(f1thenum+f2thenum)
               		
				fv1thenumsquare(jtwm1) = f1thenumsquare   ! svgrd valeurs fct f a gche du centre
               		fv2thenumsquare(jtwm1) = f2thenumsquare   ! svgrd valeurs fct f a drte du centre
			reskthenumsquare = reskthenumsquare + wgk(jtwm1)*(f1thenumsquare+f2thenumsquare)
			
				fv1the01(jtwm1) = f1the01   ! svgrd valeurs fct f a gche du centre
               		fv2the01(jtwm1) = f2the01  ! svgrd valeurs fct f a drte du centre
			reskthe01 = reskthe01 + wgk(jtwm1)*(f1the01+f2the01)
			
			fv1the02(jtwm1) = f1the02   ! svgrd valeurs fct f a gche du centre
               		fv2the02(jtwm1) = f2the02  ! svgrd valeurs fct f a drte du centre
			reskthe02 = reskthe02 + wgk(jtwm1)*(f1the02+f2the02)
			
			fv1the12(jtwm1) = f1the12   ! svgrd valeurs fct f a gche du centre
               		fv2the12(jtwm1) = f2the12  ! svgrd valeurs fct f a drte du centre
			reskthe12 = reskthe12 + wgk(jtwm1)*(f1the12+f2the12)
			
			fv1the0101(jtwm1) = f1the0101   ! svgrd valeurs fct f a gche du centre
               		fv2the0101(jtwm1) = f2the0101  ! svgrd valeurs fct f a drte du centre
			reskthe0101 = reskthe0101 + wgk(jtwm1)*(f1the0101+f2the0101)
			
			fv1the0202(jtwm1) = f1the0202   ! svgrd valeurs fct f a gche du centre
               		fv2the0202(jtwm1) = f2the0202  ! svgrd valeurs fct f a drte du centre
			reskthe0202 = reskthe0202 + wgk(jtwm1)*(f1the0202+f2the0202)
			
			fv1the1212(jtwm1) = f1the1212   ! svgrd valeurs fct f a gche du centre
               		fv2the1212(jtwm1) = f2the1212  ! svgrd valeurs fct f a drte du centre
			reskthe1212 = reskthe1212 + wgk(jtwm1)*(f1the1212+f2the1212)
			
			fv1the0102(jtwm1) = f1the0102   ! svgrd valeurs fct f a gche du centre
               		fv2the0102(jtwm1) = f2the0102  ! svgrd valeurs fct f a drte du centre
			reskthe0102 = reskthe0102 + wgk(jtwm1)*(f1the0102+f2the0102)
			
			fv1the0112(jtwm1) = f1the0112   ! svgrd valeurs fct f a gche du centre
               		fv2the0112(jtwm1) = f2the0112  ! svgrd valeurs fct f a drte du centre
			reskthe0112 = reskthe0112 + wgk(jtwm1)*(f1the0112+f2the0112)
			
			fv1the0212(jtwm1) = f1the0212   ! svgrd valeurs fct f a gche du centre
               		fv2the0212(jtwm1) = f2the0212  ! svgrd valeurs fct f a drte du centre
			reskthe0212 = reskthe0212 + wgk(jtwm1)*(f1the0212+f2the0212)
			
			
			fv1the0101square(jtwm1) = f1the0101square   ! svgrd valeurs fct f a gche du centre
               		fv2the0101square(jtwm1) = f2the0101square  ! svgrd valeurs fct f a drte du centre
			reskthe0101square = reskthe0101square + wgk(jtwm1)*(f1the0101square+f2the0101square)
			
			fv1the0202square(jtwm1) = f1the0202square   ! svgrd valeurs fct f a gche du centre
               		fv2the0202square(jtwm1) = f2the0202square  ! svgrd valeurs fct f a drte du centre
			reskthe0202square = reskthe0202square + wgk(jtwm1)*(f1the0202square+f2the0202square)
			
			fv1the1212square(jtwm1) = f1the1212square   ! svgrd valeurs fct f a gche du centre
               		fv2the1212square(jtwm1) = f2the1212square  ! svgrd valeurs fct f a drte du centre
			reskthe1212square = reskthe1212square + wgk(jtwm1)*(f1the1212square+f2the1212square)
			
			fv1the0102square(jtwm1) = f1the0102square   ! svgrd valeurs fct f a gche du centre
               		fv2the0102square(jtwm1) = f2the0102square  ! svgrd valeurs fct f a drte du centre
			reskthe0102square = reskthe0102square + wgk(jtwm1)*(f1the0102square+f2the0102square)
			
			fv1the0112square(jtwm1) = f1the0112square   ! svgrd valeurs fct f a gche du centre
               		fv2the0112square(jtwm1) = f2the0112square  ! svgrd valeurs fct f a drte du centre
			reskthe0112square = reskthe0112square + wgk(jtwm1)*(f1the0112square+f2the0112square)
			
			fv1the0212square(jtwm1) = f1the0212square   ! svgrd valeurs fct f a gche du centre
               		fv2the0212square(jtwm1) = f2the0212square  ! svgrd valeurs fct f a drte du centre
			reskthe0212square = reskthe0212square + wgk(jtwm1)*(f1the0212square+f2the0212square)
             
			fv1the0101dsquare(jtwm1) = f1the0101dsquare   ! svgrd valeurs fct f a gche du centre
               		fv2the0101dsquare(jtwm1) = f2the0101dsquare  ! svgrd valeurs fct f a drte du centre
			reskthe0101dsquare = reskthe0101dsquare + wgk(jtwm1)*(f1the0101dsquare+f2the0101dsquare)
			
			fv1the0202dsquare(jtwm1) = f1the0202dsquare   ! svgrd valeurs fct f a gche du centre
               		fv2the0202dsquare(jtwm1) = f2the0202dsquare  ! svgrd valeurs fct f a drte du centre
			reskthe0202dsquare = reskthe0202dsquare + wgk(jtwm1)*(f1the0202dsquare+f2the0202dsquare)
			
			fv1the1212dsquare(jtwm1) = f1the1212dsquare   ! svgrd valeurs fct f a gche du centre
               		fv2the1212dsquare(jtwm1) = f2the1212dsquare  ! svgrd valeurs fct f a drte du centre
			reskthe1212dsquare = reskthe1212dsquare + wgk(jtwm1)*(f1the1212dsquare+f2the1212dsquare)
			
			 !write(6,*) "reskdenum",reskdenum

         	end do
	
    !write(6,*) "end  do qgauss"
	!write(6,*) "reskdenum",reskdenum

    		resdenum = reskdenum*xr
    		res01num = resk01num*xr
    		res02num = resk02num*xr
    		res12num = resk12num*xr
		res0101num = resk0101num*xr
		res0102num = resk0102num*xr
		res0112num = resk0112num*xr
		res0202num = resk0202num*xr
		res0212num = resk0212num*xr
		res1212num = resk1212num*xr
		res0101numbis = res0101num+2*resdenum-3*res01num
		
		resthenum = reskthenum*xr
		resthenumsquare = reskthenumsquare*xr
		resthe01 = reskthe01*xr
		resthe02 = reskthe02*xr
		resthe12 = reskthe12*xr
		resthe0101 = reskthe0101*xr
		resthe0202 = reskthe0202*xr
		resthe1212 = reskthe1212*xr
		resthe0102 = reskthe0102*xr
		resthe0112 = reskthe0112*xr
		resthe0212 = reskthe0212*xr
		resthe0101square = reskthe0101square*xr
		resthe0202square = reskthe0202square*xr
		resthe1212square = reskthe1212square*xr
		resthe0102square = reskthe0102square*xr
		resthe0112square = reskthe0112square*xr
		resthe0212square = reskthe0212square*xr
		resthe0101dsquare = reskthe0101dsquare*xr
		resthe0202dsquare = reskthe0202dsquare*xr
		resthe1212dsquare = reskthe1212dsquare*xr
		!write(6,*) "end qgauss calculation"
         
              
          end subroutine qgaussweibbetaderiv
!=============================================================================================        
!================================  qgaussweibderiv  ==========================================
!=== for derivatives approximation out a 15 point Gauss-Kronrod quadrature rule for weib =====
!========================= looking only at diag hessian ======================================
!=============================================================================================  


subroutine qgaussweibfirstderiv(a,b,the01,the02,the12,resdenum,&
res01num,res02num,res12num,v01,v02,v12)

        implicit none
         double precision a,b,the01(2),the02(2),the12(2)
         double precision dx,xm,xr,reskdenum,&
         reskhdenum,resdenum,resk01num,res01num, & 
	resk02num,res02num,resk12num,res12num,v01,v02,v12
         double precision xx,f1denum,f2denum, f101num, f102num, f112num, & 
	f201num, f202num, f212num, &
	su01,ri01,ri12,su12,su02,ri02,fv1denum,fv2denum, &
	fv101num,fv102num,fv112num,fv201num,fv202num,fv212num, &  
	d1mach(5),epmach,uflow,fcdenum,fc01num,fc02num,fc12num

         double precision gl01,gl12,gl02
         integer::j,jtw,jtwm1
         double precision,dimension(8)::xgk,wgk
	 double precision,dimension(4)::wg
	save wgk,xgk

	 dimension fv1denum(7),fv2denum(7),fv101num(7),fv201num(7), &
	fv102num(7),fv202num(7), fv112num(7),fv212num(7)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)

	wg(1)=0.129484966168869693270611432679082d0
   	wg(2)=0.279705391489276667901467771423780d0
    	wg(3)=0.381830050505118944950369775488975d0
    	wg(4)=0.417959183673469387755102040816327d0

    	xgk(1)=0.991455371120812639206854697526329d0
    	xgk(2)=0.949107912342758524526189684047851d0
    	xgk(3)=0.864864423359769072789712788640926d0
    	xgk(4)=0.741531185599394439863864773280788d0
    	xgk(5)=0.586087235467691130294144838258730d0
    	xgk(6)=0.405845151377397166906606412076961d0
    	xgk(7)=0.207784955007898467600689403773245d0
    	xgk(8)=0.000000000000000000000000000000000d0

    	wgk(1)=0.022935322010529224963732008058970d0
    	wgk(2)=0.063092092629978553290700663189204d0
    	wgk(3)=0.104790010322250183839876322541518d0
    	wgk(4)=0.140653259715525918745189590510238d0
    	wgk(5)=0.169004726639267902826583426598550d0
    	wgk(6)=0.190350578064785409913256402421014d0
    	wgk(7)=0.204432940075298892414161999234649d0
    	wgk(8)=0.209482141084727828012999174891714d0

	xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)

            resdenum = 0.d0
	    res01num = 0.d0
	    res02num = 0.d0
	    res12num = 0.d0
	

        
    	    	call fonct(xm,the01,ri01,gl01,su01)
    		call fonct(xm,the02,ri02,gl02,su02)
   		call fonct(xm,the12,ri12,gl12,su12)

    		fcdenum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0
    		fc01num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
		fc02num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
		fc12num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 
		
        	reskdenum = fcdenum*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk01num = fc01num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk02num = fc02num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk12num = fc12num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	
		do j=1,3
	       		jtw = j*2
               		dx=xr*xgk(jtw)
               		xx = xm+dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)
                        
			f1denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f101num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 
			
               		xx = xm-dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)
                        f2denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f201num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

               		fv1denum(jtw) = f1denum   ! svgrd valeurs fct f a gche du centre
               		fv2denum(jtw) = f2denum   ! svgrd valeurs fct f a drte du centre

               		reskdenum = reskdenum + wgk(jtw)*(f1denum+f2denum)
               		
			fv101num(jtw) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtw) = f201num   ! svgrd valeurs fct f a drte du centre

               		resk01num = resk01num + wgk(jtw)*(f101num+f201num)
               		
			fv102num(jtw) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtw) = f202num   ! svgrd valeurs fct f a drte du centre

               		resk02num = resk02num + wgk(jtw)*(f102num+f202num)
               		
			fv112num(jtw) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtw) = f212num   ! svgrd valeurs fct f a drte du centre

               		resk12num = resk12num + wgk(jtw)*(f112num+f212num)
               		

			
         	end do
	 	do j=1,4
			jtwm1 = j*2-1
               		dx=xr*xgk(jtwm1)
               		xx = xm+dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)
      			
			f1denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f101num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 
			
               		xx = xm-dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)

      			f2denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f201num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 
			

               		fv1denum(jtwm1) = f1denum   ! svgrd valeurs fct f a gche du centre
               		fv2denum(jtwm1) = f2denum   ! svgrd valeurs fct f a drte du centre
	       		reskdenum = reskdenum + wgk(jtwm1)*(f1denum+f2denum)
               		
			fv101num(jtwm1) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtwm1) = f201num   ! svgrd valeurs fct f a drte du centre
	       		resk01num = resk01num + wgk(jtwm1)*(f101num+f201num)
               		

			fv102num(jtwm1) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtwm1) = f202num   ! svgrd valeurs fct f a drte du centre
	       		resk02num = resk02num + wgk(jtwm1)*(f102num+f202num)
               		

			fv112num(jtwm1) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtwm1) = f212num   ! svgrd valeurs fct f a drte du centre
	       		resk12num = resk12num + wgk(jtwm1)*(f112num+f212num)
               		

         	end do

    		resdenum = xr*reskdenum
    		res01num = xr*resk01num 
    		res02num = xr*resk02num
    		res12num = xr*resk12num 
	
         
              
          end subroutine qgaussweibfirstderiv
!=============================================================================================        
!================================  qgaussweibderiv  ==========================================
!=== for derivatives approximation out a 15 point Gauss-Kronrod quadrature rule for weib =====
!========================= looking only at diag hessian ======================================
!=============================================================================================  


subroutine qgaussweibderivdiag(a,b,the01,the02,the12,resdenum,&
res01num,res02num,res12num,res0101num,&
res0202num,res1212num,v01,v02,v12)

        implicit none
         double precision a,b,the01(2),the02(2),the12(2)
         double precision dx,xm,xr,reskdenum,resdenum,&
	resk01num,res01num,resk02num,res02num,resk12num,res12num, & 
	resk0101num,res0101num,resk0202num,res0202num, &
	resk1212num,res1212num,v01,v02,v12
         double precision xx,f1denum,f2denum, f101num, f102num, f112num, & 
	f10101num,f10202num,f11212num,&
	f201num, f202num, f212num, f20101num,f20202num,f21212num,&
	su01,ri01,ri12,su12,su02,ri02,fv1denum,fv2denum, &
	fv101num,fv102num,fv112num,fv201num,fv202num,fv212num, & 
	fv10101num,fv20101num, fv10202num,fv20202num, & 
	fv11212num,fv21212num, &  
	d1mach(5),epmach,uflow,fcdenum,fc01num,fc02num,fc12num, & 
	fc0101num,fc0202num,fc1212num

         double precision gl01,gl12,gl02
         integer::j,jtw,jtwm1
         double precision,dimension(8)::xgk,wgk
	 double precision,dimension(4)::wg
	save wgk,xgk

	 dimension fv1denum(7),fv2denum(7),fv101num(7),fv201num(7), &
	fv102num(7),fv202num(7), fv112num(7),fv212num(7), & 
	fv10101num(7),fv20101num(7),fv10202num(7),fv20202num(7),&
	fv11212num(7),fv21212num(7)

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)

	wg(1)=0.129484966168869693270611432679082d0
   	wg(2)=0.279705391489276667901467771423780d0
    	wg(3)=0.381830050505118944950369775488975d0
    	wg(4)=0.417959183673469387755102040816327d0

    	xgk(1)=0.991455371120812639206854697526329d0
    	xgk(2)=0.949107912342758524526189684047851d0
    	xgk(3)=0.864864423359769072789712788640926d0
    	xgk(4)=0.741531185599394439863864773280788d0
    	xgk(5)=0.586087235467691130294144838258730d0
    	xgk(6)=0.405845151377397166906606412076961d0
    	xgk(7)=0.207784955007898467600689403773245d0
    	xgk(8)=0.000000000000000000000000000000000d0

    	wgk(1)=0.022935322010529224963732008058970d0
    	wgk(2)=0.063092092629978553290700663189204d0
    	wgk(3)=0.104790010322250183839876322541518d0
    	wgk(4)=0.140653259715525918745189590510238d0
    	wgk(5)=0.169004726639267902826583426598550d0
    	wgk(6)=0.190350578064785409913256402421014d0
    	wgk(7)=0.204432940075298892414161999234649d0
    	wgk(8)=0.209482141084727828012999174891714d0

	xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)

            resdenum = 0.d0
	    res01num = 0.d0
	    res02num = 0.d0
	    res12num = 0.d0
	
	    res0101num = 0.d0
	    res0202num = 0.d0
	    res1212num = 0.d0
        
    	    	call fonct(xm,the01,ri01,gl01,su01)
    		call fonct(xm,the02,ri02,gl02,su02)
   		call fonct(xm,the12,ri12,gl12,su12)

    		fcdenum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0
    		fc01num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
		fc02num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
		fc12num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 
		fc0101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
		fc0202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
		fc1212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

        	reskdenum = fcdenum*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk01num = fc01num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk02num = fc02num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk12num = fc12num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0101num = fc0101num*wgk(8)      ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0202num = fc0202num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk1212num = fc1212num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	


		do j=1,3
	       		jtw = j*2
               		dx=xr*xgk(jtw)
               		xx = xm+dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)
                        
			f1denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f101num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f10101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f10202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f11212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

			
               		xx = xm-dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)
                        f2denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f201num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f20101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f20202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f21212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

			
               		fv1denum(jtw) = f1denum   ! svgrd valeurs fct f a gche du centre
               		fv2denum(jtw) = f2denum   ! svgrd valeurs fct f a drte du centre

               		reskdenum = reskdenum + wgk(jtw)*(f1denum+f2denum)
               		
			fv101num(jtw) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtw) = f201num   ! svgrd valeurs fct f a drte du centre

               		resk01num = resk01num + wgk(jtw)*(f101num+f201num)
               		
			fv102num(jtw) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtw) = f202num   ! svgrd valeurs fct f a drte du centre

               		resk02num = resk02num + wgk(jtw)*(f102num+f202num)
               		
			fv112num(jtw) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtw) = f212num   ! svgrd valeurs fct f a drte du centre

               		resk12num = resk12num + wgk(jtw)*(f112num+f212num)
               		
			fv10101num(jtw) = f10101num   ! svgrd valeurs fct f a gche du centre
               		fv20101num(jtw) = f20101num   ! svgrd valeurs fct f a drte du centre

               		resk0101num = resk0101num + wgk(jtw)*(f10101num+f20101num)
               		
			fv10202num(jtw) = f10202num   ! svgrd valeurs fct f a gche du centre
               		fv20202num(jtw) = f20202num   ! svgrd valeurs fct f a drte du centre

               		resk0202num = resk0202num + wgk(jtw)*(f10202num+f20202num)
               		
			fv11212num(jtw) = f11212num   ! svgrd valeurs fct f a gche du centre
               		fv21212num(jtw) = f21212num   ! svgrd valeurs fct f a drte du centre

               		resk1212num = resk1212num + wgk(jtw)*(f11212num+f21212num)
               		

			
         	end do
	 	do j=1,4
			jtwm1 = j*2-1
               		dx=xr*xgk(jtwm1)
               		xx = xm+dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)
      			
			f1denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f101num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 
			f10101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f10202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f11212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)


               		xx = xm-dx
               		call fonct(xx,the01,ri01,gl01,su01)
               		call fonct(xx,the02,ri02,gl02,su02)
	       		call fonct(xx,the12,ri12,gl12,su12)

      			f2denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f201num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
			f202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 
			f20101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f20202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f21212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)


               		fv1denum(jtwm1) = f1denum   ! svgrd valeurs fct f a gche du centre
               		fv2denum(jtwm1) = f2denum   ! svgrd valeurs fct f a drte du centre
	       		reskdenum = reskdenum + wgk(jtwm1)*(f1denum+f2denum)
               		
			fv101num(jtwm1) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtwm1) = f201num   ! svgrd valeurs fct f a drte du centre
	       		resk01num = resk01num + wgk(jtwm1)*(f101num+f201num)
               		

			fv102num(jtwm1) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtwm1) = f202num   ! svgrd valeurs fct f a drte du centre
	       		resk02num = resk02num + wgk(jtwm1)*(f102num+f202num)
               		

			fv112num(jtwm1) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtwm1) = f212num   ! svgrd valeurs fct f a drte du centre
	       		resk12num = resk12num + wgk(jtwm1)*(f112num+f212num)
               		

			fv10101num(jtwm1) = f10101num   ! svgrd valeurs fct f a gche du centre
               		fv20101num(jtwm1) = f20101num   ! svgrd valeurs fct f a drte du centre
               		resk0101num = resk0101num + wgk(jtwm1)*(f10101num+f20101num)
               		
			fv10202num(jtwm1) = f10202num   ! svgrd valeurs fct f a gche du centre
               		fv20202num(jtwm1) = f20202num   ! svgrd valeurs fct f a drte du centre
               		resk0202num = resk0202num + wgk(jtwm1)*(f10202num+f20202num)
               		
			fv11212num(jtwm1) = f11212num   ! svgrd valeurs fct f a gche du centre
               		fv21212num(jtwm1) = f21212num   ! svgrd valeurs fct f a drte du centre
               		resk1212num = resk1212num + wgk(jtwm1)*(f11212num+f21212num)
               		




			
         	end do

    		resdenum = xr*reskdenum
    		res01num = xr*resk01num 
    		res02num = xr*resk02num
    		res12num = xr*resk12num 
		res0101num = xr*resk0101num
		res0202num = xr*resk0202num
		res1212num = xr*resk1212num  
	
         
              
          end subroutine qgaussweibderivdiag

!================================  QGAUS : 1  15  ==========================
!================================ complete hessian =========================
!================================ M-spline baseline risk ===================


subroutine qgausssplinederiv(a,b,the01,the02,the12,resdenum,&
		res01num,res02num,res12num,res0101num,res0102num,res0112num,&
		res0202num,res0212num,res1212num,v01,v02,v12)

   use commun,only:zi01,zi12,zi02,nz01,nz12,nz02

        double precision a,b
         double precision dx,xm,xr,reskdenum,resdenum, & 
	resk01num,res01num,resk02num,res02num,resk12num,res12num, & 
	resk0101num,res0101num,resk0102num,res0102num, &
	resk0112num,res0112num,resk0202num,res0202num, &
	resk0212num,res0212num,resk1212num,res1212num, &
	v01,v02,v12
         double precision xx,f1denum,f2denum, f101num, f102num, f112num, & 
	f10101num,f10102num,f10112num,f10202num,f10212num,f11212num,&
	f201num, f202num, f212num, f20101num,f20102num,f20112num,&
	f20202num,f20212num,f21212num,&
	su01,ri01,ri12,su12,su02,ri02,fv1denum,fv2denum, &
	fv101num,fv102num,fv112num,fv201num,fv202num,fv212num, & 
	fv10101num,fv20101num, fv10102num,fv20102num, & 
	fv10112num,fv20112num, fv10202num,fv20202num, & 
	fv10212num,fv20212num, fv11212num,fv21212num, &  
	d1mach(5),epmach,uflow,fcdenum,fc01num,fc02num,fc12num, & 
	fc0101num,fc0102num,fc0112num,fc0202num,fc0212num,fc1212num

         double precision gl01,gl12,gl02
         integer::j,jtw,jtwm1
         double precision,dimension(8)::xgk,wgk
	 double precision,dimension(4)::wg
	save wgk,xgk

	 dimension fv1denum(7),fv2denum(7),fv101num(7),fv201num(7), &
	fv102num(7),fv202num(7), fv112num(7),fv212num(7), & 
	fv10101num(7),fv20101num(7),fv10102num(7),fv20102num(7),&
	fv10112num(7),fv20112num(7),fv10202num(7),fv20202num(7),&
	fv10212num(7),fv20212num(7),fv11212num(7),fv21212num(7)

	
         double precision,dimension(-2:(nz01-1))::the01
         double precision,dimension(-2:(nz12-1))::the12
         double precision,dimension(-2:(nz02-1))::the02

         

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)

	wg(1)=0.129484966168869693270611432679082d0
   	wg(2)=0.279705391489276667901467771423780d0
    	wg(3)=0.381830050505118944950369775488975d0
    	wg(4)=0.417959183673469387755102040816327d0

    	xgk(1)=0.991455371120812639206854697526329d0
    	xgk(2)=0.949107912342758524526189684047851d0
    	xgk(3)=0.864864423359769072789712788640926d0
    	xgk(4)=0.741531185599394439863864773280788d0
    	xgk(5)=0.586087235467691130294144838258730d0
    	xgk(6)=0.405845151377397166906606412076961d0
    	xgk(7)=0.207784955007898467600689403773245d0
    	xgk(8)=0.000000000000000000000000000000000d0

    	wgk(1)=0.022935322010529224963732008058970d0
    	wgk(2)=0.063092092629978553290700663189204d0
    	wgk(3)=0.104790010322250183839876322541518d0
    	wgk(4)=0.140653259715525918745189590510238d0
    	wgk(5)=0.169004726639267902826583426598550d0
    	wgk(6)=0.190350578064785409913256402421014d0
    	wgk(7)=0.204432940075298892414161999234649d0
    	wgk(8)=0.209482141084727828012999174891714d0

              

	
	xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)

            resdenum = 0.d0
	    res01num = 0.d0
	    res02num = 0.d0
	    res12num = 0.d0
	
	    res0101num = 0.d0
	    res0102num = 0.d0
	    res0112num = 0.d0
	    res0202num = 0.d0
	    res0212num = 0.d0
	    res1212num = 0.d0
        
    	    	call susp(xm,the01,nz01,su01,ri01,zi01,gl01)
        	call susp(xm,the02,nz02,su02,ri02,zi02,gl02)
        	call susp(xm,the12,nz12,su12,ri12,zi12,gl12)
        
    		fcdenum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0
    		fc01num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
		fc02num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
		fc12num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

		
		fc0101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
		fc0102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl01*v01/(su12**v12)
		fc0112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl01*v01/(su12**v12)
		fc0202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
		fc0212num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl12*v12/(su12**v12)
		fc1212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

    		reskdenum = fcdenum*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk01num = fc01num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk02num = fc02num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk12num = fc12num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0101num = fc0101num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0102num = fc0102num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0112num = fc0112num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0202num = fc0202num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0212num = fc0212num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk1212num = fc1212num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	
		do j=1,3
	       		jtw = j*2
               		dx=xr*xgk(jtw)
               		xx = xm+dx
               		
    	    		call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
        		call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
        		call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
                        f1denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f101num=(su01**v01)*(su02**v02)*ri01*v01*(1-gl01*v01)/(su12**v12)
			f102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f10101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f10102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl01*v01/(su12**v12)
			f10112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl01*v01/(su12**v12)
			f10202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f10212num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl12*v12/(su12**v12)
			f11212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

			
               		xx = xm-dx
               		call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
        		call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
        		call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
                        
			f2denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f201num=(su01**v01)*(su02**v02)*ri01*v01*(1-gl01*v01)/(su12**v12)
			f202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f20101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f20102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl01*v01/(su12**v12)
			f20112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl01*v01/(su12**v12)
			f20202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f20212num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl12*v12/(su12**v12)
			f21212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

			
               		fv1denum(jtw) = f1denum   ! svgrd valeurs fct f a gche du centre
               		fv2denum(jtw) = f2denum   ! svgrd valeurs fct f a drte du centre

	       		reskdenum = reskdenum + wgk(jtw)*(f1denum+f2denum)
               		
			fv101num(jtw) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtw) = f201num   ! svgrd valeurs fct f a drte du centre

	       		resk01num = resk01num + wgk(jtw)*(f101num+f201num)
               		
			fv102num(jtw) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtw) = f202num   ! svgrd valeurs fct f a drte du centre

	       		resk02num = resk02num + wgk(jtw)*(f102num+f202num)
               		
			fv112num(jtw) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtw) = f212num   ! svgrd valeurs fct f a drte du centre

	       		resk12num = resk12num + wgk(jtw)*(f112num+f212num)
               		
			fv10101num(jtw) = f10101num   ! svgrd valeurs fct f a gche du centre
               		fv20101num(jtw) = f20101num   ! svgrd valeurs fct f a drte du centre

	       		resk0101num = resk0101num + wgk(jtw)*(f10101num+f20101num)
               		
			fv10102num(jtw) = f10102num   ! svgrd valeurs fct f a gche du centre
               		fv20102num(jtw) = f20102num   ! svgrd valeurs fct f a drte du centre

	       		resk0102num = resk0102num + wgk(jtw)*(f10102num+f20102num)
               		
			fv10112num(jtw) = f10112num   ! svgrd valeurs fct f a gche du centre
               		fv20112num(jtw) = f20112num   ! svgrd valeurs fct f a drte du centre

	       		resk0112num = resk0112num + wgk(jtw)*(f10112num+f20112num)
               		
			fv10202num(jtw) = f10202num   ! svgrd valeurs fct f a gche du centre
               		fv20202num(jtw) = f20202num   ! svgrd valeurs fct f a drte du centre

	       		resk0202num = resk0202num + wgk(jtw)*(f10202num+f20202num)
               		
			fv10212num(jtw) = f10212num   ! svgrd valeurs fct f a gche du centre
               		fv20212num(jtw) = f20212num   ! svgrd valeurs fct f a drte du centre

	       		resk0212num = resk0212num + wgk(jtw)*(f10212num+f20212num)
               		
			fv11212num(jtw) = f11212num   ! svgrd valeurs fct f a gche du centre
               		fv21212num(jtw) = f21212num   ! svgrd valeurs fct f a drte du centre

	       		resk1212num = resk1212num + wgk(jtw)*(f11212num+f21212num)
               		

			
         	end do
	 	do j=1,4
			jtwm1 = j*2-1
               		dx=xr*xgk(jtwm1)
               		xx = xm+dx
               		call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
        		call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
        		call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
      			f1denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f101num=(su01**v01)*(su02**v02)*ri01*v01*(1-gl01*v01)/(su12**v12)
			f102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f10101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f10102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl01*v01/(su12**v12)
			f10112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl01*v01/(su12**v12)
			f10202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f10212num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl12*v12/(su12**v12)
			f11212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

			
               		xx = xm-dx
               		call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
        		call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
        		call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
                        f2denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f201num=(su01**v01)*(su02**v02)*ri01*v01*(1-gl01*v01)/(su12**v12)
			f202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f20101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f20102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl01*v01/(su12**v12)
			f20112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl01*v01/(su12**v12)
			f20202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f20212num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl12*v12/(su12**v12)
			f21212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)


               		fv1denum(jtwm1) = f1denum   ! svgrd valeurs fct f a gche du centre
               		fv2denum(jtwm1) = f2denum   ! svgrd valeurs fct f a drte du centre
	       		reskdenum = reskdenum + wgk(jtwm1)*(f1denum+f2denum)
               		
			fv101num(jtwm1) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtwm1) = f201num   ! svgrd valeurs fct f a drte du centre
	       		resk01num = resk01num + wgk(jtwm1)*(f101num+f201num)
               		
			fv102num(jtwm1) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtwm1) = f202num   ! svgrd valeurs fct f a drte du centre
	       		resk02num = resk02num + wgk(jtwm1)*(f102num+f202num)
               		
			fv112num(jtwm1) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtwm1) = f212num   ! svgrd valeurs fct f a drte du centre
	       		resk12num = resk12num + wgk(jtwm1)*(f112num+f212num)
               		
			fv10101num(jtwm1) = f10101num   ! svgrd valeurs fct f a gche du centre
               		fv20101num(jtwm1) = f20101num   ! svgrd valeurs fct f a drte du centre
               		resk0101num = resk0101num + wgk(jtwm1)*(f10101num+f20101num)

			fv10102num(jtwm1) = f10102num   ! svgrd valeurs fct f a gche du centre
               		fv20102num(jtwm1) = f20102num   ! svgrd valeurs fct f a drte du centre
               		resk0102num = resk0102num + wgk(jtwm1)*(f10102num+f20102num)
               		
			fv10112num(jtwm1) = f10112num   ! svgrd valeurs fct f a gche du centre
               		fv20112num(jtwm1) = f20112num   ! svgrd valeurs fct f a drte du centre
               		resk0112num = resk0112num + wgk(jtwm1)*(f10112num+f20112num)
               		
			fv10202num(jtwm1) = f10202num   ! svgrd valeurs fct f a gche du centre
               		fv20202num(jtwm1) = f20202num   ! svgrd valeurs fct f a drte du centre
               		resk0202num = resk0202num + wgk(jtwm1)*(f10202num+f20202num)
               		
			fv10212num(jtwm1) = f10212num   ! svgrd valeurs fct f a gche du centre
               		fv20212num(jtwm1) = f20212num   ! svgrd valeurs fct f a drte du centre
               		resk0212num = resk0212num + wgk(jtwm1)*(f10212num+f20212num)
               		
			fv11212num(jtwm1) = f11212num   ! svgrd valeurs fct f a gche du centre
               		fv21212num(jtwm1) = f21212num   ! svgrd valeurs fct f a drte du centre
               		resk1212num = resk1212num + wgk(jtwm1)*(f11212num+f21212num)
               		
			
         	end do
     
    
    		resdenum = reskdenum*xr
    		res01num = resk01num*xr
    		res02num = resk02num*xr
    		res12num = resk12num*xr
		res0101num = resk0101num*xr
		res0102num = resk0102num*xr
		res0112num = resk0112num*xr
		res0202num = resk0202num*xr
		res0212num = resk0212num*xr
		res1212num = resk1212num*xr
	
         
     
	
          end subroutine qgausssplinederiv

!================================  QGAUS : 1  15  ==========================
!================================ only first derivatives  ===================
!================================ M-spline baseline risk ===================

subroutine qgausssplinefirstderiv(a,b,the01,the02,the12,resdenum,&
		res01num,res02num,res12num,v01,v02,v12)

   use commun,only:zi01,zi12,zi02,nz01,nz12,nz02

        double precision a,b
         double precision dx,xm,xr,reskdenum,resdenum, & 
	resk01num,res01num,resk02num,res02num,resk12num,res12num, & 
	v01,v02,v12
         double precision xx,f1denum,f2denum, f101num, f102num, f112num, & 
	f201num, f202num, f212num, &
	su01,ri01,ri12,su12,su02,ri02,fv1denum,fv2denum, &
	fv101num,fv102num,fv112num,fv201num,fv202num,fv212num, &  
	d1mach(5),epmach,uflow,fcdenum,fc01num,fc02num,fc12num

         double precision gl01,gl12,gl02
         integer::j,jtw,jtwm1
         double precision,dimension(8)::xgk,wgk
	 double precision,dimension(4)::wg
	save wgk,xgk

	 dimension fv1denum(7),fv2denum(7),fv101num(7),fv201num(7), &
	fv102num(7),fv202num(7), fv112num(7),fv212num(7)

	
         double precision,dimension(-2:(nz01-1))::the01
         double precision,dimension(-2:(nz12-1))::the12
         double precision,dimension(-2:(nz02-1))::the02

         

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)

	wg(1)=0.129484966168869693270611432679082d0
   	wg(2)=0.279705391489276667901467771423780d0
    	wg(3)=0.381830050505118944950369775488975d0
    	wg(4)=0.417959183673469387755102040816327d0

    	xgk(1)=0.991455371120812639206854697526329d0
    	xgk(2)=0.949107912342758524526189684047851d0
    	xgk(3)=0.864864423359769072789712788640926d0
    	xgk(4)=0.741531185599394439863864773280788d0
    	xgk(5)=0.586087235467691130294144838258730d0
    	xgk(6)=0.405845151377397166906606412076961d0
    	xgk(7)=0.207784955007898467600689403773245d0
    	xgk(8)=0.000000000000000000000000000000000d0

    	wgk(1)=0.022935322010529224963732008058970d0
    	wgk(2)=0.063092092629978553290700663189204d0
    	wgk(3)=0.104790010322250183839876322541518d0
    	wgk(4)=0.140653259715525918745189590510238d0
    	wgk(5)=0.169004726639267902826583426598550d0
    	wgk(6)=0.190350578064785409913256402421014d0
    	wgk(7)=0.204432940075298892414161999234649d0
    	wgk(8)=0.209482141084727828012999174891714d0

              

	
	xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)

            resdenum = 0.d0
	    res01num = 0.d0
	    res02num = 0.d0
	    res12num = 0.d0
	
        
    	    	call susp(xm,the01,nz01,su01,ri01,zi01,gl01)
        	call susp(xm,the02,nz02,su02,ri02,zi02,gl02)
        	call susp(xm,the12,nz12,su12,ri12,zi12,gl12)
        
    		fcdenum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0
    		fc01num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
		fc02num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
		fc12num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

    		reskdenum = fcdenum*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk01num = fc01num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk02num = fc02num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk12num = fc12num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	
		do j=1,3
	       		jtw = j*2
               		dx=xr*xgk(jtw)
               		xx = xm+dx
               		
    	    		call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
        		call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
        		call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
                        f1denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f101num=(su01**v01)*(su02**v02)*ri01*v01*(1-gl01*v01)/(su12**v12)
			f102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 
			
               		xx = xm-dx
               		call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
        		call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
        		call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
                        
			f2denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f201num=(su01**v01)*(su02**v02)*ri01*v01*(1-gl01*v01)/(su12**v12)
			f202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

               		fv1denum(jtw) = f1denum   ! svgrd valeurs fct f a gche du centre
               		fv2denum(jtw) = f2denum   ! svgrd valeurs fct f a drte du centre

	       		reskdenum = reskdenum + wgk(jtw)*(f1denum+f2denum)
               		
			fv101num(jtw) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtw) = f201num   ! svgrd valeurs fct f a drte du centre

	       		resk01num = resk01num + wgk(jtw)*(f101num+f201num)
               		
			fv102num(jtw) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtw) = f202num   ! svgrd valeurs fct f a drte du centre

	       		resk02num = resk02num + wgk(jtw)*(f102num+f202num)
               		
			fv112num(jtw) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtw) = f212num   ! svgrd valeurs fct f a drte du centre

	       		resk12num = resk12num + wgk(jtw)*(f112num+f212num)

         	end do
	 	do j=1,4
			jtwm1 = j*2-1
               		dx=xr*xgk(jtwm1)
               		xx = xm+dx
               		call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
        		call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
        		call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
      			f1denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f101num=(su01**v01)*(su02**v02)*ri01*v01*(1-gl01*v01)/(su12**v12)
			f102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

               		xx = xm-dx
               		call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
        		call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
        		call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
                        f2denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f201num=(su01**v01)*(su02**v02)*ri01*v01*(1-gl01*v01)/(su12**v12)
			f202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

               		fv1denum(jtwm1) = f1denum   ! svgrd valeurs fct f a gche du centre
               		fv2denum(jtwm1) = f2denum   ! svgrd valeurs fct f a drte du centre
	       		reskdenum = reskdenum + wgk(jtwm1)*(f1denum+f2denum)
               		
			fv101num(jtwm1) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtwm1) = f201num   ! svgrd valeurs fct f a drte du centre
	       		resk01num = resk01num + wgk(jtwm1)*(f101num+f201num)
               		
			fv102num(jtwm1) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtwm1) = f202num   ! svgrd valeurs fct f a drte du centre
	       		resk02num = resk02num + wgk(jtwm1)*(f102num+f202num)
               		
			fv112num(jtwm1) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtwm1) = f212num   ! svgrd valeurs fct f a drte du centre
	       		resk12num = resk12num + wgk(jtwm1)*(f112num+f212num)
               		
			
         	end do
     
    
    		resdenum = reskdenum*xr
    		res01num = resk01num*xr
    		res02num = resk02num*xr
    		res12num = resk12num*xr
	
	
          end subroutine qgausssplinefirstderiv

!================================  QGAUS : 1  15  ==========================
!================================ only diagonal terms of hessian ===========
!================================ M-spline baseline risk ===================


subroutine qgausssplinederivdiag(a,b,the01,the02,the12,resdenum,&
		res01num,res02num,res12num,res0101num,&
		res0202num,res1212num,v01,v02,v12)

   use commun,only:zi01,zi12,zi02,nz01,nz12,nz02

        double precision a,b
         double precision dx,xm,xr,reskdenum,resdenum, & 
	resk01num,res01num,resk02num,res02num,resk12num,res12num, & 
	resk0101num,res0101num,resk0202num,res0202num, &
	resk1212num,res1212num,v01,v02,v12
         double precision xx,f1denum,f2denum, f101num, f102num, f112num, & 
	f10101num,f10202num,f11212num,&
	f201num, f202num, f212num, f20101num,&
	f20202num,f21212num,&
	su01,ri01,ri12,su12,su02,ri02,fv1denum,fv2denum, &
	fv101num,fv102num,fv112num,fv201num,fv202num,fv212num, & 
	fv10101num,fv20101num,fv10202num,fv20202num, & 
	fv11212num,fv21212num, &  
	d1mach(5),epmach,uflow,fcdenum,fc01num,fc02num,fc12num, & 
	fc0101num,fc0202num,fc1212num

         double precision gl01,gl12,gl02
         integer::j,jtw,jtwm1
         double precision,dimension(8)::xgk,wgk
	 double precision,dimension(4)::wg
	save wgk,xgk

	 dimension fv1denum(7),fv2denum(7),fv101num(7),fv201num(7), &
	fv102num(7),fv202num(7), fv112num(7),fv212num(7), & 
	fv10101num(7),fv20101num(7),fv10202num(7),fv20202num(7),&
	fv11212num(7),fv21212num(7)

	
         double precision,dimension(-2:(nz01-1))::the01
         double precision,dimension(-2:(nz12-1))::the12
         double precision,dimension(-2:(nz02-1))::the02

         

   	D1MACH(1)=2.23D-308
    	D1MACH(2)=1.79D+308
    	D1MACH(3)=1.11D-16
    	D1MACH(4)=2.22D-16
    	D1MACH(5)=0.301029995663981195D0

    	epmach = d1mach(4)
    	uflow = d1mach(1)

	wg(1)=0.129484966168869693270611432679082d0
   	wg(2)=0.279705391489276667901467771423780d0
    	wg(3)=0.381830050505118944950369775488975d0
    	wg(4)=0.417959183673469387755102040816327d0

    	xgk(1)=0.991455371120812639206854697526329d0
    	xgk(2)=0.949107912342758524526189684047851d0
    	xgk(3)=0.864864423359769072789712788640926d0
    	xgk(4)=0.741531185599394439863864773280788d0
    	xgk(5)=0.586087235467691130294144838258730d0
    	xgk(6)=0.405845151377397166906606412076961d0
    	xgk(7)=0.207784955007898467600689403773245d0
    	xgk(8)=0.000000000000000000000000000000000d0

    	wgk(1)=0.022935322010529224963732008058970d0
    	wgk(2)=0.063092092629978553290700663189204d0
    	wgk(3)=0.104790010322250183839876322541518d0
    	wgk(4)=0.140653259715525918745189590510238d0
    	wgk(5)=0.169004726639267902826583426598550d0
    	wgk(6)=0.190350578064785409913256402421014d0
    	wgk(7)=0.204432940075298892414161999234649d0
    	wgk(8)=0.209482141084727828012999174891714d0

              

	
	xm = 0.5d+00*(b+a)
        xr = 0.5d+00*(b-a)

            resdenum = 0.d0
	    res01num = 0.d0
	    res02num = 0.d0
	    res12num = 0.d0
	
	    res0101num = 0.d0
	    res0202num = 0.d0
	    res1212num = 0.d0
        
    	    	call susp(xm,the01,nz01,su01,ri01,zi01,gl01)
        	call susp(xm,the02,nz02,su02,ri02,zi02,gl02)
        	call susp(xm,the12,nz12,su12,ri12,zi12,gl12)
        
    		fcdenum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)  ! valeur fct f au milieu de intervalle (a,b), cas pnt 0
    		fc01num=(su01**v01)*(su02**v02)*ri01*v01*(1-(gl01*v01))/(su12**v12)
		fc02num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
		fc12num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

		
		fc0101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
		fc0202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
		fc1212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

    		reskdenum = fcdenum*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk01num = fc01num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk02num = fc02num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk12num = fc12num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0101num = fc0101num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk0202num = fc0202num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resk1212num = fc1212num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	
		do j=1,3
	       		jtw = j*2
               		dx=xr*xgk(jtw)
               		xx = xm+dx
               		
    	    		call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
        		call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
        		call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
                        f1denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f101num=(su01**v01)*(su02**v02)*ri01*v01*(1-gl01*v01)/(su12**v12)
			f102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f10101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f10202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f11212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

			
               		xx = xm-dx
               		call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
        		call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
        		call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
                        
			f2denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f201num=(su01**v01)*(su02**v02)*ri01*v01*(1-gl01*v01)/(su12**v12)
			f202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f20101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f20202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f21212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

			
               		fv1denum(jtw) = f1denum   ! svgrd valeurs fct f a gche du centre
               		fv2denum(jtw) = f2denum   ! svgrd valeurs fct f a drte du centre

	       		reskdenum = reskdenum + wgk(jtw)*(f1denum+f2denum)
               		
			fv101num(jtw) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtw) = f201num   ! svgrd valeurs fct f a drte du centre

	       		resk01num = resk01num + wgk(jtw)*(f101num+f201num)
               		
			fv102num(jtw) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtw) = f202num   ! svgrd valeurs fct f a drte du centre

	       		resk02num = resk02num + wgk(jtw)*(f102num+f202num)
               		
			fv112num(jtw) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtw) = f212num   ! svgrd valeurs fct f a drte du centre

	       		resk12num = resk12num + wgk(jtw)*(f112num+f212num)
               		
			fv10101num(jtw) = f10101num   ! svgrd valeurs fct f a gche du centre
               		fv20101num(jtw) = f20101num   ! svgrd valeurs fct f a drte du centre

	       		resk0101num = resk0101num + wgk(jtw)*(f10101num+f20101num)
               		
			fv10202num(jtw) = f10202num   ! svgrd valeurs fct f a gche du centre
               		fv20202num(jtw) = f20202num   ! svgrd valeurs fct f a drte du centre

	       		resk0202num = resk0202num + wgk(jtw)*(f10202num+f20202num)
               		
			fv11212num(jtw) = f11212num   ! svgrd valeurs fct f a gche du centre
               		fv21212num(jtw) = f21212num   ! svgrd valeurs fct f a drte du centre

	       		resk1212num = resk1212num + wgk(jtw)*(f11212num+f21212num)
               		

			
         	end do
	 	do j=1,4
			jtwm1 = j*2-1
               		dx=xr*xgk(jtwm1)
               		xx = xm+dx
               		call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
        		call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
        		call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
      			f1denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f101num=(su01**v01)*(su02**v02)*ri01*v01*(1-gl01*v01)/(su12**v12)
			f102num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f112num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f10101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f10202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f11212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)

			
               		xx = xm-dx
               		call susp(xx,the01,nz01,su01,ri01,zi01,gl01)
        		call susp(xx,the02,nz02,su02,ri02,zi02,gl02)
        		call susp(xx,the12,nz12,su12,ri12,zi12,gl12)
                        f2denum =(su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
			f201num=(su01**v01)*(su02**v02)*ri01*v01*(1-gl01*v01)/(su12**v12)
			f202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02/(su12**v12) 
			f212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12/(su12**v12) 

			f20101num=(su01**v01)*(su02**v02)*ri01*v01*(((1-gl01*v01)**2)-&
			gl01*v01)/(su12**v12)
			f20202num=(su01**v01)*(su02**v02)*ri01*v01*gl02*v02*gl02*v02/(su12**v12)
			f21212num=(su01**v01)*(su02**v02)*ri01*v01*gl12*v12*gl12*v12/(su12**v12)


               		fv1denum(jtwm1) = f1denum   ! svgrd valeurs fct f a gche du centre
               		fv2denum(jtwm1) = f2denum   ! svgrd valeurs fct f a drte du centre
	       		reskdenum = reskdenum + wgk(jtwm1)*(f1denum+f2denum)
               		
			fv101num(jtwm1) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtwm1) = f201num   ! svgrd valeurs fct f a drte du centre
	       		resk01num = resk01num + wgk(jtwm1)*(f101num+f201num)
               		
			fv102num(jtwm1) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtwm1) = f202num   ! svgrd valeurs fct f a drte du centre
	       		resk02num = resk02num + wgk(jtwm1)*(f102num+f202num)
               		
			fv112num(jtwm1) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtwm1) = f212num   ! svgrd valeurs fct f a drte du centre
	       		resk12num = resk12num + wgk(jtwm1)*(f112num+f212num)
               		
			fv10101num(jtwm1) = f10101num   ! svgrd valeurs fct f a gche du centre
               		fv20101num(jtwm1) = f20101num   ! svgrd valeurs fct f a drte du centre
               		resk0101num = resk0101num + wgk(jtwm1)*(f10101num+f20101num)

			fv10202num(jtwm1) = f10202num   ! svgrd valeurs fct f a gche du centre
               		fv20202num(jtwm1) = f20202num   ! svgrd valeurs fct f a drte du centre
               		resk0202num = resk0202num + wgk(jtwm1)*(f10202num+f20202num)
			
			fv11212num(jtwm1) = f11212num   ! svgrd valeurs fct f a gche du centre
               		fv21212num(jtwm1) = f21212num   ! svgrd valeurs fct f a drte du centre
               		resk1212num = resk1212num + wgk(jtwm1)*(f11212num+f21212num)
               		
			
         	end do
     
    
    		resdenum = reskdenum*xr
    		res01num = resk01num*xr
    		res02num = resk02num*xr
    		res12num = resk12num*xr
		res0101num = resk0101num*xr
		res0202num = resk0202num*xr
		res1212num = resk1212num*xr
	
         
     
	
          end subroutine qgausssplinederivdiag
!=============================================================================================  
!======================= Calculate derivatives of loglik with weibull baseline risk ==========
!======================= only beta parameters ========================================================
!=============================================================================================  


subroutine derivaweib(b0,np0,npar0,bfix0,fix0,c0,no0,ve010,ve120,ve020,&
        dimnva01,dimnva12,dimnva02,nva01,nva12,nva02,t00,&
        t10,t20,t30,troncature0,likelihood_deriv)
	
	use commun
        implicit none
         
        double precision::res2denum,res201num,res202num,res212num, &
	res20101num,res20102num,res20112num,res20202num,res20212num,&
	res21212num,vet01,vet12,vet02,resint,v,u1,u2,u3
        integer::np0,i,j,l,w,k,lfix, kfix,npar0,nva01,nva12,nva02,no0, &
	nz010,nz020,nz120,troncature0,dimnva01,dimnva02,dimnva12, & 
	nva01nofix,nva12nofix,nva02nofix,nvamax, sizespline,nva0102,&
	nvamax01,nvamax0102,nvamax0112,nvamax02,nvamax0212,nvamax12

	double precision,dimension(np0+np0*(np0+1)/2),intent(inout)::likelihood_deriv
	double precision,dimension(np0)::b0
	double precision,dimension(np0+np0*(np0+1)/2)::res,res1
        double precision,dimension(npar0)::bh
	double precision,dimension(npar0-np0)::bfix0
	integer,dimension(npar0)::fix0
	double precision,dimension(2)::the01,the12,the02
        double precision,dimension(no0,dimnva01)::ve010
	double precision,dimension(no0,dimnva02)::ve020
	double precision,dimension(no0,dimnva12)::ve120
	
	
        double precision::su01,ri01,su12,ri12,su02,ri02,gl01,gl02,gl12
	double precision,dimension(no0)::t00,t10,t20,t30
	integer,dimension(no0)::c0

	allocate(b(np0),bfix(npar0-np0),fix(npar0))
	b=b0
	bfix=bfix0
	fix=fix0
	troncature=troncature0
	
	sizespline=6


	if(nva01.gt.0) then 
	  nva01nofix=nva01-sum(fix((sizespline+1):(nva01+sizespline)))
	else 
	  nva01nofix=0
	end if 

	if(nva02.gt.0) then 
          nva02nofix=nva02-sum(fix((nva01+sizespline+1):(nva02+nva01+sizespline)))
	else 
	   nva02nofix=0
	end if 

	if(nva12.gt.0) then 
	  nva12nofix=nva12-sum(fix((nva01+nva02+sizespline+1):npar0))
	else 
	  nva12nofix=0
	end if 

	nva0102=nva01nofix+nva02nofix
	nvamax=nva01nofix+nva02nofix+nva12nofix
	nvamax01=nvamax+(nva01nofix+1)*nva01nofix/2
	nvamax0102=nvamax01+nva01nofix*nva02nofix
	nvamax0112=nvamax0102+nva01nofix*nva12nofix
	nvamax02=nvamax0112+(nva02nofix+1)*nva02nofix/2
	nvamax0212=nvamax02+nva02nofix*nva12nofix
	nvamax12=nvamax0212+(nva12nofix+1)*nva12nofix/2


	if(nva01.gt.0) then 
		allocate(ve01(no0,nva01))
		allocate(ve01nofix(no0,nva01nofix))
		allocate(ve01square(no0,nva01nofix*(nva01nofix+1)/2))
		allocate(tronc01(nva01nofix))
		allocate(tronc01square(nva01nofix*(nva01nofix+1)/2))
	else 
		allocate(ve01(no0,1))
		allocate(ve01nofix(no0,1))
		ve01nofix=0
		allocate(ve01square(no0,1))
		ve01square=0
		allocate(tronc01(1))
		allocate(tronc01square(1))
	end if 
	
	if(nva02.gt.0) then 
		allocate(ve02(no0,nva02))
		allocate(ve02nofix(no0,nva02nofix))
		allocate(ve02square(no0,nva02nofix*(nva02nofix+1)/2))
		allocate(tronc02(nva02nofix))
		allocate(tronc02square(nva02nofix*(nva02nofix+1)/2))
	else 
		allocate(ve02(no0,1))
		allocate(ve02nofix(no0,1))
		ve02nofix=0
		allocate(ve02square(no0,1))
		ve02square=0
		allocate(tronc02(1))
		allocate(tronc02square(1))
	end if 

	if(nva12.gt.0) then 
		allocate(ve12(no0,nva12))
		allocate(ve12nofix(no0,nva12nofix))
		allocate(ve12square(no0,nva12nofix*(nva12nofix+1)/2))
	else 
		allocate(ve12(no0,1))
		allocate(ve12nofix(no0,1))
		ve12nofix=0
		allocate(ve12square(no0,1))
		ve12square=0
	end if 



	ve01=ve010
	ve02=ve020
	ve12=ve120

	allocate(t0(no0),t1(no0),t2(no0),t3(no0),c(no0))
	c=c0
	t0=t00
	t1=t10
	t2=t20
	t3=t30

         
        ! we need to put bh at its original values if in posfix 


       l=0
       lfix=0
       w=0

	if(nva01.gt.0) then
	do k=1,nva01
	   if(fix((sizespline+k)).eq.0) then 
		lfix=lfix+1
		ve01nofix(:,lfix)=ve01(:,k)
	   end if 
	end do
	
	do i=1,no0
		lfix=1
		do l=1,nva01nofix
	   	ve01square(i,lfix:(lfix+nva01nofix-l))=&
		ve01nofix(i,l:nva01nofix)
		ve01square(i,lfix:(lfix+nva01nofix-l))=&
		ve01square(i,lfix:(lfix+nva01nofix-l))*ve01nofix(i,l)
		lfix=lfix+nva01nofix-l+1
		end do
	end do
	 
        end if 
	lfix=0

	if(nva02.gt.0) then 
	do k=1,nva02
	   if(fix((sizespline+nva01+k)).eq.0) then 
		lfix=lfix+1
		ve02nofix(:,lfix)=ve02(:,k)
	   end if 
	end do

	do i=1,no0
		lfix=1
		do l=1,nva02nofix
	   	ve02square(i,lfix:(lfix+nva02nofix-l))=&
		ve02nofix(i,l)*ve02nofix(i,l:nva02nofix)
		lfix=lfix+nva02nofix-l+1
		end do
	end do

	end if 
	
	lfix=0

	if(nva12.gt.0) then 
	do k=1,nva12
	   if(fix((sizespline+nva01+nva02+k)).eq.0) then 
		lfix=lfix+1
		ve12nofix(:,lfix)=ve12(:,k)
	   end if 
	end do

	do i=1,no0
		lfix=1
		do l=1,nva12nofix
	   		ve12square(i,lfix:(lfix+nva12nofix-l))=&
			ve12nofix(i,l)*ve12nofix(i,l:nva12nofix)
			lfix=lfix+nva12nofix-l+1
		end do
	end do

	end if
	
	l=0
       lfix=0
       w=0

     do k=1,npar0 
         if(fix(k).eq.0) then
            l=l+1
            bh(k)=b(l)
	 end if 
         if(fix(k).eq.1) then
            w=w+1
            bh(k)=bfix(w)
         end if
      end do

   
	


         do i=1,2
            the01(i)=(bh(i))*(bh(i))
         end do
         do i=1,2
            j = 2+i
            the02(i)=(bh(j))*(bh(j))
         end do
         do i=1,2
            j = 4+i
            the12(i)=(bh(j))*(bh(j))
         end do



	endif
	
!---------- calcul des derivees premiere ------------------   

	res = 0
        do i=1,no0


                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
				vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
				vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
				vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

		
                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)

                res1 = 0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc01 = 0
				tronc02 =  0
                        	tronc01square= 0
                        	tronc02square=0
                        else 
				call fonct(t0(i),the01,ri01,gl01,su01)
				call fonct(t0(i),the02,ri02,gl02,su02)
                                tronc01=ve01nofix(i,:)*gl01*vet01
                        	tronc02=ve02nofix(i,:)*gl02*vet02
                        	tronc01square=ve01square(i,:)*gl01*vet01
                        	tronc02square=ve02square(i,:)*gl02*vet02
                        end if
                else
                        tronc01 = 0
                  	tronc02 = 0
                  	tronc01square=0
                  	tronc02square=0
                end if
		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
			call fonct(t1(i),the01,ri01,gl01,su01)
			call fonct(t1(i),the02,ri02,gl02,su02)

			if(nva01nofix.gt.0) then 
			
			res1(1:nva01nofix)=&
			-ve01nofix(i,:)*gl01*vet01+tronc01
			res1((nvamax+1):nvamax01)=&
			-ve01square(i,:)*gl01*vet01+&
			tronc01square

			if(nva02nofix.gt.0) then 
			res1((nvamax01+1):nvamax0102)=0
			end if 
			if(nva12nofix.gt.0) then 
			res1((nvamax0102+1):nvamax0112)=0
			end if

			end if 

			if(nva02nofix.gt.0) then
			res1((nva01nofix+1):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02

			res1((nvamax0112+1):nvamax02)=&
			-ve02square(i,:)*gl02*vet02+&
			tronc02square
			if(nva12nofix.gt.0) then 
			res1((nvamax02+1):nvamax0212)=0
			end if 
			end if 

			if(nva12nofix.gt.0) then 
			res1((nva0102+1):nvamax)=0
			res1((nvamax0212+1):nvamax12)=0
			end if 
			

                else
                if(c(i).eq.2)then ! cpi 0-->1


			call fonct(t3(i),the12,ri12,gl12,su12)
			call qgaussweibderiv(t1(i),t2(i),the01,&
			the02,the12,res2denum,res201num,&
			res202num,res212num,res20101num,&
			res20102num,res20112num,res20202num,&
			res20212num,res21212num,&
			vet01,vet02,vet12)
                        
			v=res2denum*(su12**vet12)

			if(nva01nofix.gt.0) then

			u1=res201num*(su12**vet12)
			
      			res1(1:nva01nofix)=&
			ve01nofix(i,:)*u1/v
			res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01

			res1((nvamax+1):nvamax01)=&
			ve01square(i,:)*res20101num*(su12**vet12)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)/v
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)+&
			tronc01square

			end if 

			
			if(nva02nofix.gt.0) then

			u2=-res202num*(su12**vet12)
			res1((nva01nofix+1):nva0102)=&
			ve02nofix(i,:)*u2/v
			res1((nva01nofix+1):nva0102)=&
			res1((nva01nofix+1):nva0102)+&
			tronc02

			res1((nvamax0112+1):nvamax02)=&
			-ve02square(i,:)*(su12**vet12)*(-res20202num+res202num)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)/v
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)+tronc02square


			end if 

			
			if(nva12nofix.gt.0) then 

			u3=res212num-gl12*vet12*res2denum
			u3=u3*(su12**vet12)

			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v
			res1((nvamax0212+1):nvamax12)=&
			-2*gl12*vet12*res212num+res212num+&
			res21212num-gl12*vet12*(1-gl12*vet12)*res2denum
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*(su12**vet12)*ve12square(i,:)
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)/v
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)

			
			end if 
			
			

			
			if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 

			kfix=nvamax01+1
			lfix=kfix-1+nva02nofix

			do j=1,nva01nofix

			   res1(kfix:lfix)=&
			  (res20102num-res202num)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva02nofix

			end do

			end if 

			if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 

			kfix=nvamax0102+1
			lfix=kfix-1+nva12nofix

			do j=1,nva01nofix
			   res1(kfix:lfix)=-gl12*vet12*res201num-res20112num+&
			   res212num
			   res1(kfix:lfix)=res1(kfix:lfix)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do

			end if 

			
			if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 

			kfix=nvamax02+1
			lfix=kfix-1+nva12nofix

			do j=1,nva02nofix
			   res1(kfix:lfix)=-res20212num+&
			   gl12*vet12*res202num
			   res1(kfix:lfix)=res1(kfix:lfix)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve02nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u2*u3*ve02nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do

			end if 
			
                else  
                    if(c(i).eq.3)then ! obs 0-->1

			call fonct(t1(i),the01,ri01,gl01,su01)
			call fonct(t1(i),the02,ri02,gl02,su02)
			call fonct(t1(i),the12,ri12,gl12,su12)

			if(nva01nofix.gt.0) then 

			res1(1:nva01nofix)=-ve01nofix(i,:)*gl01*vet01+&
			tronc01+ve01nofix(i,:)
			res1((nvamax+1):nvamax01)=&
			-ve01square(i,:)*gl01*vet01+&
			tronc01square

			if(nva02nofix.gt.0) then 
			res1((nvamax01+1):nvamax0102)=0
			end if 

			if(nva12nofix.gt.0)then 
			res1((nvamax0102+1):nvamax0112)=0
			end if 

			end if 

			if(nva02nofix.gt.0) then 

			res1((nva01nofix+1):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02
			res1((nvamax0112+1):nvamax02)=&
			-ve02square(i,:)*gl02*vet02+&
			tronc02square

			if(nva12nofix.gt.0) then 
			res1((nvamax02+1):nvamax0212)=0
			end if 

			end if 

			if(nva12nofix.gt.0) then 

			res1((nva0102+1):nvamax)=gl12*vet12
			res1((nvamax0212+1):nvamax12)=gl12*vet12

			call fonct(t3(i),the12,ri12,gl12,su12)
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)-gl12*vet12
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)*ve12nofix(i,:)

			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)-gl12*vet12
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*ve12square(i,:)
			

			end if 
			
			

                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			
			call fonct(t3(i),the12,ri12,gl12,su12)
			call qgaussweibderiv(t1(i),t2(i),the01,the02,the12,&
			res2denum,res201num,res202num,res212num,&
			res20101num,res20102num,&
			res20112num,res20202num,res20212num,&
			res21212num,vet01,vet02,vet12)

			v=res2denum*(su12**vet12)*ri12*vet12

			if(nva01nofix.gt.0) then 
			u1=res201num*(su12**vet12)*ri12*vet12
			
			res1(1:nva01nofix)=&
			ve01nofix(i,:)*u1/v
			res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01

			res1((nvamax+1):nvamax01)=&
			res20101num*(su12**vet12)*ri12*vet12*ve01square(i,:)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)/v
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)+tronc01square

			end if 

			
			if(nva02nofix.gt.0) then 

			u2=-res202num*(su12**vet12)*ri12*vet12
			
			res1((nva01nofix+1):(nva01nofix+nva02nofix))=&
			ve02nofix(i,:)*u2/v
			res1((nva01nofix+1):nva0102)=&
			res1((nva01nofix+1):nva0102)+&
			tronc02
			
			res1((nvamax0112+1):nvamax02)=&
			(res20202num-res202num)*ri12*vet12*(su12**vet12)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)*ve02square(i,:)	
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)/v
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)+tronc02square

			end if 


			
			if(nva12nofix.gt.0) then 

			u3=res212num+&
			(1-gl12*vet12)*res2denum
			u3=u3*(su12**vet12)*ri12*vet12
			
			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v
			
			res1((nvamax0212+1):nvamax12)=&
			(1-gl12*vet12)*res212num+res21212num+&
			((1-gl12*vet12)**2)*res2denum+&
			(2-gl12*vet12)*res212num-gl12*vet12*res2denum
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*(su12**vet12)*ri12*vet12
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*ve12square(i,:)
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)/v
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)			


			end if 

			
			
			if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 
			kfix=nvamax01+1
			lfix=kfix-1+nva02nofix

			do j=1,nva01nofix
			   res1(kfix:lfix)=(-res202num+res20102num)*ri12*vet12
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva02nofix
			end do

			end if 

			if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 
			
			kfix=nvamax0102+1
			lfix=kfix-1+nva12nofix

			do j=1,nva01nofix
			   res1(kfix:lfix)=&
			   (1-gl12*vet12)*res201num+&
			   res212num-res20112num
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve12nofix(i,:)
		           res1(kfix:lfix)=&
			   res1(kfix:lfix)*ri12*vet12*(su12**vet12)
			   res1(kfix:lfix)=res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do


			end if 

			if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 
			
			kfix=nvamax02+1
			lfix=kfix-1+nva12nofix

			do j=1,nva02nofix
			   res1(kfix:lfix)=&
			   res20212num+(1-gl12*vet12)*res202num
			   res1(kfix:lfix)=&
			   -res1(kfix:lfix)*ve02nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ri12*vet12*(su12**vet12)
			   res1(kfix:lfix)=res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-ve02nofix(i,j)*ve12nofix(i,:)*u2*u3
			   res1(kfix:lfix)=res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do

			end if 


                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
				call fonct(t1(i),the01,ri01,gl01,su01)
				call fonct(t1(i),the02,ri02,gl02,su02)
				call fonct(t1(i),the12,ri12,gl12,su12)
				
				if(nva01nofix.gt.0) then 

				res1(1:nva01nofix)=&
				-ve01nofix(i,:)*gl01*vet01+&
				ve01nofix(i,:)+&
				tronc01
				res1((nvamax+1):nvamax01)=&
				-ve01square(i,:)*gl01*vet01+&
				tronc01square

				if(nva02nofix.gt.0) then 
				res1((nvamax01+1):nvamax0102)=0
				end if 

				if(nva12nofix.gt.0) then
				res1((nvamax0102+1):nvamax0112)=0
				end if 

				end if 

				if(nva02nofix.gt.0) then 

				res1((nva01nofix+1):nva0102)=&
				-ve02nofix(i,:)*gl02*vet02+tronc02
				

				res1((nvamax0112+1):nvamax02)=&
				-ve02square(i,:)*gl02*vet02+&
				tronc02square

				if(nva12nofix.gt.0) then
				res1((nvamax02+1):nvamax0212)=0
				end if 

				end if 
			
				if(nva12nofix.gt.0) then 

				res1((nvamax0212+1):nvamax12)=&
				gl12*vet12
				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*gl12*vet12+&
				ve12nofix(i,:)

				call fonct(t3(i),the12,ri12,gl12,su12)
				res1((nva0102+1):nvamax)=&
				res1((nva0102+1):nvamax)-ve12nofix(i,:)*gl12*vet12
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)-gl12*vet12
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)*ve12square(i,:)
				

				end if 
				
				
				
                         else
                            if(c(i).eq.6)then ! vivant ???


				call fonct(t3(i),the01,ri01,gl01,su01)
				call fonct(t3(i),the02,ri02,gl02,su02)
				call fonct(t3(i),the12,ri12,gl12,su12)
				call qgaussweibderiv(t1(i),t3(i),the01,the02,the12,&
				res2denum,res201num,res202num,res212num,res20101num,&
				res20102num,res20112num,res20202num,res20212num,&
				res21212num,vet01,vet02,vet12)

				
				v=(su12**vet12)*res2denum+&
				(su01**vet01)*(su02**vet02)

				if(nva01nofix.gt.0) then 

				u1=(-gl01*vet01)*(su01**vet01)*(su02**vet02)+&
				(su12**vet12)*res201num
			        res1(1:nva01nofix)=&
				ve01nofix(i,:)*u1/v
				res1(1:nva01nofix)=&
				res1(1:nva01nofix)+tronc01

				res1((nvamax+1):nvamax01)=&
				res20101num*(su12**vet12)
				resint=(su01**vet01)*(su02**vet02)
				resint=resint*gl01*vet01*(1-gl01*vet01)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-resint
				
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)*ve01square(i,:)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)/v

				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)+tronc01square

				end if 
				
				if(nva02nofix.gt.0) then 

				u2=-gl02*vet02*(su01**vet01)*(su02**vet02)
				u2=u2-(su12**vet12)*res202num
				res1((nva01nofix+1):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1):nva0102)=&
				res1((nva01nofix+1):nva0102)+&
				tronc02

				res1((nvamax0112+1):nvamax02)=&
				(su12**vet12)*(res20202num-res202num)-&
				gl02*vet02*(su01**vet01)*(su02**vet02)*(1-gl02*vet02)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)*ve02square(i,:)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)/v

				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)+tronc02square

				
				end if 

				if(nva12nofix.gt.0) then 
				
				u3=-gl12*vet12*(su12**vet12)*res2denum+&
				(su12**vet12)*res212num
				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v

				res1((nvamax0212+1):nvamax12)=&
				res21212num+(1-gl12*vet12)*res212num-&
				gl12*vet12*(1-gl12*vet12)*res2denum-&
				gl12*vet12*res212num
				resint=(su12**vet12)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)*resint*ve12square(i,:)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)/v
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)


				end if 
                        	
				
				if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 

				kfix=nvamax01+1
				lfix=kfix-1+nva02nofix

				do j=1,nva01nofix
			   		res1(kfix:lfix)=&
					-(res202num-res20102num)*(su12**vet12)
					resint=(su01**vet01)*(su02**vet02)
					res1(kfix:lfix)=res1(kfix:lfix)+&
					gl01*vet01*gl02*vet02*resint
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)*v
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva02nofix
				end do

				end if 

				if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 

				
				kfix=nvamax0102+1
				lfix=kfix-1+nva12nofix

				do j=1,nva01nofix
			  		res1(kfix:lfix)=-gl12*vet12*res201num+&
					res212num-res20112num
					resint=(su12**vet12)
			   		res1(kfix:lfix)=&
			  		res1(kfix:lfix)*resint*ve01nofix(i,j)*ve12nofix(i,:)
		           		res1(kfix:lfix)=&
			   		res1(kfix:lfix)*v
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
			   		res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva12nofix
				end do

				end if 

				
				if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 

				kfix=nvamax02+1
				lfix=kfix-1+nva12nofix

				do j=1,nva02nofix
			  		res1(kfix:lfix)=&
					-res20212num+gl12*vet12*res202num
					resint=(su12**vet12)
					res1(kfix:lfix)=&
					res1(kfix:lfix)*resint*ve02nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
					res1(kfix:lfix)*v
					res1(kfix:lfix)=&
					res1(kfix:lfix)-u2*u3*ve02nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva12nofix
				end do

				end if 


				
                            else ! passage 0-->2  
				
				call fonct(t3(i),the01,ri01,gl01,su01)
				call fonct(t3(i),the02,ri02,gl02,su02)
				call fonct(t3(i),the12,ri12,gl12,su12)
        			call qgaussweibderiv(t1(i),t3(i),the01,the02,&
        			the12,res2denum,res201num,res202num,res212num,&
        			res20101num,res20102num,&
			  	res20112num,res20202num,res20212num,res21212num,&
        			vet01,vet02,vet12)
				
				
				v=(su12**vet12)*ri12*vet12*res2denum+&
				(su01**vet01)*(su02**vet02)*ri02*vet02

				if(nva01nofix.gt.0) then 

				u1=-gl01*vet01*(su01**vet01)
				u1=u1*(su02**vet02)*ri02*vet02
				u1=u1+(su12**vet12)*ri12*vet12*res201num

				
				res1(1:nva01nofix)=&
				ve01nofix(i,:)*u1/v
				res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01

				res1((nvamax+1):nvamax01)=&
				res20101num*(su12**vet12)*ri12*vet12
				resint=(su01**vet01)*(su02**vet02)*gl01*vet01
				resint=resint*(1-gl01*vet01)*ri02*vet02
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-resint
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)*ve01square(i,:)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)/v
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)+tronc01square


				end if 

				if(nva02nofix.gt.0) then 

				u2=-(su12**vet12)*ri12*vet12*res202num+&
				(1-gl02*vet02)*(su01**vet01)*(su02**vet02)*ri02*vet02

				res1((nva01nofix+1):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1):nva0102)=&
				res1((nva01nofix+1):nva0102)+&
				tronc02

				res1((nvamax0112+1):nvamax02)=&
				-gl02*vet02*(su01**vet01)*(su02**vet02)*ri02*vet02
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)+&
				(su12**vet12)*ri12*vet12*(res20202num-res202num)+&
				((1-gl02*vet02)**2)*ri02*vet02*(su01**vet01)*(su02**vet02)


				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)*ve02square(i,:)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)/v
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)+tronc02square


				end if 

				if(nva12nofix.gt.0) then 

				u3=res212num+(1-gl12*vet12)*res2denum
				u3=u3*(su12**vet12)*ri12*vet12

				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v

				res1((nvamax0212+1):nvamax12)=&
				-gl12*vet12*res2denum+((1-gl12*vet12)**2)*res2denum+&
				(1-gl12*vet12)*2*res212num+res212num+&
				res21212num

				resint=ri12*vet12*(su12**vet12)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)*resint*ve12square(i,:)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)/v
				
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)


				end if 

				
				if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 

				kfix=nvamax01+1
				lfix=kfix-1+nva02nofix

				do j=1,nva01nofix
			    		res1(kfix:lfix)=(-res202num+&
					res20102num)*(su12**vet12)*ri12*vet12
					resint=-gl01*vet01*(1-gl02*vet02)*(su01**vet01)
					resint=resint*(su02**vet02)*ri02*vet02
					res1(kfix:lfix)=&
					res1(kfix:lfix)+resint
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			    		res1(kfix:lfix)=&
			    		res1(kfix:lfix)*v
					res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   	kfix=lfix+1
			   	lfix=lfix+nva02nofix
				end do

				end if 

				if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 

				
				kfix=nvamax0102+1
				lfix=kfix-1+nva12nofix

				do j=1,nva01nofix
			  		res1(kfix:lfix)=&
					(1-gl12*vet12)*res201num-res20112num+&
					res212num
					res1(kfix:lfix)=&
					res1(kfix:lfix)*ve01nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
					res1(kfix:lfix)*(su12**vet12)*ri12*vet12
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)*v
					res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva12nofix
				end do

				end if 


				if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 

				
				
        			kfix=nvamax02+1
				lfix=kfix-1+nva12nofix

				do j=1,nva02nofix
			  	res1(kfix:lfix)=res20212num+&
				(1-gl12*vet12)*res202num
				resint=(su12**vet12)
				res1(kfix:lfix)=&
				res1(kfix:lfix)*resint*ve02nofix(i,j)*ve12nofix(i,:)
				res1(kfix:lfix)=&
				-res1(kfix:lfix)*ri12*vet12
				res1(kfix:lfix)=&
				res1(kfix:lfix)*v
				res1(kfix:lfix)=&
				res1(kfix:lfix)-u2*u3*ve02nofix(i,j)*ve12nofix(i,:)
				res1(kfix:lfix)=&
			   	res1(kfix:lfix)/(v**2)
			   	kfix=lfix+1
			  	lfix=lfix+nva12nofix
				end do

				end if 
                            endif
                         endif                        
                      endif
                   endif  
                endif   
                endif   

                res = res + res1 


		


        end do   


        likelihood_deriv = res




123     continue 
	 
	deallocate(b,bfix,fix,ve01,ve02,ve12,ve01nofix,&
	ve02nofix,ve12nofix,tronc01,tronc02,t0,t1,t2,t3,c,&
	ve01square,ve02square,ve12square,&
	tronc01square,tronc02square)     

    end subroutine derivaweib
	
	
	!=============================================================================================  
!======================= Calculate derivatives of loglik with weibull baseline risk ==========
!=============================================================================================  


subroutine derivaweiballpara(b0,np0,npar0,bfix0,fix0,c0,no0,ve010,ve120,ve020,&
        dimnva01,dimnva12,dimnva02,nva01,nva12,nva02,t00,&
        t10,t20,t30,troncature0,likelihood_deriv)
	
		use commun
        implicit none
         
        double precision::res2denum,res201num,res202num,res212num, &
	res20101num,res20102num,res20112num,res20202num,res20212num,&
	res21212num,res20101numbis,&
	res2thenum,res2thenumsquare,res2the01,res2the02,res2the12,&
	res2the0101,res2the0202,res2the1212,res2the0102,res2the0112,&
	res2the0212,res2the0101square,res2the0202square,&
	res2the1212square,res2the0102square,res2the0112square,&
	res2the0212square,res2the0101dsquare,res2the0202dsquare,&
	res2the1212dsquare,vet01,vet12,vet02,resint,v,u1,u2,u3
	
        integer::np0,i,j,l,w,k,lfix, kfix,npar0,nva01,nva12,nva02,no0, &
	nz010,nz020,nz120,troncature0,dimnva01,dimnva02,dimnva12, & 
	nva01nofix,nva12nofix,nva02nofix,nvamax, sizespline,nva0102
	integer::nvamax01,nvamax0102,nvamax0112,nvamax02,nvamax0212,nvamax12, &
	nvaweib01,nvaweib02,nvaweib12,nvaweib,nvamax12weib12,iter,nweib

	double precision,dimension(np0+np0*(np0+1)/2),intent(inout)::likelihood_deriv
	double precision,dimension(np0)::b0
	double precision,dimension(np0+np0*(np0+1)/2)::res,res1
        double precision,dimension(npar0)::bh
	double precision,dimension(npar0-np0)::bfix0
	integer,dimension(npar0)::fix0
	double precision,dimension(2)::the01,the12,the02
        double precision,dimension(no0,dimnva01)::ve010
	double precision,dimension(no0,dimnva02)::ve020
	double precision,dimension(no0,dimnva12)::ve120
	
	
        double precision::su01,ri01,su12,ri12,su02,ri02,gl01,gl02,gl12,&
		troncweib01011,troncweib01011square, &
		troncweib010112,troncweib01012,&
		troncweib01012square,troncweib02021,troncweib02021square,&
		troncweib020212,troncweib02022,&
		troncweib02022square
	double precision,dimension(no0)::t00,t10,t20,t30
	integer,dimension(no0)::c0

	allocate(b(np0),bfix(npar0-np0),fix(npar0))
	b=b0
	bfix=bfix0
	fix=fix0
	troncature=troncature0

	sizespline=6


	nvaweib01=2-sum(fix(1:2))
	nvaweib02=2-sum(fix(3:4))
	nvaweib12=2-sum(fix(5:6))
	
	!write(6,*) "nvaweib01",nvaweib01 OK
	!write(6,*) "nvaweib02",nvaweib02 OK
	!write(6,*) "nvaweib12",nvaweib12 OK
	
	if(nva01.gt.0) then 
	  nva01nofix=nva01-sum(fix((sizespline+1):(nva01+sizespline)))
	else 
	  nva01nofix=0
	end if 

	if(nva02.gt.0) then 
          nva02nofix=nva02-sum(fix((nva01+sizespline+1):(nva02+nva01+sizespline)))
	else 
	   nva02nofix=0
	end if 

	if(nva12.gt.0) then 
	  nva12nofix=nva12-sum(fix((nva01+nva02+sizespline+1):npar0))
	else 
	  nva12nofix=0
	end if 

	nvaweib=nvaweib01+nvaweib02+nvaweib12
	nva0102=nva01nofix+nva02nofix+nvaweib
	nvamax=nva01nofix+nva02nofix+nva12nofix+nvaweib
	
	nvamax12weib12=nvamax+&
	(nvaweib+1)*nvaweib/2+&
	(nvamax-nvaweib)*nvaweib
	
	nvamax01=nvamax12weib12+(nva01nofix+1)*nva01nofix/2
	nvamax0102=nvamax01+nva01nofix*nva02nofix
	nvamax0112=nvamax0102+nva01nofix*nva12nofix
	
	nvamax02=nvamax0112+(nva02nofix+1)*nva02nofix/2
	nvamax0212=nvamax02+nva02nofix*nva12nofix
	nvamax12=nvamax0212+(nva12nofix+1)*nva12nofix/2


	if(nva01.gt.0) then 
		allocate(ve01(no0,nva01))
		allocate(ve01nofix(no0,nva01nofix))
		allocate(ve01square(no0,nva01nofix*(nva01nofix+1)/2))
		allocate(tronc01(nva01nofix))
		allocate(tronc01square(nva01nofix*(nva01nofix+1)/2))
		allocate(troncweib01011beta01(nva01nofix))
		allocate(troncweib01012beta01(nva01nofix))
	else 
		allocate(ve01(no0,1))
		allocate(ve01nofix(no0,1))
		ve01nofix=0
		allocate(ve01square(no0,1))
		ve01square=0
		allocate(tronc01(1))
		allocate(tronc01square(1))
		allocate(troncweib01011beta01(1))
		allocate(troncweib01012beta01(1))
	end if 
	
	if(nva02.gt.0) then 
		allocate(ve02(no0,nva02))
		allocate(ve02nofix(no0,nva02nofix))
		allocate(ve02square(no0,nva02nofix*(nva02nofix+1)/2))
		allocate(tronc02(nva02nofix))
		allocate(tronc02square(nva02nofix*(nva02nofix+1)/2))
		allocate(troncweib02021beta02(nva02nofix))
		allocate(troncweib02022beta02(nva02nofix))
	else 
		allocate(ve02(no0,1))
		allocate(ve02nofix(no0,1))
		ve02nofix=0
		allocate(ve02square(no0,1))
		ve02square=0
		allocate(tronc02(1))
		allocate(tronc02square(1))
		allocate(troncweib02021beta02(1))
		allocate(troncweib02022beta02(1))
	end if 

	if(nva12.gt.0) then 
		allocate(ve12(no0,nva12))
		allocate(ve12nofix(no0,nva12nofix))
		allocate(ve12square(no0,nva12nofix*(nva12nofix+1)/2))
	else 
		allocate(ve12(no0,1))
		allocate(ve12nofix(no0,1))
		ve12nofix=0
		allocate(ve12square(no0,1))
		ve12square=0
	end if 



	ve01=ve010
	ve02=ve020
	ve12=ve120

	allocate(t0(no0),t1(no0),t2(no0),t3(no0),c(no0))
	c=c0
	t0=t00
	t1=t10
	t2=t20
	t3=t30

         
        ! we need to put bh at its original values if in posfix 


       l=0
       lfix=0
       w=0

	if(nva01.gt.0) then
	do k=1,nva01
	   if(fix((sizespline+k)).eq.0) then 
		lfix=lfix+1
		ve01nofix(:,lfix)=ve01(:,k)
	   end if 
	end do
	
	do i=1,no0
		lfix=1
		do l=1,nva01nofix
	   	ve01square(i,lfix:(lfix+nva01nofix-l))=&
		ve01nofix(i,l:nva01nofix)
		ve01square(i,lfix:(lfix+nva01nofix-l))=&
		ve01square(i,lfix:(lfix+nva01nofix-l))*ve01nofix(i,l)
		lfix=lfix+nva01nofix-l+1
		end do
	end do
	 
        end if 
	lfix=0

	if(nva02.gt.0) then 
	do k=1,nva02
	   if(fix((sizespline+nva01+k)).eq.0) then 
		lfix=lfix+1
		ve02nofix(:,lfix)=ve02(:,k)
	   end if 
	end do

	do i=1,no0
		lfix=1
		do l=1,nva02nofix
	   	ve02square(i,lfix:(lfix+nva02nofix-l))=&
		ve02nofix(i,l)*ve02nofix(i,l:nva02nofix)
		lfix=lfix+nva02nofix-l+1
		end do
	end do

	end if 
	
	lfix=0

	if(nva12.gt.0) then 
	do k=1,nva12
	   if(fix((sizespline+nva01+nva02+k)).eq.0) then 
		lfix=lfix+1
		ve12nofix(:,lfix)=ve12(:,k)
	   end if 
	end do

	do i=1,no0
		lfix=1
		do l=1,nva12nofix
	   		ve12square(i,lfix:(lfix+nva12nofix-l))=&
			ve12nofix(i,l)*ve12nofix(i,l:nva12nofix)
			lfix=lfix+nva12nofix-l+1
		end do
	end do

	end if
	
	l=0
       lfix=0
       w=0

     do k=1,npar0 
         if(fix(k).eq.0) then
            l=l+1
            bh(k)=b(l)
	 end if 
         if(fix(k).eq.1) then
            w=w+1
            bh(k)=bfix(w)
         end if
      end do

   
	


	
  
    !the01(2)=dexp(bh(2))
	!the02(2)=dexp(bh(4))
	!the12(2)=dexp(bh(6))
    
	the01(2)=(bh(2))**2
	the02(2)=(bh(4))**2
	the12(2)=(bh(6))**2
	
	the01(1)=(bh(1))**2
	the02(1)=(bh(3))**2
	the12(1)=(bh(5))**2

	!write(6, *) "before loop" OK

!---------- calcul des derivees premiere ------------------   

	res = 0
        do i=1,no0

!write(6, *) "subject ",i
                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
				vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
				vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
				vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

		
                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)

                res1 = 0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                            tronc01 = 0
							tronc02 =  0
                        	tronc01square= 0
                        	tronc02square=0
							troncweib01011=0
							troncweib01011square=0
							troncweib01011beta01=0
							troncweib01012beta01=0
							troncweib010112=0
							troncweib01012=0
							troncweib01012square=0
							troncweib02021=0
							troncweib02021square=0
							troncweib020212=0
							troncweib02021beta02=0
							troncweib02022=0
							troncweib02022square=0
							troncweib02022beta02=0
                        else 
				call fonct(t0(i),the01,ri01,gl01,su01)
				call fonct(t0(i),the02,ri02,gl02,su02)
				
                            tronc01=ve01nofix(i,:)*gl01*vet01
                        	tronc02=ve02nofix(i,:)*gl02*vet02
                        	tronc01square=ve01square(i,:)*gl01*vet01
                        	tronc02square=ve02square(i,:)*gl02*vet02
							
							troncweib01011=gl01*vet01*LOG(the01(2)*t0(i))
							troncweib01011square=&
							((LOG(the01(2)*t0(i)))**2)*gl01*vet01
							
							troncweib01011beta01=LOG(the01(2)*t0(i))*tronc01
							
							troncweib01012beta01=tronc01*the01(1)/the01(2)
							
							troncweib010112=&
							(1+LOG(the01(2)*t0(i))*the01(1))*gl01*vet01/the01(2)
							
							troncweib01012=the01(1)*gl01*vet01/the01(2)
							
							troncweib01012square=&
							the01(1)*(the01(1)-1)*gl01*vet01/(the01(2)**2)
							
							troncweib02021=LOG(the02(2)*t0(i))*gl02*vet02
							troncweib02021square=&
							((LOG(the02(2)*t0(i)))**2)*gl02*vet02
							
							troncweib020212=&
							(1+LOG(the02(2)*t0(i))*the02(1))*gl02*vet02/the02(2)
							
							troncweib02021beta02=&
							tronc02*LOG(the02(2)*t0(i))
							
							troncweib02022=&
							the02(1)*gl02*vet02/the02(2)
							
							troncweib02022square=&
							the02(1)*(the02(1)-1)*gl02*vet02/(the02(2)**2)
							
							troncweib02022beta02=&
							tronc02*the02(1)/the02(2)
							
							
							
                        end if
                else
                    tronc01 = 0
							tronc02 =  0
                        	tronc01square= 0
                        	tronc02square=0
							troncweib01011=0
							troncweib01011square=0
							troncweib01011beta01=0
							troncweib01012beta01=0
							troncweib010112=0
							troncweib01012=0
							troncweib01012square=0
							troncweib02021=0
							troncweib02021square=0
							troncweib020212=0
							troncweib02021beta02=0
							troncweib02022=0
							troncweib02022square=0
							troncweib02022beta02=0
                end if
		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
				!write(6, *) "start c=1" 
			call fonct(t1(i),the01,ri01,gl01,su01)
			call fonct(t1(i),the02,ri02,gl02,su02)
			iter = 0
			nweib= 0
			if(fix(1).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				res1(nweib)=-LOG(the01(2)*t1(i))*gl01*vet01 +&
				troncweib01011
				res1((nvamax+iter))=&
				-((LOG(the01(2)*t1(i)))**2)*gl01*vet01 +&
				troncweib01011square
				if(fix(2).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=&
					-(1+the01(1)*LOG(the01(2)*t1(i)))* &
					gl01*vet01/the01(2) +&
					troncweib010112
				endif
				if(fix(3).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					-ve01nofix(i,:)*gl01*vet01*LOG(the01(2)*t1(i))+&
					troncweib01011beta01
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif 
			
			if(fix(2).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				res1(nweib)=-the01(1)*gl01*vet01/the01(2) +&
				troncweib01012
				
				res1(nvamax+iter)=&
				-the01(1)*(the01(1)-1)*gl01*vet01/(the01(2)**2)+&
				troncweib01012square
				
				if(fix(3).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					-ve01nofix(i,:)*gl01*vet01*the01(1)/the01(2)+&
					troncweib01012beta01
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif
			
			if(fix(3).eq.0)then
				iter = iter + 1
				nweib = nweib +1
				res1(nweib)=-LOG(the02(2)*t1(i))*gl02*vet02 +&
				troncweib02021
				res1((nvamax+iter))=&
				-((LOG(the02(2)*t1(i)))**2)*gl02*vet02+&
				troncweib02021square
				
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=&
					-(1+the02(1)*LOG(the02(2)*t1(i)))* &
					gl02*vet02/the02(2)+&
					troncweib020212
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					troncweib02021beta02-&
					ve02nofix(i,:)*gl02*vet02*LOG(the02(2)*t1(i))
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif 
			
			if(fix(4).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				res1(nweib)=-the02(1)*gl02*vet02/the02(2)+&
				troncweib02022
				res1(nvamax+iter)=&
				-the02(1)*(the02(1)-1)*gl02*vet02/(the02(2)**2)+&
				troncweib02022square
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					-ve02nofix(i,:)*gl02*vet02*the02(1)/the02(2)+&
					troncweib02022beta02
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif
			
			
			if(fix(5).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				res1(nweib)=0
				res1(nvamax+iter)=0
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif
			
			if(fix(6).eq.0)then
				iter = iter + 1
				nweib = nweib +1
				res1(nweib)=0
				res1(nvamax+iter)=0
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif

			if(nva01nofix.gt.0) then 
			
			res1((nvaweib+1):(nva01nofix+nvaweib))=&
			-ve01nofix(i,:)*gl01*vet01+tronc01
			res1((nvamax12weib12+1):nvamax01)=&
			-ve01square(i,:)*gl01*vet01+&
			tronc01square

			if(nva02nofix.gt.0) then 
			res1((nvamax01+1):nvamax0102)=0
			end if 
			if(nva12nofix.gt.0) then 
			res1((nvamax0102+1):nvamax0112)=0
			end if

			end if 

			if(nva02nofix.gt.0) then
			res1((nvaweib+nva01nofix+1):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02

			res1((nvamax0112+1):nvamax02)=&
			-ve02square(i,:)*gl02*vet02+&
			tronc02square
			if(nva12nofix.gt.0) then 
			res1((nvamax02+1):nvamax0212)=0
			end if 
			end if 

			if(nva12nofix.gt.0) then 
			res1((nva0102+1):nvamax)=0
			res1((nvamax0212+1):nvamax12)=0
			end if 
			
			
			

                else
                if(c(i).eq.2)then ! cpi 0-->1
			!write(6, *) "start c=2" 

			call fonct(t3(i),the12,ri12,gl12,su12)
			!call logapproxi(t1(i),t2(i),res2denum) approxi log(x) OK
			call qgaussweibbetaderiv(t1(i),t2(i),the01,&
			the02,the12,res2denum,res2the01,res2the02,&
			res2the12,res2thenum,res2thenumsquare,&
			res2the0101,res2the0102,res2the0112,res2the0202,&
			res2the0212,res2the1212,res2the0101square,&
			res2the0102square,res2the0112square,&
			res2the0202square,res2the0212square,&
			res2the1212square,res2the0101dsquare,&
			res2the0202dsquare,res2the1212dsquare,&
			res201num,res202num,res212num,res20101num,&
			res20101numbis,&
			res20102num,res20112num,&
			res20202num,res20212num,res21212num,&
			vet01,vet02,vet12)
			
			iter = 0
			nweib = 0
			!write(6,*) "res2the01",res2the01
			if(fix(1).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				!write(6,*) "nweib",nweib
				!write(6,*) "res2denum theta",res2denum
				v=(LOG(the01(2))+(1/the01(1)))*res2denum +&
				res2thenum-LOG(the01(2))*(res2denum-res201num)-&
				res2the01
				res1(nweib)=v/res2denum
				
				res1(nweib)=res1(nweib)+troncweib01011
				
				!write(6,*) "troncweib01011square",troncweib01011square
				res1((nvamax+iter))=&
				(LOG(the01(2))**2)*res2denum+ &
				LOG(the01(2))*2*res2thenum+ &
				res2thenumsquare+ &
				LOG(the01(2))*2/the01(1)*res2denum+ &
				res2thenum*2/the01(1)- &
				(LOG(the01(2))**2)*3*(res2denum-res201num)- &
				LOG(the01(2))*6*res2the01- &
				res2the0101square*3- &
				LOG(the01(2))*2*(res2denum-res201num)/the01(1)-&
				res2the01*2/the01(1)+ &
				(LOG(the01(2))**2)*res20101numbis+ &
				LOG(the01(2))*2*res2the0101+ &
				res2the0101dsquare
				
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*res2denum -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib01011square
				
				if(fix(2).eq.0)then
					
					u2=res2denum*the01(1)/the01(2)-&
					(res2denum-res201num)*the01(1)/the01(2)
					
					u1=the01(1)*LOG(the01(2))*res20101numbis/the01(2)+&
					the01(1)*res2the0101/the01(2)+&
					(the01(1)*LOG(the01(2))+2)*res2denum/the01(2)+&
					the01(1)*res2thenum/the01(2)-&
					(3*the01(1)*LOG(the01(2))+2)* &
					(res2denum-res201num)/the01(2)-3* &
					the01(1)*res2the01/the01(2)
					
					iter = iter +1
					res1((nvamax+iter))=&
					u1*res2denum-u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+troncweib010112
				endif
				if(fix(3).eq.0)then
				
					u1=-LOG(the01(2))*LOG(the02(2))*res202num-&
					LOG(the01(2))*res2the02-LOG(the02(2))*res2the02-&
					res2the0202square-LOG(the02(2))*res202num/the01(1)-&
					res2the02/the01(1)+&
					LOG(the01(2))*LOG(the02(2))*res20102num+&
					LOG(the01(2))*res2the0102+&
					LOG(the02(2))*res2the0102+&
					res2the0102square
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(fix(4).eq.0)then
				
					u1=-LOG(the01(2))*the02(1)*res202num/the02(2)-&
					the02(1)*res2the02/the02(2)-&
					the02(1)*res202num/(the01(1)*the02(2))+&
					LOG(the01(2))*the02(1)*res20102num/the02(2)+&
					the02(1)*res2the0102/the02(2)
					
					u2=-the02(1)*res202num/the02(2)
					
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(fix(5).eq.0)then
				
					u1=LOG(the01(2))*LOG(the12(2))*res212num+&
					LOG(the01(2))*res2the12+LOG(the12(2))*res2the12+&
					res2the1212square+LOG(the12(2))*res212num/the01(1)+&
					res2the12/the01(1)-&
					LOG(the01(2))*LOG(the12(2))*res20112num-&
					LOG(the01(2))*res2the0112-&
					LOG(the12(2))*res2the0112-&
					res2the0112square
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
					
					u1=&
					LOG(the01(2))*the12(1)*res212num/the12(2)+&
					the12(1)/the12(2)*res2the12+res212num* &
					the12(1)/(the12(2)*the01(1))-&
					LOG(the01(2))*the12(1)*res20112num/the12(2)-&
					the12(1)/the12(2)*res2the0112
					
					
					u2=the12(1)*res212num/the12(2)
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(nva01nofix.gt.0) then 
					
					u1=res2the0101+ &
					LOG(the01(2))*res20101numbis-&
					LOG(the01(2))*3*(res2denum-res201num)-&
					res2the01*3-(res2denum-res201num)/the01(1)+&
					(LOG(the01(2))+(1/the01(1)))*res2denum+&
					res2thenum
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01011beta01
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					u1=-(LOG(the01(2))+&
					(1/the01(1)))*res202num-&
					res2the02+&
					LOG(the01(2))*res20102num+&
					res2the0102
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					u1=(LOG(the01(2))+&
					(1/the01(1)))*res212num+&
					res2the12-&
					LOG(the01(2))*res20112num-&
					res2the0112
					
					u2=res212num
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva12nofix
				endif
			endif
			
			!write(6, *) " done x101 -c2" 
			if(fix(2).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=the01(1)*res2denum/the01(2) -&
				the01(1)*(res2denum-res201num)/the01(2)
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+troncweib01012
				
				u1=((the01(1)/the01(2))**2)*res20101numbis+&
					the01(1)*(the01(1)-1)*res2denum/((the01(2))**2)-&
					(res2denum-res201num)*(3*(the01(1)**2)-&
					the01(1))/(the01(2)**2)
				
				
				res1((nvamax+iter))=&
				u1*res2denum-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+ troncweib01012square
				if(fix(3).eq.0)then
				
					u1=-LOG(the02(2))*the01(1)* &
					res202num/the01(2)-&
					the01(1)*res2the02/the01(2)+&
					the01(1)*LOG(the02(2))*res20102num/the01(2)+&
					the01(1)*res2the0102/the01(2)
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(fix(4).eq.0)then
				
					u1=the01(1)*the02(1)*res20102num-&
					the01(1)*the02(1)*res202num
					u1=u1/(the01(2)*the02(2))
					
					u2=-the02(1)*res202num/the02(2)
				
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(fix(5).eq.0)then
				
					u1=LOG(the12(2))*the01(1)*res212num/the01(2)+&
					the01(1)*res2the12/the01(2)-&
					the01(1)*LOG(the12(2))*res20112num/the01(2)-&
					the01(1)*res2the0112/the01(2)
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
					
					u1=-the01(1)*the12(1)*res20112num+&
					the01(1)*the12(1)*res212num
					u1=u1/(the01(2)*the12(2))
					
					u2=the12(1)*res212num/the12(2)
				
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(nva01nofix.gt.0) then 
					
					u1=the01(1)*res20101numbis/the01(2)+&
					the01(1)*res2denum/the01(2)-&
					the01(1)*3* &
					(res2denum-res201num)/the01(2)
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01012beta01
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					u1=res20102num-res202num
					u1=u1*the01(1)/the01(2)
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!here 
					u1=-res20112num+res212num
					u1=u1*the01(1)/the01(2)
					
					u2=res212num
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva12nofix
				endif
			endif
			
				!write(6, *) " done x201 -c2" 	
				!here 
				
			if(fix(3).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				!write(6,*) "res2the02",res2the02
				v=-LOG(the02(2))*res202num-&
				res2the02
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+ &
				troncweib02021
				
				res1((nvamax+iter))=&
				((LOG(the02(2)))**2)* &
				res20202num+2*LOG(the02(2))*res2the0202+&
				res2the0202dsquare-&
				((LOG(the02(2)))**2)* &
				res202num-2*LOG(the02(2))*res2the02-&
				res2the0202square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*res2denum -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				res1((nvamax+iter))=&
				res1((nvamax+iter))+&
				troncweib02021square
				
				if(fix(4).eq.0)then
				
					u1=-(the02(1)*LOG(the02(2))+1)* &
					res202num/the02(2)-&
					the02(1)*res2the02/the02(2)+ &
					LOG(the02(2))*the02(1)* &
					res20202num/the02(2)+ &
					the02(1)/the02(2)*res2the0202
					
					u2=-the02(1)*res202num/the02(2)
					
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+&
					troncweib020212
				endif
				
				if(fix(5).eq.0)then
				
					u1=-LOG(the02(2))*LOG(the12(2))*res20212num-&
					LOG(the02(2))*res2the0212-LOG(the12(2))*res2the0212-&
					res2the0212square
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
				
					u1=-LOG(the02(2))*the12(1)*res20212num/the12(2)-&
					the12(1)*res2the0212/the12(2)
					
					u2=the12(1)*res212num/the12(2)
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(nva01nofix.gt.0) then 
					
					u1=LOG(the02(2))*res20102num+&
					res2the0102-res202num*LOG(the02(2))-&
					res2the02
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					u1=-LOG(the02(2))*res202num-&
					res2the02+LOG(the02(2))*res20202num+&
					res2the0202
					
					u2=-res202num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02021beta02
					iter = iter +nva02nofix
					
				endif
				if(nva12nofix.gt.0) then 
				
					u1=-LOG(the02(2))*res20212num-&
					res2the0212
					
					u2=res212num
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva12nofix
				endif
			endif
			
             
			!write(6, *) " done x102 -c2" 	
			if(fix(4).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=-the02(1)*res202num/the02(2)
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+troncweib02022
				
				u1=((the02(1)/the02(2))**2)*res20202num-&
					(the02(1)*(the02(1)-1)/((the02(2))**2))*res202num
				
				
				res1((nvamax+iter))=&
				u1*res2denum-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib02022square
				
				
				if(fix(5).eq.0)then
				
					u1=-the02(1)*LOG(the12(2))*res20212num/the02(2)-&
					the02(1)*res2the0212/the02(2)
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
					
					u1=-the02(1)*the12(1)*res20212num
					u1=u1/(the02(2)*the12(2))
					
					u2=the12(1)*res212num/the12(2)
				
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				
				if(nva01nofix.gt.0) then 
					
					u1=-the02(1)*res202num/the02(2)+&
					the02(1)*res20102num/the02(2)					
					
					u2=res201num
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				
				if(nva02nofix.gt.0) then 
					
					u1=the02(1)*res20202num/the02(2)-&
					the02(1)*res202num/the02(2)
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02022beta02
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					u1=-res20212num*the02(1)/the02(2)
					
					u2=res212num
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva12nofix
				endif
			endif	
			!write(6, *) " done x202 -c2" 	
			if(fix(5).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=LOG(the12(2))*res212num+&
				res2the12
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)-&
				LOG(the12(2)*t3(i))*gl12*vet12
				
				res1((nvamax+iter))=&
				((LOG(the12(2)))**2)* &
				res21212num+2*LOG(the12(2))*res2the1212+&
				res2the1212dsquare+&
				((LOG(the12(2)))**2)* &
				res212num+2*LOG(the12(2))*res2the12+&
				res2the1212square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*res2denum -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				res1((nvamax+iter))=&
				res1((nvamax+iter))-&
				((LOG(the12(2)*t3(i)))**2)*gl12*vet12
				

				if(fix(6).eq.0)then
			
					u1=LOG(the12(2))*the12(1)*res21212num/the12(2)+&
					the12(1)*res2the1212/the12(2)+&
					(the12(1)*LOG(the12(2))+1)*res212num/the12(2)+&
					the12(1)*res2the12/the12(2)
					
					u2=the12(1)*res212num/the12(2)
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					res1((nvamax+iter))=&
					res1((nvamax+iter))-&
					(1+LOG(the12(2)*t3(i))*the12(1))*gl12*vet12/the12(2)
				endif
				
				if(nva01nofix.gt.0) then 
					
					u1=-LOG(the12(2))*res20112num-&
					res2the0112+res212num*LOG(the12(2))+&
					res2the12
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					!here
					u1=-LOG(the12(2))*res20212num-&
					res2the0212
					
					u2=-res202num
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva02nofix
					
				endif
				if(nva12nofix.gt.0) then 
				
					u1=LOG(the12(2))*res21212num+&
					res2the1212+LOG(the12(2))*res212num+&
					res2the12
					
					u2=res212num
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*LOG(the12(2)*t3(i))*gl12*vet12
					
					iter = iter +nva12nofix
				endif
			endif
			
			  !!write(6, *) " done x112 -c2" 	

			if(fix(6).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=the12(1)*res212num/the12(2)
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)-&
				the12(1)*gl12*vet12/the12(2)
				
				u1=((the12(1)/the12(2))**2)*res21212num+&
					(the12(1)*(the12(1)-1)/((the12(2))**2))*res212num
				
				
				res1((nvamax+iter))=&
				u1*res2denum-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))-the12(1)* &
				(the12(1)-1)*gl12*vet12/(the12(2)**2)
				
				
				
				if(nva01nofix.gt.0) then 
					
					u1=the12(1)*res212num/the12(2)-&
					the12(1)*res20112num/the12(2)					
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				
				if(nva02nofix.gt.0) then 
					
					u1=-the12(1)*res20212num/the12(2)
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!here
					u1=res21212num*the12(1)/the12(2)+&
					the12(1)*res212num/the12(2)
					
					u2=res212num
							
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*the12(1)*gl12*vet12/the12(2)
					iter = iter +nva12nofix
				endif
			endif	
			
			write(6,*) "res1(1) after all",res1(1) 	
			v=res2denum*(su12**vet12)
			!write(6,*)"v",v
			!write(6,*)"nvaweib",nvaweib
			!write(6,*) "nvamax",nvamax
			!write(6,*) "nvamax01",nvamax01
			!write(6,*) "nvamax02",nvamax02
			!write(6,*) "nvamax0112",nvamax0112
			!write(6,*) "nvamax0212",nvamax0212
			!write(6,*) "nvamax12",nvamax12
			!write(6,*) "nva0102",nva0102
			!write(6,*) "iter",iter
			!write(6,*) "nvamax12weib12",nvamax12weib12
			
			if(nva01nofix.gt.0) then

			u1=res201num*(su12**vet12)
			
      		res1((nvaweib+1):(nvaweib+nva01nofix))=&
			ve01nofix(i,:)*u1/v
			res1((nvaweib+1):(nvaweib+nva01nofix))=&
			res1((nvaweib+1):(nvaweib+nva01nofix))+tronc01

			res1((nvamax12weib12+1):nvamax01)=&
			ve01square(i,:)*res20101num*(su12**vet12)
			res1((nvamax12weib12+1):nvamax01)=&
			res1((nvamax12weib12+1):nvamax01)/v
			res1((nvamax12weib12+1):nvamax01)=&
			res1((nvamax12weib12+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
			res1((nvamax12weib12+1):nvamax01)=&
			res1((nvamax12weib12+1):nvamax01)+&
			tronc01square

			end if 

			
			if(nva02nofix.gt.0) then

			u2=-res202num*(su12**vet12)
			res1((nvaweib+nva01nofix+1):nva0102)=&
			ve02nofix(i,:)*u2/v
			res1((nvaweib+nva01nofix+1):nva0102)=&
			res1((nvaweib+nva01nofix+1):nva0102)+&
			tronc02

			res1((nvamax0112+1):nvamax02)=&
			-ve02square(i,:)*(su12**vet12)*(-res20202num+res202num)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)/v
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)+tronc02square


			end if 

			
			if(nva12nofix.gt.0) then 

			u3=res212num-gl12*vet12*res2denum
			u3=u3*(su12**vet12)

			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v
			res1((nvamax0212+1):nvamax12)=&
			-2*gl12*vet12*res212num+res212num+&
			res21212num-gl12*vet12*(1-gl12*vet12)*res2denum
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*(su12**vet12)*ve12square(i,:)
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)/v
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)

			
			end if 
			
			

			
			if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 

			kfix=nvamax01+1
			lfix=kfix-1+nva02nofix

			do j=1,nva01nofix

			   res1(kfix:lfix)=&
			  (res20102num-res202num)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva02nofix

			end do

			end if 

			if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 

			kfix=nvamax0102+1
			lfix=kfix-1+nva12nofix

			do j=1,nva01nofix
			   res1(kfix:lfix)=-gl12*vet12*res201num-res20112num+&
			   res212num
			   res1(kfix:lfix)=res1(kfix:lfix)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do

			end if 

			
			if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 

			kfix=nvamax02+1
			lfix=kfix-1+nva12nofix

			do j=1,nva02nofix
			   res1(kfix:lfix)=-res20212num+&
			   gl12*vet12*res202num
			   res1(kfix:lfix)=res1(kfix:lfix)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve02nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u2*u3*ve02nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do

			end if 
			
                else  
                    if(c(i).eq.3)then ! obs 0-->1
					!here 
			!write(6, *) "start c=3" 
			call fonct(t1(i),the01,ri01,gl01,su01)
			call fonct(t1(i),the02,ri02,gl02,su02)
			call fonct(t1(i),the12,ri12,gl12,su12)

			iter = 0
			nweib= 0
			if(fix(1).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				res1(nweib)=-LOG(the01(2)*t1(i))*gl01*vet01 +&
				troncweib01011+&
				LOG(the01(2)*t1(i))+(1/the01(1))
				
				res1((nvamax+iter))=&
				-((LOG(the01(2)*t1(i)))**2)*gl01*vet01 +&
				troncweib01011square-1/(the01(1)**2)
				!here
				if(fix(2).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=&
					-(LOG(the01(2)*t1(i))* &
					the01(1)+1)*gl01*vet01/the01(2)+&
					troncweib010112+(1/the01(2))
				endif
				if(fix(3).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				!here
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					-ve01nofix(i,:)*gl01*vet01*LOG(the01(2)*t1(i))+&
					troncweib01011beta01
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif 
			
			if(fix(2).eq.0)then
			
				iter = iter + 1
				nweib = nweib + 1
				
				res1(nweib)=-the01(1)*gl01*vet01/the01(2) +&
				troncweib01012+&
				the01(1)/the01(2)

				res1(nvamax+iter)=&
				-the01(1)*(the01(1)-1)*gl01*vet01/(the01(2)**2)+&
				troncweib01012square-&
				the01(1)/(the01(2)**2)
				
				if(fix(3).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
				
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					-ve01nofix(i,:)*gl01*vet01*the01(1)/the01(2)+&
					troncweib01012beta01
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif
			!here
			if(fix(3).eq.0)then
				iter = iter + 1
				nweib = nweib +1
				res1(nweib)=-LOG(the02(2)*t1(i))*gl02*vet02 +&
				troncweib02021
				
				res1((nvamax+iter))=&
				-((LOG(the02(2)*t1(i)))**2)*gl02*vet02+&
				troncweib02021square
				!here
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=&
					-(LOG(the02(2)*t1(i))* &
					the02(1)+1)*gl02*vet02/the02(2)+&
					troncweib020212
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter):(nvamax+iter+nva02nofix))=&
					-ve02nofix(i,:)*LOG(the02(2)*t1(i))*gl12*vet12+&
					troncweib02021beta02
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif 
			!here
			if(fix(4).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				
				res1(nweib)=-the02(1)*gl02*vet02/the02(2)+&
				troncweib02022
				
				res1(nvamax+iter)=&
				-the02(1)*(the02(1)-1)*gl02*vet02/(the02(2)**2)+&
				troncweib02022square
				!here
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					-ve02nofix(i,:)*gl02*vet02*the02(1)/the02(2)+&
					troncweib02022beta02
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif
			!here
			
			if(fix(5).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				res1(nweib)=LOG(the12(2)*t1(i))*gl12*vet12
				res1(nvamax+iter)=((LOG(the12(2)*t1(i)))**2)*gl12*vet12
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=(LOG(the12(2)*t1(i))* &
					the12(1)+1)*gl12*vet12/the12(2)
				endif
				if(nva01nofix.gt.0) then 
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*LOG(the12(2)*t1(i))*gl12*vet12
					call fonct(t3(i),the12,ri12,gl12,su12)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*LOG(the12(2)*t3(i))*gl12*vet12
					
					iter = iter +nva12nofix
				endif
				
				call fonct(t3(i),the12,ri12,gl12,su12)
				
				res1(nweib)=res1(nweib)-&
				LOG(the12(2)*t3(i))*gl12*vet12
				iter=iter-nva12nofix-nva02nofix-nva01nofix
				if(fix(6).eq.0)then
					res1((nvamax+iter))=&
					res1((nvamax+iter))-&
					(LOG(the12(2)*t3(i))* &
					the12(1)+1)*gl12*vet12/the12(2)
					iter=iter-1
				endif
				
				res1(nvamax+iter)=res1(nvamax+iter)-&
				((LOG(the12(2)*t3(i)))**2)*gl12*vet12
				
				iter=iter+nva12nofix+nva02nofix+nva01nofix
				
				if(fix(6).eq.0)then
					iter=iter+1
				endif
				call fonct(t1(i),the12,ri12,gl12,su12)
			endif
			!here
			
			if(fix(6).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				res1(nweib)=the12(1)*gl12*vet12/the12(2)
				res1(nvamax+iter)=the12(1)* &
				(the12(1)-1)*gl12*vet12/(the12(2)**2)
				
				if(nva01nofix.gt.0) then 
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*the12(1)*gl12*vet12/the12(2)
					call fonct(t3(i),the12,ri12,gl12,su12)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*the12(1)*gl12*vet12/the12(2)
					
					iter = iter +nva12nofix
				endif
				
				call fonct(t3(i),the12,ri12,gl12,su12)
				
				res1(nweib)=res1(nweib)-&
				the12(1)*gl12*vet12/the12(2)
				
				iter=iter-nva12nofix-nva02nofix-nva01nofix
				
				res1(nvamax+iter)=res1(nvamax+iter)-&
				the12(1)*(the12(1)-1)*gl12*vet12/(the12(2)**2)
				
				iter=iter+nva12nofix+nva02nofix+nva01nofix
				
				call fonct(t1(i),the12,ri12,gl12,su12)
			endif
			
			
!finish
			if(nva01nofix.gt.0) then 

			res1((nvaweib+1):(nvaweib+nva01nofix))=-ve01nofix(i,:)*gl01*vet01+&
			tronc01+ve01nofix(i,:)
			res1((nvamax12weib12+1):nvamax01)=&
			-ve01square(i,:)*gl01*vet01+&
			tronc01square

			if(nva02nofix.gt.0) then 
			res1((nvamax01+1):nvamax0102)=0
			end if 

			if(nva12nofix.gt.0)then 
			res1((nvamax0102+1):nvamax0112)=0
			end if 

			end if 

			if(nva02nofix.gt.0) then 

			res1((nva01nofix+1+nvaweib):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02
			res1((nvamax0112+1):nvamax02)=&
			-ve02square(i,:)*gl02*vet02+&
			tronc02square

			if(nva12nofix.gt.0) then 
			res1((nvamax02+1):nvamax0212)=0
			end if 

			end if 

			if(nva12nofix.gt.0) then 

			res1((nva0102+1):nvamax)=gl12*vet12
			res1((nvamax0212+1):nvamax12)=gl12*vet12

			call fonct(t3(i),the12,ri12,gl12,su12)
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)-gl12*vet12
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)*ve12nofix(i,:)

			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)-gl12*vet12
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*ve12square(i,:)
			

			end if 
			
			

                    else   !here!17124
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			
			!write(6, *) "start c=4" 
			call fonct(t3(i),the12,ri12,gl12,su12)
			call qgaussweibbetaderiv(t1(i),t2(i),the01,&
			the02,the12,res2denum,res2the01,res2the02,&
			res2the12,res2thenum,res2thenumsquare,&
			res2the0101,res2the0102,res2the0112,res2the0202,&
			res2the0212,res2the1212,res2the0101square,&
			res2the0102square,res2the0112square,&
			res2the0202square,res2the0212square,&
			res2the1212square,res2the0101dsquare,&
			res2the0202dsquare,res2the1212dsquare,&
			res201num,res202num,res212num,res20101num,&
			res20101numbis,&
			res20102num,res20112num,&
			res20202num,res20212num,res21212num,&
			vet01,vet02,vet12)
			
			iter = 0
			nweib = 0
			
			
			if(fix(1).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=(LOG(the01(2))+(1/the01(1)))*res2denum +&
				res2thenum-LOG(the01(2))*(res2denum-res201num)-&
				res2the01
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+troncweib01011
				
				res1((nvamax+iter))=&
				LOG(the01(2))*2/the01(1)*res2denum+&
				res2thenum*2/the01(1)+&
				(LOG(the01(2))**2)*res2denum+&
				LOG(the01(2))*2*res2thenum+&
				res2thenumsquare-3*(LOG(the01(2))**2)* &
				(res2denum-res201num)-&
				LOG(the01(2))*6*res2the01-3*res2the0101square-&
				LOG(the01(2))*2/the01(1)*(res2denum-res201num)- &
				res2the01*2/the01(1)+ &
				(LOG(the01(2))**2)*res20101numbis+&
				LOG(the01(2))*2*res2the0101+&
				res2the0101dsquare
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*res2denum -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+&
				troncweib01011square
				
				if(fix(2).eq.0)then
					
					u2=res2denum*the01(1)/the01(2)-&
					(res2denum-res201num)*the01(1)/the01(2)
					
					u1=the01(1)*LOG(the01(2))*res20101numbis/the01(2)+&
					the01(1)*res2the0101/the01(2)+&
					(the01(1)*LOG(the01(2))+2)*res2denum/the01(2)+&
					the01(1)*res2thenum/the01(2)-&
					(3*the01(1)*LOG(the01(2))+2)* &
					(res2denum-res201num)/the01(2)-3* &
					the01(1)*res2the01/the01(2)
					
					iter = iter +1
					res1((nvamax+iter))=&
					u1*res2denum-u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+&
					troncweib010112
				endif
				if(fix(3).eq.0)then
				!here
					u1=-LOG(the01(2))*LOG(the02(2))*res202num-&
					LOG(the01(2))*res2the02-LOG(the02(2))*res2the02-&
					res2the0202square-LOG(the02(2))*res202num/the01(1)-&
					res2the02/the01(1)+&
					LOG(the01(2))*LOG(the02(2))*res20102num+&
					LOG(the01(2))*res2the0102+&
					LOG(the02(2))*res2the0102+&
					res2the0102square
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(fix(4).eq.0)then
				
					u1=-LOG(the01(2))*the02(1)*res202num/the02(2)-&
					the02(1)*res2the02/the02(2)-&
					the02(1)*res202num/(the01(1)*the02(2))+&
					LOG(the01(2))*the02(1)*res20102num/the02(2)+&
					the02(1)*res2the0102/the02(2)
					
					u2=-the02(1)*res202num/the02(2)
					
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(fix(5).eq.0)then
				!here
					u1=LOG(the01(2))*LOG(the12(2))*res212num+&
					LOG(the01(2))*res2the12+LOG(the12(2))*res2the12+&
					res2the1212square+LOG(the12(2))*res212num/the01(1)+&
					res2the12/the01(1)-&
					LOG(the01(2))*LOG(the12(2))*res20112num-&
					LOG(the01(2))*res2the0112-&
					LOG(the12(2))*res2the0112-&
					res2the0112square
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
					u1=&
					LOG(the01(2))*the12(1)*res212num/the12(2)+&
					the12(1)/the12(2)*res2the12+res212num* &
					the12(1)/(the12(2)*the01(1))-&
					LOG(the01(2))*the12(1)*res20112num/the12(2)-&
					the12(1)/the12(2)*res2the0112
					
					
					u2=the12(1)*res212num/the12(2)
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(nva01nofix.gt.0) then 
					!here
					u1=res2the0101+ &
					LOG(the01(2))*res20101numbis-&
					LOG(the01(2))*3*(res2denum-res201num)-&
					res2the01*3-(res2denum-res201num)/the01(1)+&
					(LOG(the01(2))+(1/the01(1)))*res2denum+&
					res2thenum
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01011beta01
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					u1=-(LOG(the01(2))+&
					(1/the01(1)))*res202num-&
					res2the02+&
					LOG(the01(2))*res20102num+&
					res2the0102
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					u1=(LOG(the01(2))+&
					(1/the01(1)))*res212num+&
					res2the12-&
					LOG(the01(2))*res20112num-&
					res2the0112
					
					u2=res212num
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva12nofix
				endif
			endif
			
			!write(6, *) " done x101 -c4" 	
			if(fix(2).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=the01(1)*res2denum/the01(2) -&
				the01(1)*(res2denum-res201num)/the01(2)
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+troncweib01012
				
				
					
				u1=((the01(1)/the01(2))**2)*res20101numbis+&
					the01(1)*(the01(1)-1)*res2denum/((the01(2))**2)-&
					(res2denum-res201num)*(3*(the01(1)**2)-&
					the01(1))/(the01(2)**2)
				
				res1((nvamax+iter))=&
				u1*res2denum-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib01012square
				if(fix(3).eq.0)then
				
					u1=-LOG(the02(2))*the01(1)*res202num/the01(2)-&
					the01(1)*res2the02/the01(2)+&
					the01(1)*LOG(the02(2))*res20102num/the01(2)+&
					the01(1)*res2the0102/the01(2)
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				!here
				if(fix(4).eq.0)then
				
					u1=the01(1)*the02(1)*res20102num-&
					the01(1)*the02(1)*res202num
					u1=u1/(the01(2)*the02(2))
					
					u2=-the02(1)*res202num/the02(2)
				
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(fix(5).eq.0)then
				
					u1=LOG(the12(2))*the01(1)*res212num/the01(2)+&
					the01(1)*res2the12/the01(2)-&
					the01(1)*LOG(the12(2))*res20112num/the01(2)-&
					the01(1)*res2the0112/the01(2)
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
					
					u1=-the01(1)*the12(1)*res20112num+&
					the01(1)*the12(1)*res212num
					u1=u1/(the01(2)*the12(2))
					
					u2=the12(1)*res212num/the12(2)
				
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				!here
				if(nva01nofix.gt.0) then 
					
					u1=the01(1)*res20101numbis/the01(2)+&
					the01(1)*res2denum/the01(2)-3* &
					the01(1)*(res2denum-res201num)/the01(2)
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01012beta01
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					u1=res20102num-res202num
					u1=u1*the01(1)/the01(2)
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!here 
					u1=-res20112num+res212num
					u1=u1*the01(1)/the01(2)
					
					u2=res212num
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva12nofix
				endif
			endif
			!write(6, *) " done x201 -c4" 	
					
				!!here 
				
			if(fix(3).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=-LOG(the02(2))*res202num-&
				res2the02
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+troncweib02021
				
				res1((nvamax+iter))=&
				((LOG(the02(2)))**2)*res20202num+2* &
				LOG(the02(2))*res2the0202+&
				res2the0202dsquare-&
				((LOG(the02(2)))**2)*res202num-2* &
				LOG(the02(2))*res2the02-&
				res2the0202square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*res2denum -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				res1((nvamax+iter))=&
				res1((nvamax+iter))+&
				troncweib02021square
				
				if(fix(4).eq.0)then
				
					u1=-(the02(1)*LOG(the02(2))+1)* &
					res202num/the02(2)-&
					the02(1)*res2the02/the02(2)+ &
					LOG(the02(2))*the02(1)* &
					res20202num/the02(2)+ &
					the02(1)/the02(2)*res2the0202
					
					u2=-the02(1)*res202num/the02(2)
					
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+&
					troncweib020212
				endif
				
				if(fix(5).eq.0)then
				
					u1=-LOG(the02(2))*LOG(the12(2))*res20212num-&
					LOG(the02(2))*res2the0212-LOG(the12(2))*res2the0212-&
					res2the0212square
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
				
					u1=-LOG(the02(2))*the12(1)*res20212num/the12(2)-&
					the12(1)*res2the0212/the12(2)
					
					u2=the12(1)*res212num/the12(2)
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(nva01nofix.gt.0) then 
					
					u1=LOG(the02(2))*res20102num+&
					res2the0102-res202num*LOG(the02(2))-&
					res2the02
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					
					u1=-LOG(the02(2))*res202num-&
					res2the02+LOG(the02(2))*res20202num+&
					res2the0202
					
					u2=-res202num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02021beta02
					iter = iter +nva02nofix
					
				endif
				if(nva12nofix.gt.0) then 
				
					u1=-LOG(the02(2))*res20212num-&
					res2the0212
					
					u2=res212num
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva12nofix
				endif
			endif
			
             !write(6, *) " done x102 -c4" 	
!here!
			if(fix(4).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=-the02(1)*res202num/the02(2)
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+troncweib02022
				
				
				u1=((the02(1)/the02(2))**2)*res20202num-&
					(the02(1)*(the02(1)-1)/((the02(2))**2))*res202num
				
				
				res1((nvamax+iter))=&
				u1*res2denum-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib02022square
				
				
				if(fix(5).eq.0)then
				
					u1=-the02(1)*LOG(the12(2))*res20212num/the02(2)-&
					the02(1)*res2the0212/the02(2)
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
					
					u1=-the02(1)*the12(1)*res20212num
					u1=u1/(the02(2)*the12(2))
					
					u2=the12(1)*res212num/the12(2)
				
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				
				if(nva01nofix.gt.0) then 
					
					u1=-the02(1)*res202num/the02(2)+&
					the02(1)*res20102num/the02(2)					
					
					u2=res201num
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				
				if(nva02nofix.gt.0) then 
					
					u1=the02(1)*res20202num/the02(2)-&
					the02(1)*res202num/the02(2)
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02022beta02
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					u1=-res20212num*the02(1)/the02(2)
					
					u2=res212num
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva12nofix
				endif
			endif	
!here!
!write(6, *) " done x202 -c4" 	
			if(fix(5).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=LOG(the12(2))*res212num+&
				res2the12
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+&
				LOG(the12(2)*t3(i))*(1-gl12*vet12)+1/ &
				the12(1)
				
				res1((nvamax+iter))=&
				((LOG(the12(2)))**2)*res21212num+2* &
				LOG(the12(2))*res2the1212+&
				res2the1212dsquare+&
				((LOG(the12(2)))**2)*res212num+2* &
				LOG(the12(2))*res2the12+&
				res2the1212square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*res2denum -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				res1((nvamax+iter))=&
				res1((nvamax+iter))-&
				((LOG(the12(2)*t3(i)))**2)*gl12*vet12-&
				(1/(the12(1)**2))
				
!here!
				if(fix(6).eq.0)then
			
					u1=LOG(the12(2))*the12(1)*res21212num/the12(2)+&
					the12(1)*res2the1212/the12(2)+&
					(the12(1)*LOG(the12(2))+1)*res212num/the12(2)+&
					the12(1)*res2the12/the12(2)
					
					u2=the12(1)*res212num/the12(2)
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					res1((nvamax+iter))=&
					res1((nvamax+iter))-&
					(1+LOG(the12(2)*t3(i))*the12(1))*gl12*vet12/the12(2)+&
					(1/the12(2))
				endif
				
				if(nva01nofix.gt.0) then 
					
					u1=-LOG(the12(2))*res20112num-&
					res2the0112+res212num*LOG(the12(2))+&
					res2the12
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					!here
					u1=-LOG(the12(2))*res20212num-&
					res2the0212
					
					u2=-res202num
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva02nofix
					
				endif
				!here
				if(nva12nofix.gt.0) then 
				
					u1=LOG(the12(2))*res21212num+&
					res2the1212+LOG(the12(2))*res212num+&
					res2the12
					
					u2=res212num
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*LOG(the12(2)*t3(i))*gl12*vet12
					
					iter = iter +nva12nofix
				endif
			endif
			
	!write(6, *) " done x112 -c4" 			  
!here!
			if(fix(6).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=the12(1)*res212num/the12(2)
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)-&
				the12(1)*gl12*vet12/the12(2)+&
				the12(1)/the12(2)
				
				u1=((the12(1)/the12(2))**2)*res21212num+&
					(the12(1)*(the12(1)-1)/((the12(2))**2))*res212num
				
				
				res1((nvamax+iter))=&
				u1*res2denum-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))-the12(1)*(the12(1)-1)* &
				gl12*vet12/(the12(2)**2)-&
				the12(1)/(the12(2)**2)
				
				
				
				if(nva01nofix.gt.0) then 
					
					u1=the12(1)*res212num/the12(2)-&
					the12(1)*res20112num/the12(2)					
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				
				if(nva02nofix.gt.0) then 
					
					u1=-the12(1)*res20212num/the12(2)
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!here
					u1=res21212num*the12(1)/the12(2)+&
					the12(1)*res212num/the12(2)
					
					u2=res212num
							
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*the12(1)*gl12*vet12/the12(2)
					iter = iter +nva12nofix
				endif
			endif	
!here!
!write(6, *) " done x212 -c4" 	
			v=res2denum*(su12**vet12)*ri12*vet12

			if(nva01nofix.gt.0) then 
			u1=res201num*(su12**vet12)*ri12*vet12
			
			res1((nvaweib+1):nva01nofix)=&
			ve01nofix(i,:)*u1/v
			res1((nvaweib+1):nva01nofix)=res1((nvaweib+1):nva01nofix)+tronc01

			res1((nvamax12weib12+1):nvamax01)=&
			res20101num*(su12**vet12)*ri12*vet12*ve01square(i,:)
			res1((nvamax12weib12+1):nvamax01)=&
			res1((nvamax12weib12+1):nvamax01)/v
			res1((nvamax12weib12+1):nvamax01)=&
			res1((nvamax12weib12+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
			res1((nvamax12weib12+1):nvamax01)=&
			res1((nvamax12weib12+1):nvamax01)+tronc01square

			end if 

			
			if(nva02nofix.gt.0) then 

			u2=-res202num*(su12**vet12)*ri12*vet12
			
			res1((nva01nofix+1+nvaweib):nva0102)=&
			ve02nofix(i,:)*u2/v
			res1((nva01nofix+1+nvaweib):nva0102)=&
			res1((nva01nofix+1+nvaweib):nva0102)+&
			tronc02
			
			res1((nvamax0112+1):nvamax02)=&
			(res20202num-res202num)*ri12*vet12*(su12**vet12)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)*ve02square(i,:)	
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)/v
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)+tronc02square

			end if 


			
			if(nva12nofix.gt.0) then 

			u3=res212num+&
			(1-gl12*vet12)*res2denum
			u3=u3*(su12**vet12)*ri12*vet12
			
			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v
			
			res1((nvamax0212+1):nvamax12)=&
			(1-gl12*vet12)*res212num+res21212num+&
			((1-gl12*vet12)**2)*res2denum+&
			(2-gl12*vet12)*res212num-gl12*vet12*res2denum
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*(su12**vet12)*ri12*vet12
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*ve12square(i,:)
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)/v
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)			


			end if 

			
			
			if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 
			kfix=nvamax01+1
			lfix=kfix-1+nva02nofix

			do j=1,nva01nofix
			   res1(kfix:lfix)=(-res202num+res20102num)*ri12*vet12
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva02nofix
			end do

			end if 

			if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 
			
			kfix=nvamax0102+1
			lfix=kfix-1+nva12nofix

			do j=1,nva01nofix
			   res1(kfix:lfix)=&
			   (1-gl12*vet12)*res201num+&
			   res212num-res20112num
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve12nofix(i,:)
		           res1(kfix:lfix)=&
			   res1(kfix:lfix)*ri12*vet12*(su12**vet12)
			   res1(kfix:lfix)=res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do


			end if 

			if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 
			
			kfix=nvamax02+1
			lfix=kfix-1+nva12nofix

			do j=1,nva02nofix
			   res1(kfix:lfix)=&
			   res20212num+(1-gl12*vet12)*res202num
			   res1(kfix:lfix)=&
			   -res1(kfix:lfix)*ve02nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ri12*vet12*(su12**vet12)
			   res1(kfix:lfix)=res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-ve02nofix(i,j)*ve12nofix(i,:)*u2*u3
			   res1(kfix:lfix)=res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do

			end if 

!here!
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				!write(6, *) "start c=5" 
				call fonct(t1(i),the01,ri01,gl01,su01)
				call fonct(t1(i),the02,ri02,gl02,su02)
				call fonct(t1(i),the12,ri12,gl12,su12)
				
				iter = 0
				nweib= 0
			if(fix(1).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				res1(nweib)=-LOG(the01(2)*t1(i))*gl01*vet01 +&
				troncweib01011+&
				LOG(the01(2)*t1(i))+(1/the01(1))
				
				res1((nvamax+iter))=&
				-((LOG(the01(2)*t1(i)))**2)*gl01*vet01 +&
				troncweib01011square-1/(the01(1)**2)
				!here!!
				if(fix(2).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=&
					-(LOG(the01(2)*t1(i))*the01(1)+1)* &
					gl01*vet01/the01(2)+&
					troncweib010112+(1/the01(2))
				endif
				if(fix(3).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				!here!!
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					-ve01nofix(i,:)*gl01*vet01*LOG(the01(2)*t1(i))+&
					troncweib01011beta01
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif 
			!here!!
			if(fix(2).eq.0)then
			
				iter = iter + 1
				nweib = nweib + 1
				
				res1(nweib)=-the01(1)*gl01*vet01/the01(2) +&
				troncweib01012+&
				the01(1)/the01(2)

				res1(nvamax+iter)=&
				-the01(1)*(the01(1)-1)*gl01*vet01/(the01(2)**2)+&
				troncweib01012square-&
				the01(1)/(the01(2)**2)
				
				if(fix(3).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
				
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					-ve01nofix(i,:)*gl01*vet01*the01(1)/the01(2)+&
					troncweib01012beta01
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif
			!here!!
			if(fix(3).eq.0)then
				iter = iter + 1
				nweib = nweib +1
				res1(nweib)=-LOG(the02(2)*t1(i))*gl02*vet02 +&
				troncweib02021
				
				res1((nvamax+iter))=&
				-((LOG(the02(2)*t1(i)))**2)*gl02*vet02+&
				troncweib02021square
				!here!!
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=&
					-(LOG(the02(2)*t1(i))*the02(1)+1)* &
					gl02*vet02/the02(2)+&
					troncweib02021
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter):(nvamax+iter+nva02nofix))=&
					-ve02nofix(i,:)*LOG(the02(2)*t1(i))*gl12*vet12+&
					troncweib02021beta02
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif 
			!here!!
			if(fix(4).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				
				res1(nweib)=-the02(1)*gl02*vet02/the02(2)+&
				troncweib02022
				
				res1(nvamax+iter)=&
				-the02(1)*(the02(1)-1)*gl02*vet02/(the02(2)**2)+&
				troncweib02022square
				!here!!
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					-ve02nofix(i,:)*gl02*vet02*the02(1)/the02(2)+&
					troncweib02022beta02
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif
			!here!!
			
			if(fix(5).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				res1(nweib)=LOG(the12(2)*t1(i))*gl12*vet12+&
				(1/the12(1))
				res1(nvamax+iter)=((LOG(the12(2)*t1(i)))**2)*gl12*vet12-&
				(1/(the12(1)**2))
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=&
					(LOG(the12(2)*t1(i))*the12(1)+1)* &
					gl12*vet12/the12(2)+(1/the12(2))
				endif
				if(nva01nofix.gt.0) then 
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*LOG(the12(2)*t1(i))*gl12*vet12
					call fonct(t3(i),the12,ri12,gl12,su12)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*LOG(the12(2)*t3(i))*gl12*vet12
					
					iter = iter +nva12nofix
				endif
				
				call fonct(t3(i),the12,ri12,gl12,su12)
				
				res1(nweib)=res1(nweib)-&
				LOG(the12(2)*t3(i))*gl12*vet12+&
				LOG(the12(2)*t3(i))
				iter=iter-nva12nofix-nva02nofix-nva01nofix
				if(fix(6).eq.0)then
					res1((nvamax+iter))=res1((nvamax+iter))-&
					(LOG(the12(2)*t3(i))*the12(1)+1)* &
					gl12*vet12/the12(2)
					iter=iter-1
				endif
				
				res1(nvamax+iter)=res1(nvamax+iter)-&
				((LOG(the12(2)*t3(i)))**2)*gl12*vet12
				
				iter=iter+nva12nofix+nva02nofix+nva01nofix
				
				if(fix(6).eq.0)then
					iter=iter+1
				endif
				call fonct(t1(i),the12,ri12,gl12,su12)
			endif
			!here!!
			
			if(fix(6).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				res1(nweib)=the12(1)*gl12*vet12/the12(2)+&
				the12(1)/the12(2)
				res1(nvamax+iter)=the12(1)*(-1+&
				the12(1))*gl12*vet12/(the12(2)**2)-&
				(the12(1)/(the12(2)**2))
				
				if(nva01nofix.gt.0) then 
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*the12(1)*gl12*vet12/the12(2)
					call fonct(t3(i),the12,ri12,gl12,su12)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*the12(1)*gl12*vet12/the12(2)
					
					iter = iter +nva12nofix
				endif
				
				call fonct(t3(i),the12,ri12,gl12,su12)
				
				res1(nweib)=res1(nweib)-&
				the12(1)*gl12*vet12/the12(2)
				
				iter=iter-nva12nofix-nva02nofix-nva01nofix
				
				res1(nvamax+iter)=res1(nvamax+iter)-&
				the12(1)*(the12(1)-1)*gl12*vet12/(the12(2)**2)
				
				iter=iter+nva12nofix+nva02nofix+nva01nofix
				
				call fonct(t1(i),the12,ri12,gl12,su12)
			endif
			
			
				if(nva01nofix.gt.0) then 

				res1((nvaweib+1):(nva01nofix+nvaweib))=&
				-ve01nofix(i,:)*gl01*vet01+&
				ve01nofix(i,:)+&
				tronc01
				res1((nvamax12weib12+1):nvamax01)=&
				-ve01square(i,:)*gl01*vet01+&
				tronc01square

				if(nva02nofix.gt.0) then 
				res1((nvamax01+1):nvamax0102)=0
				end if 

				if(nva12nofix.gt.0) then
				res1((nvamax0102+1):nvamax0112)=0
				end if 

				end if 

				if(nva02nofix.gt.0) then 

				res1((nva01nofix+1+nvaweib):nva0102)=&
				-ve02nofix(i,:)*gl02*vet02+tronc02
				

				res1((nvamax0112+1):nvamax02)=&
				-ve02square(i,:)*gl02*vet02+&
				tronc02square

				if(nva12nofix.gt.0) then
				res1((nvamax02+1):nvamax0212)=0
				end if 

				end if 
			
				if(nva12nofix.gt.0) then 

				res1((nvamax0212+1):nvamax12)=&
				gl12*vet12
				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*gl12*vet12+&
				ve12nofix(i,:)

				call fonct(t3(i),the12,ri12,gl12,su12)
				res1((nva0102+1):nvamax)=&
				res1((nva0102+1):nvamax)-ve12nofix(i,:)*gl12*vet12
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)-gl12*vet12
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)*ve12square(i,:)
				

				end if 
				
				
				!here!!
                         else
                            if(c(i).eq.6)then ! vivant ???

				!write(6, *) "start c=6" 
				call fonct(t3(i),the01,ri01,gl01,su01)
				call fonct(t3(i),the02,ri02,gl02,su02)
				call fonct(t3(i),the12,ri12,gl12,su12)

			call qgaussweibbetaderiv(t1(i),t3(i),the01,&
			the02,the12,res2denum,res2the01,res2the02,&
			res2the12,res2thenum,res2thenumsquare,&
			res2the0101,res2the0102,res2the0112,res2the0202,&
			res2the0212,res2the1212,res2the0101square,&
			res2the0102square,res2the0112square,&
			res2the0202square,res2the0212square,&
			res2the1212square,res2the0101dsquare,&
			res2the0202dsquare,res2the1212dsquare,&
			res201num,res202num,res212num,res20101num,&
			res20101numbis,&
			res20102num,res20112num,&
			res20202num,res20212num,res21212num,&
			vet01,vet02,vet12)
			
			iter = 0
			nweib = 0
			u1=((su12**vet12)*res2denum+(su01**vet01)*(su02**vet02))
				
				
				
			if(fix(1).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=(LOG(the01(2))+(1/the01(1)))*res2denum +&
				res2thenum-LOG(the01(2))*(res2denum-res201num)-&
				res2the01
				
				v=v*(su12**vet12)-&
				LOG(the01(2)*t3(i))*gl01*vet01*(su01**vet01)*(su02**vet02)
				
				res1(nweib)=v/u1
				
				res1(nweib)=res1(nweib)+troncweib01011
				
				res1((nvamax+iter))=&
				LOG(the01(2))*2/the01(1)*res2denum+&
				res2thenum*2/the01(1)+&
				(LOG(the01(2))**2)*res2denum+&
				LOG(the01(2))*2*res2thenum+&
				res2thenumsquare-3*(LOG(the01(2))**2)* &
				(res2denum-res201num)-&
				LOG(the01(2))*6*res2the01-3*res2the0101square-&
				LOG(the01(2))*2/the01(1)*(res2denum-res201num)- &
				res2the01*2/the01(1)+ &
				(LOG(the01(2))**2)*res20101numbis+&
				LOG(the01(2))*2*res2the0101+&
				res2the0101dsquare
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*(su12**vet12)+&
				gl01*vet01*(LOG(the01(2)*t3(i))**2)*(-1+&
				gl01*vet01)*(su01**vet01)*(su02**vet02)
				res1((nvamax+iter))=&
				(res1((nvamax+iter))*u1-v*v)/(u1*u1)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+&
				troncweib01011square
				!ici
				if(fix(2).eq.0)then
					
					u2=res2denum*the01(1)/the01(2)-&
					(res2denum-res201num)*the01(1)/the01(2)
					u2=u2*(su12**vet12)-&
					the01(1)*gl01*vet01*(su01**vet01)*(su02**vet02)/the01(2)
					
					u3=the01(1)*LOG(the01(2))*res20101numbis/the01(2)+&
					the01(1)*res2the0101/the01(2)+&
					(the01(1)*LOG(the01(2))+2)*res2denum/the01(2)+&
					the01(1)*res2thenum/the01(2)-&
					(3*the01(1)*LOG(the01(2))+2)* &
					(res2denum-res201num)/the01(2)-3* &
					the01(1)*res2the01/the01(2)
					
					u3=u3*(su12**vet12)+&
					gl01*vet01*(su01**vet01)*(su02**vet02)*(-1-&
					LOG(the01(2)*t3(i))*the01(1)+&
					gl01*vet01*LOG(the01(2)*t3(i))*the01(1))/the01(2)
					
					iter = iter +1
					res1((nvamax+iter))=&
					(u3*u1-v*u2)/(u1*u1)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+&
					troncweib010112
				endif
				!ici
				if(fix(3).eq.0)then
				
					u3=-LOG(the01(2))*LOG(the02(2))*res202num-&
					LOG(the01(2))*res2the02-LOG(the02(2))*res2the02-&
					res2the0202square-LOG(the02(2))*res202num/the01(1)-&
					res2the02/the01(1)+&
					LOG(the01(2))*LOG(the02(2))*res20102num+&
					LOG(the01(2))*res2the0102+&
					LOG(the02(2))*res2the0102+&
					res2the0102square
					u3=u3*(su12**vet12)+&
					gl01*vet01*gl02*vet02*(su01**vet01)*(su02**vet02)* &
					LOG(the01(2)*t3(i))*LOG(the02(2)*t3(i))
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					u2=u2*(su12**vet12)-&
					LOG(the02(2)*t3(i))*gl02*vet02* &
					(su01**vet01)*(su02**vet02)
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!ici
				if(fix(4).eq.0)then
				
					u3=-LOG(the01(2))*the02(1)*res202num/the02(2)-&
					the02(1)*res2the02/the02(2)-&
					the02(1)*res202num/(the01(1)*the02(2))+&
					LOG(the01(2))*the02(1)*res20102num/the02(2)+&
					the02(1)*res2the0102/the02(2)
					
					u3=u3*(su12**vet12)+&
					gl01*vet01*gl02*vet02*(su01**vet01)* &
					(su02**vet02)*the02(1)* &
					LOG(the01(2)*t3(i))/the02(2)
					
					u2=-the02(1)*res202num/the02(2)
					u2=u2*(su12**vet12)-&
					the02(1)*gl02*vet02*(su01**vet01)* &
					(su02**vet02)/the02(2)
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!ici
				if(fix(5).eq.0)then
				
				
					u3=LOG(the01(2))*LOG(the12(2))*res212num+&
					LOG(the01(2))*res2the12+LOG(the12(2))*res2the12+&
					res2the1212square+LOG(the12(2))*res212num/the01(1)+&
					res2the12/the01(1)-&
					LOG(the01(2))*LOG(the12(2))*res20112num-&
					LOG(the01(2))*res2the0112-&
					LOG(the12(2))*res2the0112-&
					res2the0112square

					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)* &
					LOG(the12(2)*t3(i))* &
					((LOG(the01(2))+(1/the01(1)))*res2denum +&
					res2thenum-LOG(the01(2))*(res2denum-res201num)-&
					res2the01)
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)-&
					LOG(the12(2)*t3(i))*gl12*vet12* &
					res2denum*(su12**vet12)
						
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					v*u2)/(u1*u1)
					
				endif
				!ici
				if(fix(6).eq.0)then
				
					
					u3=&
					LOG(the01(2))*the12(1)*res212num/the12(2)+&
					the12(1)/the12(2)*res2the12+res212num* &
					the12(1)/(the12(2)*the01(1))-&
					LOG(the01(2))*the12(1)*res20112num/the12(2)-&
					the12(1)/the12(2)*res2the0112
					
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*the12(1)* &
					((LOG(the01(2))+(1/the01(1)))*res2denum +&
					res2thenum-LOG(the01(2))*(res2denum-res201num)-&
					res2the01)/the12(2)
					
					
					u2=the12(1)*res212num/the12(2)
					
					u2=u2*(su12**vet12)-&
					the12(1)*gl12*vet12*(su12**vet12)* &
					res2denum/the12(2)
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				
				!ici
				if(nva01nofix.gt.0) then 
					
					u3=res2the0101+ &
					LOG(the01(2))*res20101numbis-&
					LOG(the01(2))*3*(res2denum-res201num)-&
					res2the01*3-(res2denum-res201num)/the01(1)+&
					(LOG(the01(2))+(1/the01(1)))*res2denum+&
					res2thenum
					
					u3= u3*(su12**vet12)-&
					LOG(the01(2)*t3(i))*gl01*vet01*(1-&
					gl01*vet01)*(su01**vet01)*(su02**vet02)
					
					u2=res201num*(su12**vet12)-&
					gl01*vet01*(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01011beta01
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					u3=-(LOG(the01(2))+&
					(1/the01(1)))*res202num-&
					res2the02+&
					LOG(the01(2))*res20102num+&
					res2the0102
					
					u3=u3*(su12**vet12)+&
					LOG(the01(2)*t3(i))*gl01*vet01* &
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
					u2=-res202num*(su12**vet12)-&
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					iter = iter +nva02nofix
				endif
				!ici
				if(nva12nofix.gt.0) then 
					u3=(LOG(the01(2))+&
					(1/the01(1)))*res212num+&
					res2the12-&
					LOG(the01(2))*res20112num-&
					res2the0112
					
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)* &
					((LOG(the01(2))+(1/the01(1)))*res2denum +&
					res2thenum-LOG(the01(2))*(res2denum-res201num)-&
					res2the01)
					
					u2=res212num*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*res2denum
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
				endif
			endif
			
			!write(6, *) " done x101 -c6" 	
			if(fix(2).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=the01(1)*res2denum/the01(2) -&
				the01(1)*(res2denum-res201num)/the01(2)
				
				v= v*(su12**vet12)-&
				the01(1)*gl01*vet01*(su01**vet01)* &
				(su02**vet02)/the01(2)
				
				res1(nweib)=v/u1
				res1(nweib)=res1(nweib)+troncweib01012
				
				
				u2=((the01(1)/the01(2))**2)*res20101numbis+&
					the01(1)*(the01(1)-1)*res2denum/((the01(2))**2)-&
					(res2denum-res201num)*(3*(the01(1)**2)-&
					the01(1))/(the01(2)**2)
					
				u2=u2*(su12**vet12)+&
				gl01*vet01*(su01**vet01)*(su02**vet02)* &
				the01(1)*(1-the01(1)+gl01*vet01*the01(1))/(the01(2)**2)
				
				res1((nvamax+iter))=&
				(u2*u1-v*v)/(u1*u1)
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib01012square
				
				if(fix(3).eq.0)then
				!ici
					u3=-LOG(the02(2))*the01(1)*res202num/the01(2)-&
					the01(1)*res2the02/the01(2)+&
					the01(1)*LOG(the02(2))*res20102num/the01(2)+&
					the01(1)*res2the0102/the01(2)
					
					u3=u3*(su12**vet12)+&
					gl01*vet01*gl02*vet02*the01(1)* &
					LOG(the02(2)*t3(i))*(su01**vet01)* &
					(su02**vet02)/the01(2)
					
					
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					u2=u2*(su12**vet12)-&
					LOG(the02(2)*t3(i))*gl02*vet02* &
					(su01**vet01)*(su02**vet02)
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!ici!
				if(fix(4).eq.0)then
				
					u3=the01(1)*the02(1)*res20102num-&
					the01(1)*the02(1)*res202num
					u3=u3/(the01(2)*the02(2))
					u3=u3*(su12**vet12)+&
					gl01*vet01*gl02*vet02* &
					(su01**vet01)*(su02**vet02)* &
					the01(1)*the02(1)/(the01(2)*the02(2))
					
					u2=-the02(1)*res202num/the02(2)
					u2=u2*(su12**vet12)-&
					the02(1)*gl02*vet02*(su01**vet01)* &
					(su02**vet02)/the02(2)
				
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!ici
				if(fix(5).eq.0)then
				
					u3=LOG(the12(2))*the01(1)*res212num/the01(2)+&
					the01(1)*res2the12/the01(2)-&
					the01(1)*LOG(the12(2))*res20112num/the01(2)-&
					the01(1)*res2the0112/the01(2)
					
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*LOG(the12(2)* &
					t3(i))*(the01(1)*res2denum/the01(2) -&
					the01(1)*(res2denum-res201num)/the01(2))
				
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)-&
					LOG(the12(2)*t3(i))*gl12*vet12* &
					(su12**vet12)*res2denum
						
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					v*u2)/(u1*u1)
					
				endif
				!ici
				if(fix(6).eq.0)then
					
					u3=-the01(1)*the12(1)*res20112num+&
					the01(1)*the12(1)*res212num
					u3=u3/(the01(2)*the12(2))
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*the12(1)* &
					(the01(1)*res2denum/the01(2) -&
					the01(1)*(res2denum-res201num)/the01(2))/ &
					the12(2)
					
					
					u2=the12(1)*res212num/the12(2)
					u2=u2*(su12**vet12)-&
					the12(1)*gl12*vet12*(su12**vet12)* &
					res2denum/the12(2)
				
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!ici
				if(nva01nofix.gt.0) then 
					
					u3=the01(1)*res20101numbis/the01(2)+&
					the01(1)*res2denum/the01(2)-&
					(3*the01(1)*(res2denum-res201num)/the01(2))
					u3=u3*(su12**vet12)+&
					the01(1)*gl01*vet01*(gl01*vet01-1)* &
					(su01**vet01)*(su02**vet02)/the01(2)
					
					u2=res201num*(su12**vet12)-&
					gl01*vet01*(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01012beta01
					
					iter = iter +nva01nofix
				endif
				!ici
				if(nva02nofix.gt.0) then 
					
					u3=res20102num-res202num
					u3=u3*the01(1)/the01(2)
					u3=u3*(su12**vet12)+&
					the01(1)*gl01*vet01*gl02*vet02* &
					(su01**vet01)*(su02**vet02)/the01(2)
					
					u2=-res202num*(su12**vet12)-&
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!ici
					u3=-res20112num+res212num
					u3=u3*the01(1)/the01(2)
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)* &
					(the01(1)*res2denum/the01(2) -&
					the01(1)*(res2denum-res201num)/the01(2))
				
					
					u2=res212num*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*res2denum
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					iter = iter +nva12nofix
				endif
			endif
			
					!write(6, *) " done x201 -c6" 	
				!ici
				
			if(fix(3).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=-LOG(the02(2))*res202num-&
				res2the02
				
				v=v*(su12**vet12)-&
				LOG(the02(2)*t3(i))*gl02*vet02* &
				(su01**vet01)*(su02**vet02)
				
				res1(nweib)=v/u1
				res1(nweib)=res1(nweib)+troncweib02021
				
				res1((nvamax+iter))=&
				((LOG(the02(2)))**2)*res20202num+&
				(2*LOG(the02(2))*res2the0202)+&
				res2the0202dsquare-&
				((LOG(the02(2)))**2)*res202num-&
				(2*LOG(the02(2))*res2the02)-&
				res2the0202square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*(su12**vet12)+&
				gl02*vet02*(su01**vet01)*(su02**vet02)* &
				(LOG(the02(2)*t3(i))**2)*(gl02*vet02-1)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*u1 -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				res1((nvamax+iter))=&
				res1((nvamax+iter))+&
				troncweib02021square
				
				!ici
				if(fix(4).eq.0)then
				
					u3=-(the02(1)*LOG(the02(2))+1)* &
					res202num/the02(2)-&
					the02(1)*res2the02/the02(2)+ &
					LOG(the02(2))*the02(1)* &
					res20202num/the02(2)+ &
					the02(1)/the02(2)*res2the0202
					
					u3=u3*(su12**vet12)+&
					gl02*vet02*(su01**vet01)*(su02**vet02)* &
					(-1-LOG(the02(2)*t3(i))*the02(1)+&
					gl02*vet02*LOG(the02(2)*t3(i))* &
					the02(1))/the02(2)
					
					u2=-the02(1)*res202num/the02(2)
					u2=u2*(su12**vet12)-&
					the02(1)*gl02*vet02*(su01**vet01)* &
					(su02**vet02)/the02(2)
					
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+&
					troncweib020212
				endif
				!ici
				if(fix(5).eq.0)then
				
					u3=-LOG(the02(2))*LOG(the12(2))*res20212num-&
					LOG(the02(2))*res2the0212-LOG(the12(2))*res2the0212-&
					res2the0212square
					
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*LOG(the12(2)* &
					t3(i))*(-LOG(the02(2))*res202num-&
					res2the02)
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)-&
					LOG(the12(2)*t3(i))*gl12*vet12* &
					(su12**vet12)*res2denum
						
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
				endif
				!ici
				if(fix(6).eq.0)then
				
					u3=-LOG(the02(2))*the12(1)*res20212num/the12(2)-&
					the12(1)*res2the0212/the12(2)
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*the12(1)* &
					(-1*LOG(the02(2))*res202num-&
					res2the02)/the12(2)
				
					
					
					u2=the12(1)*res212num/the12(2)
					u2=u2*(su12**vet12)-the12(1)*gl12*vet12* &
					(su12**vet12)*res2denum/the12(2)
					
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
				endif
				!ici
				if(nva01nofix.gt.0) then 
					
					u3=LOG(the02(2))*res20102num+&
					res2the0102-res202num*LOG(the02(2))-&
					res2the02
					
					u3=u3*(su12**vet12)+&
					LOG(the02(2)*t3(i))*gl01*vet01* &
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
					u2=res201num*(su12**vet12)-&
					gl01*vet01*(su01**vet01)* &
					(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
				endif
				!ici
				if(nva02nofix.gt.0) then 
					
					u3=-LOG(the02(2))*res202num-&
					res2the02+LOG(the02(2))*res20202num+&
					res2the0202
					u3=u3*(su12**vet12)+&
					LOG(the02(2)*t3(i))*gl02*vet02* &
					(su01**vet01)*(su02**vet02)* &
					(-1+gl02*vet02)
					
					u2=-res202num*(su12**vet12)-&
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02021beta02
					iter = iter +nva02nofix
					
				endif
				!ici
				if(nva12nofix.gt.0) then 
				
					u3=-LOG(the02(2))*res20212num-&
					res2the0212
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)* &
					(-1*LOG(the02(2))*res202num-&
					res2the02)
					
					u2=res212num*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*res2denum
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
				endif
			endif
			
             !write(6, *) " done x102 -c6" 	
			!ici
			if(fix(4).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=-the02(1)*res202num/the02(2)
				
				v=v*(su12**vet12)-&
				gl02*vet02*the02(1)*(su01**vet01)* &
				(su02**vet02)/the02(2)
				
				res1(nweib)=v/u1
				res1(nweib)=res1(nweib)+troncweib02022
				
				u3=((the02(1)/the02(2))**2)*res20202num-&
					(the02(1)*(the02(1)-1)/((the02(2))**2))*res202num

				
				u3=u3*(su12**vet12)+&
				gl02*vet02*(su01**vet01)*(su02**vet02)* &
				the02(1)*(1-the02(1)+the02(1)*gl02*vet02)/ &
				(the02(2)**2)
				
				res1((nvamax+iter))=&
				u1*u3-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib02022square
				
				!ici
				if(fix(5).eq.0)then
				
					u3=-the02(1)*LOG(the12(2))*res20212num/the02(2)-&
					the02(1)*res2the0212/the02(2)
					u3=u3*(su12**vet12)+&
					gl12*vet12*(su12**vet12)*LOG(the12(2)*t3(i))* &
					(the02(1)*res202num/the02(2))
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)-&
					LOG(the12(2)*t3(i))*gl12*vet12* &
					(su12**vet12)*res2denum
						
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
				endif
				!ici
				if(fix(6).eq.0)then
					
					u3=-the02(1)*the12(1)*res20212num
					u3=u3/(the02(2)*the12(2))
					u3=u3*(su12**vet12)+&
					gl12*vet12*(su12**vet12)*the12(1)* &
					(the02(1)*res202num/the02(2))/the12(2)
					
					u2=the12(1)*res212num/the12(2)
					u2=u2*(su12**vet12)-&
					the12(1)*gl12*vet12*(su12**vet12)* &
					res2denum/the12(2)
				
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
				endif
				
				!ici
				if(nva01nofix.gt.0) then 
					
					u3=-the02(1)*res202num/the02(2)+&
					the02(1)*res20102num/the02(2)		
					u3=u3*(su12**vet12)+&
					the02(1)*gl02*vet02*gl01*vet01* &
					(su01**vet01)*(su02**vet02)/the02(2)
					
					u2=res201num*(su12**vet12)-&
					gl01*vet01*(su01**vet01)*(su02**vet02)
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
				endif
				!ici
				if(nva02nofix.gt.0) then 
					
					u3=the02(1)*res20202num/the02(2)-&
					the02(1)*res202num/the02(2)
					u3=u3*(su12**vet12)+&
					the02(1)*gl02*vet02*(su01**vet01)* &
					(su02**vet02)*(gl02*vet02-1)/the02(2)
					
					u2=-res202num*(su12**vet12)-&
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02022beta02
					
					iter = iter +nva02nofix
				endif
				!ici
				if(nva12nofix.gt.0) then 
					u3=-res20212num*the02(1)/the02(2)
					u3=u3*(su12**vet12)+&
					gl12*vet12*(su12**vet12)* &
					(the02(1)*res202num/the02(2))
					
					u2=res212num*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*res2denum
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					iter = iter +nva12nofix
				endif
			endif	

			!ici
			!write(6, *) " done x202 -c6" 	
			if(fix(5).eq.0)then
			
			
				iter = iter + 1
				nweib = nweib + 1
				
				v=LOG(the12(2))*res212num+&
				res2the12
				v=v*(su12**vet12)-&
				LOG(the12(2)*t3(i))*gl12*vet12* &
				(su12**vet12)*res2denum
				
				res1(nweib)=v/u1
				
				res1((nvamax+iter))=&
				((LOG(the12(2)))**2)* &
				res21212num+2*LOG(the12(2))*res2the1212+&
				res2the1212dsquare+&
				((LOG(the12(2)))**2)* &
				res212num+2*LOG(the12(2))*res2the12+&
				res2the1212square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*(su12**vet12)+&
				gl12*vet12*(su12**vet12)*LOG(the12(2)* &
				t3(i))*(-2*(LOG(the12(2))*res212num+&
				res2the12)-LOG(the12(2)*t3(i))*res2denum* &
				(1-gl12*vet12))
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*u1 -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				
				!ici
				if(fix(6).eq.0)then
			
					u2=the12(1)*res212num/the12(2)
					
					u3=LOG(the12(2))*the12(1)*res21212num/the12(2)+&
					the12(1)*res2the1212/the12(2)+&
					(the12(1)*LOG(the12(2))+1)*res212num/the12(2)+&
					the12(1)*res2the12/the12(2)
					u3=u3*(su12**vet12)+&
					gl12*vet12*(su12**vet12)* &
					(-1*LOG(the12(2)*t3(i))*u2 - &
					(LOG(the12(2))*res212num+&
					res2the12)*the12(1)/the12(2)+&
					res2denum*(-1-LOG(the12(2)*t3(i))*the12(1)+&
					gl12*vet12*the12(1)*LOG(the12(2)*t3(i)))/ &
					the12(2))
					
					u2=u2*(su12**vet12)-&
					the12(1)*gl12*vet12*(su12**vet12)* &
					res2denum/the12(2)
					
					
					iter = iter +1
					
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
				endif
				
				if(nva01nofix.gt.0) then 
					!ici
					u3=-LOG(the12(2))*res20112num-&
					res2the0112+res212num*LOG(the12(2))+&
					res2the12
					u3=u3*(su12**vet12)-&
					LOG(the12(2)*t3(i))*gl12*vet12* &
					(su12**vet12)*res201num
					
					
					u2=res201num*(su12**vet12)-&
					gl01*vet01*(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					!ici
					u3=-LOG(the12(2))*res20212num-&
					res2the0212
					u3=u3*(su12**vet12)-&
					LOG(the12(2)*t3(i))*gl12*vet12* &
					(su12**vet12)*(-res202num)
					
					u2=-res202num*(su12**vet12)-&
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					iter = iter +nva02nofix
					
				endif
				if(nva12nofix.gt.0) then 
					!ici
					u3=LOG(the12(2))*res21212num+&
					res2the1212+LOG(the12(2))*res212num+&
					res2the12
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)* &
					(LOG(the12(2))*res212num+&
					res2the12)-gl12*vet12*(su12**vet12)* &
					LOG(the12(2)*t3(i))* &
					(res212num + res2denum*(1-gl12*vet12))
					
					u2=res212num*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*res2denum
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
				endif
			endif
			
			  !ici
			!write(6, *) " done x112 -c6" 	
			if(fix(6).eq.0)then
			
			
				iter = iter + 1
				nweib = nweib + 1
				
				v=the12(1)*res212num/the12(2)
				v=v*(su12**vet12)-&
				the12(1)*gl12*vet12*(su12**vet12)* &
				res2denum/the12(2)
				
				res1(nweib)=v/u1
				
				u3=((the12(1)/the12(2))**2)*res21212num+&
					(the12(1)*(the12(1)-1)/((the12(2))**2))*res212num
				u3=u3*(su12**vet12)+&
				gl12*vet12*(su12**vet12)*the12(1)* &
				(-2*(the12(1)*res212num/the12(2))-&
				res2denum*(-1+the12(1)-&
				gl12*vet12*the12(1))/the12(2))/the12(2)
				
				
				res1((nvamax+iter))=&
				u1*u3-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				
				
				if(nva01nofix.gt.0) then 
					!ici
					u3=the12(1)*res212num/the12(2)-&
					the12(1)*res20112num/the12(2)					
					u3=u3*(su12**vet12)-&
					the12(1)*gl12*vet12*res201num* &
					(su12**vet12)/the12(2)
					
					u2=res201num*(su12**vet12)-&
					gl01*vet01*(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
				endif
				
				if(nva02nofix.gt.0) then 
					!ici
					u3=-the12(1)*res20212num/the12(2)
					u3=u3*(su12**vet12)+&
					the12(1)*gl12*vet12*(su12**vet12)* &
					(res202num)/the12(2)
					
					u2=-res202num*(su12**vet12)-&
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!ici
					u3=res21212num*the12(1)/the12(2)+&
					the12(1)*res212num/the12(2)
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)* &
					(the12(1)*res212num/the12(2))-&
					gl12*vet12*(su12**vet12)*the12(1)* &
					(res212num+res2denum*(1-gl12*vet12))/ &
					the12(2)
					
					u2=res212num*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*res2denum
							
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
				endif
			endif	
			
			!write(6, *) " done x212 -c6" 	
				v=(su12**vet12)*res2denum+&
				(su01**vet01)*(su02**vet02)

				if(nva01nofix.gt.0) then 

				u1=(-gl01*vet01)*(su01**vet01)*(su02**vet02)+&
				(su12**vet12)*res201num
			        res1((nvaweib+1):(nvaweib+nva01nofix))=&
				ve01nofix(i,:)*u1/v
				res1((nvaweib+1):(nvaweib+nva01nofix))=&
				res1((nvaweib+1):(nvaweib+nva01nofix))+tronc01

				res1((nvamax12weib12+1):nvamax01)=&
				res20101num*(su12**vet12)
				resint=(su01**vet01)*(su02**vet02)
				resint=resint*gl01*vet01*(1-gl01*vet01)
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)-resint
				
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)*ve01square(i,:)
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)/v

				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)+tronc01square

				end if 
				
				if(nva02nofix.gt.0) then 

				u2=-gl02*vet02*(su01**vet01)*(su02**vet02)
				u2=u2-(su12**vet12)*res202num
				res1((nva01nofix+1+nvaweib):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1+nvaweib):nva0102)=&
				res1((nva01nofix+1+nvaweib):nva0102)+&
				tronc02

				res1((nvamax0112+1):nvamax02)=&
				(su12**vet12)*(res20202num-res202num)-&
				gl02*vet02*(su01**vet01)*(su02**vet02)*(1-gl02*vet02)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)*ve02square(i,:)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)/v

				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)+tronc02square

				
				end if 

				if(nva12nofix.gt.0) then 
				
				u3=-gl12*vet12*(su12**vet12)*res2denum+&
				(su12**vet12)*res212num
				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v

				res1((nvamax0212+1):nvamax12)=&
				res21212num+(1-gl12*vet12)*res212num-&
				gl12*vet12*(1-gl12*vet12)*res2denum-&
				gl12*vet12*res212num
				resint=(su12**vet12)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)*resint*ve12square(i,:)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)/v
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)


				end if 
                        	
				
				if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 

				kfix=nvamax01+1
				lfix=kfix-1+nva02nofix

				do j=1,nva01nofix
			   		res1(kfix:lfix)=&
					-(res202num-res20102num)*(su12**vet12)
					resint=(su01**vet01)*(su02**vet02)
					res1(kfix:lfix)=res1(kfix:lfix)+&
					gl01*vet01*gl02*vet02*resint
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)*v
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva02nofix
				end do

				end if 

				if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 

				
				kfix=nvamax0102+1
				lfix=kfix-1+nva12nofix

				do j=1,nva01nofix
			  		res1(kfix:lfix)=-gl12*vet12*res201num+&
					res212num-res20112num
					resint=(su12**vet12)
			   		res1(kfix:lfix)=&
			  		res1(kfix:lfix)*resint*ve01nofix(i,j)*ve12nofix(i,:)
		           		res1(kfix:lfix)=&
			   		res1(kfix:lfix)*v
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
			   		res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva12nofix
				end do

				end if 

				
				if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 

				kfix=nvamax02+1
				lfix=kfix-1+nva12nofix

				do j=1,nva02nofix
			  		res1(kfix:lfix)=&
					-res20212num+gl12*vet12*res202num
					resint=(su12**vet12)
					res1(kfix:lfix)=&
					res1(kfix:lfix)*resint*ve02nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
					res1(kfix:lfix)*v
					res1(kfix:lfix)=&
					res1(kfix:lfix)-u2*u3*ve02nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva12nofix
				end do

				end if 


				
                            else ! passage 0-->2  
				
				!write(6, *) "start c=7" 
				
			call qgaussweibbetaderiv(t1(i),t3(i),the01,&
			the02,the12,res2denum,res2the01,res2the02,&
			res2the12,res2thenum,res2thenumsquare,&
			res2the0101,res2the0102,res2the0112,res2the0202,&
			res2the0212,res2the1212,res2the0101square,&
			res2the0102square,res2the0112square,&
			res2the0202square,res2the0212square,&
			res2the1212square,res2the0101dsquare,&
			res2the0202dsquare,res2the1212dsquare,&
			res201num,res202num,res212num,res20101num,&
			res20101numbis,&
			res20102num,res20112num,&
			res20202num,res20212num,res21212num,&
			vet01,vet02,vet12)
			
			call fonct(t3(i),the01,ri01,gl01,su01)
			call fonct(t3(i),the02,ri02,gl02,su02)
			call fonct(t3(i),the12,ri12,gl12,su12)

			!write(6, *) "calculate integrals -c7" 
			
			iter = 0
			nweib = 0
			
			
			u1=(su12**vet12)*res2denum*ri12*vet12+&
			(su01**vet01)*(su02**vet02)*ri02*vet02
				
			
			if(fix(1).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=(LOG(the01(2))+(1/the01(1)))*res2denum +&
				res2thenum-LOG(the01(2))*(res2denum-res201num)-&
				res2the01
				!stop
				
				v=v*(su12**vet12)*ri12*vet12-&
				LOG(the01(2)*t3(i))*gl01*vet01*(su01**vet01)* &
				ri02*vet02*(su02**vet02)
				
				res1(nweib)=v/u1
				
				res1(nweib)=res1(nweib)+troncweib01011
				!ici

				
				res1((nvamax+iter))=&
				LOG(the01(2))*2/the01(1)*res2denum+&
				res2thenum*2/the01(1)+&
				(LOG(the01(2))**2)*res2denum+&
				LOG(the01(2))*2*res2thenum+&
				res2thenumsquare-3*(LOG(the01(2))**2)* &
				(res2denum-res201num)-&
				LOG(the01(2))*6*res2the01-3*res2the0101square-&
				LOG(the01(2))*2/the01(1)*(res2denum-res201num)- &
				res2the01*2/the01(1)+ &
				(LOG(the01(2))**2)*res20101numbis+&
				LOG(the01(2))*2*res2the0101+&
				res2the0101dsquare
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*(su12**vet12)*ri12*vet12+&
				gl01*vet01*(LOG(the01(2)*t3(i))**2)*(-1+&
				gl01*vet01)*(su01**vet01)*(su02**vet02)*ri02*vet02
				
				res1((nvamax+iter))=&
				(res1((nvamax+iter))*u1-v*v)/(u1*u1)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+&
				troncweib01011square
				!la
				if(fix(2).eq.0)then
					
					u2=res2denum*the01(1)/the01(2)-&
					(res2denum-res201num)*the01(1)/the01(2)
					u2=u2*(su12**vet12)*ri12*vet12-&
					ri02*vet02*the01(1)*gl01*vet01* &
					(su01**vet01)*(su02**vet02)/the01(2)
					
					u3=the01(1)*LOG(the01(2))*res20101numbis/the01(2)+&
					the01(1)*res2the0101/the01(2)+&
					(the01(1)*LOG(the01(2))+2)*res2denum/the01(2)+&
					the01(1)*res2thenum/the01(2)-&
					(3*the01(1)*LOG(the01(2))+2)* &
					(res2denum-res201num)/the01(2)-3* &
					the01(1)*res2the01/the01(2)
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri02*vet02*gl01*vet01*(su01**vet01)* &
					(su02**vet02)*(-1-&
					LOG(the01(2)*t3(i))*the01(1)+&
					gl01*vet01*LOG(the01(2)*t3(i))*the01(1))/the01(2)
					
					iter = iter +1
					res1((nvamax+iter))=&
					(u3*u1-v*u2)/(u1*u1)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+troncweib010112
				endif
				!la
				if(fix(3).eq.0)then
				
					u3=&
					-LOG(the01(2))*LOG(the02(2))*res202num-&
					LOG(the01(2))*res2the02-&
					LOG(the02(2))*res2the02-&
					res2the0202square-&
					LOG(the02(2))*res202num/the01(1)-&
					res2the02/the01(1)+&
					LOG(the01(2))*LOG(the02(2))*res20102num+&
					LOG(the01(2))*res2the0102+&
					LOG(the02(2))*res2the0102+&
					res2the0102square
					
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri02*vet02*gl01*vet01* &
					(su01**vet01)*(su02**vet02)* &
					LOG(the01(2)*t3(i))* &
					(LOG(the02(2)*t3(i))*gl02*vet02-&
					LOG(the02(2)*t3(i))-&
					(1/the02(1)))
					
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri02*vet02*(su01**vet01)*(su02**vet02)* &
					(LOG(the02(2)*t3(i))+(1/the02(1))-&
					LOG(the02(2)*t3(i))*gl02*vet02)
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!la
				if(fix(4).eq.0)then
				
					u3=-LOG(the01(2))*the02(1)*res202num/the02(2)-&
					the02(1)*res2the02/the02(2)-&
					the02(1)*res202num/(the01(1)*the02(2))+&
					LOG(the01(2))*the02(1)*res20102num/the02(2)+&
					the02(1)*res2the0102/the02(2)
					
					
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01* &
					(1-gl02*vet02)*(su01**vet01)* &
					(su02**vet02)*the02(1)* &
					LOG(the01(2)*t3(i))/the02(2)
					
					u2=-the02(1)*res202num/the02(2)
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri02*vet02*the02(1)*(1-gl02*vet02)*(su01**vet01)* &
					(su02**vet02)/the02(2)

					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!la
				if(fix(5).eq.0)then
				
					
					u3=LOG(the01(2))*LOG(the12(2))*res212num+&
					LOG(the01(2))*res2the12+LOG(the12(2))*res2the12+&
					res2the1212square+LOG(the12(2))*res212num/the01(1)+&
					res2the12/the01(1)-&
					LOG(the01(2))*LOG(the12(2))*res20112num-&
					LOG(the01(2))*res2the0112-&
					LOG(the12(2))*res2the0112-&
					res2the0112square
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)* &
					(LOG(the12(2)*t3(i))+ &
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)* &
					((LOG(the01(2))+(1/the01(1)))*res2denum +&
				res2thenum-LOG(the01(2))*(res2denum-res201num)-&
				res2the01)
					
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*res2denum* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)
						
						
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					v*u2)/(u1*u1)
					
				endif
				!la
				if(fix(6).eq.0)then
				
					u3=&
					LOG(the01(2))*the12(1)*res212num/the12(2)+&
					the12(1)/the12(2)*res2the12+res212num* &
					the12(1)/(the12(2)*the01(1))-&
					LOG(the01(2))*the12(1)*res20112num/the12(2)-&
					the12(1)/the12(2)*res2the0112
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*the12(1)* &
					(1-gl12*vet12)* &
					((LOG(the01(2))+(1/the01(1)))*res2denum +&
					res2thenum-LOG(the01(2))*(res2denum-res201num)-&
					res2the01)/the12(2)
					
					
					
					u2=the12(1)*res212num/the12(2)
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*the12(1)*(1-gl12*vet12)* &
					(su12**vet12)* &
					res2denum/the12(2)
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				
				!la
				if(nva01nofix.gt.0) then 
					
					u3=res2the0101+ &
					LOG(the01(2))*res20101numbis-&
					LOG(the01(2))*3*(res2denum-res201num)-&
					res2the01*3-(res2denum-res201num)/the01(1)+&
					(LOG(the01(2))+(1/the01(1)))*res2denum+&
					res2thenum
					
					u3= u3*(su12**vet12)*ri12*vet12-&
					ri02*vet02*LOG(the01(2)*t3(i))*gl01*vet01*(1-&
					gl01*vet01)*(su01**vet01)*(su02**vet02)
					
					u2=res201num*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01*(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01011beta01
					
					iter = iter +nva01nofix
				endif
				!la
				if(nva02nofix.gt.0) then 
					
					
					u3=-(LOG(the01(2))+&
					(1/the01(1)))*res202num-&
					res2the02+&
					LOG(the01(2))*res20102num+&
					res2the0102
					
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri02*vet02*LOG(the01(2)*t3(i))*gl01*vet01* &
					(su01**vet01)*(su02**vet02)* &
					(1-gl02*vet02)
					
					u2=-res202num*(su12**vet12)*ri12*vet12+&
					(1-gl02*vet02)*(su01**vet01)*(su02**vet02)* &
					ri02*vet02
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					iter = iter +nva02nofix
				endif
				!ici
				if(nva12nofix.gt.0) then 
				!la
					u3=(LOG(the01(2))+&
					(1/the01(1)))*res212num+&
					res2the12-&
					LOG(the01(2))*res20112num-&
					res2the0112
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1-gl12*vet12)*(su12**vet12)* &
					((LOG(the01(2))+(1/the01(1)))*res2denum +&
					res2thenum-LOG(the01(2))*(res2denum-res201num)-&
					res2the01)
					
					u2=res212num*(su12**vet12)*ri12*vet12+&
					(1-gl12*vet12)*(su12**vet12)*res2denum* &
					ri12*vet12
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
				endif
			endif
			
			!write(6, *) " done x101 -c7" 
			if(fix(2).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=the01(1)*res2denum/the01(2) -&
				the01(1)*(res2denum-res201num)/the01(2)
				
				v= v*(su12**vet12)*ri12*vet12-&
				ri02*vet02*the01(1)*gl01*vet01*(su01**vet01)* &
				(su02**vet02)/the01(2)
				
				res1(nweib)=v/u1
				res1(nweib)=res1(nweib)+troncweib01012
				!la
				
				u2=((the01(1)/the01(2))**2)*res20101numbis+&
					the01(1)*(the01(1)-1)*res2denum/((the01(2))**2)-&
					(res2denum-res201num)*(3*(the01(1)**2)-&
					the01(1))/(the01(2)**2)
					
				u2=u2*(su12**vet12)*ri12*vet12-&
				ri02*vet02*(1-gl01*vet01)*(su01**vet01)*(su02**vet02)* &
				((the01(1)/the01(2))**2)*gl01*vet01+&
				ri02*vet02*gl01*vet01*(su01**vet01)* &
				(su02**vet02)*the01(1)/(the01(2)**2)
				
				res1((nvamax+iter))=&
				(u2*u1-v*v)/(u1*u1)
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib01012square
				
				if(fix(3).eq.0)then
				!la!
					
					u3=-LOG(the02(2))*the01(1)* &
					res202num/the01(2)-&
					the01(1)*res2the02/the01(2)+&
					the01(1)*LOG(the02(2))*res20102num/the01(2)+&
					the01(1)*res2the0102/the01(2)
					
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01*the01(1)* &
					(su01**vet01)*(su02**vet02)* &
					(LOG(the02(2)*t3(i))*(1-&
					gl02*vet02)+(1/the02(1)))/the01(2)
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri02*vet02*(su01**vet01)*(su02**vet02)* &
					(LOG(the02(2)*t3(i))+&
					(1/the02(1))-&
					gl02*vet02*LOG(the02(2)*t3(i)))
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!la
				if(fix(4).eq.0)then
				
					u3=the01(1)*the02(1)*res20102num-&
					the01(1)*the02(1)*res202num
					u3=u3/(the01(2)*the02(2))
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri02*vet02*gl01*vet01*(gl02*vet02-1)* &
					(su01**vet01)*(su02**vet02)* &
					the01(1)*the02(1)/(the01(2)*the02(2))
					
					u2=-the02(1)*res202num/the02(2)
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri02*vet02*the02(1)*(1-gl02*vet02)* &
					(su01**vet01)* &
					(su02**vet02)/the02(2)
				
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!la
				if(fix(5).eq.0)then
				
					u3=LOG(the12(2))*the01(1)*res212num/the01(2)+&
					the01(1)*res2the12/the01(2)-&
					the01(1)*LOG(the12(2))*res20112num/the01(2)-&
					the01(1)*res2the0112/the01(2)
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*(LOG(the12(2)* &
					t3(i)) + &
					(1/the12(1))-LOG(the12(2)*t3(i))*gl12*vet12)* &
					(the01(1)*res2denum/the01(2) -&
					the01(1)*(res2denum-res201num)/the01(2))
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*res2denum* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)
						
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					v*u2)/(u1*u1)
					
				endif
				!la
				if(fix(6).eq.0)then
					
					u3=-the01(1)*the12(1)*res20112num+&
					the01(1)*the12(1)*res212num
					
					u3=u3/(the01(2)*the12(2))
					
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*the12(1)* &
					(the01(1)*res2denum/the01(2) -&
					the01(1)*(res2denum-res201num)/the01(2))* &
					(1-gl12*vet12)/the12(2)
					
					
					u2=the12(1)*res212num/the12(2)
					
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*the12(1)*(1-gl12*vet12)* &
					(su12**vet12)*res2denum/the12(2)
				
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!la
				if(nva01nofix.gt.0) then 
					
					u3=the01(1)*res20101numbis/the01(2)+&
					the01(1)*res2denum/the01(2)-&
					(3*the01(1)*(res2denum-res201num)/the01(2))
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri02*vet02*the01(1)*(1-gl01*vet01)* &
					gl01*vet01*(su01**vet01)*(su02**vet02)/the01(2)
					
					u2=res201num*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01*(su01**vet01)* &
					(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01012beta01
					
					iter = iter +nva01nofix
				endif
				!la
				if(nva02nofix.gt.0) then 
					
					u3=res20102num-res202num
					u3=u3*the01(1)/the01(2)
					u3=u3*(su12**vet12)*ri12*vet12-&
					the01(1)*gl01*vet01*(1-gl02*vet02)* &
					(su01**vet01)*(su02**vet02)/the01(2)* &
					ri02*vet02
					
					u2=-res202num*(su12**vet12)*ri12*vet12+&
					(1-gl02*vet02)*(su01**vet01)*(su02**vet02)* &
					ri02*vet02
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!la
					u3=-res20112num+res212num
					u3=u3*the01(1)/the01(2)
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1-gl12*vet12)*(su12**vet12)* &
					(the01(1)*res2denum/the01(2) -&
					the01(1)*(res2denum-res201num)/the01(2))
				
					
					u2=res212num*(su12**vet12)*ri12*vet12+&
					(1-gl12*vet12)*ri12*vet12*(su12**vet12)*res2denum
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					iter = iter +nva12nofix
				endif
			endif
			!write(6, *) " done x201 -c7" 
					
				!ici
		!la fin de journee : 181224
			if(fix(3).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=-LOG(the02(2))*res202num-&
				res2the02
				
				v=v*(su12**vet12)*ri12*vet12+&
				ri02*vet02*(su01**vet01)*(su02**vet02)* &
				(LOG(the02(2)*t3(i))+&
				(1/the02(1))-&
				gl02*vet02*LOG(the02(2)*t3(i)))
				
				res1(nweib)=v/u1
				res1(nweib)=res1(nweib)+troncweib02021
				res1((nvamax+iter))=&
				((LOG(the02(2)))**2)* &
				res20202num+2*LOG(the02(2))*res2the0202+&
				res2the0202dsquare-&
				((LOG(the02(2)))**2)* &
				res202num-2*LOG(the02(2))*res2the02-&
				res2the0202square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))* &
				(su12**vet12)*ri12*vet12+&
				ri02*vet02*(su01**vet01)*(su02**vet02)* &
				(((LOG(the02(2)*t3(i))+&
				(1/the02(1))-&
				gl02*vet02*LOG(the02(2)*t3(i)))**2)-&
				((1/(the02(1)**2))+&
				(LOG(the02(2)*t3(i))**2)*gl02*vet02))
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*u1 -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				res1((nvamax+iter))=&
				res1((nvamax+iter))+&
				troncweib02021square
				
				!la
				if(fix(4).eq.0)then
				
					u3=-(the02(1)*LOG(the02(2))+1)* &
					res202num/the02(2)-&
					the02(1)*res2the02/the02(2)+ &
					LOG(the02(2))*the02(1)* &
					res20202num/the02(2)+ &
					the02(1)/the02(2)*res2the0202
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri02*vet02*(1-gl02*vet02)*the02(1)* &
					(su01**vet01)*(su02**vet02)* &
					(LOG(the02(2)*t3(i))+&
					(1/the02(1))-&
					LOG(the02(2)*t3(i))*gl02*vet02)/the02(2)+ &
					ri02*vet02*(su01**vet01)*(su02**vet02)* &
					(1-gl02*vet02-LOG(the02(2)*t3(i))*the02(1)* &
					gl02*vet02)/the02(2)
					
					u2=-the02(1)*res202num/the02(2)
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri02*vet02*the02(1)*(1-gl02*vet02)* &
					(su01**vet01)* &
					(su02**vet02)/the02(2)
					
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+troncweib020212
				endif
				!la
				if(fix(5).eq.0)then
				
					u3=-LOG(the02(2))*LOG(the12(2))*res20212num-&
					LOG(the02(2))*res2the0212-LOG(the12(2))*res2the0212-&
					res2the0212square
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)* &
					(-1*LOG(the02(2))*res202num-&
					res2the02)
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*res2denum* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)
						
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
				endif
				!la
				if(fix(6).eq.0)then
				
					u3=-LOG(the02(2))*the12(1)*res20212num/the12(2)-&
					the12(1)*res2the0212/the12(2)
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1-gl12*vet12)*(su12**vet12)*the12(1)* &
					(-1*LOG(the02(2))*res202num-&
					res2the02)/the12(2)
				
					
					
					u2=the12(1)*res212num/the12(2)
					u2=u2*(su12**vet12)*ri12*vet12+ &
					ri12*vet12*the12(1)*(1-gl12*vet12)* &
					(su12**vet12)*res2denum/the12(2)
					
				
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
				endif
				!la
				if(nva01nofix.gt.0) then 
					
					
					u3=LOG(the02(2))*res20102num+&
					res2the0102-res202num*LOG(the02(2))-&
					res2the02
					
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01* &
					(su01**vet01)*(su02**vet02)* &
					(LOG(the02(2)*t3(i))+&
					(1/the02(1))-&
					LOG(the02(2)*t3(i))*gl02*vet02)
					
					u2=res201num*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01*(su01**vet01)* &
					(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
				endif
				!la
				if(nva02nofix.gt.0) then 
					
					u3=-LOG(the02(2))*res202num-&
					res2the02+LOG(the02(2))*res20202num+&
					res2the0202
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri02*vet02*(1-gl02*vet02)* &
					(su01**vet01)*(su02**vet02)* &
					(LOG(the02(2)*t3(i))+&
					(1/the02(1))-&
					LOG(the02(2)*t3(i))*gl02*vet02)-&
					ri02*vet02*gl02*vet02*(su01**vet01)* &
					(su02**vet02)*LOG(the02(2)*t3(i))
					
					u2=-res202num*(su12**vet12)*ri12*vet12+&
					(1-gl02*vet02)*(su01**vet01)* &
					(su02**vet02)*ri02*vet02
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02021beta02
					iter = iter +nva02nofix
					
				endif
				!la
				if(nva12nofix.gt.0) then 
				
					u3=-LOG(the02(2))*res20212num-&
					res2the0212
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1-gl12*vet12)*(su12**vet12)* &
					(-1*LOG(the02(2))*res202num-&
					res2the02)
					
					u2=res212num*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1-gl12*vet12)*(su12**vet12)*res2denum
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
				endif
			endif
			
             !write(6, *) " done x102 -c7" 
			!la
			if(fix(4).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=-the02(1)*res202num/the02(2)
				
				v=v*(su12**vet12)*ri12*vet12+&
				ri02*vet02*(1-gl02*vet02)* &
				the02(1)*(su01**vet01)* &
				(su02**vet02)/the02(2)
				
				res1(nweib)=v/u1
				res1(nweib)=res1(nweib)+troncweib02022
				!la
				
				u3=((the02(1)/the02(2))**2)*res20202num-&
					(the02(1)*(the02(1)-1)/((the02(2))**2))*res202num
				
				u3=u3*(su12**vet12)*ri12*vet12+&
				ri02*vet02*(1-gl02*vet02)* &
				(su01**vet01)*(su02**vet02)* &
				the02(1)*(-gl02*vet02*the02(1)/the02(2)- &
				(1/the02(2))+&
				the02(1)/the02(2))/the02(2)-&
				ri02*vet02*gl02*vet02*(the02(1)**2)* &
				(su01**vet01)*(su02**vet02)/(the02(2)**2)
				
				res1((nvamax+iter))=&
				u1*u3-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib02022square
				
				!write(6, *) " done x202/x202 -c7" 
				if(fix(5).eq.0)then
				
					u3=-the02(1)*LOG(the12(2))*res20212num/the02(2)-&
					the02(1)*res2the0212/the02(2)
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri12*vet12*(su12**vet12)* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))- &
					LOG(the12(2)*t3(i))*gl12*vet12)* &
					(the02(1)*res202num/the02(2))
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*res2denum* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)
						
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
					!write(6, *) " done x202/x112 -c7" 
					
				endif
				!la
				if(fix(6).eq.0)then
					
					u3=-the02(1)*the12(1)*res20212num
					u3=u3/(the02(2)*the12(2))
					
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri12*vet12*(su12**vet12)*the12(1)* &
					(1-gl12*vet12)* &
					(the02(1)*res202num/the02(2))/the12(2)
					
				
					u2=the12(1)*res212num/the12(2)
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*the12(1)*(1-gl12*vet12)* &
					(su12**vet12)*res2denum/the12(2)

					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
					!write(6, *) " done x202/x212 -c7" 
				endif
				
				!la
				if(nva01nofix.gt.0) then 
					
					u3=-the02(1)*res202num/the02(2)+&
					the02(1)*res20102num/the02(2)		
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri02*vet02*the02(1)*(1- &
					gl02*vet02)*gl01*vet01* &
					(su01**vet01)*(su02**vet02)/the02(2)
					
					u2=res201num*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01* &
					(su01**vet01)*(su02**vet02)
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
					
					!write(6, *) " done x202/b01 -c7" 
				endif
				!la
				if(nva02nofix.gt.0) then 
					
					u3=the02(1)*res20202num/the02(2)-&
					the02(1)*res202num/the02(2)
					u3=u3*(su12**vet12)*ri12*vet12+&
					the02(1)*((1-gl02*vet02)**2)*(su01**vet01)* &
					(su02**vet02)*ri02*vet02/the02(2)-&
					ri02*vet02*the02(1)*gl02*vet02* &
					(su01**vet01)*(su02**vet02)/the02(2)
					
					u2=-res202num*(su12**vet12)*ri12*vet12+&
					(1-gl02*vet02)*(su01**vet01)* &
					(su02**vet02)*ri02*vet02
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02022beta02
					
					iter = iter +nva02nofix
					!write(6, *) " done x202/b02 -c7" 
				endif
				!la
				if(nva12nofix.gt.0) then 
					u3=-res20212num*the02(1)/the02(2)
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri12*vet12*(1-gl12*vet12)*(su12**vet12)* &
					(the02(1)*res202num/the02(2))
					
					u2=res212num*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1- &
					gl12*vet12)*(su12**vet12)*res2denum
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					iter = iter +nva12nofix
					!write(6, *) " done x202/b12 -c7" 
				endif
			endif	
			!write(6, *) " done x202 -c7" 
			!la!
			if(fix(5).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=LOG(the12(2))*res212num+&
				res2the12
				v=v*(su12**vet12)*ri12*vet12+&
				ri12*vet12*(su12**vet12)*res2denum* &
				(LOG(the12(2)*t3(i))+&
				(1/the12(1))-&
				LOG(the12(2)*t3(i))*gl12*vet12)
				
				res1(nweib)=v/u1
				!la
				res1((nvamax+iter))=&
				((LOG(the12(2)))**2)* &
				res21212num+2*LOG(the12(2))*res2the1212+&
				res2the1212dsquare+&
				((LOG(the12(2)))**2)* &
				res212num+2*LOG(the12(2))*res2the12+&
				res2the1212square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*(su12**vet12)*ri12*vet12+&
				ri12*vet12*(su12**vet12)* &
				(LOG(the12(2)*t3(i)) +&
				(1/the12(1))-&
				LOG(the12(2)*t3(i))*gl12*vet12)* &
				(2*(LOG(the12(2))*res212num+&
				res2the12)+res2denum*(LOG(the12(2)*t3(i))+&
				(1/the12(1))-&
				LOG(the12(2)*t3(i))*gl12*vet12))-&
				ri12*vet12*(su12**vet12)*res2denum* &
				((1/(the12(1)**2))+&
				(LOG(the12(2)*t3(i))**2)*gl12*vet12)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*u1 -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				
				!write(6, *) " done x112/x112 -c7" 
				if(fix(6).eq.0)then
	
					u2=the12(1)*res212num/the12(2)
					
					u3=LOG(the12(2))*the12(1)*res21212num/the12(2)+&
					the12(1)*res2the1212/the12(2)+&
					(the12(1)*LOG(the12(2))+1)*res212num/the12(2)+&
					the12(1)*res2the12/the12(2)
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1-gl12*vet12)* &
					(su12**vet12)*the12(1)* &
					((LOG(the12(2))*res212num+&
					res2the12)+&
					res2denum*(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12))/the12(2)
					
					u3=u3+ &
					ri12*vet12*(su12**vet12)*res2denum* &
					(1-LOG(the12(2)*t3(i))*the12(1)* &
					gl12*vet12-gl12*vet12)/the12(2)+&
					ri12*vet12*(su12**vet12)*u2* &
					(LOG(the12(2)*t3(i))+(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)
					
					
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*the12(1)*(1-&
					gl12*vet12)*(su12**vet12)* &
					res2denum/the12(2)
					
					iter = iter +1
					
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
					!write(6, *) " done x112/x212 -c7" 
				endif
				!la
				if(nva01nofix.gt.0) then 
					!ici
					u3=-LOG(the12(2))*res20112num-&
					res2the0112+res212num*LOG(the12(2))+&
					res2the12
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)* &
					res201num
					
					
					u2=res201num*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01*(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
					
					!write(6, *) " done x112/b01 -c7" 
				endif
				if(nva02nofix.gt.0) then 
					!la
					u3=-LOG(the12(2))*res20212num-&
					res2the0212
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri12*vet12*(su12**vet12)* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)* &
					(res202num)
					
					u2=-res202num*(su12**vet12)*ri12*vet12+&
					(1-gl02*vet02)*(su01**vet01)* &
					(su02**vet02)*ri02*vet02
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					iter = iter +nva02nofix
					
					!write(6, *) " done x112/b02 -c7" 
					
				endif
				if(nva12nofix.gt.0) then 
					!la

					
					u3=LOG(the12(2))*res21212num+&
					res2the1212+LOG(the12(2))*res212num+&
					res2the12
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					(1-gl12*vet12)*(su12**vet12)* &
					(LOG(the12(2))*res212num+&
					res2the12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)* &
					(LOG(the12(2)*t3(i))+(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)* &
					(res212num+res2denum*(1-gl12*vet12))-&
					LOG(the12(2)*t3(i))*gl12*vet12* &
					ri12*vet12*(su12**vet12)*res2denum
					
					u2=res212num*(su12**vet12)*ri12*vet12+&
					(1-gl12*vet12)*ri12*vet12* &
					(su12**vet12)*res2denum
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
					
					!write(6, *) " done x112/b12 -c7" 
				endif
			endif
			
			  !la
			!write(6, *) " done x112 -c7" 
			if(fix(6).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=the12(1)*res212num/the12(2)
				v=v*(su12**vet12)*ri12*vet12+&
				ri12*vet12*the12(1)*(1- &
				gl12*vet12)*(su12**vet12)* &
				res2denum/the12(2)
				
				res1(nweib)=v/u1
				
				u3=((the12(1)/the12(2))**2)*res21212num+&
				(the12(1)*(the12(1)-1)/((the12(2))**2))*res212num
				

				u3=u3*(su12**vet12)*ri12*vet12+&
				ri12*vet12*(su12**vet12)*the12(1)*(1-&
				gl12*vet12)*(2*the12(1)*res212num/the12(2)+&
				res2denum*the12(1)*(1-gl12*vet12)/the12(2)- &
				res2denum/the12(2))/the12(2)-&
				(the12(1)**2)*gl12*vet12*ri12*vet12* &
				(su12**vet12)*res2denum/(the12(2)**2)
				
				
				res1((nvamax+iter))=&
				u1*u3-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				
				
				if(nva01nofix.gt.0) then 
					!la
					u3=the12(1)*res212num/the12(2)-&
					the12(1)*res20112num/the12(2)	
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*the12(1)*(1-&
					gl12*vet12)*res201num* &
					(su12**vet12)/the12(2)
					
					u2=res201num*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01* &
					(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
				endif
				
				if(nva02nofix.gt.0) then 
					!la
					
					
					u3=-the12(1)*res20212num/the12(2)
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri12*vet12*the12(1)*(1-&
					gl12*vet12)*(su12**vet12)* &
					res202num/the12(2)
					
					
					u2=-(su12**vet12)*ri12*vet12*res202num+&
				(1-gl02*vet02)*(su01**vet01)*(su02**vet02)*ri02*vet02

					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!la
					u3=res21212num*the12(1)/the12(2)+&
					the12(1)*res212num/the12(2)
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*(1-&
					gl12*vet12)*the12(1)*res212num/the12(2)+&
					ri12*vet12*(su12**vet12)*(1-&
					gl12*vet12)*the12(1)* &
					(res212num+res2denum*(1-gl12*vet12))/ &
					the12(2)-ri12*vet12*gl12*vet12* &
					(su12**vet12)*res2denum*the12(1)/ &
					the12(2)
					
					u2=res212num*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1-&
					gl12*vet12)*(su12**vet12)*res2denum
							
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
				endif
			endif	
			!write(6, *) " done x212 -c7" 
				v=(su12**vet12)*ri12*vet12*res2denum+&
				(su01**vet01)*(su02**vet02)*ri02*vet02

				if(nva01nofix.gt.0) then 

				u1=-gl01*vet01*(su01**vet01)
				u1=u1*(su02**vet02)*ri02*vet02
				u1=u1+(su12**vet12)*ri12*vet12*res201num

				
				res1((nvaweib+1):(nvaweib+nva01nofix))=&
				ve01nofix(i,:)*u1/v
				res1((nvaweib+1):(nvaweib+nva01nofix))= &
				res1((nvaweib+1):(nvaweib+nva01nofix))+tronc01

				res1((nvamax12weib12+1):nvamax01)=&
				res20101num*(su12**vet12)*ri12*vet12
				resint=(su01**vet01)*(su02**vet02)*gl01*vet01
				resint=resint*(1-gl01*vet01)*ri02*vet02
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)-resint
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)*ve01square(i,:)
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)/v
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)+tronc01square


				end if 

				if(nva02nofix.gt.0) then 

				u2=-(su12**vet12)*ri12*vet12*res202num+&
				(1-gl02*vet02)*(su01**vet01)*(su02**vet02)*ri02*vet02

				res1((nva01nofix+1+nvaweib):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1+nvaweib):nva0102)=&
				res1((nva01nofix+1+nvaweib):nva0102)+&
				tronc02

				res1((nvamax0112+1):nvamax02)=&
				-gl02*vet02*(su01**vet01)*(su02**vet02)*ri02*vet02
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)+&
				(su12**vet12)*ri12*vet12*(res20202num-res202num)+&
				((1-gl02*vet02)**2)*ri02*vet02*(su01**vet01)*(su02**vet02)


				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)*ve02square(i,:)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)/v
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)+tronc02square


				end if 

				if(nva12nofix.gt.0) then 

				u3=res212num+(1-gl12*vet12)*res2denum
				u3=u3*(su12**vet12)*ri12*vet12

				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v

				res1((nvamax0212+1):nvamax12)=&
				-gl12*vet12*res2denum+((1-gl12*vet12)**2)*res2denum+&
				(2*(1-gl12*vet12)*res212num+res212num)+&
				res21212num

				resint=ri12*vet12*(su12**vet12)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)*resint*ve12square(i,:)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)/v
				
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)


				end if 

				
				if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 

				kfix=nvamax01+1
				lfix=kfix-1+nva02nofix

				do j=1,nva01nofix
			    		res1(kfix:lfix)=(-res202num+&
					res20102num)*(su12**vet12)*ri12*vet12
					resint=-gl01*vet01*(1-gl02*vet02)*(su01**vet01)
					resint=resint*(su02**vet02)*ri02*vet02
					res1(kfix:lfix)=&
					res1(kfix:lfix)+resint
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			    		res1(kfix:lfix)=&
			    		res1(kfix:lfix)*v
					res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   	kfix=lfix+1
			   	lfix=lfix+nva02nofix
				end do

				end if 

				if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 

				
				kfix=nvamax0102+1
				lfix=kfix-1+nva12nofix

				do j=1,nva01nofix
			  		res1(kfix:lfix)=&
					(1-gl12*vet12)*res201num-res20112num+&
					res212num
					res1(kfix:lfix)=&
					res1(kfix:lfix)*ve01nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
					res1(kfix:lfix)*(su12**vet12)*ri12*vet12
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)*v
					res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva12nofix
				end do

				end if 


				if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 

				
				
        			kfix=nvamax02+1
				lfix=kfix-1+nva12nofix

				do j=1,nva02nofix
			  	res1(kfix:lfix)=res20212num+&
				(1-gl12*vet12)*res202num
				resint=(su12**vet12)
				res1(kfix:lfix)=&
				res1(kfix:lfix)*resint*ve02nofix(i,j)*ve12nofix(i,:)
				res1(kfix:lfix)=&
				-res1(kfix:lfix)*ri12*vet12
				res1(kfix:lfix)=&
				res1(kfix:lfix)*v
				res1(kfix:lfix)=&
				res1(kfix:lfix)-u2*u3*ve02nofix(i,j)*ve12nofix(i,:)
				res1(kfix:lfix)=&
			   	res1(kfix:lfix)/(v**2)
			   	kfix=lfix+1
			  	lfix=lfix+nva12nofix
				end do

				end if 
                            endif
                         endif                        
                      endif
                   endif  
                endif   
                endif   

                res = res + res1 


        end do   

!write(6, *) "end do" 
        likelihood_deriv = res




123     continue 
	 
	deallocate(b,bfix,fix,ve01,ve02,ve12,ve01nofix,&
	ve02nofix,ve12nofix,tronc01,tronc02,t0,t1,t2,t3,c,&
	ve01square,ve02square,ve12square,&
	tronc01square,tronc02square,troncweib01011beta01,&
	troncweib01012beta01,troncweib02021beta02,troncweib02022beta02)     

    end subroutine derivaweiballpara

	!=============================================================================================  
!======================= Calculate first derivatives of loglik with weibull baseline risk ==========
!=============================================================================================  


subroutine firstderivaweiballpara(b0,np0,npar0,bfix0,fix0,c0,no0,ve010,ve120,ve020,&
        dimnva01,dimnva12,dimnva02,nva01,nva12,nva02,t00,&
        t10,t20,t30,troncature0,likelihood_deriv)
	
		use commun
        implicit none
         
        double precision::res2denum,res201num,res202num,res212num, &
	res20101num,res20102num,res20112num,res20202num,res20212num,&
	res21212num,res20101numbis,&
	res2thenum,res2thenumsquare,res2the01,res2the02,res2the12,&
	res2the0101,res2the0202,res2the1212,res2the0102,res2the0112,&
	res2the0212,res2the0101square,res2the0202square,&
	res2the1212square,res2the0102square,res2the0112square,&
	res2the0212square,res2the0101dsquare,res2the0202dsquare,&
	res2the1212dsquare,vet01,vet12,vet02,resint,v,u1,u2,u3
	
        integer::np0,i,j,l,w,k,lfix, kfix,npar0,nva01,nva12,nva02,no0, &
	nz010,nz020,nz120,troncature0,dimnva01,dimnva02,dimnva12, & 
	nva01nofix,nva12nofix,nva02nofix,nvamax, sizespline,nva0102
	integer::nvamax01,nvamax0102,nvamax0112,nvamax02,nvamax0212,nvamax12, &
	nvaweib01,nvaweib02,nvaweib12,nvaweib,nvamax12weib12,iter,nweib

	double precision,dimension(np0+np0*(np0+1)/2),intent(inout)::likelihood_deriv
	double precision,dimension(np0)::b0
	double precision,dimension(np0+np0*(np0+1)/2)::res,res1
        double precision,dimension(npar0)::bh
	double precision,dimension(npar0-np0)::bfix0
	integer,dimension(npar0)::fix0
	double precision,dimension(2)::the01,the12,the02
        double precision,dimension(no0,dimnva01)::ve010
	double precision,dimension(no0,dimnva02)::ve020
	double precision,dimension(no0,dimnva12)::ve120
	
	
        double precision::su01,ri01,su12,ri12,su02,ri02,gl01,gl02,gl12,&
		troncweib01011,troncweib01011square, &
		troncweib010112,troncweib01012,&
		troncweib01012square,troncweib02021,troncweib02021square,&
		troncweib020212,troncweib02022,&
		troncweib02022square
	double precision,dimension(no0)::t00,t10,t20,t30
	integer,dimension(no0)::c0

	allocate(b(np0),bfix(npar0-np0),fix(npar0))
	b=b0
	bfix=bfix0
	fix=fix0
	troncature=troncature0

	sizespline=6


	nvaweib01=2-sum(fix(1:2))
	nvaweib02=2-sum(fix(3:4))
	nvaweib12=2-sum(fix(5:6))
	
	!write(6,*) "nvaweib01",nvaweib01 OK
	!write(6,*) "nvaweib02",nvaweib02 OK
	!write(6,*) "nvaweib12",nvaweib12 OK
	
	if(nva01.gt.0) then 
	  nva01nofix=nva01-sum(fix((sizespline+1):(nva01+sizespline)))
	else 
	  nva01nofix=0
	end if 

	if(nva02.gt.0) then 
          nva02nofix=nva02-sum(fix((nva01+sizespline+1):(nva02+nva01+sizespline)))
	else 
	   nva02nofix=0
	end if 

	if(nva12.gt.0) then 
	  nva12nofix=nva12-sum(fix((nva01+nva02+sizespline+1):npar0))
	else 
	  nva12nofix=0
	end if 

	nvaweib=nvaweib01+nvaweib02+nvaweib12
	nva0102=nva01nofix+nva02nofix+nvaweib
	nvamax=nva01nofix+nva02nofix+nva12nofix+nvaweib
	
	nvamax12weib12=nvamax+&
	(nvaweib+1)*nvaweib/2+&
	(nvamax-nvaweib)*nvaweib
	
	nvamax01=nvamax12weib12+(nva01nofix+1)*nva01nofix/2
	nvamax0102=nvamax01+nva01nofix*nva02nofix
	nvamax0112=nvamax0102+nva01nofix*nva12nofix
	
	nvamax02=nvamax0112+(nva02nofix+1)*nva02nofix/2
	nvamax0212=nvamax02+nva02nofix*nva12nofix
	nvamax12=nvamax0212+(nva12nofix+1)*nva12nofix/2


	if(nva01.gt.0) then 
		allocate(ve01(no0,nva01))
		allocate(ve01nofix(no0,nva01nofix))
		allocate(ve01square(no0,nva01nofix*(nva01nofix+1)/2))
		allocate(tronc01(nva01nofix))
		allocate(tronc01square(nva01nofix*(nva01nofix+1)/2))
		allocate(troncweib01011beta01(nva01nofix))
		allocate(troncweib01012beta01(nva01nofix))
	else 
		allocate(ve01(no0,1))
		allocate(ve01nofix(no0,1))
		ve01nofix=0
		allocate(ve01square(no0,1))
		ve01square=0
		allocate(tronc01(1))
		allocate(tronc01square(1))
		allocate(troncweib01011beta01(1))
		allocate(troncweib01012beta01(1))
	end if 
	
	if(nva02.gt.0) then 
		allocate(ve02(no0,nva02))
		allocate(ve02nofix(no0,nva02nofix))
		allocate(ve02square(no0,nva02nofix*(nva02nofix+1)/2))
		allocate(tronc02(nva02nofix))
		allocate(tronc02square(nva02nofix*(nva02nofix+1)/2))
		allocate(troncweib02021beta02(nva02nofix))
		allocate(troncweib02022beta02(nva02nofix))
	else 
		allocate(ve02(no0,1))
		allocate(ve02nofix(no0,1))
		ve02nofix=0
		allocate(ve02square(no0,1))
		ve02square=0
		allocate(tronc02(1))
		allocate(tronc02square(1))
		allocate(troncweib02021beta02(1))
		allocate(troncweib02022beta02(1))
	end if 

	if(nva12.gt.0) then 
		allocate(ve12(no0,nva12))
		allocate(ve12nofix(no0,nva12nofix))
		allocate(ve12square(no0,nva12nofix*(nva12nofix+1)/2))
	else 
		allocate(ve12(no0,1))
		allocate(ve12nofix(no0,1))
		ve12nofix=0
		allocate(ve12square(no0,1))
		ve12square=0
	end if 



	ve01=ve010
	ve02=ve020
	ve12=ve120

	allocate(t0(no0),t1(no0),t2(no0),t3(no0),c(no0))
	c=c0
	t0=t00
	t1=t10
	t2=t20
	t3=t30

         
        ! we need to put bh at its original values if in posfix 


       l=0
       lfix=0
       w=0

	if(nva01.gt.0) then
	do k=1,nva01
	   if(fix((sizespline+k)).eq.0) then 
		lfix=lfix+1
		ve01nofix(:,lfix)=ve01(:,k)
	   end if 
	end do
	
	do i=1,no0
		lfix=1
		do l=1,nva01nofix
	   	ve01square(i,lfix:(lfix+nva01nofix-l))=&
		ve01nofix(i,l:nva01nofix)
		ve01square(i,lfix:(lfix+nva01nofix-l))=&
		ve01square(i,lfix:(lfix+nva01nofix-l))*ve01nofix(i,l)
		lfix=lfix+nva01nofix-l+1
		end do
	end do
	 
        end if 
	lfix=0

	if(nva02.gt.0) then 
	do k=1,nva02
	   if(fix((sizespline+nva01+k)).eq.0) then 
		lfix=lfix+1
		ve02nofix(:,lfix)=ve02(:,k)
	   end if 
	end do

	do i=1,no0
		lfix=1
		do l=1,nva02nofix
	   	ve02square(i,lfix:(lfix+nva02nofix-l))=&
		ve02nofix(i,l)*ve02nofix(i,l:nva02nofix)
		lfix=lfix+nva02nofix-l+1
		end do
	end do

	end if 
	
	lfix=0

	if(nva12.gt.0) then 
	do k=1,nva12
	   if(fix((sizespline+nva01+nva02+k)).eq.0) then 
		lfix=lfix+1
		ve12nofix(:,lfix)=ve12(:,k)
	   end if 
	end do

	do i=1,no0
		lfix=1
		do l=1,nva12nofix
	   		ve12square(i,lfix:(lfix+nva12nofix-l))=&
			ve12nofix(i,l)*ve12nofix(i,l:nva12nofix)
			lfix=lfix+nva12nofix-l+1
		end do
	end do

	end if
	
	l=0
       lfix=0
       w=0

     do k=1,npar0 
         if(fix(k).eq.0) then
            l=l+1
            bh(k)=b(l)
	 end if 
         if(fix(k).eq.1) then
            w=w+1
            bh(k)=bfix(w)
         end if
      end do

   
	


	
  
    !the01(2)=dexp(bh(2))
	!the02(2)=dexp(bh(4))
	!the12(2)=dexp(bh(6))
    
	the01(2)=(bh(2))**2
	the02(2)=(bh(4))**2
	the12(2)=(bh(6))**2
	
	the01(1)=(bh(1))**2
	the02(1)=(bh(3))**2
	the12(1)=(bh(5))**2

	!write(6, *) "before loop" OK

!---------- calcul des derivees premiere ------------------   

	res = 0
        do i=1,no0

!write(6, *) "subject ",i
                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
				vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
				vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
				vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

		
                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)

                res1 = 0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                            tronc01 = 0
							tronc02 =  0
                        	tronc01square= 0
                        	tronc02square=0
							troncweib01011=0
							troncweib01011square=0
							troncweib01011beta01=0
							troncweib01012beta01=0
							troncweib010112=0
							troncweib01012=0
							troncweib01012square=0
							troncweib02021=0
							troncweib02021square=0
							troncweib020212=0
							troncweib02021beta02=0
							troncweib02022=0
							troncweib02022square=0
							troncweib02022beta02=0
                        else 
				call fonct(t0(i),the01,ri01,gl01,su01)
				call fonct(t0(i),the02,ri02,gl02,su02)
				
                            tronc01=ve01nofix(i,:)*gl01*vet01
                        	tronc02=ve02nofix(i,:)*gl02*vet02
                        	tronc01square=ve01square(i,:)*gl01*vet01
                        	tronc02square=ve02square(i,:)*gl02*vet02
							
							troncweib01011=gl01*vet01*LOG(the01(2)*t0(i))
							troncweib01011square=&
							((LOG(the01(2)*t0(i)))**2)*gl01*vet01
							
							troncweib01011beta01=LOG(the01(2)*t0(i))*tronc01
							
							troncweib01012beta01=tronc01*the01(1)/the01(2)
							
							troncweib010112=&
							(1+LOG(the01(2)*t0(i))*the01(1))*gl01*vet01/the01(2)
							
							troncweib01012=the01(1)*gl01*vet01/the01(2)
							
							troncweib01012square=&
							the01(1)*(the01(1)-1)*gl01*vet01/(the01(2)**2)
							
							troncweib02021=LOG(the02(2)*t0(i))*gl02*vet02
							troncweib02021square=&
							((LOG(the02(2)*t0(i)))**2)*gl02*vet02
							
							troncweib020212=&
							(1+LOG(the02(2)*t0(i))*the02(1))*gl02*vet02/the02(2)
							
							troncweib02021beta02=&
							tronc02*LOG(the02(2)*t0(i))
							
							troncweib02022=&
							the02(1)*gl02*vet02/the02(2)
							
							troncweib02022square=&
							the02(1)*(the02(1)-1)*gl02*vet02/(the02(2)**2)
							
							troncweib02022beta02=&
							tronc02*the02(1)/the02(2)
							
							
							
                        end if
                else
                    tronc01 = 0
							tronc02 =  0
                        	tronc01square= 0
                        	tronc02square=0
							troncweib01011=0
							troncweib01011square=0
							troncweib01011beta01=0
							troncweib01012beta01=0
							troncweib010112=0
							troncweib01012=0
							troncweib01012square=0
							troncweib02021=0
							troncweib02021square=0
							troncweib020212=0
							troncweib02021beta02=0
							troncweib02022=0
							troncweib02022square=0
							troncweib02022beta02=0
                end if
		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
				!write(6, *) "start c=1" 
			call fonct(t1(i),the01,ri01,gl01,su01)
			call fonct(t1(i),the02,ri02,gl02,su02)
			iter = 0
			nweib= 0
			if(fix(1).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				res1(nweib)=-LOG(the01(2)*t1(i))*gl01*vet01 +&
				troncweib01011
				res1((nvamax+iter))=&
				-((LOG(the01(2)*t1(i)))**2)*gl01*vet01 +&
				troncweib01011square
				if(fix(2).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=&
					-(1+the01(1)*LOG(the01(2)*t1(i)))* &
					gl01*vet01/the01(2) +&
					troncweib010112
				endif
				if(fix(3).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					-ve01nofix(i,:)*gl01*vet01*LOG(the01(2)*t1(i))+&
					troncweib01011beta01
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif 
			
			if(fix(2).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				res1(nweib)=-the01(1)*gl01*vet01/the01(2) +&
				troncweib01012
				
				res1(nvamax+iter)=&
				-the01(1)*(the01(1)-1)*gl01*vet01/(the01(2)**2)+&
				troncweib01012square
				
				if(fix(3).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					-ve01nofix(i,:)*gl01*vet01*the01(1)/the01(2)+&
					troncweib01012beta01
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif
			
			if(fix(3).eq.0)then
				iter = iter + 1
				nweib = nweib +1
				res1(nweib)=-LOG(the02(2)*t1(i))*gl02*vet02 +&
				troncweib02021
				res1((nvamax+iter))=&
				-((LOG(the02(2)*t1(i)))**2)*gl02*vet02+&
				troncweib02021square
				
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=&
					-(1+the02(1)*LOG(the02(2)*t1(i)))* &
					gl02*vet02/the02(2)+&
					troncweib020212
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					troncweib02021beta02-&
					ve02nofix(i,:)*gl02*vet02*LOG(the02(2)*t1(i))
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif 
			
			if(fix(4).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				res1(nweib)=-the02(1)*gl02*vet02/the02(2)+&
				troncweib02022
				res1(nvamax+iter)=&
				-the02(1)*(the02(1)-1)*gl02*vet02/(the02(2)**2)+&
				troncweib02022square
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					-ve02nofix(i,:)*gl02*vet02*the02(1)/the02(2)+&
					troncweib02022beta02
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif
			
			
			if(fix(5).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				res1(nweib)=0
				res1(nvamax+iter)=0
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif
			
			if(fix(6).eq.0)then
				iter = iter + 1
				nweib = nweib +1
				res1(nweib)=0
				res1(nvamax+iter)=0
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif

			if(nva01nofix.gt.0) then 
			
			res1((nvaweib+1):(nva01nofix+nvaweib))=&
			-ve01nofix(i,:)*gl01*vet01+tronc01
			res1((nvamax12weib12+1):nvamax01)=&
			-ve01square(i,:)*gl01*vet01+&
			tronc01square

			if(nva02nofix.gt.0) then 
			res1((nvamax01+1):nvamax0102)=0
			end if 
			if(nva12nofix.gt.0) then 
			res1((nvamax0102+1):nvamax0112)=0
			end if

			end if 

			if(nva02nofix.gt.0) then
			res1((nvaweib+nva01nofix+1):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02

			res1((nvamax0112+1):nvamax02)=&
			-ve02square(i,:)*gl02*vet02+&
			tronc02square
			if(nva12nofix.gt.0) then 
			res1((nvamax02+1):nvamax0212)=0
			end if 
			end if 

			if(nva12nofix.gt.0) then 
			res1((nva0102+1):nvamax)=0
			res1((nvamax0212+1):nvamax12)=0
			end if 
			
			
			

                else
                if(c(i).eq.2)then ! cpi 0-->1
			!write(6, *) "start c=2" 

			call fonct(t3(i),the12,ri12,gl12,su12)
			!call logapproxi(t1(i),t2(i),res2denum) approxi log(x) OK
			call qgaussweibbetaderiv(t1(i),t2(i),the01,&
			the02,the12,res2denum,res2the01,res2the02,&
			res2the12,res2thenum,res2thenumsquare,&
			res2the0101,res2the0102,res2the0112,res2the0202,&
			res2the0212,res2the1212,res2the0101square,&
			res2the0102square,res2the0112square,&
			res2the0202square,res2the0212square,&
			res2the1212square,res2the0101dsquare,&
			res2the0202dsquare,res2the1212dsquare,&
			res201num,res202num,res212num,res20101num,&
			res20101numbis,&
			res20102num,res20112num,&
			res20202num,res20212num,res21212num,&
			vet01,vet02,vet12)
			
			iter = 0
			nweib = 0
			!write(6,*) "res2the01",res2the01
			if(fix(1).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				!write(6,*) "nweib",nweib
				!write(6,*) "res2denum theta",res2denum
				v=(LOG(the01(2))+(1/the01(1)))*res2denum +&
				res2thenum-LOG(the01(2))*(res2denum-res201num)-&
				res2the01
				res1(nweib)=v/res2denum
				
				res1(nweib)=res1(nweib)+troncweib01011
				
				!write(6,*) "troncweib01011square",troncweib01011square
				res1((nvamax+iter))=&
				(LOG(the01(2))**2)*res2denum+ &
				LOG(the01(2))*2*res2thenum+ &
				res2thenumsquare+ &
				LOG(the01(2))*2/the01(1)*res2denum+ &
				res2thenum*2/the01(1)- &
				(LOG(the01(2))**2)*3*(res2denum-res201num)- &
				LOG(the01(2))*6*res2the01- &
				res2the0101square*3- &
				LOG(the01(2))*2*(res2denum-res201num)/the01(1)-&
				res2the01*2/the01(1)+ &
				(LOG(the01(2))**2)*res20101numbis+ &
				LOG(the01(2))*2*res2the0101+ &
				res2the0101dsquare
				
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*res2denum -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib01011square
				
				if(fix(2).eq.0)then
					
					u2=res2denum*the01(1)/the01(2)-&
					(res2denum-res201num)*the01(1)/the01(2)
					
					u1=the01(1)*LOG(the01(2))*res20101numbis/the01(2)+&
					the01(1)*res2the0101/the01(2)+&
					(the01(1)*LOG(the01(2))+2)*res2denum/the01(2)+&
					the01(1)*res2thenum/the01(2)-&
					(3*the01(1)*LOG(the01(2))+2)* &
					(res2denum-res201num)/the01(2)-3* &
					the01(1)*res2the01/the01(2)
					
					iter = iter +1
					res1((nvamax+iter))=&
					u1*res2denum-u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+troncweib010112
				endif
				if(fix(3).eq.0)then
				
					u1=-LOG(the01(2))*LOG(the02(2))*res202num-&
					LOG(the01(2))*res2the02-LOG(the02(2))*res2the02-&
					res2the0202square-LOG(the02(2))*res202num/the01(1)-&
					res2the02/the01(1)+&
					LOG(the01(2))*LOG(the02(2))*res20102num+&
					LOG(the01(2))*res2the0102+&
					LOG(the02(2))*res2the0102+&
					res2the0102square
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(fix(4).eq.0)then
				
					u1=-LOG(the01(2))*the02(1)*res202num/the02(2)-&
					the02(1)*res2the02/the02(2)-&
					the02(1)*res202num/(the01(1)*the02(2))+&
					LOG(the01(2))*the02(1)*res20102num/the02(2)+&
					the02(1)*res2the0102/the02(2)
					
					u2=-the02(1)*res202num/the02(2)
					
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(fix(5).eq.0)then
				
					u1=LOG(the01(2))*LOG(the12(2))*res212num+&
					LOG(the01(2))*res2the12+LOG(the12(2))*res2the12+&
					res2the1212square+LOG(the12(2))*res212num/the01(1)+&
					res2the12/the01(1)-&
					LOG(the01(2))*LOG(the12(2))*res20112num-&
					LOG(the01(2))*res2the0112-&
					LOG(the12(2))*res2the0112-&
					res2the0112square
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
					
					u1=&
					LOG(the01(2))*the12(1)*res212num/the12(2)+&
					the12(1)/the12(2)*res2the12+res212num* &
					the12(1)/(the12(2)*the01(1))-&
					LOG(the01(2))*the12(1)*res20112num/the12(2)-&
					the12(1)/the12(2)*res2the0112
					
					
					u2=the12(1)*res212num/the12(2)
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(nva01nofix.gt.0) then 
					
					u1=res2the0101+ &
					LOG(the01(2))*res20101numbis-&
					LOG(the01(2))*3*(res2denum-res201num)-&
					res2the01*3-(res2denum-res201num)/the01(1)+&
					(LOG(the01(2))+(1/the01(1)))*res2denum+&
					res2thenum
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01011beta01
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					u1=-(LOG(the01(2))+&
					(1/the01(1)))*res202num-&
					res2the02+&
					LOG(the01(2))*res20102num+&
					res2the0102
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					u1=(LOG(the01(2))+&
					(1/the01(1)))*res212num+&
					res2the12-&
					LOG(the01(2))*res20112num-&
					res2the0112
					
					u2=res212num
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva12nofix
				endif
			endif
			
			!write(6, *) " done x101 -c2" 
			if(fix(2).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=the01(1)*res2denum/the01(2) -&
				the01(1)*(res2denum-res201num)/the01(2)
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+troncweib01012
				
				u1=((the01(1)/the01(2))**2)*res20101numbis+&
					the01(1)*(the01(1)-1)*res2denum/((the01(2))**2)-&
					(res2denum-res201num)*(3*(the01(1)**2)-&
					the01(1))/(the01(2)**2)
				
				
				res1((nvamax+iter))=&
				u1*res2denum-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+ troncweib01012square
				if(fix(3).eq.0)then
				
					u1=-LOG(the02(2))*the01(1)* &
					res202num/the01(2)-&
					the01(1)*res2the02/the01(2)+&
					the01(1)*LOG(the02(2))*res20102num/the01(2)+&
					the01(1)*res2the0102/the01(2)
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(fix(4).eq.0)then
				
					u1=the01(1)*the02(1)*res20102num-&
					the01(1)*the02(1)*res202num
					u1=u1/(the01(2)*the02(2))
					
					u2=-the02(1)*res202num/the02(2)
				
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(fix(5).eq.0)then
				
					u1=LOG(the12(2))*the01(1)*res212num/the01(2)+&
					the01(1)*res2the12/the01(2)-&
					the01(1)*LOG(the12(2))*res20112num/the01(2)-&
					the01(1)*res2the0112/the01(2)
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
					
					u1=-the01(1)*the12(1)*res20112num+&
					the01(1)*the12(1)*res212num
					u1=u1/(the01(2)*the12(2))
					
					u2=the12(1)*res212num/the12(2)
				
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(nva01nofix.gt.0) then 
					
					u1=the01(1)*res20101numbis/the01(2)+&
					the01(1)*res2denum/the01(2)-&
					the01(1)*3* &
					(res2denum-res201num)/the01(2)
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01012beta01
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					u1=res20102num-res202num
					u1=u1*the01(1)/the01(2)
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!here 
					u1=-res20112num+res212num
					u1=u1*the01(1)/the01(2)
					
					u2=res212num
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva12nofix
				endif
			endif
			
				!write(6, *) " done x201 -c2" 	
				!here 
				
			if(fix(3).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				!write(6,*) "res2the02",res2the02
				v=-LOG(the02(2))*res202num-&
				res2the02
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+ &
				troncweib02021
				
				res1((nvamax+iter))=&
				((LOG(the02(2)))**2)* &
				res20202num+2*LOG(the02(2))*res2the0202+&
				res2the0202dsquare-&
				((LOG(the02(2)))**2)* &
				res202num-2*LOG(the02(2))*res2the02-&
				res2the0202square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*res2denum -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				res1((nvamax+iter))=&
				res1((nvamax+iter))+&
				troncweib02021square
				
				if(fix(4).eq.0)then
				
					u1=-(the02(1)*LOG(the02(2))+1)* &
					res202num/the02(2)-&
					the02(1)*res2the02/the02(2)+ &
					LOG(the02(2))*the02(1)* &
					res20202num/the02(2)+ &
					the02(1)/the02(2)*res2the0202
					
					u2=-the02(1)*res202num/the02(2)
					
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+&
					troncweib020212
				endif
				
				if(fix(5).eq.0)then
				
					u1=-LOG(the02(2))*LOG(the12(2))*res20212num-&
					LOG(the02(2))*res2the0212-LOG(the12(2))*res2the0212-&
					res2the0212square
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
				
					u1=-LOG(the02(2))*the12(1)*res20212num/the12(2)-&
					the12(1)*res2the0212/the12(2)
					
					u2=the12(1)*res212num/the12(2)
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(nva01nofix.gt.0) then 
					
					u1=LOG(the02(2))*res20102num+&
					res2the0102-res202num*LOG(the02(2))-&
					res2the02
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					u1=-LOG(the02(2))*res202num-&
					res2the02+LOG(the02(2))*res20202num+&
					res2the0202
					
					u2=-res202num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02021beta02
					iter = iter +nva02nofix
					
				endif
				if(nva12nofix.gt.0) then 
				
					u1=-LOG(the02(2))*res20212num-&
					res2the0212
					
					u2=res212num
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva12nofix
				endif
			endif
			
             
			!write(6, *) " done x102 -c2" 	
			if(fix(4).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=-the02(1)*res202num/the02(2)
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+troncweib02022
				
				u1=((the02(1)/the02(2))**2)*res20202num-&
					(the02(1)*(the02(1)-1)/((the02(2))**2))*res202num
				
				
				res1((nvamax+iter))=&
				u1*res2denum-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib02022square
				
				
				if(fix(5).eq.0)then
				
					u1=-the02(1)*LOG(the12(2))*res20212num/the02(2)-&
					the02(1)*res2the0212/the02(2)
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
					
					u1=-the02(1)*the12(1)*res20212num
					u1=u1/(the02(2)*the12(2))
					
					u2=the12(1)*res212num/the12(2)
				
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				
				if(nva01nofix.gt.0) then 
					
					u1=-the02(1)*res202num/the02(2)+&
					the02(1)*res20102num/the02(2)					
					
					u2=res201num
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				
				if(nva02nofix.gt.0) then 
					
					u1=the02(1)*res20202num/the02(2)-&
					the02(1)*res202num/the02(2)
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02022beta02
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					u1=-res20212num*the02(1)/the02(2)
					
					u2=res212num
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva12nofix
				endif
			endif	
			!write(6, *) " done x202 -c2" 	
			if(fix(5).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=LOG(the12(2))*res212num+&
				res2the12
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)-&
				LOG(the12(2)*t3(i))*gl12*vet12
				
				res1((nvamax+iter))=&
				((LOG(the12(2)))**2)* &
				res21212num+2*LOG(the12(2))*res2the1212+&
				res2the1212dsquare+&
				((LOG(the12(2)))**2)* &
				res212num+2*LOG(the12(2))*res2the12+&
				res2the1212square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*res2denum -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				res1((nvamax+iter))=&
				res1((nvamax+iter))-&
				((LOG(the12(2)*t3(i)))**2)*gl12*vet12
				

				if(fix(6).eq.0)then
			
					u1=LOG(the12(2))*the12(1)*res21212num/the12(2)+&
					the12(1)*res2the1212/the12(2)+&
					(the12(1)*LOG(the12(2))+1)*res212num/the12(2)+&
					the12(1)*res2the12/the12(2)
					
					u2=the12(1)*res212num/the12(2)
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					res1((nvamax+iter))=&
					res1((nvamax+iter))-&
					(1+LOG(the12(2)*t3(i))*the12(1))*gl12*vet12/the12(2)
				endif
				
				if(nva01nofix.gt.0) then 
					
					u1=-LOG(the12(2))*res20112num-&
					res2the0112+res212num*LOG(the12(2))+&
					res2the12
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					!here
					u1=-LOG(the12(2))*res20212num-&
					res2the0212
					
					u2=-res202num
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva02nofix
					
				endif
				if(nva12nofix.gt.0) then 
				
					u1=LOG(the12(2))*res21212num+&
					res2the1212+LOG(the12(2))*res212num+&
					res2the12
					
					u2=res212num
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*LOG(the12(2)*t3(i))*gl12*vet12
					
					iter = iter +nva12nofix
				endif
			endif
			
			  !!write(6, *) " done x112 -c2" 	

			if(fix(6).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=the12(1)*res212num/the12(2)
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)-&
				the12(1)*gl12*vet12/the12(2)
				
				u1=((the12(1)/the12(2))**2)*res21212num+&
					(the12(1)*(the12(1)-1)/((the12(2))**2))*res212num
				
				
				res1((nvamax+iter))=&
				u1*res2denum-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))-the12(1)* &
				(the12(1)-1)*gl12*vet12/(the12(2)**2)
				
				
				
				if(nva01nofix.gt.0) then 
					
					u1=the12(1)*res212num/the12(2)-&
					the12(1)*res20112num/the12(2)					
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				
				if(nva02nofix.gt.0) then 
					
					u1=-the12(1)*res20212num/the12(2)
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!here
					u1=res21212num*the12(1)/the12(2)+&
					the12(1)*res212num/the12(2)
					
					u2=res212num
							
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*the12(1)*gl12*vet12/the12(2)
					iter = iter +nva12nofix
				endif
			endif	
			
			write(6,*) "res1(1) after all",res1(1) 	
			v=res2denum*(su12**vet12)
			!write(6,*)"v",v
			!write(6,*)"nvaweib",nvaweib
			!write(6,*) "nvamax",nvamax
			!write(6,*) "nvamax01",nvamax01
			!write(6,*) "nvamax02",nvamax02
			!write(6,*) "nvamax0112",nvamax0112
			!write(6,*) "nvamax0212",nvamax0212
			!write(6,*) "nvamax12",nvamax12
			!write(6,*) "nva0102",nva0102
			!write(6,*) "iter",iter
			!write(6,*) "nvamax12weib12",nvamax12weib12
			
			if(nva01nofix.gt.0) then

			u1=res201num*(su12**vet12)
			
      		res1((nvaweib+1):(nvaweib+nva01nofix))=&
			ve01nofix(i,:)*u1/v
			res1((nvaweib+1):(nvaweib+nva01nofix))=&
			res1((nvaweib+1):(nvaweib+nva01nofix))+tronc01

			res1((nvamax12weib12+1):nvamax01)=&
			ve01square(i,:)*res20101num*(su12**vet12)
			res1((nvamax12weib12+1):nvamax01)=&
			res1((nvamax12weib12+1):nvamax01)/v
			res1((nvamax12weib12+1):nvamax01)=&
			res1((nvamax12weib12+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
			res1((nvamax12weib12+1):nvamax01)=&
			res1((nvamax12weib12+1):nvamax01)+&
			tronc01square

			end if 

			
			if(nva02nofix.gt.0) then

			u2=-res202num*(su12**vet12)
			res1((nvaweib+nva01nofix+1):nva0102)=&
			ve02nofix(i,:)*u2/v
			res1((nvaweib+nva01nofix+1):nva0102)=&
			res1((nvaweib+nva01nofix+1):nva0102)+&
			tronc02

			res1((nvamax0112+1):nvamax02)=&
			-ve02square(i,:)*(su12**vet12)*(-res20202num+res202num)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)/v
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)+tronc02square


			end if 

			
			if(nva12nofix.gt.0) then 

			u3=res212num-gl12*vet12*res2denum
			u3=u3*(su12**vet12)

			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v
			res1((nvamax0212+1):nvamax12)=&
			-2*gl12*vet12*res212num+res212num+&
			res21212num-gl12*vet12*(1-gl12*vet12)*res2denum
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*(su12**vet12)*ve12square(i,:)
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)/v
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)

			
			end if 
			
			

			
			if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 

			kfix=nvamax01+1
			lfix=kfix-1+nva02nofix

			do j=1,nva01nofix

			   res1(kfix:lfix)=&
			  (res20102num-res202num)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva02nofix

			end do

			end if 

			if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 

			kfix=nvamax0102+1
			lfix=kfix-1+nva12nofix

			do j=1,nva01nofix
			   res1(kfix:lfix)=-gl12*vet12*res201num-res20112num+&
			   res212num
			   res1(kfix:lfix)=res1(kfix:lfix)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do

			end if 

			
			if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 

			kfix=nvamax02+1
			lfix=kfix-1+nva12nofix

			do j=1,nva02nofix
			   res1(kfix:lfix)=-res20212num+&
			   gl12*vet12*res202num
			   res1(kfix:lfix)=res1(kfix:lfix)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve02nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u2*u3*ve02nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do

			end if 
			
                else  
                    if(c(i).eq.3)then ! obs 0-->1
					!here 
			!write(6, *) "start c=3" 
			call fonct(t1(i),the01,ri01,gl01,su01)
			call fonct(t1(i),the02,ri02,gl02,su02)
			call fonct(t1(i),the12,ri12,gl12,su12)

			iter = 0
			nweib= 0
			if(fix(1).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				res1(nweib)=-LOG(the01(2)*t1(i))*gl01*vet01 +&
				troncweib01011+&
				LOG(the01(2)*t1(i))+(1/the01(1))
				
				res1((nvamax+iter))=&
				-((LOG(the01(2)*t1(i)))**2)*gl01*vet01 +&
				troncweib01011square-1/(the01(1)**2)
				!here
				if(fix(2).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=&
					-(LOG(the01(2)*t1(i))* &
					the01(1)+1)*gl01*vet01/the01(2)+&
					troncweib010112+(1/the01(2))
				endif
				if(fix(3).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				!here
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					-ve01nofix(i,:)*gl01*vet01*LOG(the01(2)*t1(i))+&
					troncweib01011beta01
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif 
			
			if(fix(2).eq.0)then
			
				iter = iter + 1
				nweib = nweib + 1
				
				res1(nweib)=-the01(1)*gl01*vet01/the01(2) +&
				troncweib01012+&
				the01(1)/the01(2)

				res1(nvamax+iter)=&
				-the01(1)*(the01(1)-1)*gl01*vet01/(the01(2)**2)+&
				troncweib01012square-&
				the01(1)/(the01(2)**2)
				
				if(fix(3).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
				
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					-ve01nofix(i,:)*gl01*vet01*the01(1)/the01(2)+&
					troncweib01012beta01
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif
			!here
			if(fix(3).eq.0)then
				iter = iter + 1
				nweib = nweib +1
				res1(nweib)=-LOG(the02(2)*t1(i))*gl02*vet02 +&
				troncweib02021
				
				res1((nvamax+iter))=&
				-((LOG(the02(2)*t1(i)))**2)*gl02*vet02+&
				troncweib02021square
				!here
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=&
					-(LOG(the02(2)*t1(i))* &
					the02(1)+1)*gl02*vet02/the02(2)+&
					troncweib020212
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter):(nvamax+iter+nva02nofix))=&
					-ve02nofix(i,:)*LOG(the02(2)*t1(i))*gl12*vet12+&
					troncweib02021beta02
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif 
			!here
			if(fix(4).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				
				res1(nweib)=-the02(1)*gl02*vet02/the02(2)+&
				troncweib02022
				
				res1(nvamax+iter)=&
				-the02(1)*(the02(1)-1)*gl02*vet02/(the02(2)**2)+&
				troncweib02022square
				!here
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					-ve02nofix(i,:)*gl02*vet02*the02(1)/the02(2)+&
					troncweib02022beta02
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif
			!here
			
			if(fix(5).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				res1(nweib)=LOG(the12(2)*t1(i))*gl12*vet12
				res1(nvamax+iter)=((LOG(the12(2)*t1(i)))**2)*gl12*vet12
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=(LOG(the12(2)*t1(i))* &
					the12(1)+1)*gl12*vet12/the12(2)
				endif
				if(nva01nofix.gt.0) then 
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*LOG(the12(2)*t1(i))*gl12*vet12
					call fonct(t3(i),the12,ri12,gl12,su12)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*LOG(the12(2)*t3(i))*gl12*vet12
					
					iter = iter +nva12nofix
				endif
				
				call fonct(t3(i),the12,ri12,gl12,su12)
				
				res1(nweib)=res1(nweib)-&
				LOG(the12(2)*t3(i))*gl12*vet12
				iter=iter-nva12nofix-nva02nofix-nva01nofix
				if(fix(6).eq.0)then
					res1((nvamax+iter))=&
					res1((nvamax+iter))-&
					(LOG(the12(2)*t3(i))* &
					the12(1)+1)*gl12*vet12/the12(2)
					iter=iter-1
				endif
				
				res1(nvamax+iter)=res1(nvamax+iter)-&
				((LOG(the12(2)*t3(i)))**2)*gl12*vet12
				
				iter=iter+nva12nofix+nva02nofix+nva01nofix
				
				if(fix(6).eq.0)then
					iter=iter+1
				endif
				call fonct(t1(i),the12,ri12,gl12,su12)
			endif
			!here
			
			if(fix(6).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				res1(nweib)=the12(1)*gl12*vet12/the12(2)
				res1(nvamax+iter)=the12(1)* &
				(the12(1)-1)*gl12*vet12/(the12(2)**2)
				
				if(nva01nofix.gt.0) then 
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*the12(1)*gl12*vet12/the12(2)
					call fonct(t3(i),the12,ri12,gl12,su12)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*the12(1)*gl12*vet12/the12(2)
					
					iter = iter +nva12nofix
				endif
				
				call fonct(t3(i),the12,ri12,gl12,su12)
				
				res1(nweib)=res1(nweib)-&
				the12(1)*gl12*vet12/the12(2)
				
				iter=iter-nva12nofix-nva02nofix-nva01nofix
				
				res1(nvamax+iter)=res1(nvamax+iter)-&
				the12(1)*(the12(1)-1)*gl12*vet12/(the12(2)**2)
				
				iter=iter+nva12nofix+nva02nofix+nva01nofix
				
				call fonct(t1(i),the12,ri12,gl12,su12)
			endif
			
			
!finish
			if(nva01nofix.gt.0) then 

			res1((nvaweib+1):(nvaweib+nva01nofix))=-ve01nofix(i,:)*gl01*vet01+&
			tronc01+ve01nofix(i,:)
			res1((nvamax12weib12+1):nvamax01)=&
			-ve01square(i,:)*gl01*vet01+&
			tronc01square

			if(nva02nofix.gt.0) then 
			res1((nvamax01+1):nvamax0102)=0
			end if 

			if(nva12nofix.gt.0)then 
			res1((nvamax0102+1):nvamax0112)=0
			end if 

			end if 

			if(nva02nofix.gt.0) then 

			res1((nva01nofix+1+nvaweib):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02
			res1((nvamax0112+1):nvamax02)=&
			-ve02square(i,:)*gl02*vet02+&
			tronc02square

			if(nva12nofix.gt.0) then 
			res1((nvamax02+1):nvamax0212)=0
			end if 

			end if 

			if(nva12nofix.gt.0) then 

			res1((nva0102+1):nvamax)=gl12*vet12
			res1((nvamax0212+1):nvamax12)=gl12*vet12

			call fonct(t3(i),the12,ri12,gl12,su12)
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)-gl12*vet12
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)*ve12nofix(i,:)

			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)-gl12*vet12
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*ve12square(i,:)
			

			end if 
			
			

                    else   !here!17124
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			
			!write(6, *) "start c=4" 
			call fonct(t3(i),the12,ri12,gl12,su12)
			call qgaussweibbetaderiv(t1(i),t2(i),the01,&
			the02,the12,res2denum,res2the01,res2the02,&
			res2the12,res2thenum,res2thenumsquare,&
			res2the0101,res2the0102,res2the0112,res2the0202,&
			res2the0212,res2the1212,res2the0101square,&
			res2the0102square,res2the0112square,&
			res2the0202square,res2the0212square,&
			res2the1212square,res2the0101dsquare,&
			res2the0202dsquare,res2the1212dsquare,&
			res201num,res202num,res212num,res20101num,&
			res20101numbis,&
			res20102num,res20112num,&
			res20202num,res20212num,res21212num,&
			vet01,vet02,vet12)
			
			iter = 0
			nweib = 0
			
			
			if(fix(1).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=(LOG(the01(2))+(1/the01(1)))*res2denum +&
				res2thenum-LOG(the01(2))*(res2denum-res201num)-&
				res2the01
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+troncweib01011
				
				res1((nvamax+iter))=&
				LOG(the01(2))*2/the01(1)*res2denum+&
				res2thenum*2/the01(1)+&
				(LOG(the01(2))**2)*res2denum+&
				LOG(the01(2))*2*res2thenum+&
				res2thenumsquare-3*(LOG(the01(2))**2)* &
				(res2denum-res201num)-&
				LOG(the01(2))*6*res2the01-3*res2the0101square-&
				LOG(the01(2))*2/the01(1)*(res2denum-res201num)- &
				res2the01*2/the01(1)+ &
				(LOG(the01(2))**2)*res20101numbis+&
				LOG(the01(2))*2*res2the0101+&
				res2the0101dsquare
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*res2denum -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+&
				troncweib01011square
				
				if(fix(2).eq.0)then
					
					u2=res2denum*the01(1)/the01(2)-&
					(res2denum-res201num)*the01(1)/the01(2)
					
					u1=the01(1)*LOG(the01(2))*res20101numbis/the01(2)+&
					the01(1)*res2the0101/the01(2)+&
					(the01(1)*LOG(the01(2))+2)*res2denum/the01(2)+&
					the01(1)*res2thenum/the01(2)-&
					(3*the01(1)*LOG(the01(2))+2)* &
					(res2denum-res201num)/the01(2)-3* &
					the01(1)*res2the01/the01(2)
					
					iter = iter +1
					res1((nvamax+iter))=&
					u1*res2denum-u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+&
					troncweib010112
				endif
				if(fix(3).eq.0)then
				!here
					u1=-LOG(the01(2))*LOG(the02(2))*res202num-&
					LOG(the01(2))*res2the02-LOG(the02(2))*res2the02-&
					res2the0202square-LOG(the02(2))*res202num/the01(1)-&
					res2the02/the01(1)+&
					LOG(the01(2))*LOG(the02(2))*res20102num+&
					LOG(the01(2))*res2the0102+&
					LOG(the02(2))*res2the0102+&
					res2the0102square
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(fix(4).eq.0)then
				
					u1=-LOG(the01(2))*the02(1)*res202num/the02(2)-&
					the02(1)*res2the02/the02(2)-&
					the02(1)*res202num/(the01(1)*the02(2))+&
					LOG(the01(2))*the02(1)*res20102num/the02(2)+&
					the02(1)*res2the0102/the02(2)
					
					u2=-the02(1)*res202num/the02(2)
					
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(fix(5).eq.0)then
				!here
					u1=LOG(the01(2))*LOG(the12(2))*res212num+&
					LOG(the01(2))*res2the12+LOG(the12(2))*res2the12+&
					res2the1212square+LOG(the12(2))*res212num/the01(1)+&
					res2the12/the01(1)-&
					LOG(the01(2))*LOG(the12(2))*res20112num-&
					LOG(the01(2))*res2the0112-&
					LOG(the12(2))*res2the0112-&
					res2the0112square
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
					u1=&
					LOG(the01(2))*the12(1)*res212num/the12(2)+&
					the12(1)/the12(2)*res2the12+res212num* &
					the12(1)/(the12(2)*the01(1))-&
					LOG(the01(2))*the12(1)*res20112num/the12(2)-&
					the12(1)/the12(2)*res2the0112
					
					
					u2=the12(1)*res212num/the12(2)
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(nva01nofix.gt.0) then 
					!here
					u1=res2the0101+ &
					LOG(the01(2))*res20101numbis-&
					LOG(the01(2))*3*(res2denum-res201num)-&
					res2the01*3-(res2denum-res201num)/the01(1)+&
					(LOG(the01(2))+(1/the01(1)))*res2denum+&
					res2thenum
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01011beta01
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					u1=-(LOG(the01(2))+&
					(1/the01(1)))*res202num-&
					res2the02+&
					LOG(the01(2))*res20102num+&
					res2the0102
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					u1=(LOG(the01(2))+&
					(1/the01(1)))*res212num+&
					res2the12-&
					LOG(the01(2))*res20112num-&
					res2the0112
					
					u2=res212num
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva12nofix
				endif
			endif
			
			!write(6, *) " done x101 -c4" 	
			if(fix(2).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=the01(1)*res2denum/the01(2) -&
				the01(1)*(res2denum-res201num)/the01(2)
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+troncweib01012
				
				
					
				u1=((the01(1)/the01(2))**2)*res20101numbis+&
					the01(1)*(the01(1)-1)*res2denum/((the01(2))**2)-&
					(res2denum-res201num)*(3*(the01(1)**2)-&
					the01(1))/(the01(2)**2)
				
				res1((nvamax+iter))=&
				u1*res2denum-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib01012square
				if(fix(3).eq.0)then
				
					u1=-LOG(the02(2))*the01(1)*res202num/the01(2)-&
					the01(1)*res2the02/the01(2)+&
					the01(1)*LOG(the02(2))*res20102num/the01(2)+&
					the01(1)*res2the0102/the01(2)
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				!here
				if(fix(4).eq.0)then
				
					u1=the01(1)*the02(1)*res20102num-&
					the01(1)*the02(1)*res202num
					u1=u1/(the01(2)*the02(2))
					
					u2=-the02(1)*res202num/the02(2)
				
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(fix(5).eq.0)then
				
					u1=LOG(the12(2))*the01(1)*res212num/the01(2)+&
					the01(1)*res2the12/the01(2)-&
					the01(1)*LOG(the12(2))*res20112num/the01(2)-&
					the01(1)*res2the0112/the01(2)
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
					
					u1=-the01(1)*the12(1)*res20112num+&
					the01(1)*the12(1)*res212num
					u1=u1/(the01(2)*the12(2))
					
					u2=the12(1)*res212num/the12(2)
				
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				!here
				if(nva01nofix.gt.0) then 
					
					u1=the01(1)*res20101numbis/the01(2)+&
					the01(1)*res2denum/the01(2)-3* &
					the01(1)*(res2denum-res201num)/the01(2)
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01012beta01
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					u1=res20102num-res202num
					u1=u1*the01(1)/the01(2)
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!here 
					u1=-res20112num+res212num
					u1=u1*the01(1)/the01(2)
					
					u2=res212num
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva12nofix
				endif
			endif
			!write(6, *) " done x201 -c4" 	
					
				!!here 
				
			if(fix(3).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=-LOG(the02(2))*res202num-&
				res2the02
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+troncweib02021
				
				res1((nvamax+iter))=&
				((LOG(the02(2)))**2)*res20202num+2* &
				LOG(the02(2))*res2the0202+&
				res2the0202dsquare-&
				((LOG(the02(2)))**2)*res202num-2* &
				LOG(the02(2))*res2the02-&
				res2the0202square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*res2denum -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				res1((nvamax+iter))=&
				res1((nvamax+iter))+&
				troncweib02021square
				
				if(fix(4).eq.0)then
				
					u1=-(the02(1)*LOG(the02(2))+1)* &
					res202num/the02(2)-&
					the02(1)*res2the02/the02(2)+ &
					LOG(the02(2))*the02(1)* &
					res20202num/the02(2)+ &
					the02(1)/the02(2)*res2the0202
					
					u2=-the02(1)*res202num/the02(2)
					
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+&
					troncweib020212
				endif
				
				if(fix(5).eq.0)then
				
					u1=-LOG(the02(2))*LOG(the12(2))*res20212num-&
					LOG(the02(2))*res2the0212-LOG(the12(2))*res2the0212-&
					res2the0212square
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
				
					u1=-LOG(the02(2))*the12(1)*res20212num/the12(2)-&
					the12(1)*res2the0212/the12(2)
					
					u2=the12(1)*res212num/the12(2)
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				if(nva01nofix.gt.0) then 
					
					u1=LOG(the02(2))*res20102num+&
					res2the0102-res202num*LOG(the02(2))-&
					res2the02
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					
					u1=-LOG(the02(2))*res202num-&
					res2the02+LOG(the02(2))*res20202num+&
					res2the0202
					
					u2=-res202num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02021beta02
					iter = iter +nva02nofix
					
				endif
				if(nva12nofix.gt.0) then 
				
					u1=-LOG(the02(2))*res20212num-&
					res2the0212
					
					u2=res212num
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva12nofix
				endif
			endif
			
             !write(6, *) " done x102 -c4" 	
!here!
			if(fix(4).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=-the02(1)*res202num/the02(2)
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+troncweib02022
				
				
				u1=((the02(1)/the02(2))**2)*res20202num-&
					(the02(1)*(the02(1)-1)/((the02(2))**2))*res202num
				
				
				res1((nvamax+iter))=&
				u1*res2denum-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib02022square
				
				
				if(fix(5).eq.0)then
				
					u1=-the02(1)*LOG(the12(2))*res20212num/the02(2)-&
					the02(1)*res2the0212/the02(2)
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
						
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					
				endif
				
				if(fix(6).eq.0)then
					
					u1=-the02(1)*the12(1)*res20212num
					u1=u1/(the02(2)*the12(2))
					
					u2=the12(1)*res212num/the12(2)
				
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
				endif
				
				
				if(nva01nofix.gt.0) then 
					
					u1=-the02(1)*res202num/the02(2)+&
					the02(1)*res20102num/the02(2)					
					
					u2=res201num
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				
				if(nva02nofix.gt.0) then 
					
					u1=the02(1)*res20202num/the02(2)-&
					the02(1)*res202num/the02(2)
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02022beta02
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					u1=-res20212num*the02(1)/the02(2)
					
					u2=res212num
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva12nofix
				endif
			endif	
!here!
!write(6, *) " done x202 -c4" 	
			if(fix(5).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=LOG(the12(2))*res212num+&
				res2the12
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)+&
				LOG(the12(2)*t3(i))*(1-gl12*vet12)+1/ &
				the12(1)
				
				res1((nvamax+iter))=&
				((LOG(the12(2)))**2)*res21212num+2* &
				LOG(the12(2))*res2the1212+&
				res2the1212dsquare+&
				((LOG(the12(2)))**2)*res212num+2* &
				LOG(the12(2))*res2the12+&
				res2the1212square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*res2denum -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				res1((nvamax+iter))=&
				res1((nvamax+iter))-&
				((LOG(the12(2)*t3(i)))**2)*gl12*vet12-&
				(1/(the12(1)**2))
				
!here!
				if(fix(6).eq.0)then
			
					u1=LOG(the12(2))*the12(1)*res21212num/the12(2)+&
					the12(1)*res2the1212/the12(2)+&
					(the12(1)*LOG(the12(2))+1)*res212num/the12(2)+&
					the12(1)*res2the12/the12(2)
					
					u2=the12(1)*res212num/the12(2)
					iter = iter +1
					res1((nvamax+iter))=u1*res2denum-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(res2denum*res2denum)
					res1((nvamax+iter))=&
					res1((nvamax+iter))-&
					(1+LOG(the12(2)*t3(i))*the12(1))*gl12*vet12/the12(2)+&
					(1/the12(2))
				endif
				
				if(nva01nofix.gt.0) then 
					
					u1=-LOG(the12(2))*res20112num-&
					res2the0112+res212num*LOG(the12(2))+&
					res2the12
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					!here
					u1=-LOG(the12(2))*res20212num-&
					res2the0212
					
					u2=-res202num
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					iter = iter +nva02nofix
					
				endif
				!here
				if(nva12nofix.gt.0) then 
				
					u1=LOG(the12(2))*res21212num+&
					res2the1212+LOG(the12(2))*res212num+&
					res2the12
					
					u2=res212num
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*LOG(the12(2)*t3(i))*gl12*vet12
					
					iter = iter +nva12nofix
				endif
			endif
			
	!write(6, *) " done x112 -c4" 			  
!here!
			if(fix(6).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=the12(1)*res212num/the12(2)
				res1(nweib)=v/res2denum
				res1(nweib)=res1(nweib)-&
				the12(1)*gl12*vet12/the12(2)+&
				the12(1)/the12(2)
				
				u1=((the12(1)/the12(2))**2)*res21212num+&
					(the12(1)*(the12(1)-1)/((the12(2))**2))*res212num
				
				
				res1((nvamax+iter))=&
				u1*res2denum-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(res2denum*res2denum)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))-the12(1)*(the12(1)-1)* &
				gl12*vet12/(the12(2)**2)-&
				the12(1)/(the12(2)**2)
				
				
				
				if(nva01nofix.gt.0) then 
					
					u1=the12(1)*res212num/the12(2)-&
					the12(1)*res20112num/the12(2)					
					
					u2=res201num
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva01nofix
				endif
				
				if(nva02nofix.gt.0) then 
					
					u1=-the12(1)*res20212num/the12(2)
					
					u2=-res202num
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!here
					u1=res21212num*the12(1)/the12(2)+&
					the12(1)*res212num/the12(2)
					
					u2=res212num
							
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*res2denum-&
					v*u2)/(res2denum*res2denum)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*the12(1)*gl12*vet12/the12(2)
					iter = iter +nva12nofix
				endif
			endif	
!here!
!write(6, *) " done x212 -c4" 	
			v=res2denum*(su12**vet12)*ri12*vet12

			if(nva01nofix.gt.0) then 
			u1=res201num*(su12**vet12)*ri12*vet12
			
			res1((nvaweib+1):nva01nofix)=&
			ve01nofix(i,:)*u1/v
			res1((nvaweib+1):nva01nofix)=res1((nvaweib+1):nva01nofix)+tronc01

			res1((nvamax12weib12+1):nvamax01)=&
			res20101num*(su12**vet12)*ri12*vet12*ve01square(i,:)
			res1((nvamax12weib12+1):nvamax01)=&
			res1((nvamax12weib12+1):nvamax01)/v
			res1((nvamax12weib12+1):nvamax01)=&
			res1((nvamax12weib12+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
			res1((nvamax12weib12+1):nvamax01)=&
			res1((nvamax12weib12+1):nvamax01)+tronc01square

			end if 

			
			if(nva02nofix.gt.0) then 

			u2=-res202num*(su12**vet12)*ri12*vet12
			
			res1((nva01nofix+1+nvaweib):nva0102)=&
			ve02nofix(i,:)*u2/v
			res1((nva01nofix+1+nvaweib):nva0102)=&
			res1((nva01nofix+1+nvaweib):nva0102)+&
			tronc02
			
			res1((nvamax0112+1):nvamax02)=&
			(res20202num-res202num)*ri12*vet12*(su12**vet12)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)*ve02square(i,:)	
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)/v
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)+tronc02square

			end if 


			
			if(nva12nofix.gt.0) then 

			u3=res212num+&
			(1-gl12*vet12)*res2denum
			u3=u3*(su12**vet12)*ri12*vet12
			
			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v
			
			res1((nvamax0212+1):nvamax12)=&
			(1-gl12*vet12)*res212num+res21212num+&
			((1-gl12*vet12)**2)*res2denum+&
			(2-gl12*vet12)*res212num-gl12*vet12*res2denum
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*(su12**vet12)*ri12*vet12
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*ve12square(i,:)
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)/v
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)			


			end if 

			
			
			if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 
			kfix=nvamax01+1
			lfix=kfix-1+nva02nofix

			do j=1,nva01nofix
			   res1(kfix:lfix)=(-res202num+res20102num)*ri12*vet12
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva02nofix
			end do

			end if 

			if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 
			
			kfix=nvamax0102+1
			lfix=kfix-1+nva12nofix

			do j=1,nva01nofix
			   res1(kfix:lfix)=&
			   (1-gl12*vet12)*res201num+&
			   res212num-res20112num
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve12nofix(i,:)
		           res1(kfix:lfix)=&
			   res1(kfix:lfix)*ri12*vet12*(su12**vet12)
			   res1(kfix:lfix)=res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do


			end if 

			if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 
			
			kfix=nvamax02+1
			lfix=kfix-1+nva12nofix

			do j=1,nva02nofix
			   res1(kfix:lfix)=&
			   res20212num+(1-gl12*vet12)*res202num
			   res1(kfix:lfix)=&
			   -res1(kfix:lfix)*ve02nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ri12*vet12*(su12**vet12)
			   res1(kfix:lfix)=res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-ve02nofix(i,j)*ve12nofix(i,:)*u2*u3
			   res1(kfix:lfix)=res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do

			end if 

!here!
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				!write(6, *) "start c=5" 
				call fonct(t1(i),the01,ri01,gl01,su01)
				call fonct(t1(i),the02,ri02,gl02,su02)
				call fonct(t1(i),the12,ri12,gl12,su12)
				
				iter = 0
				nweib= 0
			if(fix(1).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				res1(nweib)=-LOG(the01(2)*t1(i))*gl01*vet01 +&
				troncweib01011+&
				LOG(the01(2)*t1(i))+(1/the01(1))
				
				res1((nvamax+iter))=&
				-((LOG(the01(2)*t1(i)))**2)*gl01*vet01 +&
				troncweib01011square-1/(the01(1)**2)
				!here!!
				if(fix(2).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=&
					-(LOG(the01(2)*t1(i))*the01(1)+1)* &
					gl01*vet01/the01(2)+&
					troncweib010112+(1/the01(2))
				endif
				if(fix(3).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				!here!!
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					-ve01nofix(i,:)*gl01*vet01*LOG(the01(2)*t1(i))+&
					troncweib01011beta01
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif 
			!here!!
			if(fix(2).eq.0)then
			
				iter = iter + 1
				nweib = nweib + 1
				
				res1(nweib)=-the01(1)*gl01*vet01/the01(2) +&
				troncweib01012+&
				the01(1)/the01(2)

				res1(nvamax+iter)=&
				-the01(1)*(the01(1)-1)*gl01*vet01/(the01(2)**2)+&
				troncweib01012square-&
				the01(1)/(the01(2)**2)
				
				if(fix(3).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
				
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					-ve01nofix(i,:)*gl01*vet01*the01(1)/the01(2)+&
					troncweib01012beta01
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif
			!here!!
			if(fix(3).eq.0)then
				iter = iter + 1
				nweib = nweib +1
				res1(nweib)=-LOG(the02(2)*t1(i))*gl02*vet02 +&
				troncweib02021
				
				res1((nvamax+iter))=&
				-((LOG(the02(2)*t1(i)))**2)*gl02*vet02+&
				troncweib02021square
				!here!!
				if(fix(4).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=&
					-(LOG(the02(2)*t1(i))*the02(1)+1)* &
					gl02*vet02/the02(2)+&
					troncweib02021
				endif
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter):(nvamax+iter+nva02nofix))=&
					-ve02nofix(i,:)*LOG(the02(2)*t1(i))*gl12*vet12+&
					troncweib02021beta02
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif 
			!here!!
			if(fix(4).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				
				res1(nweib)=-the02(1)*gl02*vet02/the02(2)+&
				troncweib02022
				
				res1(nvamax+iter)=&
				-the02(1)*(the02(1)-1)*gl02*vet02/(the02(2)**2)+&
				troncweib02022square
				!here!!
				if(fix(5).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=0
				endif
				if(nva01nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					-ve02nofix(i,:)*gl02*vet02*the02(1)/the02(2)+&
					troncweib02022beta02
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=0.d0
					iter = iter +nva12nofix
				endif
			endif
			!here!!
			
			if(fix(5).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				res1(nweib)=LOG(the12(2)*t1(i))*gl12*vet12+&
				(1/the12(1))
				res1(nvamax+iter)=((LOG(the12(2)*t1(i)))**2)*gl12*vet12-&
				(1/(the12(1)**2))
				if(fix(6).eq.0)then
					iter = iter +1
					res1((nvamax+iter))=&
					(LOG(the12(2)*t1(i))*the12(1)+1)* &
					gl12*vet12/the12(2)+(1/the12(2))
				endif
				if(nva01nofix.gt.0) then 
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*LOG(the12(2)*t1(i))*gl12*vet12
					call fonct(t3(i),the12,ri12,gl12,su12)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*LOG(the12(2)*t3(i))*gl12*vet12
					
					iter = iter +nva12nofix
				endif
				
				call fonct(t3(i),the12,ri12,gl12,su12)
				
				res1(nweib)=res1(nweib)-&
				LOG(the12(2)*t3(i))*gl12*vet12+&
				LOG(the12(2)*t3(i))
				iter=iter-nva12nofix-nva02nofix-nva01nofix
				if(fix(6).eq.0)then
					res1((nvamax+iter))=res1((nvamax+iter))-&
					(LOG(the12(2)*t3(i))*the12(1)+1)* &
					gl12*vet12/the12(2)
					iter=iter-1
				endif
				
				res1(nvamax+iter)=res1(nvamax+iter)-&
				((LOG(the12(2)*t3(i)))**2)*gl12*vet12
				
				iter=iter+nva12nofix+nva02nofix+nva01nofix
				
				if(fix(6).eq.0)then
					iter=iter+1
				endif
				call fonct(t1(i),the12,ri12,gl12,su12)
			endif
			!here!!
			
			if(fix(6).eq.0)then
				iter = iter + 1
				nweib= nweib +1
				res1(nweib)=the12(1)*gl12*vet12/the12(2)+&
				the12(1)/the12(2)
				res1(nvamax+iter)=the12(1)*(-1+&
				the12(1))*gl12*vet12/(the12(2)**2)-&
				(the12(1)/(the12(2)**2))
				
				if(nva01nofix.gt.0) then 
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=0.d0
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=0.d0
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*the12(1)*gl12*vet12/the12(2)
					call fonct(t3(i),the12,ri12,gl12,su12)
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))-&
					ve12nofix(i,:)*the12(1)*gl12*vet12/the12(2)
					
					iter = iter +nva12nofix
				endif
				
				call fonct(t3(i),the12,ri12,gl12,su12)
				
				res1(nweib)=res1(nweib)-&
				the12(1)*gl12*vet12/the12(2)
				
				iter=iter-nva12nofix-nva02nofix-nva01nofix
				
				res1(nvamax+iter)=res1(nvamax+iter)-&
				the12(1)*(the12(1)-1)*gl12*vet12/(the12(2)**2)
				
				iter=iter+nva12nofix+nva02nofix+nva01nofix
				
				call fonct(t1(i),the12,ri12,gl12,su12)
			endif
			
			
				if(nva01nofix.gt.0) then 

				res1((nvaweib+1):(nva01nofix+nvaweib))=&
				-ve01nofix(i,:)*gl01*vet01+&
				ve01nofix(i,:)+&
				tronc01
				res1((nvamax12weib12+1):nvamax01)=&
				-ve01square(i,:)*gl01*vet01+&
				tronc01square

				if(nva02nofix.gt.0) then 
				res1((nvamax01+1):nvamax0102)=0
				end if 

				if(nva12nofix.gt.0) then
				res1((nvamax0102+1):nvamax0112)=0
				end if 

				end if 

				if(nva02nofix.gt.0) then 

				res1((nva01nofix+1+nvaweib):nva0102)=&
				-ve02nofix(i,:)*gl02*vet02+tronc02
				

				res1((nvamax0112+1):nvamax02)=&
				-ve02square(i,:)*gl02*vet02+&
				tronc02square

				if(nva12nofix.gt.0) then
				res1((nvamax02+1):nvamax0212)=0
				end if 

				end if 
			
				if(nva12nofix.gt.0) then 

				res1((nvamax0212+1):nvamax12)=&
				gl12*vet12
				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*gl12*vet12+&
				ve12nofix(i,:)

				call fonct(t3(i),the12,ri12,gl12,su12)
				res1((nva0102+1):nvamax)=&
				res1((nva0102+1):nvamax)-ve12nofix(i,:)*gl12*vet12
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)-gl12*vet12
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)*ve12square(i,:)
				

				end if 
				
				
				!here!!
                         else
                            if(c(i).eq.6)then ! vivant ???

				!write(6, *) "start c=6" 
				call fonct(t3(i),the01,ri01,gl01,su01)
				call fonct(t3(i),the02,ri02,gl02,su02)
				call fonct(t3(i),the12,ri12,gl12,su12)

			call qgaussweibbetaderiv(t1(i),t3(i),the01,&
			the02,the12,res2denum,res2the01,res2the02,&
			res2the12,res2thenum,res2thenumsquare,&
			res2the0101,res2the0102,res2the0112,res2the0202,&
			res2the0212,res2the1212,res2the0101square,&
			res2the0102square,res2the0112square,&
			res2the0202square,res2the0212square,&
			res2the1212square,res2the0101dsquare,&
			res2the0202dsquare,res2the1212dsquare,&
			res201num,res202num,res212num,res20101num,&
			res20101numbis,&
			res20102num,res20112num,&
			res20202num,res20212num,res21212num,&
			vet01,vet02,vet12)
			
			iter = 0
			nweib = 0
			u1=((su12**vet12)*res2denum+(su01**vet01)*(su02**vet02))
				
				
				
			if(fix(1).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=(LOG(the01(2))+(1/the01(1)))*res2denum +&
				res2thenum-LOG(the01(2))*(res2denum-res201num)-&
				res2the01
				
				v=v*(su12**vet12)-&
				LOG(the01(2)*t3(i))*gl01*vet01*(su01**vet01)*(su02**vet02)
				
				res1(nweib)=v/u1
				
				res1(nweib)=res1(nweib)+troncweib01011
				
				res1((nvamax+iter))=&
				LOG(the01(2))*2/the01(1)*res2denum+&
				res2thenum*2/the01(1)+&
				(LOG(the01(2))**2)*res2denum+&
				LOG(the01(2))*2*res2thenum+&
				res2thenumsquare-3*(LOG(the01(2))**2)* &
				(res2denum-res201num)-&
				LOG(the01(2))*6*res2the01-3*res2the0101square-&
				LOG(the01(2))*2/the01(1)*(res2denum-res201num)- &
				res2the01*2/the01(1)+ &
				(LOG(the01(2))**2)*res20101numbis+&
				LOG(the01(2))*2*res2the0101+&
				res2the0101dsquare
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*(su12**vet12)+&
				gl01*vet01*(LOG(the01(2)*t3(i))**2)*(-1+&
				gl01*vet01)*(su01**vet01)*(su02**vet02)
				res1((nvamax+iter))=&
				(res1((nvamax+iter))*u1-v*v)/(u1*u1)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+&
				troncweib01011square
				!ici
				if(fix(2).eq.0)then
					
					u2=res2denum*the01(1)/the01(2)-&
					(res2denum-res201num)*the01(1)/the01(2)
					u2=u2*(su12**vet12)-&
					the01(1)*gl01*vet01*(su01**vet01)*(su02**vet02)/the01(2)
					
					u3=the01(1)*LOG(the01(2))*res20101numbis/the01(2)+&
					the01(1)*res2the0101/the01(2)+&
					(the01(1)*LOG(the01(2))+2)*res2denum/the01(2)+&
					the01(1)*res2thenum/the01(2)-&
					(3*the01(1)*LOG(the01(2))+2)* &
					(res2denum-res201num)/the01(2)-3* &
					the01(1)*res2the01/the01(2)
					
					u3=u3*(su12**vet12)+&
					gl01*vet01*(su01**vet01)*(su02**vet02)*(-1-&
					LOG(the01(2)*t3(i))*the01(1)+&
					gl01*vet01*LOG(the01(2)*t3(i))*the01(1))/the01(2)
					
					iter = iter +1
					res1((nvamax+iter))=&
					(u3*u1-v*u2)/(u1*u1)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+&
					troncweib010112
				endif
				!ici
				if(fix(3).eq.0)then
				
					u3=-LOG(the01(2))*LOG(the02(2))*res202num-&
					LOG(the01(2))*res2the02-LOG(the02(2))*res2the02-&
					res2the0202square-LOG(the02(2))*res202num/the01(1)-&
					res2the02/the01(1)+&
					LOG(the01(2))*LOG(the02(2))*res20102num+&
					LOG(the01(2))*res2the0102+&
					LOG(the02(2))*res2the0102+&
					res2the0102square
					u3=u3*(su12**vet12)+&
					gl01*vet01*gl02*vet02*(su01**vet01)*(su02**vet02)* &
					LOG(the01(2)*t3(i))*LOG(the02(2)*t3(i))
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					u2=u2*(su12**vet12)-&
					LOG(the02(2)*t3(i))*gl02*vet02* &
					(su01**vet01)*(su02**vet02)
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!ici
				if(fix(4).eq.0)then
				
					u3=-LOG(the01(2))*the02(1)*res202num/the02(2)-&
					the02(1)*res2the02/the02(2)-&
					the02(1)*res202num/(the01(1)*the02(2))+&
					LOG(the01(2))*the02(1)*res20102num/the02(2)+&
					the02(1)*res2the0102/the02(2)
					
					u3=u3*(su12**vet12)+&
					gl01*vet01*gl02*vet02*(su01**vet01)* &
					(su02**vet02)*the02(1)* &
					LOG(the01(2)*t3(i))/the02(2)
					
					u2=-the02(1)*res202num/the02(2)
					u2=u2*(su12**vet12)-&
					the02(1)*gl02*vet02*(su01**vet01)* &
					(su02**vet02)/the02(2)
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!ici
				if(fix(5).eq.0)then
				
				
					u3=LOG(the01(2))*LOG(the12(2))*res212num+&
					LOG(the01(2))*res2the12+LOG(the12(2))*res2the12+&
					res2the1212square+LOG(the12(2))*res212num/the01(1)+&
					res2the12/the01(1)-&
					LOG(the01(2))*LOG(the12(2))*res20112num-&
					LOG(the01(2))*res2the0112-&
					LOG(the12(2))*res2the0112-&
					res2the0112square

					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)* &
					LOG(the12(2)*t3(i))* &
					((LOG(the01(2))+(1/the01(1)))*res2denum +&
					res2thenum-LOG(the01(2))*(res2denum-res201num)-&
					res2the01)
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)-&
					LOG(the12(2)*t3(i))*gl12*vet12* &
					res2denum*(su12**vet12)
						
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					v*u2)/(u1*u1)
					
				endif
				!ici
				if(fix(6).eq.0)then
				
					
					u3=&
					LOG(the01(2))*the12(1)*res212num/the12(2)+&
					the12(1)/the12(2)*res2the12+res212num* &
					the12(1)/(the12(2)*the01(1))-&
					LOG(the01(2))*the12(1)*res20112num/the12(2)-&
					the12(1)/the12(2)*res2the0112
					
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*the12(1)* &
					((LOG(the01(2))+(1/the01(1)))*res2denum +&
					res2thenum-LOG(the01(2))*(res2denum-res201num)-&
					res2the01)/the12(2)
					
					
					u2=the12(1)*res212num/the12(2)
					
					u2=u2*(su12**vet12)-&
					the12(1)*gl12*vet12*(su12**vet12)* &
					res2denum/the12(2)
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				
				!ici
				if(nva01nofix.gt.0) then 
					
					u3=res2the0101+ &
					LOG(the01(2))*res20101numbis-&
					LOG(the01(2))*3*(res2denum-res201num)-&
					res2the01*3-(res2denum-res201num)/the01(1)+&
					(LOG(the01(2))+(1/the01(1)))*res2denum+&
					res2thenum
					
					u3= u3*(su12**vet12)-&
					LOG(the01(2)*t3(i))*gl01*vet01*(1-&
					gl01*vet01)*(su01**vet01)*(su02**vet02)
					
					u2=res201num*(su12**vet12)-&
					gl01*vet01*(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01011beta01
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					
					u3=-(LOG(the01(2))+&
					(1/the01(1)))*res202num-&
					res2the02+&
					LOG(the01(2))*res20102num+&
					res2the0102
					
					u3=u3*(su12**vet12)+&
					LOG(the01(2)*t3(i))*gl01*vet01* &
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
					u2=-res202num*(su12**vet12)-&
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					iter = iter +nva02nofix
				endif
				!ici
				if(nva12nofix.gt.0) then 
					u3=(LOG(the01(2))+&
					(1/the01(1)))*res212num+&
					res2the12-&
					LOG(the01(2))*res20112num-&
					res2the0112
					
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)* &
					((LOG(the01(2))+(1/the01(1)))*res2denum +&
					res2thenum-LOG(the01(2))*(res2denum-res201num)-&
					res2the01)
					
					u2=res212num*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*res2denum
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
				endif
			endif
			
			!write(6, *) " done x101 -c6" 	
			if(fix(2).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=the01(1)*res2denum/the01(2) -&
				the01(1)*(res2denum-res201num)/the01(2)
				
				v= v*(su12**vet12)-&
				the01(1)*gl01*vet01*(su01**vet01)* &
				(su02**vet02)/the01(2)
				
				res1(nweib)=v/u1
				res1(nweib)=res1(nweib)+troncweib01012
				
				
				u2=((the01(1)/the01(2))**2)*res20101numbis+&
					the01(1)*(the01(1)-1)*res2denum/((the01(2))**2)-&
					(res2denum-res201num)*(3*(the01(1)**2)-&
					the01(1))/(the01(2)**2)
					
				u2=u2*(su12**vet12)+&
				gl01*vet01*(su01**vet01)*(su02**vet02)* &
				the01(1)*(1-the01(1)+gl01*vet01*the01(1))/(the01(2)**2)
				
				res1((nvamax+iter))=&
				(u2*u1-v*v)/(u1*u1)
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib01012square
				
				if(fix(3).eq.0)then
				!ici
					u3=-LOG(the02(2))*the01(1)*res202num/the01(2)-&
					the01(1)*res2the02/the01(2)+&
					the01(1)*LOG(the02(2))*res20102num/the01(2)+&
					the01(1)*res2the0102/the01(2)
					
					u3=u3*(su12**vet12)+&
					gl01*vet01*gl02*vet02*the01(1)* &
					LOG(the02(2)*t3(i))*(su01**vet01)* &
					(su02**vet02)/the01(2)
					
					
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					u2=u2*(su12**vet12)-&
					LOG(the02(2)*t3(i))*gl02*vet02* &
					(su01**vet01)*(su02**vet02)
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!ici!
				if(fix(4).eq.0)then
				
					u3=the01(1)*the02(1)*res20102num-&
					the01(1)*the02(1)*res202num
					u3=u3/(the01(2)*the02(2))
					u3=u3*(su12**vet12)+&
					gl01*vet01*gl02*vet02* &
					(su01**vet01)*(su02**vet02)* &
					the01(1)*the02(1)/(the01(2)*the02(2))
					
					u2=-the02(1)*res202num/the02(2)
					u2=u2*(su12**vet12)-&
					the02(1)*gl02*vet02*(su01**vet01)* &
					(su02**vet02)/the02(2)
				
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!ici
				if(fix(5).eq.0)then
				
					u3=LOG(the12(2))*the01(1)*res212num/the01(2)+&
					the01(1)*res2the12/the01(2)-&
					the01(1)*LOG(the12(2))*res20112num/the01(2)-&
					the01(1)*res2the0112/the01(2)
					
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*LOG(the12(2)* &
					t3(i))*(the01(1)*res2denum/the01(2) -&
					the01(1)*(res2denum-res201num)/the01(2))
				
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)-&
					LOG(the12(2)*t3(i))*gl12*vet12* &
					(su12**vet12)*res2denum
						
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					v*u2)/(u1*u1)
					
				endif
				!ici
				if(fix(6).eq.0)then
					
					u3=-the01(1)*the12(1)*res20112num+&
					the01(1)*the12(1)*res212num
					u3=u3/(the01(2)*the12(2))
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*the12(1)* &
					(the01(1)*res2denum/the01(2) -&
					the01(1)*(res2denum-res201num)/the01(2))/ &
					the12(2)
					
					
					u2=the12(1)*res212num/the12(2)
					u2=u2*(su12**vet12)-&
					the12(1)*gl12*vet12*(su12**vet12)* &
					res2denum/the12(2)
				
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!ici
				if(nva01nofix.gt.0) then 
					
					u3=the01(1)*res20101numbis/the01(2)+&
					the01(1)*res2denum/the01(2)-&
					(3*the01(1)*(res2denum-res201num)/the01(2))
					u3=u3*(su12**vet12)+&
					the01(1)*gl01*vet01*(gl01*vet01-1)* &
					(su01**vet01)*(su02**vet02)/the01(2)
					
					u2=res201num*(su12**vet12)-&
					gl01*vet01*(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01012beta01
					
					iter = iter +nva01nofix
				endif
				!ici
				if(nva02nofix.gt.0) then 
					
					u3=res20102num-res202num
					u3=u3*the01(1)/the01(2)
					u3=u3*(su12**vet12)+&
					the01(1)*gl01*vet01*gl02*vet02* &
					(su01**vet01)*(su02**vet02)/the01(2)
					
					u2=-res202num*(su12**vet12)-&
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!ici
					u3=-res20112num+res212num
					u3=u3*the01(1)/the01(2)
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)* &
					(the01(1)*res2denum/the01(2) -&
					the01(1)*(res2denum-res201num)/the01(2))
				
					
					u2=res212num*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*res2denum
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					iter = iter +nva12nofix
				endif
			endif
			
					!write(6, *) " done x201 -c6" 	
				!ici
				
			if(fix(3).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=-LOG(the02(2))*res202num-&
				res2the02
				
				v=v*(su12**vet12)-&
				LOG(the02(2)*t3(i))*gl02*vet02* &
				(su01**vet01)*(su02**vet02)
				
				res1(nweib)=v/u1
				res1(nweib)=res1(nweib)+troncweib02021
				
				res1((nvamax+iter))=&
				((LOG(the02(2)))**2)*res20202num+&
				(2*LOG(the02(2))*res2the0202)+&
				res2the0202dsquare-&
				((LOG(the02(2)))**2)*res202num-&
				(2*LOG(the02(2))*res2the02)-&
				res2the0202square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*(su12**vet12)+&
				gl02*vet02*(su01**vet01)*(su02**vet02)* &
				(LOG(the02(2)*t3(i))**2)*(gl02*vet02-1)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*u1 -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				res1((nvamax+iter))=&
				res1((nvamax+iter))+&
				troncweib02021square
				
				!ici
				if(fix(4).eq.0)then
				
					u3=-(the02(1)*LOG(the02(2))+1)* &
					res202num/the02(2)-&
					the02(1)*res2the02/the02(2)+ &
					LOG(the02(2))*the02(1)* &
					res20202num/the02(2)+ &
					the02(1)/the02(2)*res2the0202
					
					u3=u3*(su12**vet12)+&
					gl02*vet02*(su01**vet01)*(su02**vet02)* &
					(-1-LOG(the02(2)*t3(i))*the02(1)+&
					gl02*vet02*LOG(the02(2)*t3(i))* &
					the02(1))/the02(2)
					
					u2=-the02(1)*res202num/the02(2)
					u2=u2*(su12**vet12)-&
					the02(1)*gl02*vet02*(su01**vet01)* &
					(su02**vet02)/the02(2)
					
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+&
					troncweib020212
				endif
				!ici
				if(fix(5).eq.0)then
				
					u3=-LOG(the02(2))*LOG(the12(2))*res20212num-&
					LOG(the02(2))*res2the0212-LOG(the12(2))*res2the0212-&
					res2the0212square
					
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*LOG(the12(2)* &
					t3(i))*(-LOG(the02(2))*res202num-&
					res2the02)
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)-&
					LOG(the12(2)*t3(i))*gl12*vet12* &
					(su12**vet12)*res2denum
						
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
				endif
				!ici
				if(fix(6).eq.0)then
				
					u3=-LOG(the02(2))*the12(1)*res20212num/the12(2)-&
					the12(1)*res2the0212/the12(2)
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*the12(1)* &
					(-1*LOG(the02(2))*res202num-&
					res2the02)/the12(2)
				
					
					
					u2=the12(1)*res212num/the12(2)
					u2=u2*(su12**vet12)-the12(1)*gl12*vet12* &
					(su12**vet12)*res2denum/the12(2)
					
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
				endif
				!ici
				if(nva01nofix.gt.0) then 
					
					u3=LOG(the02(2))*res20102num+&
					res2the0102-res202num*LOG(the02(2))-&
					res2the02
					
					u3=u3*(su12**vet12)+&
					LOG(the02(2)*t3(i))*gl01*vet01* &
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
					u2=res201num*(su12**vet12)-&
					gl01*vet01*(su01**vet01)* &
					(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
				endif
				!ici
				if(nva02nofix.gt.0) then 
					
					u3=-LOG(the02(2))*res202num-&
					res2the02+LOG(the02(2))*res20202num+&
					res2the0202
					u3=u3*(su12**vet12)+&
					LOG(the02(2)*t3(i))*gl02*vet02* &
					(su01**vet01)*(su02**vet02)* &
					(-1+gl02*vet02)
					
					u2=-res202num*(su12**vet12)-&
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02021beta02
					iter = iter +nva02nofix
					
				endif
				!ici
				if(nva12nofix.gt.0) then 
				
					u3=-LOG(the02(2))*res20212num-&
					res2the0212
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)* &
					(-1*LOG(the02(2))*res202num-&
					res2the02)
					
					u2=res212num*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*res2denum
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
				endif
			endif
			
             !write(6, *) " done x102 -c6" 	
			!ici
			if(fix(4).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=-the02(1)*res202num/the02(2)
				
				v=v*(su12**vet12)-&
				gl02*vet02*the02(1)*(su01**vet01)* &
				(su02**vet02)/the02(2)
				
				res1(nweib)=v/u1
				res1(nweib)=res1(nweib)+troncweib02022
				
				u3=((the02(1)/the02(2))**2)*res20202num-&
					(the02(1)*(the02(1)-1)/((the02(2))**2))*res202num

				
				u3=u3*(su12**vet12)+&
				gl02*vet02*(su01**vet01)*(su02**vet02)* &
				the02(1)*(1-the02(1)+the02(1)*gl02*vet02)/ &
				(the02(2)**2)
				
				res1((nvamax+iter))=&
				u1*u3-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib02022square
				
				!ici
				if(fix(5).eq.0)then
				
					u3=-the02(1)*LOG(the12(2))*res20212num/the02(2)-&
					the02(1)*res2the0212/the02(2)
					u3=u3*(su12**vet12)+&
					gl12*vet12*(su12**vet12)*LOG(the12(2)*t3(i))* &
					(the02(1)*res202num/the02(2))
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)-&
					LOG(the12(2)*t3(i))*gl12*vet12* &
					(su12**vet12)*res2denum
						
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
				endif
				!ici
				if(fix(6).eq.0)then
					
					u3=-the02(1)*the12(1)*res20212num
					u3=u3/(the02(2)*the12(2))
					u3=u3*(su12**vet12)+&
					gl12*vet12*(su12**vet12)*the12(1)* &
					(the02(1)*res202num/the02(2))/the12(2)
					
					u2=the12(1)*res212num/the12(2)
					u2=u2*(su12**vet12)-&
					the12(1)*gl12*vet12*(su12**vet12)* &
					res2denum/the12(2)
				
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
				endif
				
				!ici
				if(nva01nofix.gt.0) then 
					
					u3=-the02(1)*res202num/the02(2)+&
					the02(1)*res20102num/the02(2)		
					u3=u3*(su12**vet12)+&
					the02(1)*gl02*vet02*gl01*vet01* &
					(su01**vet01)*(su02**vet02)/the02(2)
					
					u2=res201num*(su12**vet12)-&
					gl01*vet01*(su01**vet01)*(su02**vet02)
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
				endif
				!ici
				if(nva02nofix.gt.0) then 
					
					u3=the02(1)*res20202num/the02(2)-&
					the02(1)*res202num/the02(2)
					u3=u3*(su12**vet12)+&
					the02(1)*gl02*vet02*(su01**vet01)* &
					(su02**vet02)*(gl02*vet02-1)/the02(2)
					
					u2=-res202num*(su12**vet12)-&
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02022beta02
					
					iter = iter +nva02nofix
				endif
				!ici
				if(nva12nofix.gt.0) then 
					u3=-res20212num*the02(1)/the02(2)
					u3=u3*(su12**vet12)+&
					gl12*vet12*(su12**vet12)* &
					(the02(1)*res202num/the02(2))
					
					u2=res212num*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*res2denum
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					iter = iter +nva12nofix
				endif
			endif	

			!ici
			!write(6, *) " done x202 -c6" 	
			if(fix(5).eq.0)then
			
			
				iter = iter + 1
				nweib = nweib + 1
				
				v=LOG(the12(2))*res212num+&
				res2the12
				v=v*(su12**vet12)-&
				LOG(the12(2)*t3(i))*gl12*vet12* &
				(su12**vet12)*res2denum
				
				res1(nweib)=v/u1
				
				res1((nvamax+iter))=&
				((LOG(the12(2)))**2)* &
				res21212num+2*LOG(the12(2))*res2the1212+&
				res2the1212dsquare+&
				((LOG(the12(2)))**2)* &
				res212num+2*LOG(the12(2))*res2the12+&
				res2the1212square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*(su12**vet12)+&
				gl12*vet12*(su12**vet12)*LOG(the12(2)* &
				t3(i))*(-2*(LOG(the12(2))*res212num+&
				res2the12)-LOG(the12(2)*t3(i))*res2denum* &
				(1-gl12*vet12))
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*u1 -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				
				!ici
				if(fix(6).eq.0)then
			
					u2=the12(1)*res212num/the12(2)
					
					u3=LOG(the12(2))*the12(1)*res21212num/the12(2)+&
					the12(1)*res2the1212/the12(2)+&
					(the12(1)*LOG(the12(2))+1)*res212num/the12(2)+&
					the12(1)*res2the12/the12(2)
					u3=u3*(su12**vet12)+&
					gl12*vet12*(su12**vet12)* &
					(-1*LOG(the12(2)*t3(i))*u2 - &
					(LOG(the12(2))*res212num+&
					res2the12)*the12(1)/the12(2)+&
					res2denum*(-1-LOG(the12(2)*t3(i))*the12(1)+&
					gl12*vet12*the12(1)*LOG(the12(2)*t3(i)))/ &
					the12(2))
					
					u2=u2*(su12**vet12)-&
					the12(1)*gl12*vet12*(su12**vet12)* &
					res2denum/the12(2)
					
					
					iter = iter +1
					
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
				endif
				
				if(nva01nofix.gt.0) then 
					!ici
					u3=-LOG(the12(2))*res20112num-&
					res2the0112+res212num*LOG(the12(2))+&
					res2the12
					u3=u3*(su12**vet12)-&
					LOG(the12(2)*t3(i))*gl12*vet12* &
					(su12**vet12)*res201num
					
					
					u2=res201num*(su12**vet12)-&
					gl01*vet01*(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
				endif
				if(nva02nofix.gt.0) then 
					!ici
					u3=-LOG(the12(2))*res20212num-&
					res2the0212
					u3=u3*(su12**vet12)-&
					LOG(the12(2)*t3(i))*gl12*vet12* &
					(su12**vet12)*(-res202num)
					
					u2=-res202num*(su12**vet12)-&
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					iter = iter +nva02nofix
					
				endif
				if(nva12nofix.gt.0) then 
					!ici
					u3=LOG(the12(2))*res21212num+&
					res2the1212+LOG(the12(2))*res212num+&
					res2the12
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)* &
					(LOG(the12(2))*res212num+&
					res2the12)-gl12*vet12*(su12**vet12)* &
					LOG(the12(2)*t3(i))* &
					(res212num + res2denum*(1-gl12*vet12))
					
					u2=res212num*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*res2denum
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
				endif
			endif
			
			  !ici
			!write(6, *) " done x112 -c6" 	
			if(fix(6).eq.0)then
			
			
				iter = iter + 1
				nweib = nweib + 1
				
				v=the12(1)*res212num/the12(2)
				v=v*(su12**vet12)-&
				the12(1)*gl12*vet12*(su12**vet12)* &
				res2denum/the12(2)
				
				res1(nweib)=v/u1
				
				u3=((the12(1)/the12(2))**2)*res21212num+&
					(the12(1)*(the12(1)-1)/((the12(2))**2))*res212num
				u3=u3*(su12**vet12)+&
				gl12*vet12*(su12**vet12)*the12(1)* &
				(-2*(the12(1)*res212num/the12(2))-&
				res2denum*(-1+the12(1)-&
				gl12*vet12*the12(1))/the12(2))/the12(2)
				
				
				res1((nvamax+iter))=&
				u1*u3-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				
				
				if(nva01nofix.gt.0) then 
					!ici
					u3=the12(1)*res212num/the12(2)-&
					the12(1)*res20112num/the12(2)					
					u3=u3*(su12**vet12)-&
					the12(1)*gl12*vet12*res201num* &
					(su12**vet12)/the12(2)
					
					u2=res201num*(su12**vet12)-&
					gl01*vet01*(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
				endif
				
				if(nva02nofix.gt.0) then 
					!ici
					u3=-the12(1)*res20212num/the12(2)
					u3=u3*(su12**vet12)+&
					the12(1)*gl12*vet12*(su12**vet12)* &
					(res202num)/the12(2)
					
					u2=-res202num*(su12**vet12)-&
					gl02*vet02*(su01**vet01)*(su02**vet02)
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!ici
					u3=res21212num*the12(1)/the12(2)+&
					the12(1)*res212num/the12(2)
					u3=u3*(su12**vet12)-&
					gl12*vet12*(su12**vet12)* &
					(the12(1)*res212num/the12(2))-&
					gl12*vet12*(su12**vet12)*the12(1)* &
					(res212num+res2denum*(1-gl12*vet12))/ &
					the12(2)
					
					u2=res212num*(su12**vet12)-&
					gl12*vet12*(su12**vet12)*res2denum
							
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
				endif
			endif	
			
			!write(6, *) " done x212 -c6" 	
				v=(su12**vet12)*res2denum+&
				(su01**vet01)*(su02**vet02)

				if(nva01nofix.gt.0) then 

				u1=(-gl01*vet01)*(su01**vet01)*(su02**vet02)+&
				(su12**vet12)*res201num
			        res1((nvaweib+1):(nvaweib+nva01nofix))=&
				ve01nofix(i,:)*u1/v
				res1((nvaweib+1):(nvaweib+nva01nofix))=&
				res1((nvaweib+1):(nvaweib+nva01nofix))+tronc01

				res1((nvamax12weib12+1):nvamax01)=&
				res20101num*(su12**vet12)
				resint=(su01**vet01)*(su02**vet02)
				resint=resint*gl01*vet01*(1-gl01*vet01)
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)-resint
				
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)*ve01square(i,:)
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)/v

				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)+tronc01square

				end if 
				
				if(nva02nofix.gt.0) then 

				u2=-gl02*vet02*(su01**vet01)*(su02**vet02)
				u2=u2-(su12**vet12)*res202num
				res1((nva01nofix+1+nvaweib):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1+nvaweib):nva0102)=&
				res1((nva01nofix+1+nvaweib):nva0102)+&
				tronc02

				res1((nvamax0112+1):nvamax02)=&
				(su12**vet12)*(res20202num-res202num)-&
				gl02*vet02*(su01**vet01)*(su02**vet02)*(1-gl02*vet02)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)*ve02square(i,:)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)/v

				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)+tronc02square

				
				end if 

				if(nva12nofix.gt.0) then 
				
				u3=-gl12*vet12*(su12**vet12)*res2denum+&
				(su12**vet12)*res212num
				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v

				res1((nvamax0212+1):nvamax12)=&
				res21212num+(1-gl12*vet12)*res212num-&
				gl12*vet12*(1-gl12*vet12)*res2denum-&
				gl12*vet12*res212num
				resint=(su12**vet12)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)*resint*ve12square(i,:)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)/v
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)


				end if 
                        	
				
				if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 

				kfix=nvamax01+1
				lfix=kfix-1+nva02nofix

				do j=1,nva01nofix
			   		res1(kfix:lfix)=&
					-(res202num-res20102num)*(su12**vet12)
					resint=(su01**vet01)*(su02**vet02)
					res1(kfix:lfix)=res1(kfix:lfix)+&
					gl01*vet01*gl02*vet02*resint
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)*v
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva02nofix
				end do

				end if 

				if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 

				
				kfix=nvamax0102+1
				lfix=kfix-1+nva12nofix

				do j=1,nva01nofix
			  		res1(kfix:lfix)=-gl12*vet12*res201num+&
					res212num-res20112num
					resint=(su12**vet12)
			   		res1(kfix:lfix)=&
			  		res1(kfix:lfix)*resint*ve01nofix(i,j)*ve12nofix(i,:)
		           		res1(kfix:lfix)=&
			   		res1(kfix:lfix)*v
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
			   		res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva12nofix
				end do

				end if 

				
				if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 

				kfix=nvamax02+1
				lfix=kfix-1+nva12nofix

				do j=1,nva02nofix
			  		res1(kfix:lfix)=&
					-res20212num+gl12*vet12*res202num
					resint=(su12**vet12)
					res1(kfix:lfix)=&
					res1(kfix:lfix)*resint*ve02nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
					res1(kfix:lfix)*v
					res1(kfix:lfix)=&
					res1(kfix:lfix)-u2*u3*ve02nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva12nofix
				end do

				end if 


				
                            else ! passage 0-->2  
				
				!write(6, *) "start c=7" 
				
			call qgaussweibbetaderiv(t1(i),t3(i),the01,&
			the02,the12,res2denum,res2the01,res2the02,&
			res2the12,res2thenum,res2thenumsquare,&
			res2the0101,res2the0102,res2the0112,res2the0202,&
			res2the0212,res2the1212,res2the0101square,&
			res2the0102square,res2the0112square,&
			res2the0202square,res2the0212square,&
			res2the1212square,res2the0101dsquare,&
			res2the0202dsquare,res2the1212dsquare,&
			res201num,res202num,res212num,res20101num,&
			res20101numbis,&
			res20102num,res20112num,&
			res20202num,res20212num,res21212num,&
			vet01,vet02,vet12)
			
			call fonct(t3(i),the01,ri01,gl01,su01)
			call fonct(t3(i),the02,ri02,gl02,su02)
			call fonct(t3(i),the12,ri12,gl12,su12)

			!write(6, *) "calculate integrals -c7" 
			
			iter = 0
			nweib = 0
			
			
			u1=(su12**vet12)*res2denum*ri12*vet12+&
			(su01**vet01)*(su02**vet02)*ri02*vet02
				
			
			if(fix(1).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=(LOG(the01(2))+(1/the01(1)))*res2denum +&
				res2thenum-LOG(the01(2))*(res2denum-res201num)-&
				res2the01
				!stop
				
				v=v*(su12**vet12)*ri12*vet12-&
				LOG(the01(2)*t3(i))*gl01*vet01*(su01**vet01)* &
				ri02*vet02*(su02**vet02)
				
				res1(nweib)=v/u1
				
				res1(nweib)=res1(nweib)+troncweib01011
				!ici

				
				res1((nvamax+iter))=&
				LOG(the01(2))*2/the01(1)*res2denum+&
				res2thenum*2/the01(1)+&
				(LOG(the01(2))**2)*res2denum+&
				LOG(the01(2))*2*res2thenum+&
				res2thenumsquare-3*(LOG(the01(2))**2)* &
				(res2denum-res201num)-&
				LOG(the01(2))*6*res2the01-3*res2the0101square-&
				LOG(the01(2))*2/the01(1)*(res2denum-res201num)- &
				res2the01*2/the01(1)+ &
				(LOG(the01(2))**2)*res20101numbis+&
				LOG(the01(2))*2*res2the0101+&
				res2the0101dsquare
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*(su12**vet12)*ri12*vet12+&
				gl01*vet01*(LOG(the01(2)*t3(i))**2)*(-1+&
				gl01*vet01)*(su01**vet01)*(su02**vet02)*ri02*vet02
				
				res1((nvamax+iter))=&
				(res1((nvamax+iter))*u1-v*v)/(u1*u1)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+&
				troncweib01011square
				!la
				if(fix(2).eq.0)then
					
					u2=res2denum*the01(1)/the01(2)-&
					(res2denum-res201num)*the01(1)/the01(2)
					u2=u2*(su12**vet12)*ri12*vet12-&
					ri02*vet02*the01(1)*gl01*vet01* &
					(su01**vet01)*(su02**vet02)/the01(2)
					
					u3=the01(1)*LOG(the01(2))*res20101numbis/the01(2)+&
					the01(1)*res2the0101/the01(2)+&
					(the01(1)*LOG(the01(2))+2)*res2denum/the01(2)+&
					the01(1)*res2thenum/the01(2)-&
					(3*the01(1)*LOG(the01(2))+2)* &
					(res2denum-res201num)/the01(2)-3* &
					the01(1)*res2the01/the01(2)
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri02*vet02*gl01*vet01*(su01**vet01)* &
					(su02**vet02)*(-1-&
					LOG(the01(2)*t3(i))*the01(1)+&
					gl01*vet01*LOG(the01(2)*t3(i))*the01(1))/the01(2)
					
					iter = iter +1
					res1((nvamax+iter))=&
					(u3*u1-v*u2)/(u1*u1)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+troncweib010112
				endif
				!la
				if(fix(3).eq.0)then
				
					u3=&
					-LOG(the01(2))*LOG(the02(2))*res202num-&
					LOG(the01(2))*res2the02-&
					LOG(the02(2))*res2the02-&
					res2the0202square-&
					LOG(the02(2))*res202num/the01(1)-&
					res2the02/the01(1)+&
					LOG(the01(2))*LOG(the02(2))*res20102num+&
					LOG(the01(2))*res2the0102+&
					LOG(the02(2))*res2the0102+&
					res2the0102square
					
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri02*vet02*gl01*vet01* &
					(su01**vet01)*(su02**vet02)* &
					LOG(the01(2)*t3(i))* &
					(LOG(the02(2)*t3(i))*gl02*vet02-&
					LOG(the02(2)*t3(i))-&
					(1/the02(1)))
					
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri02*vet02*(su01**vet01)*(su02**vet02)* &
					(LOG(the02(2)*t3(i))+(1/the02(1))-&
					LOG(the02(2)*t3(i))*gl02*vet02)
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!la
				if(fix(4).eq.0)then
				
					u3=-LOG(the01(2))*the02(1)*res202num/the02(2)-&
					the02(1)*res2the02/the02(2)-&
					the02(1)*res202num/(the01(1)*the02(2))+&
					LOG(the01(2))*the02(1)*res20102num/the02(2)+&
					the02(1)*res2the0102/the02(2)
					
					
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01* &
					(1-gl02*vet02)*(su01**vet01)* &
					(su02**vet02)*the02(1)* &
					LOG(the01(2)*t3(i))/the02(2)
					
					u2=-the02(1)*res202num/the02(2)
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri02*vet02*the02(1)*(1-gl02*vet02)*(su01**vet01)* &
					(su02**vet02)/the02(2)

					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!la
				if(fix(5).eq.0)then
				
					
					u3=LOG(the01(2))*LOG(the12(2))*res212num+&
					LOG(the01(2))*res2the12+LOG(the12(2))*res2the12+&
					res2the1212square+LOG(the12(2))*res212num/the01(1)+&
					res2the12/the01(1)-&
					LOG(the01(2))*LOG(the12(2))*res20112num-&
					LOG(the01(2))*res2the0112-&
					LOG(the12(2))*res2the0112-&
					res2the0112square
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)* &
					(LOG(the12(2)*t3(i))+ &
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)* &
					((LOG(the01(2))+(1/the01(1)))*res2denum +&
				res2thenum-LOG(the01(2))*(res2denum-res201num)-&
				res2the01)
					
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*res2denum* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)
						
						
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					v*u2)/(u1*u1)
					
				endif
				!la
				if(fix(6).eq.0)then
				
					u3=&
					LOG(the01(2))*the12(1)*res212num/the12(2)+&
					the12(1)/the12(2)*res2the12+res212num* &
					the12(1)/(the12(2)*the01(1))-&
					LOG(the01(2))*the12(1)*res20112num/the12(2)-&
					the12(1)/the12(2)*res2the0112
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*the12(1)* &
					(1-gl12*vet12)* &
					((LOG(the01(2))+(1/the01(1)))*res2denum +&
					res2thenum-LOG(the01(2))*(res2denum-res201num)-&
					res2the01)/the12(2)
					
					
					
					u2=the12(1)*res212num/the12(2)
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*the12(1)*(1-gl12*vet12)* &
					(su12**vet12)* &
					res2denum/the12(2)
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				
				!la
				if(nva01nofix.gt.0) then 
					
					u3=res2the0101+ &
					LOG(the01(2))*res20101numbis-&
					LOG(the01(2))*3*(res2denum-res201num)-&
					res2the01*3-(res2denum-res201num)/the01(1)+&
					(LOG(the01(2))+(1/the01(1)))*res2denum+&
					res2thenum
					
					u3= u3*(su12**vet12)*ri12*vet12-&
					ri02*vet02*LOG(the01(2)*t3(i))*gl01*vet01*(1-&
					gl01*vet01)*(su01**vet01)*(su02**vet02)
					
					u2=res201num*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01*(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01011beta01
					
					iter = iter +nva01nofix
				endif
				!la
				if(nva02nofix.gt.0) then 
					
					
					u3=-(LOG(the01(2))+&
					(1/the01(1)))*res202num-&
					res2the02+&
					LOG(the01(2))*res20102num+&
					res2the0102
					
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri02*vet02*LOG(the01(2)*t3(i))*gl01*vet01* &
					(su01**vet01)*(su02**vet02)* &
					(1-gl02*vet02)
					
					u2=-res202num*(su12**vet12)*ri12*vet12+&
					(1-gl02*vet02)*(su01**vet01)*(su02**vet02)* &
					ri02*vet02
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					iter = iter +nva02nofix
				endif
				!ici
				if(nva12nofix.gt.0) then 
				!la
					u3=(LOG(the01(2))+&
					(1/the01(1)))*res212num+&
					res2the12-&
					LOG(the01(2))*res20112num-&
					res2the0112
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1-gl12*vet12)*(su12**vet12)* &
					((LOG(the01(2))+(1/the01(1)))*res2denum +&
					res2thenum-LOG(the01(2))*(res2denum-res201num)-&
					res2the01)
					
					u2=res212num*(su12**vet12)*ri12*vet12+&
					(1-gl12*vet12)*(su12**vet12)*res2denum* &
					ri12*vet12
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
				endif
			endif
			
			!write(6, *) " done x101 -c7" 
			if(fix(2).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=the01(1)*res2denum/the01(2) -&
				the01(1)*(res2denum-res201num)/the01(2)
				
				v= v*(su12**vet12)*ri12*vet12-&
				ri02*vet02*the01(1)*gl01*vet01*(su01**vet01)* &
				(su02**vet02)/the01(2)
				
				res1(nweib)=v/u1
				res1(nweib)=res1(nweib)+troncweib01012
				!la
				
				u2=((the01(1)/the01(2))**2)*res20101numbis+&
					the01(1)*(the01(1)-1)*res2denum/((the01(2))**2)-&
					(res2denum-res201num)*(3*(the01(1)**2)-&
					the01(1))/(the01(2)**2)
					
				u2=u2*(su12**vet12)*ri12*vet12-&
				ri02*vet02*(1-gl01*vet01)*(su01**vet01)*(su02**vet02)* &
				((the01(1)/the01(2))**2)*gl01*vet01+&
				ri02*vet02*gl01*vet01*(su01**vet01)* &
				(su02**vet02)*the01(1)/(the01(2)**2)
				
				res1((nvamax+iter))=&
				(u2*u1-v*v)/(u1*u1)
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib01012square
				
				if(fix(3).eq.0)then
				!la!
					
					u3=-LOG(the02(2))*the01(1)* &
					res202num/the01(2)-&
					the01(1)*res2the02/the01(2)+&
					the01(1)*LOG(the02(2))*res20102num/the01(2)+&
					the01(1)*res2the0102/the01(2)
					
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01*the01(1)* &
					(su01**vet01)*(su02**vet02)* &
					(LOG(the02(2)*t3(i))*(1-&
					gl02*vet02)+(1/the02(1)))/the01(2)
					
					u2=-LOG(the02(2))*(res202num)-&
						res2the02
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri02*vet02*(su01**vet01)*(su02**vet02)* &
					(LOG(the02(2)*t3(i))+&
					(1/the02(1))-&
					gl02*vet02*LOG(the02(2)*t3(i)))
					
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!la
				if(fix(4).eq.0)then
				
					u3=the01(1)*the02(1)*res20102num-&
					the01(1)*the02(1)*res202num
					u3=u3/(the01(2)*the02(2))
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri02*vet02*gl01*vet01*(gl02*vet02-1)* &
					(su01**vet01)*(su02**vet02)* &
					the01(1)*the02(1)/(the01(2)*the02(2))
					
					u2=-the02(1)*res202num/the02(2)
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri02*vet02*the02(1)*(1-gl02*vet02)* &
					(su01**vet01)* &
					(su02**vet02)/the02(2)
				
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!la
				if(fix(5).eq.0)then
				
					u3=LOG(the12(2))*the01(1)*res212num/the01(2)+&
					the01(1)*res2the12/the01(2)-&
					the01(1)*LOG(the12(2))*res20112num/the01(2)-&
					the01(1)*res2the0112/the01(2)
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*(LOG(the12(2)* &
					t3(i)) + &
					(1/the12(1))-LOG(the12(2)*t3(i))*gl12*vet12)* &
					(the01(1)*res2denum/the01(2) -&
					the01(1)*(res2denum-res201num)/the01(2))
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*res2denum* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)
						
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					v*u2)/(u1*u1)
					
				endif
				!la
				if(fix(6).eq.0)then
					
					u3=-the01(1)*the12(1)*res20112num+&
					the01(1)*the12(1)*res212num
					
					u3=u3/(the01(2)*the12(2))
					
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*the12(1)* &
					(the01(1)*res2denum/the01(2) -&
					the01(1)*(res2denum-res201num)/the01(2))* &
					(1-gl12*vet12)/the12(2)
					
					
					u2=the12(1)*res212num/the12(2)
					
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*the12(1)*(1-gl12*vet12)* &
					(su12**vet12)*res2denum/the12(2)
				
					iter = iter +1
					res1((nvamax+iter))=(u3*u1-&
					u2*v)/(u1*u1)
				endif
				!la
				if(nva01nofix.gt.0) then 
					
					u3=the01(1)*res20101numbis/the01(2)+&
					the01(1)*res2denum/the01(2)-&
					(3*the01(1)*(res2denum-res201num)/the01(2))
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri02*vet02*the01(1)*(1-gl01*vet01)* &
					gl01*vet01*(su01**vet01)*(su02**vet02)/the01(2)
					
					u2=res201num*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01*(su01**vet01)* &
					(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))+&
					troncweib01012beta01
					
					iter = iter +nva01nofix
				endif
				!la
				if(nva02nofix.gt.0) then 
					
					u3=res20102num-res202num
					u3=u3*the01(1)/the01(2)
					u3=u3*(su12**vet12)*ri12*vet12-&
					the01(1)*gl01*vet01*(1-gl02*vet02)* &
					(su01**vet01)*(su02**vet02)/the01(2)* &
					ri02*vet02
					
					u2=-res202num*(su12**vet12)*ri12*vet12+&
					(1-gl02*vet02)*(su01**vet01)*(su02**vet02)* &
					ri02*vet02
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!la
					u3=-res20112num+res212num
					u3=u3*the01(1)/the01(2)
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1-gl12*vet12)*(su12**vet12)* &
					(the01(1)*res2denum/the01(2) -&
					the01(1)*(res2denum-res201num)/the01(2))
				
					
					u2=res212num*(su12**vet12)*ri12*vet12+&
					(1-gl12*vet12)*ri12*vet12*(su12**vet12)*res2denum
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u3*u1-&
					v*u2)/(u1*u1)
					iter = iter +nva12nofix
				endif
			endif
			!write(6, *) " done x201 -c7" 
					
				!ici
		!la fin de journee : 181224
			if(fix(3).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=-LOG(the02(2))*res202num-&
				res2the02
				
				v=v*(su12**vet12)*ri12*vet12+&
				ri02*vet02*(su01**vet01)*(su02**vet02)* &
				(LOG(the02(2)*t3(i))+&
				(1/the02(1))-&
				gl02*vet02*LOG(the02(2)*t3(i)))
				
				res1(nweib)=v/u1
				res1(nweib)=res1(nweib)+troncweib02021
				res1((nvamax+iter))=&
				((LOG(the02(2)))**2)* &
				res20202num+2*LOG(the02(2))*res2the0202+&
				res2the0202dsquare-&
				((LOG(the02(2)))**2)* &
				res202num-2*LOG(the02(2))*res2the02-&
				res2the0202square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))* &
				(su12**vet12)*ri12*vet12+&
				ri02*vet02*(su01**vet01)*(su02**vet02)* &
				(((LOG(the02(2)*t3(i))+&
				(1/the02(1))-&
				gl02*vet02*LOG(the02(2)*t3(i)))**2)-&
				((1/(the02(1)**2))+&
				(LOG(the02(2)*t3(i))**2)*gl02*vet02))
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*u1 -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				res1((nvamax+iter))=&
				res1((nvamax+iter))+&
				troncweib02021square
				
				!la
				if(fix(4).eq.0)then
				
					u3=-(the02(1)*LOG(the02(2))+1)* &
					res202num/the02(2)-&
					the02(1)*res2the02/the02(2)+ &
					LOG(the02(2))*the02(1)* &
					res20202num/the02(2)+ &
					the02(1)/the02(2)*res2the0202
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri02*vet02*(1-gl02*vet02)*the02(1)* &
					(su01**vet01)*(su02**vet02)* &
					(LOG(the02(2)*t3(i))+&
					(1/the02(1))-&
					LOG(the02(2)*t3(i))*gl02*vet02)/the02(2)+ &
					ri02*vet02*(su01**vet01)*(su02**vet02)* &
					(1-gl02*vet02-LOG(the02(2)*t3(i))*the02(1)* &
					gl02*vet02)/the02(2)
					
					u2=-the02(1)*res202num/the02(2)
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri02*vet02*the02(1)*(1-gl02*vet02)* &
					(su01**vet01)* &
					(su02**vet02)/the02(2)
					
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
					res1((nvamax+iter))=&
					res1((nvamax+iter))+troncweib020212
				endif
				!la
				if(fix(5).eq.0)then
				
					u3=-LOG(the02(2))*LOG(the12(2))*res20212num-&
					LOG(the02(2))*res2the0212-LOG(the12(2))*res2the0212-&
					res2the0212square
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)* &
					(-1*LOG(the02(2))*res202num-&
					res2the02)
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*res2denum* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)
						
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
				endif
				!la
				if(fix(6).eq.0)then
				
					u3=-LOG(the02(2))*the12(1)*res20212num/the12(2)-&
					the12(1)*res2the0212/the12(2)
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1-gl12*vet12)*(su12**vet12)*the12(1)* &
					(-1*LOG(the02(2))*res202num-&
					res2the02)/the12(2)
				
					
					
					u2=the12(1)*res212num/the12(2)
					u2=u2*(su12**vet12)*ri12*vet12+ &
					ri12*vet12*the12(1)*(1-gl12*vet12)* &
					(su12**vet12)*res2denum/the12(2)
					
				
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
				endif
				!la
				if(nva01nofix.gt.0) then 
					
					
					u3=LOG(the02(2))*res20102num+&
					res2the0102-res202num*LOG(the02(2))-&
					res2the02
					
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01* &
					(su01**vet01)*(su02**vet02)* &
					(LOG(the02(2)*t3(i))+&
					(1/the02(1))-&
					LOG(the02(2)*t3(i))*gl02*vet02)
					
					u2=res201num*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01*(su01**vet01)* &
					(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
				endif
				!la
				if(nva02nofix.gt.0) then 
					
					u3=-LOG(the02(2))*res202num-&
					res2the02+LOG(the02(2))*res20202num+&
					res2the0202
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri02*vet02*(1-gl02*vet02)* &
					(su01**vet01)*(su02**vet02)* &
					(LOG(the02(2)*t3(i))+&
					(1/the02(1))-&
					LOG(the02(2)*t3(i))*gl02*vet02)-&
					ri02*vet02*gl02*vet02*(su01**vet01)* &
					(su02**vet02)*LOG(the02(2)*t3(i))
					
					u2=-res202num*(su12**vet12)*ri12*vet12+&
					(1-gl02*vet02)*(su01**vet01)* &
					(su02**vet02)*ri02*vet02
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02021beta02
					iter = iter +nva02nofix
					
				endif
				!la
				if(nva12nofix.gt.0) then 
				
					u3=-LOG(the02(2))*res20212num-&
					res2the0212
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1-gl12*vet12)*(su12**vet12)* &
					(-1*LOG(the02(2))*res202num-&
					res2the02)
					
					u2=res212num*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1-gl12*vet12)*(su12**vet12)*res2denum
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
				endif
			endif
			
             !write(6, *) " done x102 -c7" 
			!la
			if(fix(4).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				v=-the02(1)*res202num/the02(2)
				
				v=v*(su12**vet12)*ri12*vet12+&
				ri02*vet02*(1-gl02*vet02)* &
				the02(1)*(su01**vet01)* &
				(su02**vet02)/the02(2)
				
				res1(nweib)=v/u1
				res1(nweib)=res1(nweib)+troncweib02022
				!la
				
				u3=((the02(1)/the02(2))**2)*res20202num-&
					(the02(1)*(the02(1)-1)/((the02(2))**2))*res202num
				
				u3=u3*(su12**vet12)*ri12*vet12+&
				ri02*vet02*(1-gl02*vet02)* &
				(su01**vet01)*(su02**vet02)* &
				the02(1)*(-gl02*vet02*the02(1)/the02(2)- &
				(1/the02(2))+&
				the02(1)/the02(2))/the02(2)-&
				ri02*vet02*gl02*vet02*(the02(1)**2)* &
				(su01**vet01)*(su02**vet02)/(the02(2)**2)
				
				res1((nvamax+iter))=&
				u1*u3-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))+troncweib02022square
				
				!write(6, *) " done x202/x202 -c7" 
				if(fix(5).eq.0)then
				
					u3=-the02(1)*LOG(the12(2))*res20212num/the02(2)-&
					the02(1)*res2the0212/the02(2)
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri12*vet12*(su12**vet12)* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))- &
					LOG(the12(2)*t3(i))*gl12*vet12)* &
					(the02(1)*res202num/the02(2))
					
					u2=LOG(the12(2))*(res212num)+&
						res2the12
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*res2denum* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)
						
					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					v*u2
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
					!write(6, *) " done x202/x112 -c7" 
					
				endif
				!la
				if(fix(6).eq.0)then
					
					u3=-the02(1)*the12(1)*res20212num
					u3=u3/(the02(2)*the12(2))
					
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri12*vet12*(su12**vet12)*the12(1)* &
					(1-gl12*vet12)* &
					(the02(1)*res202num/the02(2))/the12(2)
					
				
					u2=the12(1)*res212num/the12(2)
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*the12(1)*(1-gl12*vet12)* &
					(su12**vet12)*res2denum/the12(2)

					iter = iter +1
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
					!write(6, *) " done x202/x212 -c7" 
				endif
				
				!la
				if(nva01nofix.gt.0) then 
					
					u3=-the02(1)*res202num/the02(2)+&
					the02(1)*res20102num/the02(2)		
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri02*vet02*the02(1)*(1- &
					gl02*vet02)*gl01*vet01* &
					(su01**vet01)*(su02**vet02)/the02(2)
					
					u2=res201num*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01* &
					(su01**vet01)*(su02**vet02)
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
					
					!write(6, *) " done x202/b01 -c7" 
				endif
				!la
				if(nva02nofix.gt.0) then 
					
					u3=the02(1)*res20202num/the02(2)-&
					the02(1)*res202num/the02(2)
					u3=u3*(su12**vet12)*ri12*vet12+&
					the02(1)*((1-gl02*vet02)**2)*(su01**vet01)* &
					(su02**vet02)*ri02*vet02/the02(2)-&
					ri02*vet02*the02(1)*gl02*vet02* &
					(su01**vet01)*(su02**vet02)/the02(2)
					
					u2=-res202num*(su12**vet12)*ri12*vet12+&
					(1-gl02*vet02)*(su01**vet01)* &
					(su02**vet02)*ri02*vet02
					
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))+&
					troncweib02022beta02
					
					iter = iter +nva02nofix
					!write(6, *) " done x202/b02 -c7" 
				endif
				!la
				if(nva12nofix.gt.0) then 
					u3=-res20212num*the02(1)/the02(2)
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri12*vet12*(1-gl12*vet12)*(su12**vet12)* &
					(the02(1)*res202num/the02(2))
					
					u2=res212num*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1- &
					gl12*vet12)*(su12**vet12)*res2denum
							
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					iter = iter +nva12nofix
					!write(6, *) " done x202/b12 -c7" 
				endif
			endif	
			!write(6, *) " done x202 -c7" 
			!la!
			if(fix(5).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=LOG(the12(2))*res212num+&
				res2the12
				v=v*(su12**vet12)*ri12*vet12+&
				ri12*vet12*(su12**vet12)*res2denum* &
				(LOG(the12(2)*t3(i))+&
				(1/the12(1))-&
				LOG(the12(2)*t3(i))*gl12*vet12)
				
				res1(nweib)=v/u1
				!la
				res1((nvamax+iter))=&
				((LOG(the12(2)))**2)* &
				res21212num+2*LOG(the12(2))*res2the1212+&
				res2the1212dsquare+&
				((LOG(the12(2)))**2)* &
				res212num+2*LOG(the12(2))*res2the12+&
				res2the1212square
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*(su12**vet12)*ri12*vet12+&
				ri12*vet12*(su12**vet12)* &
				(LOG(the12(2)*t3(i)) +&
				(1/the12(1))-&
				LOG(the12(2)*t3(i))*gl12*vet12)* &
				(2*(LOG(the12(2))*res212num+&
				res2the12)+res2denum*(LOG(the12(2)*t3(i))+&
				(1/the12(1))-&
				LOG(the12(2)*t3(i))*gl12*vet12))-&
				ri12*vet12*(su12**vet12)*res2denum* &
				((1/(the12(1)**2))+&
				(LOG(the12(2)*t3(i))**2)*gl12*vet12)
				
				res1((nvamax+iter))=&
				res1((nvamax+iter))*u1 -v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				
				!write(6, *) " done x112/x112 -c7" 
				if(fix(6).eq.0)then
	
					u2=the12(1)*res212num/the12(2)
					
					u3=LOG(the12(2))*the12(1)*res21212num/the12(2)+&
					the12(1)*res2the1212/the12(2)+&
					(the12(1)*LOG(the12(2))+1)*res212num/the12(2)+&
					the12(1)*res2the12/the12(2)
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1-gl12*vet12)* &
					(su12**vet12)*the12(1)* &
					((LOG(the12(2))*res212num+&
					res2the12)+&
					res2denum*(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12))/the12(2)
					
					u3=u3+ &
					ri12*vet12*(su12**vet12)*res2denum* &
					(1-LOG(the12(2)*t3(i))*the12(1)* &
					gl12*vet12-gl12*vet12)/the12(2)+&
					ri12*vet12*(su12**vet12)*u2* &
					(LOG(the12(2)*t3(i))+(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)
					
					
					u2=u2*(su12**vet12)*ri12*vet12+&
					ri12*vet12*the12(1)*(1-&
					gl12*vet12)*(su12**vet12)* &
					res2denum/the12(2)
					
					iter = iter +1
					
					res1((nvamax+iter))=u1*u3-&
					u2*v
					res1((nvamax+iter))=&
					res1((nvamax+iter))/(u1*u1)
					
					!write(6, *) " done x112/x212 -c7" 
				endif
				!la
				if(nva01nofix.gt.0) then 
					!ici
					u3=-LOG(the12(2))*res20112num-&
					res2the0112+res212num*LOG(the12(2))+&
					res2the12
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)* &
					res201num
					
					
					u2=res201num*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01*(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
					
					!write(6, *) " done x112/b01 -c7" 
				endif
				if(nva02nofix.gt.0) then 
					!la
					u3=-LOG(the12(2))*res20212num-&
					res2the0212
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri12*vet12*(su12**vet12)* &
					(LOG(the12(2)*t3(i))+&
					(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)* &
					(res202num)
					
					u2=-res202num*(su12**vet12)*ri12*vet12+&
					(1-gl02*vet02)*(su01**vet01)* &
					(su02**vet02)*ri02*vet02
					
				
					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					iter = iter +nva02nofix
					
					!write(6, *) " done x112/b02 -c7" 
					
				endif
				if(nva12nofix.gt.0) then 
					!la

					
					u3=LOG(the12(2))*res21212num+&
					res2the1212+LOG(the12(2))*res212num+&
					res2the12
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					(1-gl12*vet12)*(su12**vet12)* &
					(LOG(the12(2))*res212num+&
					res2the12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)* &
					(LOG(the12(2)*t3(i))+(1/the12(1))-&
					LOG(the12(2)*t3(i))*gl12*vet12)* &
					(res212num+res2denum*(1-gl12*vet12))-&
					LOG(the12(2)*t3(i))*gl12*vet12* &
					ri12*vet12*(su12**vet12)*res2denum
					
					u2=res212num*(su12**vet12)*ri12*vet12+&
					(1-gl12*vet12)*ri12*vet12* &
					(su12**vet12)*res2denum
				
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
					
					!write(6, *) " done x112/b12 -c7" 
				endif
			endif
			
			  !la
			!write(6, *) " done x112 -c7" 
			if(fix(6).eq.0)then
				iter = iter + 1
				nweib = nweib + 1
				
				v=the12(1)*res212num/the12(2)
				v=v*(su12**vet12)*ri12*vet12+&
				ri12*vet12*the12(1)*(1- &
				gl12*vet12)*(su12**vet12)* &
				res2denum/the12(2)
				
				res1(nweib)=v/u1
				
				u3=((the12(1)/the12(2))**2)*res21212num+&
				(the12(1)*(the12(1)-1)/((the12(2))**2))*res212num
				

				u3=u3*(su12**vet12)*ri12*vet12+&
				ri12*vet12*(su12**vet12)*the12(1)*(1-&
				gl12*vet12)*(2*the12(1)*res212num/the12(2)+&
				res2denum*the12(1)*(1-gl12*vet12)/the12(2)- &
				res2denum/the12(2))/the12(2)-&
				(the12(1)**2)*gl12*vet12*ri12*vet12* &
				(su12**vet12)*res2denum/(the12(2)**2)
				
				
				res1((nvamax+iter))=&
				u1*u3-v*v
				res1((nvamax+iter))=&
				res1((nvamax+iter))/(u1*u1)
				
				
				if(nva01nofix.gt.0) then 
					!la
					u3=the12(1)*res212num/the12(2)-&
					the12(1)*res20112num/the12(2)	
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*the12(1)*(1-&
					gl12*vet12)*res201num* &
					(su12**vet12)/the12(2)
					
					u2=res201num*(su12**vet12)*ri12*vet12-&
					ri02*vet02*gl01*vet01* &
					(su01**vet01)*(su02**vet02)
					
					
					res1((nvamax+iter+1):(nvamax+iter+nva01nofix))=&
					ve01nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva01nofix
				endif
				
				if(nva02nofix.gt.0) then 
					!la
					
					
					u3=-the12(1)*res20212num/the12(2)
					u3=u3*(su12**vet12)*ri12*vet12-&
					ri12*vet12*the12(1)*(1-&
					gl12*vet12)*(su12**vet12)* &
					res202num/the12(2)
					
					
					u2=-(su12**vet12)*ri12*vet12*res202num+&
				(1-gl02*vet02)*(su01**vet01)*(su02**vet02)*ri02*vet02

					res1((nvamax+iter+1):(nvamax+iter+nva02nofix))=&
					ve02nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva02nofix
				endif
				if(nva12nofix.gt.0) then 
				!la
					u3=res21212num*the12(1)/the12(2)+&
					the12(1)*res212num/the12(2)
					
					u3=u3*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(su12**vet12)*(1-&
					gl12*vet12)*the12(1)*res212num/the12(2)+&
					ri12*vet12*(su12**vet12)*(1-&
					gl12*vet12)*the12(1)* &
					(res212num+res2denum*(1-gl12*vet12))/ &
					the12(2)-ri12*vet12*gl12*vet12* &
					(su12**vet12)*res2denum*the12(1)/ &
					the12(2)
					
					u2=res212num*(su12**vet12)*ri12*vet12+&
					ri12*vet12*(1-&
					gl12*vet12)*(su12**vet12)*res2denum
							
					
					res1((nvamax+iter+1):(nvamax+iter+nva12nofix))=&
					ve12nofix(i,:)*(u1*u3-&
					v*u2)/(u1*u1)
					
					iter = iter +nva12nofix
				endif
			endif	
			!write(6, *) " done x212 -c7" 
				v=(su12**vet12)*ri12*vet12*res2denum+&
				(su01**vet01)*(su02**vet02)*ri02*vet02

				if(nva01nofix.gt.0) then 

				u1=-gl01*vet01*(su01**vet01)
				u1=u1*(su02**vet02)*ri02*vet02
				u1=u1+(su12**vet12)*ri12*vet12*res201num

				
				res1((nvaweib+1):(nvaweib+nva01nofix))=&
				ve01nofix(i,:)*u1/v
				res1((nvaweib+1):(nvaweib+nva01nofix))= &
				res1((nvaweib+1):(nvaweib+nva01nofix))+tronc01

				res1((nvamax12weib12+1):nvamax01)=&
				res20101num*(su12**vet12)*ri12*vet12
				resint=(su01**vet01)*(su02**vet02)*gl01*vet01
				resint=resint*(1-gl01*vet01)*ri02*vet02
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)-resint
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)*ve01square(i,:)
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)/v
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
				res1((nvamax12weib12+1):nvamax01)=&
				res1((nvamax12weib12+1):nvamax01)+tronc01square


				end if 

				if(nva02nofix.gt.0) then 

				u2=-(su12**vet12)*ri12*vet12*res202num+&
				(1-gl02*vet02)*(su01**vet01)*(su02**vet02)*ri02*vet02

				res1((nva01nofix+1+nvaweib):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1+nvaweib):nva0102)=&
				res1((nva01nofix+1+nvaweib):nva0102)+&
				tronc02

				res1((nvamax0112+1):nvamax02)=&
				-gl02*vet02*(su01**vet01)*(su02**vet02)*ri02*vet02
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)+&
				(su12**vet12)*ri12*vet12*(res20202num-res202num)+&
				((1-gl02*vet02)**2)*ri02*vet02*(su01**vet01)*(su02**vet02)


				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)*ve02square(i,:)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)/v
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)+tronc02square


				end if 

				if(nva12nofix.gt.0) then 

				u3=res212num+(1-gl12*vet12)*res2denum
				u3=u3*(su12**vet12)*ri12*vet12

				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v

				res1((nvamax0212+1):nvamax12)=&
				-gl12*vet12*res2denum+((1-gl12*vet12)**2)*res2denum+&
				(2*(1-gl12*vet12)*res212num+res212num)+&
				res21212num

				resint=ri12*vet12*(su12**vet12)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)*resint*ve12square(i,:)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)/v
				
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)


				end if 

				
				if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 

				kfix=nvamax01+1
				lfix=kfix-1+nva02nofix

				do j=1,nva01nofix
			    		res1(kfix:lfix)=(-res202num+&
					res20102num)*(su12**vet12)*ri12*vet12
					resint=-gl01*vet01*(1-gl02*vet02)*(su01**vet01)
					resint=resint*(su02**vet02)*ri02*vet02
					res1(kfix:lfix)=&
					res1(kfix:lfix)+resint
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			    		res1(kfix:lfix)=&
			    		res1(kfix:lfix)*v
					res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   	kfix=lfix+1
			   	lfix=lfix+nva02nofix
				end do

				end if 

				if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 

				
				kfix=nvamax0102+1
				lfix=kfix-1+nva12nofix

				do j=1,nva01nofix
			  		res1(kfix:lfix)=&
					(1-gl12*vet12)*res201num-res20112num+&
					res212num
					res1(kfix:lfix)=&
					res1(kfix:lfix)*ve01nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
					res1(kfix:lfix)*(su12**vet12)*ri12*vet12
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)*v
					res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva12nofix
				end do

				end if 


				if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 

				
				
        			kfix=nvamax02+1
				lfix=kfix-1+nva12nofix

				do j=1,nva02nofix
			  	res1(kfix:lfix)=res20212num+&
				(1-gl12*vet12)*res202num
				resint=(su12**vet12)
				res1(kfix:lfix)=&
				res1(kfix:lfix)*resint*ve02nofix(i,j)*ve12nofix(i,:)
				res1(kfix:lfix)=&
				-res1(kfix:lfix)*ri12*vet12
				res1(kfix:lfix)=&
				res1(kfix:lfix)*v
				res1(kfix:lfix)=&
				res1(kfix:lfix)-u2*u3*ve02nofix(i,j)*ve12nofix(i,:)
				res1(kfix:lfix)=&
			   	res1(kfix:lfix)/(v**2)
			   	kfix=lfix+1
			  	lfix=lfix+nva12nofix
				end do

				end if 
                            endif
                         endif                        
                      endif
                   endif  
                endif   
                endif   

                res = res + res1 


        end do   

!write(6, *) "end do" 
        likelihood_deriv = res




123     continue 
	 
	deallocate(b,bfix,fix,ve01,ve02,ve12,ve01nofix,&
	ve02nofix,ve12nofix,tronc01,tronc02,t0,t1,t2,t3,c,&
	ve01square,ve02square,ve12square,&
	tronc01square,tronc02square,troncweib01011beta01,&
	troncweib01012beta01,troncweib02021beta02,troncweib02022beta02)     

    end subroutine firstderivaweiballpara

!=============================================================================================  
!======================= Calculate first derivatives of loglik with weibull baseline risk ==========
!======================= only beta parameters ========================================================
!=============================================================================================  


subroutine firstderivaweib(b0,np0,npar0,bfix0,fix0,c0,no0,ve010,ve120,ve020,&
        dimnva01,dimnva12,dimnva02,nva01,nva12,nva02,t00,&
        t10,t20,t30,troncature0,likelihood_deriv)
	
	use commun
        implicit none
         
        double precision::res2denum,res201num,res202num,res212num,&
	vet01,vet12,vet02,resint,v,u1,u2,u3
        integer::np0,i,j,l,w,k,lfix, kfix,npar0,nva01,nva12,nva02,no0, &
	nz010,nz020,nz120,troncature0,dimnva01,dimnva02,dimnva12, & 
	nva01nofix,nva12nofix,nva02nofix,nvamax,sizespline,nva0102

	double precision,dimension(np0),intent(inout)::likelihood_deriv
	double precision,dimension(np0)::b0
	double precision,dimension(np0)::res,res1
        double precision,dimension(npar0)::bh
	double precision,dimension(npar0-np0)::bfix0
	integer,dimension(npar0)::fix0
	double precision,dimension(2)::the01,the12,the02
        double precision,dimension(no0,dimnva01)::ve010
	double precision,dimension(no0,dimnva02)::ve020
	double precision,dimension(no0,dimnva12)::ve120
	
	
        double precision::su01,ri01,su12,ri12,su02,ri02,gl01,gl02,gl12
	double precision,dimension(no0)::t00,t10,t20,t30
	integer,dimension(no0)::c0

	allocate(b(np0),bfix(npar0-np0),fix(npar0))
	b=b0
	bfix=bfix0
	fix=fix0
	troncature=troncature0

	sizespline=6


	if(nva01.gt.0) then 
	  nva01nofix=nva01-sum(fix((sizespline+1):(nva01+sizespline)))
	else 
	  nva01nofix=0
	end if 

	if(nva02.gt.0) then 
          nva02nofix=nva02-sum(fix((nva01+sizespline+1):(nva02+nva01+sizespline)))
	else 
	   nva02nofix=0
	end if 

	if(nva12.gt.0) then 
	  nva12nofix=nva12-sum(fix((nva01+nva02+sizespline+1):npar0))
	else 
	  nva12nofix=0
	end if 

	nva0102=nva01nofix+nva02nofix
	nvamax=nva01nofix+nva02nofix+nva12nofix
	

	if(nva01.gt.0) then 
		allocate(ve01(no0,nva01))
		allocate(ve01nofix(no0,nva01nofix))
		allocate(tronc01(nva01nofix))
		
	else 
		allocate(ve01(no0,1))
		allocate(ve01nofix(no0,1))
		ve01nofix=0
		allocate(tronc01(1))
	end if 
	
	if(nva02.gt.0) then 
		allocate(ve02(no0,nva02))
		allocate(ve02nofix(no0,nva02nofix))
		allocate(tronc02(nva02nofix))
		
	else 
		allocate(ve02(no0,1))
		allocate(ve02nofix(no0,1))
		ve02nofix=0
		allocate(tronc02(1))
	end if 

	if(nva12.gt.0) then 
		allocate(ve12(no0,nva12))
		allocate(ve12nofix(no0,nva12nofix))
		
	else 
		allocate(ve12(no0,1))
		allocate(ve12nofix(no0,1))
		ve12nofix=0
	end if 



	ve01=ve010
	ve02=ve020
	ve12=ve120

	allocate(t0(no0),t1(no0),t2(no0),t3(no0),c(no0))
	c=c0
	t0=t00
	t1=t10
	t2=t20
	t3=t30

         
        ! we need to put bh at its original values if in posfix 


       l=0
       lfix=0
       w=0

	if(nva01.gt.0) then
	do k=1,nva01
	   if(fix((sizespline+k)).eq.0) then 
		lfix=lfix+1
		ve01nofix(:,lfix)=ve01(:,k)
	   end if 
	end do
	end if 

	lfix=0

	if(nva02.gt.0) then 
	do k=1,nva02
	   if(fix((sizespline+nva01+k)).eq.0) then 
		lfix=lfix+1
		ve02nofix(:,lfix)=ve02(:,k)
	   end if 
	end do
	end if 
	
	lfix=0

	if(nva12.gt.0) then 
	do k=1,nva12
	   if(fix((sizespline+nva01+nva02+k)).eq.0) then 
		lfix=lfix+1
		ve12nofix(:,lfix)=ve12(:,k)
	   end if 
	end do
	end if
	
	l=0
       lfix=0
       w=0

     do k=1,npar0 
         if(fix(k).eq.0) then
            l=l+1
            bh(k)=b(l)
	 end if 
         if(fix(k).eq.1) then
            w=w+1
            bh(k)=bfix(w)
         end if
      end do

   
	



         do i=1,2
            the01(i)=(bh(i))*(bh(i))
         end do
         do i=1,2
            j = 2+i
            the02(i)=(bh(j))*(bh(j))
         end do
         do i=1,2
            j = 4+i
            the12(i)=(bh(j))*(bh(j))
         end do

!---------- calcul des derivees premiere ------------------   

	res = 0
        do i=1,no0


                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
				vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
				vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
				vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

		
                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)

                res1 = 0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc01 = 0
				tronc02 =  0
                        else 
				call fonct(t0(i),the01,ri01,gl01,su01)
				call fonct(t0(i),the02,ri02,gl02,su02)
                                tronc01=ve01nofix(i,:)*gl01*vet01
                        	tronc02=ve02nofix(i,:)*gl02*vet02
                        end if
                else
                        tronc01 = 0
                  	tronc02 = 0
                end if
		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
			call fonct(t1(i),the01,ri01,gl01,su01)
			call fonct(t1(i),the02,ri02,gl02,su02)

			if(nva01nofix.gt.0) then 
			
			res1(1:nva01nofix)=&
			-ve01nofix(i,:)*gl01*vet01+tronc01
			end if 

			if(nva02nofix.gt.0) then
			res1((nva01nofix+1):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02
			end if 

			if(nva12nofix.gt.0) then 
			res1((nva0102+1):nvamax)=0
			end if 
			

                else
                if(c(i).eq.2)then ! cpi 0-->1


			call fonct(t3(i),the12,ri12,gl12,su12)
			call qgaussweibfirstderiv(t1(i),t2(i),the01,&
			the02,the12,res2denum,res201num,&
			res202num,res212num,&
			vet01,vet02,vet12)
                        
			v=res2denum*(su12**vet12)

			if(nva01nofix.gt.0) then

			u1=res201num*(su12**vet12)
			
      			res1(1:nva01nofix)=&
			ve01nofix(i,:)*u1/v
			res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01
			end if 

			
			if(nva02nofix.gt.0) then

			u2=-res202num*(su12**vet12)
			res1((nva01nofix+1):nva0102)=&
			ve02nofix(i,:)*u2/v
			res1((nva01nofix+1):nva0102)=&
			res1((nva01nofix+1):nva0102)+&
			tronc02
			end if 

			
			if(nva12nofix.gt.0) then 

			u3=res212num-gl12*vet12*res2denum
			u3=u3*(su12**vet12)

			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v
			end if 
			
			
                else  
                    if(c(i).eq.3)then ! obs 0-->1

			call fonct(t1(i),the01,ri01,gl01,su01)
			call fonct(t1(i),the02,ri02,gl02,su02)
			call fonct(t1(i),the12,ri12,gl12,su12)

			if(nva01nofix.gt.0) then 

			res1(1:nva01nofix)=-ve01nofix(i,:)*gl01*vet01+&
			tronc01+ve01nofix(i,:)
			end if 

			if(nva02nofix.gt.0) then 

			res1((nva01nofix+1):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02
			end if 

			if(nva12nofix.gt.0) then 

			res1((nva0102+1):nvamax)=gl12*vet12
			call fonct(t3(i),the12,ri12,gl12,su12)
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)-gl12*vet12
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)*ve12nofix(i,:)

			end if 
			
			

                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			
			call fonct(t3(i),the12,ri12,gl12,su12)
			call qgaussweibfirstderiv(t1(i),t2(i),the01,the02,the12,&
			res2denum,res201num,res202num,res212num,&
			vet01,vet02,vet12)

			v=res2denum*(su12**vet12)*ri12*vet12

			if(nva01nofix.gt.0) then 
			u1=res201num*(su12**vet12)*ri12*vet12
			
			res1(1:nva01nofix)=&
			ve01nofix(i,:)*u1/v
			res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01
			end if 

			
			if(nva02nofix.gt.0) then 

			u2=-res202num*(su12**vet12)*ri12*vet12
			
			res1((nva01nofix+1):(nva01nofix+nva02nofix))=&
			ve02nofix(i,:)*u2/v
			res1((nva01nofix+1):nva0102)=&
			res1((nva01nofix+1):nva0102)+&
			tronc02
			end if 


			
			if(nva12nofix.gt.0) then 

			u3=res212num+&
			(1-gl12*vet12)*res2denum
			u3=u3*(su12**vet12)*ri12*vet12
			
			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v

			end if 

                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
				call fonct(t1(i),the01,ri01,gl01,su01)
				call fonct(t1(i),the02,ri02,gl02,su02)
				call fonct(t1(i),the12,ri12,gl12,su12)
				
				if(nva01nofix.gt.0) then 

				res1(1:nva01nofix)=&
				-ve01nofix(i,:)*gl01*vet01+&
				ve01nofix(i,:)+&
				tronc01
				end if 

				if(nva02nofix.gt.0) then 

				res1((nva01nofix+1):nva0102)=&
				-ve02nofix(i,:)*gl02*vet02+tronc02
				end if 
			
				if(nva12nofix.gt.0) then 
				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*gl12*vet12+&
				ve12nofix(i,:)

				call fonct(t3(i),the12,ri12,gl12,su12)
				res1((nva0102+1):nvamax)=&
				res1((nva0102+1):nvamax)-ve12nofix(i,:)*gl12*vet12
				end if 
				
				
				
                         else
                            if(c(i).eq.6)then ! vivant ???


				call fonct(t3(i),the01,ri01,gl01,su01)
				call fonct(t3(i),the02,ri02,gl02,su02)
				call fonct(t3(i),the12,ri12,gl12,su12)
				call qgaussweibfirstderiv(t1(i),t3(i),the01,the02,the12,&
				res2denum,res201num,res202num,res212num,&
				vet01,vet02,vet12)

				
				v=(su12**vet12)*res2denum+&
				(su01**vet01)*(su02**vet02)

				if(nva01nofix.gt.0) then 

				u1=(-gl01*vet01)*(su01**vet01)*(su02**vet02)+&
				(su12**vet12)*res201num
			        res1(1:nva01nofix)=&
				ve01nofix(i,:)*u1/v
				res1(1:nva01nofix)=&
				res1(1:nva01nofix)+tronc01
				end if 
				
				if(nva02nofix.gt.0) then 

				u2=-gl02*vet02*(su01**vet01)*(su02**vet02)
				u2=u2-(su12**vet12)*res202num
				res1((nva01nofix+1):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1):nva0102)=&
				res1((nva01nofix+1):nva0102)+&
				tronc02
				end if 

				if(nva12nofix.gt.0) then 
				
				u3=-gl12*vet12*(su12**vet12)*res2denum+&
				(su12**vet12)*res212num
				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v
				end if 
                        	

				
                            else ! passage 0-->2  
				
				call fonct(t3(i),the01,ri01,gl01,su01)
				call fonct(t3(i),the02,ri02,gl02,su02)
				call fonct(t3(i),the12,ri12,gl12,su12)
        			call qgaussweibfirstderiv(t1(i),t3(i),the01,the02,&
        			the12,res2denum,res201num,res202num,res212num,&
        			vet01,vet02,vet12)
				
				v=(su12**vet12)*ri12*vet12*res2denum+&
				(su01**vet01)*(su02**vet02)*ri02*vet02

				if(nva01nofix.gt.0) then 

				u1=-gl01*vet01*(su01**vet01)
				u1=u1*(su02**vet02)*ri02*vet02
				u1=u1+(su12**vet12)*ri12*vet12*res201num

				
				res1(1:nva01nofix)=&
				ve01nofix(i,:)*u1/v
				res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01
				end if 

				if(nva02nofix.gt.0) then 

				u2=-(su12**vet12)*ri12*vet12*res202num+&
				(1-gl02*vet02)*(su01**vet01)*(su02**vet02)*ri02*vet02

				res1((nva01nofix+1):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1):nva0102)=&
				res1((nva01nofix+1):nva0102)+&
				tronc02

				end if 

				if(nva12nofix.gt.0) then 

				u3=res212num+(1-gl12*vet12)*res2denum
				u3=u3*(su12**vet12)*ri12*vet12

				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v

				end if 

                            endif
                         endif                        
                      endif
                   endif  
                endif   
                endif   

                res = res + res1 


		


        end do   


        likelihood_deriv = res




123     continue 
	 
	deallocate(b,bfix,fix,ve01,ve02,ve12,ve01nofix,&
	ve02nofix,ve12nofix,tronc01,tronc02,t0,t1,t2,t3,c)    

    end subroutine firstderivaweib

!=============================================================================================  
!======================= Calculate derivatives of loglik with weibull baseline risk ==========
!======================= only diagnola terms for second derivatives of beta parameters             ==========
!=============================================================================================  


subroutine derivaweibdiag(b0,np0,npar0,bfix0,fix0,c0,no0,ve010,ve120,ve020,&
        dimnva01,dimnva12,dimnva02,nva01,nva12,nva02,t00,&
        t10,t20,t30,troncature0,likelihood_deriv)
	
	use commun
        implicit none
         
        double precision::res2denum,res201num,res202num,res212num, &
	res20101num,res20202num,res21212num,vet01,vet12,vet02, & 
	resint,v,u1,u2,u3
        integer::np0,i,j,l,w,k,lfix, kfix,npar0,nva01,nva12,nva02,no0, &
	nz010,nz020,nz120,troncature0,dimnva01,dimnva02,dimnva12, & 
	nva01nofix,nva12nofix,nva02nofix,nvamax, sizespline,nva0102,&
	nvamax01,nvamax02,nvamax12

	double precision,dimension(2*np0),intent(inout)::likelihood_deriv
	double precision,dimension(np0)::b0
	double precision,dimension(2*np0)::res,res1
        double precision,dimension(npar0)::bh
	double precision,dimension(npar0-np0)::bfix0
	integer,dimension(npar0)::fix0
	double precision,dimension(2)::the01,the12,the02
        double precision,dimension(no0,dimnva01)::ve010
	double precision,dimension(no0,dimnva02)::ve020
	double precision,dimension(no0,dimnva12)::ve120
	
	
        double precision::su01,ri01,su12,ri12,su02,ri02,gl01,gl02,gl12
	double precision,dimension(no0)::t00,t10,t20,t30
	integer,dimension(no0)::c0

!	PRINT *,'allocate'
	allocate(b(np0),bfix(npar0-np0),fix(npar0))
	b=b0
	bfix=bfix0
	fix=fix0
	troncature=troncature0

	sizespline=6

    
!	 PRINT *,'start'
	if(nva01.gt.0) then 
	  nva01nofix=nva01-sum(fix((sizespline+1):(nva01+sizespline)))
	else 
	  nva01nofix=0
	end if 

	if(nva02.gt.0) then 
          nva02nofix=nva02-sum(fix((nva01+sizespline+1):(nva02+nva01+sizespline)))
	else 
	   nva02nofix=0
	end if 

	if(nva12.gt.0) then 
	  nva12nofix=nva12-sum(fix((nva01+nva02+sizespline+1):npar0))
	else 
	  nva12nofix=0
	end if 

	nva0102=nva01nofix+nva02nofix
	nvamax=nva01nofix+nva02nofix+nva12nofix
	nvamax01=nvamax+nva01nofix
	nvamax02=nvamax01+nva02nofix
	nvamax12=nvamax02+nva12nofix


	if(nva01.gt.0) then 
		allocate(ve01(no0,nva01))
		allocate(ve01nofix(no0,nva01nofix))
		allocate(ve01square(no0,nva01nofix))
		allocate(tronc01(nva01nofix))
		allocate(tronc01square(nva01nofix))
	else 
		allocate(ve01(no0,1))
		allocate(ve01nofix(no0,1))
		ve01nofix=0
		allocate(ve01square(no0,1))
		ve01square=0
		allocate(tronc01(1))
		allocate(tronc01square(1))
	end if 
	
	if(nva02.gt.0) then 
		allocate(ve02(no0,nva02))
		allocate(ve02nofix(no0,nva02nofix))
		allocate(ve02square(no0,nva02nofix))
		allocate(tronc02(nva02nofix))
		allocate(tronc02square(nva02nofix))
	else 
		allocate(ve02(no0,1))
		allocate(ve02nofix(no0,1))
		ve02nofix=0
		allocate(ve02square(no0,1))
		ve02square=0
		allocate(tronc02(1))
		allocate(tronc02square(1))
	end if 

	if(nva12.gt.0) then 
		allocate(ve12(no0,nva12))
		allocate(ve12nofix(no0,nva12nofix))
		allocate(ve12square(no0,nva12nofix))
	else 
		allocate(ve12(no0,1))
		allocate(ve12nofix(no0,1))
		ve12nofix=0
		allocate(ve12square(no0,1))
		ve12square=0
	end if 



	ve01=ve010
	ve02=ve020
	ve12=ve120

	allocate(t0(no0),t1(no0),t2(no0),t3(no0),c(no0))
	c=c0
	t0=t00
	t1=t10
	t2=t20
	t3=t30

         
        ! we need to put bh at its original values if in posfix 


       l=0
       lfix=0
       w=0

	if(nva01.gt.0) then
	do k=1,nva01
	   if(fix((sizespline+k)).eq.0) then 
		lfix=lfix+1
		ve01nofix(:,lfix)=ve01(:,k)
	   end if 
	end do
	ve01square=ve01nofix*ve01nofix
	end if 

	lfix=0

	if(nva02.gt.0) then 
	do k=1,nva02
	   if(fix((sizespline+nva01+k)).eq.0) then 
		lfix=lfix+1
		ve02nofix(:,lfix)=ve02(:,k)
	   end if 
	end do
	ve02square=ve02nofix*ve02nofix
	end if 
	
	lfix=0

	if(nva12.gt.0) then 
	do k=1,nva12
	   if(fix((sizespline+nva01+nva02+k)).eq.0) then 
		lfix=lfix+1
		ve12nofix(:,lfix)=ve12(:,k)
	   end if 
	end do
	ve12square=ve12nofix*ve12nofix
	end if
	
	l=0
       lfix=0
       w=0

     do k=1,npar0 
         if(fix(k).eq.0) then
            l=l+1
            bh(k)=b(l)
	 end if 
         if(fix(k).eq.1) then
            w=w+1
            bh(k)=bfix(w)
         end if
      end do

!   PRINT *,'update ve'
	



         do i=1,2
            the01(i)=(bh(i))*(bh(i))
         end do
         do i=1,2
            j = 2+i
            the02(i)=(bh(j))*(bh(j))
         end do
         do i=1,2
            j = 4+i
            the12(i)=(bh(j))*(bh(j))
         end do




 !	PRINT *,'start loop'
	
!---------- calcul des derivees premiere ------------------   

	res = 0
        do i=1,no0


                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
				vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
				vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
				vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

		
                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)

                res1 = 0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc01 = 0
				tronc02 =  0
                        	tronc01square= 0
                        	tronc02square=0
                        else 
				call fonct(t0(i),the01,ri01,gl01,su01)
				call fonct(t0(i),the02,ri02,gl02,su02)
                                tronc01=ve01nofix(i,:)*gl01*vet01
                        	tronc02=ve02nofix(i,:)*gl02*vet02
                        	tronc01square=ve01square(i,:)*gl01*vet01
                        	tronc02square=ve02square(i,:)*gl02*vet02
                        end if
                else
                        tronc01 = 0
                  	tronc02 = 0
                  	tronc01square=0
                  	tronc02square=0
                end if
		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
!			 PRINT *,'profile 1'
			call fonct(t1(i),the01,ri01,gl01,su01)
			call fonct(t1(i),the02,ri02,gl02,su02)

			if(nva01nofix.gt.0) then 
			
			res1(1:nva01nofix)=&
			-ve01nofix(i,:)*gl01*vet01+tronc01
			res1((nvamax+1):nvamax01)=&
			-ve01square(i,:)*gl01*vet01+&
			tronc01square

			end if 

			if(nva02nofix.gt.0) then
			res1((nva01nofix+1):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02

			res1((nvamax01+1):nvamax02)=&
			-ve02square(i,:)*gl02*vet02+&
			tronc02square
			end if 

			if(nva12nofix.gt.0) then 
			res1((nva0102+1):nvamax)=0
			res1((nvamax02+1):nvamax12)=0
			end if 
			

                else
                if(c(i).eq.2)then ! cpi 0-->1

!			PRINT *,'profile 2'
			call fonct(t3(i),the12,ri12,gl12,su12)
			call qgaussweibderivdiag(t1(i),t2(i),the01,&
			the02,the12,res2denum,res201num,&
			res202num,res212num,res20101num,&
			res20202num,res21212num,&
			vet01,vet02,vet12)
                        
			v=res2denum*(su12**vet12)

			if(nva01nofix.gt.0) then

			u1=res201num*(su12**vet12)
			
      			res1(1:nva01nofix)=&
			ve01nofix(i,:)*u1/v
			res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01

			res1((nvamax+1):nvamax01)=&
			ve01square(i,:)*res20101num*(su12**vet12)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)/v
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)+&
			tronc01square

			end if 

			
			if(nva02nofix.gt.0) then

			u2=-res202num*(su12**vet12)
			res1((nva01nofix+1):nva0102)=&
			ve02nofix(i,:)*u2/v
			res1((nva01nofix+1):nva0102)=&
			res1((nva01nofix+1):nva0102)+&
			tronc02

			res1((nvamax01+1):nvamax02)=&
			-ve02square(i,:)*(su12**vet12)*(-res20202num+res202num)
			res1((nvamax01+1):nvamax02)=&
			res1((nvamax01+1):nvamax02)/v
			res1((nvamax01+1):nvamax02)=&
			res1((nvamax01+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
			res1((nvamax01+1):nvamax02)=&
			res1((nvamax01+1):nvamax02)+tronc02square


			end if 

			
			if(nva12nofix.gt.0) then 

			u3=res212num-gl12*vet12*res2denum
			u3=u3*(su12**vet12)

			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v
			res1((nvamax02+1):nvamax12)=&
			-2*gl12*vet12*res212num+res212num+&
			res21212num-gl12*vet12*(1-gl12*vet12)*res2denum
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)*(su12**vet12)*ve12square(i,:)
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)/v
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)-ve12square(i,:)*((u3/v)**2)

			
			end if 
			
                else  
                    if(c(i).eq.3)then ! obs 0-->1
!			 PRINT *,'profile 3'
			call fonct(t1(i),the01,ri01,gl01,su01)
			call fonct(t1(i),the02,ri02,gl02,su02)
			call fonct(t1(i),the12,ri12,gl12,su12)

			if(nva01nofix.gt.0) then 

			res1(1:nva01nofix)=-ve01nofix(i,:)*gl01*vet01+&
			tronc01+ve01nofix(i,:)
			res1((nvamax+1):nvamax01)=&
			-ve01square(i,:)*gl01*vet01+&
			tronc01square

			end if 

			if(nva02nofix.gt.0) then 

			res1((nva01nofix+1):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02
			res1((nvamax01+1):nvamax02)=&
			-ve02square(i,:)*gl02*vet02+&
			tronc02square

			end if 

			if(nva12nofix.gt.0) then 

			res1((nva0102+1):nvamax)=gl12*vet12
			res1((nvamax02+1):nvamax12)=gl12*vet12

			call fonct(t3(i),the12,ri12,gl12,su12)
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)-gl12*vet12
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)*ve12nofix(i,:)

			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)-gl12*vet12
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)*ve12square(i,:)
			

			end if 
			
			

                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
!			PRINT *,'profile 4'
			call fonct(t3(i),the12,ri12,gl12,su12)
			call qgaussweibderivdiag(t1(i),t2(i),the01,&
			the02,the12,res2denum,res201num,&
			res202num,res212num,&
			res20101num,res20202num,&
			res21212num,vet01,vet02,vet12)

			v=res2denum*(su12**vet12)*ri12*vet12

			if(nva01nofix.gt.0) then 
			u1=res201num*(su12**vet12)*ri12*vet12
			
			res1(1:nva01nofix)=&
			ve01nofix(i,:)*u1/v
			res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01

			res1((nvamax+1):nvamax01)=&
			res20101num*(su12**vet12)*ri12*vet12*ve01square(i,:)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)/v
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)+tronc01square

			end if 

			
			if(nva02nofix.gt.0) then 

			u2=-res202num*(su12**vet12)*ri12*vet12
			
			res1((nva01nofix+1):(nva01nofix+nva02nofix))=&
			ve02nofix(i,:)*u2/v
			res1((nva01nofix+1):nva0102)=&
			res1((nva01nofix+1):nva0102)+&
			tronc02
			
			res1((nvamax01+1):nvamax02)=&
			(res20202num-res202num)*ri12*vet12*(su12**vet12)
			res1((nvamax01+1):nvamax02)=&
			res1((nvamax01+1):nvamax02)*ve02square(i,:)	
			res1((nvamax01+1):nvamax02)=&
			res1((nvamax01+1):nvamax02)/v
			res1((nvamax01+1):nvamax02)=&
			res1((nvamax01+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
			res1((nvamax01+1):nvamax02)=&
			res1((nvamax01+1):nvamax02)+tronc02square

			end if 


			
			if(nva12nofix.gt.0) then 

			u3=res212num+&
			(1-gl12*vet12)*res2denum
			u3=u3*(su12**vet12)*ri12*vet12
			
			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v
			
			res1((nvamax02+1):nvamax12)=&
			(1-gl12*vet12)*res212num+res21212num+&
			((1-gl12*vet12)**2)*res2denum+&
			(2-gl12*vet12)*res212num-gl12*vet12*res2denum
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)*(su12**vet12)*ri12*vet12
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)*ve12square(i,:)
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)/v
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)-ve12square(i,:)*((u3/v)**2)			


			end if 


                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
!				PRINT *,'profile 5'
				call fonct(t1(i),the01,ri01,gl01,su01)
				call fonct(t1(i),the02,ri02,gl02,su02)
				call fonct(t1(i),the12,ri12,gl12,su12)
				
				if(nva01nofix.gt.0) then 

				res1(1:nva01nofix)=&
				-ve01nofix(i,:)*gl01*vet01+&
				ve01nofix(i,:)+&
				tronc01
				res1((nvamax+1):nvamax01)=&
				-ve01square(i,:)*gl01*vet01+&
				tronc01square

				end if 

				if(nva02nofix.gt.0) then 

				res1((nva01nofix+1):nva0102)=&
				-ve02nofix(i,:)*gl02*vet02+tronc02

				res1((nvamax01+1):nvamax02)=&
				-ve02square(i,:)*gl02*vet02+&
				tronc02square

				end if 
			
				if(nva12nofix.gt.0) then 

				res1((nvamax02+1):nvamax12)=&
				gl12*vet12
				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*gl12*vet12+&
				ve12nofix(i,:)

				call fonct(t3(i),the12,ri12,gl12,su12)
				res1((nva0102+1):nvamax)=&
				res1((nva0102+1):nvamax)-ve12nofix(i,:)*gl12*vet12
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)-gl12*vet12
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)*ve12square(i,:)
				

				end if 
				
				
				
                         else
                            if(c(i).eq.6)then ! vivant ???

				!write(6,*)'profile 6'
				call fonct(t3(i),the01,ri01,gl01,su01)
				call fonct(t3(i),the02,ri02,gl02,su02)
				call fonct(t3(i),the12,ri12,gl12,su12)
				call qgaussweibderivdiag(t1(i),t3(i),the01,&
				the02,the12,res2denum,res201num,res202num,&
				res212num,res20101num,res20202num,&
				res21212num,vet01,vet02,vet12)

				
				v=(su12**vet12)*res2denum+&
				(su01**vet01)*(su02**vet02)

				if(nva01nofix.gt.0) then 

				u1=(-gl01*vet01)*(su01**vet01)*(su02**vet02)+&
				(su12**vet12)*res201num
			        res1(1:nva01nofix)=&
				ve01nofix(i,:)*u1/v
				res1(1:nva01nofix)=&
				res1(1:nva01nofix)+tronc01

				res1((nvamax+1):nvamax01)=&
				res20101num*(su12**vet12)
				resint=(su01**vet01)*(su02**vet02)
				resint=resint*gl01*vet01*(1-gl01*vet01)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-resint
				
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)*ve01square(i,:)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)/v

				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)+tronc01square

				end if 
				
				if(nva02nofix.gt.0) then 

				u2=-gl02*vet02*(su01**vet01)*(su02**vet02)
				u2=u2-(su12**vet12)*res202num
				res1((nva01nofix+1):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1):nva0102)=&
				res1((nva01nofix+1):nva0102)+&
				tronc02

				res1((nvamax01+1):nvamax02)=&
				(su12**vet12)*(res20202num-res202num)-&
				gl02*vet02*(su01**vet01)*(su02**vet02)*(1-gl02*vet02)
				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)*ve02square(i,:)
				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)/v

				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)+tronc02square

				
				end if 

				if(nva12nofix.gt.0) then 
				
				u3=-gl12*vet12*(su12**vet12)*res2denum+&
				(su12**vet12)*res212num
				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v

				res1((nvamax02+1):nvamax12)=&
				res21212num+(1-gl12*vet12)*res212num-&
				gl12*vet12*(1-gl12*vet12)*res2denum-&
				gl12*vet12*res212num
				resint=(su12**vet12)
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)*resint*ve12square(i,:)
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)/v
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)-ve12square(i,:)*((u3/v)**2)


				end if 
                        	
				
                            else ! passage 0-->2  
!				PRINT *,'profile 7'
				call fonct(t3(i),the01,ri01,gl01,su01)
				call fonct(t3(i),the02,ri02,gl02,su02)
				call fonct(t3(i),the12,ri12,gl12,su12)
        			call qgaussweibderivdiag(t1(i),t3(i),&
				the01,the02,&
        			the12,res2denum,res201num,res202num,res212num,&
        			res20101num,res20202num,res21212num,&
        			vet01,vet02,vet12)
				
				
				v=(su12**vet12)*ri12*vet12*res2denum+&
				(su01**vet01)*(su02**vet02)*ri02*vet02

				if(nva01nofix.gt.0) then 

				u1=-gl01*vet01*(su01**vet01)
				u1=u1*(su02**vet02)*ri02*vet02
				u1=u1+(su12**vet12)*ri12*vet12*res201num

				
				res1(1:nva01nofix)=&
				ve01nofix(i,:)*u1/v
				res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01

				res1((nvamax+1):nvamax01)=&
				res20101num*(su12**vet12)*ri12*vet12
				resint=(su01**vet01)*(su02**vet02)*gl01*vet01
				resint=resint*(1-gl01*vet01)*ri02*vet02
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-resint
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)*ve01square(i,:)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)/v
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)+tronc01square


				end if 

				if(nva02nofix.gt.0) then 

				u2=-(su12**vet12)*ri12*vet12*res202num+&
				(1-gl02*vet02)*(su01**vet01)*(su02**vet02)*ri02*vet02

				res1((nva01nofix+1):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1):nva0102)=&
				res1((nva01nofix+1):nva0102)+&
				tronc02

				res1((nvamax01+1):nvamax02)=&
				-gl02*vet02*(su01**vet01)*(su02**vet02)*ri02*vet02
				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)+&
				(su12**vet12)*ri12*vet12*(res20202num-res202num)+&
				((1-gl02*vet02)**2)*ri02*vet02*(su01**vet01)*(su02**vet02)


				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)*ve02square(i,:)
				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)/v
				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)+tronc02square


				end if 

				if(nva12nofix.gt.0) then 

				u3=res212num+(1-gl12*vet12)*res2denum
				u3=u3*(su12**vet12)*ri12*vet12

				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v

				res1((nvamax02+1):nvamax12)=&
				-gl12*vet12*res2denum+((1-gl12*vet12)**2)*res2denum+&
				2*(1-gl12*vet12)*res212num+res212num+&
				res21212num

				resint=ri12*vet12*(su12**vet12)
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)*resint*ve12square(i,:)
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)/v
				
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)-ve12square(i,:)*((u3/v)**2)


				end if 

				
				
                            endif
                         endif                        
                      endif
                   endif  
                endif   
                endif   

                res = res + res1 


		


        end do   


        likelihood_deriv = res




123     continue 
	 
	deallocate(b,bfix,fix,ve01,ve02,ve12,ve01nofix,&
	ve02nofix,ve12nofix,tronc01,tronc02,t0,t1,t2,t3,c,&
	ve01square,ve02square,ve12square,&
	tronc01square,tronc02square)     

    end subroutine derivaweibdiag

!attention commencer une ligne par (- non possible 

!=============================================================================================  
!======================= Calculate derivatives of loglik with M-splines baseline risk ========
 !======================= only beta parameters ========================================================
!============================================================================================= 


 subroutine derivaspline(b0,np0,npar0,bfix0,fix0,zi010,zi120,&
      zi020,c0,no0,nz010,nz120,nz020,ve010,ve120,ve020,&
        dimnva01,dimnva12,dimnva02,nva01,nva12,nva02,t00,&
        t10,t20,t30,troncature0,likelihood_deriv)
	
	use commun
        implicit none
         
        double precision::res2denum,res201num,res202num,res212num, &
	res20101num,res20102num,res20112num,res20202num,res20212num,&
	res21212num,vet01,vet12,vet02,resint,v,u1,u2,u3
        integer::np0,i,j,l,w,k,lfix, kfix,npar0,nva01,nva12,nva02,no0, &
	nz010,nz020,nz120,troncature0,dimnva01,dimnva02,dimnva12, & 
	nva01nofix,nva12nofix,nva02nofix,nvamax, sizespline,nva0102,&
	nvamax01,nvamax0102,nvamax0112,nvamax02,nvamax0212,nvamax12

	double precision,dimension(np0+np0*(np0+1)/2),intent(inout)::likelihood_deriv
	double precision,dimension(np0)::b0
	double precision,dimension(np0+np0*(np0+1)/2)::res,res1
        double precision,dimension(npar0)::bh
	double precision,dimension(npar0-np0)::bfix0
	integer,dimension(npar0)::fix0
	double precision,dimension(-2:(nz010+3))::zi010
	double precision,dimension(-2:(nz020+3))::zi020
	double precision,dimension(-2:(nz120+3))::zi120
	double precision,dimension(-2:(nz010-1))::the01
	double precision,dimension(-2:(nz120-1))::the12
	double precision,dimension(-2:(nz020-1))::the02
        double precision,dimension(no0,dimnva01)::ve010
	double precision,dimension(no0,dimnva02)::ve020
	double precision,dimension(no0,dimnva12)::ve120
	
	
        double precision::su01,ri01,su12,ri12,su02,ri02,gl01,gl02,gl12
	double precision,dimension(no0)::t00,t10,t20,t30
	integer,dimension(no0)::c0

	allocate(b(np0),bfix(npar0-np0),fix(npar0))
	b=b0
	bfix=bfix0
	fix=fix0
	allocate(zi01(-2:(nz01+3)),zi12(-2:(nz12+3)),zi02(-2:(nz02+3)))
	zi01=zi010
	zi02=zi020
	zi12=zi120

	
	nz01=nz010
	nz02=nz020
	nz12=nz120
	sizespline=nz01+nz02+nz12+6
	troncature=troncature0


	if(nva01.gt.0) then 
	  nva01nofix=nva01-sum(fix((sizespline+1):(nva01+sizespline)))
	else 
	  nva01nofix=0
	end if 

	if(nva02.gt.0) then 
          nva02nofix=nva02-sum(fix((nva01+sizespline+1):(nva02+nva01+sizespline)))
	else 
	   nva02nofix=0
	end if 

	if(nva12.gt.0) then 
	  nva12nofix=nva12-sum(fix((nva01+nva02+sizespline+1):npar0))
	else 
	  nva12nofix=0
	end if 

	nva0102=nva01nofix+nva02nofix
	nvamax=nva01nofix+nva02nofix+nva12nofix
	nvamax01=nvamax+(nva01nofix+1)*nva01nofix/2
	nvamax0102=nvamax01+nva01nofix*nva02nofix
	nvamax0112=nvamax0102+nva01nofix*nva12nofix
	nvamax02=nvamax0112+(nva02nofix+1)*nva02nofix/2
	nvamax0212=nvamax02+nva02nofix*nva12nofix
	nvamax12=nvamax0212+(nva12nofix+1)*nva12nofix/2


	if(nva01.gt.0) then 
		allocate(ve01(no0,nva01))
		allocate(ve01nofix(no0,nva01nofix))
		allocate(ve01square(no0,nva01nofix*(nva01nofix+1)/2))
		allocate(tronc01(nva01nofix))
		allocate(tronc01square(nva01nofix*(nva01nofix+1)/2))
	else 
		allocate(ve01(no0,1))
		allocate(ve01nofix(no0,1))
		allocate(ve01square(no0,1))
		ve01square=0
		ve01nofix=0
		allocate(tronc01(1))
		allocate(tronc01square(1))
	end if 
	
	if(nva02.gt.0) then 
		allocate(ve02(no0,nva02))
		allocate(ve02nofix(no0,nva02nofix))
		allocate(ve02square(no0,nva02nofix*(nva02nofix+1)/2))
		allocate(tronc02(nva02nofix))
		allocate(tronc02square(nva02nofix*(nva02nofix+1)/2))
	else 
		allocate(ve02(no0,1))
		allocate(ve02nofix(no0,1))
		allocate(ve02square(no0,1))
		ve02square=0
		ve02nofix=0
		allocate(tronc02(1))
		allocate(tronc02square(1))
	end if 

	if(nva12.gt.0) then 
		allocate(ve12(no0,nva12))
		allocate(ve12nofix(no0,nva12nofix))
		allocate(ve12square(no0,nva12nofix*(nva12nofix+1)/2))
	else 
		allocate(ve12(no0,1))
		allocate(ve12nofix(no0,1))
		allocate(ve12square(no0,1))
	end if 


	ve01=ve010
	ve02=ve020
	ve12=ve120

	allocate(t0(no0),t1(no0),t2(no0),t3(no0),c(no0))
	c=c0
	t0=t00
	t1=t10
	t2=t20
	t3=t30

         
        ! we need to put bh at its original values if in posfix 


       l=0
       lfix=0
       w=0


	if(nva01.gt.0) then 
	do k=1,nva01
	   if(fix((sizespline+k)).eq.0) then 
		lfix=lfix+1
		ve01nofix(:,lfix)=ve01(:,k)
	   end if 
	end do

	
	do i=1,no0
		lfix=1
		do l=1,nva01nofix
	   	ve01square(i,lfix:(lfix+nva01nofix-l))=&
		ve01nofix(i,l:nva01nofix)
		ve01square(i,lfix:(lfix+nva01nofix-l))=&
		ve01square(i,lfix:(lfix+nva01nofix-l))*ve01nofix(i,l)
		lfix=lfix+nva01nofix-l+1
		end do
	end do
	 
	end if 


	lfix=0
	if(nva02.gt.0) then 

	do k=1,nva02
	   if(fix((sizespline+nva01+k)).eq.0) then 
		lfix=lfix+1
		ve02nofix(:,lfix)=ve02(:,k)
	   end if 
	end do

	do i=1,no0
		lfix=1
		do l=1,nva02nofix
	   	ve02square(i,lfix:(lfix+nva02nofix-l))=&
		ve02nofix(i,l)*ve02nofix(i,l:nva02nofix)
		lfix=lfix+nva02nofix-l+1
		end do
	end do

	end if 

	
	lfix=0

	if(nva12.gt.0) then 
	do k=1,nva12
	   if(fix((sizespline+nva01+nva02+k)).eq.0) then 
		lfix=lfix+1
		ve12nofix(:,lfix)=ve12(:,k)
	   end if 
	end do

	do i=1,no0
		lfix=1
		do l=1,nva12nofix
	   		ve12square(i,lfix:(lfix+nva12nofix-l))=&
			ve12nofix(i,l)*ve12nofix(i,l:nva12nofix)
			lfix=lfix+nva12nofix-l+1
		end do
	end do

	end if 
	
	l=0
       lfix=0
       w=0

     do k=1,npar0 
         if(fix(k).eq.0) then
            l=l+1
            bh(k)=b(l)
	 end if 
         if(fix(k).eq.1) then
            w=w+1
            bh(k)=bfix(w)
         end if
      end do
   
	



         do i=1,nz01+2
            the01(i-3)=(bh(i))*(bh(i))
!       the01(i-3)=dexp(bh(i))
         end do
         do i=1,nz02+2
            j = nz01+2+i
            the02(i-3)=(bh(j))*(bh(j))
!       the12(i-3)=dexp(bh(j))
         end do
         do i=1,nz12+2
            j = nz02+2+nz01+2+i
            the12(i-3)=(bh(j))*(bh(j))
!       the02(i-3)=dexp(bh(j))
         end do


	
!---------- calcul des derivees premiere ------------------   

	res = 0
        do i=1,no0


                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
				vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
				vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
				vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

		
                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)

                res1 = 0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc01 = 0
				tronc02 =  0
                        	tronc01square= 0
                        	tronc02square=0
                        else 
				
				call susp(t0(i),the01,nz01,su01,ri01,zi01,gl01)
				call susp(t0(i),the02,nz02,su02,ri02,zi02,gl02)
                                tronc01=ve01nofix(i,:)*(gl01*vet01)
                        	tronc02=ve02nofix(i,:)*(gl02*vet02)
                        	tronc01square=ve01square(i,:)*(gl01*vet01)
                        	tronc02square=ve02square(i,:)*(gl02*vet02)
                        end if
                else
                        tronc01 = 0
                  	tronc02 = 0
                  	tronc01square=0
                  	tronc02square=0
                end if
		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
			call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
			call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)

			if(nva01nofix.gt.0) then 

			res1(1:nva01nofix)=&
			-ve01nofix(i,:)*gl01*vet01+tronc01
			res1((nvamax+1):nvamax01)=&
			-ve01square(i,:)*gl01*vet01+&
			tronc01square

			if(nva02nofix.gt.0) then 
			res1((nvamax01+1):nvamax0102)=0
			end if 

			if(nva12nofix.gt.0) then 
			res1((nvamax0102+1):nvamax0112)=0
			end if 

			end if

			if(nva02nofix.gt.0) then 

			res1((nva01nofix+1):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02
			res1((nvamax0112+1):nvamax02)=&
			-ve02square(i,:)*gl02*vet02+&
			tronc02square
			
			if(nva12nofix.gt.0) then 
			res1((nvamax02+1):nvamax0212)=0
			end if 

			end if 

			if(nva12nofix.gt.0) then 

			res1((nva0102+1):nvamax)=0
			res1((nvamax0212+1):nvamax12)=0

			end if 

		
			

                else
                if(c(i).eq.2)then ! cpi 0-->1
			call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
			call qgausssplinederiv(t1(i),t2(i),the01,&
			the02,the12,res2denum,res201num,&
			res202num,res212num,res20101num,&
			res20102num,res20112num,res20202num,&
			res20212num,res21212num,&
			vet01,vet02,vet12)
                        
			
			v=res2denum*(su12**vet12)

			if(nva01nofix.gt.0) then 

			u1=res201num*(su12**vet12)
      			res1(1:nva01nofix)=&
			ve01nofix(i,:)*u1/v
			res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01
			res1((nvamax+1):nvamax01)=&
			ve01square(i,:)*res20101num*(su12**vet12)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)/v
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)+&
			tronc01square

			end if 

			if(nva02nofix.gt.0) then 

			u2=-res202num*(su12**vet12)
			res1((nva01nofix+1):nva0102)=&
			ve02nofix(i,:)*u2/v
			res1((nva01nofix+1):nva0102)=&
			res1((nva01nofix+1):nva0102)+&
			tronc02

			res1((nvamax0112+1):nvamax02)=&
			-ve02square(i,:)*(su12**vet12)*(-res20202num+res202num)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)/v
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)+tronc02square

			end if 

			if(nva12nofix.gt.0) then

			u3=res212num-gl12*vet12*res2denum
			u3=u3*(su12**vet12)
			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v

			res1((nvamax0212+1):nvamax12)=&
			-2*gl12*vet12*res212num+res212num+&
			res21212num-gl12*vet12*(1-gl12*vet12)*res2denum
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*(su12**vet12)*ve12square(i,:)
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)/v
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)

			end if 
			

			
			if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 

			kfix=nvamax01+1
			lfix=kfix-1+nva02nofix

			do j=1,nva01nofix

			   res1(kfix:lfix)=&
			  (res20102num-res202num)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva02nofix

			end do

			end if 

			if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 

			kfix=nvamax0102+1
			lfix=kfix-1+nva12nofix

			do j=1,nva01nofix
			   res1(kfix:lfix)=-gl12*vet12*res201num-res20112num+&
			   res212num
			   res1(kfix:lfix)=res1(kfix:lfix)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do

			end if 

			
			if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 

			kfix=nvamax02+1
			lfix=kfix-1+nva12nofix

			do j=1,nva02nofix
			   res1(kfix:lfix)=-res20212num+&
			   gl12*vet12*res202num
			   res1(kfix:lfix)=res1(kfix:lfix)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve02nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u2*u3*ve02nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do

			end if 

			
                else  
                    if(c(i).eq.3)then ! obs 0-->1

			call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
			call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
			call susp(t1(i),the12,nz12,su12,ri12,zi12,gl12)


			if(nva01nofix.gt.0) then 

			res1(1:nva01nofix)=-ve01nofix(i,:)*gl01*vet01+&
			tronc01+ve01nofix(i,:)
			res1((nvamax+1):nvamax01)=&
			-ve01square(i,:)*gl01*vet01+&
			tronc01square

			if(nva02nofix.gt.0) then 
			res1((nvamax01+1):nvamax0102)=0
			end if 
	
			if(nva12nofix.gt.0) then 
			res1((nvamax0102+1):nvamax0112)=0
			end if 

			end if 

			if(nva02nofix.gt.0) then 

			res1((nva01nofix+1):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02
			res1((nvamax0112+1):nvamax02)=&
			-ve02square(i,:)*gl02*vet02+&
			tronc02square

			if(nva12nofix.gt.0) then 
			res1((nvamax02+1):nvamax0212)=0
			end if 

			end if 

			if(nva12nofix.gt.0) then 

			res1((nva0102+1):nvamax)=gl12*vet12
			res1((nvamax0212+1):nvamax12)=gl12*vet12
			call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)-gl12*vet12
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)*ve12nofix(i,:)
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)-gl12*vet12
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*ve12square(i,:)
			

			end if 

			

                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			
			call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
			call qgausssplinederiv(t1(i),t2(i),the01,the02,the12,&
			res2denum,res201num,res202num,res212num,&
			res20101num,res20102num,&
			res20112num,res20202num,res20212num,&
			res21212num,vet01,vet02,vet12)

			
			v=res2denum*(su12**vet12)*ri12*vet12

			if(nva01nofix.gt.0) then 

			u1=res201num*(su12**vet12)*ri12*vet12
			res1(1:nva01nofix)=&
			ve01nofix(i,:)*u1/v
			res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01

			res1((nvamax+1):nvamax01)=&
			res20101num*(su12**vet12)*ri12*vet12*ve01square(i,:)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)/v
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)+tronc01square
			

			end if 

			if(nva02nofix.gt.0) then

			u2=-res202num*(su12**vet12)*ri12*vet12
			res1((nva01nofix+1):(nva01nofix+nva02nofix))=&
			ve02nofix(i,:)*u2/v
			res1((nva01nofix+1):nva0102)=&
			res1((nva01nofix+1):nva0102)+&
			tronc02

			res1((nvamax0112+1):nvamax02)=&
			(res20202num-res202num)*ri12*vet12*(su12**vet12)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)*ve02square(i,:)	
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)/v
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
			res1((nvamax0112+1):nvamax02)=&
			res1((nvamax0112+1):nvamax02)+tronc02square


			end if 

			if(nva12nofix.gt.0) then

			u3=res212num+&
			(1-gl12*vet12)*res2denum
			u3=u3*(su12**vet12)*ri12*vet12
			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v

			res1((nvamax0212+1):nvamax12)=&
			(1-gl12*vet12)*res212num+res21212num+&
			((1-gl12*vet12)**2)*res2denum+&
			(2-gl12*vet12)*res212num-gl12*vet12*res2denum
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*(su12**vet12)*ri12*vet12
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)*ve12square(i,:)
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)/v
			res1((nvamax0212+1):nvamax12)=&
			res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)			


			end if 

			

			if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 

			kfix=nvamax01+1
			lfix=kfix-1+nva02nofix

			do j=1,nva01nofix
			   res1(kfix:lfix)=(-res202num+res20102num)*ri12*vet12
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*(su12**vet12)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva02nofix
			end do

			end if 


			if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 

			kfix=nvamax0102+1
			lfix=kfix-1+nva12nofix

			do j=1,nva01nofix
			   res1(kfix:lfix)=&
			   (1-gl12*vet12)*res201num+&
			   res212num-res20112num
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ve01nofix(i,j)*ve12nofix(i,:)
		           res1(kfix:lfix)=&
			   res1(kfix:lfix)*ri12*vet12*(su12**vet12)
			   res1(kfix:lfix)=res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do


			end if 

			
			if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 

			kfix=nvamax02+1
			lfix=kfix-1+nva12nofix

			do j=1,nva02nofix
			   res1(kfix:lfix)=&
			   res20212num+(1-gl12*vet12)*res202num
			   res1(kfix:lfix)=&
			   -res1(kfix:lfix)*ve02nofix(i,j)*ve12nofix(i,:)
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)*ri12*vet12*(su12**vet12)
			   res1(kfix:lfix)=res1(kfix:lfix)*v
			   res1(kfix:lfix)=&
			   res1(kfix:lfix)-ve02nofix(i,j)*ve12nofix(i,:)*u2*u3
			   res1(kfix:lfix)=res1(kfix:lfix)/(v**2)
			   kfix=lfix+1
			   lfix=lfix+nva12nofix
			end do

			end if 

		
                        
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
				
			         call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
			         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
			         call susp(t1(i),the12,nz12,su12,ri12,zi12,gl12)

				if(nva01nofix.gt.0) then 

				res1(1:nva01nofix)=&
				-ve01nofix(i,:)*gl01*vet01+&
				ve01nofix(i,:)+&
				tronc01
				res1((nvamax+1):nvamax01)=&
				-ve01square(i,:)*gl01*vet01+&
				tronc01square
			
				if(nva02nofix.gt.0) then
				res1((nvamax01+1):nvamax0102)=0
				end if 

				if(nva12nofix.gt.0) then 
				res1((nvamax0102+1):nvamax0112)=0
				end if 

				end if 

				if(nva02nofix.gt.0) then 

				res1((nva01nofix+1):nva0102)=&
				-ve02nofix(i,:)*gl02*vet02+tronc02
				res1((nvamax0112+1):nvamax02)=&
				-ve02square(i,:)*gl02*vet02+&
				tronc02square
				
				if(nva12nofix.gt.0) then 
				res1((nvamax02+1):nvamax0212)=0
				end if 

				end if 

				if(nva12nofix.gt.0) then

				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*gl12*vet12+&
				ve12nofix(i,:)
				res1((nvamax0212+1):nvamax12)=&
				gl12*vet12
				call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                res1((nva0102+1):nvamax)=&
				res1((nva0102+1):nvamax)-ve12nofix(i,:)*gl12*vet12
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)-gl12*vet12
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)*ve12square(i,:)

				end if 
				
                         else
                            if(c(i).eq.6)then ! vivant ???
				
				call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
				call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
			        call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
				call qgausssplinederiv(t1(i),t3(i),the01,the02,the12,&
				res2denum,res201num,res202num,res212num,res20101num,&
				res20102num,res20112num,res20202num,res20212num,&
				res21212num,vet01,vet02,vet12)

				
				v=(su12**vet12)*res2denum+&
				(su01**vet01)*(su02**vet02)
			
				if(nva01nofix.gt.0) then 

				u1=(-gl01*vet01)*(su01**vet01)*(su02**vet02)+&
				(su12**vet12)*res201num
			        res1(1:nva01nofix)=&
				ve01nofix(i,:)*u1/v
				res1(1:nva01nofix)=&
				res1(1:nva01nofix)+tronc01

				res1((nvamax+1):nvamax01)=&
				res20101num*(su12**vet12)
				resint=(su01**vet01)*(su02**vet02)
				resint=resint*gl01*vet01*(1-gl01*vet01)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-resint
				
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)*ve01square(i,:)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)/v

				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)+tronc01square

				end if 


				if(nva02nofix.gt.0) then 

				u2=-gl02*vet02*(su01**vet01)*(su02**vet02)
				u2=u2-(su12**vet12)*res202num
				res1((nva01nofix+1):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1):nva0102)=&
				res1((nva01nofix+1):nva0102)+&
				tronc02

				res1((nvamax0112+1):nvamax02)=&
				(su12**vet12)*(res20202num-res202num)-&
				gl02*vet02*(su01**vet01)*(su02**vet02)*(1-gl02*vet02)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)*ve02square(i,:)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)/v

				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)+tronc02square

				end if 

				if(nva12nofix.gt.0) then 
				
				u3=-gl12*vet12*(su12**vet12)*res2denum+&
				(su12**vet12)*res212num
				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v

				res1((nvamax0212+1):nvamax12)=&
				res21212num+(1-gl12*vet12)*res212num-&
				gl12*vet12*(1-gl12*vet12)*res2denum-&
				gl12*vet12*res212num
				resint=(su12**vet12)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)*resint*ve12square(i,:)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)/v
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)
			
				end if 
                        	
				

				if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 

				kfix=nvamax01+1
				lfix=kfix-1+nva02nofix

				do j=1,nva01nofix
			   		res1(kfix:lfix)=&
					-(res202num-res20102num)*(su12**vet12)
					resint=(su01**vet01)*(su02**vet02)
					res1(kfix:lfix)=res1(kfix:lfix)+&
					gl01*vet01*gl02*vet02*resint
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)*v
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva02nofix
				end do

				end if 

				if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 

				kfix=nvamax0102+1
				lfix=kfix-1+nva12nofix

				do j=1,nva01nofix
			  		res1(kfix:lfix)=-gl12*vet12*res201num+&
					res212num-res20112num
					resint=(su12**vet12)
			   		res1(kfix:lfix)=&
			  		res1(kfix:lfix)*resint*ve01nofix(i,j)*ve12nofix(i,:)
		           		res1(kfix:lfix)=&
			   		res1(kfix:lfix)*v
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
			   		res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva12nofix
				end do

				end if 

			


				if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 

				kfix=nvamax02+1
				lfix=kfix-1+nva12nofix

				do j=1,nva02nofix
			  		res1(kfix:lfix)=&
					-res20212num+gl12*vet12*res202num
					resint=(su12**vet12)
					res1(kfix:lfix)=&
					res1(kfix:lfix)*resint*ve02nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
					res1(kfix:lfix)*v
					res1(kfix:lfix)=&
					res1(kfix:lfix)-u2*u3*ve02nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva12nofix
				end do

				end if 

				

                            else ! passage 0-->2  
				
				call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
				call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
			        call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
        			call qgausssplinederiv(t1(i),t3(i),the01,the02,&
        			the12,res2denum,res201num,res202num,res212num,&
        			res20101num,res20102num,&
			  	res20112num,res20202num,res20212num,res21212num,&
        			vet01,vet02,vet12)


				
				v=(su12**vet12)*ri12*vet12*res2denum+&
				(su01**vet01)*(su02**vet02)*ri02*vet02


				if(nva01nofix.gt.0) then 

				u1=-gl01*vet01*(su01**vet01)
				u1=u1*(su02**vet02)*ri02*vet02
				u1=u1+(su12**vet12)*ri12*vet12*res201num

				res1(1:nva01nofix)=&
				ve01nofix(i,:)*u1/v
				res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01

				res1((nvamax+1):nvamax01)=&
				res20101num*(su12**vet12)*ri12*vet12
				resint=(su01**vet01)*(su02**vet02)*gl01*vet01
				resint=resint*(1-gl01*vet01)*ri02*vet02
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-resint
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)*ve01square(i,:)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)/v
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)+tronc01square

				end if 


				if(nva02nofix.gt.0) then 

				u2=-(su12**vet12)*ri12*vet12*res202num+&
				(1-gl02*vet02)*(su01**vet01)*(su02**vet02)*ri02*vet02

				res1((nva01nofix+1):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1):nva0102)=&
				res1((nva01nofix+1):nva0102)+&
				tronc02

				res1((nvamax0112+1):nvamax02)=&
				-gl02*vet02*(su01**vet01)*(su02**vet02)*ri02*vet02
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)+&
				(su12**vet12)*ri12*vet12*(res20202num-res202num)+&
				((1-gl02*vet02)**2)*ri02*vet02*(su01**vet01)*(su02**vet02)


				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)*ve02square(i,:)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)/v
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
				res1((nvamax0112+1):nvamax02)=&
				res1((nvamax0112+1):nvamax02)+tronc02square

				end if 

				if(nva12nofix.gt.0) then 

				u3=res212num+(1-gl12*vet12)*res2denum
				u3=u3*(su12**vet12)*ri12*vet12

				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v

				res1((nvamax0212+1):nvamax12)=&
				-gl12*vet12*res2denum+((1-gl12*vet12)**2)*res2denum+&
				2*(1-gl12*vet12)*res212num+res212num+&
				res21212num

				resint=ri12*vet12*(su12**vet12)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)*resint*ve12square(i,:)
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)/v
				
				res1((nvamax0212+1):nvamax12)=&
				res1((nvamax0212+1):nvamax12)-ve12square(i,:)*((u3/v)**2)

				end if 


				

				if(nva01nofix.gt.0 .AND. nva02nofix.gt.0) then 

				kfix=nvamax01+1
				lfix=kfix-1+nva02nofix

				do j=1,nva01nofix
			    		res1(kfix:lfix)=(-res202num+&
					res20102num)*(su12**vet12)*ri12*vet12
					resint=-gl01*vet01*(1-gl02*vet02)*(su01**vet01)
					resint=resint*(su02**vet02)*ri02*vet02
					res1(kfix:lfix)=&
					res1(kfix:lfix)+resint
			   		res1(kfix:lfix)=&
					res1(kfix:lfix)*ve01nofix(i,j)*ve02nofix(i,:)
			    		res1(kfix:lfix)=&
			    		res1(kfix:lfix)*v
					res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u2*ve01nofix(i,j)*ve02nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   	kfix=lfix+1
			   	lfix=lfix+nva02nofix
				end do

				end if 


				if(nva01nofix.gt.0 .AND. nva12nofix.gt.0) then 

				kfix=nvamax0102+1
				lfix=kfix-1+nva12nofix

				do j=1,nva01nofix
			  		res1(kfix:lfix)=&
					(1-gl12*vet12)*res201num-res20112num+&
					res212num
					res1(kfix:lfix)=&
					res1(kfix:lfix)*ve01nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
					res1(kfix:lfix)*(su12**vet12)*ri12*vet12
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)*v
					res1(kfix:lfix)=&
					res1(kfix:lfix)-u1*u3*ve01nofix(i,j)*ve12nofix(i,:)
					res1(kfix:lfix)=&
			   		res1(kfix:lfix)/(v**2)
			   		kfix=lfix+1
			   		lfix=lfix+nva12nofix
				end do

				end if 



				if(nva12nofix.gt.0 .AND. nva02nofix.gt.0) then 

        			kfix=nvamax02+1
				lfix=kfix-1+nva12nofix

				do j=1,nva02nofix
			  	res1(kfix:lfix)=res20212num+&
				(1-gl12*vet12)*res202num
				resint=(su12**vet12)
				res1(kfix:lfix)=&
				res1(kfix:lfix)*resint*ve02nofix(i,j)*ve12nofix(i,:)
				res1(kfix:lfix)=&
				-res1(kfix:lfix)*ri12*vet12
				res1(kfix:lfix)=&
				res1(kfix:lfix)*v
				res1(kfix:lfix)=&
				res1(kfix:lfix)-u2*u3*ve02nofix(i,j)*ve12nofix(i,:)
				res1(kfix:lfix)=&
			   	res1(kfix:lfix)/(v**2)
			   	kfix=lfix+1
			  	lfix=lfix+nva12nofix
				end do

				end if 

                            endif
                         endif                        
                      endif
                   endif  
                endif   
                endif   

                res = res + res1 

        end do   

        likelihood_deriv = res




123     continue 
	 
	deallocate(b,bfix,fix,zi01,zi02,zi12,ve01,ve02,ve12,ve01nofix,&
	ve02nofix,ve12nofix,tronc01,tronc02,t0,t1,t2,t3,c,&
	ve01square,ve02square,ve12square,&
	tronc01square,tronc02square)     

    end subroutine derivaspline



!=============================================================================================  
!======================= Calculate first derivatives of loglik with M-splines baseline risk ==
!======================= only beta parameters ========================================================
!============================================================================================= 


 subroutine firstderivaspline(b0,np0,npar0,bfix0,fix0,zi010,zi120,&
      zi020,c0,no0,nz010,nz120,nz020,ve010,ve120,ve020,&
        dimnva01,dimnva12,dimnva02,nva01,nva12,nva02,t00,&
        t10,t20,t30,troncature0,likelihood_deriv)
	
	use commun
        implicit none
         
        double precision::res2denum,res201num,res202num,res212num,&
	vet01,vet12,vet02,resint,v,u1,u2,u3
        integer::np0,i,j,l,w,k,lfix, kfix,npar0,nva01,nva12,nva02,no0, &
	nz010,nz020,nz120,troncature0,dimnva01,dimnva02,dimnva12, & 
	nva01nofix,nva12nofix,nva02nofix,nvamax,sizespline,nva0102

	double precision,dimension(np0),intent(inout)::likelihood_deriv
	double precision,dimension(np0)::b0
	double precision,dimension(np0)::res,res1
        double precision,dimension(npar0)::bh
	double precision,dimension(npar0-np0)::bfix0
	integer,dimension(npar0)::fix0
	double precision,dimension(-2:(nz010+3))::zi010
	double precision,dimension(-2:(nz020+3))::zi020
	double precision,dimension(-2:(nz120+3))::zi120
	double precision,dimension(-2:(nz010-1))::the01
	double precision,dimension(-2:(nz120-1))::the12
	double precision,dimension(-2:(nz020-1))::the02
        double precision,dimension(no0,dimnva01)::ve010
	double precision,dimension(no0,dimnva02)::ve020
	double precision,dimension(no0,dimnva12)::ve120
	
	
        double precision::su01,ri01,su12,ri12,su02,ri02,gl01,gl02,gl12
	double precision,dimension(no0)::t00,t10,t20,t30
	integer,dimension(no0)::c0

	allocate(b(np0),bfix(npar0-np0),fix(npar0))
	b=b0
	bfix=bfix0
	fix=fix0
	allocate(zi01(-2:(nz01+3)),zi12(-2:(nz12+3)),zi02(-2:(nz02+3)))
	zi01=zi010
	zi02=zi020
	zi12=zi120

	
	nz01=nz010
	nz02=nz020
	nz12=nz120
	sizespline=nz01+nz02+nz12+6
	troncature=troncature0


	if(nva01.gt.0) then 
	  nva01nofix=nva01-sum(fix((sizespline+1):(nva01+sizespline)))
	else 
	  nva01nofix=0
	end if 

	if(nva02.gt.0) then 
          nva02nofix=nva02-sum(fix((nva01+sizespline+1):(nva02+nva01+sizespline)))
	else 
	   nva02nofix=0
	end if 

	if(nva12.gt.0) then 
	  nva12nofix=nva12-sum(fix((nva01+nva02+sizespline+1):npar0))
	else 
	  nva12nofix=0
	end if 

	nva0102=nva01nofix+nva02nofix
	nvamax=nva01nofix+nva02nofix+nva12nofix
	

	if(nva01.gt.0) then 
		allocate(ve01(no0,nva01))
		allocate(ve01nofix(no0,nva01nofix))
		allocate(tronc01(nva01nofix))
	else 
		allocate(ve01(no0,1))
		allocate(ve01nofix(no0,1))
		ve01nofix=0
		allocate(tronc01(1))
	end if 
	
	if(nva02.gt.0) then 
		allocate(ve02(no0,nva02))
		allocate(ve02nofix(no0,nva02nofix))
		allocate(tronc02(nva02nofix))
	else 
		allocate(ve02(no0,1))
		allocate(ve02nofix(no0,1))
		ve02nofix=0
		allocate(tronc02(1))
	end if 

	if(nva12.gt.0) then 
		allocate(ve12(no0,nva12))
		allocate(ve12nofix(no0,nva12nofix))
	else 
		allocate(ve12(no0,1))
		allocate(ve12nofix(no0,1))
	end if 


	ve01=ve010
	ve02=ve020
	ve12=ve120

	allocate(t0(no0),t1(no0),t2(no0),t3(no0),c(no0))
	c=c0
	t0=t00
	t1=t10
	t2=t20
	t3=t30

         
        ! we need to put bh at its original values if in posfix 


       l=0
       lfix=0
       w=0


	if(nva01.gt.0) then 
	do k=1,nva01
	   if(fix((sizespline+k)).eq.0) then 
		lfix=lfix+1
		ve01nofix(:,lfix)=ve01(:,k)
	   end if 
	end do
	 
	end if 


	lfix=0
	if(nva02.gt.0) then 

	do k=1,nva02
	   if(fix((sizespline+nva01+k)).eq.0) then 
		lfix=lfix+1
		ve02nofix(:,lfix)=ve02(:,k)
	   end if 
	end do
	end if 

	
	lfix=0

	if(nva12.gt.0) then 
	do k=1,nva12
	   if(fix((sizespline+nva01+nva02+k)).eq.0) then 
		lfix=lfix+1
		ve12nofix(:,lfix)=ve12(:,k)
	   end if 
	end do

	end if 
	
	l=0
       lfix=0
       w=0

     do k=1,npar0 
         if(fix(k).eq.0) then
            l=l+1
            bh(k)=b(l)
	 end if 
         if(fix(k).eq.1) then
            w=w+1
            bh(k)=bfix(w)
         end if
      end do
   
	



         do i=1,nz01+2
            the01(i-3)=(bh(i))*(bh(i))
!       the01(i-3)=dexp(bh(i))
         end do
         do i=1,nz02+2
            j = nz01+2+i
            the02(i-3)=(bh(j))*(bh(j))
!       the12(i-3)=dexp(bh(j))
         end do
         do i=1,nz12+2
            j = nz02+2+nz01+2+i
            the12(i-3)=(bh(j))*(bh(j))
!       the02(i-3)=dexp(bh(j))
         end do


	
!---------- calcul des derivees premiere ------------------   

	res = 0
        do i=1,no0


                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
				vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
				vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
				vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

		
                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)

                res1 = 0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc01 = 0
				tronc02 =  0
                        else 
				
				call susp(t0(i),the01,nz01,su01,ri01,zi01,gl01)
				call susp(t0(i),the02,nz02,su02,ri02,zi02,gl02)
                                tronc01=ve01nofix(i,:)*(gl01*vet01)
                        	tronc02=ve02nofix(i,:)*(gl02*vet02)
                        end if
                else
                        tronc01 = 0
                  	tronc02 = 0
                end if
		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
			call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
			call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)

			if(nva01nofix.gt.0) then 

			res1(1:nva01nofix)=&
			-ve01nofix(i,:)*gl01*vet01+tronc01

			end if

			if(nva02nofix.gt.0) then 

			res1((nva01nofix+1):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02
			
			end if 

			if(nva12nofix.gt.0) then 

			res1((nva0102+1):nvamax)=0

			end if 

		
			

                else
                if(c(i).eq.2)then ! cpi 0-->1
			call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
			call qgausssplinefirstderiv(t1(i),t2(i),the01,&
			the02,the12,res2denum,res201num,&
			res202num,res212num,&
			vet01,vet02,vet12)
                        
			
			v=res2denum*(su12**vet12)

			if(nva01nofix.gt.0) then 

			u1=res201num*(su12**vet12)
      			res1(1:nva01nofix)=&
			ve01nofix(i,:)*u1/v
			res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01
			end if 

			if(nva02nofix.gt.0) then 

			u2=-res202num*(su12**vet12)
			res1((nva01nofix+1):nva0102)=&
			ve02nofix(i,:)*u2/v
			res1((nva01nofix+1):nva0102)=&
			res1((nva01nofix+1):nva0102)+&
			tronc02
			end if 

			if(nva12nofix.gt.0) then

			u3=res212num-gl12*vet12*res2denum
			u3=u3*(su12**vet12)
			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v
			end if 
			

			
                else  
                    if(c(i).eq.3)then ! obs 0-->1

			call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
			call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
			call susp(t1(i),the12,nz12,su12,ri12,zi12,gl12)


			if(nva01nofix.gt.0) then 

			res1(1:nva01nofix)=-ve01nofix(i,:)*gl01*vet01+&
			tronc01+ve01nofix(i,:)
			end if 

			if(nva02nofix.gt.0) then 

			res1((nva01nofix+1):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02
			end if 

			if(nva12nofix.gt.0) then 

			res1((nva0102+1):nvamax)=gl12*vet12
			call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)-gl12*vet12
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)*ve12nofix(i,:)
			end if 

			

                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			
			call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
			call qgausssplinefirstderiv(t1(i),t2(i),the01,the02,the12,&
			res2denum,res201num,res202num,res212num,&
			vet01,vet02,vet12)

			
			v=res2denum*(su12**vet12)*ri12*vet12

			if(nva01nofix.gt.0) then 

			u1=res201num*(su12**vet12)*ri12*vet12
			res1(1:nva01nofix)=&
			ve01nofix(i,:)*u1/v
			res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01
			end if 

			if(nva02nofix.gt.0) then

			u2=-res202num*(su12**vet12)*ri12*vet12
			res1((nva01nofix+1):(nva01nofix+nva02nofix))=&
			ve02nofix(i,:)*u2/v
			res1((nva01nofix+1):nva0102)=&
			res1((nva01nofix+1):nva0102)+&
			tronc02
			end if 

			if(nva12nofix.gt.0) then

			u3=res212num+&
			(1-gl12*vet12)*res2denum
			u3=u3*(su12**vet12)*ri12*vet12
			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v
			end if 

                        
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
				
			         call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
			         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
			         call susp(t1(i),the12,nz12,su12,ri12,zi12,gl12)

				if(nva01nofix.gt.0) then 

				res1(1:nva01nofix)=&
				-ve01nofix(i,:)*gl01*vet01+&
				ve01nofix(i,:)+&
				tronc01
				end if 

				if(nva02nofix.gt.0) then 

				res1((nva01nofix+1):nva0102)=&
				-ve02nofix(i,:)*gl02*vet02+tronc02
				end if 

				if(nva12nofix.gt.0) then

				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*gl12*vet12+&
				ve12nofix(i,:)
				call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                res1((nva0102+1):nvamax)=&
				res1((nva0102+1):nvamax)-ve12nofix(i,:)*gl12*vet12
				end if 
				
                         else
                            if(c(i).eq.6)then ! vivant ???
				
				call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
				call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
			        call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
				call qgausssplinefirstderiv(t1(i),t3(i),the01,the02,the12,&
				res2denum,res201num,res202num,res212num,&
				vet01,vet02,vet12)

				
				v=(su12**vet12)*res2denum+&
				(su01**vet01)*(su02**vet02)
			
				if(nva01nofix.gt.0) then 

				u1=(-gl01*vet01)*(su01**vet01)*(su02**vet02)+&
				(su12**vet12)*res201num
			        res1(1:nva01nofix)=&
				ve01nofix(i,:)*u1/v
				res1(1:nva01nofix)=&
				res1(1:nva01nofix)+tronc01

				end if 


				if(nva02nofix.gt.0) then 

				u2=-gl02*vet02*(su01**vet01)*(su02**vet02)
				u2=u2-(su12**vet12)*res202num
				res1((nva01nofix+1):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1):nva0102)=&
				res1((nva01nofix+1):nva0102)+&
				tronc02

				end if 

				if(nva12nofix.gt.0) then 
				
				u3=-gl12*vet12*(su12**vet12)*res2denum+&
				(su12**vet12)*res212num
				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v

				end if 
                        	
                            else ! passage 0-->2  
				
				call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
				call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
			        call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
        			call qgausssplinefirstderiv(t1(i),t3(i),the01,the02,&
        			the12,res2denum,res201num,res202num,res212num,&
        			vet01,vet02,vet12)

				v=(su12**vet12)*ri12*vet12*res2denum+&
				(su01**vet01)*(su02**vet02)*ri02*vet02


				if(nva01nofix.gt.0) then 

				u1=-gl01*vet01*(su01**vet01)
				u1=u1*(su02**vet02)*ri02*vet02
				u1=u1+(su12**vet12)*ri12*vet12*res201num

				res1(1:nva01nofix)=&
				ve01nofix(i,:)*u1/v
				res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01

				end if 


				if(nva02nofix.gt.0) then 

				u2=-(su12**vet12)*ri12*vet12*res202num+&
				(1-gl02*vet02)*(su01**vet01)*(su02**vet02)*ri02*vet02

				res1((nva01nofix+1):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1):nva0102)=&
				res1((nva01nofix+1):nva0102)+&
				tronc02
				end if 

				if(nva12nofix.gt.0) then 

				u3=res212num+(1-gl12*vet12)*res2denum
				u3=u3*(su12**vet12)*ri12*vet12

				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v
				end if 


                            endif
                         endif                        
                      endif
                   endif  
                endif   
                endif   

                res = res + res1 

        end do   

        likelihood_deriv = res




123     continue 
	 
	deallocate(b,bfix,fix,zi01,zi02,zi12,ve01,ve02,ve12,ve01nofix,&
	ve02nofix,ve12nofix,tronc01,tronc02,t0,t1,t2,t3,c)     

    end subroutine firstderivaspline



!=============================================================================================  
!======================= Calculate derivatives of loglik with M-splines baseline risk ========
!======================= only diagonal terms of hessian of beta parameters ======================================
!============================================================================================= 


 subroutine derivasplinediag(b0,np0,npar0,bfix0,fix0,zi010,zi120,&
      zi020,c0,no0,nz010,nz120,nz020,ve010,ve120,ve020,&
        dimnva01,dimnva12,dimnva02,nva01,nva12,nva02,t00,&
        t10,t20,t30,troncature0,likelihood_deriv)
	
	use commun
        implicit none
         
        double precision::res2denum,res201num,res202num,res212num, &
	res20101num,res20202num,res21212num,vet01,vet12,vet02,resint,v,u1,u2,u3
        integer::np0,i,j,l,w,k,lfix, kfix,npar0,nva01,nva12,nva02,no0, &
	nz010,nz020,nz120,troncature0,dimnva01,dimnva02,dimnva12, & 
	nva01nofix,nva12nofix,nva02nofix,nvamax, sizespline,nva0102,&
	nvamax01,nvamax0102,nvamax0112,nvamax02,nvamax0212,nvamax12

	double precision,dimension(np0+np0),intent(inout)::likelihood_deriv
	double precision,dimension(np0)::b0
	double precision,dimension(np0+np0)::res,res1
        double precision,dimension(npar0)::bh
	double precision,dimension(npar0-np0)::bfix0
	integer,dimension(npar0)::fix0
	double precision,dimension(-2:(nz010+3))::zi010
	double precision,dimension(-2:(nz020+3))::zi020
	double precision,dimension(-2:(nz120+3))::zi120
	double precision,dimension(-2:(nz010-1))::the01
	double precision,dimension(-2:(nz120-1))::the12
	double precision,dimension(-2:(nz020-1))::the02
        double precision,dimension(no0,dimnva01)::ve010
	double precision,dimension(no0,dimnva02)::ve020
	double precision,dimension(no0,dimnva12)::ve120
	
	
        double precision::su01,ri01,su12,ri12,su02,ri02,gl01,gl02,gl12
	double precision,dimension(no0)::t00,t10,t20,t30
	integer,dimension(no0)::c0

	allocate(b(np0),bfix(npar0-np0),fix(npar0))
	b=b0
	bfix=bfix0
	fix=fix0
	allocate(zi01(-2:(nz01+3)),zi12(-2:(nz12+3)),zi02(-2:(nz02+3)))
	zi01=zi010
	zi02=zi020
	zi12=zi120

	
	nz01=nz010
	nz02=nz020
	nz12=nz120
	sizespline=nz01+nz02+nz12+6
	troncature=troncature0


	if(nva01.gt.0) then 
	  nva01nofix=nva01-sum(fix((sizespline+1):(nva01+sizespline)))
	else 
	  nva01nofix=0
	end if 

	if(nva02.gt.0) then 
          nva02nofix=nva02-sum(fix((nva01+sizespline+1):(nva02+nva01+sizespline)))
	else 
	   nva02nofix=0
	end if 

	if(nva12.gt.0) then 
	  nva12nofix=nva12-sum(fix((nva01+nva02+sizespline+1):npar0))
	else 
	  nva12nofix=0
	end if 

	nva0102=nva01nofix+nva02nofix
	nvamax=nva01nofix+nva02nofix+nva12nofix
	nvamax01=nvamax+nva01nofix
	nvamax02=nvamax01+nva02nofix
	nvamax12=nvamax02+nva12nofix


	if(nva01.gt.0) then 
		allocate(ve01(no0,nva01))
		allocate(ve01nofix(no0,nva01nofix))
		allocate(ve01square(no0,nva01nofix))
		allocate(tronc01(nva01nofix))
		allocate(tronc01square(nva01nofix))
	else 
		allocate(ve01(no0,1))
		allocate(ve01nofix(no0,1))
		allocate(ve01square(no0,1))
		ve01square=0
		ve01nofix=0
		allocate(tronc01(1))
		allocate(tronc01square(1))
	end if 
	
	if(nva02.gt.0) then 
		allocate(ve02(no0,nva02))
		allocate(ve02nofix(no0,nva02nofix))
		allocate(ve02square(no0,nva02nofix))
		allocate(tronc02(nva02nofix))
		allocate(tronc02square(nva02nofix))
	else 
		allocate(ve02(no0,1))
		allocate(ve02nofix(no0,1))
		allocate(ve02square(no0,1))
		ve02square=0
		ve02nofix=0
		allocate(tronc02(1))
		allocate(tronc02square(1))
	end if 

	if(nva12.gt.0) then 
		allocate(ve12(no0,nva12))
		allocate(ve12nofix(no0,nva12nofix))
		allocate(ve12square(no0,nva12nofix))
	else 
		allocate(ve12(no0,1))
		allocate(ve12nofix(no0,1))
		allocate(ve12square(no0,1))
	end if 


	ve01=ve010
	ve02=ve020
	ve12=ve120

	allocate(t0(no0),t1(no0),t2(no0),t3(no0),c(no0))
	c=c0
	t0=t00
	t1=t10
	t2=t20
	t3=t30

         
        ! we need to put bh at its original values if in posfix 


       l=0
       lfix=0
       w=0


	if(nva01.gt.0) then 
	do k=1,nva01
	   if(fix((sizespline+k)).eq.0) then 
		lfix=lfix+1
		ve01nofix(:,lfix)=ve01(:,k)
		ve01square(:,lfix)=ve01nofix(:,lfix)*ve01nofix(:,lfix)
	   end if 
	end do
	 
	end if 


	lfix=0
	if(nva02.gt.0) then 

	do k=1,nva02
	   if(fix((sizespline+nva01+k)).eq.0) then 
		lfix=lfix+1
		ve02nofix(:,lfix)=ve02(:,k)
		ve02square(:,lfix)=ve02nofix(:,lfix)*ve02nofix(:,lfix)
	   end if 
	end do
	end if 

	
	lfix=0

	if(nva12.gt.0) then 
	do k=1,nva12
	   if(fix((sizespline+nva01+nva02+k)).eq.0) then 
		lfix=lfix+1
		ve12nofix(:,lfix)=ve12(:,k)
		ve12square(:,lfix)=ve12nofix(:,lfix)*ve12nofix(:,lfix)
	   end if 
	end do
	end if 
	
	l=0
       lfix=0
       w=0

     do k=1,npar0 
         if(fix(k).eq.0) then
            l=l+1
            bh(k)=b(l)
	 end if 
         if(fix(k).eq.1) then
            w=w+1
            bh(k)=bfix(w)
         end if
      end do
   
	



         do i=1,nz01+2
            the01(i-3)=(bh(i))*(bh(i))
!       the01(i-3)=dexp(bh(i))
         end do
         do i=1,nz02+2
            j = nz01+2+i
            the02(i-3)=(bh(j))*(bh(j))
!       the12(i-3)=dexp(bh(j))
         end do
         do i=1,nz12+2
            j = nz02+2+nz01+2+i
            the12(i-3)=(bh(j))*(bh(j))
!       the02(i-3)=dexp(bh(j))
         end do


	
!---------- calcul des derivees premiere ------------------   

	res = 0
        do i=1,no0


                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0


                if(nva01.gt.0)then
                        do j=1,nva01
				vet01 =vet01 +&
                                bh(npar0-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
				vet02 =vet02 +&
                                bh(npar0-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
				vet12 =vet12 +&
                                bh(npar0-nva12+j)*dble(ve12(i,j))
                        end do
                endif

		
                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)

                res1 = 0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc01 = 0
				tronc02 =  0
                        	tronc01square= 0
                        	tronc02square=0
                        else 
				
				call susp(t0(i),the01,nz01,su01,ri01,zi01,gl01)
				call susp(t0(i),the02,nz02,su02,ri02,zi02,gl02)
                                tronc01=ve01nofix(i,:)*(gl01*vet01)
                        	tronc02=ve02nofix(i,:)*(gl02*vet02)
                        	tronc01square=ve01square(i,:)*(gl01*vet01)
                        	tronc02square=ve02square(i,:)*(gl02*vet02)
                        end if
                else
                        tronc01 = 0
                  	tronc02 = 0
                  	tronc01square=0
                  	tronc02square=0
                end if
		
                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
			call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
			call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)

			if(nva01nofix.gt.0) then 

			res1(1:nva01nofix)=&
			-ve01nofix(i,:)*gl01*vet01+tronc01
			res1((nvamax+1):nvamax01)=&
			-ve01square(i,:)*gl01*vet01+&
			tronc01square


			end if

			if(nva02nofix.gt.0) then 

			res1((nva01nofix+1):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02
			res1((nvamax01+1):nvamax02)=&
			-ve02square(i,:)*gl02*vet02+&
			tronc02square
			

			end if 

			if(nva12nofix.gt.0) then 

			res1((nva0102+1):nvamax)=0
			res1((nvamax02+1):nvamax12)=0

			end if 

		
			

                else
                if(c(i).eq.2)then ! cpi 0-->1
			call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
			call qgausssplinederivdiag(t1(i),t2(i),the01,&
			the02,the12,res2denum,res201num,&
			res202num,res212num,res20101num,&
			res20202num,res21212num,&
			vet01,vet02,vet12)
                        
			
			v=res2denum*(su12**vet12)

			if(nva01nofix.gt.0) then 

			u1=res201num*(su12**vet12)
      			res1(1:nva01nofix)=&
			ve01nofix(i,:)*u1/v
			res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01
			res1((nvamax+1):nvamax01)=&
			ve01square(i,:)*res20101num*(su12**vet12)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)/v
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)+&
			tronc01square

			end if 

			if(nva02nofix.gt.0) then 

			u2=-res202num*(su12**vet12)
			res1((nva01nofix+1):nva0102)=&
			ve02nofix(i,:)*u2/v
			res1((nva01nofix+1):nva0102)=&
			res1((nva01nofix+1):nva0102)+&
			tronc02

			res1((nvamax01+1):nvamax02)=&
			-ve02square(i,:)*(su12**vet12)*(-res20202num+res202num)
			res1((nvamax01+1):nvamax02)=&
			res1((nvamax01+1):nvamax02)/v
			res1((nvamax01+1):nvamax02)=&
			res1((nvamax01+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
			res1((nvamax01+1):nvamax02)=&
			res1((nvamax01+1):nvamax02)+tronc02square

			end if 

			if(nva12nofix.gt.0) then

			u3=res212num-gl12*vet12*res2denum
			u3=u3*(su12**vet12)
			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v

			res1((nvamax02+1):nvamax12)=&
			-2*gl12*vet12*res212num+res212num+&
			res21212num-gl12*vet12*(1-gl12*vet12)*res2denum
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)*(su12**vet12)*ve12square(i,:)
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)/v
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)-ve12square(i,:)*((u3/v)**2)

			end if 
			
			
                else  
                    if(c(i).eq.3)then ! obs 0-->1

			call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
			call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
			call susp(t1(i),the12,nz12,su12,ri12,zi12,gl12)


			if(nva01nofix.gt.0) then 

			res1(1:nva01nofix)=-ve01nofix(i,:)*gl01*vet01+&
			tronc01+ve01nofix(i,:)
			res1((nvamax+1):nvamax01)=&
			-ve01square(i,:)*gl01*vet01+&
			tronc01square

			end if 

			if(nva02nofix.gt.0) then 

			res1((nva01nofix+1):nva0102)=&
			-ve02nofix(i,:)*gl02*vet02+&
			tronc02
			res1((nvamax01+1):nvamax02)=&
			-ve02square(i,:)*gl02*vet02+&
			tronc02square

			end if 

			if(nva12nofix.gt.0) then 

			res1((nva0102+1):nvamax)=gl12*vet12
			res1((nvamax02+1):nvamax12)=gl12*vet12
			call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)-gl12*vet12
			res1((nva0102+1):nvamax)=&
			res1((nva0102+1):nvamax)*ve12nofix(i,:)
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)-gl12*vet12
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)*ve12square(i,:)
			

			end if 

			

                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			
			call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
			call qgausssplinederivdiag(t1(i),t2(i),the01,the02,the12,&
			res2denum,res201num,res202num,res212num,&
			res20101num,res20202num,&
			res21212num,vet01,vet02,vet12)

			
			v=res2denum*(su12**vet12)*ri12*vet12

			if(nva01nofix.gt.0) then 

			u1=res201num*(su12**vet12)*ri12*vet12
			res1(1:nva01nofix)=&
			ve01nofix(i,:)*u1/v
			res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01

			res1((nvamax+1):nvamax01)=&
			res20101num*(su12**vet12)*ri12*vet12*ve01square(i,:)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)/v
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
			res1((nvamax+1):nvamax01)=&
			res1((nvamax+1):nvamax01)+tronc01square
			

			end if 

			if(nva02nofix.gt.0) then

			u2=-res202num*(su12**vet12)*ri12*vet12
			res1((nva01nofix+1):(nva01nofix+nva02nofix))=&
			ve02nofix(i,:)*u2/v
			res1((nva01nofix+1):nva0102)=&
			res1((nva01nofix+1):nva0102)+&
			tronc02

			res1((nvamax01+1):nvamax02)=&
			(res20202num-res202num)*ri12*vet12*(su12**vet12)
			res1((nvamax01+1):nvamax02)=&
			res1((nvamax01+1):nvamax02)*ve02square(i,:)	
			res1((nvamax01+1):nvamax02)=&
			res1((nvamax01+1):nvamax02)/v
			res1((nvamax01+1):nvamax02)=&
			res1((nvamax01+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
			res1((nvamax01+1):nvamax02)=&
			res1((nvamax01+1):nvamax02)+tronc02square


			end if 

			if(nva12nofix.gt.0) then

			u3=res212num+&
			(1-gl12*vet12)*res2denum
			u3=u3*(su12**vet12)*ri12*vet12
			res1((nva0102+1):nvamax)=&
			ve12nofix(i,:)*u3/v

			res1((nvamax02+1):nvamax12)=&
			(1-gl12*vet12)*res212num+res21212num+&
			((1-gl12*vet12)**2)*res2denum+&
			(2-gl12*vet12)*res212num-gl12*vet12*res2denum
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)*(su12**vet12)*ri12*vet12
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)*ve12square(i,:)
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)/v
			res1((nvamax02+1):nvamax12)=&
			res1((nvamax02+1):nvamax12)-ve12square(i,:)*((u3/v)**2)			


			end if 

			
                        
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
				
			         call susp(t1(i),the01,nz01,su01,ri01,zi01,gl01)
			         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
			         call susp(t1(i),the12,nz12,su12,ri12,zi12,gl12)

				if(nva01nofix.gt.0) then 

				res1(1:nva01nofix)=&
				-ve01nofix(i,:)*gl01*vet01+&
				ve01nofix(i,:)+&
				tronc01
				res1((nvamax+1):nvamax01)=&
				-ve01square(i,:)*gl01*vet01+&
				tronc01square

				end if 

				if(nva02nofix.gt.0) then 

				res1((nva01nofix+1):nva0102)=&
				-ve02nofix(i,:)*gl02*vet02+tronc02
				res1((nvamax01+1):nvamax02)=&
				-ve02square(i,:)*gl02*vet02+&
				tronc02square
				

				end if 

				if(nva12nofix.gt.0) then

				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*gl12*vet12+&
				ve12nofix(i,:)
				res1((nvamax02+1):nvamax12)=&
				gl12*vet12
				call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                res1((nva0102+1):nvamax)=&
				res1((nva0102+1):nvamax)-ve12nofix(i,:)*gl12*vet12
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)-gl12*vet12
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)*ve12square(i,:)

				end if 
				
                         else
                            if(c(i).eq.6)then ! vivant ???
				
				call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
				call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
			        call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
				call qgausssplinederivdiag(t1(i),t3(i),the01,the02,the12,&
				res2denum,res201num,res202num,res212num,res20101num,&
				res20202num,res21212num,vet01,vet02,vet12)

				
				v=(su12**vet12)*res2denum+&
				(su01**vet01)*(su02**vet02)
			
				if(nva01nofix.gt.0) then 

				u1=(-gl01*vet01)*(su01**vet01)*(su02**vet02)+&
				(su12**vet12)*res201num
			        res1(1:nva01nofix)=&
				ve01nofix(i,:)*u1/v
				res1(1:nva01nofix)=&
				res1(1:nva01nofix)+tronc01

				res1((nvamax+1):nvamax01)=&
				res20101num*(su12**vet12)
				resint=(su01**vet01)*(su02**vet02)
				resint=resint*gl01*vet01*(1-gl01*vet01)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-resint
				
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)*ve01square(i,:)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)/v

				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)+tronc01square

				end if 


				if(nva02nofix.gt.0) then 

				u2=-gl02*vet02*(su01**vet01)*(su02**vet02)
				u2=u2-(su12**vet12)*res202num
				res1((nva01nofix+1):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1):nva0102)=&
				res1((nva01nofix+1):nva0102)+&
				tronc02

				res1((nvamax01+1):nvamax02)=&
				(su12**vet12)*(res20202num-res202num)-&
				gl02*vet02*(su01**vet01)*(su02**vet02)*(1-gl02*vet02)
				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)*ve02square(i,:)
				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)/v

				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)+tronc02square

				end if 

				if(nva12nofix.gt.0) then 
				
				u3=-gl12*vet12*(su12**vet12)*res2denum+&
				(su12**vet12)*res212num
				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v

				res1((nvamax02+1):nvamax12)=&
				res21212num+(1-gl12*vet12)*res212num-&
				gl12*vet12*(1-gl12*vet12)*res2denum-&
				gl12*vet12*res212num
				resint=(su12**vet12)
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)*resint*ve12square(i,:)
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)/v
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)-ve12square(i,:)*((u3/v)**2)
			
				end if 
                        	

                            else ! passage 0-->2  
				
				call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
				call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
			        call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
        			call qgausssplinederivdiag(t1(i),t3(i),the01,the02,&
        			the12,res2denum,res201num,res202num,res212num,&
        			res20101num,res20202num,res21212num,&
        			vet01,vet02,vet12)

		
				v=(su12**vet12)*ri12*vet12*res2denum+&
				(su01**vet01)*(su02**vet02)*ri02*vet02


				if(nva01nofix.gt.0) then 

				u1=-gl01*vet01*(su01**vet01)
				u1=u1*(su02**vet02)*ri02*vet02
				u1=u1+(su12**vet12)*ri12*vet12*res201num

				res1(1:nva01nofix)=&
				ve01nofix(i,:)*u1/v
				res1(1:nva01nofix)=res1(1:nva01nofix)+tronc01

				res1((nvamax+1):nvamax01)=&
				res20101num*(su12**vet12)*ri12*vet12
				resint=(su01**vet01)*(su02**vet02)*gl01*vet01
				resint=resint*(1-gl01*vet01)*ri02*vet02
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-resint
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)*ve01square(i,:)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)/v
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)-ve01square(i,:)*((u1/v)**2)
				res1((nvamax+1):nvamax01)=&
				res1((nvamax+1):nvamax01)+tronc01square

				end if 


				if(nva02nofix.gt.0) then 

				u2=-(su12**vet12)*ri12*vet12*res202num+&
				(1-gl02*vet02)*(su01**vet01)*(su02**vet02)*ri02*vet02

				res1((nva01nofix+1):nva0102)=&
				ve02nofix(i,:)*u2/v
				res1((nva01nofix+1):nva0102)=&
				res1((nva01nofix+1):nva0102)+&
				tronc02

				res1((nvamax01+1):nvamax02)=&
				-gl02*vet02*(su01**vet01)*(su02**vet02)*ri02*vet02
				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)+&
				(su12**vet12)*ri12*vet12*(res20202num-res202num)+&
				((1-gl02*vet02)**2)*ri02*vet02*(su01**vet01)*(su02**vet02)


				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)*ve02square(i,:)
				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)/v
				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)-ve02square(i,:)*((u2/v)**2)
				res1((nvamax01+1):nvamax02)=&
				res1((nvamax01+1):nvamax02)+tronc02square

				end if 

				if(nva12nofix.gt.0) then 

				u3=res212num+(1-gl12*vet12)*res2denum
				u3=u3*(su12**vet12)*ri12*vet12

				res1((nva0102+1):nvamax)=&
				ve12nofix(i,:)*u3/v

				res1((nvamax02+1):nvamax12)=&
				-gl12*vet12*res2denum+((1-gl12*vet12)**2)*res2denum+&
				2*(1-gl12*vet12)*res212num+res212num+&
				res21212num

				resint=ri12*vet12*(su12**vet12)
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)*resint*ve12square(i,:)
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)/v
				
				res1((nvamax02+1):nvamax12)=&
				res1((nvamax02+1):nvamax12)-ve12square(i,:)*((u3/v)**2)

				end if 


				
                            endif
                         endif                        
                      endif
                   endif  
                endif   
                endif   

                res = res + res1 

        end do   

        likelihood_deriv = res




123     continue 
	 
	deallocate(b,bfix,fix,zi01,zi02,zi12,ve01,ve02,ve12,ve01nofix,&
	ve02nofix,ve12nofix,tronc01,tronc02,t0,t1,t2,t3,c,&
	ve01square,ve02square,ve12square,&
	tronc01square,tronc02square)     

    end subroutine derivasplinediag




