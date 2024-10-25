 



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
!               res1=dlog(res2)-gl12*vet12 (autre ecriture)
		                      res1=dlog(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
                    call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                    call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                    call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                    res2=-gl01*vet01+dlog(ri01*vet01) -&
                    gl12*vet12-gl02*vet02
                    call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                    res1=res2 +gl12*vet12
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                          call qgaussPL(t1(i),t2(i),the01,the12,the02,&
                        	res2,vet01,vet12,vet02)
                          res1=dlog(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
                         call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                         call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                         res2=-gl01*vet01 + dlog(ri01*vet01) -&
                         gl12*vet12 + dlog(ri12*vet12) -gl02*vet02
                         call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                         res1=res2+gl12*vet12

                         else
                            if(c(i).eq.6)then ! vivant ???
                          		call qgaussPL(t1(i),t2(i),the01,the12,the02,&
                              res2,vet01,vet12,vet02)
                              call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                              call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                              call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                              res1=dlog(res2*(su12**vet12) +&
                            (su01**vet01)*(su02**vet02))  

                            else ! passage 0-->2  

                                call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                                call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                        	    	call qgaussPL(t1(i),t3(i),the01,the12,&
                        		the02,res2,vet01,vet12,vet02)
                                res1=dlog(res2*ri12*vet12*(su12**vet12) +&
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
!               res1=dlog(res2)-gl12*vet12 (autre ecriture)
		                      res1=dlog(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
                    call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                    call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                    call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                    res2=-gl01*vet01+dlog(ri01*vet01) -&
                    gl12*vet12-gl02*vet02
                    call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                    res1=res2 +gl12*vet12
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                          call qgaussPL15(t1(i),t2(i),the01,the12,the02,&
                        	res2,vet01,vet12,vet02)
                          res1=dlog(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
                         call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                         call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                         res2=-gl01*vet01 + dlog(ri01*vet01) -&
                         gl12*vet12 + dlog(ri12*vet12) -gl02*vet02
                         call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                         res1=res2+gl12*vet12

                         else
                            if(c(i).eq.6)then ! vivant ???
                          		call qgaussPL15(t1(i),t2(i),the01,the12,the02,&
                              res2,vet01,vet12,vet02)
                              call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                              call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                              call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                              res1=dlog(res2*(su12**vet12) +&
                            (su01**vet01)*(su02**vet02))  

                            else ! passage 0-->2  

                                call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                                call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                        	    	call qgaussPL15(t1(i),t3(i),the01,the12,&
                        		the02,res2,vet01,vet12,vet02)
                                res1=dlog(res2*ri12*vet12*(su12**vet12) +&
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
!               res1=dlog(res2)-gl12*vet12 (autre ecriture)
		                      res1=dlog(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
                    call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                    call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                    call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                    res2=-gl01*vet01+dlog(ri01*vet01) -&
                    gl12*vet12-gl02*vet02
                    call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                    res1=res2 +gl12*vet12
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                          call qgaussPL21(t1(i),t2(i),the01,the12,the02,&
                        	res2,vet01,vet12,vet02)
                          res1=dlog(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
                         call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                         call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                         res2=-gl01*vet01 + dlog(ri01*vet01) -&
                         gl12*vet12 + dlog(ri12*vet12) -gl02*vet02
                         call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                         res1=res2+gl12*vet12

                         else
                            if(c(i).eq.6)then ! vivant ???
                          		call qgaussPL21(t1(i),t2(i),the01,the12,the02,&
                              res2,vet01,vet12,vet02)
                              call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                              call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                              call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                              res1=dlog(res2*(su12**vet12) +&
                            (su01**vet01)*(su02**vet02))  

                            else ! passage 0-->2  

                                call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                                call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                        	    	call qgaussPL21(t1(i),t3(i),the01,the12,&
                        		the02,res2,vet01,vet12,vet02)
                                res1=dlog(res2*ri12*vet12*(su12**vet12) +&
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
!               res1=dlog(res2)-gl12*vet12 (autre ecriture)
		                      res1=dlog(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
                    call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                    call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                    call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                    res2=-gl01*vet01+dlog(ri01*vet01) -&
                    gl12*vet12-gl02*vet02
                    call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                    res1=res2 +gl12*vet12
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                          call qgaussPL31(t1(i),t2(i),the01,the12,the02,&
                        	res2,vet01,vet12,vet02)
                          res1=dlog(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
                         call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                         call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                         res2=-gl01*vet01 + dlog(ri01*vet01) -&
                         gl12*vet12 + dlog(ri12*vet12) -gl02*vet02
                         call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                         res1=res2+gl12*vet12

                         else
                            if(c(i).eq.6)then ! vivant ???
                          		call qgaussPL31(t1(i),t2(i),the01,the12,the02,&
                              res2,vet01,vet12,vet02)
                              call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                              call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                              call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                              res1=dlog(res2*(su12**vet12) +&
                            (su01**vet01)*(su02**vet02))  

                            else ! passage 0-->2  

                                call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                                call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                        	    	call qgaussPL31(t1(i),t3(i),the01,the12,&
                        		the02,res2,vet01,vet12,vet02)
                                res1=dlog(res2*ri12*vet12*(su12**vet12) +&
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
!               res1=dlog(res2)-gl12*vet12 (autre ecriture)
		                      res1=dlog(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
                    call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                    call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                    call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                    res2=-gl01*vet01+dlog(ri01*vet01) -&
                    gl12*vet12-gl02*vet02
                    call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                    res1=res2 +gl12*vet12
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                          call qgaussPL41(t1(i),t2(i),the01,the12,the02,&
                        	res2,vet01,vet12,vet02)
                          res1=dlog(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
                         call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                         call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                         res2=-gl01*vet01 + dlog(ri01*vet01) -&
                         gl12*vet12 + dlog(ri12*vet12) -gl02*vet02
                         call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                         res1=res2+gl12*vet12

                         else
                            if(c(i).eq.6)then ! vivant ???
                          		call qgaussPL41(t1(i),t2(i),the01,the12,the02,&
                              res2,vet01,vet12,vet02)
                              call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                              call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                              call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                              res1=dlog(res2*(su12**vet12) +&
                            (su01**vet01)*(su02**vet02))  

                            else ! passage 0-->2  

                                call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                                call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                        	    	call qgaussPL41(t1(i),t3(i),the01,the12,&
                        		the02,res2,vet01,vet12,vet02)
                                res1=dlog(res2*ri12*vet12*(su12**vet12) +&
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
!               res1=dlog(res2)-gl12*vet12 (autre ecriture)
		                      res1=dlog(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
                    call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                    call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                    call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                    res2=-gl01*vet01+dlog(ri01*vet01) -&
                    gl12*vet12-gl02*vet02
                    call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                    res1=res2 +gl12*vet12
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                          call qgaussPL51(t1(i),t2(i),the01,the12,the02,&
                        	res2,vet01,vet12,vet02)
                          res1=dlog(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
                         call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                         call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                         res2=-gl01*vet01 + dlog(ri01*vet01) -&
                         gl12*vet12 + dlog(ri12*vet12) -gl02*vet02
                         call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                         res1=res2+gl12*vet12

                         else
                            if(c(i).eq.6)then ! vivant ???
                          		call qgaussPL51(t1(i),t2(i),the01,the12,the02,&
                              res2,vet01,vet12,vet02)
                              call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                              call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                              call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                              res1=dlog(res2*(su12**vet12) +&
                            (su01**vet01)*(su02**vet02))  

                            else ! passage 0-->2  

                                call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                                call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                        	    	call qgaussPL51(t1(i),t3(i),the01,the12,&
                        		the02,res2,vet01,vet12,vet02)
                                res1=dlog(res2*ri12*vet12*(su12**vet12) +&
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
!               res1=dlog(res2)-gl12*vet12 (autre ecriture)
		                      res1=dlog(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
                    call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                    call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                    call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                    res2=-gl01*vet01+dlog(ri01*vet01) -&
                    gl12*vet12-gl02*vet02
                    call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                    res1=res2 +gl12*vet12
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
                          call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                          call qgaussPL61(t1(i),t2(i),the01,the12,the02,&
                        	res2,vet01,vet12,vet02)
                          res1=dlog(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				
                         call susp(t2(i),the01,nz01,su01,ri01,zi01,gl01)
                         call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                         call susp(t1(i),the02,nz02,su02,ri02,zi02,gl02)
                         res2=-gl01*vet01 + dlog(ri01*vet01) -&
                         gl12*vet12 + dlog(ri12*vet12) -gl02*vet02
                         call susp(t2(i),the12,nz12,su12,ri12,zi12,gl12)
                         res1=res2+gl12*vet12

                         else
                            if(c(i).eq.6)then ! vivant ???
                          		call qgaussPL61(t1(i),t2(i),the01,the12,the02,&
                              res2,vet01,vet12,vet02)
                              call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                              call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                              call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                              res1=dlog(res2*(su12**vet12) +&
                            (su01**vet01)*(su02**vet02))  

                            else ! passage 0-->2  

                                call susp(t3(i),the12,nz12,su12,ri12,zi12,gl12)
                                call susp(t3(i),the02,nz02,su02,ri02,zi02,gl02)
                                call susp(t3(i),the01,nz01,su01,ri01,zi01,gl01)
                        	    	call qgaussPL61(t1(i),t3(i),the01,the12,&
                        		the02,res2,vet01,vet12,vet02)
                                res1=dlog(res2*ri12*vet12*(su12**vet12) +&
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
      nva12,nva02,t00,t10,t20,t30,troncature0,gausspoint0,weib0,likelihood_res)

	use commun
        implicit none
         
        double precision::res,res1,res2,tronc, &
        vet01,vet12,vet02
	double precision, intent(inout)::likelihood_res
        integer::np0,i,j,l,w,k,npar0,nva01,nva12,nva02,no0, &
	weib0,troncature0,dimnva01,dimnva02,dimnva12,gausspoint0

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
       weib=weib0

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
 
	

	if(weib.eq.1)then
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
	else 
	 do i=1,2
            the01(i)=dexp(bh(i))
         end do
         do i=1,2
            j = 2+i
            the02(i)=dexp(bh(j))
         end do
         do i=1,2
            j = 4+i
            the12(i)=dexp(bh(j))
         end do
	endif
	
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
                        res1=dlog(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t1(i),the02,ri02,gl02,su02)
                        call fonct(t1(i),the12,ri12,gl12,su12)
                        res1 = -(gl01*vet01)-(gl02*vet02)+&
                        dlog(ri01*vet01)+(gl12*vet12)
                        call fonct(t3(i),the12,ri12,gl12,su12)
                        res1 = res1 -(gl12*vet12)
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			call fonct(t3(i),the12,ri12,gl12,su12)
                        call  qgaussPLweib(t1(i),t2(i),the01,the02,the12,&
                        res2,vet01,vet02,vet12)
                        res1=dlog(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				call fonct(t1(i),the01,ri01,gl01,su01)
                                call fonct(t1(i),the02,ri02,gl02,su02)
                                call fonct(t1(i),the12,ri12,gl12,su12)
                                res1 = -(gl01*vet01)-(gl02*vet02)+&
                                dlog(ri01*vet01)+(gl12*vet12)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                res1 = res1 -(gl12*vet12) + dlog(ri12*vet12)
                         else
                            if(c(i).eq.6)then ! vivant ???
				 call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPLweib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12))+&
                                ((su01**vet01)*(su02**vet02))
                                res1 = dlog(res1)

                            else ! passage 0-->2  

				call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPLweib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12)*ri12*vet12)+&
                                ((su01**vet01)*(su02**vet02)*ri02*vet02)
                                res1 = dlog(res1)
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
                        res1=dlog(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t1(i),the02,ri02,gl02,su02)
                        call fonct(t1(i),the12,ri12,gl12,su12)
                        res1 = -(gl01*vet01)-(gl02*vet02)+&
                        dlog(ri01*vet01)+(gl12*vet12)
                        call fonct(t3(i),the12,ri12,gl12,su12)
                        res1 = res1 -(gl12*vet12)
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			call fonct(t3(i),the12,ri12,gl12,su12)
                        call  qgaussPL15weib(t1(i),t2(i),the01,the02,the12,&
                        res2,vet01,vet02,vet12)
                        res1=dlog(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				call fonct(t1(i),the01,ri01,gl01,su01)
                                call fonct(t1(i),the02,ri02,gl02,su02)
                                call fonct(t1(i),the12,ri12,gl12,su12)
                                res1 = -(gl01*vet01)-(gl02*vet02)+&
                                dlog(ri01*vet01)+(gl12*vet12)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                res1 = res1 -(gl12*vet12) + dlog(ri12*vet12)
                         else
                            if(c(i).eq.6)then ! vivant ???
				 call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL15weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12))+&
                                ((su01**vet01)*(su02**vet02))
                                res1 = dlog(res1)

                            else ! passage 0-->2  

				call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL15weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12)*ri12*vet12)+&
                                ((su01**vet01)*(su02**vet02)*ri02*vet02)
                                res1 = dlog(res1)
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
                        res1=dlog(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t1(i),the02,ri02,gl02,su02)
                        call fonct(t1(i),the12,ri12,gl12,su12)
                        res1 = -(gl01*vet01)-(gl02*vet02)+&
                        dlog(ri01*vet01)+(gl12*vet12)
                        call fonct(t3(i),the12,ri12,gl12,su12)
                        res1 = res1 -(gl12*vet12)
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			call fonct(t3(i),the12,ri12,gl12,su12)
                        call  qgaussPL21weib(t1(i),t2(i),the01,the02,the12,&
                        res2,vet01,vet02,vet12)
                        res1=dlog(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				call fonct(t1(i),the01,ri01,gl01,su01)
                                call fonct(t1(i),the02,ri02,gl02,su02)
                                call fonct(t1(i),the12,ri12,gl12,su12)
                                res1 = -(gl01*vet01)-(gl02*vet02)+&
                                dlog(ri01*vet01)+(gl12*vet12)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                res1 = res1 -(gl12*vet12) + dlog(ri12*vet12)
                         else
                            if(c(i).eq.6)then ! vivant ???
				 call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL21weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12))+&
                                ((su01**vet01)*(su02**vet02))
                                res1 = dlog(res1)

                            else ! passage 0-->2  

				call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL21weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12)*ri12*vet12)+&
                                ((su01**vet01)*(su02**vet02)*ri02*vet02)
                                res1 = dlog(res1)
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
                        res1=dlog(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t1(i),the02,ri02,gl02,su02)
                        call fonct(t1(i),the12,ri12,gl12,su12)
                        res1 = -(gl01*vet01)-(gl02*vet02)+&
                        dlog(ri01*vet01)+(gl12*vet12)
                        call fonct(t3(i),the12,ri12,gl12,su12)
                        res1 = res1 -(gl12*vet12)
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			call fonct(t3(i),the12,ri12,gl12,su12)
                        call  qgaussPL31weib(t1(i),t2(i),the01,the02,the12,&
                        res2,vet01,vet02,vet12)
                        res1=dlog(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				call fonct(t1(i),the01,ri01,gl01,su01)
                                call fonct(t1(i),the02,ri02,gl02,su02)
                                call fonct(t1(i),the12,ri12,gl12,su12)
                                res1 = -(gl01*vet01)-(gl02*vet02)+&
                                dlog(ri01*vet01)+(gl12*vet12)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                res1 = res1 -(gl12*vet12) + dlog(ri12*vet12)
                         else
                            if(c(i).eq.6)then ! vivant ???
				 call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL31weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12))+&
                                ((su01**vet01)*(su02**vet02))
                                res1 = dlog(res1)

                            else ! passage 0-->2  

				call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL31weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12)*ri12*vet12)+&
                                ((su01**vet01)*(su02**vet02)*ri02*vet02)
                                res1 = dlog(res1)
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
                        res1=dlog(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t1(i),the02,ri02,gl02,su02)
                        call fonct(t1(i),the12,ri12,gl12,su12)
                        res1 = -(gl01*vet01)-(gl02*vet02)+&
                        dlog(ri01*vet01)+(gl12*vet12)
                        call fonct(t3(i),the12,ri12,gl12,su12)
                        res1 = res1 -(gl12*vet12)
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			call fonct(t3(i),the12,ri12,gl12,su12)
                        call  qgaussPL41weib(t1(i),t2(i),the01,the02,the12,&
                        res2,vet01,vet02,vet12)
                        res1=dlog(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				call fonct(t1(i),the01,ri01,gl01,su01)
                                call fonct(t1(i),the02,ri02,gl02,su02)
                                call fonct(t1(i),the12,ri12,gl12,su12)
                                res1 = -(gl01*vet01)-(gl02*vet02)+&
                                dlog(ri01*vet01)+(gl12*vet12)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                res1 = res1 -(gl12*vet12) + dlog(ri12*vet12)
                         else
                            if(c(i).eq.6)then ! vivant ???
				 call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL41weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12))+&
                                ((su01**vet01)*(su02**vet02))
                                res1 = dlog(res1)

                            else ! passage 0-->2  

				call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL41weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12)*ri12*vet12)+&
                                ((su01**vet01)*(su02**vet02)*ri02*vet02)
                                res1 = dlog(res1)
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
                        res1=dlog(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t1(i),the02,ri02,gl02,su02)
                        call fonct(t1(i),the12,ri12,gl12,su12)
                        res1 = -(gl01*vet01)-(gl02*vet02)+&
                        dlog(ri01*vet01)+(gl12*vet12)
                        call fonct(t3(i),the12,ri12,gl12,su12)
                        res1 = res1 -(gl12*vet12)
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			call fonct(t3(i),the12,ri12,gl12,su12)
                        call  qgaussPL51weib(t1(i),t2(i),the01,the02,the12,&
                        res2,vet01,vet02,vet12)
                        res1=dlog(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				call fonct(t1(i),the01,ri01,gl01,su01)
                                call fonct(t1(i),the02,ri02,gl02,su02)
                                call fonct(t1(i),the12,ri12,gl12,su12)
                                res1 = -(gl01*vet01)-(gl02*vet02)+&
                                dlog(ri01*vet01)+(gl12*vet12)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                res1 = res1 -(gl12*vet12) + dlog(ri12*vet12)
                         else
                            if(c(i).eq.6)then ! vivant ???
				 call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL51weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12))+&
                                ((su01**vet01)*(su02**vet02))
                                res1 = dlog(res1)

                            else ! passage 0-->2  

				call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL51weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12)*ri12*vet12)+&
                                ((su01**vet01)*(su02**vet02)*ri02*vet02)
                                res1 = dlog(res1)
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
                        res1=dlog(res2*(su12**vet12))

                else  
                    if(c(i).eq.3)then ! obs 0-->1
			call fonct(t1(i),the01,ri01,gl01,su01)
                        call fonct(t1(i),the02,ri02,gl02,su02)
                        call fonct(t1(i),the12,ri12,gl12,su12)
                        res1 = -(gl01*vet01)-(gl02*vet02)+&
                        dlog(ri01*vet01)+(gl12*vet12)
                        call fonct(t3(i),the12,ri12,gl12,su12)
                        res1 = res1 -(gl12*vet12)
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2
			call fonct(t3(i),the12,ri12,gl12,su12)
                        call  qgaussPL61weib(t1(i),t2(i),the01,the02,the12,&
                        res2,vet01,vet02,vet12)
                        res1=dlog(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
				call fonct(t1(i),the01,ri01,gl01,su01)
                                call fonct(t1(i),the02,ri02,gl02,su02)
                                call fonct(t1(i),the12,ri12,gl12,su12)
                                res1 = -(gl01*vet01)-(gl02*vet02)+&
                                dlog(ri01*vet01)+(gl12*vet12)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                res1 = res1 -(gl12*vet12) + dlog(ri12*vet12)
                         else
                            if(c(i).eq.6)then ! vivant ???
				 call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL61weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12))+&
                                ((su01**vet01)*(su02**vet02))
                                res1 = dlog(res1)

                            else ! passage 0-->2  

				call fonct(t3(i),the01,ri01,gl01,su01)
                                call fonct(t3(i),the02,ri02,gl02,su02)
                                call fonct(t3(i),the12,ri12,gl12,su12)
                                call  qgaussPL61weib(t1(i),t3(i),the01,the02,&
                                the12,res2,vet01,vet02,vet12)
                                res1 = (res2*(su12**vet12)*ri12*vet12)+&
                                ((su01**vet01)*(su02**vet02)*ri02*vet02)
                                res1 = dlog(res1)
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
         double precision::a,b,dx,xm,xr,res,resabs,resk,reskh,resasc,resg,v01,v02,v12,&
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

    	resg = fc*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        resk = fc*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        resabs = dabs(resk)    ! init res absolue Kronrod   
         
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
	       	resg = resg+wg(j)*(f1+f2)
               	resk = resk + wgk(jtw)*(f1+f2)
               	resabs = resabs + wgk(jtw)*(dabs(f1)+dabs(f2))
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
               resabs = resabs + wgk(jtwm1)*(dabs(f1)+dabs(f2))
              end do
	    

	reskh = resk*0.5d+00         ! res Kronrod / 2
        resasc = wgk(8)*dabs(fc-reskh)  ! init  ! w8 * abs(fc - res Kronrod / 2 )
         do j=1,7
        	resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))  ! ite
    	 end do
    
    	res = resk*xr
	endif
    
          end subroutine qgaussPL15weib

!=============================================================================================  
!==== QGAUS15 out a 15 point Gauss-Kronrod quadrature rule for splines   =====================
!=============================================================================================  

      subroutine qgaussPL15(a,b,the01,the12,the02,res,v1,v2,v3)

 
	use commun,only:zi01,zi12,zi02,nz01,nz12,nz02

         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resabs,resk,reskh,resasc,resg,v1,v2,v3,&
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

    	resg = fc*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        resk = fc*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        resabs = dabs(resk)    ! init res absolue Kronrod   
         

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
	       resg = resg+wg(j)*(f1+f2)
               resk = resk + wgk(jtw)*(f1+f2)
               resabs = resabs + wgk(jtw)*(dabs(f1)+dabs(f2))
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
               resabs = resabs + wgk(jtwm1)*(dabs(f1)+dabs(f2))
         end do
     ! pr calcul erreur 
         reskh = resk*0.5d+00         ! res Kronrod / 2
         resasc = wgk(8)*dabs(fc-reskh)  ! init  ! w8 * abs(fc - res Kronrod / 2 )
         do j=1,7
        	resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))  ! ite
    	 end do
    
    res = resk*xr
    !resabs = resabs*dabs(xr)
    !resasc = resasc*dabs(xr)
    !abserr = dabs((resk-resg)*xr)   ! estimation erreur Kronrod pra Gauss
    !if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00) abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00) ! si erreurs non nulles
    !if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1((epmach*0.5d+02)*resabs,abserr)

         end subroutine qgaussPL15

!=============================================================================================  
!=====QGAUS21 out a 21 point Gauss-Kronrod quadrature rule for weib ==========================
!=============================================================================================  

subroutine qgaussPL21weib(a,b,the01,the02,the12,res,v01,v02,v12)
         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resabs,resk,reskh,resasc,resg,v01,v02,v12,&
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

    	resg = 0        ! init 0
        resk = fc*wgk(11)       ! init res Kronrod   ! fc * 8e poids Kronrod
        resabs = dabs(resk)    ! init res absolue Kronrod   
         
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
	       	resg = resg+wg(j)*(f1+f2)
                resk = resk + wgk(jtw)*(f1+f2)
                resabs = resabs + wgk(jtw)*(dabs(f1)+dabs(f2))
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
               resabs = resabs + wgk(jtwm1)*(dabs(f1)+dabs(f2))
              end do
	    

	 reskh = resk*0.5d+00         ! res Kronrod / 2
         resasc = wgk(11)*dabs(fc-reskh)  ! init  ! w8 * abs(fc - res Kronrod / 2 )
         do j=1,10
        	resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))  ! ite
    	 end do
    
    	res = resk*xr
	endif
    
          end subroutine qgaussPL21weib

!=============================================================================================  
!=====QGAUS21 out a 21 point Gauss-Kronrod quadrature rule for splines =======================
!=============================================================================================  


      subroutine qgaussPL21(a,b,the01,the12,the02,res,v1,v2,v3)

	use commun,only:zi01,zi12,zi02,nz01,nz12,nz02

         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resabs,resk,reskh,resasc,resg,v1,v2,v3,&
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

    	resg = 0        ! init 0
        resk = fc*wgk(11)       ! init res Kronrod   ! fc * 8e poids Kronrod
        resabs = dabs(resk)    ! init res absolue Kronrod   
         

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
	       resg = resg+wg(j)*(f1+f2)
               resk = resk + wgk(jtw)*(f1+f2)
               resabs = resabs + wgk(jtw)*(dabs(f1)+dabs(f2))
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
               resabs = resabs + wgk(jtwm1)*(dabs(f1)+dabs(f2))
         end do
     ! pr calcul erreur 
         reskh = resk*0.5d+00         ! res Kronrod / 2
         resasc = wgk(11)*dabs(fc-reskh)  ! init  ! w8 * abs(fc - res Kronrod / 2 )
         do j=1,10
        	resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))  ! ite
    	 end do
    
    res = resk*xr
    !resabs = resabs*dabs(xr)
    !resasc = resasc*dabs(xr)
    !abserr = dabs((resk-resg)*xr)   ! estimation erreur Kronrod pra Gauss
    !if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00) abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00) ! si erreurs non nulles
    !if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1((epmach*0.5d+02)*resabs,abserr)

         end subroutine qgaussPL21

!=============================================================================================  
!=====QGAUS31 out a 31 point Gauss-Kronrod quadrature rule for weib ==========================
!=============================================================================================  


subroutine qgaussPL31weib(a,b,the01,the02,the12,res,v01,v02,v12)
         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resabs,resk,reskh,resasc,resg,v01,v02,v12,&
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

    	resg = wg(8)*fc        ! init at wg(8)
        resk = fc*wgk(16)       ! init res Kronrod   ! fc * 8e poids Kronrod
        resabs = dabs(resk)    ! init res absolue Kronrod   
         
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
	       	resg = resg+wg(j)*(f1+f2)
                resk = resk + wgk(jtw)*(f1+f2)
                resabs = resabs + wgk(jtw)*(dabs(f1)+dabs(f2))
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
               resabs = resabs + wgk(jtwm1)*(dabs(f1)+dabs(f2))
              end do
	    

	 reskh = resk*0.5d+00         ! res Kronrod / 2
         resasc = wgk(16)*dabs(fc-reskh)  ! init  ! w8 * abs(fc - res Kronrod / 2 )
         do j=1,15
        	resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))  ! ite
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
         double precision::a,b,dx,xm,xr,res,resabs,resk,reskh,resasc,resg,v1,v2,v3,&
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

    	resg = wg(8)*fc        ! init at wg(8)
        resk = fc*wgk(16)       ! init res Kronrod   ! fc * 8e poids Kronrod
        resabs = dabs(resk)    ! init res absolue Kronrod   
         

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
	       resg = resg+wg(j)*(f1+f2)
               resk = resk + wgk(jtw)*(f1+f2)
               resabs = resabs + wgk(jtw)*(dabs(f1)+dabs(f2))
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
               resabs = resabs + wgk(jtwm1)*(dabs(f1)+dabs(f2))
         end do
     ! pr calcul erreur 
         reskh = resk*0.5d+00         ! res Kronrod / 2
         resasc = wgk(16)*dabs(fc-reskh)  ! init  ! w8 * abs(fc - res Kronrod / 2 )
         do j=1,15
        	resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))  ! ite
    	 end do
    
    res = resk*xr
    !resabs = resabs*dabs(xr)
    !resasc = resasc*dabs(xr)
    !abserr = dabs((resk-resg)*xr)   ! estimation erreur Kronrod pra Gauss
    !if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00) abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00) ! si erreurs non nulles
    !if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1((epmach*0.5d+02)*resabs,abserr)

         end subroutine qgaussPL31

!=============================================================================================  
!=====QGAUS41 out a 41 point Gauss-Kronrod quadrature rule for weib ==========================
!=============================================================================================  

subroutine qgaussPL41weib(a,b,the01,the02,the12,res,v01,v02,v12)
         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resabs,resk,reskh,resasc,resg,v01,v02,v12,&
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

    	resg = 0        ! init at wg(8)
        resk = fc*wgk(21)       ! init res Kronrod   ! fc * 8e poids Kronrod
        resabs = dabs(resk)    ! init res absolue Kronrod   
         
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
	       	resg = resg+wg(j)*(f1+f2)
                resk = resk + wgk(jtw)*(f1+f2)
                resabs = resabs + wgk(jtw)*(dabs(f1)+dabs(f2))
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
               resabs = resabs + wgk(jtwm1)*(dabs(f1)+dabs(f2))
              end do
	    

	  reskh = resk*0.5d+00         ! res Kronrod / 2
         resasc = wgk(21)*dabs(fc-reskh)  ! init  ! w8 * abs(fc - res Kronrod / 2 )
         do j=1,20
        	resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))  ! ite
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
         double precision::a,b,dx,xm,xr,res,resabs,resk,reskh,resasc,resg,v1,v2,v3,&
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

    	resg = 0        ! init at 0
        resk = fc*wgk(21)       ! init res Kronrod   ! fc * 8e poids Kronrod
        resabs = dabs(resk)    ! init res absolue Kronrod   
         

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
	       resg = resg+wg(j)*(f1+f2)
               resk = resk + wgk(jtw)*(f1+f2)
               resabs = resabs + wgk(jtw)*(dabs(f1)+dabs(f2))
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
               resabs = resabs + wgk(jtwm1)*(dabs(f1)+dabs(f2))
         end do
     ! pr calcul erreur 
         reskh = resk*0.5d+00         ! res Kronrod / 2
         resasc = wgk(21)*dabs(fc-reskh)  ! init  ! w8 * abs(fc - res Kronrod / 2 )
         do j=1,20
        	resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))  ! ite
    	 end do
    
    res = resk*xr
    !resabs = resabs*dabs(xr)
    !resasc = resasc*dabs(xr)
    !abserr = dabs((resk-resg)*xr)   ! estimation erreur Kronrod pra Gauss
    !if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00) abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00) ! si erreurs non nulles
    !if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1((epmach*0.5d+02)*resabs,abserr)

         end subroutine qgaussPL41

!=============================================================================================  
!=====QGAUS51 out a 51 point Gauss-Kronrod quadrature rule for weib ==========================
!=============================================================================================  


subroutine qgaussPL51weib(a,b,the01,the02,the12,res,v01,v02,v12)
         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resabs,resk,reskh,resasc,resg,v01,v02,v12,&
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

    	resg = fc*wgk(13)        ! init at wg(13)
        resk = fc*wgk(26)       ! init res Kronrod   ! fc * 8e poids Kronrod
        resabs = dabs(resk)    ! init res absolue Kronrod   
         
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
	       	resg = resg+wg(j)*(f1+f2)
                resk = resk + wgk(jtw)*(f1+f2)
                resabs = resabs + wgk(jtw)*(dabs(f1)+dabs(f2))
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
               resabs = resabs + wgk(jtwm1)*(dabs(f1)+dabs(f2))
              end do
	    

	  reskh = resk*0.5d+00         ! res Kronrod / 2
         resasc = wgk(21)*dabs(fc-reskh)  ! init  ! w8 * abs(fc - res Kronrod / 2 )
         do j=1,25
        	resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))  ! ite
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
         double precision::a,b,dx,xm,xr,res,resabs,resk,reskh,resasc,resg,v1,v2,v3,&
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

    	resg = fc*wgk(13)        ! init at wg(13)
        resk = fc*wgk(26)       ! init res Kronrod   ! fc * 8e poids Kronrod
        resabs = dabs(resk)    ! init res absolue Kronrod   
         

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
	       resg = resg+wg(j)*(f1+f2)
               resk = resk + wgk(jtw)*(f1+f2)
               resabs = resabs + wgk(jtw)*(dabs(f1)+dabs(f2))
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
               resabs = resabs + wgk(jtwm1)*(dabs(f1)+dabs(f2))
         end do
     ! pr calcul erreur 
         reskh = resk*0.5d+00         ! res Kronrod / 2
         resasc = wgk(21)*dabs(fc-reskh)  ! init  ! w8 * abs(fc - res Kronrod / 2 )
         do j=1,25
        	resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))  ! ite
    	 end do
    
    res = resk*xr
    !resabs = resabs*dabs(xr)
    !resasc = resasc*dabs(xr)
    !abserr = dabs((resk-resg)*xr)   ! estimation erreur Kronrod pra Gauss
    !if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00) abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00) ! si erreurs non nulles
    !if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1((epmach*0.5d+02)*resabs,abserr)

         end subroutine qgaussPL51

!=============================================================================================  
!=====QGAUS61 out a 61 point Gauss-Kronrod quadrature rule for weib ==========================
!=============================================================================================  


subroutine qgaussPL61weib(a,b,the01,the02,the12,res,v01,v02,v12)
         implicit none
         
         integer::j,jtw,jtwm1
         double precision::a,b,dx,xm,xr,res,resabs,resk,reskh,resasc,resg,v01,v02,v12,&
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

    	resg = 0       ! init at 0
        resk = fc*wgk(31)       ! init res Kronrod   ! fc * 8e poids Kronrod
        resabs = dabs(resk)    ! init res absolue Kronrod   
         
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
	       	resg = resg+wg(j)*(f1+f2)
                resk = resk + wgk(jtw)*(f1+f2)
                resabs = resabs + wgk(jtw)*(dabs(f1)+dabs(f2))
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
               resabs = resabs + wgk(jtwm1)*(dabs(f1)+dabs(f2))
              end do
	    

	  reskh = resk*0.5d+00         ! res Kronrod / 2
         resasc = wgk(31)*dabs(fc-reskh)  ! init  ! w8 * abs(fc - res Kronrod / 2 )
         do j=1,30
        	resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))  ! ite
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
         double precision::a,b,dx,xm,xr,res,resabs,resk,reskh,resasc,resg,v1,v2,v3,&
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

    	resg = 0        ! init at 0
        resk = fc*wgk(31)       ! init res Kronrod   ! fc * 8e poids Kronrod
        resabs = dabs(resk)    ! init res absolue Kronrod   
         

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
	       resg = resg+wg(j)*(f1+f2)
               resk = resk + wgk(jtw)*(f1+f2)
               resabs = resabs + wgk(jtw)*(dabs(f1)+dabs(f2))
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
               resabs = resabs + wgk(jtwm1)*(dabs(f1)+dabs(f2))
         end do
     ! pr calcul erreur 
         reskh = resk*0.5d+00         ! res Kronrod / 2
         resasc = wgk(31)*dabs(fc-reskh)  ! init  ! w8 * abs(fc - res Kronrod / 2 )
         do j=1,30
        	resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))  ! ite
    	 end do
    
    res = resk*xr
    !resabs = resabs*dabs(xr)
    !resasc = resasc*dabs(xr)
    !abserr = dabs((resk-resg)*xr)   ! estimation erreur Kronrod pra Gauss
    !if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00) abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00) ! si erreurs non nulles
    !if(resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1((epmach*0.5d+02)*resabs,abserr)

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
         double precision dx,xm,xr,resabsdenum,reskdenum,&
         reskhdenum,resascdenum,resgdenum,resdenum, & 
	resabs01num,resk01num,reskh01num,resasc01num,resg01num,res01num, & 
	resabs02num,resk02num,reskh02num,resasc02num,resg02num,res02num, & 
	resabs12num,resk12num,reskh12num,resasc12num,resg12num,res12num, & 
	resabs0101num,resk0101num,reskh0101num,resasc0101num,resg0101num,res0101num, &
	resabs0102num,resk0102num,reskh0102num,resasc0102num,resg0102num,res0102num, &
	resabs0112num,resk0112num,reskh0112num,resasc0112num,resg0112num,res0112num, &
	resabs0202num,resk0202num,reskh0202num,resasc0202num,resg0202num,res0202num, &
	resabs0212num,resk0212num,reskh0212num,resasc0212num,resg0212num,res0212num, &
	resabs1212num,resk1212num,reskh1212num,resasc1212num,resg1212num,res1212num, &
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

    		resgdenum = fcdenum*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	reskdenum = fcdenum*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabsdenum = dabs(reskdenum)    ! init res absolue Kronrod  

		resg01num = fc01num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk01num = fc01num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs01num = dabs(resk01num)    ! init res absolue Kronrod   

		resg02num = fc02num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk02num = fc02num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs02num = dabs(resk02num)    ! init res absolue Kronrod   

		resg12num = fc12num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk12num = fc12num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs12num = dabs(resk12num)    ! init res absolue Kronrod  

		resg0101num = fc0101num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk0101num = fc0101num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs0101num = dabs(resk0101num)    ! init res absolue Kronrod   
		
		resg0102num = fc0102num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk0102num = fc0102num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs0102num = dabs(resk0102num)    ! init res absolue Kronrod   

		resg0112num = fc0112num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk0112num = fc0112num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs0112num = dabs(resk0112num)    ! init res absolue Kronrod   

		resg0202num = fc0202num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk0202num = fc0202num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs0202num = dabs(resk0202num)    ! init res absolue Kronrod   

		resg0212num = fc0212num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk0212num = fc0212num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs0212num = dabs(resk0212num)    ! init res absolue Kronrod   


		resg1212num = fc1212num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk1212num = fc1212num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs1212num = dabs(resk1212num)    ! init res absolue Kronrod   






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

	       		resgdenum = resgdenum+wg(j)*(f1denum+f2denum)
               		reskdenum = reskdenum + wgk(jtw)*(f1denum+f2denum)
               		resabsdenum = resabsdenum + wgk(jtw)*(dabs(f1denum)+&
               		dabs(f2denum))

			fv101num(jtw) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtw) = f201num   ! svgrd valeurs fct f a drte du centre

	       		resg01num = resg01num+wg(j)*(f101num+f201num)
               		resk01num = resk01num + wgk(jtw)*(f101num+f201num)
               		resabs01num = resabs01num + wgk(jtw)*(dabs(f101num)+&
               		dabs(f201num))

			fv102num(jtw) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtw) = f202num   ! svgrd valeurs fct f a drte du centre

	       		resg02num = resg02num+wg(j)*(f102num+f202num)
               		resk02num = resk02num + wgk(jtw)*(f102num+f202num)
               		resabs02num = resabs02num + wgk(jtw)*(dabs(f102num)+&
               		dabs(f202num))

			fv112num(jtw) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtw) = f212num   ! svgrd valeurs fct f a drte du centre

	       		resg12num = resg12num+wg(j)*(f112num+f212num)
               		resk12num = resk12num + wgk(jtw)*(f112num+f212num)
               		resabs12num = resabs12num + wgk(jtw)*(dabs(f112num)+&
               		dabs(f212num))

			fv10101num(jtw) = f10101num   ! svgrd valeurs fct f a gche du centre
               		fv20101num(jtw) = f20101num   ! svgrd valeurs fct f a drte du centre

	       		resg0101num = resg0101num+wg(j)*(f10101num+f20101num)
               		resk0101num = resk0101num + wgk(jtw)*(f10101num+f20101num)
               		resabs0101num = resabs0101num + wgk(jtw)*(dabs(f10101num)+&
               		dabs(f20101num))

			fv10102num(jtw) = f10102num   ! svgrd valeurs fct f a gche du centre
               		fv20102num(jtw) = f20102num   ! svgrd valeurs fct f a drte du centre

	       		resg0102num = resg0102num+wg(j)*(f10102num+f20102num)
               		resk0102num = resk0102num + wgk(jtw)*(f10102num+f20102num)
               		resabs0102num = resabs0102num + wgk(jtw)*(dabs(f10102num)+&
               		dabs(f20102num))

			fv10112num(jtw) = f10112num   ! svgrd valeurs fct f a gche du centre
               		fv20112num(jtw) = f20112num   ! svgrd valeurs fct f a drte du centre

	       		resg0112num = resg0112num+wg(j)*(f10112num+f20112num)
               		resk0112num = resk0112num + wgk(jtw)*(f10112num+f20112num)
               		resabs0112num = resabs0112num + wgk(jtw)*(dabs(f10112num)+&
               		dabs(f20112num))

			fv10202num(jtw) = f10202num   ! svgrd valeurs fct f a gche du centre
               		fv20202num(jtw) = f20202num   ! svgrd valeurs fct f a drte du centre

	       		resg0202num = resg0202num+wg(j)*(f10202num+f20202num)
               		resk0202num = resk0202num + wgk(jtw)*(f10202num+f20202num)
               		resabs0202num = resabs0202num + wgk(jtw)*(dabs(f10202num)+&
               		dabs(f20202num))

			fv10212num(jtw) = f10212num   ! svgrd valeurs fct f a gche du centre
               		fv20212num(jtw) = f20212num   ! svgrd valeurs fct f a drte du centre

	       		resg0212num = resg0212num+wg(j)*(f10212num+f20212num)
               		resk0212num = resk0212num + wgk(jtw)*(f10212num+f20212num)
               		resabs0212num = resabs0212num + wgk(jtw)*(dabs(f10212num)+&
               		dabs(f20212num))

			fv11212num(jtw) = f11212num   ! svgrd valeurs fct f a gche du centre
               		fv21212num(jtw) = f21212num   ! svgrd valeurs fct f a drte du centre

	       		resg1212num = resg1212num+wg(j)*(f11212num+f21212num)
               		resk1212num = resk1212num + wgk(jtw)*(f11212num+f21212num)
               		resabs1212num = resabs1212num + wgk(jtw)*(dabs(f11212num)+&
               		dabs(f21212num))

			
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
               		resabsdenum = resabsdenum + wgk(jtwm1)*(dabs(f1denum)+&
               		dabs(f2denum))

			fv101num(jtwm1) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtwm1) = f201num   ! svgrd valeurs fct f a drte du centre
	       		resk01num = resk01num + wgk(jtwm1)*(f101num+f201num)
               		resabs01num = resabs01num + wgk(jtwm1)*(dabs(f101num)+&
               		dabs(f201num))

			fv102num(jtwm1) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtwm1) = f202num   ! svgrd valeurs fct f a drte du centre
	       		resk02num = resk02num + wgk(jtwm1)*(f102num+f202num)
               		resabs02num = resabs02num + wgk(jtwm1)*(dabs(f102num)+&
               		dabs(f202num))

			fv112num(jtwm1) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtwm1) = f212num   ! svgrd valeurs fct f a drte du centre
	       		resk12num = resk12num + wgk(jtwm1)*(f112num+f212num)
               		resabs12num = resabs12num + wgk(jtwm1)*(dabs(f112num)+&
               		dabs(f212num))

			fv10101num(jtwm1) = f10101num   ! svgrd valeurs fct f a gche du centre
               		fv20101num(jtwm1) = f20101num   ! svgrd valeurs fct f a drte du centre
               		resk0101num = resk0101num + wgk(jtwm1)*(f10101num+f20101num)
               		resabs0101num = resabs0101num + wgk(jtwm1)*(dabs(f10101num)+&
               		dabs(f20101num))


			fv10102num(jtwm1) = f10102num   ! svgrd valeurs fct f a gche du centre
               		fv20102num(jtwm1) = f20102num   ! svgrd valeurs fct f a drte du centre
               		resk0102num = resk0102num + wgk(jtwm1)*(f10102num+f20102num)
               		resabs0102num = resabs0102num + wgk(jtwm1)*(dabs(f10102num)+&
               		dabs(f20102num))

			fv10112num(jtwm1) = f10112num   ! svgrd valeurs fct f a gche du centre
               		fv20112num(jtwm1) = f20112num   ! svgrd valeurs fct f a drte du centre
               		resk0112num = resk0112num + wgk(jtwm1)*(f10112num+f20112num)
               		resabs0112num = resabs0112num + wgk(jtwm1)*(dabs(f10112num)+&
               		dabs(f20112num))

			fv10202num(jtwm1) = f10202num   ! svgrd valeurs fct f a gche du centre
               		fv20202num(jtwm1) = f20202num   ! svgrd valeurs fct f a drte du centre
               		resk0202num = resk0202num + wgk(jtwm1)*(f10202num+f20202num)
               		resabs0202num = resabs0202num + wgk(jtwm1)*(dabs(f10202num)+&
               		dabs(f20202num))

			fv10212num(jtwm1) = f10212num   ! svgrd valeurs fct f a gche du centre
               		fv20212num(jtwm1) = f20212num   ! svgrd valeurs fct f a drte du centre
               		resk0212num = resk0212num + wgk(jtwm1)*(f10212num+f20212num)
               		resabs0212num = resabs0212num + wgk(jtwm1)*(dabs(f10212num)+&
               		dabs(f20212num))

			fv11212num(jtwm1) = f11212num   ! svgrd valeurs fct f a gche du centre
               		fv21212num(jtwm1) = f21212num   ! svgrd valeurs fct f a drte du centre
               		resk1212num = resk1212num + wgk(jtwm1)*(f11212num+f21212num)
               		resabs1212num = resabs1212num + wgk(jtwm1)*(dabs(f11212num)+&
               		dabs(f21212num))




			
         	end do
     ! pr calcul erreur 
         	reskhdenum = reskdenum*0.5d+00         ! res Kronrod / 2
         	resascdenum = wgk(8)*dabs(fcdenum-reskhdenum)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

     ! pr calcul erreur 
         	reskh01num = resk01num*0.5d+00         ! res Kronrod / 2
         	resasc01num = wgk(8)*dabs(fc01num-reskh01num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

		     ! pr calcul erreur 
         	reskh02num = resk02num*0.5d+00         ! res Kronrod / 2
         	resasc02num = wgk(8)*dabs(fc02num-reskh02num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

    ! pr calcul erreur 
         	reskh12num = resk12num*0.5d+00         ! res Kronrod / 2
         	resasc12num = wgk(8)*dabs(fc12num-reskh12num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

		reskh0101num = resk0101num*0.5d+00         ! res Kronrod / 2
         	resasc0101num = wgk(8)*dabs(fc0101num-reskh0101num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

		reskh0102num = resk0102num*0.5d+00         ! res Kronrod / 2
         	resasc0102num = wgk(8)*dabs(fc0102num-reskh0102num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

		reskh0112num = resk0112num*0.5d+00         ! res Kronrod / 2
         	resasc0112num = wgk(8)*dabs(fc0112num-reskh0112num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

		reskh0202num = resk0202num*0.5d+00         ! res Kronrod / 2
         	resasc0202num = wgk(8)*dabs(fc0202num-reskh0202num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

		reskh0212num = resk0212num*0.5d+00         ! res Kronrod / 2
         	resasc0212num = wgk(8)*dabs(fc0212num-reskh0212num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )


		reskh1212num = resk1212num*0.5d+00         ! res Kronrod / 2
         	resasc1212num = wgk(8)*dabs(fc1212num-reskh1212num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

		
		do j=1,7
        		resascdenum = resascdenum+wgk(j)*(dabs(fv1denum(j)-&
        		reskhdenum)+dabs(fv2denum(j)-reskhdenum))  ! ite
			resascdenum = resascdenum+wgk(j)*(dabs(fv1denum(j)-&
			reskhdenum)+dabs(fv2denum(j)-reskhdenum))  ! ite

        		resasc01num = resasc01num+wgk(j)*(dabs(fv101num(j)-&
        		reskh01num)+dabs(fv201num(j)-reskh01num))  ! ite
			resasc01num = resasc01num+wgk(j)*(dabs(fv101num(j)-&
			reskh01num)+dabs(fv201num(j)-reskh01num))  ! ite

        		resasc02num = resasc02num+wgk(j)*(dabs(fv102num(j)-&
        		reskh02num)+dabs(fv202num(j)-reskh02num))  ! ite
			resasc02num = resasc02num+wgk(j)*(dabs(fv102num(j)-&
			reskh02num)+dabs(fv202num(j)-reskh02num))  ! ite

        		resasc12num = resasc12num+wgk(j)*(dabs(fv112num(j)-&
        		reskh12num)+dabs(fv212num(j)-reskh12num))  ! ite
			resasc12num = resasc12num+wgk(j)*(dabs(fv112num(j)-&
			reskh12num)+dabs(fv212num(j)-reskh12num))  ! ite

			resasc0101num = resasc0101num+wgk(j)*(dabs(fv10101num(j)-&
        		reskh0101num)+dabs(fv20101num(j)-reskh0101num))  ! ite
			resasc0101num = resasc0101num+wgk(j)*(dabs(fv10101num(j)-&
			reskh0101num)+dabs(fv20101num(j)-reskh0101num))  ! ite

			resasc0102num = resasc0102num+wgk(j)*(dabs(fv10102num(j)-&
        		reskh0102num)+dabs(fv20102num(j)-reskh0102num))  ! ite
			resasc0102num = resasc0102num+wgk(j)*(dabs(fv10102num(j)-&
			reskh0102num)+dabs(fv20102num(j)-reskh0102num))  ! ite

			resasc0112num = resasc0112num+wgk(j)*(dabs(fv10112num(j)-&
        		reskh0112num)+dabs(fv20112num(j)-reskh0112num))  ! ite
			resasc0112num = resasc0112num+wgk(j)*(dabs(fv10112num(j)-&
			reskh0112num)+dabs(fv20112num(j)-reskh0112num))  ! ite


			resasc0202num = resasc0202num+wgk(j)*(dabs(fv10202num(j)-&
        		reskh0202num)+dabs(fv20202num(j)-reskh0202num))  ! ite
			resasc0202num = resasc0202num+wgk(j)*(dabs(fv10202num(j)-&
			reskh0202num)+dabs(fv20202num(j)-reskh0202num))  ! ite


			resasc0212num = resasc0212num+wgk(j)*(dabs(fv10212num(j)-&
        		reskh0212num)+dabs(fv20212num(j)-reskh0212num))  ! ite
			resasc0212num = resasc0212num+wgk(j)*(dabs(fv10212num(j)-&
			reskh0212num)+dabs(fv20212num(j)-reskh0212num))  ! ite


        		
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


!================================  QGAUS : 1  15  ==========================


subroutine qgausssplinederiv(a,b,the01,the02,the12,resdenum,&
		res01num,res02num,res12num,res0101num,res0102num,res0112num,&
		res0202num,res0212num,res1212num,v01,v02,v12)

   use commun,only:zi01,zi12,zi02,nz01,nz12,nz02

        double precision a,b
         double precision dx,xm,xr,resabsdenum,reskdenum,&
         reskhdenum,resascdenum,resgdenum,resdenum, & 
	resabs01num,resk01num,reskh01num,resasc01num,resg01num,res01num, & 
	resabs02num,resk02num,reskh02num,resasc02num,resg02num,res02num, & 
	resabs12num,resk12num,reskh12num,resasc12num,resg12num,res12num, & 
	resabs0101num,resk0101num,reskh0101num,resasc0101num,resg0101num,res0101num, &
	resabs0102num,resk0102num,reskh0102num,resasc0102num,resg0102num,res0102num, &
	resabs0112num,resk0112num,reskh0112num,resasc0112num,resg0112num,res0112num, &
	resabs0202num,resk0202num,reskh0202num,resasc0202num,resg0202num,res0202num, &
	resabs0212num,resk0212num,reskh0212num,resasc0212num,resg0212num,res0212num, &
	resabs1212num,resk1212num,reskh1212num,resasc1212num,resg1212num,res1212num, &
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

    		resgdenum = fcdenum*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	reskdenum = fcdenum*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabsdenum = dabs(reskdenum)    ! init res absolue Kronrod  

		resg01num = fc01num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk01num = fc01num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs01num = dabs(resk01num)    ! init res absolue Kronrod   

		resg02num = fc02num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk02num = fc02num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs02num = dabs(resk02num)    ! init res absolue Kronrod   

		resg12num = fc12num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk12num = fc12num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs12num = dabs(resk12num)    ! init res absolue Kronrod  

		resg0101num = fc0101num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk0101num = fc0101num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs0101num = dabs(resk0101num)    ! init res absolue Kronrod   
		
		resg0102num = fc0102num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk0102num = fc0102num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs0102num = dabs(resk0102num)    ! init res absolue Kronrod   

		resg0112num = fc0112num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk0112num = fc0112num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs0112num = dabs(resk0112num)    ! init res absolue Kronrod   

		resg0202num = fc0202num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk0202num = fc0202num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs0202num = dabs(resk0202num)    ! init res absolue Kronrod   

		resg0212num = fc0212num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk0212num = fc0212num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs0212num = dabs(resk0212num)    ! init res absolue Kronrod   


		resg1212num = fc1212num*wg(4)        ! init res Gauss  ! fc * 4e poids Gauss
        	resk1212num = fc1212num*wgk(8)       ! init res Kronrod   ! fc * 8e poids Kronrod
        	resabs1212num = dabs(resk1212num)    ! init res absolue Kronrod   






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

	       		resgdenum = resgdenum+wg(j)*(f1denum+f2denum)
               		reskdenum = reskdenum + wgk(jtw)*(f1denum+f2denum)
               		resabsdenum = resabsdenum + wgk(jtw)*(dabs(f1denum)+&
               		dabs(f2denum))

			fv101num(jtw) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtw) = f201num   ! svgrd valeurs fct f a drte du centre

	       		resg01num = resg01num+wg(j)*(f101num+f201num)
               		resk01num = resk01num + wgk(jtw)*(f101num+f201num)
               		resabs01num = resabs01num + wgk(jtw)*(dabs(f101num)+&
               		dabs(f201num))

			fv102num(jtw) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtw) = f202num   ! svgrd valeurs fct f a drte du centre

	       		resg02num = resg02num+wg(j)*(f102num+f202num)
               		resk02num = resk02num + wgk(jtw)*(f102num+f202num)
               		resabs02num = resabs02num + wgk(jtw)*(dabs(f102num)+&
               		dabs(f202num))

			fv112num(jtw) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtw) = f212num   ! svgrd valeurs fct f a drte du centre

	       		resg12num = resg12num+wg(j)*(f112num+f212num)
               		resk12num = resk12num + wgk(jtw)*(f112num+f212num)
               		resabs12num = resabs12num + wgk(jtw)*(dabs(f112num)+&
               		dabs(f212num))

			fv10101num(jtw) = f10101num   ! svgrd valeurs fct f a gche du centre
               		fv20101num(jtw) = f20101num   ! svgrd valeurs fct f a drte du centre

	       		resg0101num = resg0101num+wg(j)*(f10101num+f20101num)
               		resk0101num = resk0101num + wgk(jtw)*(f10101num+f20101num)
               		resabs0101num = resabs0101num + wgk(jtw)*(dabs(f10101num)+&
               		dabs(f20101num))

			fv10102num(jtw) = f10102num   ! svgrd valeurs fct f a gche du centre
               		fv20102num(jtw) = f20102num   ! svgrd valeurs fct f a drte du centre

	       		resg0102num = resg0102num+wg(j)*(f10102num+f20102num)
               		resk0102num = resk0102num + wgk(jtw)*(f10102num+f20102num)
               		resabs0102num = resabs0102num + wgk(jtw)*(dabs(f10102num)+&
               		dabs(f20102num))

			fv10112num(jtw) = f10112num   ! svgrd valeurs fct f a gche du centre
               		fv20112num(jtw) = f20112num   ! svgrd valeurs fct f a drte du centre

	       		resg0112num = resg0112num+wg(j)*(f10112num+f20112num)
               		resk0112num = resk0112num + wgk(jtw)*(f10112num+f20112num)
               		resabs0112num = resabs0112num + wgk(jtw)*(dabs(f10112num)+&
               		dabs(f20112num))

			fv10202num(jtw) = f10202num   ! svgrd valeurs fct f a gche du centre
               		fv20202num(jtw) = f20202num   ! svgrd valeurs fct f a drte du centre

	       		resg0202num = resg0202num+wg(j)*(f10202num+f20202num)
               		resk0202num = resk0202num + wgk(jtw)*(f10202num+f20202num)
               		resabs0202num = resabs0202num + wgk(jtw)*(dabs(f10202num)+&
               		dabs(f20202num))

			fv10212num(jtw) = f10212num   ! svgrd valeurs fct f a gche du centre
               		fv20212num(jtw) = f20212num   ! svgrd valeurs fct f a drte du centre

	       		resg0212num = resg0212num+wg(j)*(f10212num+f20212num)
               		resk0212num = resk0212num + wgk(jtw)*(f10212num+f20212num)
               		resabs0212num = resabs0212num + wgk(jtw)*(dabs(f10212num)+&
               		dabs(f20212num))

			fv11212num(jtw) = f11212num   ! svgrd valeurs fct f a gche du centre
               		fv21212num(jtw) = f21212num   ! svgrd valeurs fct f a drte du centre

	       		resg1212num = resg1212num+wg(j)*(f11212num+f21212num)
               		resk1212num = resk1212num + wgk(jtw)*(f11212num+f21212num)
               		resabs1212num = resabs1212num + wgk(jtw)*(dabs(f11212num)+&
               		dabs(f21212num))

			
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
               		resabsdenum = resabsdenum + wgk(jtwm1)*(dabs(f1denum)+&
               		dabs(f2denum))

			fv101num(jtwm1) = f101num   ! svgrd valeurs fct f a gche du centre
               		fv201num(jtwm1) = f201num   ! svgrd valeurs fct f a drte du centre
	       		resk01num = resk01num + wgk(jtwm1)*(f101num+f201num)
               		resabs01num = resabs01num + wgk(jtwm1)*(dabs(f101num)+&
               		dabs(f201num))

			fv102num(jtwm1) = f102num   ! svgrd valeurs fct f a gche du centre
               		fv202num(jtwm1) = f202num   ! svgrd valeurs fct f a drte du centre
	       		resk02num = resk02num + wgk(jtwm1)*(f102num+f202num)
               		resabs02num = resabs02num + wgk(jtwm1)*(dabs(f102num)+&
               		dabs(f202num))

			fv112num(jtwm1) = f112num   ! svgrd valeurs fct f a gche du centre
               		fv212num(jtwm1) = f212num   ! svgrd valeurs fct f a drte du centre
	       		resk12num = resk12num + wgk(jtwm1)*(f112num+f212num)
               		resabs12num = resabs12num + wgk(jtwm1)*(dabs(f112num)+&
               		dabs(f212num))

			fv10101num(jtwm1) = f10101num   ! svgrd valeurs fct f a gche du centre
               		fv20101num(jtwm1) = f20101num   ! svgrd valeurs fct f a drte du centre
               		resk0101num = resk0101num + wgk(jtwm1)*(f10101num+f20101num)
               		resabs0101num = resabs0101num + wgk(jtwm1)*(dabs(f10101num)+&
               		dabs(f20101num))


			fv10102num(jtwm1) = f10102num   ! svgrd valeurs fct f a gche du centre
               		fv20102num(jtwm1) = f20102num   ! svgrd valeurs fct f a drte du centre
               		resk0102num = resk0102num + wgk(jtwm1)*(f10102num+f20102num)
               		resabs0102num = resabs0102num + wgk(jtwm1)*(dabs(f10102num)+&
               		dabs(f20102num))

			fv10112num(jtwm1) = f10112num   ! svgrd valeurs fct f a gche du centre
               		fv20112num(jtwm1) = f20112num   ! svgrd valeurs fct f a drte du centre
               		resk0112num = resk0112num + wgk(jtwm1)*(f10112num+f20112num)
               		resabs0112num = resabs0112num + wgk(jtwm1)*(dabs(f10112num)+&
               		dabs(f20112num))

			fv10202num(jtwm1) = f10202num   ! svgrd valeurs fct f a gche du centre
               		fv20202num(jtwm1) = f20202num   ! svgrd valeurs fct f a drte du centre
               		resk0202num = resk0202num + wgk(jtwm1)*(f10202num+f20202num)
               		resabs0202num = resabs0202num + wgk(jtwm1)*(dabs(f10202num)+&
               		dabs(f20202num))

			fv10212num(jtwm1) = f10212num   ! svgrd valeurs fct f a gche du centre
               		fv20212num(jtwm1) = f20212num   ! svgrd valeurs fct f a drte du centre
               		resk0212num = resk0212num + wgk(jtwm1)*(f10212num+f20212num)
               		resabs0212num = resabs0212num + wgk(jtwm1)*(dabs(f10212num)+&
               		dabs(f20212num))

			fv11212num(jtwm1) = f11212num   ! svgrd valeurs fct f a gche du centre
               		fv21212num(jtwm1) = f21212num   ! svgrd valeurs fct f a drte du centre
               		resk1212num = resk1212num + wgk(jtwm1)*(f11212num+f21212num)
               		resabs1212num = resabs1212num + wgk(jtwm1)*(dabs(f11212num)+&
               		dabs(f21212num))




			
         	end do
     ! pr calcul erreur 
         	reskhdenum = reskdenum*0.5d+00         ! res Kronrod / 2
         	resascdenum = wgk(8)*dabs(fcdenum-reskhdenum)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

     ! pr calcul erreur 
         	reskh01num = resk01num*0.5d+00         ! res Kronrod / 2
         	resasc01num = wgk(8)*dabs(fc01num-reskh01num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

		     ! pr calcul erreur 
         	reskh02num = resk02num*0.5d+00         ! res Kronrod / 2
         	resasc02num = wgk(8)*dabs(fc02num-reskh02num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

    ! pr calcul erreur 
         	reskh12num = resk12num*0.5d+00         ! res Kronrod / 2
         	resasc12num = wgk(8)*dabs(fc12num-reskh12num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

		reskh0101num = resk0101num*0.5d+00         ! res Kronrod / 2
         	resasc0101num = wgk(8)*dabs(fc0101num-reskh0101num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

		reskh0102num = resk0102num*0.5d+00         ! res Kronrod / 2
         	resasc0102num = wgk(8)*dabs(fc0102num-reskh0102num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

		reskh0112num = resk0112num*0.5d+00         ! res Kronrod / 2
         	resasc0112num = wgk(8)*dabs(fc0112num-reskh0112num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

		reskh0202num = resk0202num*0.5d+00         ! res Kronrod / 2
         	resasc0202num = wgk(8)*dabs(fc0202num-reskh0202num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

		reskh0212num = resk0212num*0.5d+00         ! res Kronrod / 2
         	resasc0212num = wgk(8)*dabs(fc0212num-reskh0212num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )


		reskh1212num = resk1212num*0.5d+00         ! res Kronrod / 2
         	resasc1212num = wgk(8)*dabs(fc1212num-reskh1212num)  ! init  ! w8 * abs(fc - res Kronrod / 2 )

		
		do j=1,7
        		resascdenum = resascdenum+wgk(j)*(dabs(fv1denum(j)-&
        		reskhdenum)+dabs(fv2denum(j)-reskhdenum))  ! ite
			resascdenum = resascdenum+wgk(j)*(dabs(fv1denum(j)-&
			reskhdenum)+dabs(fv2denum(j)-reskhdenum))  ! ite

        		resasc01num = resasc01num+wgk(j)*(dabs(fv101num(j)-&
        		reskh01num)+dabs(fv201num(j)-reskh01num))  ! ite
			resasc01num = resasc01num+wgk(j)*(dabs(fv101num(j)-&
			reskh01num)+dabs(fv201num(j)-reskh01num))  ! ite

        		resasc02num = resasc02num+wgk(j)*(dabs(fv102num(j)-&
        		reskh02num)+dabs(fv202num(j)-reskh02num))  ! ite
			resasc02num = resasc02num+wgk(j)*(dabs(fv102num(j)-&
			reskh02num)+dabs(fv202num(j)-reskh02num))  ! ite

        		resasc12num = resasc12num+wgk(j)*(dabs(fv112num(j)-&
        		reskh12num)+dabs(fv212num(j)-reskh12num))  ! ite
			resasc12num = resasc12num+wgk(j)*(dabs(fv112num(j)-&
			reskh12num)+dabs(fv212num(j)-reskh12num))  ! ite

			resasc0101num = resasc0101num+wgk(j)*(dabs(fv10101num(j)-&
        		reskh0101num)+dabs(fv20101num(j)-reskh0101num))  ! ite
			resasc0101num = resasc0101num+wgk(j)*(dabs(fv10101num(j)-&
			reskh0101num)+dabs(fv20101num(j)-reskh0101num))  ! ite

			resasc0102num = resasc0102num+wgk(j)*(dabs(fv10102num(j)-&
        		reskh0102num)+dabs(fv20102num(j)-reskh0102num))  ! ite
			resasc0102num = resasc0102num+wgk(j)*(dabs(fv10102num(j)-&
			reskh0102num)+dabs(fv20102num(j)-reskh0102num))  ! ite

			resasc0112num = resasc0112num+wgk(j)*(dabs(fv10112num(j)-&
        		reskh0112num)+dabs(fv20112num(j)-reskh0112num))  ! ite
			resasc0112num = resasc0112num+wgk(j)*(dabs(fv10112num(j)-&
			reskh0112num)+dabs(fv20112num(j)-reskh0112num))  ! ite


			resasc0202num = resasc0202num+wgk(j)*(dabs(fv10202num(j)-&
        		reskh0202num)+dabs(fv20202num(j)-reskh0202num))  ! ite
			resasc0202num = resasc0202num+wgk(j)*(dabs(fv10202num(j)-&
			reskh0202num)+dabs(fv20202num(j)-reskh0202num))  ! ite


			resasc0212num = resasc0212num+wgk(j)*(dabs(fv10212num(j)-&
        		reskh0212num)+dabs(fv20212num(j)-reskh0212num))  ! ite
			resasc0212num = resasc0212num+wgk(j)*(dabs(fv10212num(j)-&
			reskh0212num)+dabs(fv20212num(j)-reskh0212num))  ! ite


        		
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



!=============================================================================================  
!======================= Calculate derivatives of loglik with weibull baseline risk ==========
!=============================================================================================  


subroutine derivaweib(b0,np0,npar0,bfix0,fix0,c0,no0,ve010,ve120,ve020,&
        dimnva01,dimnva12,dimnva02,nva01,nva12,nva02,t00,&
        t10,t20,t30,troncature0,weib0,likelihood_deriv)
	
	use commun
        implicit none
         
        double precision::res2denum,res201num,res202num,res212num, &
	res20101num,res20102num,res20112num,res20202num,res20212num,&
	res21212num,vet01,vet12,vet02,resint,v,u1,u2,u3
        integer::np0,i,j,l,w,k,lfix, kfix,npar0,nva01,nva12,nva02,no0, &
	weib0,nz010,nz020,nz120,troncature0,dimnva01,dimnva02,dimnva12, & 
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
	weib=weib0
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

   
	


	if(weib.eq.1)then
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
	else 
	do i=1,2
            the01(i)=dexp(bh(i))
         end do
         do i=1,2
            j = 2+i
            the02(i)=dexp(bh(j))
         end do
         do i=1,2
            j = 4+i
            the12(i)=dexp(bh(j))
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
	 
	deallocate(b,bfix,fix,ve01,ve02,ve12,ve01nofix,&
	ve02nofix,ve12nofix,tronc01,tronc02,t0,t1,t2,t3,c,&
	ve01square,ve02square,ve12square,&
	tronc01square,tronc02square)     

    end subroutine derivaweib


!=============================================================================================  
!======================= Calculate derivatives of loglik with M-splines baseline risk ========
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









