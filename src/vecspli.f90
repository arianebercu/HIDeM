!========================== VECSPLI ==============================

        subroutine vecspli(n,no,zi,t,mm3,mm2,mm1,mm,im3,im2,im1,im)
        
        implicit none
        
        integer::no,n
        double precision,dimension(no)::t        
        double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2, & 
        h3,h4,h3m,h2n,hn,hh3,hh2
        integer::i,j,k
        double precision,dimension(-2:(n+1)),intent(inout)::zi
        double precision,dimension(no),intent(inout)::mm3,mm2,mm1,mm,im3,im2,im1,im

!----------  calcul de u(ti) ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)

!       write(*,*)'taille de zi dans vecspli',size(zi)
        i=0
        k=0
        j=0     
        do i=1,no
                do k = 2,n-2
                        if ((t(i).ge.zi(k-1)).and.(t(i).lt.zi(k))) then
                                j = k-1
                        endif
                end do 
        
                ht = t(i)-zi(j)
                htm= t(i)-zi(j-1)
                h2t= t(i)-zi(j+2)       
                ht2 = zi(j+1)-t(i)
                ht3 = zi(j+3)-t(i)
                hht = t(i)-zi(j-2)
        
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
                mm3(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
                mm2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
                *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
                mm1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
                h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
                mm(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
        
                im3(i) = (0.25d0*(t(i)-zi(j-3))*mm3(i))+(0.25d0*hh2*mm2(i))+(0.25d0*h3m*mm1(i))+(0.25d0*h4*mm(i))
        
                im2(i) = (0.25d0*hht*mm2(i))+(h3m*mm1(i)*0.25d0) &
                        +(h4*mm(i)*0.25d0)
                im1(i) = (htm*mm1(i)*0.25d0)+(h4*mm(i)*0.25d0)
                im(i)  = ht*mm(i)*0.25d0
        end do

        end subroutine vecspli