module commun
        integer,save::nn
!-------------------mem1 & mem2----------------------------------------
         double precision,dimension(:),allocatable,save::mm3a,mm2a,mm1a,mma, &
         im3a,im2a,im1a,ima
!-------------------  mem3 & mem4 ----------------------------------
         double precision,dimension(:),allocatable,save::mm3b,mm2b,mm1b,mmb, &
         im3b,im2b,im1b,imb
!-----------------  mem5 & mem6-----------------------------------------
         double precision,dimension(:),allocatable,save::mm3c,mm2c,mm1c,mmc, &
         im3c,im2c,im1c,imc
!-------------------  dace1 -------------------------------------------
        double precision,dimension(:),allocatable,save::zi01,zi12,zi02
!-------------------  dace2 -------------------------------------------
        integer::no
!-------------------  dace3 -------------------------------------------

        integer::nz01,nz12,nz02
!-------------------  ve1 -------------------------------------------
        double precision,dimension(:,:),allocatable,save::ve01,ve12,ve02,ve01nofix, &
	ve02nofix,ve12nofix,ve01square,ve02square,ve12square
	double precision,dimension(:),allocatable,save::tronc01, tronc02, & 
	tronc01square, tronc02square

!-------------------  dace1new -------------------------------------------
        double precision,dimension(:),allocatable,save::t0,t1,t2,t3,t4
        integer,dimension(:),allocatable,save::c
!-------------------  dace1new -------------------------------------------
        integer,save::troncature,ind_hess
        double precision,dimension(:,:),allocatable,save::hessienne     
!-------------------  dace1 -------------------------------------------
        double precision,dimension(:),allocatable,save::zi
!-------------------  dace -------------------------------------------
        integer::nz,verSurv
!-------------------  ve1 -------------------------------------------
        double precision,dimension(:,:),allocatable,save::ve
        double precision::k0surv
!---------------------- ca, bd and dd -----------------------------
	integer,dimension(:),allocatable,save::fix
        double precision,dimension(:),allocatable,save::b,bfix,bfix1
        double precision,dimension(:),allocatable,save::eta
!-------------------- gausspoints ------------------------------
	integer::gausspoint
!-------------------- weibullpara ------------------------------
	integer::weib


end module commun

module tailles
        integer,save::np
end module tailles

