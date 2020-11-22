       subroutine chis2_min_2020(min_th,min_nosc,min_data,min_chi2_pull)
       implicit none
!this subroutine uses the shifted evets given by 
! N_bin_bar= evbinOSC12*(1+sum(i))
! where sum = Xpull(j)*errmatk2k12(n,j)
      include 'MINOSPULL.inc'

        real*8  min_th(40),min_nosc(40),min_data(40),min_chi2_pull

       real*8 Xpull(7),F(7),PAR(1),FNORM,WK(5000) 
       real*8 suma,sumpull, sumatix1,sumatix2	
       real*8 sumcoef,sumcoef2
       real*8 N_th_tot_bar, N_bin_bar(40),N_data_tot
       real*8 mat(7,7),sigma_minos(7)
       real*8 chi2_shape

       integer i,j,k,h,n
       
       external FCN_minos_20

       
       
! errores systematicos obtenidos de sep 2008  PRL 101, 131802 (2008)
c (a) Absolute hadronic E scale (10:3%) 
c (b) Relative hadronic E scale (3:3%)
c (c) Normalization (4%)
c (d) NC contamination (50%)
c (e) " muon momentum (range 2%, curvature 3%)
c (f ) cross section   < 10 GeVÞ (12%)
c (g) Beam flux    free

c        sigma_minos(1)=0.103                !  a) Absolute hadronic E scale (10:3%) 
c        sigma_minos(2)=0.033                !  b) Relative hadronic E scale (3:3%)
c        sigma_minos(3)=0.04                 !  c) Normalization (4%)
c        sigma_minos(4)=0.50                 !  d) NC contamination (50%)
c        sigma_minos(5)=0.03                 !  e) " muon momentum (range 2%, curvature 3%)
c        sigma_minos(6)=0.12                 !  f ) cross section   < 10 GeVÞ (12%)
       
 ! errores systematicos en analisis de 2020 estan en PHYSICAL REVIEW LETTERS week ending 
 !PRL 106, 181801 (2011)
 
c (a) Hadronic energy                                      anterior    10.3%
c (b) "muon  energy (range 2%, curv. 3%)     anterior  (range 2%, curvature 3%)
c (c) Relative normalization (1.6%)                  anterior   (4%)
c (d) NC contamination (20%)                         anterior (50%)
c (e) Relative hadronic energy (2.2%)              anterior  (3.3%)
c (f)' cross section  <10GeVÞ                           anterior (12%)
c (g) Beam flux   free
c (h) Neutrino-antineutrino separation
c (i) Partially reconstructed events
       
       sigma_minos(1)=0.103                       !anterior    a) Absolute hadronic E scale (10:3%) 
       sigma_minos(2)=0.022                        !anterior    b) Relative hadronic E scale (3:3%)
       sigma_minos(3)=0.016                          !anterior    c) Normalization (4%)
       sigma_minos(4)=0.20                            !anterior    d) NC contamination (50%)
       sigma_minos(5)=0.03                            !anterior    e) " muon momentum (range 2%, curvature 3%)
       sigma_minos(6)=0.12                            !anterior    f ) cross section   < 10 GeVÞ (12%)
       


	do i=1,7
	do j=1,7
	if (i.eq.j) then
	invmat(i,j)=1.0d0/sigma_minos(i)**2
	else
	invmat(i,j)=0.0d0
	endif
	enddo
	enddo 
	invmat(7,7)=0.0d0



	do i=1,6
	do n=1,39	
	fik_mat(n,i)= sigma_minos(i)/min_th(n) 
	fik_mat(n,7)= 1.0d0

 	if (n.gt.27) then
 	fik_mat(n,6)=0.0	
 	endif 


	min_th_1(n)=min_th(n)
	min_data_1(n)=min_data(n)
! 	print*,min_th_1(n),min_data_1(n)
	
	
	enddo
	Xpull(i)=0.00d0
	Xpull(7)=0.00d0

	enddo

!  	fik_mat(1,7)= 1.0d0
!  	fik_mat(2,7)= 1.0d0
!  	fik_mat(3,7)= 1.0d0
!    	fik_mat(4,7)= 1.0d0
!  	fik_mat(5,7)= 1.0d0
!  	fik_mat(6,7)= 1.0d0





        call ZSPOW(FCN_minos_20,5,7,200,PAR,Xpull,FNORM,WK)


!   	write(*,*) (Xpull(n),n=1,7)

	N_th_tot_bar=0.0
	N_data_tot=0.0d0

	do n=1,39
	sumcoef=0.0d0

	do j=1,7
	sumcoef=sumcoef+Xpull(j)*(fik_mat(n,j))
	enddo !j=1,39
	
	N_bin_bar(n)=min_th(n)*(1.0d0+sumcoef)
	N_th_tot_bar=N_th_tot_bar+N_bin_bar(n)
	N_data_tot=N_data_tot+ min_data(n)	

!    	print*,N_bin_bar(n),min_th(n),min_data(n)

	enddo !n=1,8

!   	print*,N_th_tot_bar,N_data_tot


	chi2_shape=0.0d0
 	do n=1,39

	if (min_data(n).ne.0.0d0) then
	
 	chi2_shape=chi2_shape  + 
     c	2.0d0*( N_bin_bar(n)- min_data(n) -
     c	min_data(n)*log(N_bin_bar(n)/min_data(n) ) )
	
	endif

 	enddo !n=1,8


	sumcoef2=0.0d0
	do i=1,7
	do j=1,7
	sumcoef2=sumcoef2+Xpull(j)*invmat(i,j)*Xpull(i)
	enddo
	enddo
	
! 	print*,chi2_shape,sumcoef2


 	min_chi2_pull=chi2_shape+sumcoef2

! 	print*, min_chi2_pull,Xpull(7)


	end 



	SUBROUTINE FCN_minos_20(Xpull,F,M,PAR)
!========================================================
!=======================================================
!this subroutine tests the set of non-linear equations	
!that have been solved by linearization F(X)=0.0
!in case of not being relatively close to zero 
!this same  subroutine is used as a external function to be solved
!iteratively by ZSPOW
!REAL*8 X(N),F(N),PAR(1)
!F(1)=
!F(N)=
!RETURN
!END
!GIVEN X(1)...X(N), FCN MUST EVALUATE THE
!FUNCTIONS F(1)...F(N) WHICH ARE TO BE MADE
!ZERO. X SHOULD NOT BE ALTERED BY FCN.
!=========================================================
!======================================================
        IMPLICIT NONE

 	include 'MINOSPULL.inc'

	integer j,k,n,i
	real*8 sumcoef,sumcoef2
         REAL*8 Xpull(7),F(7),PAR(1)
	integer M
	real*8 A,D,N_data_tot
	real*8 N_th_tot_bar, N_bin_bar(40)




	N_th_tot_bar=0.0d0
	N_data_tot=0.0
	do n=1,39
	sumcoef=0.0d0
	do j=1,7
	sumcoef=sumcoef+Xpull(j)*(fik_mat(n,j))
	enddo ! j=1,39

	N_bin_bar(n)=min_th_1(n)*(1.0d0+sumcoef)

	N_th_tot_bar=N_th_tot_bar+N_bin_bar(n)
	N_data_tot=N_data_tot+min_data_1(n)
	
	enddo !n=1,8

	

	do j=1,7
	A=0.0d0

	do n=1,39
  	A=A+min_th_1(n)*fik_mat(n,j)*
     c	(1.0d0-(min_data_1(n)/N_bin_bar(n)))

!  	A=A+ fik_mat(n,j)*
!      c	(1.0d0-(min_data_1(n)/N_bin_bar(n)))



	enddo !n=1,39



	D=0.0d0
	do i=1,7
	D=D+invmat(i,j)*Xpull(i)
! 	print*,Xpull(i),'X'
	enddo !i=1,39

	F(j)= A+D


	enddo !	 j=1,3


	end




!=============================================================================
!==========================================================================
!  shape + norm analysis
! 

	subroutine chis2_min_2020_SN(min_th,min_nosc,
     c		min_data,min_chi2_pull)
	implicit none
!this subroutine uses the shifted evets given by 
! N_bin_bar= evbinOSC12*(1+sum(i))
! where sum = Xpull(j)*errmatk2k12(n,j)

 	
 	include 'MINOSPULL.inc'

        real*8  min_th(40),min_nosc(40),min_data(40),min_chi2_pull

        real*8 Xpull(7),F(7),PAR(1),FNORM,WK(5000) 
        real*8 suma,sumpull, sumatix1,sumatix2	
        real*8 sumcoef,sumcoef2
        real*8 N_th_tot_bar, N_bin_bar(40),N_data_tot
        real*8 mat(7,7),sigma_minos(7)
        real*8 chi2_shape,chi2_norm,AA,chi2_i(40)
        real*8 turn_off(40)
	integer i,j,k,h,n

	external FCN_minos_20_SN
	
       
       
! errores systematicos obtenidos de sep 2008  PRL 101, 131802 (2008)
c (a) Absolute hadronic E scale (10:3%) 
c (b) Relative hadronic E scale (3:3%)
c (c) Normalization (4%)
c (d) NC contamination (50%)
c (e) " muon momentum (range 2%, curvature 3%)
c (f ) cross section   < 10 GeVÞ (12%)
c (g) Beam flux    free

c        sigma_minos(1)=0.103                !  a) Absolute hadronic E scale (10:3%) 
c        sigma_minos(2)=0.033                !  b) Relative hadronic E scale (3:3%)
c        sigma_minos(3)=0.04                 !  c) Normalization (4%)
c        sigma_minos(4)=0.50                 !  d) NC contamination (50%)
c        sigma_minos(5)=0.03                 !  e) " muon momentum (range 2%, curvature 3%)
c        sigma_minos(6)=0.12                 !  f ) cross section   < 10 GeVÞ (12%)
       
 ! errores systematicos en analisis de 2020 estan en PHYSICAL REVIEW LETTERS week ending 
 !PRL 106, 181801 (2011)
 
c (a) Hadronic energy                                      anterior    10.3%
c (b) "muon  energy (range 2%, curv. 3%)     anterior  (range 2%, curvature 3%)
c (c) Relative normalization (1.6%)                  anterior   (4%)
c (d) NC contamination (20%)                         anterior (50%)
c (e) Relative hadronic energy (2.2%)              anterior  (3.3%)
c (f)' cross section  <10GeVÞ                           anterior (12%)
c (g) Beam flux   free
c (h) Neutrino-antineutrino separation
c (i) Partially reconstructed events

        turn_off(:)=1.0d0
       
       sigma_minos(1)=   0.103          !anterior    a) Absolute hadronic E scale (10:3%) 
       sigma_minos(2)=   0.022           !anterior    b) Relative hadronic E scale (3:3%)
       sigma_minos(3)=   0.016             !anterior    c) Normalization (4%)
       sigma_minos(4)=   0.50               !anterior    d) NC contamination (50%)
       sigma_minos(5)=   0.03               !anterior    e) " muon momentum (range 2%, curvature 3%)
       sigma_minos(6)=   0.12               !anterior    f ) cross section   < 10 GeVÞ (12%)
       
       	
	

	
          do i=1,7
          do j=1,7
          if (i.eq.j) then
          invmat(i,j)=1.0d0/sigma_minos(i)**2
          else
          invmat(i,j)=0.0d0
          endif
          enddo
          enddo 
          
          invmat(7,7)=1.0d0


       do i=1,6
       do n=1,39	
       fik_mat(n,i)= sigma_minos(i)/min_th(n) 
       fik_mat(n,7)= 0.0d0
       
       
       if (n.gt.27) then
       fik_mat(n,6)=0.0	
       endif 
c        print*,fik_mat(n,i)
       
       min_th_1(n)=min_th(n)
       min_data_1(n)=min_data(n)
c        print*,min_th_1(n),min_data_1(n)
       
       
       enddo
       Xpull(i)=0.00d0
       Xpull(7)=0.00d0
       
       enddo



       call ZSPOW(FCN_minos_20_SN,5,7,200,PAR,Xpull,FNORM,WK)


c        	write(*,*) (Xpull(n),n=1,7)


        N_th_tot_bar=0.0
        N_data_tot=0.0d0
        
        do n=1,39
        
        sumcoef=0.0d0
        
        do j=1,7
        sumcoef=sumcoef+fik_mat(n,j)*Xpull(j)
        enddo !j=1,39
        
        
        N_bin_bar(n)=min_th(n)*(1.0d0+sumcoef)
        N_th_tot_bar=N_th_tot_bar+N_bin_bar(n)
        N_data_tot=N_data_tot+ min_data(n)	
        
        
        enddo !n=1,8

        turn_off(3)=1.0d0
        turn_off(4)=1.0d0 
        turn_off(5)=1.0d0 
c          turn_off(6)=.0d0 
c          turn_off(7)=0.0d0 
c         turn_off(n)=0.0d0         

       chi2_shape=0.0d0
       do n=1,39
       
       if (min_data(n).ne.0.0d0) then
       
          chi2_i(n)=2.0d0  * ( N_bin_bar(n)- min_data(n) -
     c	min_data(n)*log(N_bin_bar(n)/min_data(n) ) )
       
          chi2_shape=chi2_shape  +  chi2_i(n)*turn_off(n)
       
c             print*,N_bin_bar(n),min_data(n),chi2_i(n)
       
       endif
       enddo !n=1,8


       chi2_norm=2.0d0*(N_th_tot_bar-N_data_tot-
     c	     (N_data_tot*log(N_th_tot_bar/N_data_tot)) )

c        print*, 'N_th_tot_bar,N_data_tot,chi2_norm'     
c        print*, N_th_tot_bar,N_data_tot,chi2_norm

         AA=0.0d0
         do i=1,7
         do j=1,7
         AA=AA+ 
     c	 Xpull(j)*invmat(i,j)*Xpull(i)
         
         enddo
         enddo




       min_chi2_pull=chi2_shape+AA+chi2_norm



        end 



	SUBROUTINE FCN_minos_20_SN(Xpull,F,M,PAR)
!========================================================
!=======================================================
!this subroutine tests the set of non-linear equations	
!that have been solved by linearization F(X)=0.0
!in case of not being relatively close to zero 
!this same  subroutine is used as a external function to be solved
!iteratively by ZSPOW
!REAL*8 X(N),F(N),PAR(1)
!F(1)=
!F(N)=
!RETURN
!END
!GIVEN X(1)...X(N), FCN MUST EVALUATE THE
!FUNCTIONS F(1)...F(N) WHICH ARE TO BE MADE
!ZERO. X SHOULD NOT BE ALTERED BY FCN.
!=========================================================
!======================================================
        IMPLICIT NONE

 	include 'MINOSPULL.inc'

	integer j,k,n,i
	real*8 sumcoef,sumcoef2
         REAL*8 Xpull(7),F(7),PAR(1)
	integer M
	real*8 A,D,N_data_tot
	real*8 N_th_tot_bar, N_bin_bar(40)




	N_th_tot_bar=0.0d0
	N_data_tot=0.0

	do n=1,39
	sumcoef=0.0d0
	do j=1,7
	sumcoef=sumcoef+Xpull(j)*(fik_mat(n,j))
	enddo ! j=1,39

	N_bin_bar(n)=min_th_1(n)*(1.0d0+sumcoef)

	N_th_tot_bar=N_th_tot_bar+N_bin_bar(n)
	N_data_tot=N_data_tot+min_data_1(n)
	
	enddo !n=1,8

	

	do j=1,7
	A=0.0d0

	do n=1,39

  	A=A+fik_mat(n,j)*min_th_1(n)*
     c	(2.0d0-(min_data_1(n)/N_bin_bar(n))-
     c	(N_data_tot/N_th_tot_bar ))

!  	A=A+ fik_mat(n,j)*
!      c	(1.0d0-(min_data_1(n)/N_bin_bar(n)))

	enddo !n=1,39


	D=0.0d0
	do i=1,7
	D=D+invmat(i,j)*Xpull(i)
! 	print*,Xpull(i),'X'
	enddo !i=1,39

	F(j)= A+D
	enddo !	 j=1,3

	end

	
