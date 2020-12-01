        program minos_grid
        
        implicit none
        integer  i,j,points,k
        !parameter (points=44800 )        
        parameter (points=30 )

        real*8 chi_dayabay, chi_reno, chi_doubleCHOOZ
        real*8 chi_min_global,db_chi_min
        real*8 db_chi_square_spectral_analysis2_period
        real*8 DC_FUNC    
        real*8 reno_chi_square_spectral_analysis
        real*8 Y(13)      !Y=( t12, t13,t14,t23,t24,t34,d13,d24,d34,dm21,dm31,dm41)        
        real*8 delta_dm23
        real*8  var_dm23max,var_the13max,var_the23max
        real*8 var_the23min,var_the13min,var_dm23min
        real*8 sumth13, sumth23, sumdm23,  delta_th13,delta_th23
        real*8 var_th13(points), var_th23(points),var_dm23(points)
        real*8 chi2_grid(points,points,points)
        real*8 min_chi2,pi,sin2_th
        real*8 grid(points,3),minos_ji(points)

     
        CHARACTER(50)  names        
       
        names='3D.nh'        
        
         call Read_MinData20                
         
       var_the23min=0.5
       var_the23max=1.0 

       var_the13min=-0.45
       var_the13max=0.45
 
       var_dm23max=2.70d-3
       var_dm23min=2.20d-3

       !var_dm23max=9.0d-3
       !var_dm23min=0.90d-3
       
       delta_th13=abs( var_the13max - var_the13min) /float(points-1)     
       delta_th23=abs(var_the23max - var_the23min)/float(points-1)
       delta_dm23=abs( var_dm23max - var_dm23min) /float(points-1)     

       
       
       sumth13=var_the13min!initial parameters
       sumth23=var_the23min
       sumdm23=var_dm23min
       	
       do i=1,points
       var_th13(i)=sumth13
       var_th23(i)=sumth23 
       
        var_dm23(i)=sumdm23
       
!        	print*,var_th13(i),var_th23(i),var_dm23(i)
       
       sumth13=sumth13+delta_th13
       sumth23=sumth23+delta_th23
       sumdm23=sumdm23+delta_dm23
       
       
       enddo
       
!c        Y(1)=0.0d0*7.650d-5 ! dm2_12   
!c        !Y(2)=2.5d-3 !dm2_23
!c        Y(3)=0.0d0 !dm2_34
!c        
!c        Y(4)=0.0d0*0.584d0 !theta_12
!c        !Y(5)=0.1  !theta_13
!c       ! Y(6)=0.78  !theta_23
!c        
!c        Y(7)=0.0d0 !theta_14
!c        Y(8)=0.0d0 !theta_24
!c        Y(9)=0.0d0 !theta_34
!c        Y(10)=0.0     
       
          Y(1) =7.54d-5
          Y(2) =2.480d-3
          Y(3) =0.0d0
          Y(4) =0.5873
          Y(5) =0.1454
          Y(6) =0.6642
          Y(7) =0.0d0
          Y(8) =0.0d0
          Y(9) =0.0d0
          Y(10)=0.0
         
!=======================================================       
         pi= 3.14159265
!c          Y(2)=2.5d-3
!c          Y(5)=0.0
!c          Y(6)=pi/4.0d0             
!c          call minos_2020_print( Y,min_chi2)
!c          print*, min_chi2
!=======================================================       
       !######################################
       !
       !   Grid de sk-nh con Delta = 0
       !
       !######################################
       !
!c       open(30, file='MINOS_data/sk-nh.dat')
!c       do i=1,points
!c         read(30,*) grid(i,:)
!c       enddo
!c       close(30)
!cc
!c       names='3D.grid.of.sk.dm32.s2t13.s2t23.nh'
!c       min_chi2=0.0d0
!c       do i=1,points
!c        Y(2)  = grid(i,1)            !dm32
!c        Y(5)  = asin(sqrt(grid(i,2))) !t13
!c        Y(6)  = asin(sqrt(grid(i,3))) !t23
!c        call minos_2020( Y,min_chi2)
!c        minos_ji(i)=min_chi2
!c       enddo
!c
!c       open(30,file='minos.'//names)
!c       do i=1,points
!c        write(30, '(4F20.8)') grid(i,1),
!c     c  grid(i,2), grid(i,3), minos_ji(i)
!c       enddo            
!c       close(30)
      
!c       names='3D.grid.of.sk.dm32.s2t13.s2t23.ih'
!c       min_chi2=0.0d0
!c       do i=1,points
!c        Y(2)  = -grid(i,1)            !dm32
!c        Y(5)  = asin(sqrt(grid(i,2))) !t13
!c        Y(6)  = asin(sqrt(grid(i,3))) !t23
!c        call minos_2020( Y,min_chi2)
!c        minos_ji(i)=min_chi2
!c       enddo
!cc
!c       open(30,file='minos.'//names)
!c       do i=1,points
!c        write(30, '(4F20.8)') grid(i,1),
!c     c  grid(i,2), grid(i,3), minos_ji(i)
!c       enddo            
!c       close(30)

       !#######################################
       !#######################################
       !#######################################

         min_chi2=0.0d0
!$omp parallel do private (Y)
         do i=1,points

!!$omp  do
         do k=1,points
         !do j=1,points         
         Y(2)=var_dm23(i)
         !Y(5)=var_th13(j)
         Y(6)=var_th23(k)
               call minos_2020( Y,min_chi2)
               chi2_grid(i,1,k) = min_chi2 
        enddo
!!$omp  end do

        !enddo
        enddo
!$omp  end parallel do


       !open(30,file='g.dm.s2th13.s2th23.'//names)
       open(30,file='g.dm.s2th23.'//names)
       do k=1,points
       !do j=1,points
       do i=1,points
          sin2_th=sin(var_th23(k))**2

!c          write(30, '(4F20.8)') var_dm23(i),
!c     c  sin(var_th13(j))**2,sin2_th, chi2_grid(i,j,k)

          write(30, *)  sin2_th, var_dm23(i),chi2_grid(i,1,k)
          !write(*, *)  sin2_th, var_dm23(i),chi2_grid(i,1,k)
        enddo
        !enddo
        enddo            
        close(30)

       end