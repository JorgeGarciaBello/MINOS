program main_grid

    implicit none
    integer  i,j,points,k
    parameter (points=6)
    integer, parameter :: n_bins=44800
    !integer, parameter :: n_bins=851200    
    real*8 reno_chi_min,db_chi_min
    real*8 Y(13)                                          !Y=( t12 , t13 , t14 , t23 , t24 , t34 , d13 , d24 , d34 , dm21 , dm31 , dm41 )            
    real*8 var_dm23max, var_the13max, var_the23max
    real*8 var_the23min, var_the13min, var_dm23min
    real*8 sumth13, sumth23, sumdm23, delta_th13, delta_th23, delta_dm23
     
    real*8 var_th13(points), var_th23(points),var_dm23(points)
    real*8 chi2_grid(points,points,points)
    real*8 grid(n_bins,3)
    integer :: u     
    CHARACTER(30)  names    

    names='3D.ih'
    !names='prueba4D'        
    call db_read_data()      ! Lee datos de Dayabay
    call readRENOData()    ! Lee datos de RENO    

    !######################################
    !
    !   Grid de sk-nh con Delta = 0
    !
    !######################################
    !
     open(newunit=u, file='data/sk-nh.dat')
       do i=1,n_bins
         read(u,*) grid(i,:)
       enddo
     close(u)
    
    !######################################
    !open(newunit=u, file='data/sk-complete-nh.dat')
    !  do i=1,n_bins
    !    read(u,*) grid(i,:)
    !  enddo
    !close(u)

    !grid =>  Delta m^2_{32}    sin^2 theta_{13}  sin^2 theta_{23}

    
    Y(1)=7.650d-5 ! dm2_12   
    !Y(2)=2.5d-3 !dm2_23
    Y(3)=0.0d0 !dm2_34
    Y(4)=0.584d0 !theta_12
    !Y(5)=0.1  !theta_13
    !Y(6)=0.78d0  !theta_23
    Y(7)=0.0d0 !theta_14
    Y(8)=0.0d0 !theta_24
    Y(9)=0.0d0 !theta_34
    Y(10)=0.0
    Y(11)=0.0       
    Y(12)=0.0       
    Y(13)=0.0

    
   
    select case(1)
        case(1)
    !########################################################################################
    !
    !   Ieración para un grid de 3D = [dm32,s2t13,s2t23,]
    !
    !########################################################################################
            print*, 'ih'
            print*, 'RENO'
            open(newunit=u,file='reno.grid.dm32.s2t13.s2t23.of.SK.'//trim(names)//'.dat')
              do i=1,n_bins   !s44800
                Y(2)  = -grid(i,1)             !dm32
                Y(5)  = asin(sqrt(grid(i,2))) !t13
                Y(6)  = asin(sqrt(grid(i,3))) !t23                

                call renoChi2(Y,reno_chi_min)
                write(u,*) grid(i,1), grid(i,2), grid(i,3), reno_chi_min                
                print*, i
              enddo
            close(u)
            print*, 'ih'
            print*, 'Daya Bay '
            open(newunit=u,file='dbreno.grid.dm32.s2t13.s2t23.of.SK.'//trim(names)//'.dat')
              do i=1,n_bins !44800
                Y(2)  = -grid(i,1)             !dm32
                Y(5)  = asin(sqrt(grid(i,2))) !t13
                Y(6)  = asin(sqrt(grid(i,3))) !t23               

                call daya_bay_cov(Y,db_chi_min)
                write(u,*) grid(i,1), grid(i,2), grid(i,3), db_chi_min                
                print*, 'DB: ', i
              enddo
            close(u) 
        case(2)
    !########################################################################################
    !
    !   Ieración para un grid de 4D = [delta_CP,dm32,s2t13,s2t23,]
    !
    !########################################################################################
            print*, 'RENO'
            open(newunit=u,file='reno.full.sk.grid.delta.dm32.s2t13.s2t23.'//trim(names)//'.dat')
              do i=1,n_bins 
                Y(2)  = grid(i,2)             !dm32
                Y(5)  = asin(sqrt(grid(i,3))) !t13
                Y(6)  = asin(sqrt(grid(i,4))) !t23
                Y(10) = grid(i,1)

                call renoChi2(Y,reno_chi_min)                
                write(u,*) grid(i,1), grid(i,2), grid(i,3), grid(i,4), reno_chi_min
                print*, i
              enddo
            close(u)

            print*, 'Daya Bay '
            open(newunit=u,file='db.full.sk.grid.delta.dm32.s2t13.s2t23.'//trim(names)//'.dat')
              do i=1,n_bins
                Y(2)  = grid(i,2)             !dm32
                Y(5)  = asin(sqrt(grid(i,3))) !t13
                Y(6)  = asin(sqrt(grid(i,4))) !t23
                Y(10) = grid(i,1)

                call daya_bay_cov(Y,db_chi_min)                
                write(u,*) grid(i,1), grid(i,2), grid(i,3), grid(i,4), db_chi_min
                print*, 'DB: ', i
              enddo
            close(u)
    end select        
end program 
