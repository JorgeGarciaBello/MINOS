


       subroutine  minos_2020( Y,min_chi2)
       implicit none
       include 'MINOS_20.inc'       
       !	include 'UMAT4.inc'
       !	include 'GLOB.inc'
       	
       	
       integer i,j,k
       real*8 Y(13)
       real*8  dm223,theta13,theta23,norm
       real*8  w12,w13,w14,w23,w24,w34
       real*8  dm212,dm234	
       real*8 theta12,theta14,theta24,theta34
       real*8 min_tot_th,min_tot_nosc,min_tot_data       
       real*8 pi,epsilon1, Ener
       real*8 int_min_data,int_min_nosc, int_min_th
       real*8 P_aver,LoE,Pmumu
       	real*8 min_th(40),min_nosc(40),min_data(40)
        real*8 min_chi2_pull,min_chi2,Prob_k(40)
        real*8 NOS_Minos(40),err_th(40),Ener_cent
        real*8 min_osc(40),min_mc(40),min_bf(40),
     c sigma(40),Enu(40),bin_min(40),bin_max(40),
     c mc_tot, osc_tot
        real*8 sigma_k, sigma_minos(6)
        real*8 res

c i -> MinData20(i,2,40)      
c i=1 ->MINOS-POT-numu-minos-mas.dat
c i=2 ->MINOS-POT-numu-numubar.dat
c i=3 ->MINOS-best-fit.dat
c i=4 ->MINOS-bines.dat
c i=5 ->MINOS-bottom_sigma_limits.dat
c i=6 ->MINOS-data-points.dat
c i=7 ->MINOS-no-oscillation.dat
c i=8 -> MINOS-upper-sigma-limits.dat

       sigma_minos(1)=   0.103          !anterior    a) Absolute hadronic E scale (10:3%) 
       sigma_minos(2)=   0.022           !anterior    b) Relative hadronic E scale (3:3%)
       sigma_minos(3)=   0.016             !anterior    c) Normalization (4%)
       sigma_minos(4)=   0.50               !anterior    d) NC contamination (50%)
       sigma_minos(5)=   0.03               !anterior    e) " muon momentum (range 2%, curvature 3%)
       sigma_minos(6)=   0.12               !anterior    f ) cross section   < 10 GeVÞ (12%)

        do i=1,39
        min_bf(i)=MinData20(3,2,i) 
        bin_min(i)=MinData20(4,1,i) 
        bin_max(i)=MinData20(4,2,i) 
        
        min_osc(i)=MinData20(6,2,i) 
        min_mc(i)=MinData20(7,2,i)  
        sigma(i)= abs(MinData20(8,2,i)-MinData20(5,2,i))/2.0
        sigma(i)=sigma(i)/min_osc(i)        
c        write(*,*) sigma(i),sigma(i)/min_osc(i),i    
         enddo
         pi=3.14159265    
         
       dm212=  7.54d-5
       dm223=  2.480d-3
       dm234=  0.0d0
       theta12= 0.5873
       theta13= 0.1454
       theta23= 0.6642
       theta14= 0.0d0
       theta24= 0.0d0
       theta34= 0.0d0
! 	  norm=Y(10)
 	    norm=1.0d0 

c       call umat(theta12,theta13,theta23,
c     c		theta14,theta24,theta34)
         
         epsilon1=0.005d0

          do k=1,39
         Ener=bin_min(k)         
         P_aver=0.0d0 
 	   do  while (Ener.le.bin_max(k))         
 	   
 	    LoE=735.0d0/Ener
c           w12=1.27d0*dm212*LoE 
c           w23=1.27d0*dm223*LoE 
c             w34=1.27d0*dm234*LoE 
c             w13=w12+w23      
c             w14=w13+w34      
c             w24=w23+w34      
                       
c            Pmumu= 1.0d0  -4.d0*       
c     c        (cmm(1)*(sin(w12)**2)      
c     c       +cmm(2)*(sin(w13)**2)      
c     c       +cmm(3)*(sin(w23)**2)      
c     c       +cmm(4)*(sin(w14)**2)       
c     c       +cmm(5)*(sin(w24)**2)             
c     c       +cmm(6)*(sin(w34)**2) )   
       call minosoneslab(2,2,
     c theta12,theta23,theta13,0.0d0,dm212,dm223,Ener,1,res)
       

         Pmumu=res
         P_aver=P_aver+Pmumu*(epsilon1)  
         Ener=Ener+epsilon1
  	   
  	   enddo
     
        P_aver=P_aver/(bin_max(k)-bin_min(k))
        
        Prob_k(k)=P_aver
        
         NOS_Minos(k)=  min_bf(k)/Prob_k(k)    ! min_osc(k)/Prob_k(k)       
          err_th(k)=abs(NOS_Minos(k)-min_mc(k)) /min_mc(k)
         enddo       
       !################################################
       !
       !  Calibrations to fix the confidence region       !
       !
       !################################################       
       do k=1,5
        NOS_Minos(k)=NOS_Minos(k)*0.9d0  ! FIxed to 0.7
       enddo
       do k=6,10
        NOS_Minos(k)=NOS_Minos(k)*0.8d0  ! FIxed to 1.2 !0.6,0.8 casi funciona, manditne minimo y valores de masa y angulo, pero neceitamosaplanar la curva
       enddo
       do k=11,15
        NOS_Minos(k)=NOS_Minos(k)*1.07d0  ! Fixed 1.07
       enddo

       do k=21,25
        NOS_Minos(k)=NOS_Minos(k)*0.98d0  ! Fixed 0.98 ! 0.89
       enddo
       do k=26,30
        NOS_Minos(k)=NOS_Minos(k)*0.98d0  ! Fixed 0.98 ! 0.89
       enddo
       do k=31,35
        NOS_Minos(k)=NOS_Minos(k)*0.98d0  ! Fixed 0.98 ! 0.89
       enddo
       do k=36,39
        NOS_Minos(k)=NOS_Minos(k)*0.98d0  ! Fixed 0.98 ! 0.89
       enddo
       !################################################
              
       dm212=Y(1)
       dm223=1.0d0*Y(2)
       dm234=Y(3)
       theta12=Y(4)
       theta13=Y(5)
       theta23=Y(6)
       theta14=Y(7)
       theta24=Y(8)
       theta34=Y(9)
c! 	norm=Y(10)
 	    norm=1.0d0

c
c       call umat(theta12,theta13,theta23,
c     c		theta14,theta24,theta34)               
     
     
     
        epsilon1=0.005d0

       min_tot_th=0.0d0
       min_tot_nosc=0.0d0
       min_tot_data=0.0d0        
       
        

         do k=1,39

        Ener=bin_min(k)
        int_min_th=0.0d0
        int_min_nosc=0.0d0
        int_min_data=0.0d0
        P_aver=0.0d0
        
          
           if  (  k.eq.12 ) then 
           Ener= 0.999* Ener
           endif 
           if  (  k.eq.11) then 
           Ener= 0.996* Ener
           endif 
           if  (  k.eq.10 ) then 
           Ener= 0.995* Ener
           endif 
          if  (  k.eq.9) then 
           Ener= 0.99* Ener
           endif          
           if  (  k.eq.8) then 
           Ener= 0.98* Ener
           endif                
          if  (  k.eq.7) then 
           Ener= 0.97* Ener
           endif           
           if  (  k.eq.6) then 
           Ener= 0.93 * Ener
           endif                  
           if  (  k.eq.5) then 
           Ener= 0.80 * Ener
           endif                     
           if  (  k.eq.4) then 
           Ener= 0.78 * Ener
           endif                              
           if  (  k.eq.3) then 
           Ener= 0.90* Ener
           endif                              
      
        
        Ener_cent=bin_min(k) + ((bin_max(k)-bin_min(k))/2.0d0) 	   

	      LoE=735.0d0/Ener_cent
	    
 
c          w12=1.27d0*dm212*LoE 
c          w23=1.27d0*dm223*LoE 
c            w34=1.27d0*dm234*LoE 
c            w13=w12+w23      
c            w14=w13+w34      
c            w24=w23+w34              
       !####################################################
       !##########                       ###################
       !####################################################    
       call minosoneslab(2,2,
     c theta12,theta23,theta13,0.0d0,dm212,dm223,Ener_cent,1,res)
       Prob_k(k) = res
c               Prob_k(k) = 1.0d0  -4.d0*       
c     c        (cmm(1)*(sin(w12)**2)      
c     c       +cmm(2)*(sin(w13)**2)      
c     c       +cmm(3)*(sin(w23)**2)      
c     c       +cmm(4)*(sin(w14)**2)        
c     c       +cmm(5)*(sin(w24)**2)             
c     c       +cmm(6)*(sin(w34)**2) )    
       !####################################################s           
            
        
        
      
	   do  while (Ener.le.bin_max(k))      
      
	    LoE=735.0d0/Ener     
 
c          w12=1.27d0*dm212*LoE 
c          w23=1.27d0*dm223*LoE 
c            w34=1.27d0*dm234*LoE 
c            w13=w12+w23      
c            w14=w13+w34      
c            w24=w23+w34      
           
           
c           Pmumu= 1.0d0  -4.d0*       
c     c        (cmm(1)*(sin(w12)**2)      
c     c       +cmm(2)*(sin(w13)**2)      
c     c       +cmm(3)*(sin(w23)**2)      
c     c       +cmm(4)*(sin(w14)**2)        
c     c       +cmm(5)*(sin(w24)**2)             
c     c       +cmm(6)*(sin(w34)**2) )    

       !####################################################
       !##########                       ###################
       !####################################################
       call minosoneslab(2,2,
     c theta12,theta23,theta13,0.0d0,dm212,dm223,Ener,1,res)
       Pmumu=res
       !####################################################
        int_min_data=int_min_data + min_osc(k)*(epsilon1)
        int_min_th= int_min_th +  NOS_Minos(k)*Pmumu*(epsilon1)  
        int_min_nosc= int_min_nosc +  NOS_Minos(k)*(epsilon1)        
        
        P_aver=P_aver+Pmumu*(epsilon1)  
        Ener=Ener+epsilon1

 	   enddo
       P_aver=P_aver/(bin_max(k)-bin_min(k))
      
       
              min_th(k)=int_min_th          !  / (bin_max(k)-bin_min(k))
              min_nosc(k)=int_min_nosc  ! /(bin_max(k)-bin_min(k))
              min_data(k)=int_min_data   !/(bin_max(k)-bin_min(k))
              Prob_k(k)=P_aver
         
         
         min_tot_th=min_tot_th+min_th(k)
         min_tot_nosc=min_tot_nosc+min_nosc(k)
         min_tot_data=min_tot_data+min_data(k)
        enddo 
        
       min_chi2_pull=0.0d0
       do k=1,39                  
          sigma_k= min_data(k) ! * sigma(k)  !min_th(k)  !   !*   
          if(k<=10) sigma_k=sigma_k*13.0   ! FIxed value 10.0 && 12.0 en 20
          
          min_chi2_pull=min_chi2_pull+
     c    (( min_th(k)*0.99d0 -  min_data(k))  )**2   / (sigma_k)
       enddo
        min_chi2=min_chi2_pull
            print*, min_chi2       
         
       end 
        

