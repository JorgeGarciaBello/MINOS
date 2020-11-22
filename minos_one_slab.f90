subroutine minosoneslab(flvr1,flvr2,t12,t23,t13,delta,sm,aM,P,nu,result)
 implicit none
 real(8), parameter :: PI= acos(-1.0d0)
 real(8), parameter :: N_A= 6.0221415D23  ! N_A is the Avogadro's number [1/mol]
 integer, parameter :: n=1     ! number of slabs
 real(8) :: t12,t23,t13,delta  ! angle mixing and phase
 real(8) :: sm,aM              ! sm,aM are the squared mass difference m=m_21 y M=m_32
 real(8) :: P                  ! Energy
 integer :: nu                 ! nu is 1 for neutrinos an 2 for antineutrino
 integer :: flvr1              ! flavor intial 1 2 3 = e, m, t
 integer :: flvr2              ! flavor final 1 2 3 = e, m, t
 real(8), dimension(n) :: L    ! L is the length of slabs in serie. unitis [km]
 real(8), dimension(n) :: Ne   ! Ne is the electron density. unitis [mol/cm^{3}]

 real(8) :: probabilityOfTransitionAB    ! fuction
 integer :: k                  ! counter
 real(8) :: result


! real(8),parameter :: scalaFactor=2.5399811853d10  ! scalaFactor is the scala factor to obtein length in [eV^{-1}]
! real(8) :: Pt                                     ! Probabilidad de transición
! real(8) :: eta,rEarth,r1,r2,r3                    ! other variables atmospheric
! real(8) :: matterDensity
! real(8) :: jump

! inner parameters    
! t12=PI / 4.0d0       ! equiv to 45 degrees
! t23=PI *33.0d0 / 180 ! equiv to 45 degrees
! t13=PI * 5.0d0 / 180 ! equiv to 5 degrees
! delta= 0.0d0         ! fase
! sm= 0.0d0            ! [eV^2]
! aM= 3.2d-3           ! [eV^2]
! P=  1.0d-1           ! Energía en [GeV]
! eta=cos(0.0d0)       ! angulo nadir
! rEarth= 6378.0d0     ! Longitud en [Km]
! nu= 1                ! 1 neutrino 2 antineutrino
! r1= 1220.0d0         ! [km]
! r2= 3470             ! [km]
! r3= 6336.0d0         ! [km]

! slab longitude

 L(1)= 735           ! [km] far detector minos for de source

! electron density

 Ne(1)= 1.36d0  ! rho1 Z1/A1 [mol/cm^{3}] acelerator minos
 !Ne(1)= 0.0d0  ! rho1 Z1/A1 [mol/cm^{3}] acelerator minos


! jump=1.0d-2

 

! function
! probabilityOfTransitionAB(flvr1,flvr2,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
! flvr2 and flvr1 = 1,2,3. 1 is flavor electron, 2 flavor muon, 3 flavor tau
! L longitud in km
! tij Mixing angle
! delta fase
! sm, suaqre solar mass diference m_sol ^{2} = m^{2} - m^{1}
! aM, squared atmosferic mass diference m_ATM ^{2} = m^{2} - m^{1}
! P, energy
! nu 1 neutrino, 2 antineutrino
! Ne eledtron density 

 
  
  result=0.0d0
  result=probabilityOfTransitionAB(flvr1,flvr2,L(1),t12,t23,t13,delta,sm,aM,P,nu,Ne(1))

  !probabilityOfTransitionAB(1,2,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)

  !probabilityOfTransitionAB(2,3,L,t12,t23,t13,delta,sm,aM,P,nu,Ne)
 


 return

end subroutine minosoneslab
