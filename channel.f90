!============================================!
!                                            !
!     Direct Numerical Simulation (DNS)      !
!       of a turbulent channel flow          !
!                                            !
!============================================!
! 
! This program has been written following the
! KISS (Keep it Simple and Stupid) philosophy
!           
! Author: Dr.-Ing. Davide Gatti
! Date  : 28/Jul/2015
! 

PROGRAM channel

  USE dnsdata
  REAL timei,timee 

  ! Read simulation params
  CALL read_dnsin()
  CALL init_memory()

  ! Init various subroutines
  CALL init_fft(VVd,rVVd,nxd,nzd)
  CALL setup_derivatives()
  CALL setup_boundary_conditions()
  CALL read_restart_file()

  ! Output DNS.in
  WRITE(*,*) " "
  WRITE(*,*) "!====================================================!"
  WRITE(*,*) "!                     D   N   S                      !"
  WRITE(*,*) "!====================================================!"
  WRITE(*,*) " "
  WRITE(*,"(A,I5,A,I5,A,I5)") "   nx =",nx,"   ny =",ny,"   nz =",nz
  WRITE(*,"(A,I5,A,I5)") "   nxd =",nxd,"  nzd =",nzd
  WRITE(*,"(A,F6.4,A,F6.4,A,F8.6)") "   alfa0 =",alfa0,"       beta0 =",beta0,"   ni =",ni
  WRITE(*,"(A,F6.4,A,F6.4)") "   meanpx =",meanpx,"      meanpz =",meanpz
  WRITE(*,"(A,F6.4,A,F6.4)") "   meanflowx =",meanflowx, "   meanflowz =", meanflowz
  WRITE(*,*) " "

  ! Compute CFL
  DO iy=0,ny
   CALL convolutions(iy,1,.TRUE.)
  END DO
  ! Time loop
  CALL outstats()
  timeloop: DO WHILE (time<t_max-deltat/2.0) 
    CALL CPU_TIME(timei)
    time=time+2.0/RK1_rai(1)*deltat
    CALL buildrhs(RK1_rai,.TRUE. ); CALL linsolve(RK1_rai(1)/deltat)
    time=time+2.0/RK2_rai(1)*deltat
    CALL buildrhs(RK2_rai,.FALSE.); CALL linsolve(RK2_rai(1)/deltat)
    time=time+2.0/RK3_rai(1)*deltat
    CALL buildrhs(RK3_rai,.FALSE.); CALL linsolve(RK3_rai(1)/deltat)
    CALL outstats()
    CALL CPU_TIME(timee)
    WRITE(*,*) timee-timei
  END DO timeloop

  ! Realease memory
  CALL free_fft()
  CALL free_memory()

END PROGRAM channel
