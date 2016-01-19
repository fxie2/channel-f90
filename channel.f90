!============================================!
!                                            !
!     Direct Numerical Simulation (DNS)      !
!       of a turbulent channel flow          !
!                                            !
!============================================!
! 
! Author: M.Sc. Davide Gatti
! Date  : 28/Jul/2015
! 

PROGRAM channel

  USE dnsdata
 
  WRITE(*,*) " "
  WRITE(*,*) "!====================================================!"
  WRITE(*,*) "!                     D   N   S                      !"
  WRITE(*,*) "!====================================================!"
  WRITE(*,*) " "

  CALL init_fft(rVd,Vd,rVVd,VVd,nxd,nzd)
  CALL setup_derivatives()
  CALL setup_boundary_conditions()
  CALL read_restart_file(V)

  ! Output DNS.in
  WRITE(*,"(A,I5,A,I5,A,I5)") "   nx =",nx,"   ny =",ny,"   nz =",nz
  WRITE(*,"(A,I5,A,I5)") "   nxd =",nxd,"  nzd =",nzd
  WRITE(*,"(A,F6.4,A,F6.4,A,F8.6)") "   alfa0 =",alfa0,"       beta0 =",beta0,"   ni =",ni
  WRITE(*,"(A,F6.4,A,F6.4)") "   meanpx =",meanpx,"      meanpz =",meanpz
  WRITE(*,"(A,F6.4,A,F6.4)") "   meanflowx =",meanflowx, "   meanflowz =", meanflowx
  WRITE(*,*) " "

  !Time loop
  WRITE(*,*) time,sum(d140(-2:2)*dreal(V(0,0,-1:3)%u)),-sum(d14n(-2:2)*dreal(V(0,0,ny-3:ny+1)%u))
  timeloop: DO WHILE (time<t_max-deltat/2.0) 
    time=time+2.0/RK1_rai_coeff*deltat
    CALL buildrhs(RK1_rai); CALL linsolve(RK1_rai_coeff/deltat)
    time=time+2.0/RK2_rai_coeff*deltat
    CALL buildrhs(RK2_rai); CALL linsolve(RK2_rai_coeff/deltat)
    time=time+2.0/RK3_rai_coeff*deltat
    CALL buildrhs(RK3_rai); CALL linsolve(RK3_rai_coeff/deltat)
    !Outstats
    WRITE(*,*) time,sum(d140(-2:2)*dreal(V(0,0,-1:3)%u)),-sum(d14n(-2:2)*dreal(V(0,0,ny-3:ny+1)%u))
    !Write output field
    IF (FLOOR((time+0.5*deltat)/dt_save) > FLOOR((time-0.5*deltat)/dt_save)) THEN
      WRITE(*,*) "Writing Dati.cart.out at time ", time
      CALL save_restart_file(V)
    END IF
  END DO timeloop

  CALL free_fft()

END PROGRAM channel
