!============================================!
!                                            !
!           Fast Fourier Transforms          !
!                  for the                   !
!      Direct Numerical Simulation (DNS)     !
!        of a turbulent channel flow         !
!                                            !
!============================================!
! 
! Author: Dr. Davide Gatti
! Date  : 28/Jul/2015
! 

MODULE ffts

  USE, intrinsic :: iso_c_binding
  USE typedef
  IMPLICIT NONE
  INCLUDE 'fftw3.f03'

  integer, save        :: plan_type=FFTW_PATIENT
  TYPE(C_PTR), save    :: pFFT,pIFT,ptrVVd

CONTAINS

  !--------------------------------------------------------------!
  !------------------ init aligned FFT vectors ------------------!
  SUBROUTINE init_fft(VVd,rVVd,nxd,nzd)
    integer(C_INT), intent(in) :: nxd,nzd
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:), intent(out) :: VVd
    real(C_DOUBLE), pointer, dimension(:,:,:,:), intent(out) :: rVVd
    integer(C_INT), dimension(2) :: n
    n=[nzd,2*nxd];
    !Allocate aligned memory
    ptrVVd=fftw_alloc_complex(int((nxd+1)*nzd*6*5, C_SIZE_T))
    !Convert C to F pointer 
    CALL c_f_pointer(ptrVVd, VVd, [nxd+1,nzd,6,5]);   CALL c_f_pointer(ptrVVd,  rVVd, [2*(nxd+1),nzd,6,5])  
                                                     
    !FFTs plans
    pIFT=fftw_plan_many_dft_c2r(2, n, 3,   VVd, [nzd,nxd+1], 1, (nxd+1)*nzd, &
                                              rVVd, [nzd,2*(nxd+1)], 1, 2*(nxd+1)*nzd, plan_type)
    pFFT=fftw_plan_many_dft_r2c(2, n, 6,  rVVd, [nzd,2*(nxd+1)], 1, 2*(nxd+1)*nzd, &
                                               VVd, [nzd,nxd+1], 1, (nxd+1)*nzd, plan_type)
  END SUBROUTINE init_fft

  LOGICAL FUNCTION fftFIT(i) result(isFIT)
    integer(C_INT), intent(in) :: i
    integer(C_INT) :: j
    j=i
    DO WHILE ( MOD(j,2)==0 )
      j=SHIFTA(j,1)
    END DO
    isFIT=( (j==1) .OR. (j==3) )
    !WRITE(*,*) i,j,isFIT
  END FUNCTION fftFIT

  SUBROUTINE IFT(x,rx) 
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:,:,:)
    real(C_DOUBLE), intent(inout) :: rx(:,:,:)
    CALL fftw_execute_dft_c2r(pIFT,x,rx)
  END SUBROUTINE IFT

  SUBROUTINE FFT(rx,x) 
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:,:,:)
    real(C_DOUBLE),intent(inout) :: rx(:,:,:)
    CALL fftw_execute_dft_r2c(pFFT,rx,x)
  END SUBROUTINE FFT

  SUBROUTINE free_fft()
    CALL fftw_free(ptrVVd);
  END SUBROUTINE free_fft


END MODULE ffts
