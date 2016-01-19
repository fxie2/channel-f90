!============================================!
!                                            !
!           Fast Fourier Transforms          !
!                  for the                   !
!      Direct Numerical Simulation (DNS)     !
!        of a turbulent channel flow         !
!                                            !
!============================================!
! 
! Author: M.Sc. Davide Gatti
! Date  : 28/Jul/2015
! 

MODULE ffts

  USE, intrinsic :: iso_c_binding
  USE typedef
  IMPLICIT NONE
  INCLUDE 'fftw3.f03'

  integer, save        :: plan_type=FFTW_ESTIMATE
  TYPE(C_PTR), save    :: pFFT,pIFT,pRFT,pHFT,pFFT6,pIFT3,ptrVd,ptrVVd

CONTAINS

  !--------------------------------------------------------------!
  !------------------ init aligned FFT vectors ------------------!
  SUBROUTINE init_fft(rVd,Vd,rVVd,VVd,nxd,nzd)
    integer(C_INT), intent(in) :: nxd,nzd
    complex(C_DOUBLE_COMPLEX), pointer, intent(out) :: Vd(:,:,:)
    complex(C_DOUBLE_COMPLEX), pointer,  intent(out) :: VVd(:,:,:)
    real(C_DOUBLE), pointer, intent(out) :: rVd(:,:,:)
    real(C_DOUBLE), pointer, intent(out) :: rVVd(:,:,:)
    integer(C_INT), dimension(1) :: n_z, n_x, rn_x
    integer(C_INT), dimension(2) :: n
    n_z=[nzd]; n_x=[nxd]; rn_x=[2*nxd]; n=[nzd, 2*nxd]
    ptrVd =fftw_alloc_complex(int((nxd+1)*nzd*3, C_SIZE_T))
    ptrVVd=fftw_alloc_complex(int((nxd+1)*nzd*6, C_SIZE_T))
    CALL c_f_pointer(ptrVd, Vd, [nxd+1,nzd,3]);      CALL c_f_pointer(ptrVVd, VVd, [nxd+1,nzd,6])
    CALL c_f_pointer(ptrVd, rVd, [2*(nxd+1),nzd,3]); CALL c_f_pointer(ptrVVd, rVVd, [2*(nxd+1),nzd,6])    
    pFFT=fftw_plan_many_dft(1, n_z, 1, Vd(1,1,:), n_z, 1, 0, Vd(1,1,:), n_z, 1, 0, FFTW_FORWARD,  plan_type)
    pIFT=fftw_plan_many_dft(1, n_z, 1, Vd(1,1,:), n_z, 1, 0, Vd(1,1,:), n_z, 1, 0, FFTW_BACKWARD, plan_type)
    pRFT=fftw_plan_many_dft_c2r(1, rn_x, 1, Vd(1,:,1),  n_x+1,  1, 0, rVd(:,1,1), rn_x+2, 1, 0,   plan_type)
    pHFT=fftw_plan_many_dft_r2c(1, rn_x, 1, rVd(:,1,1), rn_x+2, 1, 0, Vd(:,1,1),     n_x+1, 1, 0, plan_type)
    pFFT6=fftw_plan_many_dft_r2c(2, n, 6, rVd, [nzd,2*(nxd+1)], 1, 2*(nxd+1)*nzd,  Vd, [nzd,nxd+1], 1,   (nxd+1)*nzd, plan_type)
    pIFT3=fftw_plan_many_dft_c2r(2, n, 3, Vd,  [nzd,nxd+1], 1,   (nxd+1)*nzd, rVd, [nzd,2*(nxd+1)], 1, 2*(nxd+1)*nzd, plan_type)    
  END SUBROUTINE init_fft

  SUBROUTINE FFT(x)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:)
    CALL fftw_execute_dft(pFFT,x,x)
  END SUBROUTINE FFT

  SUBROUTINE IFT(x)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:)
    CALL fftw_execute_dft(pIFT,x,x)
  END SUBROUTINE IFT

  SUBROUTINE RFT(x,rx) !Ritestare
    complex(C_DOUBLE_COMPLEX) :: x(:)
    real(C_DOUBLE) :: rx(:)
    CALL fftw_execute_dft_c2r(pRFT,x,rx)
  END SUBROUTINE RFT

  SUBROUTINE HFT(x,rx) !Ritestare
    complex(C_DOUBLE_COMPLEX) :: x(:)
    real(C_DOUBLE) :: rx(:)
    CALL fftw_execute_dft_r2c(pHFT,rx,x)
  END SUBROUTINE HFT

  SUBROUTINE free_fft()
    CALL fftw_free(ptrVd); CALL fftw_free(ptrVVd);
  END SUBROUTINE free_fft


END MODULE ffts
