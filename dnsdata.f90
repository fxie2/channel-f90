!============================================!
!                                            !
!    Data Structures, Definitions and I/O    !
!                  for the                   !
!     Direct Numerical Simulation (DNS)      !
!        of a turbulent channel flow         !
!                                            !
!============================================!
! 
! Author: M.Sc. Davide Gatti
! Date  : 28/Jul/2015
! 

MODULE dnsdata

  USE, intrinsic :: iso_c_binding
  USE rbmat
  USE ffts
  IMPLICIT NONE

  !Simulation parameters
  real(C_DOUBLE), parameter :: PI=3.1415926535897932384626433832795028841971
  integer(C_INT), parameter :: nx=64, ny=100, nz=64, nxd=3*nx/2, nzd=3*nz
  real(C_DOUBLE), parameter :: alfa0=1.0d0,beta0=2.0d0,ni=1.0d0/4760.0d0,tmax=100.0,a=1.6d0, ymin=0.0d0, ymax=2.0d0
  integer(C_INT), private :: iy
  real(C_DOUBLE) :: cfl=0.0, deltat=1.0d-2,cflmax=0,time=0,dt_field=0,dt_save=10,t_max=5
  real(C_DOUBLE) :: meanpx=0.0d0,meanpz=0.0d0,meanflowx=1.3333d0,meanflowz=0.0d0
  !Grid
  real(C_DOUBLE), dimension(-1:ny+1) :: &
                  y=(/(ymin+0.5d0*(ymax-ymin)*(tanh(a*(2.0d0*real(iy)/real(ny)-1))/tanh(a)+0.5d0*(ymax-ymin)), iy=-1, ny+1)/)
  real(C_DOUBLE), parameter :: dx=PI/(alfa0*nxd), dz=2.0d0*PI/(beta0*nzd),  factor=1.0d0/(2.0d0*nxd*nzd)
  !Derivatives
  TYPE(Di), dimension(1:ny-1) :: der
  real(C_DOUBLE), dimension(-2:2) :: d040,d140,d14m1,d04n,d14n,d24n,d14np1
  real(C_DOUBLE), dimension(1:ny-1,-2:2) :: D0mat, etamat, D2vmat
  !Fourier-transformable arrays (allocated in ffts.f90)
  complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) ::  Vd   
  real(C_DOUBLE), pointer, dimension(:,:,:) :: rVd              
  complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: VVd
  real(C_DOUBLE), pointer, dimension(:,:,:) :: rVVd
  !Solution
  TYPE(RHSTYPE)  :: newrhs(-2:0,0:nx,-nz:nz), oldrhs(-1:ny+1,0:nx,-nz:nz)
  TYPE(MOMFLUX)  :: VV(0:nx,-nz:nz,-2:2)
  TYPE(VELOCITY) :: V(0:nx,-nz:nz,-1:ny+1)
  !Boundary conditions
  real(C_DOUBLE), dimension(-2:2) :: v0bc,v0m1bc,vnbc,vnp1bc,eta0bc,eta0m1bc,etanbc,etanp1bc
  !ODE library coefficients
  real(C_DOUBLE), parameter :: CN_AB_coeff=2.0d0
  real(C_DOUBLE), parameter :: RK1_rai_coeff=120.0d0/32.0d0, RK2_rai_coeff=120.0d0/8.0d0, RK3_rai_coeff=120.0d0/20.0d0

  CONTAINS

  !--------------------------------------------------------------!
  !--------------- Set-up the compact derivatives ---------------!
  SUBROUTINE setup_derivatives()
    real(C_DOUBLE) :: M(0:4,0:4), t(0:4)
    integer(C_INT) :: iy,i,j
    DO iy=1,ny-1
      FORALL (i=0:4, j=0:4) M(i,j)=(y(iy-2+j)-y(iy))**(4.0d0-i); CALL LUdecomp(M)
      t=0; t(0)=24
      der(iy)%d4(-2:2)=M.bs.t
      FORALL (i=0:4, j=0:4) M(i,j)=(5.0d0-i)*(6.0d0-i)*(7.0d0-i)*(8.0d0-i)*(y(iy-2+j)-y(iy))**(4.0d0-i); CALL LUdecomp(M)
      FORALL (i=0:4) t(i)=sum( der(iy)%d4(-2:2)*(y(iy-2:iy+2)-y(iy))**(8.0d0-i) )
      der(iy)%d0(-2:2)=M.bs.t
      FORALL (i=0:4, j=0:4) M(i,j)=(y(iy-2+j)-y(iy))**(4.0d0-i); CALL LUdecomp(M)
      t=0; FORALL (i=0:2) t(i)=sum( der(iy)%d0(-2:2)*(4.0d0-i)*(3.0d0-i)*(y(iy-2:iy+2)-y(iy))**(2.0d0-i) )
      der(iy)%d2(-2:2)=M.bs.t
      t=0; FORALL (i=0:3) t(i)=sum( der(iy)%d0(-2:2)*(4.0d0-i)*(y(iy-2:iy+2)-y(iy))**(3.0d0-i) )
      der(iy)%d1(-2:2)=M.bs.t
    END DO
    FORALL (i=0:4, j=0:4) M(i,j)=(y(-1+j)-y(0))**(4.0d0-i); CALL LUdecomp(M)
    t=0; t(3)=1.0; d140(-2:2)=M.bs.t
    FORALL (i=0:4, j=0:4) M(i,j)=(y(-1+j)-y(-1))**(4.0d0-i); CALL LUdecomp(M)
    t=0; t(3)=1.0; d14m1(-2:2)=M.bs.t
    d04n=0; d04n(1)=1; d040=0; d040(-1)=1
    FORALL (i=0:4, j=0:4) M(i,j)=(y(ny-3+j)-y(ny))**(4.0d0-i); CALL LUdecomp(M)
    t=0; t(3)=1; d14n(-2:2)=M.bs.t
    t=0; t(2)=2; d24n(-2:2)=M.bs.t
    FORALL (i=0:4, j=0:4) M(i,j)=(y(ny-3+j)-y(ny+1))**(4.0d0-i); CALL LUdecomp(M)
    t=0; t(3)=1; d14np1(-2:2)=M.bs.t
    FORALL (iy=1:ny-1) D0mat(iy,-2:2)=der(iy)%d0(-2:2); CALL LU5decomp(D0mat)
  END SUBROUTINE setup_derivatives

  !--------------------------------------------------------------!
  !--------------- Set-up the boundary conditions ---------------!
  SUBROUTINE setup_boundary_conditions()
    v0bc=d040; v0m1bc=d140; eta0bc=d040
    vnbc=d04n; vnp1bc=d14n; etanbc=d04n
    etanp1bc=der(ny-1)%d4
    eta0m1bc=der(1)%d4
    v0bc(-1:2)=v0bc(-1:2)-v0bc(-2)*v0m1bc(-1:2)/v0m1bc(-2)
    eta0bc(-1:2)=eta0bc(-1:2)-eta0bc(-2)*eta0m1bc(-1:2)/eta0m1bc(-2)
    vnbc(-2:1)=vnbc(-2:1)-vnbc(2)*vnp1bc(-2:1)/vnp1bc(2)
    etanbc(-2:1)=etanbc(-2:1)-etanbc(2)*etanp1bc(-2:1)/etanp1bc(2)
  END SUBROUTINE setup_boundary_conditions

  !--------------------------------------------------------------!
  !---------------- integral in the y-direction -----------------!
  PURE FUNCTION yintegr(f) result(II)
    real(C_DOUBLE), intent(in) :: f(-1:ny+1)
    real(C_DOUBLE) :: II, yp1, ym1, a1, a2, a3
    integer(C_INT) :: iy
    II=0.0d0
    DO iy=1,ny-1,2
      yp1=y(iy+1)-y(iy); ym1=y(iy-1)-y(iy)
      a1=-1.0d0/3.0d0*ym1+1.0d0/6.0d0*yp1+1.0d0/6.0d0*yp1*yp1/ym1
      a3=+1.0d0/3.0d0*yp1-1.0d0/6.0d0*ym1-1.0d0/6.0d0*ym1*ym1/yp1
      a2=yp1-ym1-a1-a3
      II=II+a1*f(iy-1)+a2*f(iy)+a3*f(iy+1)
    END DO
  END FUNCTION yintegr

#define rD0(f,g,k) dcmplx(sum(der(iy)%d0(-2:2)*dreal(f(ix,iz,iy-2:iy+2)%g)) ,sum(der(iy)%d0(-2:2)*dreal(f(ix,iz,iy-2:iy+2)%k)))
#define rD1(f,g,k) dcmplx(sum(der(iy)%d1(-2:2)*dreal(f(ix,iz,iy-2:iy+2)%g)) ,sum(der(iy)%d1(-2:2)*dreal(f(ix,iz,iy-2:iy+2)%k)))
#define rD2(f,g,k) dcmplx(sum(der(iy)%d2(-2:2)*dreal(f(ix,iz,iy-2:iy+2)%g)) ,sum(der(iy)%d2(-2:2)*dreal(f(ix,iz,iy-2:iy+2)%k)))
#define rD4(f,g,k) dcmplx(sum(der(iy)%d4(-2:2)*dreal(f(ix,iz,iy-2:iy+2)%g)) ,sum(der(iy)%d4(-2:2)*dreal(f(ix,iz,iy-2:iy+2)%k)))
#define D0(f,g) dcmplx(sum(der(iy)%d0(-2:2)*dreal(f(ix,iz,iy-2:iy+2)%g)) ,sum(der(iy)%d0(-2:2)*dimag(f(ix,iz,iy-2:iy+2)%g)))
#define D1(f,g) dcmplx(sum(der(iy)%d1(-2:2)*dreal(f(ix,iz,iy-2:iy+2)%g)) ,sum(der(iy)%d1(-2:2)*dimag(f(ix,iz,iy-2:iy+2)%g)))
#define D2(f,g) dcmplx(sum(der(iy)%d2(-2:2)*dreal(f(ix,iz,iy-2:iy+2)%g)) ,sum(der(iy)%d2(-2:2)*dimag(f(ix,iz,iy-2:iy+2)%g)))
#define D4(f,g) dcmplx(sum(der(iy)%d4(-2:2)*dreal(f(ix,iz,iy-2:iy+2)%g)) ,sum(der(iy)%d4(-2:2)*dimag(f(ix,iz,iy-2:iy+2)%g)))
#define iD0(g) dcmplx(sum(der(iy)%d0(-2:2)*dreal(VV(ix,iz,-2:2)%g)) ,sum(der(iy)%d0(-2:2)*dimag(VV(ix,iz,-2:2)%g)))
#define iD1(g) dcmplx(sum(der(iy)%d1(-2:2)*dreal(VV(ix,iz,-2:2)%g)) ,sum(der(iy)%d1(-2:2)*dimag(VV(ix,iz,-2:2)%g)))
#define iD2(g) dcmplx(sum(der(iy)%d2(-2:2)*dreal(VV(ix,iz,-2:2)%g)) ,sum(der(iy)%d2(-2:2)*dimag(VV(ix,iz,-2:2)%g)))
#define iD4(g) dcmplx(sum(der(iy)%d4(-2:2)*dreal(VV(ix,iz,-2:2)%g)) ,sum(der(iy)%d4(-2:2)*dimag(VV(ix,iz,-2:2)%g)))

  !--------------------------------------------------------------!
  !---COMPLEX----- derivative in the y-direction ----------------!
  SUBROUTINE COMPLEXderiv(f0,f1)
    complex(C_DOUBLE_COMPLEX), intent(in)  :: f0(-1:ny+1)
    complex(C_DOUBLE_COMPLEX), intent(out) :: f1(-1:ny+1)
    f1(0)=sum(d140(-2:2)*f0(-1:3))
    f1(-1)=sum(d14m1(-2:2)*f0(-1:3))
    f1(ny)=sum(d14n(-2:2)*f0(ny-3:ny+1))
    f1(ny+1)=sum(d14np1(-2:2)*f0(ny-3:ny+1))
    DO iy=1,ny-1
      f1(iy)=dcmplx(sum(der(iy)%d1(-2:2)*dreal(f0(iy-2:iy+2))) ,sum(der(iy)%d1(-2:2)*dimag(f0(iy-2:iy+2))))
    END DO
    f1(1)=f1(1)-(der(1)%d0(-1)*f1(0)+der(1)%d0(-2)*f1(-1))
    f1(2)=f1(2)-der(2)%d0(-2)*f1(0)
    f1(ny-1)=f1(ny-1)-(der(ny-1)%d0(1)*f1(ny)+der(ny-1)%d0(2)*f1(ny+1))
    f1(ny-2)=f1(ny-2)-der(ny-2)%d0(2)*f1(ny)
    f1(1:ny-1)=dcmplx(D0mat.bsr.dreal(f1(1:ny-1)),D0mat.bsr.dimag(f1(1:ny-1)))
  END SUBROUTINE COMPLEXderiv

  !--------------------------------------------------------------!
  !----------------- apply the boundary conditions --------------!
  SUBROUTINE applybc_0(EQ,bc0,bc0m1,RHS,rhs0,rhs0m1)
    real(C_DOUBLE), intent(inout) :: EQ(1:ny-1,-2:2)
    real(C_DOUBLE), intent(in) :: bc0(-2:2),bc0m1(-2:2)
    complex(C_DOUBLE_COMPLEX), intent(in) :: rhs0,rhs0m1
    complex(C_DOUBLE_COMPLEX), intent(inout) :: RHS(-1:ny+1)
    EQ(1,-1:2)=EQ(1,-1:2)-EQ(1,-2)*bc0m1(-1:2)/bc0m1(-2)
    EQ(1, 0:2)=EQ(1, 0:2)-EQ(1,-1)*bc0(0:2)/bc0(-1)
    EQ(2,-1:1)=EQ(2,-1:1)-EQ(2,-2)*bc0(0:2)/bc0(-1)
    RHS(1)=RHS(1)-EQ(1,-2)*rhs0m1/bc0m1(-2)
    RHS(1)=RHS(1)-EQ(1,-1)*rhs0/bc0(-1)
    RHS(2)=RHS(2)-EQ(2,-2)*rhs0/bc0(-1)
  END SUBROUTINE applybc_0

  SUBROUTINE applybc_n(EQ,bcn,bcnp1,RHS,rhsn,rhsnp1)
    real(C_DOUBLE), intent(inout) :: EQ(1:ny-1,-2:2)
    real(C_DOUBLE), intent(in) :: bcn(-2:2),bcnp1(-2:2)
    complex(C_DOUBLE_COMPLEX), intent(in) :: rhsn,rhsnp1
    complex(C_DOUBLE_COMPLEX), intent(inout) :: RHS(-1:ny+1)
    EQ(ny-1,-2:1)=EQ(ny-1,-2:1)-EQ(ny-1,2)*bcnp1(-2:1)/bcnp1(2)
    EQ(ny-1,-2:0)=EQ(ny-1,-2:0)-EQ(ny-1,1)*bcn(-2:0)/bcn(1)
    EQ(ny-2,-1:1)=EQ(ny-2,-1:1)-EQ(ny-2,2)*bcn(-2:0)/bcn(1)
    RHS(ny-1)=RHS(ny-1)-EQ(ny-1,2)*rhsnp1/bcnp1(2)
    RHS(ny-1)=RHS(ny-1)-EQ(ny-1,1)*rhsn/bcn(1)
    RHS(ny-2)=RHS(ny-2)-EQ(ny-2,2)*rhsn/bcn(1)
  END SUBROUTINE applybc_n


#define OS(iy,j) (ni*(der(iy)%d4(j)-2.0d0*k2*der(iy)%d2(j)+k2*k2*der(iy)%d0(j)))
#define SQ(iy,j) (ni*(der(iy)%d2(j)-k2*der(iy)%d0(j)))
  !--------------------------------------------------------------!
  !------------------- solve the linear system  -----------------!
  SUBROUTINE linsolve(lambda)
    real(C_DOUBLE), intent(in) :: lambda
    complex(C_DOUBLE_COMPLEX) :: A0,B0,An,Bn
    integer(C_INT) :: ix,iz
    complex(C_DOUBLE_COMPLEX) :: ialfa,ibeta,temp(-1:ny+1)
    real(C_DOUBLE) :: k2
    real(C_DOUBLE) :: ucor(-1:ny+1)
    DO ix=0,nx
      ialfa=dcmplx(0.0d0,alfa0*ix)
      DO iz=-nz,nz
        A0=0; An=A0; B0=0; Bn=0; 
        A0=A0-v0bc(-2)*B0/v0m1bc(-2); An=An-vnbc(2)*Bn/vnp1bc(2);
        ibeta=dcmplx(0.0d0,beta0*iz); k2=(alfa0*ix)**2.0d0+(beta0*iz)**2.0d0
        DO iy=1,ny-1
          D2vmat(iy,-2:2)=lambda*(der(iy)%d2(-2:2)-k2*der(iy)%d0(-2:2))-OS(iy,-2:2)
          etamat(iy,-2:2)=lambda*der(iy)%d0(-2:2)-SQ(iy,-2:2) 
        END DO
        CALL applybc_0(D2vmat,v0bc,v0m1bc,V(ix,iz,:)%v,A0,B0)
        CALL applybc_n(D2vmat,vnbc,vnp1bc,V(ix,iz,:)%v,An,Bn)
        CALL applybc_0(etamat,eta0bc,eta0m1bc,V(ix,iz,:)%u,(0.0D0,0.0D0),(0.0D0,0.0D0))
        CALL applybc_n(etamat,etanbc,etanp1bc,V(ix,iz,:)%u,(0.0D0,0.0D0),(0.0D0,0.0D0))
        CALL LU5decomp(D2vmat); CALL LU5decomp(etamat)
        V(ix,iz,1:ny-1)%v=dcmplx(D2vmat.bsr.dreal(V(ix,iz,1:ny-1)%v), D2vmat.bsr.dimag(V(ix,iz,1:ny-1)%v))
        V(ix,iz,0)%v=(A0-sum(V(ix,iz,1:3)%v*v0bc(0:2)))/v0bc(-1)
        V(ix,iz,-1)%v=(B0-sum(V(ix,iz,0:3)%v*v0m1bc(-1:2)))/v0m1bc(-2)
        V(ix,iz,ny)%v=(An-sum(V(ix,iz,ny-3:ny-1)%v*vnbc(-2:0)))/vnbc(1)
        V(ix,iz,ny+1)%v=(Bn-sum(V(ix,iz,ny-3:ny)%v*vnp1bc(-2:1)))/vnp1bc(2)
        V(ix,iz,1:ny-1)%u=dcmplx(etamat.bsr.dreal(V(ix,iz,1:ny-1)%u), etamat.bsr.dimag(V(ix,iz,1:ny-1)%u))
        V(ix,iz,0)%u=-sum(V(ix,iz,1:3)%u*eta0bc(0:2))/eta0bc(-1)
        V(ix,iz,-1)%u=-sum(V(ix,iz,0:3)%u*eta0m1bc(-1:2))/eta0m1bc(-2)
        V(ix,iz,ny)%u=-sum(V(ix,iz,ny-3:ny-1)%u*etanbc(-2:0))/etanbc(1)
        V(ix,iz,ny+1)%u=-sum(V(ix,iz,ny-3:ny)%u*etanp1bc(-2:1))/etanp1bc(2)
        IF (ix==0 .AND. iz==0) THEN 
          IF (abs(meanflowx)>0) THEN
            ucor=1
            ucor(1:ny-1)=etamat.bsr.ucor(1:ny-1)
            ucor(0)=-sum(ucor(1:3)*eta0bc(0:2))/eta0bc(-1)
            ucor(-1)=-sum(ucor(0:3)*eta0m1bc(-1:2))/eta0m1bc(-2)
            ucor(ny)=-sum(ucor(ny-3:ny-1)*etanbc(-2:0))/etanbc(1)
            ucor(ny+1)=-sum(ucor(ny-3:ny)*etanp1bc(-2:1))/etanp1bc(2)
            V(0,0,:)%u=dcmplx(dreal(V(0,0,:)%u)+(meanflowx-yintegr(dreal(V(0,0,:)%u)))/yintegr(ucor)*ucor,dimag(V(0,0,:)%u))
            V(0,0,:)%w=dcmplx(dreal(V(0,0,:)%w)+(meanflowz-yintegr(dreal(V(0,0,:)%w)))/yintegr(ucor)*ucor,dimag(V(0,0,:)%w))
          END IF
        ELSE
            CALL COMPLEXderiv(V(ix,iz,:)%v,V(ix,iz,:)%w)
            temp=(ialfa*V(ix,iz,:)%w-ibeta*V(ix,iz,:)%u)/k2
            V(ix,iz,:)%w=(ibeta*V(ix,iz,:)%w+ialfa*V(ix,iz,:)%u)/k2
            V(ix,iz,:)%u=temp
        END IF
      END DO
    END DO
  END SUBROUTINE linsolve

  !--------------------------------------------------------------!
  !------------------------ convolutions ------------------------!
  SUBROUTINE convolutions(V,VV)
    TYPE(MOMFLUX), intent(inout):: VV(0:nx,-nz:nz)
    TYPE(VELOCITY),intent(in) ::  V(0:nx,-nz:nz)
    Vd=0 ! XXX This must be improved and can be reduced to about three lines XXX
    Vd(1:nx+1,1:nz+1,1)=V(0:nx,0:nz)%u;           Vd(1:nx+1,1:nz+1,2)=V(0:nx,0:nz)%v;          Vd(1:nx+1,1:nz+1,3)=V(0:nx,0:nz)%w
    Vd(1:nx+1,nzd+1-nz:nzd,1)=V(0:nx,-nz:-1)%u;   Vd(1:nx+1,nzd+1-nz:nzd,2)=V(0:nx,-nz:-1)%v;  Vd(1:nx+1,nzd+1-nz:nzd,3)=V(0:nx,-nz:-1)%w
    CALL fftw_execute_dft_c2r(pIFT3,Vd,rVd)
    rVVd(1:2*nxd,1:nzd,1:3)=rVd(1:2*nxd,1:nzd,1:3)*rVd(1:2*nxd,1:nzd,1:3)*factor
    rVVd(1:2*nxd,1:nzd,4  )=rVd(1:2*nxd,1:nzd,1  )*rVd(1:2*nxd,1:nzd,2)*factor
    rVVd(1:2*nxd,1:nzd,5  )=rVd(1:2*nxd,1:nzd,2  )*rVd(1:2*nxd,1:nzd,3)*factor
    rVVd(1:2*nxd,1:nzd,6  )=rVd(1:2*nxd,1:nzd,1  )*rVd(1:2*nxd,1:nzd,3)*factor
    CALL fftw_execute_dft_r2c(pFFT6,rVVd,VVd)
    VV(0:nx,0:nz)%uu=VVd(1:nx+1,1:nz+1,1); VV(0:nx,-nz:-1)%uu=VVd(1:nx+1,nzd+1-nz:nzd,1)
    VV(0:nx,0:nz)%vv=VVd(1:nx+1,1:nz+1,2); VV(0:nx,-nz:-1)%vv=VVd(1:nx+1,nzd+1-nz:nzd,2)
    VV(0:nx,0:nz)%ww=VVd(1:nx+1,1:nz+1,3); VV(0:nx,-nz:-1)%ww=VVd(1:nx+1,nzd+1-nz:nzd,3)
    VV(0:nx,0:nz)%uv=VVd(1:nx+1,1:nz+1,4); VV(0:nx,-nz:-1)%uv=VVd(1:nx+1,nzd+1-nz:nzd,4)
    VV(0:nx,0:nz)%vw=VVd(1:nx+1,1:nz+1,5); VV(0:nx,-nz:-1)%vw=VVd(1:nx+1,nzd+1-nz:nzd,5)
    VV(0:nx,0:nz)%uw=VVd(1:nx+1,1:nz+1,6); VV(0:nx,-nz:-1)%uw=VVd(1:nx+1,nzd+1-nz:nzd,6)
  END SUBROUTINE convolutions


  !--------------------------------------------------------------!
  !-------------------------- buildRHS --------------------------! 
  SUBROUTINE buildrhs(timescheme)
    integer(C_INT) :: iy,ix,iz,i
    complex(C_DOUBLE_COMPLEX) :: ialfa,ibeta,rhsu,rhsv,rhsw
    TYPE(rhstype) :: temp(0:nx,-nz:nz)
    real(C_DOUBLE) :: k2
    DO iy=-1,2
      CALL convolutions(V(:,:,iy),VV(:,:,iy))
    END DO
    DO iy=1,ny-1
      DO i=-2,1
        VV(:,:,i)=VV(:,:,i+1)
      END DO
      CALL convolutions(V(:,:,iy+2),VV(:,:,2))
      DO ix=0,nx
        ialfa=dcmplx(0.0d0,alfa0*ix)
        DO iz=-nz,nz
          ibeta=dcmplx(0.0d0,beta0*iz)
          k2=(alfa0*ix)**2.0d0 + (beta0*iz)**2.0d0
          rhsu=-ialfa*iD0(uu)-iD1(uv)-ibeta*iD0(uw)
          rhsv=-ialfa*iD0(uv)-iD1(vv)-ibeta*iD0(vw)
          rhsw=-ialfa*iD0(uw)-iD1(vw)-ibeta*iD0(ww)
          CALL timescheme(newrhs(0,ix,iz)%D2v, oldrhs(iy,ix,iz)%D2v, D2(V,v)-k2*D0(V,v),&
                          sum(OS(iy,-2:2)*V(ix,iz,iy-2:iy+2)%v),&
                          ialfa*(ialfa*iD1(uu)+iD2(uv)+ibeta*iD1(uw))+&
                          ibeta*(ialfa*iD1(uw)+iD2(vw)+ibeta*iD1(ww))-k2*rhsv)   !D2v
          IF (ix==0 .AND. iz==0) THEN
            CALL timescheme(newrhs(0,0,0)%eta,oldrhs(iy,0,0)%eta,rD0(V,u,w),&
                            ni*rD2(V,u,w),&
                            dcmplx(dreal(rhsu)+meanpx,dreal(rhsw)+meanpz)) ! (Ubar, Wbar)
          ELSE
            CALL timescheme(newrhs(0,ix,iz)%eta, oldrhs(iy,ix,iz)%eta,ibeta*D0(V,u)-ialfa*D0(V,w),&
                            sum(SQ(iy,-2:2)*[ibeta*V(ix,iz,iy-2:iy+2)%u-ialfa*V(ix,iz,iy-2:iy+2)%w]),&    
                            ibeta*rhsu-ialfa*rhsw) !eta
          END IF
          V(ix,iz,iy-2)%u=newrhs(-2,ix,iz)%eta; V(ix,iz,iy-2)%v=newrhs(-2,ix,iz)%D2v
        END DO
      END DO
      temp=newrhs(-2,:,:); newrhs(-2,:,:)=newrhs(-1,:,:); newrhs(-1,:,:)=newrhs(0,:,:); newrhs(0,:,:)=temp
    END DO
    DO i=-2,-1
      V(:,:,ny+i)%u=newrhs(i,:,:)%eta; V(:,:,ny+i)%v=newrhs(i,:,:)%D2v
    END DO
  END SUBROUTINE buildrhs


  !--------------------------------------------------------------!
  !-----is it necessary in Fortran to keep (ix,iz) separated?----!
  !--------------------- ODE schemes library --------------------!
  SUBROUTINE CN_AB(rhs,old,unkn,impl,expl)
    complex(C_DOUBLE_COMPLEX), intent(out)   :: rhs
    complex(C_DOUBLE_COMPLEX), intent(inout) :: old
    complex(C_DOUBLE_COMPLEX), intent(in) :: unkn,impl,expl
    rhs=2.0d0/deltat*unkn+impl+3.0d0*expl-old
    old=expl
  END SUBROUTINE CN_AB

  SUBROUTINE RK1_rai(rhs,old,unkn,impl,expl)
    complex(C_DOUBLE_COMPLEX), intent(out)   :: rhs
    complex(C_DOUBLE_COMPLEX), intent(inout) :: old
    complex(C_DOUBLE_COMPLEX), intent(in) :: unkn,impl,expl
    rhs=120.0d0/32.0d0/deltat*unkn+impl+2.0d0*expl
    old=expl
  END SUBROUTINE RK1_rai

  SUBROUTINE RK2_rai(rhs,old,unkn,impl,expl)
    complex(C_DOUBLE_COMPLEX), intent(out)   :: rhs
    complex(C_DOUBLE_COMPLEX), intent(inout) :: old
    complex(C_DOUBLE_COMPLEX), intent(in) :: unkn,impl,expl
    rhs=120.0d0/(8.0d0*deltat)*unkn+impl+50.0d0/8.0d0*expl-34.0d0/8.0d0*old
    old=expl
  END SUBROUTINE RK2_rai

  SUBROUTINE RK3_rai(rhs,old,unkn,impl,expl)
    complex(C_DOUBLE_COMPLEX), intent(out)   :: rhs
    complex(C_DOUBLE_COMPLEX), intent(inout) :: old
    complex(C_DOUBLE_COMPLEX), intent(in) :: unkn,impl,expl
    rhs=120.0d0/(20.0d0*deltat)*unkn+impl+90.0d0/20.0d0*expl-50.0d0/20.0d0*old
  END SUBROUTINE RK3_rai

  !--------------------------------------------------------------!
  !-------------------- read_restart_file -----------------------! 
  SUBROUTINE read_restart_file(V)
    TYPE(VELOCITY), intent(out) :: V(0:nx,-nz:nz,-1:ny+1)
    integer :: ix,iy,iz
    OPEN( unit=75, file="Dati.cart.out", status='old', access='stream', form='unformatted', action='read' )
    READ (75) V(0:nx,-nz:nz,-1:ny+1)
    CLOSE(75)
  END SUBROUTINE read_restart_file

  !--------------------------------------------------------------!
  !-------------------- save_restart_file -----------------------! 
  SUBROUTINE save_restart_file(V)
    TYPE(VELOCITY), intent(out) :: V(0:nx,-nz:nz,-1:ny+1)
    OPEN( unit=85, file="Dati.cart.out", access='stream', form='unformatted', action='write' )
    WRITE (85) V(0:nx,-nz:nz,-1:ny+1)
    CLOSE(85)
  END SUBROUTINE save_restart_file

END MODULE dnsdata
