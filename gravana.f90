!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,f,dx,ncell)
  use amr_parameters
  use poisson_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::f ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  !================================================================
  integer::idim,i
  real(dp)::gmass,emass,xmass,ymass,zmass,rr,rx,ry,rz

  ! Constant vector
  if(gravity_type==1)then
     do idim=1,ndim
        do i=1,ncell
           f(i,idim)=gravity_params(idim)
        end do
     end do
  end if

  ! Point mass
  if(gravity_type==2)then
     gmass=gravity_params(1) ! GM
     emass=dx
     emass=gravity_params(2) ! Softening length
     xmass=gravity_params(3) ! Point mass coordinates
     ymass=gravity_params(4)
     zmass=gravity_params(5)
     do i=1,ncell
        rx=0.0d0; ry=0.0d0; rz=0.0d0
        rx=x(i,1)-xmass
#if NDIM>1
        ry=x(i,2)-ymass
#endif
#if NDIM>2
        rz=x(i,3)-zmass
#endif
        rr=sqrt(rx**2+ry**2+rz**2+emass**2)
        f(i,1)=-gmass*rx/rr**3
#if NDIM>1
        f(i,2)=-gmass*ry/rr**3
#endif
#if NDIM>2
        f(i,3)=-gmass*rz/rr**3
#endif
     end do
  end if

end subroutine gravana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
#ifdef EFE
subroutine phi_ana(rr,pp,ngrid,xx)
#else
subroutine phi_ana(rr,pp,ngrid)
#endif
  use amr_commons
  use poisson_commons
  use constants, only: twopi
#ifdef EFE
  use mond_commons, only: mond_a0, &
      & g_ext, g_ext_hat_x, g_ext_hat_y, g_ext_hat_z, beta
#else
  use mond_commons, only: mond_a0
#endif
  implicit none
  integer::ngrid
  real(dp),dimension(1:nvector)::rr,pp
#ifdef EFE
  real(dp),dimension(1:nvector,1:ndim)::xx
  real(dp)::r_inv,w,L_0,mu_ext,mu_ext_inv,cos_theta
#endif
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------

  integer :: i
  real(dp):: fourpi

  fourpi=2*twopi

#if NDIM==1
  do i=1,ngrid
     pp(i)=multipole(1)*fourpi/2*rr(i)
  end do
#endif
#if NDIM==2
  do i=1,ngrid
     pp(i)=multipole(1)*2*log(rr(i))
  end do
#endif
#if NDIM==3
  do i=1,ngrid
#ifdef EFE
    r_inv = 1.0d0/rr(i)
    w = g_ext/mond_a0
    L_0 = 1.0d0/(1.0d0 + w)
    ! beta is 1/(mond_a0 + g_ext)
    mu_ext = g_ext*beta
    mu_ext_inv = 1.0d0/mu_ext
    cos_theta = (xx(i,1)*g_ext_hat_x + xx(i,2)*g_ext_hat_y + xx(i,3)*g_ext_hat_z)*r_inv
    pp(i) = -multipole(1)*r_inv*mu_ext_inv/dsqrt(1.0d0 + L_0*(1.0d0 - cos_theta*cos_theta))
#else
     pp(i)=SQRT(multipole(1)*mond_a0)*DLOG(rr(i))
     !pp(i)=-multipole(1)/rr(i)
#endif
  end do
#endif
end subroutine phi_ana
