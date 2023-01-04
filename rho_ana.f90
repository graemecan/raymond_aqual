!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine rho_ana(x,d,dx,ncell)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  use constants, only: twopi
  implicit none
  integer ::ncell                         ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector)       ::d ! Density
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates analytical Poisson source term.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! d(i) is the density field in user units.
  !================================================================
  integer::i
  real(dp)::rr,rx,ry,rz
  real(dp)::rp,xmass,ymass,zmass,M0,fourpi

  rp=gravity_params(1) ! Softening length
  xmass=gravity_params(2) ! Point mass coordinates
  ymass=gravity_params(3)
  zmass=gravity_params(4)

  !Parameters for Plummer
  M0 = 1.0d0 !Mass in 10**10 M_sol
  fourpi = 2*twopi

  do i=1,ncell
     rx=0.0d0; ry=0.0d0; rz=0.0d0
     rx=x(i,1)-xmass
#if NDIM>1
     ry=x(i,2)-ymass
#endif
#if NDIM>2
     rz=x(i,3)-zmass
#endif
     rr=sqrt(rx**2+ry**2+rz**2)
     d(i)=(3.0d0*M0/(fourpi*rp**3))*(1.0d0+(rr/rp)**2)**(-5.0/2.0)
  end do

end subroutine rho_ana
