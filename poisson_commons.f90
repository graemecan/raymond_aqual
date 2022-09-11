module poisson_commons
  use amr_commons
  use poisson_parameters

  real(dp),allocatable,dimension(:)  ::phi,phi_old       ! Potential
  real(dp),allocatable,dimension(:)  ::rho               ! Density
  real(dp),allocatable,dimension(:,:)::f                 ! 3-force

  real(dp),allocatable,dimension(:)  ::rho_top   ! Density at last CIC level

  ! Multigrid lookup table for amr -> mg index mapping
  integer, allocatable, dimension(:) :: lookup_mg   ! Lookup table

  ! Communicator arrays for multigrid levels
  type(communicator), allocatable, dimension(:,:) :: active_mg
  type(communicator), allocatable, dimension(:,:) :: emission_mg

  ! Minimum MG level
  integer :: levelmin_mg

  ! Multigrid safety switch
  logical, allocatable, dimension(:) :: safe_mode

  ! Multipole coefficients
  real(dp),dimension(1:ndim+1)::multipole

end module poisson_commons

module mond_commons
   use amr_parameters
   implicit none

   integer :: mond_n = 2         !Read in from namelist
   real(dp) :: imond_a0 = 1.2d-20 !Read in from namelist
   integer  :: maxstoredcells = 0   !Read in from namelist
   integer  :: diaginterplevel = 7 !Read in from namelist

#ifdef EFE
   real(dp) :: g_ext_x, g_ext_y, g_ext_z, g_ext, beta
   real(dp) :: g_ext_hat_x, g_ext_hat_y, g_ext_hat_z
#endif

   integer :: tnbors = twondim*ndim

   real(dp), allocatable, dimension(:,:) :: amr_interp_pts
   real(dp), dimension(1:3,1:3) :: hhh
   real(dp) :: mond_a0, mond_a0_sqd, oneover_mond_a0

   integer, dimension(1:3,1:6,1:8) :: kkk, lll, mmm
   integer, dimension(1:6,1:8) :: uuu, yyy, ttt
   integer, dimension(1:ndim,1:2,1:ndim,1:4) :: sss
   integer, dimension(1:8) :: gs_ordering
   integer, dimension(1:6) :: bbb

contains   

subroutine init_mond_arrays
   use amr_parameters
   implicit none

   integer :: nx_loc
   real(dp) :: scale
   real(dp) :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

   kkk(1,1,1:8) = (/1,0,1,0,1,0,1,0/) ; lll(1,1,1:8) = (/0,0,0,0,0,0,0,0/)
   kkk(1,2,1:8) = (/0,2,0,2,0,2,0,2/) ; lll(1,2,1:8) = (/0,0,0,0,0,0,0,0/)
   kkk(1,3,1:8) = (/1,3,1,0,1,3,1,0/) ; lll(1,3,1:8) = (/3,0,0,0,3,0,0,0/)
   kkk(1,4,1:8) = (/0,2,4,2,0,2,4,2/) ; lll(1,4,1:8) = (/0,0,0,4,0,0,0,4/)
   kkk(1,5,1:8) = (/1,5,1,5,1,0,1,0/) ; lll(1,5,1:8) = (/5,0,5,0,0,0,0,0/)
   kkk(1,6,1:8) = (/0,2,0,2,6,2,6,2/) ; lll(1,6,1:8) = (/0,0,0,0,0,6,0,6/)

   kkk(2,1,1:8) = (/3,3,0,0,3,3,0,0/) ; lll(2,1,1:8) = (/0,0,0,0,0,0,0,0/)
   kkk(2,2,1:8) = (/0,0,4,4,0,0,4,4/) ; lll(2,2,1:8) = (/0,0,0,0,0,0,0,0/)
   kkk(2,3,1:8) = (/3,2,0,2,3,2,0,2/) ; lll(2,3,1:8) = (/0,3,0,0,0,3,0,0/)
   kkk(2,4,1:8) = (/1,0,1,4,1,0,1,4/) ; lll(2,4,1:8) = (/0,0,4,0,0,0,4,0/)
   kkk(2,5,1:8) = (/3,3,5,5,3,3,0,0/) ; lll(2,5,1:8) = (/5,5,0,0,0,0,0,0/)
   kkk(2,6,1:8) = (/0,0,4,4,6,6,4,4/) ; lll(2,6,1:8) = (/0,0,0,0,0,0,6,6/)

   kkk(3,1,1:8) = (/5,5,5,5,0,0,0,0/) ; lll(3,1,1:8) = (/0,0,0,0,0,0,0,0/)
   kkk(3,2,1:8) = (/0,0,0,0,6,6,6,6/) ; lll(3,2,1:8) = (/0,0,0,0,0,0,0,0/)
   kkk(3,3,1:8) = (/5,2,5,2,0,2,0,2/) ; lll(3,3,1:8) = (/0,5,0,5,0,0,0,0/)
   kkk(3,4,1:8) = (/1,0,1,0,1,6,1,6/) ; lll(3,4,1:8) = (/0,0,0,0,6,0,6,0/)
   kkk(3,5,1:8) = (/5,5,5,5,0,0,4,4/) ; lll(3,5,1:8) = (/0,0,4,4,0,0,0,0/)
   kkk(3,6,1:8) = (/3,3,0,0,6,6,6,6/) ; lll(3,6,1:8) = (/0,0,0,0,3,3,0,0/)

   mmm(1,1,1:8) = (/2,1,4,3,6,5,8,7/) ; uuu(1,1:8) = (/1,0,1,0,1,0,1,0/)
   mmm(1,2,1:8) = (/2,1,4,3,6,5,8,7/) ; uuu(2,1:8) = (/0,1,0,1,0,1,0,1/)
   mmm(1,3,1:8) = (/4,3,2,1,8,7,6,5/) ; uuu(3,1:8) = (/1,1,0,0,1,1,0,0/)
   mmm(1,4,1:8) = (/4,3,2,1,8,7,6,5/) ; uuu(4,1:8) = (/0,0,1,1,0,0,1,1/)
   mmm(1,5,1:8) = (/6,5,8,7,2,1,4,3/) ; uuu(5,1:8) = (/1,1,1,1,0,0,0,0/)
   mmm(1,6,1:8) = (/6,5,8,7,2,1,4,3/) ; uuu(6,1:8) = (/0,0,0,0,1,1,1,1/)

   mmm(2,1,1:8) = (/3,4,1,2,7,8,5,6/) ; yyy(1,1:8) = (/2,1,4,3,6,5,8,7/)
   mmm(2,2,1:8) = (/3,4,1,2,7,8,5,6/) ; yyy(2,1:8) = (/2,1,4,3,6,5,8,7/)
   mmm(2,3,1:8) = (/4,3,2,1,8,7,6,5/) ; yyy(3,1:8) = (/3,4,1,2,7,8,5,6/)
   mmm(2,4,1:8) = (/4,3,2,1,8,7,6,5/) ; yyy(4,1:8) = (/3,4,1,2,7,8,5,6/)
   mmm(2,5,1:8) = (/7,8,5,6,3,4,1,2/) ; yyy(5,1:8) = (/5,6,7,8,1,2,3,4/)
   mmm(2,6,1:8) = (/7,8,5,6,3,4,1,2/) ; yyy(6,1:8) = (/5,6,7,8,1,2,3,4/)

   mmm(3,1,1:8) = (/5,6,7,8,1,2,3,4/)
   mmm(3,2,1:8) = (/5,6,7,8,1,2,3,4/)
   mmm(3,3,1:8) = (/6,5,8,7,2,1,4,3/)
   mmm(3,4,1:8) = (/6,5,8,7,2,1,4,3/)
   mmm(3,5,1:8) = (/7,8,5,6,3,4,1,2/)
   mmm(3,6,1:8) = (/7,8,5,6,3,4,1,2/)

   ttt(1,1:8) = (/-2,2,-2,2,-2,2,-2,2/)
   ttt(2,1:8) = (/1,-1,1,-1,1,-1,1,-1/)
   ttt(3,1:8) = (/3,3,-3,-3,3,3,-3,-3/)
   ttt(4,1:8) = (/-4,-4,4,4,-4,-4,4,4/)
   ttt(5,1:8) = (/5,5,5,5,-5,-5,-5,-5/)
   ttt(6,1:8) = (/-6,-6,-6,-6,6,6,6,6/)

   bbb(1:6) = (/1,1,2,2,3,3/)

#if NDIM == 1
   sss(1,1,1,1:4) = (/ tnbors+1,-1, 1,-1/)
   sss(1,2,1,1:4) = (/ 2,-1, tnbors+1,-1/)
#endif
#if NDIM == 2
   sss(1,1,1,1:4) = (/ tnbors+1,-1, 1,-1/)
   sss(1,1,2,1:4) = (/ 8, 4, 5, 2/)

   sss(1,2,1,1:4) = (/ 3,-1, tnbors+1,-1/)
   sss(1,2,2,1:4) = (/ 4, 7, 2, 6/)

   sss(2,1,1,1:4) = (/ 3, 6, 1, 5/)
   sss(2,1,2,1:4) = (/ tnbors+1,-1, 2,-1/)

   sss(2,2,1,1:4) = (/ 7, 3, 8, 1/)
   sss(2,2,2,1:4) = (/ 4,-1, tnbors+1,-1/)
#endif
#if NDIM == 3
   sss(1,1,1,1:4) = (/ tnbors+1,-1, 1,-1/)
   sss(1,1,2,1:4) = (/ 5,11, 2, 7/)
   sss(1,1,3,1:4) = (/ 6,12, 3,13/)

   sss(1,2,1,1:4) = (/ 4,-1, tnbors+1,-1/)
   sss(1,2,2,1:4) = (/10, 5, 8, 2/)
   sss(1,2,3,1:4) = (/16, 6, 9, 3/)

   sss(2,1,1,1:4) = (/ 4, 8, 1, 7/)
   sss(2,1,2,1:4) = (/ tnbors+1,-1, 2,-1/)
   sss(2,1,3,1:4) = (/18, 6,14, 3/)

   sss(2,2,1,1:4) = (/10, 4,11, 1/)
   sss(2,2,2,1:4) = (/ 5,-1, tnbors+1,-1/)
   sss(2,2,3,1:4) = (/ 6,17, 3,15/)

   sss(3,1,1,1:4) = (/ 4, 9, 1,13/)
   sss(3,1,2,1:4) = (/ 5,15, 2,14/)
   sss(3,1,3,1:4) = (/ tnbors+1,-1, 3,-1/)

   sss(3,2,1,1:4) = (/16, 4,12, 1/)
   sss(3,2,2,1:4) = (/17, 5,18, 2/)
   sss(3,2,3,1:4) = (/ 6,-1, tnbors+1,-1/)
#endif

   hhh(1,1:3) = (/1.0d0,4.0d0,4.0d0/)
   hhh(2,1:3) = (/4.0d0,1.0d0,4.0d0/)
   hhh(3,1:3) = (/4.0d0,4.0d0,1.0d0/)

   gs_ordering = (/1,4,6,7,2,3,5,8/)

   nx_loc = icoarse_max-icoarse_min+1
   scale  = boxlen/dble(nx_loc)

   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   mond_a0 = (imond_a0*100.0d0*scale_t**2)/scale_l
   mond_a0_sqd = mond_a0**2
   oneover_mond_a0 = 1.0d0/mond_a0

#ifdef EFE
   ! Assume that external field acceleration is given in same physical units as a0
   g_ext_x = (g_ext_x*100.0d0*scale_t**2)/scale_l
   g_ext_y = (g_ext_y*100.0d0*scale_t**2)/scale_l
   g_ext_z = (g_ext_z*100.0d0*scale_t**2)/scale_l

   g_ext = dsqrt(g_ext_x**2 + g_ext_y**2 + g_ext_z**2)
   g_ext_hat_x = g_ext_x/g_ext
   g_ext_hat_y = g_ext_y/g_ext
   g_ext_hat_z = g_ext_z/g_ext

   beta = 1.0d0/(mond_a0 + g_ext)
#endif

   allocate(amr_interp_pts(1:maxstoredcells,1:twotondim))

end subroutine init_mond_arrays

!subroutine calculate_mu(ilevel,phi_nbor,mu)
!   use amr_parameters
!   implicit none

!   integer, intent(in)  :: ilevel
!   real(dp), dimension(1:tnbors+1), intent(in)  :: phi_nbor
!   real(dp), dimension(1:ndim,1:2), intent(out) :: mu

!   real(dp) :: dx, scale, gradx, grady, gradz, x
!   integer  :: inbor, idim, c, nx_loc

!   dx  = 0.5d0**ilevel
!   c   = tnbors+1

!   nx_loc = icoarse_max-icoarse_min+1
!   scale  = boxlen/dble(nx_loc)

!#if NDIM == 1

         !Face 1
!         gradx = (phi_nbor(c) - phi_nbor(1))/dx
!         x = ABS(gradx)/(mond_a0*scale)
!         mu(1,1) = x/((1.0d0+x**mond_n)**(1.0/mond_n))
         !Face 2
!         gradx = (phi_nbor(2) - phi_nbor(c))/dx
!         x = ABS(gradx)/(mond_a0*scale)
!         mu(1,2) = x/((1.0d0+x**mond_n)**(1.0/mond_n))

!#endif
!#if NDIM == 2

         !Face 1
!         gradx = (phi_nbor(c) - phi_nbor(1))/dx
!         grady = (phi_nbor(8)+phi_nbor(4)-phi_nbor(5)-phi_nbor(2))/(4.0d0*dx)
!         x = DSQRT(gradx**2 + grady**2)/(mond_a0*scale)
!         mu(1,1) = x/((1.0d0+x**mond_n)**(1.0/mond_n))
         !Face 2
!         gradx = (phi_nbor(3) - phi_nbor(c))/dx
!         grady = (phi_nbor(4)+phi_nbor(7)-phi_nbor(2)-phi_nbor(6))/(4.0d0*dx)
!         x = DSQRT(gradx**2 + grady**2)/(mond_a0*scale)
!         mu(1,2) = x/((1.0d0+x**mond_n)**(1.0/mond_n))
         !Face 3
!         gradx = (phi_nbor(3)+phi_nbor(6)-phi_nbor(1)-phi_nbor(5))/(4.0d0*dx)
!         grady = (phi_nbor(c)-phi_nbor(2))/dx
!         x = DSQRT(gradx**2 + grady**2)/(mond_a0*scale)
!         mu(2,1) = x/((1.0d0+x**mond_n)**(1.0/mond_n))
         !Face 4
!         gradx = (phi_nbor(7)+phi_nbor(3)-phi_nbor(8)-phi_nbor(1))/(4.0d0*dx)
!         grady = (phi_nbor(4)-phi_nbor(c))/dx
!         x = DSQRT(gradx**2 + grady**2)/(mond_a0*scale)
!         mu(2,2) = x/((1.0d0+x**mond_n)**(1.0/mond_n))

!#endif
!#if NDIM == 3

         !Face 1
!         gradx = (phi_nbor(c) - phi_nbor(1))/dx
!         grady = (phi_nbor(5)+phi_nbor(11)-phi_nbor(2)-phi_nbor(7))/(4.0d0*dx)
!         gradz = (phi_nbor(6)+phi_nbor(12)-phi_nbor(3)-phi_nbor(13))/(4.0d0*dx)
!         x = SQRT(gradx**2 + grady**2 + gradz**2)/(mond_a0*scale)
!         mu(1,1) = x/((1.0d0+x**mond_n)**(1.0/mond_n))
         !Face 2
!         gradx = (phi_nbor(4) - phi_nbor(c))/dx
!         grady = (phi_nbor(10)+phi_nbor(5)-phi_nbor(8)-phi_nbor(2))/(4.0d0*dx)
!         gradz = (phi_nbor(16)+phi_nbor(6)-phi_nbor(9)-phi_nbor(3))/(4.0d0*dx)
!         x = SQRT(gradx**2 + grady**2 + gradz**2)/(mond_a0*scale)
!         mu(1,2) = x/((1.0d0+x**mond_n)**(1.0/mond_n))
         !Face 3
!         gradx = (phi_nbor(4)+phi_nbor(8)-phi_nbor(1)-phi_nbor(7))/(4.0d0*dx)
!         grady = (phi_nbor(c)-phi_nbor(2))/dx
!         gradz = (phi_nbor(18)+phi_nbor(6)-phi_nbor(14)-phi_nbor(3))/(4.0d0*dx)
!         x = SQRT(gradx**2 + grady**2 + gradz**2)/(mond_a0*scale)
!         mu(2,1) = x/((1.0d0+x**mond_n)**(1.0/mond_n))
         !Face 4
!         gradx = (phi_nbor(10)+phi_nbor(4)-phi_nbor(11)-phi_nbor(1))/(4.0d0*dx)
!         grady = (phi_nbor(5)-phi_nbor(c))/dx
!         gradz = (phi_nbor(6)+phi_nbor(17)-phi_nbor(3)-phi_nbor(15))/(4.0d0*dx)
!         x = SQRT(gradx**2 + grady**2 + gradz**2)/(mond_a0*scale)
!         mu(2,2) = x/((1.0d0+x**mond_n)**(1.0/mond_n))
         !Face 5
!         gradx = (phi_nbor(4)+phi_nbor(9)-phi_nbor(1)-phi_nbor(13))/(4.0d0*dx)
!         grady = (phi_nbor(5)+phi_nbor(15)-phi_nbor(2)-phi_nbor(14))/(4.0d0*dx)
!         gradz = (phi_nbor(c)-phi_nbor(3))/dx
!         x = SQRT(gradx**2 + grady**2 + gradz**2)/(mond_a0*scale)
!         mu(3,1) = x/((1.0d0+x**mond_n)**(1.0/mond_n))
         !Face 6
!         gradx = (phi_nbor(16)+phi_nbor(4)-phi_nbor(12)-phi_nbor(1))/(4.0d0*dx)
!         grady = (phi_nbor(17)+phi_nbor(5)-phi_nbor(18)-phi_nbor(2))/(4.0d0*dx)
!         gradz = (phi_nbor(6)-phi_nbor(c))/dx
!         x = SQRT(gradx**2 + grady**2 + gradz**2)/(mond_a0*scale)
!         mu(3,2) = x/((1.0d0+x**mond_n)**(1.0/mond_n))
 
!#endif

!end subroutine calculate_mu

!subroutine calculate_mu(ilevel,phi_nbor,mu)
!   use amr_parameters
!   implicit none

!   integer, intent(in)  :: ilevel
!   real(dp), dimension(1:tnbors+1), intent(in)  :: phi_nbor
!   real(dp), dimension(1:ndim,1:2), intent(out) :: mu

!   real(dp) :: dx, grad_sqrd, phi_sum, x
!   integer  :: inbor, idim, gdim, i

!   dx  = 0.5d0**ilevel

   !Loop over cell faces
   !do inbor=1,2
   !   do idim=1,ndim
         !Loop over cells in stencil
!         grad_sqrd = 0.0d0
!         do gdim=1,ndim
!            phi_sum = 0.0d0
!            do i=1,2
!               if (sss(idim,inbor,gdim,i) < 0) cycle
!               phi_sum = phi_sum + phi_nbor(sss(idim,inbor,gdim,i))
!            enddo
!            do i=3,4
!               if (sss(idim,inbor,gdim,i) < 0) cycle
!               phi_sum = phi_sum - phi_nbor(sss(idim,inbor,gdim,i))
!            enddo
!            grad_sqrd = grad_sqrd + (phi_sum/(hhh(idim,gdim)*dx))**2
!         enddo
         !x = SQRT(grad_sqrd)/(mond_a0*scale)
!         mu(idim,inbor) = SQRT(grad_sqrd/(mond_a0_L_sqd + grad_sqrd))
    !  enddo
   !enddo

!end subroutine calculate_mu

real(dp) function mu_function(gsqd)
   implicit none

   real(dp), intent(in) :: gsqd
   real(dp) :: g

   !x = SQRT(gsqd)*oneover_mond_a0
   !mu_function = x/(1+x**mond_n)**(1.0d0/mond_n)
   g = SQRT(gsqd)
   mu_function = g/(mond_a0+g)
   !mu_function = 1.0d0

end function mu_function

real(dp) function inverse_mu_function(y)
   implicit none

   real(dp), intent(in) :: y

   inverse_mu_function = (0.5d0 + 0.5d0*DSQRT(1.0d0 + 4.0d0/(y**mond_n)))**(1.0d0/mond_n)

end function inverse_mu_function

end module mond_commons
