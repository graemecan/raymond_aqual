! ------------------------------------------------------------------------
! Multigrid Poisson solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all MG-fine-level related routines
!
! Used variables:
!                       finest(AMR)level     coarse(MG)levels
!     -----------------------------------------------------------------
!     potential            phi            active_mg(myid,ilevel)%u(:,1)
!     physical RHS         rho            active_mg(myid,ilevel)%u(:,2)
!     residual             f(:,1)         active_mg(myid,ilevel)%u(:,3)
!     BC-modified RHS      f(:,2)                  N/A
!     mask                 f(:,3)         active_mg(myid,ilevel)%u(:,4)
!
! ------------------------------------------------------------------------


! ------------------------------------------------------------------------
! Mask restriction (top-down, OBSOLETE, UNUSED)
! ------------------------------------------------------------------------

subroutine restrict_mask_fine(ifinelevel,allmasked)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in)  :: ifinelevel
   logical, intent(out) :: allmasked

   integer :: ind_c_cell,ind_f_cell

   integer :: iskip_f_amr, iskip_c_amr, iskip_c_mg
   integer :: igrid_f_amr, igrid_c_amr, igrid_c_mg
   integer :: icell_f_amr, icell_c_amr, icell_c_mg

   real(dp) :: ngpmask
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
   icoarselevel=ifinelevel-1
   allmasked=.true.

   if(ifinelevel==1) return

   ! Loop over coarse cells of the coarse active comm for myid
   do ind_c_cell=1,twotondim
      iskip_c_amr=ncoarse+(ind_c_cell-1)*ngridmax
      iskip_c_mg =(ind_c_cell-1)*active_mg(myid,icoarselevel)%ngrid

      ! Loop over coarse grids
      do igrid_c_mg=1,active_mg(myid,icoarselevel)%ngrid
         igrid_c_amr=active_mg(myid,icoarselevel)%igrid(igrid_c_mg)
         icell_c_amr=iskip_c_amr+igrid_c_amr
         icell_c_mg =iskip_c_mg +igrid_c_mg

         if(son(icell_c_amr)==0) then
            ! Cell is not refined
            ngpmask      = -1.0d0
         else
            igrid_f_amr=son(icell_c_amr)
            ngpmask = 0.0d0
            do ind_f_cell=1,twotondim
               iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax
               icell_f_amr=iskip_f_amr+igrid_f_amr
               ngpmask=ngpmask+f(icell_f_amr,3)
            end do
            ngpmask=ngpmask/dtwotondim
         end if
         ! Store cell mask
         active_mg(myid,icoarselevel)%u(icell_c_mg,4)=ngpmask
         allmasked=allmasked .and. (ngpmask<=0.0)
      end do
   end do

end subroutine restrict_mask_fine

! ------------------------------------------------------------------------
! Mask restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_mask_fine_reverse(ifinelevel)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr

   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_amr
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_amr

   real(dp) :: ngpmask
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
   icoarselevel=ifinelevel-1

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg)
         icell_f_amr=igrid_f_amr+iskip_f_amr

         ! Get coarse grid AMR index and CPU id
         icell_c_amr=father(igrid_f_amr)
         ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax
         cpu_amr=cpu_map(father(igrid_c_amr))

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg

         ! Stack cell volume fraction in coarse cell
         ngpmask=(1d0+f(icell_f_amr,3))/2/dtwotondim
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)+ngpmask
      end do
   end do
end subroutine restrict_mask_fine_reverse

! ------------------------------------------------------------------------
! Residual computation
! ------------------------------------------------------------------------
! RAyMOND - substantially modified for AQUAL solver
subroutine cmp_residual_mg_fine(ilevel)
   ! Computes the residual the fine (AMR) level, and stores it into f(:,1)
   use amr_commons
   use poisson_commons
   use mond_commons
   implicit none
   integer, intent(in) :: ilevel

   real(dp), dimension(1:tnbors+1) :: phi_nbor
   real(dp), dimension(1:twondim) :: phi_vals
   real(dp), dimension(1:twondim) :: mu

   real(dp) :: dx, oneoverdx2, nb_sum, factor, grad_sqrd, phi_sum
   real(dp) :: oneoverdx, oneoverfourdx
   integer  :: ngrid, ind, gdim, i, interp_cnt
   integer  :: igrid_mg, idim, inbor, nbors_ind
   integer  :: igrid_amr, icell_amr, iskip_amr
   integer  :: kgshift, lgshift, igrid_nbor_amr
   integer  :: cpu_nbor_amr, icell_nbor_amr, ifathercell_nbor_amr
   integer  :: upper_cell_index, upper_cell_amr, upper_grid_amr
   integer  :: nbor_grid_amr, nbor_cell_index, nbor_cell_amr

   ! Set constants
   dx  = 0.5d0**ilevel
   oneoverdx2 = 1.0d0/(dx*dx)
   oneoverdx = 1.0d0/dx
   oneoverfourdx = 1.0d0/(4.0d0*dx)

   ngrid=active(ilevel)%ngrid

   ! Loop over cells
   do ind=1,twotondim
      interp_cnt=0
      iskip_amr = ncoarse+(ind-1)*ngridmax

      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr

         phi_nbor(tnbors+1) = phi(icell_amr)  ! Value of potential on center cell

         ! SCAN FLAG TEST
         if(flag2(icell_amr)/ngridmax==0) then
            nbors_ind = 0
            do inbor=1,twondim
               do idim=1,ndim
                  nbors_ind = nbors_ind+1
                  ! Get neighbour gridshift
                  kgshift = kkk(idim,inbor,ind)
                  lgshift = lll(idim,inbor,ind)
                  ! Get neighbor grid
                  if(kgshift==0) then
                     igrid_nbor_amr = igrid_amr
                  else if (lgshift==0) then
                     igrid_nbor_amr = son(nbor(igrid_amr,kgshift))
                  else
                     igrid_nbor_amr = son(nbor(son(nbor(igrid_amr,kgshift)),lgshift))
                  end if
                  if(igrid_nbor_amr==0) then
                     ! No DIAGONAL neighbor
                     ! Use central phi value
                     !phi_nbor(nbors_ind) = phi(icell_amr)
                     interp_cnt=interp_cnt+1
                     if(ilevel>=diaginterplevel)then
                        phi_nbor(nbors_ind)=amr_interp_pts(interp_cnt,ind)
                     else
                        phi_nbor(nbors_ind)=phi(icell_amr)
                     endif
                  else
                     ! Get phi values on neighbouring cells
                     icell_nbor_amr  = igrid_nbor_amr + &
                        (ncoarse + (mmm(idim,inbor,ind)-1)*ngridmax)
                     phi_nbor(nbors_ind) = phi(icell_nbor_amr)
                  endif
               end do
            end do
         else ! PERFORM SCAN
            if(f(icell_amr,3)<=0.0) then
               f(icell_amr,1)=0.0d0
               cycle
            end if
            nbors_ind = 0
            do inbor=1,twondim
               do idim=1,ndim
                  nbors_ind = nbors_ind+1
                  ! Get neighbor grid shift
                  kgshift = kkk(idim,inbor,ind)
                  lgshift = lll(idim,inbor,ind)
                  ! Get neighbor grid and its parent cell
                  if (kgshift==0) then
                     igrid_nbor_amr = igrid_amr
                     ifathercell_nbor_amr = father(igrid_nbor_amr)
                  else if (lgshift==0) then
                     igrid_nbor_amr = son(nbor(igrid_amr,kgshift))
                     ifathercell_nbor_amr = nbor(igrid_amr,kgshift)
                  else
                     if (son(nbor(igrid_amr,kgshift)) == 0) then
                     ! special case when we can't "leapfrog" to the diagonal cell because
                     ! the first shift takes us out of the AMR mesh.
                     ! First, integer division to find which number of cell we have in upper level grid
                     ! and then we calculate the amr index of that upper level grid
                        upper_cell_amr = nbor(igrid_amr,kgshift)
                        upper_cell_index = (upper_cell_amr/ngridmax)+1
                        ! Check if we need to move to the neighbour grid at the upper level or not
                        if (uuu(lgshift,upper_cell_index)==0) then
                            nbor_grid_amr = upper_cell_amr - ncoarse - (upper_cell_index - 1)*ngridmax
                        else
                            upper_grid_amr = upper_cell_amr - ncoarse - (upper_cell_index - 1)*ngridmax
                            nbor_grid_amr = son(nbor(upper_grid_amr,lgshift))
                        endif
                        ! Determine cell index of neighbour cell, depending on direction of 2nd shift
                        nbor_cell_index = yyy(lgshift,upper_cell_index)
                        ! Calculate amr index of neighbour cell
                        nbor_cell_amr = nbor_grid_amr + ncoarse + (nbor_cell_index - 1)*ngridmax
                        ! Find son grid in AMR mesh, if it exists
                        igrid_nbor_amr = son(nbor_cell_amr)
                        ifathercell_nbor_amr = nbor_cell_amr
                     else
                        igrid_nbor_amr = son(nbor(son(nbor(igrid_amr,kgshift)),lgshift))
                        ifathercell_nbor_amr = nbor(son(nbor(igrid_amr,kgshift)),lgshift)
                     endif
                  end if
                  if(igrid_nbor_amr==0) then
                     ! No neighbor
                     ! Use previously stored interpolated values from upper level
                     interp_cnt=interp_cnt+1
                     phi_nbor(nbors_ind) = amr_interp_pts(interp_cnt,ind)
                  else
                     ! Fetch neighbor cell id
                     icell_nbor_amr = igrid_nbor_amr + (ncoarse + (mmm(idim,inbor,ind)-1)*ngridmax)
                     phi_nbor(nbors_ind)  = phi(icell_nbor_amr)
                  end if
               end do
            end do
         end if ! END SCAN TEST

         grad_sqrd = ((phi_nbor(tnbors+1)-phi_nbor(1))*oneoverdx)**2 + &
                     ((phi_nbor( 5)+phi_nbor(11)-phi_nbor( 2)-phi_nbor( 7))*oneoverfourdx)**2 + &
                     ((phi_nbor( 6)+phi_nbor(12)-phi_nbor( 3)-phi_nbor(13))*oneoverfourdx)**2
         mu(1) = mu_function(grad_sqrd)
         grad_sqrd = ((phi_nbor(4)-phi_nbor(tnbors+1))*oneoverdx)**2 + &
                     ((phi_nbor(10)+phi_nbor( 5)-phi_nbor( 8)-phi_nbor( 2))*oneoverfourdx)**2 + &
                     ((phi_nbor(16)+phi_nbor( 6)-phi_nbor( 9)-phi_nbor( 3))*oneoverfourdx)**2
         mu(4) = mu_function(grad_sqrd)
         grad_sqrd = ((phi_nbor( 4)+phi_nbor( 8)-phi_nbor( 1)-phi_nbor( 7))*oneoverfourdx)**2 + &
                     ((phi_nbor(tnbors+1)-phi_nbor(2))*oneoverdx)**2 + &
                     ((phi_nbor(18)+phi_nbor( 6)-phi_nbor(14)-phi_nbor( 3))*oneoverfourdx)**2
         mu(2) = mu_function(grad_sqrd)
         grad_sqrd = ((phi_nbor(10)+phi_nbor( 4)-phi_nbor(11)-phi_nbor( 1))*oneoverfourdx)**2 + &
                     ((phi_nbor( 5)-phi_nbor(tnbors+1))*oneoverdx)**2 + &
                     ((phi_nbor( 6)+phi_nbor(17)-phi_nbor( 3)-phi_nbor(15))*oneoverfourdx)**2
         mu(5) = mu_function(grad_sqrd)
         grad_sqrd = ((phi_nbor( 4)+phi_nbor( 9)-phi_nbor( 1)-phi_nbor(13))*oneoverfourdx)**2 + &
                     ((phi_nbor( 5)+phi_nbor(15)-phi_nbor( 2)-phi_nbor(14))*oneoverfourdx)**2 + &
                     ((phi_nbor(tnbors+1)-phi_nbor( 3))*oneoverdx)**2
         mu(3) = mu_function(grad_sqrd)
         grad_sqrd = ((phi_nbor(16)+phi_nbor( 4)-phi_nbor(12)-phi_nbor( 1))*oneoverfourdx)**2 + &
                     ((phi_nbor(17)+phi_nbor( 5)-phi_nbor(18)-phi_nbor( 2))*oneoverfourdx)**2 + &
                     ((phi_nbor( 6)-phi_nbor(tnbors+1))*oneoverdx)**2
         mu(6) = mu_function(grad_sqrd)

         phi_vals = phi_nbor(1:twondim)
         nb_sum = mu(1)*phi_vals(1)+mu(2)*phi_vals(2)+mu(3)*phi_vals(3) + &
                  mu(4)*phi_vals(4)+mu(5)*phi_vals(5)+mu(6)*phi_vals(6)
         factor = mu(1)+mu(2)+mu(3)+mu(4)+mu(5)+mu(6)

         ! Store ***MINUS THE RESIDUAL*** in f(:,1), using BC-modified RHS
         ! Now includes central phi terms
         f(icell_amr,1) = -oneoverdx2*( nb_sum - factor*phi(icell_amr) )+f(icell_amr,2)
      end do
   end do

end subroutine cmp_residual_mg_fine

! ##################################################################
! ##################################################################

subroutine cmp_residual_norm2_fine(ilevel, norm2)
   use amr_commons
   use poisson_commons
   implicit none

   integer,  intent(in)  :: ilevel
   real(kind=8), intent(out) :: norm2

   real(kind=8) :: dx2
   integer  :: ngrid
   integer  :: ind, igrid_mg
   integer  :: igrid_amr, icell_amr, iskip_amr

   ! Set constants
   dx2  = (0.5d0**ilevel)**ndim
   ngrid=active(ilevel)%ngrid

   norm2 = 0.0d0
   ! Loop over cells
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax
      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr
         if(f(icell_amr,3)<=0.0) then      ! Do not count masked cells
            cycle
         end if
         norm2 = norm2 + f(icell_amr,1)**2
      end do
   end do
   norm2 = dx2*norm2

end subroutine cmp_residual_norm2_fine

subroutine cmp_ivar_norm2_fine(ilevel, ivar, norm2)
   use amr_commons
   use poisson_commons
   implicit none

   integer,  intent(in)  :: ilevel, ivar
   real(kind=8), intent(out) :: norm2

   real(kind=8) :: dx2
   integer  :: ngrid
   integer  :: ind, igrid_mg
   integer  :: igrid_amr, icell_amr, iskip_amr

   ! Set constants
   dx2  = (0.5d0**ilevel)**ndim
   ngrid=active(ilevel)%ngrid

   norm2 = 0.0d0
   ! Loop over cells
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax
      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr
         if(f(icell_amr,3)<=0.0) then      ! Do not count masked cells
            cycle
         end if
         if(ivar>0)then
            norm2 = norm2 + f(icell_amr,ivar)**2
         else
            norm2 = norm2 + phi(icell_amr)**2
         endif
      end do
   end do
   norm2 = dx2*norm2

end subroutine cmp_ivar_norm2_fine

! ------------------------------------------------------------------------
! Gauss-Seidel smoothing
! ------------------------------------------------------------------------
! RAyMOND - substantially modified for AQUAL solver
subroutine gauss_seidel_mg_fine(ilevel,ind)
   use amr_commons
   use pm_commons
   use poisson_commons
   use mond_commons
   implicit none
   integer, intent(in) :: ilevel
   integer, intent(in) :: ind

   real(dp), dimension(1:tnbors+1) :: phi_nbor
   real(dp), dimension(1:twondim) :: phi_vals
   real(dp), dimension(1:twondim) :: mu

   real(dp) :: dx, dx2, nb_sum, factor, grad_sqrd, phi_sum, dterm
   real(dp) :: oneoverdx, oneoverfourdx, factor_tol, omega, oneoverdx2
   integer  :: ngrid, ind_loop, gdim, i, interp_cnt, ix
   integer  :: igrid_mg, idim, inbor, nbors_ind
   integer  :: igrid_amr, icell_amr, iskip_amr
   integer  :: kgshift, lgshift, igrid_nbor_amr
   integer  :: cpu_nbor_amr, icell_nbor_amr, ifathercell_nbor_amr
   integer  :: upper_cell_index, upper_cell_amr, upper_grid_amr
   integer  :: nbor_grid_amr, nbor_cell_index, nbor_cell_amr

   ! Set constants
   dx2  = (0.5d0**ilevel)**2
   dx   = 0.5d0**ilevel
   oneoverdx = 1.0d0/dx
   oneoverdx2 = 1.0d0/dx2
   oneoverfourdx = 1.0d0/(4.0d0*dx)
   factor_tol = 1.0d-8

   !RAyMOND - use a more aggressive SOR method for levelmin
   !due to boundary condition issues.
   if (ilevel == levelmin) then
      omega = 1.4d0
   else
      omega = 1.0d0
   endif

   ngrid=active(ilevel)%ngrid
   interp_cnt=0

      iskip_amr = ncoarse+(ind-1)*ngridmax

      ! Loop over active grids
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr

         phi_nbor(tnbors+1) = phi(icell_amr)

         ! Read scan flag
         if(flag2(icell_amr)/ngridmax==0) then
            ! Use max-speed "dumb" Gauss-Seidel for "inner" cells
            ! Those cells are active, have all their neighbors active
            ! and all neighbors are in the AMR+MG trees
            nbors_ind = 0
            do inbor=1,twondim
               do idim=1,ndim
                  nbors_ind = nbors_ind+1
                  ! Get neighbour gridshift
                  kgshift = kkk(idim,inbor,ind)
                  lgshift = lll(idim,inbor,ind)
                  ! Get neighbor grid
                  if(kgshift==0) then
                     igrid_nbor_amr = igrid_amr
                  else if (lgshift==0) then
                     igrid_nbor_amr = son(nbor(igrid_amr,kgshift))
                  else
                     igrid_nbor_amr = son(nbor(son(nbor(igrid_amr,kgshift)),lgshift))
                  end if
                  if(igrid_nbor_amr==0) then
                     ! No DIAGONAL neighbor
                     ! Use central phi value
                     !phi_nbor(nbors_ind) = phi(icell_amr)
                     interp_cnt=interp_cnt+1
                     if(ilevel>=diaginterplevel)then
                        phi_nbor(nbors_ind)=amr_interp_pts(interp_cnt,ind)
                     else
                        phi_nbor(nbors_ind)=phi(icell_amr)
                     end if
                  else
                     ! Get phi values on neighbouring cells
                     icell_nbor_amr  = igrid_nbor_amr + &
                        (ncoarse + (mmm(idim,inbor,ind)-1)*ngridmax)
                     phi_nbor(nbors_ind) = phi(icell_nbor_amr)
                  endif
               end do
            end do 
         else
            ! Use the finer "solve" Gauss-Seidel near boundaries,
            ! with all necessary checks
            if (f(icell_amr,3)<=0.0) cycle
            if (safe_mode(ilevel) .and. f(icell_amr,3)<1.0) cycle
            ! Use more complex mu calculation valid at boundaries

            nbors_ind = 0
            do inbor=1,twondim
               do idim=1,ndim
                  nbors_ind = nbors_ind+1
                  ! Get neighbor grid shift
                  kgshift = kkk(idim,inbor,ind)
                  lgshift = lll(idim,inbor,ind)
                  ! Get neighbor grid and its parent cell
                  if (kgshift==0) then
                     igrid_nbor_amr = igrid_amr
                     ifathercell_nbor_amr = father(igrid_nbor_amr)
                  else if (lgshift==0) then
                     igrid_nbor_amr = son(nbor(igrid_amr,kgshift))
                     ifathercell_nbor_amr = nbor(igrid_amr,kgshift)
                  else
                     if (son(nbor(igrid_amr,kgshift)) == 0) then
                     ! special case when we can't "leapfrog" to the diagonal cell because
                     ! the first shift takes us out of the AMR mesh.
                     ! First, integer division to find which number of cell we have in upper level grid
                     ! and then we calculate the amr index of that upper level grid
                        upper_cell_amr = nbor(igrid_amr,kgshift)
                        upper_cell_index = (upper_cell_amr/ngridmax)+1
                        ! Check if we need to move to the neighbour grid at the upper level or not
                        if (uuu(lgshift,upper_cell_index)==0) then
                            nbor_grid_amr = upper_cell_amr - ncoarse - (upper_cell_index - 1)*ngridmax
                        else
                            upper_grid_amr = upper_cell_amr - ncoarse - (upper_cell_index - 1)*ngridmax
                            nbor_grid_amr = son(nbor(upper_grid_amr,lgshift))
                        endif
                        ! Determine cell index of neighbour cell, depending on direction of 2nd shift
                        nbor_cell_index = yyy(lgshift,upper_cell_index)
                        ! Calculate amr index of neighbour cell
                        nbor_cell_amr = nbor_grid_amr + ncoarse + (nbor_cell_index - 1)*ngridmax
                        ! Find son grid in AMR mesh, if it exists
                        igrid_nbor_amr = son(nbor_cell_amr)
                        ifathercell_nbor_amr = nbor_cell_amr
                     else
                        igrid_nbor_amr = son(nbor(son(nbor(igrid_amr,kgshift)),lgshift))
                        ifathercell_nbor_amr = nbor(son(nbor(igrid_amr,kgshift)),lgshift)
                     endif
                  end if
                  if(igrid_nbor_amr==0) then
                     ! No neighbor
                     ! Use previously stored interpolated values from upper level
                     interp_cnt=interp_cnt+1
                     phi_nbor(nbors_ind) = amr_interp_pts(interp_cnt,ind)
                  else
                     icell_nbor_amr = igrid_nbor_amr + (ncoarse + (mmm(idim,inbor,ind)-1)*ngridmax)
                     phi_nbor(nbors_ind)  = phi(icell_nbor_amr)
                 end if
               end do
            end do
         end if
         ! Update the potential, solving for potential on icell_amr

         grad_sqrd = ((phi_nbor(tnbors+1)-phi_nbor(1))*oneoverdx)**2 + &
                     ((phi_nbor( 5)+phi_nbor(11)-phi_nbor( 2)-phi_nbor( 7))*oneoverfourdx)**2 + &
                     ((phi_nbor( 6)+phi_nbor(12)-phi_nbor( 3)-phi_nbor(13))*oneoverfourdx)**2
         mu(1) = mu_function(grad_sqrd)
         grad_sqrd = ((phi_nbor(4)-phi_nbor(tnbors+1))*oneoverdx)**2 + &
                     ((phi_nbor(10)+phi_nbor( 5)-phi_nbor( 8)-phi_nbor( 2))*oneoverfourdx)**2 + &
                     ((phi_nbor(16)+phi_nbor( 6)-phi_nbor( 9)-phi_nbor( 3))*oneoverfourdx)**2
         mu(4) = mu_function(grad_sqrd)
         grad_sqrd = ((phi_nbor( 4)+phi_nbor( 8)-phi_nbor( 1)-phi_nbor( 7))*oneoverfourdx)**2 + &
                     ((phi_nbor(tnbors+1)-phi_nbor(2))*oneoverdx)**2 + &
                     ((phi_nbor(18)+phi_nbor( 6)-phi_nbor(14)-phi_nbor( 3))*oneoverfourdx)**2
         mu(2) = mu_function(grad_sqrd)
         grad_sqrd = ((phi_nbor(10)+phi_nbor( 4)-phi_nbor(11)-phi_nbor( 1))*oneoverfourdx)**2 + &
                     ((phi_nbor( 5)-phi_nbor(tnbors+1))*oneoverdx)**2 + &
                     ((phi_nbor( 6)+phi_nbor(17)-phi_nbor( 3)-phi_nbor(15))*oneoverfourdx)**2
         mu(5) = mu_function(grad_sqrd)
         grad_sqrd = ((phi_nbor( 4)+phi_nbor( 9)-phi_nbor( 1)-phi_nbor(13))*oneoverfourdx)**2 + &
                     ((phi_nbor( 5)+phi_nbor(15)-phi_nbor( 2)-phi_nbor(14))*oneoverfourdx)**2 + &
                     ((phi_nbor(tnbors+1)-phi_nbor( 3))*oneoverdx)**2
         mu(3) = mu_function(grad_sqrd)
         grad_sqrd = ((phi_nbor(16)+phi_nbor( 4)-phi_nbor(12)-phi_nbor( 1))*oneoverfourdx)**2 + &
                     ((phi_nbor(17)+phi_nbor( 5)-phi_nbor(18)-phi_nbor( 2))*oneoverfourdx)**2 + &
                     ((phi_nbor( 6)-phi_nbor(tnbors+1))*oneoverdx)**2
         mu(6) = mu_function(grad_sqrd)

         phi_vals = phi_nbor(1:twondim)
         nb_sum = mu(1)*phi_vals(1)+mu(2)*phi_vals(2)+mu(3)*phi_vals(3) + &
                  mu(4)*phi_vals(4)+mu(5)*phi_vals(5)+mu(6)*phi_vals(6)
         factor = mu(1)+mu(2)+mu(3)+mu(4)+mu(5)+mu(6)

         ! Phi is not updated if "factor" is zero, which is usually more likely in cosmological runs
         if (factor .ne. 0.0d0) then
            phi(icell_amr) = (1.0d0-omega)*phi(icell_amr) + omega*(nb_sum - dx2*f(icell_amr,2)) / factor
         end if

      end do

end subroutine gauss_seidel_mg_fine

! ------------------------------------------------------------------------
! Residual restriction (top-down, OBSOLETE, UNUSED)
! ------------------------------------------------------------------------

subroutine restrict_residual_fine(ifinelevel)
   ! Restrict fine (AMR) residual at level ifinelevel using injection
   ! into coarser residual at level ifinelevel-1
   ! Restricted residual is stored into the RHS at the coarser level
   use amr_commons
   use pm_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   real(dp) :: val
   real(dp) :: dtwotondim = (twotondim)

   integer  :: icoarselevel
   integer  :: ngrid_c, ind_c, iskip_c_amr, iskip_c_mg
   integer  :: igrid_c_amr, icell_c_amr, icell_c_mg, igrid_c_mg
   integer  :: ind_f, igrid_f_amr, iskip_f_amr, icell_f_amr

   icoarselevel=ifinelevel-1

   ! Loop over coarse MG cells
   ngrid_c=active_mg(myid,icoarselevel)%ngrid
   do ind_c=1,twotondim
      iskip_c_amr = ncoarse + (ind_c-1)*ngridmax
      iskip_c_mg  = (ind_c-1)*ngrid_c

      do igrid_c_mg=1,ngrid_c
         igrid_c_amr = active_mg(myid,icoarselevel)%igrid(igrid_c_mg)
         icell_c_amr = igrid_c_amr + iskip_c_amr
         icell_c_mg  = igrid_c_mg  + iskip_c_mg

         ! Get AMR child grid
         igrid_f_amr = son(icell_c_amr)
         if(igrid_f_amr==0) then
            ! Nullify residual (coarser RHS)
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = 0.0d0
            cycle
         end if

         val = 0.0d0
         ! Loop over child (fine MG) cells
         do ind_f=1,twotondim
            iskip_f_amr = ncoarse + (ind_f-1)*ngridmax
            icell_f_amr = igrid_f_amr + iskip_f_amr

            if (f(icell_f_amr,3)<=0.0) cycle
            val = val + f(icell_f_amr,1)
         end do
         ! Store restricted residual into RHS of coarse level
         active_mg(myid,icoarselevel)%u(icell_c_mg,2) = val/dtwotondim
      end do
   end do
end subroutine restrict_residual_fine


! ------------------------------------------------------------------------
! Residual restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_residual_fine_reverse(ifinelevel)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr

   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_amr
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_amr

   real(dp) :: res
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
   icoarselevel=ifinelevel-1

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg)
         icell_f_amr=igrid_f_amr+iskip_f_amr
         ! Is fine cell masked?
         if(f(icell_f_amr,3)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         icell_c_amr=father(igrid_f_amr)
         ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax
         cpu_amr=cpu_map(father(igrid_c_amr))

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg

         ! Is coarse cell masked?
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0d0) cycle

         ! Stack fine cell residual in coarse cell rhs
         res=f(icell_f_amr,1)/dtwotondim
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,2)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,2)+res
      end do
   end do
end subroutine restrict_residual_fine_reverse

! ------------------------------------------------------------------------
! RAyMOND - Phi restriction (bottom-up, added for AQUAL solver)
! ------------------------------------------------------------------------

subroutine restrict_phi_fine_reverse(ifinelevel)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr
   
   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_amr
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_amr

   real(dp) :: res
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
   icoarselevel=ifinelevel-1

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg)
         icell_f_amr=igrid_f_amr+iskip_f_amr
         ! Is fine cell masked?
         if(f(icell_f_amr,3)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         icell_c_amr=father(igrid_f_amr)
         ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax
         cpu_amr=cpu_map(father(igrid_c_amr))

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg

         ! Is coarse cell masked?
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0d0) cycle

         ! Stack fine cell phi in coarse cell (new MG field)
         ! Only difference with residual restriction is in these
         ! two lines...
         res=phi(icell_f_amr)/dtwotondim
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,5)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,5)+res
      end do
   end do
end subroutine restrict_phi_fine_reverse

! ------------------------------------------------------------------------
! Interpolation and correction
! ------------------------------------------------------------------------

subroutine interpolate_and_correct_fine(ifinelevel)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   integer  :: i, ind_father, ind_average, ind_f, iskip_f_amr
   integer  :: ngrid_f, istart, nbatch
   integer  :: icell_c_amr, igrid_c_amr, igrid_c_mg, icell_c_mg
   integer  :: icoarselevel, ind_c, cpu_amr

   real(dp) :: a, b, c, d, coeff
   real(dp), dimension(1:8)     :: bbb
   integer,  dimension(1:8,1:8) :: ccc

   integer,  dimension(1:nvector), save               :: igrid_f_amr, icell_amr
   integer,  dimension(1:nvector,1:threetondim), save :: nbors_father_cells
   integer,  dimension(1:nvector,1:twotondim), save   :: nbors_father_grids
   real(dp), dimension(1:nvector), save               :: corr

   ! Local constants
   a = 1d0/4d0**ndim
   b = 3*a
   c = 9*a
   d = 27*a
   icoarselevel=ifinelevel-1

   bbb(:)  =(/a ,b ,b ,c ,b ,c ,c ,d/)

   ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
   ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
   ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
   ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
   ccc(:,5)=(/19,20,22,23,10,11,13,14/)
   ccc(:,6)=(/21,20,24,23,12,11,15,14/)
   ccc(:,7)=(/25,26,22,23,16,17,13,14/)
   ccc(:,8)=(/27,26,24,23,18,17,15,14/)

   ! Loop over fine grids by vector sweeps
   ngrid_f=active(ifinelevel)%ngrid
   do istart=1,ngrid_f,nvector

      ! Gather nvector grids
      nbatch=MIN(nvector,ngrid_f-istart+1)
      do i=1,nbatch
         igrid_f_amr(i)=active(ifinelevel)%igrid(istart+i-1)
      end do

      ! Compute father (coarse) cell index
      do i=1,nbatch
         icell_amr(i)=father(igrid_f_amr(i))
      end do

      ! Gather 3x3x3 neighboring parent cells
      call get3cubefather(icell_amr,nbors_father_cells,nbors_father_grids, &
              nbatch,ifinelevel)

      ! Update solution for fine grid cells
      do ind_f=1,twotondim
         iskip_f_amr = ncoarse+(ind_f-1)*ngridmax

         do i=1,nbatch
            ! Compute fine cell indices
            icell_amr(i) = iskip_f_amr + igrid_f_amr(i)
         end do
         corr=0.0d0

         ! Loop over relevant parent cells
         do ind_average=1,twotondim
            ind_father = ccc(ind_average,ind_f)
            coeff      = bbb(ind_average)
            do i=1,nbatch
               if(f(icell_amr(i),3)<=0.0) then
                  corr(i)=0.0d0        ! Fine cell is masked : no correction
                  cycle
               end if
               icell_c_amr = nbors_father_cells(i,ind_father)
               ind_c       = (icell_c_amr-ncoarse-1)/ngridmax + 1
               igrid_c_amr = icell_c_amr - ncoarse - (ind_c-1)*ngridmax
               cpu_amr     = cpu_map(father(igrid_c_amr))
               igrid_c_mg  = lookup_mg(igrid_c_amr)
               if(igrid_c_mg<=0) cycle

               icell_c_mg=(ind_c-1)*active_mg(cpu_amr,icoarselevel)%ngrid+igrid_c_mg
               corr(i)=corr(i)+coeff*active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,1)
            end do
         end do

         ! Correct potential
         do i=1,nbatch
            phi(icell_amr(i))=phi(icell_amr(i))+corr(i)
         end do

      end do
      ! End loop over cells

   end do
   ! End loop over grids
end subroutine interpolate_and_correct_fine


! ------------------------------------------------------------------------
! Flag setting
! ------------------------------------------------------------------------

subroutine set_scan_flag_fine(ilevel)
   use amr_commons
   use poisson_commons
   implicit none

   integer, intent(in) :: ilevel

   integer :: ind, ngrid, scan_flag
   integer :: igrid_mg, inbor, idim, igshift
   integer :: igrid_amr, igrid_nbor_amr

   integer :: iskip_amr, icell_amr, icell_nbor_amr

   integer, dimension(1:3,1:2,1:8) :: iii, jjj

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ngrid = active(ilevel)%ngrid

   ! Loop over cells and set fine SCAN flag
   do ind=1,twotondim
      iskip_amr = ncoarse+(ind-1)*ngridmax
      do igrid_mg=1,ngrid
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr

         if(f(icell_amr,3)==1.0) then
            scan_flag=0       ! Init flag to 'no scan needed'
            scan_flag_loop: do inbor=1,2
               do idim=1,ndim
                  igshift = iii(idim,inbor,ind)
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  end if

                  if(igrid_nbor_amr==0) then
                     scan_flag=1
                     exit scan_flag_loop
                  else
                     icell_nbor_amr = igrid_nbor_amr + &
                           ncoarse+(jjj(idim,inbor,ind)-1)*ngridmax
                     if(f(icell_nbor_amr,3)<=0.0) then
                        scan_flag=1
                        exit scan_flag_loop
                     end if
                  end if
               end do
            end do scan_flag_loop
         else
            scan_flag=1
         end if
         ! Update flag2 with scan flag,
         ! BEWARE as lookup_mg backups are stored in flag2
         ! Safety init:
         if(flag2(icell_amr)>ngridmax .or. flag2(icell_amr)<0) flag2(icell_amr)=0
         ! Do NOT overwrite flag2 !
         flag2(icell_amr)=flag2(icell_amr)+ngridmax*scan_flag
      end do
   end do
end subroutine set_scan_flag_fine
