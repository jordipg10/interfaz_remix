! !> Computes aqueous species concentrations after iteration of WMA using Euler fully implicit in chemical reactions for kinetic system
! !! We assume all primary species are aqueous
! !! The Jacobians are computed analytically
! subroutine water_mixing_iter_EfI_eq_kin_anal(this,c1_old,c2nc_ig,c_tilde,conc_nc,porosity,Delta_t)
!     use aqueous_chemistry_m, only: aqueous_chemistry_c, initialise_iterative_method, Newton_EfI_rk_eq_kin_aq_anal
!     implicit none
! !> Arguments
!     class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
!     real(kind=8), intent(in) :: c1_old(:)
!     real(kind=8), intent(in) :: c2nc_ig(:)
!     real(kind=8), intent(in) :: c_tilde(:)
!     real(kind=8), intent(out) :: conc_nc(:)
!     !real(kind=8), intent(out) :: conc_comp(:)
!     real(kind=8), intent(in), optional :: porosity !> at this target
!     real(kind=8), intent(in), optional :: Delta_t !> time step
! !> Variables
!     real(kind=8), allocatable :: c1(:) !> aqueous primary concentrations
!     integer(kind=4) :: n_p !> number of primary species
!     integer(kind=4) :: niter !> number of Newton iterations
!     real(kind=8) :: mu !> Newton initialisation parameter
!     logical :: CV_flag !> convergence flag
! !> Pre-process
!     n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
!     !> We assign aqueous primary species concentrations
!         c1=this%get_c1()
!     !> Newton initialisation parameter     
!         mu=0d0 
!     ! !> Update old kinetic reaction rates
!     ! call this%update_rk_old() !> we update old values of kinetic reaction rates
!     ! call this%solid_chemistry%update_rk_old() !> we update old values of kinetic reaction rates
! !> Process
!         do !> loop until convergence is reached
!         !> Initial guess aquoeus primary concentrations
!             call initialise_iterative_method(c1_old,c1,mu,this%concentrations(1:n_p))
!         !> We apply Newton method to compute aqueous concentrations
!             call this%Newton_EfI_rk_eq_kin_aq_anal(c2nc_ig,c_tilde,Delta_t,conc_nc,niter,CV_flag)
!         !> We check convergence
!             if (.not. CV_flag) then !> no CV
!                 if (mu<1d0) then
!                     mu=mu+0.25 !> we increase Newton initialisation parameter
!                 else
!                     error stop "Iterative method does not CV"
!                 end if
!             else 
!                 exit
!             end if
!         end do
!     !> We compute component concentrations
!         !call this%compute_conc_comp_aq()
!         !conc_comp=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,conc_nc)
! !< Post-process
!     deallocate(c1)
! end subroutine 