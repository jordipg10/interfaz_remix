! !> Computes aqueous species concentrations after iteration of WMA using Euler fully implicit in chemical reactions for kinetic system
! !! We assume all primary species are aqueous
! !! The Jacobians are computed analytically
! !! We DO NOT apply lumping to the mixing ratios of kinetic reaction rates
! !! We use option 2 (see documentation)
! subroutine water_mixing_iter_EfI_eq_kin_anal_ideal_opt2(this,c1_old,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,conc_nc)
!     use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, initialise_c1_aq_iterative_method, &
!         Newton_EfI_rk_eq_kin_aq_anal_ideal_opt2
!     implicit none
! !> Arguments
!     class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
!     real(kind=8), intent(in) :: c1_old(:) !> aqueous primary concentrations at previous time step
!     real(kind=8), intent(in) :: c_tilde(:) !> variable activity species concentrations after mixing
!     !real(kind=8), intent(in) :: rk_tilde_old(:) !> old kinetic reaction rates after mixing (de momento no se usa)
!     real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rates after mixing, ordered in upstream & downstream respectively
!     !real(kind=8), intent(in) :: rk_tilde_new(:) !> new kinetic reaction rates after mixing (de momento no se usa)
!     real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target

!     real(kind=8), intent(out) :: conc_nc(:) !> variable activity species concentrations (already allocated)
!     !real(kind=8), intent(in), optional :: porosity !> at this target
!     !real(kind=8), intent(in), optional :: Delta_t !> time step
! !> Variables
!     real(kind=8), allocatable :: c1(:) !> aqueous primary concentrations
!     integer(kind=4) :: n_p !> number of primary species
!     integer(kind=4) :: k_div !> counter of time step divisions
!     integer(kind=4) :: k !> counter of time steps
!     integer(kind=4) :: niter !> number of Newton iterations
!     real(kind=8) :: mu !> Newton initialisation parameter
!     real(kind=8) :: Delta_t_bis !> reduced time step
!     logical :: CV_flag !> convergence flag
! !> Pre-process
!     n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
!     !> We get aqueous primary species concentrations
!         c1=this%get_c1()
!     !> Newton initialisation parameter
!         mu=0d0 
!     !> Update old kinetic reaction rates
!     call this%update_rk_old() !> we update old values of kinetic reaction rates
!     call this%solid_chemistry%update_rk_old() !> we update old values of kinetic reaction rates
! !> Process
!     !> Initial guess aquoeus primary concentrations
!         call this%initialise_c1_aq_iterative_method(c1_old,c1,mu)
!         k=0
!         k_div=0
!         Delta_t_bis=Delta_t
!         !call this%compute_molalities() !> we change units to compute activities
!         do
!             do !> loop until convergence is reached
!             !> We apply Newton method to compute aqueous concentrations
!                 call this%Newton_EfI_rk_eq_kin_aq_anal_ideal_opt2(c_tilde,rk_tilde,mix_ratio_Rk,Delta_t_bis,conc_nc,niter,CV_flag)
!             !> We check convergence
!                 if (CV_flag.eqv..false.) then !> no CV
!                     if (mu<1d0) then
!                         mu=mu+0.25 !> we increase Newton initialisation parameter
!                     else
!                         mu=0d0 !> we reset Newton initialisation parameter
!                         k_div=k_div+1 !> we increase time step division counter
!                         Delta_t_bis=Delta_t_bis/2d0 !> we reduce time step
!                     end if
!                     call this%initialise_c1_aq_iterative_method(c1_old,c1,mu) !> we re-initialise Newton method
!                 else 
!                     k=k+1
!                     exit
!                 end if
!             end do
!             if (abs(2d0**(k_div)-k)<this%solid_chemistry%reactive_zone%CV_params%zero) then
!                 exit
!             end if
!         end do
!     !> We compute component concentrations
!         !conc_comp=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,conc_nc)
! !< Post-process
!     deallocate(c1)
! end subroutine