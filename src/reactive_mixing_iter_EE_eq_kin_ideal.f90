!> Computes aqueous variable activity species concentrations and aqueous component concentrations after iteration of WMA-EE method for equilibrium-kinetic chemical system
subroutine reactive_mixing_iter_EE_eq_kin_ideal(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, initialise_iterative_method
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    real(kind=8), intent(in) :: c1_old(:) !> primary species concentrations before previous time step
    !real(kind=8), intent(in) :: c2v_ig(:)
    !real(kind=8), intent(in) :: c_tilde(:) !> variable activity species concentrations after mixing
    !real(kind=8), intent(in) :: R_tilde(:) !> reaction contributions after mixing
    real(kind=8), intent(in) :: c_hat(:) !> mixing contribution to variable activity species concentrations
    real(kind=8), intent(in) :: mix_ratio_r !> mixing ratio for reaction contributions
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(inout) :: conc_nc(:) !> variable activity species concentrations after time step
    real(kind=8), intent(inout) :: conc_comp(:) !> component concentrations after time step
    !real(kind=8), intent(in) :: Delta_t !> time step
!> Variables
    real(kind=8), allocatable :: u_tilde(:)
    real(kind=8), allocatable :: conc_comp_react(:)
    !real(kind=8), allocatable :: conc_comp(:)
    real(kind=8), allocatable :: c1(:)
    real(kind=8), allocatable :: c1_ig(:)
    real(kind=8), allocatable :: rk_new(:)
    integer(kind=4) :: niter !> number of iterations in Newton speciation
    integer(kind=4) :: k !> counter of time steps
    integer(kind=4) :: k_div !> counter of time step divisions
    logical :: CV_flag !> convergence flag
    real(kind=8) :: mu=0d0 !> Newton initialistaion parameter
    real(kind=8) :: Delta_t_aux !> time step (auxiliary)
!> Pre-process
    allocate(conc_comp_react(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    !allocate(conc_comp(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    allocate(c1_ig(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    allocate(rk_new(this%indices_rk%num_cols))
    c1=conc_nc(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species) !> primary species concentrations after previous time step
    Delta_t_aux=Delta_t !> we set default value for time step
    k=0 !> we initialise counter of time steps
    k_div=0 !> we initialise counter of time step divisions
!> Process
    !> Components concentrations
    !u_tilde=this%compute_u_mix(c_tilde) !> transport part: concentrations of variable activity species after mixing
    conc_comp=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,c_hat) !> we compute mixing term
    ! !> Update old kinetic reaction rates
    ! call this%update_rk_old() !> we update old values of kinetic reaction rates
    ! call this%solid_chemistry%update_rk_old() !> we update old values of kinetic reaction rates
!> We initialise primary concentrations for Newton speciation
    call initialise_iterative_method(c1_old,c1,mu,c1_ig)
    !call this%set_conc_aq_prim_species(conc_nc(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)) !> we set initial guess aqueous primary species concentrations
    !call this%compute_log_act_coeffs_aq_chem()
    !log_gamma_var_act=this%get_log_gamma_var_act()
    do
        !k=0
        !k_div=0
        !Delta_t_bis=Delta_t
        call this%compute_rk_new(rk_new) !> we compute new values of kinetic reaction rates
        !call this%compute_react_term_EE_eq_kin(Delta_t_aux,rk_tilde,conc_comp_react) !> chemical part
        !conc_comp_react=Delta_t*rk_tilde,theta) !> chemical part
        conc_comp_react=Delta_t*mix_ratio_r*matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,&
            rk_new) !> we add reactions contribution
        conc_comp=conc_comp+conc_comp_react !> we add the reactions contribution
    !> Speciation for variable activity species concentrations
        !> Loop until speciation converges
        do
        !> We compute variable activity concentrations from component concentrations
            call this%compute_c_nc_from_u_Newton_ideal(c1_ig,conc_comp,conc_nc,niter,CV_flag)
        !> We check convergence
            if (.not. CV_flag) then !> NO CV
                if (mu<1d0) then
                    mu=mu+0.25
                    !call this%initialise_iterative_method(c1_old,c1,mu,c1_ig)
                else
                    mu=0d0
                    k_div=k_div+1
                    if (k_div>this%solid_chemistry%reactive_zone%CV_params%k_div_max) error stop "Too many time step divisions"
                    Delta_t_aux=Delta_t_aux/2d0 !> we reduce time step
                    conc_comp=conc_comp-conc_comp_react !> we reset component concentrations
                    !call this%initialise_iterative_method(c1_old,c1,mu,c1_ig)
                    exit
                    !error stop "Newton speciation does not converge"
                end if
            else
                c1_ig=conc_nc(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species) !> we update primary species concentrations
                k=k+1 !> we update counter of time steps
                mu=0d0 !> we reset mu
                exit
            end if
        end do
        if (abs(2d0**(k_div)-k)<this%solid_chemistry%reactive_zone%CV_params%zero) then
            exit
        end if
    end do
        !> We check results
            !call this%check_conc_var_act_species(conc_nc,conc_comp)
            !call this%check_act_aq_species()
!> Post-process
    deallocate(conc_comp_react,c1)
end subroutine