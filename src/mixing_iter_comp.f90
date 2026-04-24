!> Computes component and variable activity concentrations after iteration of WMA method for equilibrium chemical system
!> We assume all primary species are aqueous
subroutine mixing_iter_comp(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
    use aqueous_chemistry_m,  only: aqueous_chemistry_c, initialise_iterative_method
    use vectors_m, only: inf_norm_vec_real
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
    real(kind=8), intent(in) :: c1_old(:)
    !real(kind=8), intent(in) :: c2v_ig(:)
    real(kind=8), intent(in) :: c_hat(:)
    !real(kind=8), intent(in) :: rk_tilde(:)
    real(kind=8), intent(in) :: mix_ratio_r
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta !> weighting factor for implicitness
    real(kind=8), intent(inout) :: conc_nc(:)
    real(kind=8), intent(inout) :: conc_comp(:) !> concentration primary species for iteration
    !real(kind=8), intent(in), optional :: porosity
    
!> Variables
    integer(kind=4) :: niter !> number of iterations in Newton speciation
    integer(kind=4) :: n_p_aq !> number of aqueous primary species
    logical :: CV_flag !> convergence flag
    real(kind=8) :: mu=0d0 !> Newton initialistaion parameter
    real(kind=8), allocatable :: c1(:) !> concentration primary species
    real(kind=8), allocatable :: c1_ig(:) !> initial guess primary species
    real(kind=8), allocatable :: c2v_ig(:) !> initial guess secondary variable activity species
    real(kind=8), allocatable :: log_act_coeffs(:) !> logarithm (base 10) of activity coefficients of variable activity species
!> Pre-process
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
    c1=this%get_c1()
    allocate(c1_ig(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
    allocate(c2v_ig(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
!> Process  
    !> We compute component concentrations after mixing
        conc_comp=MATMUL(THIS%solid_chemistry%reactive_zone%speciation_alg%comp_mat,c_hat)
        if (inf_norm_vec_real(c_hat(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)-c1)<&
            this%solid_chemistry%reactive_zone%CV_params%zero) then !> chapuza
            conc_nc=this%get_conc_nc()
        else
            !> Loop until speciation converges
            do
            !> We initialise primary concentrations for Newton speciation
                call initialise_iterative_method(c1_old,c1,mu,c1_ig)
                call this%set_conc_aq_prim_species(c1_ig(1:n_p_aq))
                call this%compute_log_act_coeffs_aq_chem()
                log_act_coeffs=this%get_log_gamma_var_act()
                call this%compute_c2v_from_c1_ideal(c1_ig,log_act_coeffs,c2v_ig)
            !> We compute variable activity concentrations from component concentrations
                call this%compute_c_nc_from_u_Newton(c1_ig,c2v_ig,conc_comp,conc_nc,niter,CV_flag)
            !> We check convergence
                if (CV_flag.eqv..false.) then !> NO CV
                    if (mu<1d0) then
                        mu=mu+0.25
                    else
                        error stop "Newton speciation does not converge"
                    end if
                else
                    exit
                end if
            end do
        end if
!> Post-process
    !> We check concentrations
        !call this%check_conc_aq_var_act_species(conc_comp)
        !call this%check_act_aq_species()
        deallocate(c1_ig,c1)
end subroutine