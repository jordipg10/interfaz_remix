!> Computes aqueous species concentrations after iteration of WMA using Euler implicit in chemical reactions for kinetic system
!! We assume all primary species are aqueous
!! The Jacobians are computed analytically
!! We DO NOT apply lumping to the mixing ratios of kinetic reaction rates
!! We average the kinetic reaction rates of this target water

subroutine reactive_mixing_iter_EI_kin_anal_ideal_opt2(this,c1_old,c_hat,&
    mix_ratio_r_old,mix_ratio_r_new,Delta_t,theta,conc_nc,conc_comp,c1_downstream)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, initialise_iterative_method
    use vectors_m, only: inf_norm_vec_real
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
    real(kind=8), intent(in) :: c1_old(:) !> primary concentrations at previous time step
    !real(kind=8), intent(in) :: c_tilde(:) !> variable activity species concentrations after mixing
    !real(kind=8), intent(in) :: rk_tilde_old(:) !> old kinetic reaction rates after mixing (de momento no se usa)
    !real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rate contributions after mixing, ordered in upstream & downstream respectively
    real(kind=8), intent(in) :: c_hat(:) !> new kinetic reaction rates after mixing (de momento no se usa)
    real(kind=8), intent(in) :: mix_ratio_r_old !> mixing ratio of kinetic reaction amount in this target
    real(kind=8), intent(in) :: mix_ratio_r_new !> mixing ratio of kinetic reaction amount in this target
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(inout) :: conc_nc(:) !> variable activity species concentrations (already allocated)
    real(kind=8), intent(inout) :: conc_comp(:) !> component concentrations (already allocated)
    real(kind=8), intent(in), optional :: c1_downstream(:) !> primary concentrations of closest downstream target water
!> Variables
    real(kind=8), allocatable :: c1(:) !> primary concentrations
    !real(kind=8), allocatable :: u_tilde(:) !> component concentrations after mixing
    real(kind=8), allocatable :: u_hat(:) !> component concentrations after mixing
    !real(kind=8), allocatable :: u_rk_tilde(:) !> complementary reaction part of components after mixing (llamalo de otra forma)
    integer(kind=4) :: n_p !> number of primary species
    integer(kind=4) :: n_p_aq !> number of primary aqueous species
    integer(kind=4) :: n_v !> number of variable activity species
    integer(kind=4) :: k_div !> counter of time step divisions
    integer(kind=4) :: k !> counter of time steps
    integer(kind=4) :: niter !> number of Newton iterations
    real(kind=8) :: mu !> Newton initialisation parameter
    real(kind=8) :: Delta_t_bis !> reduced time step
    logical :: CV_flag !> convergence flag
!> Pre-process
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
    n_v=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    !> We get primary species concentrations
    allocate(c1(n_p))
    c1=conc_nc(1:n_p)
    !> Newton initialisation parameter
        mu=0d0 
    !> Update old kinetic reaction rates
    !call this%update_rk_old() !> we update old values of kinetic reaction rates
    !call this%solid_chemistry%update_rk_old() !> we update old values of kinetic reaction rates
    !> We compute component concentrations after mixing
        !u_tilde=this%compute_u_mix(c_tilde)
    !> We compute complementary reaction part of components after mixing (llamalo de otra forma)
        !u_r_tilde=Delta_t*matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,rk_tilde)
    !c_hat=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,c_hat)
!> Process
    !> Initial guess variable activity concentrations
         if (present(c1_downstream)) then
            conc_nc(1:n_p)=c1_downstream(1:n_p)                                !< Use downstream target water concentrations as initial guess
        else
            call initialise_iterative_method(c1_old,c1,mu,conc_nc(1:n_p)) !> we compute initial guess for primary concentrations
        end if
        call this%set_conc_var_act_species(conc_nc) !> we set initial guess aqueous primary species concentrations
        !call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p),conc_nc(n_p+1:n_v)) !> we compute initial guess secondary variable activity concentrations
        call this%set_act_aq_species() !> we set activity of aqueous species
        k=0
        k_div=0
        Delta_t_bis=Delta_t
        !R_tilde=rk_tilde*Delta_t
        !call this%compute_molalities() !> we change units to compute activities
        do
            do !> loop until convergence is reached
            !> We apply Newton method to compute aqueous concentrations
                call this%Newton_EI_kin_anal_ideal_opt2(c_hat,mix_ratio_r_old,&
                    mix_ratio_r_new,Delta_t_bis,theta,conc_nc,niter,&
                    CV_flag)
            !> We check convergence
                if (CV_flag.eqv..false.) then !> no CV
                    if (mu<1d0) then
                        mu=mu+0.25 !> we increase Newton initialisation parameter
                    else
                        mu=0d0 !> we reset Newton initialisation parameter
                        k_div=k_div+1 !> we increase time step division counter
                        if (k_div>this%solid_chemistry%reactive_zone%CV_params%k_div_max) then
                            error stop "Too many time step divisions"
                        end if
                        Delta_t_bis=Delta_t_bis/2d0 !> we reduce time step
                    end if
                    call initialise_iterative_method(c1_old,c1,mu,conc_nc(1:n_p)) !> we compute initial guess for primary concentrations
                    call this%set_conc_var_act_species(conc_nc) !> we set initial guess aqueous primary species concentrations
                    !call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p),conc_nc(n_p+1:n_v)) !> we compute initial guess secondary variable activity concentrations
                    call this%set_act_aq_species() !> we set activity of aqueous species
                else 
                    !call this%update_rk_old() !< Update aqueous kinetic rates
                    !call this%solid_chemistry%update_rk_old() !< Update solid kinetic rates
                    !!> faltan gases
                    k=k+1
                    exit
                end if
            end do
            if (abs(2d0**(k_div)-k)<this%solid_chemistry%reactive_zone%CV_params%zero) then
                exit
            end if
        end do
        !> We compute component concentrations (trivial without equilibrium reactions)
        conc_comp=conc_nc
!< Post-process
    deallocate(c1)
end subroutine