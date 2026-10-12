!> \file equilibrate_water.f90
!> \brief Speciate (equilibrate) a water with its reactive zone (non-ideal solution).
!> \details Copy of mixing_iter_comp WITHOUT the "unchanged primaries" early-exit: the Newton
!> speciation ALWAYS runs. Used by the per-reactive-zone mixing (mix_within_zone), where the water
!> handed in is a displaced (transported) water whose current state is NOT yet at equilibrium with
!> the (new) reactive zone, so the early-exit of mixing_iter_comp would wrongly skip speciation.
subroutine equilibrate_water(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, initialise_iterative_method
    implicit none
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object (speciated in place)
    real(kind=8), intent(in) :: c1_old(:) !> primary species concentrations from previous time step (Newton initial guess)
    real(kind=8), intent(in) :: c_hat(:) !> mixed variable activity species concentrations
    real(kind=8), intent(in) :: mix_ratio_r !> mixing ratio for kinetic reactions (unused; kept for signature parity)
    real(kind=8), intent(in) :: Delta_t !> time step (unused; kept for signature parity)
    real(kind=8), intent(in) :: theta !> time weighting factor (unused; kept for signature parity)
    real(kind=8), intent(inout) :: conc_nc(:) !> variable activity species concentrations (output, already allocated)
    real(kind=8), intent(inout) :: conc_comp(:) !> component concentrations (already allocated)

    integer(kind=4) :: niter !> Newton iterations performed
    integer(kind=4) :: n_p_aq !> number of aqueous primary species
    integer(kind=4) :: n_p !> number of primary species
    integer(kind=4) :: n_eq !> number of equilibrium reactions
    logical :: CV_flag !> Newton convergence flag
    real(kind=8) :: mu !> Newton initialisation parameter in [0,1]
    real(kind=8), allocatable :: c1(:) !> current primary species concentrations
    real(kind=8), allocatable :: c1_ig(:) !> Newton initial guess for primary species
    real(kind=8), allocatable :: c2v_ig(:) !> initial guess secondary variable activity species
    real(kind=8), allocatable :: log_act_coeffs(:) !> log10 activity coefficients of variable activity species

    mu=0d0
    associate(spec_alg => this%solid_chemistry%reactive_zone%speciation_alg)
    n_p_aq=spec_alg%num_aq_prim_species
    n_p=spec_alg%num_prim_species
    n_eq=spec_alg%num_eq_reactions
    c1=this%get_c1()
    !> Component concentrations from the mixed species: u = U*c_hat
    conc_comp=matmul(spec_alg%comp_mat,c_hat)
    !> ALWAYS speciate (no early-exit): Newton with adaptive initialisation mu = 0 -> 0.25 -> ... -> 1
    allocate(c1_ig(n_p))
    allocate(c2v_ig(n_eq))
    allocate(log_act_coeffs(spec_alg%num_var_act_species))
    do
        call initialise_iterative_method(c1_old,c1,mu,c1_ig)
        call this%set_conc_aq_prim_species(c1_ig(1:n_p_aq))
        call this%compute_log_act_coeffs_aq_chem()
        log_act_coeffs=this%get_log_gamma_var_act()
        call this%compute_c2v_from_c1_ideal(c1_ig,log_act_coeffs,c2v_ig)
        call this%compute_c_nc_from_u_Newton(c1_ig,c2v_ig,conc_comp,conc_nc,niter,CV_flag)
        if (CV_flag.eqv..false.) then
            if (mu<1d0) then
                mu=mu+0.25d0 !> retry with more weight on the current mixed state
            else
                error stop "Newton speciation does not converge"
            end if
        else
            exit
        end if
    end do
    call this%compute_pH()
    call this%compute_salinity()
    call this%compute_ionic_strength()
    end associate
end subroutine
