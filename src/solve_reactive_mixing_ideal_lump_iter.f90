!> @brief Performs one time iteration of the reactive mixing loop (ideal, lumped)
!> @details Extracted from the main time loop in solve_reactive_mixing_ideal_lump.
!> Handles updating old attributes, conservative mixing, reactive mixing solver,
!> derived quantities, solid/gas chemistry, output writing, and chemistry snapshots
!> for a single time step k.
subroutine solve_reactive_mixing_ideal_lump_iter(this, k, Delta_t, theta_r, &
    num_tar_wat, ind_tar_wat, &
    n_p_cache, n_v_cache, n_eq_cache, &
    has_minerals_flag, has_gases_flag, &
    mix_conc_indices, mixing_ratios_conc, &
    lumped_lambdas, &
    all_conc_old, all_ind_aq_sp, &
    conc_comp, &
    p_solver)
    use chemistry_m, only: chemistry_c
    use aqueous_chemistry_m, only: aqueous_chemistry_c, &
        reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2, &
        compute_c_mix_global
    use arrays_m, only: real_array_c, int_array_c
    implicit none

!> Arguments
    class(chemistry_c) :: this
    integer(kind=4), intent(in) :: k                        !< Current time step index
    real(kind=8), intent(in) :: Delta_t                     !< Time step size
    real(kind=8), intent(in) :: theta_r                     !< Time-weighting parameter
    integer(kind=4), intent(in) :: num_tar_wat              !< Number of target waters
    integer(kind=4), intent(in) :: ind_tar_wat(:)           !< Indices of target waters
    integer(kind=4), intent(in) :: n_p_cache(:)             !< Cached num_prim_species per target water
    integer(kind=4), intent(in) :: n_v_cache(:)             !< Cached num_var_act_species per target water
    integer(kind=4), intent(in) :: n_eq_cache(:)            !< Cached num_eq_reactions per target water
    logical, intent(in) :: has_minerals_flag(:)             !< Whether each target water has minerals
    logical, intent(in) :: has_gases_flag(:)                !< Whether each target water has gases
    class(int_array_c), intent(in) :: mix_conc_indices      !< Mixing concentration indices
    class(real_array_c), intent(in) :: mixing_ratios_conc   !< Mixing ratios for concentrations
    real(kind=8), intent(in) :: lumped_lambdas(:)           !< Pre-computed lumped lambdas
    real(kind=8), intent(in) :: all_conc_old(:,:)           !< Global cache of old concentrations
    integer(kind=4), intent(in) :: all_ind_aq_sp(:,:)       !< Global cache of aq species indices
    real(kind=8), intent(inout) :: conc_comp(:)             !< Work buffer for component concentrations
    procedure(reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2) :: p_solver  !< Reactive mixing solver

!> Local variables
    integer(kind=4) :: i, j
    integer(kind=4) :: idx
    integer(kind=4) :: n_p, n_v, n_eq
    integer(kind=4) :: num_mix_loc
    real(kind=8), allocatable :: c_mix(:)
    real(kind=8), allocatable :: conc_nc(:)
    real(kind=8), allocatable :: mal_residual(:)
    real(kind=8) :: mal_max

    !> Update old attributes
    do i = 1, this%num_target_waters
        call this%waters(this%tar_wat_indices(i))%update_old_attributes()
    end do
    !> Target waters loop
    do i = 1, num_tar_wat
        idx = this%tar_wat_indices(ind_tar_wat(i))
        n_p = n_p_cache(i)
        n_v = n_v_cache(i)
        n_eq = n_eq_cache(i)
        !> Data preparation
        if (mix_conc_indices%cols(ind_tar_wat(i))%col_1(1) /= idx) then
            error stop "Target waters not in the right order"
        end if
        num_mix_loc = mix_conc_indices%cols(ind_tar_wat(i))%dim - 3
        !> Conservative mixing
        call compute_c_mix_global(this%waters(mix_conc_indices%cols(ind_tar_wat(i))%col_1(1)), &
            all_conc_old, all_ind_aq_sp, &
            mix_conc_indices%cols(ind_tar_wat(i))%col_1(2:mix_conc_indices%cols(ind_tar_wat(i))%dim-2), num_mix_loc, &
            mixing_ratios_conc%cols(ind_tar_wat(i))%col_1, c_mix)
        !> Getters (get_conc_nc)
        conc_nc = this%waters(idx)%get_conc_nc()
        !> Reactive mixing solver
        call p_solver(this%waters(idx), this%waters(idx)%get_c1_old_old(), c_mix(1:n_v), &
            lumped_lambdas(ind_tar_wat(i)), &
            lumped_lambdas(ind_tar_wat(i)), &
            Delta_t, theta_r, conc_nc, conc_comp(1:n_p))
        !> Mass action law check for equilibrium reactions (ideal solution)
        if (n_eq > 0) then
            mal_residual = matmul(this%waters(idx)%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star, &
                log10(conc_nc(1:n_p))) + &
                this%waters(idx)%solid_chemistry%reactive_zone%speciation_alg%logK_star - &
                log10(conc_nc(n_p+1:n_p+n_eq))
            mal_max = maxval(abs(mal_residual))
            if (mal_max > 1d-5) then
                print *, "WARNING: mass action law violated after p_solver, water", idx, &
                    " time step", k, " max|residual|=", mal_max
                print *, "  n_p=", n_p, " n_eq=", n_eq, " n_v=", n_v
                print *, "  c1 (primary)=", conc_nc(1:n_p)
                print *, "  c2 (secondary)=", conc_nc(n_p+1:n_p+n_eq)
                print *, "  logK_star=", this%waters(idx)%solid_chemistry%reactive_zone%speciation_alg%logK_star
                do j = 1, n_eq
                    print *, "  eq reaction", j, " residual=", mal_residual(j), &
                        " log10(c2)=", log10(conc_nc(n_p+j)), &
                        " expected=", log10(conc_nc(n_p+j)) - mal_residual(j)
                    print *, "    Se_nc_1_star row=", &
                        this%waters(idx)%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star(j,:)
                end do
                error stop "Mass action law violation after p_solver"
            end if
            deallocate(mal_residual)
        end if
        !> Derived quantities (pH, salinity, ionic strength)
        call this%waters(idx)%compute_pH()
        call this%waters(idx)%compute_salinity()
        call this%waters(idx)%compute_ionic_strength()
        !> Equilibrium reaction amounts
        if (n_eq > 0 .and. this%waters(idx)%indices_rk%num_cols > 0) then
            call this%waters(idx)%compute_Re_kin(c_mix(n_p+1:n_v), Delta_t, &
                lumped_lambdas(ind_tar_wat(i)))
        else if (n_eq > 0) then
            call this%waters(idx)%compute_Re(c_mix(1:n_v), Delta_t, &
                lumped_lambdas(ind_tar_wat(i)))
        end if
        !> Solid/gas chemistry
        if (has_minerals_flag(i)) then
            call this%waters(idx)%solid_chemistry%compute_mass_bal_mins(Delta_t)
            call this%waters(idx)%solid_chemistry%compute_conc_minerals_iter(Delta_t)
        end if
        if (has_gases_flag(i)) then
            call this%waters(idx)%gas_chemistry%compute_conc_gases_iter(Delta_t, &
                this%waters(ind_tar_wat(i))%volume, &
                [this%waters(ind_tar_wat(i))%re_mean, this%waters(ind_tar_wat(i))%rk_mean])
            call this%waters(idx)%gas_chemistry%compute_vol_gas_species_conc()
            call this%waters(idx)%gas_chemistry%compute_log_act_coeffs_gases()
        end if
        !> Accumulate reaction amounts
        this%waters(idx)%Rk_accum = this%waters(idx)%Rk_accum + this%waters(idx)%Rk
        this%waters(idx)%solid_chemistry%Rk_accum = &
            this%waters(idx)%solid_chemistry%Rk_accum + this%waters(idx)%solid_chemistry%Rk
        deallocate(c_mix, conc_nc)
    end do
end subroutine solve_reactive_mixing_ideal_lump_iter
