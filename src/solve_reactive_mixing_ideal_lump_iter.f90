!> @brief Performs one time iteration of the reactive mixing loop (ideal, lumped)
!> @details Extracted from the main time loop in solve_reactive_mixing_ideal_lump.
!> Handles updating old attributes, conservative mixing, reactive mixing solver,
!> derived quantities, solid/gas chemistry, output writing, and chemistry snapshots
!> for a single time step k.
subroutine solve_reactive_mixing_ideal_lump_iter(this, k, Delta_t, theta_r, time, &
    num_tar_wat, ind_tar_wat, &
    n_p_cache, n_v_cache, n_eq_cache, &
    has_minerals_flag, has_gases_flag, &
    mix_conc_indices, mixing_ratios_conc, &
    mixing_ratios_R, &
    lumped_lambdas, &
    conc_comp, iunit_mb, iunit_celldbg, &
    p_solver, p_mixer)
    use chemistry_m, only: chemistry_c, reset_wat_indices_zones
    use aqueous_chemistry_m, only: aqueous_chemistry_c, write_celldbg_cell, &
        reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2, &
        compute_c_mix_by_zones_lump, within_zone_dbg_step
    use arrays_m, only: real_array_c, int_array_c
    implicit none

!> Arguments
    class(chemistry_c) :: this
    integer(kind=4), intent(in) :: k                        !< Current time step index
    real(kind=8), intent(in) :: Delta_t                     !< Time step size
    real(kind=8), intent(in) :: theta_r                     !< Time-weighting parameter
    real(kind=8), intent(in) :: time                        !< Current simulation time (mass-balance header)
    integer(kind=4), intent(in) :: num_tar_wat              !< Number of target waters
    integer(kind=4), intent(in) :: ind_tar_wat(:)           !< Indices of target waters
    integer(kind=4), intent(in) :: n_p_cache(:)             !< Cached num_prim_species per target water
    integer(kind=4), intent(in) :: n_v_cache(:)             !< Cached num_var_act_species per target water
    integer(kind=4), intent(in) :: n_eq_cache(:)            !< Cached num_eq_reactions per target water
    logical, intent(in) :: has_minerals_flag(:)             !< Whether each target water has minerals
    logical, intent(in) :: has_gases_flag(:)                !< Whether each target water has gases
    class(int_array_c), intent(in) :: mix_conc_indices      !< Mixing concentration indices
    class(real_array_c), intent(in) :: mixing_ratios_conc   !< Mixing ratios for concentrations
    class(real_array_c), intent(in) :: mixing_ratios_R      !< Reaction mixing ratios (mix_conc_indices layout)
    real(kind=8), intent(in) :: lumped_lambdas(:)           !< Pre-computed lumped lambdas
    real(kind=8), intent(inout) :: conc_comp(:)             !< Work buffer for component concentrations
    integer(kind=4), intent(in) :: iunit_mb                 !< Mass-balance output file unit (already open)
    integer(kind=4), intent(in) :: iunit_celldbg            !< Per-cell [CELL-DBG] trace file unit (already open)
    procedure(reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2) :: p_solver  !< Reactive mixing solver
    procedure(compute_c_mix_by_zones_lump), pointer, intent(in) :: p_mixer

!> Local variables
    integer(kind=4) :: i, j
    integer(kind=4) :: idx
    integer(kind=4) :: n_p, n_v, n_eq
    real(kind=8), allocatable :: c_mix(:)
    real(kind=8), allocatable :: c_hat(:)
    real(kind=8), allocatable :: conc_nc(:)
    real(kind=8), allocatable :: mal_residual(:)
    real(kind=8) :: mal_max
    real(kind=8), allocatable :: mb_lhs(:), mb_transp(:), mb_react(:)   !< domain-summed mass-balance terms
    real(kind=8), allocatable :: mb_lhs_c(:), mb_transp_c(:), mb_react_c(:) !< per-cell mass-balance terms
    integer(kind=4) :: n_v_ref, idx_ref                                !< reference cell dimension/index
    integer(kind=4) :: celldbg_sel(10), n_celldbg_sel, isel            !< [CELL-DBG] 10 equidistant target waters

    !> Update old attributes (copies new conc from previous step into conc_old)
    within_zone_dbg_step = k !< current step for the mixer debug prints (step-gated traces)
    do i = 1, this%num_waters
        call this%waters(i)%update_attributes() !< Shift all old attributes
    end do
    !> Domain-integrated mass-balance accumulators (reference cell = first target water)
    n_v_ref = n_v_cache(1)
    idx_ref = this%tar_wat_indices(ind_tar_wat(1))
    allocate(mb_lhs(n_v_ref), mb_transp(n_v_ref), mb_react(n_v_ref))
    mb_lhs = 0d0; mb_transp = 0d0; mb_react = 0d0
    !> [CELL-DBG] pick up to 10 target waters equidistant along the mesh (endpoints included)
    n_celldbg_sel = min(10, this%num_target_waters)
    do isel = 1, n_celldbg_sel
        if (n_celldbg_sel == 1) then
            celldbg_sel(isel) = this%tar_wat_indices(1)
        else
            celldbg_sel(isel) = this%tar_wat_indices( &
                nint(1d0 + dble(isel-1)*dble(this%num_target_waters-1)/dble(n_celldbg_sel-1)) )
        end if
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
        !> Conservative mixing using the selected mixing procedure
        associate(mc => mix_conc_indices%cols(ind_tar_wat(i)), &
                  rc_col => mixing_ratios_conc%cols(ind_tar_wat(i)))
            call p_mixer(this%waters, this%zone_waters(ind_tar_wat(i),:), &
                mc%col_1(1:mc%dim-2), rc_col%col_1(1:mc%dim-2), &
                lumped_lambdas(ind_tar_wat(i)), &
                this%wat_indices_zones, this%ext_waters_indices, Delta_t, theta_r, c_hat, c_mix, .true.)
        end associate
        !> [NAN-DBG] one-shot diagnostic: catch the first cell whose mixed concentration is
        !> non-finite (NaN/Inf) so the failing time step and cell index are identified.
        block
            logical, save :: nan_dbg_done=.false.
            if ((.not. nan_dbg_done) .and. any(c_mix(1:n_v) /= c_mix(1:n_v))) then
                nan_dbg_done=.true.
                print '(A,I0,A,I0,A,I0)', "[NAN-DBG] step=", k, " cell idx=", idx, " n_v=", n_v
                print '(A,*(ES12.4))', "[NAN-DBG] c_mix=", c_mix(1:n_v)
                print '(A,*(ES12.4))', "[NAN-DBG] conc_nc_old=", this%waters(idx)%get_conc_old_var_act_new_rz()
                print '(A,*(ES12.4))', "[NAN-DBG] conc_nc_cur=", this%waters(idx)%get_conc_var_act()
            end if
        end block
        !> Getters (get_conc_var_act)
        conc_nc = this%waters(idx)%get_conc_var_act()
        !> Reactive mixing solver
        call this%waters(idx)%set_conc_hat(c_hat) !> aqueous blocks of the Newton speciation input, for the conc_hat attribute
        call p_solver(this%waters(idx), this%waters(idx)%get_conc_old_old_prim_new_rz(), c_hat, &
            lumped_lambdas(ind_tar_wat(i)), &
            Delta_t, theta_r, conc_nc, conc_comp)
        !> Mass action law check for equilibrium reactions (ideal solution)
        if (n_eq > 0) then
            block
                real(kind=8), allocatable :: log_gamma_mal(:)
                log_gamma_mal = this%waters(idx)%compute_log_act_coeffs_var_act_ideal()
                mal_residual = matmul(this%waters(idx)%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star, &
                    log10(conc_nc(1:n_p)) + log_gamma_mal(1:n_p)) + &
                    this%waters(idx)%logK_star - &
                    log10(conc_nc(n_p+1:n_p+n_eq)) - log_gamma_mal(n_p+1:n_p+n_eq)
                deallocate(log_gamma_mal)
            end block
            mal_max = maxval(abs(mal_residual))
            if (mal_max > 1d-5) then
                print *, "WARNING: mass action law violated after p_solver, water", idx, &
                    " time step", k, " max|residual|=", mal_max
                print *, "  n_p=", n_p, " n_eq=", n_eq, " n_v=", n_v
                print *, "  c1 (primary)=", conc_nc(1:n_p)
                print *, "  c2 (secondary)=", conc_nc(n_p+1:n_p+n_eq)
                print *, "  logK_star=", this%waters(idx)%logK_star
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
        call this%waters(idx)%compute_saturation_mins() !< saturation indices of all minerals at the end-of-step state
        call this%waters(idx)%compute_salinity()
        call this%waters(idx)%compute_ionic_strength()
        !> Equilibrium reaction amounts
        if (n_eq > 0) then
            call this%waters(idx)%compute_Re(c_hat(1:n_v), Delta_t, &
                lumped_lambdas(ind_tar_wat(i)), theta_r)
        end if
        !> [CELL-DBG] intermediate results for 10 target waters equidistant along the mesh
        if (k<=10 .and. any(celldbg_sel(1:n_celldbg_sel)==idx)) then
            call write_celldbg_cell(this%waters(idx), iunit_celldbg, k, idx, n_v, c_mix, c_mix, c_hat)
        end if
        !> Domain-integrated mass-balance contribution of this cell (lumped: c_hat = c_mix)
        call this%waters(idx)%get_react_mix_mass_bal_terms(c_mix(1:n_v), c_hat(1:n_v), &
            lumped_lambdas(ind_tar_wat(i)), Delta_t, theta_r, mb_lhs_c, mb_transp_c, mb_react_c)
        if (size(mb_lhs_c) == n_v_ref) then
            mb_lhs = mb_lhs + mb_lhs_c
            mb_transp = mb_transp + mb_transp_c
            mb_react = mb_react + mb_react_c
        end if
        deallocate(mb_lhs_c, mb_transp_c, mb_react_c)
        !> Solid/gas chemistry
        if (has_minerals_flag(i)) then
            call this%waters(idx)%solid_chemistry%compute_mass_bal_mins(Delta_t)
            call this%waters(idx)%solid_chemistry%compute_conc_minerals_iter(Delta_t,lumped_lambdas(ind_tar_wat(i)))
        end if
        if (has_gases_flag(i)) then
            !> gas update of the TARGET water idx (ind_tar_wat(i) is its position in the target list, not an index into this%waters)
            call this%waters(idx)%gas_chemistry%compute_conc_gases_iter(Delta_t, &
                this%waters(idx)%volume, &
                [this%waters(idx)%rxn_rate_eq, this%waters(idx)%rxn_rate_kin_avg])
            call this%waters(idx)%gas_chemistry%compute_vol_gas_species_conc()
            call this%waters(idx)%gas_chemistry%compute_log_act_coeffs_gases()
        end if
        !> Accumulate reaction amounts: kinetic AND equilibrium (update_cum_rxn_amnt adds rxn_amnt_kin and
        !> the rxn_amnt_eq that compute_Re set this step); the abs_cum_react_amnt output reports
        !> cum_rxn_amnt_eq for the aqueous equilibrium reactions and the equilibrium minerals
        call this%waters(idx)%update_cum_rxn_amnt()
        call this%waters(idx)%solid_chemistry%update_cum_rxn_amnt()
        if (has_gases_flag(i)) call this%waters(idx)%gas_chemistry%update_cum_rxn_amnt()
        deallocate(c_mix, c_hat, conc_nc)
    end do
    !> Write this time step's domain-integrated mass balance
    call this%write_react_mix_mass_bal_dom(idx_ref, iunit_mb, time, n_v_ref, mb_lhs, mb_transp, mb_react)
    deallocate(mb_lhs, mb_transp, mb_react)
    !> Sum equilibrium reaction rates at interfaces if mix_by_zones_flag is set
    if (this%mix_by_zones_flag) then
        do i=1,num_tar_wat
            !print '(A,I0)', "[SUMEQ-DBG] Summing equilibrium reaction rates for target water ", this%tar_wat_indices(ind_tar_wat(i))
            call this%waters(this%tar_wat_indices(ind_tar_wat(i)))%sum_eq_rxn_rates(this%zone_waters(&
                this%tar_wat_indices(ind_tar_wat(i)),:))
            !print '(A,I0)', "[SUMEQ-DBG] Finished summing equilibrium reaction rates for target water ", this%tar_wat_indices(ind_tar_wat(i))
        end do
    end if
    !> Relink target waters to reactive zones if the zone set changed this step
    call reset_wat_indices_zones(this)
end subroutine solve_reactive_mixing_ideal_lump_iter
