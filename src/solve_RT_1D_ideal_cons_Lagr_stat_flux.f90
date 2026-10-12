!> \file solve_RT_1D_ideal_cons_Lagr_stat_flux.f90
!> \brief Solve 1D reactive transport with ideal speciation, consistent WMA, and Lagrangian particle tracking
!> \details
!>   Solves the 1D reactive transport problem using:
!>   - Ideal speciation (gamma = 1)
!>   - Consistent Water Mixing Approach (WMA)
!>   - Lagrangian particle tracking with stationary flux
!>
!>   Operator splitting decouples transport (advection-dispersion via particle mixing)
!>   from reactions (equilibrium/kinetic chemistry at each target).
!>
!> \param[inout] this Reactive transport object (transient)
!> \param[in] dir  Directory for output files
!> \param[in] root Root name for output file
!> \author Jordi Petchame-Guerrero
!> \date November 2025
subroutine solve_RT_ideal_cons_Lagr_stat_flux_1D(this, dir, root)
    use RT_m, only: RT_1D_transient_c, move_particles_stat_flux_1D, move_particles_stat_flux_EC_1D
    use chemistry_m, only: write_out_step_vars, reset_wat_indices_zones !< Shared per-output-step writer (conc/react_rate/vol_fract); post-displacement zone regrouping
    use aqueous_chemistry_m, only: &
        aqueous_chemistry_c, &
        reactive_mixing_iter_EI_eq_anal_ideal, &
        reactive_mixing_iter_EI_kin_anal_ideal_opt2, &
        compute_c_mix_global_cons, compute_c_mix_by_zones_cons, &
        compute_r_tilde_impl_opt1_idx, compute_r_tilde_impl_opt2_idx, &
        compute_r_tilde_impl_opt3_idx, compute_r_tilde_impl_opt4_idx, &
        reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2, add_limited_increment
    use arrays_m, only: int_array_c
    implicit none

    ! Arguments
    class(RT_1D_transient_c) :: this
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: root

    ! Loop counters
    integer(kind=4) :: i, k, kk

    ! Problem dimensions
    integer(kind=4) :: n_p, n_v, n_eq, n_v_aq
    integer(kind=4) :: num_can_vec, num_non_can_vec
    integer(kind=4) :: out_unit

    ! Index arrays
    integer(kind=4), allocatable :: ind_can_vec(:), ind_non_can_vec(:)
    integer(kind=4), allocatable :: mix_ext_indices(:) !< mixer externals: boundary/recharge waters (minus the Neumann outflow water) + the inflow cell (its parcel is the boundary water each step)
    integer(kind=4), allocatable :: exit_wat_indices(:) !< Neumann outflow boundary water: holds the parcel that left the last cell, grouped in that cell's zone column (empty for a Dirichlet outflow)
    integer(kind=4) :: q !< loop index over the external waters

    ! Time and mixing
    real(kind=8) :: time, Delta_t, lambda_R

    ! Concentration arrays
    real(kind=8), allocatable :: c_hat(:), c_mix(:), conc_comp(:)
    real(kind=8), allocatable :: c_hat_old(:)
    real(kind=8), allocatable :: conc_nc(:), r_tilde(:)
    ! Domain-integrated reactive-mixing mass balance (per variable-activity species)
    integer(kind=4) :: iunit_mb                   !< mass-balance output file unit
    integer(kind=4) :: idx_ref                    !< reference target-water index (species names/dimension)
    integer(kind=4) :: n_v_ref                    !< reference number of variable-activity species
    real(kind=8), allocatable :: mb_lhs(:), mb_transp(:), mb_react(:) !< domain-summed terms
    real(kind=8), allocatable :: mb_lhs_c(:), mb_transp_c(:), mb_react_c(:) !< per-cell terms
    ! Procedure pointers
    procedure(reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2), pointer :: p_solver => null()
    procedure(compute_r_tilde_impl_opt4_idx), pointer :: p_r_tilde_idx => null()
    procedure(compute_c_mix_by_zones_cons), pointer :: p_mixer
    procedure(move_particles_stat_flux_1D), pointer :: p_displace_particles => null()

    ! Cached sizes for allocation optimisation
    integer(kind=4) :: prev_n_v, prev_n_p

    ! CPU profiling
    integer(kind=4) :: iunit_cpu                  !< I/O unit for the CPU profiling report
    integer(kind=4) :: num_profiled_steps         !< Number of profiled time steps (= all of them here)
    character(len=:), allocatable :: cpu_root     !< root + transport/WMA tags for the CPU profiling filename
    real(kind=8) :: cpu_total_t0, cpu_total_t1, cpu_total_elapsed
    real(kind=8) :: cpu_prof_t0, cpu_prof_t1
    real(kind=8) :: cpu_update_old_total, cpu_displace_total, cpu_mix_total
    real(kind=8) :: cpu_data_prep_total, cpu_getters_total, cpu_solver_total
    real(kind=8) :: cpu_derived_total, cpu_eq_rates_total, cpu_solid_gas_total
    real(kind=8) :: cpu_output_total, cpu_snapshot_total

    ! Count equilibrium reactions across all target waters
    if (this%chemistry%mix_by_zones_flag) then
        p_mixer => compute_c_mix_by_zones_cons
    else
        p_mixer => compute_c_mix_global_cons
    end if
    n_eq = 0
    do i = 1, this%chemistry%num_target_waters
        n_eq = n_eq + this%target_waters(i)%aq_chem%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    end do

    select type (this)
    type is (RT_1D_transient_c)
        time = 0d0
        kk = 2
        block
            character(len=:), allocatable :: rtag
            if (this%chemistry%lump_flag) then
                rtag = '_Lagr_lumped'
            else
                rtag = '_Lagr_consistent'
            end if
            open(newunit=out_unit, file=dir//trim(root)//rtag//'_output.out', form="formatted", access="sequential", status="unknown")
        end block

        !--- Solver selection based on chemical system type ---
        if (this%chemistry%chem_syst%num_kin_reacts > 0 .and. n_eq > 0) then
            ! Equilibrium + kinetic reactions
            if (this%chemistry%rk_avg_opt == 1) then
                error stop "rk average option 1 not implemented yet"
            else if (this%chemistry%rk_avg_opt == 2) then
                p_solver => reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2
            else
                error stop "rk average option not implemented yet"
            end if
        else if (n_eq > 0) then
            ! Equilibrium-only reactions
            ! int_method 1 (explicit), 2 (implicit) and 3 (Crank-Nicolson) share this solver: the
            ! scheme enters only through theta_r (r_tilde/c_hat), not the algebraic equilibrium speciation.
            if (this%transport%time_discr%int_method_chem >= 1 .and. this%transport%time_discr%int_method_chem <= 3) then
                p_solver => reactive_mixing_iter_EI_eq_anal_ideal
            else
                error stop "Integration method for equilibrium-only reactions not implemented yet"
            end if
        else
            ! Kinetic-only reactions
            if (this%chemistry%rk_avg_opt == 1) then
                error stop "rk average option 1 not implemented yet"
            else if (this%chemistry%rk_avg_opt == 2) then
                p_solver => reactive_mixing_iter_EI_kin_anal_ideal_opt2
            else
                error stop "rk average option not implemented yet"
            end if
        end if

        ! Downstream reaction rate estimation (common to all chemical system types)
        select case (this%chemistry%r_down_opt)
        case (1)
            p_r_tilde_idx => compute_r_tilde_impl_opt1_idx
        case (2)
            p_r_tilde_idx => compute_r_tilde_impl_opt2_idx
        case (3)
            p_r_tilde_idx => compute_r_tilde_impl_opt3_idx
        case (4)
            p_r_tilde_idx => compute_r_tilde_impl_opt4_idx
        case default
            error stop "Downstream reaction rate estimation option not implemented yet"
        end select

        ! Particle displacement method
        if (this%transport%spatial_discr%targets_flag == 0) then
            p_displace_particles => move_particles_stat_flux_1D
        else if (this%transport%spatial_discr%targets_flag == 1) then
            p_displace_particles => move_particles_stat_flux_EC_1D
        else
            error stop "Stationary flux Lagrangian scheme not implemented yet"
        end if

        ! Extract canonical vectors from mixing ratios matrix
        call this%transport%mixing_ratios_conc%get_can_vec( &
            this%chemistry%CV_params%abs_tol, num_can_vec, ind_can_vec, num_non_can_vec, ind_non_can_vec)

        ! Initialise concentration history (and the previous-step solid chemistry + index arrays of
        ! EVERY water, external ones included: the per-step update_attributes only visits targets)
        do i = 1, this%chemistry%num_waters
            call this%chemistry%waters(i)%set_conc_old()
            call this%chemistry%waters(i)%set_conc_old_old()
            call this%chemistry%waters(i)%solid_chemistry%set_conc_old()
            call this%chemistry%waters(i)%solid_chemistry%set_conc_old_old()
            !call this%chemistry%waters(i)%set_solid_chemistry_old()
        end do

        ! Output headers/blocks are written per step by write_out_step_vars (honours chem_out_options%vars)

        ! Tagged root: transport (Lagr/Euler) + WMA (lumped/consistent), shared by the
        ! mass-balance and CPU-profile filenames ({root}{_Lagr|_Euler}{_lumped|_consistent}_*.out).
        if (this%transport%Lagr_flag) then
            cpu_root = trim(root)//'_Lagr'
        else
            cpu_root = trim(root)//'_Euler'
        end if
        if (this%chemistry%lump_flag) then
            cpu_root = cpu_root//'_lumped'
        else
            cpu_root = cpu_root//'_consistent'
        end if

        ! Domain-integrated mass balance: reference cell (homogeneous reactive zone) + output file
        idx_ref = this%chemistry%tar_wat_indices(ind_non_can_vec(1))
        n_v_ref = this%chemistry%waters(idx_ref)%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
        allocate(mb_lhs(n_v_ref), mb_transp(n_v_ref), mb_react(n_v_ref))
        open(newunit=iunit_mb, file=dir//cpu_root//'_mass_balance.out', status='replace', action='write')
        write(iunit_mb,*) 'Domain-integrated reactive-mixing mass balance (per variable-activity species)'
        write(iunit_mb,*) 'Per time step & species: LHS = sum V*(c^{k+1}-c^k)/Dt, RHS = transport + reaction, residual = LHS - RHS'
        ! Initialise allocation tracking for inner loop
        prev_n_v = -1
        prev_n_p = -1

        ! --- CPU profiling setup ---
        call cpu_time(cpu_total_t0)
        cpu_update_old_total = 0d0; cpu_displace_total = 0d0; cpu_mix_total = 0d0
        cpu_data_prep_total = 0d0;  cpu_getters_total = 0d0;  cpu_solver_total = 0d0
        cpu_derived_total = 0d0;    cpu_eq_rates_total = 0d0; cpu_solid_gas_total = 0d0
        cpu_output_total = 0d0;     cpu_snapshot_total = 0d0
        num_profiled_steps = this%transport%time_discr%Num_time  !< every step is profiled (no sampling/scaling)
        open(newunit=iunit_cpu, file=dir//cpu_root//'_cpu_profile.out', status='replace', action='write')

        !> Write initial state (t=0) before the first time step, honouring chem_out_options%vars
        if (this%chemistry%chem_out_options%time_steps(1)==0) then
            call write_out_step_vars(this%chemistry, out_unit, 0d0)
        end if
        !> Mixer grouping of the waters that are not plain target waters. The inflow cell's parcel is
        !> the boundary water (zone not in the domain): external, mixed directly. With a Neumann
        !> outflow, the outflow boundary water is the parcel that left the last cell: its previous-step
        !> concentrations were equilibrated in that cell's zone, so it joins that zone's within-zone
        !> mix (classified by reset_wat_indices_zones) and is NOT an external. A Dirichlet
        !> outflow water stays fixed and external.
        if (this%transport%BCs%labels(2)==2) then
            exit_wat_indices=[this%chemistry%bd_waters_indices(2)]
        else
            allocate(exit_wat_indices(0))
        end if
        ! mix_ext_indices=[pack(this%chemistry%ext_waters_indices, &
        !     [(all(this%chemistry%ext_waters_indices(q)/=exit_wat_indices), q=1,this%chemistry%num_ext_waters)]), &
        !     this%chemistry%tar_wat_indices(1)]
        !====================================================================
        ! Main time loop
        !====================================================================
        do k = 1, this%transport%time_discr%Num_time
            Delta_t = this%transport%time_discr%get_Delta_t(k)
            time = time + Delta_t

            ! Reset the domain-integrated mass-balance accumulators for this time step
            mb_lhs = 0d0
            mb_transp = 0d0
            mb_react = 0d0

            ! Shift time levels: current -> old, old -> old_old
            call cpu_time(cpu_prof_t0)
            do i = 1, this%chemistry%num_waters
                call this%chemistry%waters(i)%update_attributes() !< Shift all old attributes
            end do
            call cpu_time(cpu_prof_t1)
            cpu_update_old_total = cpu_update_old_total + (cpu_prof_t1 - cpu_prof_t0)

            call cpu_time(cpu_prof_t0)
            call p_displace_particles(this, k)
            call this%introduce_particle(k)
            !> Regroup the target waters by the reactive zone of their previous-step state (the zone of
            !> the cell each parcel comes from, solid_chemistry_old): wat_indices_zones is the
            !> classification compute_c_mix_by_zones uses for the within-zone mixing. The Neumann
            !> outflow water (exited parcel) is classified with them.
            call this%chemistry%set_wat_indices_zones()
            call cpu_time(cpu_prof_t1)
            cpu_displace_total = cpu_displace_total + (cpu_prof_t1 - cpu_prof_t0)

            !----------------------------------------------------------------
            ! Target waters loop (non-canonical vectors only)
            !----------------------------------------------------------------
            do i = 1, num_non_can_vec
                associate( &
                    idx   => ind_non_can_vec(i), &
                    aq    => this%target_waters(ind_non_can_vec(i))%aq_chem, &
                    sc    => this%target_waters(ind_non_can_vec(i))%aq_chem%solid_chemistry, &
                    rz    => this%target_waters(ind_non_can_vec(i))%aq_chem%solid_chemistry%reactive_zone, &
                    sa    => this%target_waters(ind_non_can_vec(i))%aq_chem%solid_chemistry%reactive_zone%speciation_alg, &
                    mc    => this%transport%mix_conc_indices%cols(ind_non_can_vec(i)), &
                    mr_c  => this%transport%mixing_ratios_conc%cols(ind_non_can_vec(i)), &
                    mr_R  => this%transport%mixing_ratios_R%cols(ind_non_can_vec(i)) &
                )

                n_p = sa%num_prim_species
                n_v = sa%num_var_act_species
                n_v_aq = sa%num_aq_var_act_species
                ! Reallocate only when sizes change
                if (n_v /= prev_n_v .or. n_p /= prev_n_p) then
                    if (allocated(c_hat)) deallocate(c_hat, conc_comp, r_tilde)
                    allocate(c_hat(n_v), conc_comp(n_p), r_tilde(n_v))
                    prev_n_v = n_v
                    prev_n_p = n_p
                end if

                ! Validate target water ordering
                if (mc%col_1(1) /= this%chemistry%tar_wat_indices(idx)) then
                    print *, mc%col_1(1), this%chemistry%tar_wat_indices(idx)
                    print *, "Problematic target water", idx
                    error stop "Target waters not in the right order in mixing waters indices"
                end if

                ! Transport mixing: page-4 two-stage per-reactive-zone consistent mix -> c_mix
                ! (stage 2 within-zone mixing + speciation, stage 3 cross-zone aqueous mixing).
                ! The reaction mixing ratios share the concentration stencil (0 for the external waters).
                call cpu_time(cpu_prof_t0)
                call p_mixer(this%chemistry%waters, this%chemistry%zone_waters(idx,:), &
                    mc%col_1(1:mc%dim-2), mr_c%col_1(1:mc%dim-2), mr_R%col_1(1:mc%dim-2), &
                    this%chemistry%wat_indices_zones, this%chemistry%ext_waters_indices, Delta_t, &
                    this%transport%time_discr%theta_r, c_hat_old, c_mix, .true.)
                call cpu_time(cpu_prof_t1)
                cpu_mix_total = cpu_mix_total + (cpu_prof_t1 - cpu_prof_t0)

                ! Mixed reaction rates: r_tilde is the NEW-STEP (theta) term only; the previous-step
                ! (1-theta) term of self and neighbours is already folded into c_hat_old by p_mixer, so
                ! c_hat = c_hat_old + Delta_t*r_tilde counts it exactly once.
                cpu_prof_t0 = cpu_prof_t1
                ! All downstream options (1/2/3/4) index the full waters array directly via
                ! p_r_tilde_idx (no vector-subscript gather/deep-copy): pass the whole array
                ! plus the mixing-water index list mc%col_1(2:mc%dim-2) (externals carry a zero ratio).
                call p_r_tilde_idx( &
                    this%chemistry%waters, mc%col_1(2:mc%dim - 2), &
                    rz%ind_var_act_species, &
                    mr_R%col_1(2:), &
                    mc%col_1(mc%dim - 1), mc%col_1(mc%dim), &
                    this%transport%time_discr%theta_r, r_tilde, &
                    curr_water_idx=mc%col_1(1), &
                    lambda_self=mr_R%col_1(1), &
                    wat_indices_zones=this%chemistry%wat_indices_zones)

                ! Positivity-preserving reactive update (same uniform-omega limiter as the mixers):
                ! keeps c_hat >= 0 WITHOUT modifying the reaction mixing ratios (dispersion unchanged).
                call add_limited_increment(c_hat_old, Delta_t*r_tilde, c_hat)
                call aq%set_conc_hat(c_hat(1:n_v_aq)) ! aqueous blocks of the Newton speciation input, for the conc_hat attribute
                !c_mix = c_hat  ! Update c_mix to reflect the new concentrations after the reactive update
                conc_nc = aq%get_conc_var_act()

                ! Compute reaction mixing ratio
                if (this%chemistry%r_down_opt == 4) then
                    lambda_R = this%transport%compute_mix_ratio_R_opt4(idx)
                else
                    lambda_R = mr_R%col_1(1)
                end if
                call cpu_time(cpu_prof_t1)
                cpu_data_prep_total = cpu_data_prep_total + (cpu_prof_t1 - cpu_prof_t0)

                ! Reactive mixing iteration
                cpu_prof_t0 = cpu_prof_t1
                call p_solver(aq, aq%get_conc_old_old_prim_new_rz(), c_hat, &
                    lambda_R, Delta_t, &
                    this%transport%time_discr%theta_r, conc_nc, conc_comp)
                call cpu_time(cpu_prof_t1)
                cpu_solver_total = cpu_solver_total + (cpu_prof_t1 - cpu_prof_t0)

                ! Update derived quantities
                cpu_prof_t0 = cpu_prof_t1
                call aq%compute_pH()
                call aq%compute_saturation_mins() !< saturation indices of all minerals at the end-of-step state
                call aq%compute_salinity()
                call aq%compute_ionic_strength()
                call cpu_time(cpu_prof_t1)
                cpu_derived_total = cpu_derived_total + (cpu_prof_t1 - cpu_prof_t0)

                ! Accumulate reaction amounts
                aq%cum_rxn_amnt_kin = aq%cum_rxn_amnt_kin + aq%rxn_amnt_kin
                sc%cum_rxn_amnt_kin = sc%cum_rxn_amnt_kin + sc%rxn_amnt_kin

                ! Compute equilibrium reaction rates
                cpu_prof_t0 = cpu_prof_t1
                if (sa%num_eq_reactions > 0) then
                    call aq%compute_Re(c_hat(1:n_v), Delta_t, lambda_R, this%transport%time_discr%theta_r)
                end if
                call cpu_time(cpu_prof_t1)
                cpu_eq_rates_total = cpu_eq_rates_total + (cpu_prof_t1 - cpu_prof_t0)

                ! Domain-integrated mass-balance contribution of this cell
                call aq%get_react_mix_mass_bal_terms(c_mix(1:n_v), c_hat(1:n_v), lambda_R, &
                    Delta_t, this%transport%time_discr%theta_r, mb_lhs_c, mb_transp_c, mb_react_c)
                if (size(mb_lhs_c) == n_v_ref) then
                    mb_lhs = mb_lhs + mb_lhs_c
                    mb_transp = mb_transp + mb_transp_c
                    mb_react = mb_react + mb_react_c
                end if
                deallocate(mb_lhs_c, mb_transp_c, mb_react_c)

                ! Update mineral / gas state
                cpu_prof_t0 = cpu_prof_t1
                if (associated(sc%mineral_zone)) then
                    call sc%compute_mass_bal_mins(Delta_t)
                    call sc%compute_conc_minerals_iter(Delta_t,lambda_R)
                end if

                ! Update gas state
                if (associated(aq%gas_chemistry)) then
                    call aq%gas_chemistry%compute_conc_gases_iter( &
                        Delta_t, aq%volume, [aq%rxn_rate_eq, aq%rxn_rate_kin_avg])
                    call aq%gas_chemistry%compute_vol_gas_species_conc()
                    call aq%gas_chemistry%compute_log_act_coeffs_gases()
                end if
                call cpu_time(cpu_prof_t1)
                cpu_solid_gas_total = cpu_solid_gas_total + (cpu_prof_t1 - cpu_prof_t0)

                !> per-cell arrays only: c_mix and c_hat_old are allocated by p_mixer, conc_nc by the
                !> get_conc_var_act assignment. c_hat, conc_comp and r_tilde are the size-cached work arrays
                !> (reallocated above only when n_v/n_p change) and must survive to the next cell.
                deallocate(c_mix, c_hat_old, conc_nc)
                end associate
            end do ! target waters
            ! Domain-integrated mass balance for this time step (per variable-activity species)
            call this%chemistry%write_react_mix_mass_bal_dom(idx_ref,iunit_mb,time,n_v_ref,mb_lhs,mb_transp,mb_react)
            !> Sum equilibrium reaction rates at interfaces if mix_by_zones_flag is set
            if (this%chemistry%mix_by_zones_flag) then
                do i=1,num_non_can_vec
                    call this%target_waters(ind_non_can_vec(i))%aq_chem%sum_eq_rxn_rates(this%chemistry%zone_waters(&
                        ind_non_can_vec(i),:))
                end do
            end if
            if (k == this%chemistry%chem_out_options%time_steps(kk)) then
                call cpu_time(cpu_prof_t0)
                call write_out_step_vars(this%chemistry, out_unit, time)
                if (kk < this%chemistry%chem_out_options%num_time_steps) kk=kk+1
                call cpu_time(cpu_prof_t1)
                cpu_output_total = cpu_output_total + (cpu_prof_t1 - cpu_prof_t0)
            end if
        end do ! time loop

        if (allocated(c_hat)) deallocate(c_hat, conc_comp, r_tilde)
        if (allocated(c_hat_old)) deallocate(c_hat_old)
        close(out_unit)
        close(iunit_mb)
        deallocate(mb_lhs, mb_transp, mb_react)

        ! --- Write the CPU profiling report (tagged filename) ---
        call cpu_time(cpu_total_t1)
        cpu_total_elapsed = cpu_total_t1 - cpu_total_t0
        call write_cpu_profile_report( &
            iunit_cpu, &
            num_profiled_steps, this%transport%time_discr%Num_time, num_non_can_vec, &
            cpu_total_elapsed, &
            cpu_update_old_total, cpu_displace_total, &
            cpu_data_prep_total, cpu_mix_total, cpu_getters_total, &
            cpu_solver_total, cpu_derived_total, cpu_eq_rates_total, &
            cpu_solid_gas_total, cpu_output_total, cpu_snapshot_total, &
            dir, cpu_root, 'solve_RT_ideal_cons_Lagr_stat_flux_1D')
    end select

    !print *, "Number of lumpings performed: ", this%chemistry%num_lump
end subroutine
