!> \file solve_RT_2D_ideal_cons_Lagr_stat_flux.f90
!> \brief Solve 2D reactive transport with ideal speciation, consistent WMA, and Lagrangian particle tracking
!> \details
!>   2D analogue of solve_RT_1D_ideal_cons_Lagr_stat_flux. Solves the 2D reactive
!>   transport problem using:
!>   - Ideal speciation (gamma = 1)
!>   - Consistent (non-lumped) Water Mixing Approach (WMA)
!>   - Lagrangian particle tracking with stationary flux
!>
!>   Operator splitting decouples transport (advection-dispersion via particle mixing)
!>   from reactions (equilibrium/kinetic chemistry at each target).
!>
!>   **Column-major processing (flow along +x).**
!>   The consistent WMA estimates the implicit upstream reaction rates from the
!>   *new* (this-step) reaction rates of the upstream cells, so every cell must be
!>   solved *after* all of its x-upstream cells. The cells are stored row-major
!>   (index = (row-1)*Nx + col), but a cell's upstream domain neighbours (any row,
!>   strictly smaller column) only all precede it under a **column-major** sweep
!>   (top->bottom within a column, columns left->right). This subroutine therefore
!>   visits the non-canonical target waters in column-major order. The upstream /
!>   downstream split itself is encoded in mix_react_indices, which is built by
!>   compute_mixing_ratios_Delta_t_homog_2D using the same physical-column rule
!>   (upstream <=> strictly smaller column; same-x neighbours are downstream).
!>
!> \param[inout] this Reactive transport object (2D transient)
!> \param[in] dir  Directory for output files
!> \param[in] root Root name for output file
!> \author Jordi Petchame-Guerrero
!> \date June 2026
subroutine solve_RT_ideal_cons_Lagr_stat_flux_2D(this, dir, root)
    use RT_m, only: RT_2D_transient_c, move_particles_stat_flux_2D
    use chemistry_m, only: write_out_step_vars, reset_wat_indices_zones !< Shared per-output-step writer (conc/react_rate/vol_fract); post-displacement zone regrouping
    use aqueous_chemistry_m, only: &
        aqueous_chemistry_c, &
        reactive_mixing_iter_EI_eq_anal_ideal, &
        reactive_mixing_iter_EI_kin_anal_ideal_opt2, &
        reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2, &
        compute_r_tilde_impl_opt1, compute_r_tilde_impl_opt2, &
        compute_r_tilde_impl_opt3, compute_r_tilde_impl_opt4, &
        compute_r_tilde_impl_opt1_idx, compute_r_tilde_impl_opt2_idx, &
        compute_r_tilde_impl_opt3_idx, compute_r_tilde_impl_opt4_idx, &
        compute_c_mix_global_cons, compute_c_mix_by_zones_cons, compute_c_hat_patankar
    use spatial_discr_2D_m, only: mesh_2D_Euler_homog_c
    use arrays_m, only: int_array_c
    implicit none

    ! Arguments
    class(RT_2D_transient_c) :: this
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: root

    ! Loop counters
    integer(kind=4) :: i, j, k, kk, mloop, io, idx, idx_o, num_lump

    ! Problem dimensions
    integer(kind=4) :: n_p, n_v, n_eq
    integer(kind=4) :: num_can_vec, num_non_can_vec
    integer(kind=4) :: num_mix_loc, out_unit
    integer(kind=4) :: Nx, Ny, c_col, r_row, m
    integer(kind=4) :: num_out_ts, num_out_waters, num_out_aq_sp
    integer(kind=4) :: num_waters_loc, w, max_n_conc, max_n_idx

    ! Index arrays
    integer(kind=4), allocatable :: ind_can_vec(:), ind_non_can_vec(:)
    integer(kind=4), allocatable :: mix_ext_indices(:) !< mixer externals: boundary waters (minus the Neumann outflow ones) + the inflow (first-column) cells, whose parcels are the boundary waters each step
    integer(kind=4), allocatable :: exit_wat_indices(:) !< Neumann outflow (right) boundary waters: hold the parcels that left the last column, grouped in those cells' zone columns (empty for a Dirichlet outflow)
    integer(kind=4) :: q !< loop index over the external waters
    integer(kind=4), allocatable :: colmaj_order(:) !< processing order: positions into ind_non_can_vec, column-major
    integer(kind=4), allocatable :: pos_of_idx(:)   !< reverse map: global target index -> position in ind_non_can_vec (0 if canonical)

    ! Time and mixing
    real(kind=8) :: time, Delta_t, lambda_R

    ! Concentration arrays
    real(kind=8), allocatable :: c_hat(:), c_mix(:), conc_comp(:)
    real(kind=8), allocatable :: c_tilde(:)
    real(kind=8), allocatable :: conc_nc(:), r_tilde(:)
    ! Domain-integrated reactive-mixing mass balance (per variable-activity species)
    integer(kind=4) :: iunit_mb                   !< mass-balance output file unit
    integer(kind=4) :: idx_ref                    !< reference target-water index (species names/dimension)
    integer(kind=4) :: n_v_ref                    !< reference number of variable-activity species
    real(kind=8), allocatable :: mb_lhs(:), mb_transp(:), mb_react(:) !< domain-summed terms
    real(kind=8), allocatable :: mb_lhs_c(:), mb_transp_c(:), mb_react_c(:) !< per-cell terms
    ! Global concentration caches for cross-zone conservative mixing (compute_c_mix_global)
    real(kind=8), allocatable :: all_conc_old(:,:)
    integer(kind=4), allocatable :: all_ind_aq_sp(:,:)
    ! Procedure pointers
    procedure(reactive_mixing_iter_EI_kin_anal_ideal_opt2), pointer :: p_solver => null()
    procedure(compute_r_tilde_impl_opt4_idx), pointer :: p_r_tilde_idx => null()
    procedure(compute_c_mix_by_zones_cons), pointer :: p_mixer

    ! By-reaction-zone previous-step map (rebuilt each time step): cols(jz)%col_1 = global target-water indices in old reaction zone jz
    type(int_array_c) :: wat_indices_zones
    integer(kind=4) :: num_rz

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

    ! Extract 2D grid dimensions (needed for the column-major processing order)
    if (this%chemistry%mix_by_zones_flag) then
        p_mixer => compute_c_mix_by_zones_cons
    else
        p_mixer => compute_c_mix_global_cons
    end if
    select type (mesh => this%transport%spatial_discr)
    type is (mesh_2D_Euler_homog_c)
        Nx = mesh%Num_cells_x
        Ny = mesh%Num_cells_y
    class default
        error stop "2D homogeneous mesh expected in solve_RT_ideal_cons_Lagr_stat_flux_2D"
    end select

    ! Count equilibrium reactions across all target waters
    n_eq = 0
    do i = 1, this%chemistry%num_target_waters
        n_eq = n_eq + this%target_waters(i)%aq_chem%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    end do

    time = 0d0
    kk = 2
    num_out_ts = this%chemistry%chem_out_options%num_time_steps
    num_out_waters = this%chemistry%chem_out_options%num_waters
    num_out_aq_sp = this%chemistry%chem_out_options%num_aq_species
    this%chemistry%chem_out_options%time_steps(num_out_ts) = this%transport%time_discr%Num_time

    out_unit = 7
    block
        character(len=:), allocatable :: rtag
        if (this%chemistry%lump_flag) then
            rtag = '_Lagr_lumped'
        else
            rtag = '_Lagr_consistent'
        end if
        open(out_unit, file=dir//trim(root)//rtag//'_output.out', form="formatted", access="sequential", status="unknown")
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
        ! Equilibrium-only reactions: the speciation does not depend on the reactions integration method (Euler/RKF45)
        p_solver => reactive_mixing_iter_EI_eq_anal_ideal
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

    ! Particle displacement method (only the standard cell-centred 2D mover is supported)
    if (this%transport%spatial_discr%targets_flag /= 0) then
        error stop "Stationary flux Lagrangian scheme not implemented yet"
    end if

    ! Extract canonical vectors from mixing ratios matrix
    call this%transport%mixing_ratios_conc%get_can_vec( &
        this%chemistry%CV_params%abs_tol, num_can_vec, ind_can_vec, num_non_can_vec, ind_non_can_vec)

    ! Build the column-major processing order of the non-canonical target waters.
    ! pos_of_idx(idx) = position in ind_non_can_vec (0 if idx is canonical/skipped).
    allocate(pos_of_idx(this%transport%spatial_discr%Num_targets))
    pos_of_idx = 0
    do i = 1, num_non_can_vec
        pos_of_idx(ind_non_can_vec(i)) = i
    end do
    allocate(colmaj_order(num_non_can_vec))
    m = 0
    do c_col = 1, Nx
        do r_row = 1, Ny
            idx = (r_row - 1) * Nx + c_col   !< row-major storage index of physical cell (r_row, c_col)
            if (pos_of_idx(idx) > 0) then
                m = m + 1
                colmaj_order(m) = pos_of_idx(idx)
            end if
        end do
    end do

    ! Initialise concentration history
    do i = 1, this%chemistry%num_waters
        call this%chemistry%waters(i)%set_conc_old()
        call this%chemistry%waters(i)%set_conc_old_old()
        call this%chemistry%waters(i)%solid_chemistry%set_conc_old()
        call this%chemistry%waters(i)%set_solid_chemistry_old()
    end do

    ! Output headers/blocks are written per step by write_out_step_vars (honours chem_out_options%vars)

    ! Pre-compute global caches for cross-zone conservative mixing (see compute_c_mix_global).
    num_waters_loc = this%chemistry%num_waters
    max_n_conc = 0
    do w = 1, num_waters_loc
        j = size(this%chemistry%waters(w)%conc_old)
        if (j > max_n_conc) max_n_conc = j
    end do
    allocate(all_conc_old(max_n_conc, num_waters_loc))
    max_n_idx = 0
    do w = 1, num_waters_loc
        j = this%chemistry%waters(w)%aq_phase%num_species
        if (j > max_n_idx) max_n_idx = j
    end do
    allocate(all_ind_aq_sp(max_n_idx, num_waters_loc))
    do w = 1, num_waters_loc
        all_ind_aq_sp(1:size(this%chemistry%waters(w)%indices_aq_species), w) = &
            this%chemistry%waters(w)%indices_aq_species
    end do

    ! Domain-integrated mass balance: reference cell (homogeneous reactive zone) + output file
    idx_ref = this%chemistry%tar_wat_indices(ind_non_can_vec(1))
    n_v_ref = this%chemistry%waters(idx_ref)%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    allocate(mb_lhs(n_v_ref), mb_transp(n_v_ref), mb_react(n_v_ref))
    !> Tagged root shared by mass-balance and CPU-profile filenames.
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
    num_profiled_steps = this%transport%time_discr%Num_time   !< every step is profiled (no sampling/scaling)
    ! CPU profiling filename carries the transport (Lagr) and WMA (lumped/consistent) tags,
    ! matching the {root}{_Lagr|_Euler}{_lumped|_consistent}_*.out scheme.
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
    open(newunit=iunit_cpu, file=dir//cpu_root//'_cpu_profile.out', status='replace', action='write')

    !> Write initial state (t=0) before the first time step, honouring chem_out_options%vars
    if (this%chemistry%chem_out_options%time_steps(1)==0) then
        call write_out_step_vars(this%chemistry, out_unit, 0d0)
    end if
    !> Mixer grouping of the waters that are not plain target waters. The inflow (first-column) cells'
    !> parcels are the boundary waters (zone not in the domain): external, mixed directly. With a
    !> Neumann outflow, the right boundary waters are the parcels that left the last column: their
    !> previous-step concentrations were equilibrated in those cells' zones, so they join those zones'
    !> within-zone mix (classified by reset_wat_indices_zones) and are NOT externals.
    if (this%transport%BCs%labels(2)==2) then
        exit_wat_indices=this%chemistry%bd_waters_indices(this%chemistry%num_left_bd_waters+1:&
            this%chemistry%num_left_bd_waters+this%chemistry%num_right_bd_waters)
    else
        allocate(exit_wat_indices(0))
    end if
    mix_ext_indices=[pack(this%chemistry%ext_waters_indices, &
        [(all(this%chemistry%ext_waters_indices(q)/=exit_wat_indices), q=1,this%chemistry%num_ext_waters)]), &
        (this%chemistry%tar_wat_indices(1+(r_row-1)*Nx), r_row=1,this%chemistry%num_left_bd_waters)]
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
        do i = 1, this%chemistry%num_target_waters
            call this%target_waters(i)%update_attributes()
        end do
        call cpu_time(cpu_prof_t1)
        cpu_update_old_total = cpu_update_old_total + (cpu_prof_t1 - cpu_prof_t0)

        ! Transport: displace particles (advection along +x) and inject inflow
        call cpu_time(cpu_prof_t0)
        call move_particles_stat_flux_2D(this, k)
        call this%introduce_particles(k)
        !> Regroup the target waters by the reactive zone of their previous-step state (the zone of
        !> the cell each parcel comes from, solid_chemistry_old): wat_indices_zones is the
        !> classification compute_c_mix_by_zones uses for the within-zone mixing. The Neumann
        !> outflow waters (exited parcels) are classified with them.
        call reset_wat_indices_zones(this%chemistry)

        ! Refresh global conc_old cache AFTER displacement (post-displacement conc_old)
        do w = 1, num_waters_loc
            all_conc_old(1:size(this%chemistry%waters(w)%conc_old), w) = &
                this%chemistry%waters(w)%conc_old
        end do
        call cpu_time(cpu_prof_t1)
        cpu_displace_total = cpu_displace_total + (cpu_prof_t1 - cpu_prof_t0)

        ! Build the target-water-indices-per-old-reaction-zone map used by the by-zone
        ! previous-time-step contribution to r_tilde (rebuilt each step from solid_chemistry_old).
        call this%chemistry%link_target_waters_old_reaction_zones(num_rz, wat_indices_zones)

        !----------------------------------------------------------------
        ! Target waters loop (non-canonical vectors only), COLUMN-MAJOR order
        ! so that every cell's x-upstream cells are solved before it.
        !----------------------------------------------------------------
        do mloop = 1, num_non_can_vec
            i = colmaj_order(mloop)
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
            ! (within-zone mixing + speciation, then cross-zone aqueous mixing). This folds the
            ! (1-theta) previous-step reaction term into c_mix, so p_r_tilde_idx is new-step only.
            ! The reaction mixing ratios share the concentration stencil (0 for the external waters).
            call cpu_time(cpu_prof_t0)
            call p_mixer(this%chemistry%waters, this%chemistry%zone_waters(idx,:), &
                mc%col_1(1:mc%dim-2), mr_c%col_1(1:mc%dim-2), mr_R%col_1(1:mc%dim-2), &
                this%chemistry%wat_indices_zones, mix_ext_indices, Delta_t, &
                this%transport%time_discr%theta_r, c_mix, c_tilde, .true.)
            call cpu_time(cpu_prof_t1)
            cpu_mix_total = cpu_mix_total + (cpu_prof_t1 - cpu_prof_t0)

            ! Mixed reaction rates: implicit upstream, estimated downstream. curr_water lets opt4
            ! take solid/gas rates from the (fixed) current cell in the theta-average.
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
                wat_indices_zones=wat_indices_zones)

            ! Positivity-preserving reactive update (conservative reaction-extent limiter):
            ! keeps c_hat >= 0 WITHOUT modifying the reaction mixing ratios (dispersion unchanged).
            call compute_c_hat_patankar(c_mix, Delta_t, r_tilde, c_hat)
            call aq%set_conc_hat(c_hat) ! aqueous blocks of the Newton speciation input, for the conc_hat attribute
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

            deallocate(c_mix, c_tilde, conc_nc)
            end associate
        end do ! target waters

        ! Domain-integrated mass balance for this time step (per variable-activity species)
        call this%chemistry%write_react_mix_mass_bal_dom(idx_ref,iunit_mb,time,n_v_ref,mb_lhs,mb_transp,mb_react)

        !----------------------------------------------------------------
        ! Output selected waters at output time steps. Decoupled from the
        ! (column-major) solve order so the listed output waters are written
        ! in their natural order regardless of processing order.
        !----------------------------------------------------------------
        call cpu_time(cpu_prof_t0)
        if (k == this%chemistry%chem_out_options%time_steps(kk)) then
            call write_out_step_vars(this%chemistry, out_unit, time)
            if (kk < num_out_ts) kk = kk + 1
        end if
        call cpu_time(cpu_prof_t1)
        cpu_output_total = cpu_output_total + (cpu_prof_t1 - cpu_prof_t0)
    end do ! time loop

    if (allocated(c_hat)) deallocate(c_hat, conc_comp, r_tilde)
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
        dir, cpu_root, 'solve_RT_ideal_cons_Lagr_stat_flux_2D')
end subroutine
