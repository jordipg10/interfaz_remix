!> @file solve_reactive_mixing_lump.f90
!> @brief Solves the reactive mixing problem assuming NON-ideal conditions with a lumped formulation.
!> @details Non-ideal counterpart of solve_reactive_mixing_ideal_lump. It shares the same structure
!> (canonical/non-canonical vector loop, lumped mixing ratios, global conservative mixing via
!> compute_c_mix_global, per-target-water dimension caches and CPU profiling) with the single
!> time-step iteration inlined in the main time loop. The reactive step uses the non-ideal
!> iterations and non-ideal activity coefficients: mixing_iter_comp for equilibrium-only systems
!> and, for equilibrium + kinetic systems, the globalised Euler-implicit Newton
!> reactive_mixing_iter_EI_eq_kin_anal_opt2 (or the adaptive RKF45 variant when the reaction
!> integration method is 3), as the non-ideal Lagrangian lumped solver does. The lumped mixer folds
!> the previous-step (1-theta_r) kinetic term of the target water into the speciation input, so the
!> kinetic iteration only carries the implicit self term. Kinetic-only systems are not implemented.
!>
!> @author Jordi Petchamé-Guerrero
!> @date November 2025

subroutine solve_reactive_mixing_lump(this,dir,root,mixing_ratios_conc,mixing_ratios_R,mix_conc_indices,&
    time_discr,theta_r)
    use chemistry_m, only: chemistry_c, write_out_step_vars
    !> Non-ideal reactive-mixing iterations and global conservative mixing
    use aqueous_chemistry_m, only: aqueous_chemistry_c, write_celldbg_cell, mixing_iter_comp, compute_c_mix_global_lump, &
        compute_c_mix_by_zones_lump, reactive_mixing_iter_EI_eq_kin_anal_opt2, reactive_mixing_iter_RKF45_eq_kin_anal_opt2
    use time_discr_m, only: time_discr_c
    use arrays_m, only: real_array_c, int_array_c
    implicit none
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: root
    class(real_array_c), intent(in) :: mixing_ratios_conc !> mixing ratios concentration matrix
    class(real_array_c), intent(in) :: mixing_ratios_R !> mixing ratios reaction rates matrix (mix_conc_indices layout)
    class(int_array_c), intent(in) :: mix_conc_indices !> indices of target waters that mix with each target water
    class(time_discr_c), intent(in) :: time_discr !> time discretisation object
    real(kind=8), intent(in) :: theta_r !> time weighting factor for chemical reactions (0=explicit, 1=implicit)
!> Loop counters and indices
    integer(kind=4) :: i, k, kk
    integer(kind=4) :: idx !< Absolute target-water index tar_wat_indices(ind_non_can_vec(i))
!> Problem dimensions
    integer(kind=4) :: n_p, n_v, n_eq
    integer(kind=4) :: num_can_vec, num_non_can_vec
    integer(kind=4) :: unit !< Output file unit
    integer(kind=4) :: iunit_cpu !< CPU-profiling file unit
    character(len=:), allocatable :: cpu_root !< root + Euler/lumped tag for the CPU-profiling filename
    integer(kind=4) :: num_time_loc !< Cached total number of time steps
!> Index arrays
    integer(kind=4), allocatable :: ind_can_vec(:)
    integer(kind=4), allocatable :: ind_non_can_vec(:)
!> Time variables
    real(kind=8) :: time, Delta_t
!> Concentration and reaction-rate arrays
    real(kind=8), allocatable :: conc_nc(:)
    real(kind=8), allocatable :: conc_comp(:)
    real(kind=8), allocatable :: lumped_lambdas(:)
!> Per-step working arrays
    real(kind=8), allocatable :: c_mix(:)
    real(kind=8), allocatable :: c_tilde(:)
    real(kind=8), allocatable :: mal_residual(:)
    real(kind=8) :: mal_max
!> Domain-integrated reactive-mixing mass balance (per variable-activity species)
    integer(kind=4) :: iunit_mb !< mass-balance output file unit
    integer(kind=4) :: iunit_celldbg !< per-cell [CELL-DBG] trace file unit
    integer(kind=4) :: idx_ref !< reference target-water index (species names/dimension)
    integer(kind=4) :: n_v_ref !< reference number of variable-activity species
    real(kind=8), allocatable :: mb_lhs(:), mb_transp(:), mb_react(:) !< domain-summed terms
    real(kind=8), allocatable :: mb_lhs_c(:), mb_transp_c(:), mb_react_c(:) !< per-cell terms
!> Pre-computed per-target-water caches
    logical, allocatable :: has_minerals_flag(:)
    logical, allocatable :: has_gases_flag(:)
    integer(kind=4), allocatable :: n_p_cache(:)
    integer(kind=4), allocatable :: n_v_cache(:)
    integer(kind=4), allocatable :: n_eq_cache(:)
    integer(kind=4) :: max_n_p, max_n_v
!> CPU profiling variables
    real(kind=8) :: cpu_total_t0, cpu_total_t1, cpu_total_elapsed
    real(kind=8) :: cpu_prof_t0, cpu_prof_t1
    real(kind=8) :: cpu_solver_total
    real(kind=8) :: cpu_output_total
!> Procedure pointers
    logical :: kin_flag !> equilibrium + kinetic system: the reactive step is p_solver_kin, otherwise mixing_iter_comp
    procedure(reactive_mixing_iter_EI_eq_kin_anal_opt2), pointer :: p_solver_kin=>null() !> non-ideal equilibrium + kinetic iteration (EI or RKF45)
    procedure(compute_c_mix_by_zones_lump), pointer :: p_mixer

!> -----------------------------------------------------------------------
!> Initialization
    call cpu_time(cpu_total_t0)
    if (this%mix_by_zones_flag) then
        p_mixer => compute_c_mix_by_zones_lump
    else
        p_mixer => compute_c_mix_global_lump
    end if
    n_eq=0
    do i=1,this%num_target_waters
        n_eq=n_eq+this%waters(this%tar_wat_indices(i))%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    end do

    time=0d0
    kk=2
    num_time_loc = time_discr%Num_time
    do i=2,this%chem_out_options%num_time_steps-1
        if (this%chem_out_options%time_steps(i) > num_time_loc) then
            print *, "Output time step index", this%chem_out_options%time_steps(i), &
                "exceeds total number of time steps", num_time_loc
            error stop "Mismatch between out_opts.dat and discr_temp.dat: output time step exceeds total time steps"
        end if
    end do
    unit=7
    block
        character(len=:), allocatable :: rtag
        if (this%lump_flag) then
            rtag = '_Euler_lumped'
        else
            rtag = '_Euler_consistent'
        end if
        open(unit,file=dir//trim(root)//rtag//'_output.out',form="formatted",access="sequential",status="unknown")
    end block
    !> Output headers/blocks are written per step by write_out_step_vars (honours chem_out_options%vars)

!> -----------------------------------------------------------------------
!> Solver selection (as solve_RT_1D_lump_Lagr_stat_flux, the non-ideal Lagrangian lumped solver)
    kin_flag=.false.
    if (this%chem_syst%num_kin_reacts>0 .and. n_eq>0) then
        if (this%Jac_opt/=1) error stop "Non-ideal lumped solver: only analytical Jacobians (Jac_opt = 1) are implemented"
        if (this%rk_avg_opt/=2) error stop "Non-ideal lumped solver: rk average option not implemented yet (only option 2)"
        kin_flag=.true.
        if (time_discr%int_method_chem==3) then
            p_solver_kin=>reactive_mixing_iter_RKF45_eq_kin_anal_opt2 !> adaptive Runge-Kutta-Fehlberg 4(5) kinetics
        else
            p_solver_kin=>reactive_mixing_iter_EI_eq_kin_anal_opt2 !> Euler-implicit (theta-method) kinetics, globalised Newton
        end if
    else if (n_eq>0) then
        continue !> equilibrium-only: mixing_iter_comp (called directly; its interface has no c1_downstream argument)
    else
        error stop "Non-ideal lumped solver: kinetic-only reactions not implemented yet"
    end if

!> -----------------------------------------------------------------------
!> Pre-compute canonical/non-canonical vectors and lumped lambdas
    call mixing_ratios_conc%get_can_vec(this%CV_params%abs_tol,num_can_vec,ind_can_vec,num_non_can_vec,ind_non_can_vec)
    allocate(lumped_lambdas(mixing_ratios_conc%num_cols))
    do i=1,mixing_ratios_conc%num_cols
        lumped_lambdas(i)=sum(mixing_ratios_conc%cols(i)%col_1) !< lumped lambda = sum of the CONCENTRATION mixing ratios (same definition in every lumped solver)
    end do

!> -----------------------------------------------------------------------
!> Initial conditions setup
    do i=1,this%num_waters
        call this%waters(i)%set_conc_old()
        call this%waters(i)%set_conc_old_old()
        call this%waters(i)%solid_chemistry%set_conc_old()
        call this%waters(i)%solid_chemistry%set_conc_old_old()
    end do

!> -----------------------------------------------------------------------
!> Pre-compute per-target-water dimensions and flags
    allocate(n_p_cache(num_non_can_vec))
    allocate(n_v_cache(num_non_can_vec))
    allocate(n_eq_cache(num_non_can_vec))
    allocate(has_minerals_flag(num_non_can_vec))
    allocate(has_gases_flag(num_non_can_vec))
    max_n_p = 0
    max_n_v = 0
    do i=1,num_non_can_vec
        idx = this%tar_wat_indices(ind_non_can_vec(i))
        associate(spec_alg => this%waters(idx)%solid_chemistry%reactive_zone%speciation_alg)
            n_p_cache(i) = spec_alg%num_prim_species
            n_v_cache(i) = spec_alg%num_var_act_species
            n_eq_cache(i) = spec_alg%num_eq_reactions
        end associate
        has_minerals_flag(i) = associated(this%waters(idx)%solid_chemistry%mineral_zone)
        has_gases_flag(i) = associated(this%waters(idx)%gas_chemistry)
        if (n_p_cache(i) > max_n_p) max_n_p = n_p_cache(i)
        if (n_v_cache(i) > max_n_v) max_n_v = n_v_cache(i)
    end do

!> Pre-allocate reusable buffers
    allocate(conc_comp(max_n_p))

!> Domain-integrated mass balance: reference cell (homogeneous reactive zone) + output file
    idx_ref = this%tar_wat_indices(ind_non_can_vec(1))
    n_v_ref = n_v_cache(1)
    allocate(mb_lhs(n_v_ref), mb_transp(n_v_ref), mb_react(n_v_ref))
    !> Tagged root shared by mass-balance and CPU-profile filenames (Eulerian reactive-mixing path).
    if (this%lump_flag) then
        cpu_root = root//'_Euler_lumped'
    else
        cpu_root = root//'_Euler_consistent'
    end if
    open(newunit=iunit_mb, file=dir//cpu_root//'_mass_balance.out', status='replace', action='write')
    write(iunit_mb,*) 'Domain-integrated reactive-mixing mass balance (per variable-activity species)'
    write(iunit_mb,*) 'Per time step & species: LHS = sum V*(c^{k+1}-c^k)/Dt, RHS = transport + reaction, residual = LHS - RHS'
    !> Per-cell [CELL-DBG] trace file (bulk calcite cells 3/8 + interface cells 15/16/17)
    open(newunit=iunit_celldbg, file=dir//cpu_root//'_cell_dbg.out', status='replace', action='write')

!> -----------------------------------------------------------------------
!> Initialize profiling accumulators
    cpu_solver_total = 0d0
    cpu_output_total = 0d0

!> Open CPU profiling output file (Eulerian reactive-mixing path; same tagged root as the mass balance)
    open(newunit=iunit_cpu, file=dir//cpu_root//'_cpu_profile.out', status='replace', action='write')

!> -----------------------------------------------------------------------
!> Write initial state (t=0) before the first time step, honouring chem_out_options%vars
    if (this%chem_out_options%time_steps(1)==0) then
        call write_out_step_vars(this, unit, 0d0)
    end if

!> Main time loop
    do k=1,num_time_loc
        Delta_t=time_discr%get_Delta_t(k)
        time=time+Delta_t
        !> Reset the domain-integrated mass-balance accumulators for this time step
        mb_lhs = 0d0
        mb_transp = 0d0
        mb_react = 0d0
        !> ---- Single time-step reactive-mixing iteration (inlined) ----
        call cpu_time(cpu_prof_t0)
        !> Update old attributes (copies previous step's new conc into conc_old)
        do i = 1, this%num_target_waters
            call this%waters(this%tar_wat_indices(i))%update_attributes() !< Shift all old attributes
        end do
        !> Target waters loop (non-canonical columns only)
        do i = 1, num_non_can_vec
            idx = this%tar_wat_indices(ind_non_can_vec(i))
            n_p = n_p_cache(i)
            n_v = n_v_cache(i)
            n_eq = n_eq_cache(i)
            if (mix_conc_indices%cols(ind_non_can_vec(i))%col_1(1) /= idx) then
                error stop "Target waters not in the right order"
            end if
            !> Conservative mixing (per-reactive-zone, non-ideal speciation). zone_waters is indexed by
            !> TARGET-WATER position (num_target_solids x num_reactive_zones), i.e. ind_non_can_vec(i), not
            !> by the position i in the non-canonical list (they differ as soon as a target water is canonical)
            call p_mixer(this%waters, this%zone_waters(ind_non_can_vec(i), :), &
                mix_conc_indices%cols(ind_non_can_vec(i))%col_1(1:mix_conc_indices%cols(ind_non_can_vec(i))%dim-2), &
                mixing_ratios_conc%cols(ind_non_can_vec(i))%col_1(1:mix_conc_indices%cols(ind_non_can_vec(i))%dim-2), &
                lumped_lambdas(ind_non_can_vec(i)), &
                this%wat_indices_zones, this%ext_waters_indices, Delta_t, theta_r, c_mix, c_tilde, ideal=.false.)
            !> Non-ideal reactive mixing solver: equilibrium speciation with activity coefficients, plus the
            !> implicit kinetic self term for equilibrium + kinetic systems (the explicit kinetic terms are
            !> already in c_mix, the first mixer output). Both iterations finalise pH, salinity, ionic
            !> strength / activity coefficients and saturation indices of the water themselves.
            conc_nc = this%waters(idx)%get_conc_var_act()
            call this%waters(idx)%set_conc_hat(c_mix(1:n_v)) !> aqueous blocks of the Newton speciation input (the lumped mixer's output), for the conc_hat attribute
            if (kin_flag) then
                call p_solver_kin(this%waters(idx), this%waters(idx)%get_conc_old_old_prim_new_rz(), c_mix(1:n_v), &
                    lumped_lambdas(ind_non_can_vec(i)), Delta_t, theta_r, conc_nc, conc_comp(1:n_p))
            else
                call mixing_iter_comp(this%waters(idx), this%waters(idx)%get_conc_old_old_prim_new_rz(), c_mix(1:n_v), &
                    lumped_lambdas(ind_non_can_vec(i)), Delta_t, theta_r, conc_nc, conc_comp(1:n_p))
            end if
            !> Mass action law check for equilibrium reactions (non-ideal activity coefficients)
            if (n_eq > 0) then
                block
                    real(kind=8), allocatable :: log_gamma_mal(:)
                    log_gamma_mal = this%waters(idx)%get_log_gamma_var_act()
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
                    error stop "Mass action law violation after p_solver"
                end if
                deallocate(mal_residual)
            end if
            !> Equilibrium reaction amounts
            if (n_eq > 0) then
                call this%waters(idx)%compute_Re(c_mix(1:n_v), Delta_t, &
                    lumped_lambdas(ind_non_can_vec(i)), theta_r)
            end if
            !> [CELL-DBG] intermediate results for bulk calcite cells 3/8 and interface cells 15/16/17
            if (k<=10 .and. (this%num_target_waters<=10 .or. idx==3 .or. idx==8 .or. idx==15 .or. idx==16 .or. idx==17)) then
                call write_celldbg_cell(this%waters(idx), iunit_celldbg, k, idx, n_v, c_mix, c_tilde, c_mix)
            end if
            !> Solid/gas chemistry (of the TARGET water idx: ind_non_can_vec(i) is its position in the
            !> target list, not an index into this%waters)
            if (has_minerals_flag(i)) then
                call this%waters(idx)%solid_chemistry%compute_mass_bal_mins(Delta_t)
                call this%waters(idx)%solid_chemistry%compute_conc_minerals_iter(Delta_t,lumped_lambdas(ind_non_can_vec(i)))
            end if
            if (has_gases_flag(i)) then
                call this%waters(idx)%gas_chemistry%compute_conc_gases_iter(Delta_t, &
                    this%waters(idx)%volume, &
                    [this%waters(idx)%rxn_rate_eq, this%waters(idx)%rxn_rate_kin_avg])
                call this%waters(idx)%gas_chemistry%compute_vol_gas_species_conc()
                !call this%waters(idx)%gas_chemistry%compute_log_act_coeffs_gases()
            end if
            !> Accumulate reaction amounts: kinetic AND equilibrium (update_cum_rxn_amnt adds rxn_amnt_kin
            !> and the rxn_amnt_eq that compute_Re set this step); the abs_cum_react_amnt output reports
            !> cum_rxn_amnt_eq for the aqueous equilibrium reactions and the equilibrium minerals
            call this%waters(idx)%update_cum_rxn_amnt()
            call this%waters(idx)%solid_chemistry%update_cum_rxn_amnt()
            if (has_gases_flag(i)) call this%waters(idx)%gas_chemistry%update_cum_rxn_amnt()
            !> Domain-integrated mass-balance contribution of this cell (as solve_reactive_mixing_ideal_lump_iter):
            !> transport-only vector = c_tilde (second mixer output, the conservative mix, aqueous var-act layout),
            !> reacted vector = c_mix (first mixer output: conservative mix + the previous-step (1-theta_r) kinetic
            !> term, i.e. the speciation input, full var-act layout). The two layouts coincide when every variable
            !> activity species is aqueous; otherwise c_mix stands in for both (the split between the two vectors
            !> only moves the previous-step kinetic term between the transport and reaction columns, it does not
            !> affect the residual). The same lumped lambda as compute_Re is used, so the reaction term closes
            !> with the fitted Re.
            if (size(c_tilde)==n_v) then
                call this%waters(idx)%get_react_mix_mass_bal_terms(c_tilde(1:n_v), c_mix(1:n_v), &
                    lumped_lambdas(ind_non_can_vec(i)), Delta_t, theta_r, &
                    mb_lhs_c, mb_transp_c, mb_react_c)
            else
                call this%waters(idx)%get_react_mix_mass_bal_terms(c_mix(1:n_v), c_mix(1:n_v), &
                    lumped_lambdas(ind_non_can_vec(i)), Delta_t, theta_r, &
                    mb_lhs_c, mb_transp_c, mb_react_c)
            end if
            if (size(mb_lhs_c) == n_v_ref) then
                mb_lhs = mb_lhs + mb_lhs_c
                mb_transp = mb_transp + mb_transp_c
                mb_react = mb_react + mb_react_c
            end if
            deallocate(mb_lhs_c, mb_transp_c, mb_react_c)
            deallocate(c_mix, c_tilde, conc_nc)
        end do
        call cpu_time(cpu_prof_t1)
        cpu_solver_total = cpu_solver_total + (cpu_prof_t1 - cpu_prof_t0)
        !> Output writing honouring chem_out_options%vars (conc / react_rate / vol_fract)
        call cpu_time(cpu_prof_t0)
        if (k==this%chem_out_options%time_steps(kk)) then
            call write_out_step_vars(this, unit, time)
            if (kk<this%chem_out_options%num_time_steps) kk=kk+1
        end if
        call cpu_time(cpu_prof_t1)
        cpu_output_total = cpu_output_total + (cpu_prof_t1 - cpu_prof_t0)
        !> Domain-integrated mass balance for this time step (per variable-activity species)
        call this%write_react_mix_mass_bal_dom(idx_ref,iunit_mb,time,n_v_ref,mb_lhs,mb_transp,mb_react)
    end do

!> -----------------------------------------------------------------------
!> CPU profiling report
    call cpu_time(cpu_total_t1)
    cpu_total_elapsed = cpu_total_t1 - cpu_total_t0
    write(iunit_cpu,*) '====================================================================================='
    write(iunit_cpu,*) '  CPU PROFILING REPORT: solve_reactive_mixing_lump'
    write(iunit_cpu,*) '====================================================================================='
    write(iunit_cpu,'(A,I10)')     '  Total time steps:       ', num_time_loc
    write(iunit_cpu,'(A,I10)')     '  Non-canonical waters:   ', num_non_can_vec
    write(iunit_cpu,*) '-------------------------------------------------------------------------------------'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Iteration (total):      ', cpu_solver_total, ' s  (', cpu_solver_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Output writing:         ', cpu_output_total, ' s  (', cpu_output_total/60d0, ' min)'
    write(iunit_cpu,*) '-------------------------------------------------------------------------------------'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  TOTAL SUBROUTINE TIME:  ', cpu_total_elapsed, ' s  (', cpu_total_elapsed/60d0, ' min)'
    write(iunit_cpu,*) '====================================================================================='
    close(iunit_cpu)
    print *, 'CPU profiling written to: ', dir//cpu_root//'_cpu_profile.out'

!> -----------------------------------------------------------------------
!> Close output file and clean up
    close(unit)
    close(iunit_mb)
    close(iunit_celldbg)
    deallocate(mb_lhs, mb_transp, mb_react)
    deallocate(conc_comp)
    deallocate(lumped_lambdas)
    deallocate(n_p_cache, n_v_cache, n_eq_cache)
    deallocate(has_minerals_flag, has_gases_flag)
end subroutine solve_reactive_mixing_lump
