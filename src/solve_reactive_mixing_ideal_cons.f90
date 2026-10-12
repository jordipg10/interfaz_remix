!> @brief Solves reactive mixing problem assuming ideal conditions with conservative formulation
!> @details Computes concentrations, activities, activity coefficients, reaction rates
!> and volumetric fractions of minerals (if present) using operator splitting.
subroutine solve_reactive_mixing_ideal_cons(this,dir,root,mixing_ratios_conc,mixing_ratios_R,&
    mix_conc_indices,time_discr,theta_r)
    use chemistry_m, only: chemistry_c, write_out_step_vars, check_surf_compl_act_sum
    use aqueous_chemistry_m, only: aqueous_chemistry_c, &
        reactive_mixing_iter_EI_eq_anal_ideal, &
        compute_c_mix_by_zones_cons, compute_c_mix_global_cons, &
        reactive_mixing_iter_EI_kin_anal_ideal_opt2, &
        reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2, &
        compute_c_hat_patankar, add_limited_increment, &
        compute_r_tilde_impl_opt1_idx, compute_r_tilde_impl_opt2_idx, &
        compute_r_tilde_impl_opt3_idx, compute_r_tilde_impl_opt4_idx
    use arrays_m, only: real_array_c, int_array_c
    use time_discr_m, only: time_discr_c
    implicit none
!> Arguments
    class(chemistry_c) :: this
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: root
    class(real_array_c), intent(in) :: mixing_ratios_conc
    class(real_array_c), intent(in) :: mixing_ratios_R !> mix_conc_indices layout (0 for conc-only waters)
    class(int_array_c), intent(in) :: mix_conc_indices
    class(time_discr_c), intent(in) :: time_discr
    real(kind=8), intent(in) :: theta_r
!> Loop counters and indices
    integer(kind=4) :: i, j, k, kk, ii, m
    integer(kind=4) :: idx          !< Cached target water index: tar_wat_indices(ind_non_can_vec(i))
    integer(kind=4) :: j_mix        !< Loop index for mixing waters
!> Problem dimensions
    integer(kind=4) :: n_p, n_v, n_eq
    integer(kind=4) :: num_can_vec, num_non_can_vec
    integer(kind=4) :: unit
    integer(kind=4) :: iunit_cpu
    character(len=:), allocatable :: cpu_root !< root + transport/WMA tags for the CPU profiling filename
    integer(kind=4) :: num_time_loc
    integer(kind=4) :: kk_snap      !< Counter for chemistry snapshot output time steps
    integer(kind=4) :: num_lump
!> Index arrays
    integer(kind=4), allocatable :: ind_can_vec(:)
    integer(kind=4), allocatable :: ind_non_can_vec(:)
!> Time variables
    real(kind=8) :: time, Delta_t
!> Concentration and reaction rate arrays
    real(kind=8), allocatable :: c_hat(:) !< reacted mix: c_hat_old + Delta_t*r_tilde (limited)
    real(kind=8), allocatable :: r_tilde(:)
    real(kind=8), allocatable :: c_hat_old(:) !< p_mixer output: transported mix with the previous-step (1-theta) reaction term folded in
    real(kind=8), allocatable :: c_mix(:) !< p_mixer output: purely conservative mix (aqueous var-act species, no reaction)
    real(kind=8), allocatable :: conc_nc(:)
    real(kind=8), allocatable :: mal_residual(:)
    real(kind=8) :: mal_max
    real(kind=8), allocatable :: conc_comp(:)
    real(kind=8), allocatable :: lambdas_R(:)
!> Domain-integrated reactive-mixing mass balance (per variable-activity species)
    integer(kind=4) :: iunit_mb !< mass-balance output file unit
    integer(kind=4) :: idx_ref !< reference target-water index (species names/dimension)
    integer(kind=4) :: n_v_ref !< reference number of variable-activity species
    real(kind=8), allocatable :: mb_lhs(:), mb_transp(:), mb_react(:) !< domain-summed terms
    real(kind=8), allocatable :: mb_lhs_c(:), mb_transp_c(:), mb_react_c(:) !< per-cell terms
!> Pre-computed caches
    logical, allocatable :: has_minerals_flag(:)
    logical, allocatable :: has_gases_flag(:)
    integer(kind=4), allocatable :: n_p_cache(:)
    integer(kind=4), allocatable :: n_v_cache(:)
    integer(kind=4), allocatable :: n_eq_cache(:)
    integer(kind=4) :: max_n_p
!> CPU profiling variables
    real(kind=8) :: cpu_total_t0, cpu_total_t1, cpu_total_elapsed
    real(kind=8) :: cpu_prof_t0, cpu_prof_t1
    real(kind=8) :: cpu_update_old_total
    real(kind=8) :: cpu_mix_total
    real(kind=8) :: cpu_r_tilde_total
    real(kind=8) :: cpu_solver_total
    real(kind=8) :: cpu_derived_total
    real(kind=8) :: cpu_eq_rates_total
    real(kind=8) :: cpu_solid_gas_total
    real(kind=8) :: cpu_output_total
    real(kind=8) :: cpu_snapshot_total
    real(kind=8) :: cpu_data_prep_total
    real(kind=8) :: cpu_getters_total
    real(kind=8) :: cpu_cum_rxn_amnt_kin_total
    real(kind=8) :: cpu_c_hat_total
!> Procedure pointers
    procedure(reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2), pointer :: p_solver=>null()
    procedure(compute_r_tilde_impl_opt4_idx), pointer :: p_r_tilde=>null()
    procedure(compute_c_mix_by_zones_cons), pointer :: p_mixer

!> -----------------------------------------------------------------------
!> Initialization
    write(0,*) '[SOLVER] ideal_cons entered'; flush(0)
    call cpu_time(cpu_total_t0)
    if (this%mix_by_zones_flag) then
        p_mixer => compute_c_mix_by_zones_cons
    else
        p_mixer => compute_c_mix_global_cons
    end if
    n_eq=0
    do i=1,this%num_target_waters
        !> Only count eq reactions from waters that also have kinetic reactions (true domain waters).
        !> Boundary-adjacent waters with eq-only minerals should not influence solver selection.
        if (this%waters(this%tar_wat_indices(i))%solid_chemistry%mineral_zone%num_minerals_kin>0) then
            n_eq=n_eq+this%waters(this%tar_wat_indices(i))%solid_chemistry%mineral_zone%num_minerals_eq
        end if
    end do

    time=0d0
    kk=2
    ii=1
    kk_snap=2
    num_time_loc=time_discr%Num_time
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
!> Solver selection
!> The kin-only solver assumes n_p == n_v (no secondary species).
!> When aqueous speciation equilibria create secondary species (n_p < n_v),
!> we must use the eq_kin solver even if there are no mineral equilibria.
    n_p=this%waters(this%tar_wat_indices(1))%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_v=this%waters(this%tar_wat_indices(1))%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    if (this%chem_syst%num_kin_reacts>0 .and. (n_eq>0 .or. n_p<n_v)) then
        if (this%Jac_opt==1) then
            if (this%rk_avg_opt==2) then
                p_solver=>reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2
            else
                error stop "rk average option not implemented yet"
            end if
        end if
    else if (n_eq>0 .or. n_p<n_v) then
        p_solver=>reactive_mixing_iter_EI_eq_anal_ideal
    else
        if (this%Jac_opt==1) then
            if (this%rk_avg_opt==2) then
                p_solver=>reactive_mixing_iter_EI_kin_anal_ideal_opt2
            else
                error stop "rk average option not implemented yet"
            end if
        end if
    end if
    !> New-step reaction-mixing term method (downstream estimation option), selected once
    if (this%r_down_opt==1) then
        p_r_tilde=>compute_r_tilde_impl_opt1_idx
    else if (this%r_down_opt==2) then
        p_r_tilde=>compute_r_tilde_impl_opt2_idx
    else if (this%r_down_opt==3) then
        p_r_tilde=>compute_r_tilde_impl_opt3_idx
    else if (this%r_down_opt==4) then
        p_r_tilde=>compute_r_tilde_impl_opt4_idx
    else
        error stop "r down option not implemented yet"
    end if

!> -----------------------------------------------------------------------
!> Pre-compute canonical/non-canonical vectors and lambdas
    call mixing_ratios_conc%get_can_vec(this%CV_params%abs_tol,num_can_vec,ind_can_vec,num_non_can_vec,ind_non_can_vec)
    if (this%r_down_opt==4) then
        lambdas_R=mixing_ratios_R%compute_mix_ratios_R_opt4(mix_conc_indices)
    else
        allocate(lambdas_R(mixing_ratios_R%num_cols))
        do m=1,mixing_ratios_R%num_cols
            lambdas_R(m)=mixing_ratios_R%cols(m)%col_1(1)
        end do
    end if

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
    end do

!> Pre-allocate reusable buffer
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

!> -----------------------------------------------------------------------
!> Initialize profiling accumulators
    cpu_update_old_total = 0d0
    cpu_mix_total = 0d0
    cpu_r_tilde_total = 0d0
    cpu_solver_total = 0d0
    cpu_derived_total = 0d0
    cpu_eq_rates_total = 0d0
    cpu_solid_gas_total = 0d0
    cpu_output_total = 0d0
    cpu_snapshot_total = 0d0
    cpu_data_prep_total = 0d0
    cpu_getters_total = 0d0
    cpu_cum_rxn_amnt_kin_total = 0d0
    cpu_c_hat_total = 0d0

!> Open CPU profiling output file
    if (this%lump_flag) then
        cpu_root = root//'_Euler_lumped'
    else
        cpu_root = root//'_Euler_consistent'
    end if
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
        !> Update old attributes
        call cpu_time(cpu_prof_t0)
        do i=1,this%num_target_waters
            call this%waters(this%tar_wat_indices(i))%update_attributes() !< Shift all old attributes
        end do
        call cpu_time(cpu_prof_t1)
        cpu_update_old_total = cpu_update_old_total + (cpu_prof_t1 - cpu_prof_t0)
        !> Target waters loop
        do i=1,num_non_can_vec
            idx = this%tar_wat_indices(ind_non_can_vec(i))
            n_p = n_p_cache(i)
            n_v = n_v_cache(i)
            n_eq = n_eq_cache(i)
            !> Data preparation
            call cpu_time(cpu_prof_t0)
            if (mix_conc_indices%cols(ind_non_can_vec(i))%col_1(1)/=idx) then
                print *, mix_conc_indices%cols(ind_non_can_vec(i))%col_1(1), idx
                error stop "Target waters not in the right order"
            end if
            allocate(c_hat(n_v))
            allocate(r_tilde(n_v))
            call cpu_time(cpu_prof_t1)
            cpu_data_prep_total = cpu_data_prep_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Conservative mixing: page-4 two-stage consistent mix -> c_hat_old (stage 2 within-zone
            !> mixing + speciation, stage 3 cross-zone mixing in aqueous-phase space) and its purely
            !> conservative part c_mix, also used by the domain mass-balance diagnostic.
            call cpu_time(cpu_prof_t0)
            associate(mc => mix_conc_indices%cols(ind_non_can_vec(i)), &
                      rc_col => mixing_ratios_conc%cols(ind_non_can_vec(i)), &
                      rr_col => mixing_ratios_R%cols(ind_non_can_vec(i)))
                call p_mixer(this%waters, this%zone_waters(ind_non_can_vec(i),:), &
                    mc%col_1(1:mc%dim-2), rc_col%col_1(1:mc%dim-2), rr_col%col_1(1:mc%dim-2), &
                    this%wat_indices_zones, this%ext_waters_indices, Delta_t, theta_r, c_hat_old, c_mix, .true.)
            end associate
            call cpu_time(cpu_prof_t1)
            cpu_mix_total = cpu_mix_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Reacted concentration c_hat = c_hat_old + Delta_t*r_tilde. r_tilde is the NEW-STEP (theta)
            !> term only; the previous-step (1-theta) term is already folded into c_hat_old by p_mixer.
            call cpu_time(cpu_prof_t0)
            associate(mc => mix_conc_indices%cols(ind_non_can_vec(i)), &
                      rr_col => mixing_ratios_R%cols(ind_non_can_vec(i)))
                call p_r_tilde(this%waters, mc%col_1(2:mc%dim-2), &
                    this%waters(idx)%solid_chemistry%reactive_zone%ind_var_act_species, &
                    rr_col%col_1(2:mc%dim-2), mc%col_1(mc%dim-1), mc%col_1(mc%dim), theta_r, r_tilde, &
                    curr_water_idx=mc%col_1(1), lambda_self=rr_col%col_1(1), &
                    wat_indices_zones=this%wat_indices_zones)
                !> Positivity-preserving update (same uniform-omega limiter as the mixers)
                call add_limited_increment(c_hat_old,Delta_t*r_tilde,c_hat)
            end associate
            call this%waters(idx)%set_conc_hat(c_hat) !> aqueous blocks of the Newton speciation input, for the conc_hat attribute
            call cpu_time(cpu_prof_t1)
            cpu_r_tilde_total = cpu_r_tilde_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Getters (get_conc_var_act)
            call cpu_time(cpu_prof_t0)
            conc_nc=this%waters(idx)%get_conc_var_act()
            call cpu_time(cpu_prof_t1)
            cpu_getters_total = cpu_getters_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Reactive mixing solver
            call cpu_time(cpu_prof_t0)
            call p_solver(this%waters(idx),this%waters(idx)%get_conc_old_old_prim_new_rz(),c_hat,&
                lambdas_R(ind_non_can_vec(i)),Delta_t,theta_r,conc_nc,conc_comp(1:n_p))
            call cpu_time(cpu_prof_t1)
            cpu_solver_total = cpu_solver_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Mass action law check for equilibrium reactions (ideal solution)
            if (n_eq > 0) then
                mal_residual = matmul(this%waters(idx)%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star, &
                    log10(conc_nc(1:n_p))) + &
                    this%waters(idx)%logK_star - &
                    log10(conc_nc(n_p+1:n_p+n_eq))
                mal_max = maxval(abs(mal_residual))
                if (mal_max > 1d-6) then
                    print *, "WARNING: mass action law violated after p_solver, water", idx, &
                        " time step", k, " max|residual|=", mal_max
                    do j = 1, n_eq
                        if (abs(mal_residual(j)) > 1d-6) then
                            print *, "  eq reaction", j, " residual=", mal_residual(j), &
                                " log10(c2)=", log10(conc_nc(n_p+j)), &
                                " expected=", log10(conc_nc(n_p+j)) - mal_residual(j)
                        end if
                    end do
                end if
                deallocate(mal_residual)
            end if
            !> Derived quantities (pH, salinity, ionic strength)
            call cpu_time(cpu_prof_t0)
            call this%waters(idx)%compute_pH()
            call this%waters(idx)%compute_saturation_mins() !< saturation indices of all minerals at the end-of-step state
            call this%waters(idx)%compute_salinity()
            call this%waters(idx)%compute_ionic_strength()
            call cpu_time(cpu_prof_t1)
            cpu_derived_total = cpu_derived_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Equilibrium reaction amounts
            call cpu_time(cpu_prof_t0)
            if (n_eq>0) then
                call this%waters(idx)%compute_Re(c_hat(1:n_v),Delta_t,&
                    lambdas_R(ind_non_can_vec(i)),theta_r)
            end if
            call cpu_time(cpu_prof_t1)
            cpu_eq_rates_total = cpu_eq_rates_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Domain-integrated mass-balance contribution of this cell.
            !> Must use the SAME self reaction mixing ratio as compute_Re/p_solver (lambdas_R),
            !> otherwise the reaction term is reconstructed with a mismatched lambda and the
            !> secondary species (so4-2, co3-2) residual blows up at the reaction front.
            call this%waters(idx)%get_react_mix_mass_bal_terms(c_hat_old(1:n_v), c_hat(1:n_v), &
                lambdas_R(ind_non_can_vec(i)), Delta_t, theta_r, &
                mb_lhs_c, mb_transp_c, mb_react_c)
            if (size(mb_lhs_c) == n_v_ref) then
                mb_lhs = mb_lhs + mb_lhs_c
                mb_transp = mb_transp + mb_transp_c
                mb_react = mb_react + mb_react_c
            end if
            deallocate(mb_lhs_c, mb_transp_c, mb_react_c)
            !> Solid/gas chemistry
            call cpu_time(cpu_prof_t0)
            if (has_minerals_flag(i)) then
                call this%waters(idx)%solid_chemistry%compute_mass_bal_mins(Delta_t)
                call this%waters(idx)%solid_chemistry%compute_conc_minerals_iter(Delta_t,lambdas_R(ind_non_can_vec(i)))
            end if
            if (has_gases_flag(i)) then
                call this%waters(idx)%gas_chemistry%compute_conc_gases_iter(Delta_t,&
                    this%waters(idx)%volume,&
                    [this%waters(idx)%rxn_rate_eq,this%waters(idx)%rxn_rate_kin_avg])
                call this%waters(idx)%gas_chemistry%compute_vol_gas_species_conc()
                call this%waters(idx)%gas_chemistry%compute_log_act_coeffs_gases()
            end if
            call cpu_time(cpu_prof_t1)
            cpu_solid_gas_total = cpu_solid_gas_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Accumulate reaction amounts
            call cpu_time(cpu_prof_t0)
            this%waters(idx)%cum_rxn_amnt_kin=this%waters(idx)%cum_rxn_amnt_kin+this%waters(idx)%rxn_amnt_kin
            this%waters(idx)%solid_chemistry%cum_rxn_amnt_kin=&
                this%waters(idx)%solid_chemistry%cum_rxn_amnt_kin+this%waters(idx)%solid_chemistry%rxn_amnt_kin
            call cpu_time(cpu_prof_t1)
            cpu_cum_rxn_amnt_kin_total = cpu_cum_rxn_amnt_kin_total + (cpu_prof_t1 - cpu_prof_t0)
            deallocate(c_hat, r_tilde, c_hat_old, c_mix, conc_nc)
        end do
        !> Verify the surface-complex site-balance invariants at every cell this time step
        call check_surf_compl_act_sum(this, k, time)
        !> Domain-integrated mass balance for this time step (per variable-activity species)
        call this%write_react_mix_mass_bal_dom(idx_ref,iunit_mb,time,n_v_ref,mb_lhs,mb_transp,mb_react)
        !> Output writing honouring chem_out_options%vars (conc / react_rate / vol_fract)
        call cpu_time(cpu_prof_t0)
        !> Sum equilibrium reaction rates at interfaces if mix_by_zones_flag is set
        if (this%mix_by_zones_flag) then
            do i=1,num_non_can_vec
                call this%waters(this%tar_wat_indices(ind_non_can_vec(i)))%sum_eq_rxn_rates(this%zone_waters(&
                    this%tar_wat_indices(ind_non_can_vec(i)),:))
            end do
        end if
        if (k==this%chem_out_options%time_steps(kk)) then
            call write_out_step_vars(this, unit, time)
            if (kk<this%chem_out_options%num_time_steps) kk=kk+1
        end if
        call cpu_time(cpu_prof_t1)
        cpu_output_total = cpu_output_total + (cpu_prof_t1 - cpu_prof_t0)
        !> Write chemistry snapshot at output time steps
        !> (Snapshots disabled — write_chemistry is called once at the end via write_RT)
    end do

!> -----------------------------------------------------------------------
!> CPU profiling report
    call cpu_time(cpu_total_t1)
    cpu_total_elapsed = cpu_total_t1 - cpu_total_t0
    write(iunit_cpu,*) '====================================================================================='
    write(iunit_cpu,*) '  CPU PROFILING REPORT: solve_reactive_mixing_ideal_cons'
    write(iunit_cpu,*) '====================================================================================='
    write(iunit_cpu,'(A,I10)')     '  Total time steps:       ', num_time_loc
    write(iunit_cpu,'(A,I10)')     '  Non-canonical waters:   ', num_non_can_vec
    write(iunit_cpu,*) '-------------------------------------------------------------------------------------'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Update old attributes:  ', cpu_update_old_total, ' s  (', cpu_update_old_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Data preparation:       ', cpu_data_prep_total,  ' s  (', cpu_data_prep_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Conservative mixing:    ', cpu_mix_total,        ' s  (', cpu_mix_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Reaction rate (r_tilde):', cpu_r_tilde_total,    ' s  (', cpu_r_tilde_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Reactive mixing solver: ', cpu_solver_total,     ' s  (', cpu_solver_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Derived quantities:     ', cpu_derived_total,    ' s  (', cpu_derived_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Equilibrium rates:      ', cpu_eq_rates_total,   ' s  (', cpu_eq_rates_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Solid/gas chemistry:    ', cpu_solid_gas_total,  ' s  (', cpu_solid_gas_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Output writing:         ', cpu_output_total,     ' s  (', cpu_output_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Snapshot (write_chem):  ', cpu_snapshot_total,   ' s  (', cpu_snapshot_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Getters (conc_nc):      ', cpu_getters_total,    ' s  (', cpu_getters_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  c_hat limiter:          ', cpu_c_hat_total, ' s  (', cpu_c_hat_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Rk accumulation:        ', cpu_cum_rxn_amnt_kin_total,   ' s  (', cpu_cum_rxn_amnt_kin_total/60d0, ' min)'
    write(iunit_cpu,*) '-------------------------------------------------------------------------------------'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  TOTAL SUBROUTINE TIME:  ', cpu_total_elapsed,    ' s  (', cpu_total_elapsed/60d0, ' min)'
    write(iunit_cpu,*) '====================================================================================='
    if (cpu_total_elapsed > 0.0d0) then
        write(iunit_cpu,*) ''
        write(iunit_cpu,*) '  Breakdown (% of total):'
        write(iunit_cpu,'(A,F8.2,A)') '    Update old attributes:   ', 100d0*cpu_update_old_total/cpu_total_elapsed, ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Data preparation:        ', 100d0*cpu_data_prep_total/cpu_total_elapsed,  ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Conservative mixing:     ', 100d0*cpu_mix_total/cpu_total_elapsed,        ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Reaction rate (r_tilde): ', 100d0*cpu_r_tilde_total/cpu_total_elapsed,    ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Reactive mixing solver:  ', 100d0*cpu_solver_total/cpu_total_elapsed,     ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Derived quantities:      ', 100d0*cpu_derived_total/cpu_total_elapsed,    ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Equilibrium rates:       ', 100d0*cpu_eq_rates_total/cpu_total_elapsed,   ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Solid/gas chemistry:     ', 100d0*cpu_solid_gas_total/cpu_total_elapsed,  ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Output writing:          ', 100d0*cpu_output_total/cpu_total_elapsed,     ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Snapshot (write_chem):   ', 100d0*cpu_snapshot_total/cpu_total_elapsed,   ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Getters (conc_nc):       ', 100d0*cpu_getters_total/cpu_total_elapsed,    ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    c_hat limiter:           ', 100d0*cpu_c_hat_total/cpu_total_elapsed, ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Rk accumulation:         ', 100d0*cpu_cum_rxn_amnt_kin_total/cpu_total_elapsed,   ' %'
        write(iunit_cpu,*) '====================================================================================='
    end if
    write(iunit_cpu,*) ''
    write(iunit_cpu,'(A,F12.6,A)') '  Avg conservative mixing per time step:       ', cpu_mix_total / dble(num_time_loc), ' s'
    write(iunit_cpu,'(A,F12.6,A)') '  Avg reactive solver per time step:           ', cpu_solver_total / dble(num_time_loc), ' s'
    write(iunit_cpu,'(A,F12.6,A)') '  Avg r_tilde per time step:                   ', cpu_r_tilde_total / dble(num_time_loc), ' s'
    if (num_non_can_vec > 0) then
        write(iunit_cpu,'(A,F12.6,A)') '  Avg conservative mixing per step/water:  ', cpu_mix_total / dble(num_time_loc * num_non_can_vec), ' s'
        write(iunit_cpu,'(A,F12.6,A)') '  Avg reactive solver per step/water:      ', cpu_solver_total / dble(num_time_loc * num_non_can_vec), ' s'
        write(iunit_cpu,'(A,F12.6,A)') '  Avg r_tilde per step/water:              ', cpu_r_tilde_total / dble(num_time_loc * num_non_can_vec), ' s'
    end if
    write(iunit_cpu,*) ''
    write(iunit_cpu,'(A,F12.6,A)') '  Avg getters per time step:                   ', cpu_getters_total / dble(num_time_loc), ' s'
    write(iunit_cpu,'(A,F12.6,A)') '  Avg c_hat limiter per time step:             ', cpu_c_hat_total / dble(num_time_loc), ' s'
    write(iunit_cpu,'(A,F12.6,A)') '  Avg Rk accum per time step:                  ', cpu_cum_rxn_amnt_kin_total / dble(num_time_loc), ' s'
    close(iunit_cpu)
    print *, 'CPU profiling written to: ', dir//cpu_root//'_cpu_profile.out'

!> Close output file and clean up
    close(unit)
    close(iunit_mb)
    deallocate(mb_lhs, mb_transp, mb_react)
    deallocate(conc_comp)
    deallocate(lambdas_R)
    deallocate(n_p_cache, n_v_cache, n_eq_cache)
    deallocate(has_minerals_flag, has_gases_flag)
end subroutine solve_reactive_mixing_ideal_cons
