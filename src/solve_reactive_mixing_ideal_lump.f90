!> @brief Solves reactive mixing problem assuming ideal conditions with lumped formulation
!> @details Computes concentrations, activities, activity coefficients, reaction rates
!> and volumetric fractions of minerals (if present) using operator splitting.
subroutine solve_reactive_mixing_ideal_lump(this,dir,root,mixing_ratios_conc,mixing_ratios_R,&
    mix_conc_indices,time_discr,theta_r)
    use chemistry_m, only: chemistry_c, write_out_step_vars, check_surf_compl_act_sum
    use aqueous_chemistry_m, only: aqueous_chemistry_c, &
        reactive_mixing_iter_EI_eq_anal_ideal, &
        compute_c_mix_by_zones_lump, compute_c_mix_global_lump, reactive_mixing_iter_EI_kin_anal_ideal_opt2, &
        reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2
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
    integer(kind=4) :: i, j, k, kk, ii
    integer(kind=4) :: idx          !< Cached target water index: tar_wat_indices(ind_non_can_vec(i))
!> Problem dimensions
    integer(kind=4) :: n_p, n_v, n_eq
    integer(kind=4) :: num_can_vec, num_non_can_vec
    integer(kind=4) :: unit         !< Fortran I/O unit number for output file
    integer(kind=4) :: iunit_cpu    !< Fortran I/O unit number for CPU profiling output
    character(len=:), allocatable :: cpu_root !< root + transport/WMA tags for the CPU profiling filename
    integer(kind=4) :: num_time_loc !< Cached total number of time steps
    integer(kind=4) :: kk_snap      !< Counter for chemistry snapshot output time steps
!> Index arrays
    integer(kind=4), allocatable :: ind_can_vec(:)
    integer(kind=4), allocatable :: ind_non_can_vec(:)
!> Time variables
    real(kind=8) :: time, Delta_t
!> Concentration and reaction rate arrays
    real(kind=8), allocatable :: conc_nc(:)
    real(kind=8), allocatable :: conc_comp(:)
    real(kind=8), allocatable :: lumped_lambdas(:)
!> Domain-integrated reactive-mixing mass balance (per variable-activity species)
    integer(kind=4) :: iunit_mb !< mass-balance output file unit
    integer(kind=4) :: idx_ref !< reference target-water index (species names/dimension)
    integer(kind=4) :: n_v_ref !< reference number of variable-activity species
    real(kind=8), allocatable :: c_mix(:) !< per-zone transport mixed concentrations
    real(kind=8), allocatable :: mb_lhs(:), mb_transp(:), mb_react(:) !< domain-summed terms
    real(kind=8), allocatable :: mb_lhs_c(:), mb_transp_c(:), mb_react_c(:) !< per-cell terms
!> Pre-computed caches
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
    real(kind=8) :: cpu_snapshot_total
!> Procedure pointer
    procedure(reactive_mixing_iter_EI_kin_anal_ideal_opt2), pointer :: p_solver=>null()
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
    ii=1
    kk_snap=2
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
!> Solver selection
    if (this%chem_syst%num_kin_reacts>0 .and. n_eq>0) then
        if (this%Jac_opt==1) then
            if (this%rk_avg_opt==2) then
                p_solver=>reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2
            else
                error stop "rk average option not implemented yet"
            end if
        end if
    else if (n_eq>0) then
        if (theta_r>0d0) then
            p_solver=>reactive_mixing_iter_EI_eq_anal_ideal
        end if
    else
        if (this%Jac_opt==1) then
            if (this%rk_avg_opt==2) then
                p_solver=>reactive_mixing_iter_EI_kin_anal_ideal_opt2
            else
                error stop "rk average option not implemented yet"
            end if
        end if
    end if

!> -----------------------------------------------------------------------
!> Pre-compute canonical/non-canonical vectors and lumped lambdas
    call mixing_ratios_conc%get_can_vec(this%CV_params%abs_tol,num_can_vec,ind_can_vec,num_non_can_vec,ind_non_can_vec)
    allocate(lumped_lambdas(mixing_ratios_conc%num_cols))
    do i=1,mixing_ratios_conc%num_cols
        lumped_lambdas(i)=sum(mixing_ratios_conc%cols(i)%col_1) !< Concentration mixing ratios sum to 1 (mass conservation)
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
    cpu_solver_total = 0d0
    cpu_output_total = 0d0
    cpu_snapshot_total = 0d0

!> Open CPU profiling output file. The transport approach (Euler — this is the Eulerian
!> reactive-mixing path) and WMA option (lumped/consistent) are encoded in the filename,
!> matching the {root}{_Lagr|_Euler}{_lumped|_consistent}_*.out scheme.
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
        !> Call reactive mixing iteration subroutine with the selected mixing procedure
        call cpu_time(cpu_prof_t0)
        call this%solve_reactive_mixing_ideal_lump_iter(k, Delta_t, theta_r, time, &
            num_non_can_vec, ind_non_can_vec, &
            n_p_cache, n_v_cache, n_eq_cache, &
            has_minerals_flag, has_gases_flag, &
            mix_conc_indices, mixing_ratios_conc, &
            mixing_ratios_R, &
            lumped_lambdas, &
            conc_comp, iunit_mb, &
            p_solver, p_mixer)
        call cpu_time(cpu_prof_t1)
        cpu_solver_total = cpu_solver_total + (cpu_prof_t1 - cpu_prof_t0)
        !> Verify the surface-complex site-balance invariants at every cell this time step
        call check_surf_compl_act_sum(this, k, time)
        !> Output writing honouring chem_out_options%vars (conc / react_rate / vol_fract)
        call cpu_time(cpu_prof_t0)
        if (k==this%chem_out_options%time_steps(kk)) then
            call write_out_step_vars(this, unit, time)
            if (kk<this%chem_out_options%num_time_steps) kk=kk+1
        end if
        call cpu_time(cpu_prof_t1)
        cpu_output_total = cpu_output_total + (cpu_prof_t1 - cpu_prof_t0)
        !> (Domain-integrated mass balance is accumulated and written inside the iter routine,
        !>  where the per-reactive-zone c_mix is already computed — avoids recomputing it here.)
        !> Write chemistry snapshot at output time steps
        !> (Snapshots disabled — write_chemistry is called once at the end via write_RT)
    end do

!> -----------------------------------------------------------------------
!> CPU profiling report
    call cpu_time(cpu_total_t1)
    cpu_total_elapsed = cpu_total_t1 - cpu_total_t0
    write(iunit_cpu,*) '====================================================================================='
    write(iunit_cpu,*) '  CPU PROFILING REPORT: solve_reactive_mixing_ideal_lump'
    write(iunit_cpu,*) '====================================================================================='
    write(iunit_cpu,'(A,I10)')     '  Total time steps:       ', num_time_loc
    write(iunit_cpu,'(A,I10)')     '  Non-canonical waters:   ', num_non_can_vec
    write(iunit_cpu,*) '-------------------------------------------------------------------------------------'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Iteration (total):      ', cpu_solver_total,     ' s  (', cpu_solver_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Output writing:         ', cpu_output_total,     ' s  (', cpu_output_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Snapshot (write_chem):  ', cpu_snapshot_total,   ' s  (', cpu_snapshot_total/60d0, ' min)'
    write(iunit_cpu,*) '-------------------------------------------------------------------------------------'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  TOTAL SUBROUTINE TIME:  ', cpu_total_elapsed,    ' s  (', cpu_total_elapsed/60d0, ' min)'
    write(iunit_cpu,*) '====================================================================================='
    if (cpu_total_elapsed > 0.0d0) then
        write(iunit_cpu,*) ''
        write(iunit_cpu,*) '  Breakdown (% of total):'
        write(iunit_cpu,'(A,F8.2,A)') '    Iteration (total):       ', 100d0*cpu_solver_total/cpu_total_elapsed,     ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Output writing:          ', 100d0*cpu_output_total/cpu_total_elapsed,     ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Snapshot (write_chem):   ', 100d0*cpu_snapshot_total/cpu_total_elapsed,   ' %'
        write(iunit_cpu,*) '====================================================================================='
    end if
    write(iunit_cpu,*) ''
    write(iunit_cpu,'(A,F12.6,A)') '  Avg iteration per time step:                 ', cpu_solver_total / dble(num_time_loc), ' s'
    if (num_non_can_vec > 0) then
        write(iunit_cpu,'(A,F12.6,A)') '  Avg iteration per step/water:            ', cpu_solver_total / dble(num_time_loc * num_non_can_vec), ' s'
    end if
    close(iunit_cpu)
    print *, 'CPU profiling written to: ', dir//cpu_root//'_cpu_profile.out'

!> Close output file and clean up
    close(unit)
    close(iunit_mb)
    deallocate(conc_comp)
    deallocate(lumped_lambdas)
    deallocate(n_p_cache, n_v_cache, n_eq_cache)
    deallocate(has_minerals_flag, has_gases_flag)
end subroutine solve_reactive_mixing_ideal_lump