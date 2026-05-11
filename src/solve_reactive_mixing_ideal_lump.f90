!> @brief Solves reactive mixing problem assuming ideal conditions with lumped formulation
!> @details Computes concentrations, activities, activity coefficients, reaction rates
!> and volumetric fractions of minerals (if present) using operator splitting.
subroutine solve_reactive_mixing_ideal_lump(this,dir,root,mixing_ratios_conc,mixing_ratios_R,&
    mix_conc_indices,mix_react_indices,time_discr,theta_r)
    use chemistry_m, only: chemistry_c
    use aqueous_chemistry_m, only: aqueous_chemistry_c, &
        reactive_mixing_iter_EI_eq_anal_ideal, &
        compute_c_mix_global, reactive_mixing_iter_EI_kin_anal_ideal_opt2, &
        reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2
    use arrays_m, only: real_array_c, int_array_c
    use time_discr_m, only: time_discr_c
    implicit none
!> Arguments
    class(chemistry_c) :: this
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: root
    class(real_array_c), intent(in) :: mixing_ratios_conc
    class(real_array_c), intent(in) :: mixing_ratios_R
    class(int_array_c), intent(in) :: mix_conc_indices
    class(int_array_c), intent(in) :: mix_react_indices
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
!> Global caches for compute_c_mix_global
    real(kind=8), allocatable :: all_conc_old(:,:)
    integer(kind=4), allocatable :: all_ind_aq_sp(:,:)
    integer(kind=4) :: w, max_n_conc, max_n_idx
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

!> -----------------------------------------------------------------------
!> Initialization
    call cpu_time(cpu_total_t0)
    n_eq=0
    do i=1,this%num_target_waters
        n_eq=n_eq+this%waters(this%tar_wat_indices(i))%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    end do

    time=0d0
    kk=2
    ii=1
    kk_snap=2
    num_time_loc = time_discr%Num_time
    this%chem_out_options%time_steps(this%chem_out_options%num_time_steps)=num_time_loc
    do i=2,this%chem_out_options%num_time_steps-1
        if (this%chem_out_options%time_steps(i) > num_time_loc) then
            print *, "Output time step index", this%chem_out_options%time_steps(i), &
                "exceeds total number of time steps", num_time_loc
            error stop "Mismatch between out_opts.dat and discr_temp.dat: output time step exceeds total time steps"
        end if
    end do
    unit=7
    open(unit,file=dir//root//'.output',form="formatted",access="sequential",status="unknown")
    idx = this%chem_out_options%ind_waters(1)
    write(unit,"(A10)", advance='no') 'Water'
    do j=1,this%waters(idx)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
        if (this%chem_syst%species(this%waters(idx)%ind_var_act_species(j))%cst_act_flag .eqv. .false.) then
            write(unit,"(A15)", advance='no') this%chem_syst%species(this%waters(idx)%ind_var_act_species(j))%name
        end if
    end do
    write(unit,*)

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
    allocate(lumped_lambdas(mixing_ratios_R%num_cols))
    do i=1,mixing_ratios_R%num_cols
        lumped_lambdas(i)=sum(mixing_ratios_R%cols(i)%col_1)
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

!> Pre-allocate global caches for compute_c_mix_global
    max_n_conc = 0
    max_n_idx = 0
    do w = 1, this%num_waters
        j = size(this%waters(w)%conc_old)
        if (j > max_n_conc) max_n_conc = j
        j = this%waters(w)%aq_phase%num_species
        if (j > max_n_idx) max_n_idx = j
    end do
    allocate(all_conc_old(max_n_conc, this%num_waters))
    allocate(all_ind_aq_sp(max_n_idx, this%num_waters))
    do w = 1, this%num_waters
        all_ind_aq_sp(1:size(this%waters(w)%indices_aq_species), w) = &
            this%waters(w)%indices_aq_species
    end do

!> -----------------------------------------------------------------------
!> Initialize profiling accumulators
    cpu_solver_total = 0d0
    cpu_output_total = 0d0
    cpu_snapshot_total = 0d0

!> Open CPU profiling output file
    open(newunit=iunit_cpu, file=dir//root//'_cpu_profile.out', status='replace', action='write')

!> -----------------------------------------------------------------------
!> Main time loop
    do k=1,num_time_loc
        Delta_t=time_discr%get_Delta_t(k)
        time=time+Delta_t
        if (k==this%chem_out_options%time_steps(kk)) then
           write(unit,"(2x,'t = ',*(ES15.5))") time
        end if
        !> Note: all_conc_old is refreshed inside the iter routine, after
        !> update_old_attributes has copied the previous step's new conc into conc_old.
        !> Call reactive mixing iteration subroutine
        call cpu_time(cpu_prof_t0)
        call this%solve_reactive_mixing_ideal_lump_iter(k, Delta_t, theta_r, &
            num_non_can_vec, ind_non_can_vec, &
            n_p_cache, n_v_cache, n_eq_cache, &
            has_minerals_flag, has_gases_flag, &
            mix_conc_indices, mixing_ratios_conc, &
            lumped_lambdas, &
            all_conc_old, all_ind_aq_sp, &
            conc_comp, &
            p_solver)
        call cpu_time(cpu_prof_t1)
        cpu_solver_total = cpu_solver_total + (cpu_prof_t1 - cpu_prof_t0)
        !> Output writing
        call cpu_time(cpu_prof_t0)
        if (k==this%chem_out_options%time_steps(kk)) then
            do i=1,num_non_can_vec
                idx = this%tar_wat_indices(ind_non_can_vec(i))
                n_v = n_v_cache(i)
                if (idx==this%chem_out_options%ind_waters(ii)) then
                    conc_nc = this%waters(idx)%get_conc_nc()
                    write(unit,"(I10)", advance='no') this%chem_out_options%ind_waters(ii)
                    do j=1,n_v
                        !> conc_nc holds variable-activity species only — all entries are non-cst-act by construction
                        write(unit,"(ES15.5)", advance='no') conc_nc(j)
                    end do
                    write(unit,*)
                    deallocate(conc_nc)
                    if (ii<this%chem_out_options%num_waters) then
                        ii=ii+1
                    else if (kk<this%chem_out_options%num_time_steps) then
                        kk=kk+1
                        ii=1
                    else
                        exit
                    end if
                end if
            end do
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
    print *, 'CPU profiling written to: ', dir//root//'_cpu_profile.out'

!> Close output file and clean up
    close(unit)
    deallocate(conc_comp)
    deallocate(lumped_lambdas)
    deallocate(all_conc_old, all_ind_aq_sp)
    deallocate(n_p_cache, n_v_cache, n_eq_cache)
    deallocate(has_minerals_flag, has_gases_flag)
end subroutine solve_reactive_mixing_ideal_lump