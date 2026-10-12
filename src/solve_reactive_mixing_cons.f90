!> @file solve_reactive_mixing_cons.f90
!> @brief Solves the NON-IDEAL conservative (consistent) reactive-mixing problem.
!> @details Non-ideal counterpart of solve_reactive_mixing_ideal_cons. It shares the same
!> optimized structure: global conservative mixing via compute_c_mix_global (with a
!> per-step cache of every water's conc_old), per-target-water dimension/flag caches,
!> NO full deep copies of the target-water arrays (state is advanced in place on
!> this%waters via update_attributes), consistent reaction-rate mixing via
!> compute_r_tilde, a positivity-preserving reactive update via compute_c_hat_patankar,
!> and CPU profiling. The equilibrium speciation uses the NON-IDEAL mixing_iter_comp
!> iteration (which also finalises pH/salinity/ionic strength and non-ideal activity
!> coefficients) and a non-ideal (activity-coefficient-aware) mass-action-law check.
!> Only equilibrium-only systems are implemented; kinetic reactions trigger an error stop.
!>
!> @author Jordi Petchamé-Guerrero
!> @date November 2025
subroutine solve_reactive_mixing_cons(this,dir,root,mixing_ratios_conc,mixing_ratios_R,mix_conc_indices,&
        time_discr,theta_r)
    use chemistry_m, only: chemistry_c, write_out_step_vars
    !> Non-ideal reactive-mixing iteration, global conservative mixing, consistent r_tilde and c_hat limiter
    use aqueous_chemistry_m, only: aqueous_chemistry_c, mixing_iter_comp, compute_c_mix_global_cons, &
        compute_c_mix_by_zones_cons, &
        compute_r_tilde_impl_opt1_idx, compute_r_tilde_impl_opt2_idx, compute_r_tilde_impl_opt3_idx, &
        compute_r_tilde_impl_opt4_idx, compute_c_hat_patankar
    use arrays_m, only: real_array_c, int_array_c
    use time_discr_m, only: time_discr_c
    implicit none
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: dir !> directory for output files
    character(len=*), intent(in) :: root !> root name for output files
    class(real_array_c), intent(in) :: mixing_ratios_conc !> mixing ratios matrix for concentrations
    class(real_array_c), intent(in) :: mixing_ratios_R !> mixing ratios matrix for reaction amounts (mix_conc_indices layout, 0 for conc-only waters)
    class(int_array_c), intent(in) :: mix_conc_indices !> indices of the waters that mix with each target water (concentrations and reaction rates)
    class(time_discr_c), intent(in) :: time_discr !> time discretisation object
    real(kind=8), intent(in) :: theta_r !> time weighting factor for chemical reactions (0=explicit, 1=implicit)
!> Loop counters and indices
    integer(kind=4) :: i, j, k, kk, ii, m
    integer(kind=4) :: idx !< Absolute target-water index tar_wat_indices(ind_non_can_vec(i))
!> Problem dimensions
    integer(kind=4) :: n_p, n_v, n_eq
    integer(kind=4) :: num_can_vec, num_non_can_vec
    integer(kind=4) :: num_mix_loc
    integer(kind=4) :: unit !< Output file unit
    integer(kind=4) :: iunit_cpu !< CPU-profiling file unit
    character(len=:), allocatable :: cpu_root !< root + Euler/consistent tag for the CPU-profiling filename
    integer(kind=4) :: num_time_loc !< Cached total number of time steps
!> Index arrays
    integer(kind=4), allocatable :: ind_can_vec(:)
    integer(kind=4), allocatable :: ind_non_can_vec(:)
!> Time variables
    real(kind=8) :: time, Delta_t
!> Concentration and reaction-rate arrays
    real(kind=8), allocatable :: c_hat(:)
    real(kind=8), allocatable :: r_tilde(:)
    real(kind=8), allocatable :: c_mix(:)
    real(kind=8), allocatable :: c_tilde(:)
    real(kind=8), allocatable :: conc_nc(:)
    real(kind=8), allocatable :: conc_comp(:)
    real(kind=8), allocatable :: lambdas_R(:)
    real(kind=8), allocatable :: mal_residual(:)
    real(kind=8) :: mal_max
!> Domain-integrated reactive-mixing mass balance (per variable-activity species)
    integer(kind=4) :: iunit_mb !< mass-balance output file unit
    integer(kind=4) :: idx_ref !< reference target-water index (species names/dimension)
    integer(kind=4) :: n_v_ref !< reference number of variable-activity species
    real(kind=8), allocatable :: mb_lhs(:), mb_transp(:), mb_react(:) !< domain-summed terms
    real(kind=8), allocatable :: mb_lhs_c(:), mb_transp_c(:), mb_react_c(:) !< per-cell terms
!> Global caches for compute_c_mix_global
    real(kind=8), allocatable :: all_conc_old(:,:)
    integer(kind=4), allocatable :: all_ind_aq_sp(:,:)
    integer(kind=4) :: w, max_n_conc, max_n_idx
!> Pre-computed per-target-water caches
    logical, allocatable :: has_minerals_flag(:)
    logical, allocatable :: has_gases_flag(:)
    integer(kind=4), allocatable :: n_p_cache(:)
    integer(kind=4), allocatable :: n_v_cache(:)
    integer(kind=4), allocatable :: n_eq_cache(:)
    integer(kind=4) :: max_n_p
!> CPU profiling variables
    real(kind=8) :: cpu_total_t0, cpu_total_t1, cpu_total_elapsed
    real(kind=8) :: cpu_prof_t0, cpu_prof_t1
    real(kind=8) :: cpu_mix_total
    real(kind=8) :: cpu_r_tilde_total
    real(kind=8) :: cpu_solver_total
    real(kind=8) :: cpu_solid_gas_total
    real(kind=8) :: cpu_output_total
!> Procedure pointers (non-ideal iteration + consistent reaction-rate mixing)
    procedure(mixing_iter_comp), pointer :: p_solver=>null()
    procedure(compute_r_tilde_impl_opt4_idx), pointer :: p_r_tilde=>null()
    procedure(compute_c_mix_by_zones_cons), pointer :: p_mixer

!> -----------------------------------------------------------------------
!> Initialization
    call cpu_time(cpu_total_t0)
    if (this%mix_by_zones_flag) then
        p_mixer => compute_c_mix_by_zones_cons
    else
        p_mixer => compute_c_mix_global_cons
    end if
    write(0,*) '[cons] Z0 entered'; flush(0)
    n_eq=0
    do i=1,this%num_target_waters
        n_eq=n_eq+this%waters(this%tar_wat_indices(i))%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    end do

    time=0d0
    kk=2
    ii=1
    num_time_loc=time_discr%Num_time
    do i=2,this%chem_out_options%num_time_steps-1
        if (this%chem_out_options%time_steps(i) > num_time_loc) then
            print *, "Output time step index", this%chem_out_options%time_steps(i), &
                "exceeds total number of time steps", num_time_loc
            error stop "Mismatch between out_opts.dat and discr_temp.dat: output time step exceeds total time steps"
        end if
    end do
    unit=7
    open(unit,file=dir//trim(root)//'_Euler_consistent_output.out',form="formatted",access="sequential",status="unknown")
    !> Output headers/blocks are written per step by write_out_step_vars (honours chem_out_options%vars)

!> -----------------------------------------------------------------------
!> Solver selection (only equilibrium-only is implemented for the non-ideal consistent solver)
    if (this%chem_syst%num_kin_reacts>0 .and. n_eq>0) then
        error stop "Non-ideal consistent solver: equilibrium + kinetic reactions not implemented yet"
    else if (n_eq>0) then
        p_solver=>mixing_iter_comp
    else
        error stop "Non-ideal consistent solver: kinetic-only reactions not implemented yet"
    end if
    !> Consistent reaction-rate mixing (r_tilde) method
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
!> Pre-compute canonical/non-canonical vectors and reaction-rate mixing ratios
    call mixing_ratios_conc%get_can_vec(this%CV_params%abs_tol,num_can_vec,ind_can_vec,num_non_can_vec,ind_non_can_vec)
    if (this%r_down_opt==4) then
        lambdas_R=mixing_ratios_R%compute_mix_ratios_R_opt4(mix_conc_indices)
    else
        allocate(lambdas_R(mixing_ratios_R%num_cols))
        do m=1,mixing_ratios_R%num_cols
            lambdas_R(m)=mixing_ratios_R%cols(m)%col_1(1)
        end do
    end if
    write(0,*) '[cons] Z2 lambdas_R done'; flush(0)

!> -----------------------------------------------------------------------
!> Initial conditions setup
    do i=1,this%num_waters
        call this%waters(i)%set_conc_old()
        call this%waters(i)%set_conc_old_old()
        call this%waters(i)%solid_chemistry%set_conc_old()
        call this%waters(i)%solid_chemistry%set_conc_old_old()
    end do
    write(0,*) '[cons] Z3 initial-conditions setup done'; flush(0)

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
    cpu_root = root//'_Euler_consistent'
    open(newunit=iunit_mb, file=dir//cpu_root//'_mass_balance.out', status='replace', action='write')
    write(iunit_mb,*) 'Domain-integrated reactive-mixing mass balance (per variable-activity species)'
    write(iunit_mb,*) 'Per time step & species: LHS = sum V*(c^{k+1}-c^k)/Dt, RHS = transport + reaction, residual = LHS - RHS'

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
    cpu_mix_total = 0d0
    cpu_r_tilde_total = 0d0
    cpu_solver_total = 0d0
    cpu_solid_gas_total = 0d0
    cpu_output_total = 0d0

!> Open CPU profiling output file (Eulerian reactive-mixing path; consistent tag)
    cpu_root = root//'_Euler_consistent'
    open(newunit=iunit_cpu, file=dir//cpu_root//'_cpu_profile.out', status='replace', action='write')
    write(0,*) '[cons] Z4 caches done, entering time loop'; flush(0)

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
        !> Advance state in place
        do i=1,this%num_target_waters
            call this%waters(this%tar_wat_indices(i))%update_attributes() !< Shift all old attributes
        end do
        write(0,*) '[cons] A: update loop done, k=',k; flush(0)
        !> Refresh the global conc_old cache after the update
        do w = 1, this%num_waters
            all_conc_old(1:size(this%waters(w)%conc_old), w) = this%waters(w)%conc_old
        end do
        !> Target waters loop (non-canonical columns only)
        do i=1,num_non_can_vec
            idx = this%tar_wat_indices(ind_non_can_vec(i))
            n_p = n_p_cache(i)
            n_v = n_v_cache(i)
            n_eq = n_eq_cache(i)
            if (mix_conc_indices%cols(ind_non_can_vec(i))%col_1(1)/=idx) then
                print *, mix_conc_indices%cols(ind_non_can_vec(i))%col_1(1), idx
                error stop "Target waters not in the right order"
            end if
            allocate(c_hat(n_v))
            allocate(r_tilde(n_v))
            num_mix_loc=mix_conc_indices%cols(ind_non_can_vec(i))%dim-3
            !> Conservative mixing (per-reactive-zone, non-ideal speciation); the reaction mixing
            !> ratios share the concentration stencil (0 for the external waters)
            call cpu_time(cpu_prof_t0)
            call p_mixer(this%waters,this%zone_waters(i,:),&
                mix_conc_indices%cols(ind_non_can_vec(i))%col_1(1:mix_conc_indices%cols(ind_non_can_vec(i))%dim-2),&
                mixing_ratios_conc%cols(ind_non_can_vec(i))%col_1(1:mix_conc_indices%cols(ind_non_can_vec(i))%dim-2),&
                mixing_ratios_R%cols(ind_non_can_vec(i))%col_1(1:mix_conc_indices%cols(ind_non_can_vec(i))%dim-2),&
                this%wat_indices_zones,this%ext_waters_indices,Delta_t,theta_r,c_mix,c_tilde,ideal=.false.)
            call cpu_time(cpu_prof_t1)
            cpu_mix_total = cpu_mix_total + (cpu_prof_t1 - cpu_prof_t0)
            write(0,*) '[cons] C: compute_c_mix_global done, i=',i; flush(0)
            !> Consistent reaction-rate mixing
            call cpu_time(cpu_prof_t0)
            call p_r_tilde(this%waters,&
                mix_conc_indices%cols(ind_non_can_vec(i))%col_1(&
                2:mix_conc_indices%cols(ind_non_can_vec(i))%dim-2),&
                this%waters(idx)%solid_chemistry%reactive_zone%ind_var_act_species,&
                mixing_ratios_R%cols(ind_non_can_vec(i))%col_1(2:),&
                mix_conc_indices%cols(ind_non_can_vec(i))%col_1(&
                mix_conc_indices%cols(ind_non_can_vec(i))%dim-1),&
                mix_conc_indices%cols(ind_non_can_vec(i))%col_1(&
                mix_conc_indices%cols(ind_non_can_vec(i))%dim),&
                theta_r,r_tilde,&
                curr_water_idx=idx,&
                lambda_self=mixing_ratios_R%cols(ind_non_can_vec(i))%col_1(1),&
                wat_indices_zones=this%wat_indices_zones)
            call cpu_time(cpu_prof_t1)
            cpu_r_tilde_total = cpu_r_tilde_total + (cpu_prof_t1 - cpu_prof_t0)
            write(0,*) '[cons] D: p_r_tilde done, i=',i; flush(0)
            !> Positivity-preserving, component-conservative reactive update
            call compute_c_hat_patankar(c_mix(1:n_v),Delta_t,r_tilde,c_hat)
            call this%waters(idx)%set_conc_hat(c_hat) !> aqueous blocks of the Newton speciation input, for the conc_hat attribute
            !> Non-ideal equilibrium speciation (also finalises pH, salinity, ionic strength
            !> and non-ideal activity coefficients inside mixing_iter_comp)
            conc_nc=this%waters(idx)%get_conc_var_act()
            call cpu_time(cpu_prof_t0)
            call p_solver(this%waters(idx),this%waters(idx)%get_conc_old_old_prim_new_rz(),c_hat,&
                mixing_ratios_R%cols(ind_non_can_vec(i))%col_1(1),Delta_t,theta_r,conc_nc,conc_comp(1:n_p))
            call cpu_time(cpu_prof_t1)
            cpu_solver_total = cpu_solver_total + (cpu_prof_t1 - cpu_prof_t0)
            write(0,*) '[cons] E: p_solver done, i=',i; flush(0)
            write(0,*) '[cons] F i=',i,' p_solver done'; flush(0)
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
                end if
                deallocate(mal_residual)
            end if
            !> Equilibrium reaction amounts
            if (n_eq>0) then
                call this%waters(idx)%compute_Re(c_hat(1:n_v),Delta_t,&
                    lambdas_R(ind_non_can_vec(i)),theta_r)
            end if
            !> Domain-integrated mass-balance contribution of this cell
            call this%waters(idx)%get_react_mix_mass_bal_terms(c_mix(1:n_v), c_hat(1:n_v), &
                mixing_ratios_R%cols(ind_non_can_vec(i))%col_1(1), Delta_t, theta_r, &
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
            end if
            call cpu_time(cpu_prof_t1)
            cpu_solid_gas_total = cpu_solid_gas_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Accumulate reaction amounts
            this%waters(idx)%cum_rxn_amnt_kin=this%waters(idx)%cum_rxn_amnt_kin+this%waters(idx)%rxn_amnt_kin
            this%waters(idx)%solid_chemistry%cum_rxn_amnt_kin=&
                this%waters(idx)%solid_chemistry%cum_rxn_amnt_kin+this%waters(idx)%solid_chemistry%rxn_amnt_kin
            deallocate(c_hat, r_tilde, c_mix, c_tilde, conc_nc)
        end do
        !> Domain-integrated mass balance for this time step (per variable-activity species)
        call this%write_react_mix_mass_bal_dom(idx_ref,iunit_mb,time,n_v_ref,mb_lhs,mb_transp,mb_react)
        !> Output writing honouring chem_out_options%vars (conc / react_rate / vol_fract)
        call cpu_time(cpu_prof_t0)
        if (k==this%chem_out_options%time_steps(kk)) then
            call write_out_step_vars(this, unit, time)
            if (kk<this%chem_out_options%num_time_steps) kk=kk+1
        end if
        call cpu_time(cpu_prof_t1)
        cpu_output_total = cpu_output_total + (cpu_prof_t1 - cpu_prof_t0)
    end do

!> -----------------------------------------------------------------------
!> CPU profiling report
    call cpu_time(cpu_total_t1)
    cpu_total_elapsed = cpu_total_t1 - cpu_total_t0
    write(iunit_cpu,*) '====================================================================================='
    write(iunit_cpu,*) '  CPU PROFILING REPORT: solve_reactive_mixing_cons (non-ideal consistent)'
    write(iunit_cpu,*) '====================================================================================='
    write(iunit_cpu,'(A,I10)')     '  Total time steps:       ', num_time_loc
    write(iunit_cpu,'(A,I10)')     '  Non-canonical waters:   ', num_non_can_vec
    write(iunit_cpu,*) '-------------------------------------------------------------------------------------'
    write(iunit_cpu,'(A,F14.4,A)') '  Conservative mixing:    ', cpu_mix_total, ' s'
    write(iunit_cpu,'(A,F14.4,A)') '  Reaction-rate mixing:   ', cpu_r_tilde_total, ' s'
    write(iunit_cpu,'(A,F14.4,A)') '  Reactive iteration:     ', cpu_solver_total, ' s'
    write(iunit_cpu,'(A,F14.4,A)') '  Solid/gas chemistry:    ', cpu_solid_gas_total, ' s'
    write(iunit_cpu,'(A,F14.4,A)') '  Output writing:         ', cpu_output_total, ' s'
    write(iunit_cpu,*) '-------------------------------------------------------------------------------------'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  TOTAL SUBROUTINE TIME:  ', cpu_total_elapsed, ' s  (', cpu_total_elapsed/60d0, ' min)'
    write(iunit_cpu,*) '====================================================================================='
    close(iunit_cpu)
    print *, 'CPU profiling written to: ', dir//cpu_root//'_cpu_profile.out'

!> -----------------------------------------------------------------------
!> Close output file and clean up (state already advanced in place on this%waters)
    close(unit)
    close(iunit_mb)
    deallocate(mb_lhs, mb_transp, mb_react)
    deallocate(conc_comp)
    deallocate(lambdas_R)
    deallocate(all_conc_old, all_ind_aq_sp)
    deallocate(n_p_cache, n_v_cache, n_eq_cache)
    deallocate(has_minerals_flag, has_gases_flag)
end subroutine solve_reactive_mixing_cons
