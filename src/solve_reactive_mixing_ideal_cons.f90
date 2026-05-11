!> @brief Solves reactive mixing problem assuming ideal conditions with conservative formulation
!> @details Computes concentrations, activities, activity coefficients, reaction rates
!> and volumetric fractions of minerals (if present) using operator splitting.
subroutine solve_reactive_mixing_ideal_cons(this,dir,root,mixing_ratios_conc,mixing_ratios_R,&
    mix_conc_indices,mix_react_indices,time_discr,theta_r)
    use chemistry_m, only: chemistry_c
    use aqueous_chemistry_m, only: aqueous_chemistry_c, &
        reactive_mixing_iter_EI_eq_anal_ideal, &
        compute_r_tilde_impl_opt1, compute_c_mix_global, &
        reactive_mixing_iter_EI_kin_anal_ideal_opt2, &
        compute_r_tilde_impl_opt2, compute_r_tilde_impl_opt3, &
        compute_r_tilde_impl_opt4, reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2
    use arrays_m, only: real_array_c, int_array_c, copy_real_array
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
    integer(kind=4) :: i, j, k, kk, ii, m
    integer(kind=4) :: idx          !< Cached target water index: tar_wat_indices(ind_non_can_vec(i))
    integer(kind=4) :: j_mix        !< Loop index for mixing waters
!> Problem dimensions
    integer(kind=4) :: n_p, n_v, n_eq
    integer(kind=4) :: num_can_vec, num_non_can_vec
    integer(kind=4) :: num_mix_loc
    integer(kind=4) :: unit
    integer(kind=4) :: iunit_cpu
    integer(kind=4) :: num_time_loc
    integer(kind=4) :: kk_snap      !< Counter for chemistry snapshot output time steps
    integer(kind=4) :: num_lump
!> Index arrays
    integer(kind=4), allocatable :: ind_can_vec(:)
    integer(kind=4), allocatable :: ind_non_can_vec(:)
!> Time variables
    real(kind=8) :: time, Delta_t
!> Concentration and reaction rate arrays
    real(kind=8), allocatable :: c_hat(:)
    real(kind=8), allocatable :: r_tilde(:)
    real(kind=8), allocatable :: c_mix(:)
    real(kind=8), allocatable :: conc_nc(:)
    real(kind=8), allocatable :: mal_residual(:)
    real(kind=8) :: mal_max
    real(kind=8), allocatable :: conc_comp(:)
    real(kind=8), allocatable :: lambdas_R(:)
    type(real_array_c) :: mixing_ratios_R_modif
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
    real(kind=8) :: cpu_cache_total
    real(kind=8) :: cpu_getters_total
    real(kind=8) :: cpu_Rk_accum_total
    real(kind=8) :: cpu_modify_mix_total
!> Procedure pointers
    procedure(reactive_mixing_iter_EI_kin_anal_ideal_opt2), pointer :: p_solver=>null()
    procedure(compute_r_tilde_impl_opt1), pointer :: p_r_tilde=>null()

!> -----------------------------------------------------------------------
!> Initialization
    call cpu_time(cpu_total_t0)
    n_eq=0
    do i=1,this%num_target_waters
        !> Only count eq reactions from waters that also have kinetic reactions (true domain waters).
        !> Boundary-adjacent waters with eq-only minerals should not influence solver selection.
        if (this%waters(this%tar_wat_indices(i))%solid_chemistry%mineral_zone%num_minerals_kin>0) then
            n_eq=n_eq+this%waters(this%tar_wat_indices(i))%solid_chemistry%mineral_zone%num_minerals_eq
        end if
    end do
    call copy_real_array(mixing_ratios_R,mixing_ratios_R_modif)

    time=0d0
    kk=2
    ii=1
    kk_snap=2
    num_time_loc=time_discr%Num_time
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
    do j=1,this%waters(idx)%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
        if (this%chem_syst%species(this%waters(idx)%ind_var_act_species(j))%cst_act_flag .eqv. .false.) then
            write(unit,"(A15)", advance='no') this%chem_syst%species(this%waters(idx)%ind_var_act_species(j))%name
        end if
    end do
    write(unit,*)

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
    !> r_tilde computation method
    if (this%r_down_opt==1) then
        p_r_tilde=>compute_r_tilde_impl_opt1
    else if (this%r_down_opt==2) then
        p_r_tilde=>compute_r_tilde_impl_opt2
    else if (this%r_down_opt==3) then
        p_r_tilde=>compute_r_tilde_impl_opt3
    else if (this%r_down_opt==4) then
        p_r_tilde=>compute_r_tilde_impl_opt4
    else
        error stop "r down option not implemented yet"
    end if

!> -----------------------------------------------------------------------
!> Pre-compute canonical/non-canonical vectors and lambdas
    call mixing_ratios_conc%get_can_vec(this%CV_params%abs_tol,num_can_vec,ind_can_vec,num_non_can_vec,ind_non_can_vec)
    if (this%r_down_opt==4) then
        lambdas_R=mixing_ratios_R%compute_mix_ratios_R_opt4()
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
    cpu_cache_total = 0d0
    cpu_getters_total = 0d0
    cpu_Rk_accum_total = 0d0
    cpu_modify_mix_total = 0d0

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
        !> Update old attributes
        call cpu_time(cpu_prof_t0)
        do i=1,this%num_target_waters
            call this%waters(this%tar_wat_indices(i))%update_old_attributes()
        end do
        call cpu_time(cpu_prof_t1)
        cpu_update_old_total = cpu_update_old_total + (cpu_prof_t1 - cpu_prof_t0)
        !> Refresh global conc_old cache
        call cpu_time(cpu_prof_t0)
        do w = 1, this%num_waters
            all_conc_old(1:size(this%waters(w)%conc_old), w) = &
                this%waters(w)%conc_old
        end do
        call cpu_time(cpu_prof_t1)
        cpu_cache_total = cpu_cache_total + (cpu_prof_t1 - cpu_prof_t0)
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
            num_mix_loc=mix_conc_indices%cols(ind_non_can_vec(i))%dim-3
            call cpu_time(cpu_prof_t1)
            cpu_data_prep_total = cpu_data_prep_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Conservative mixing
            call cpu_time(cpu_prof_t0)
            call compute_c_mix_global(this%waters(mix_conc_indices%cols(ind_non_can_vec(i))%col_1(1)),&
                all_conc_old, all_ind_aq_sp, &
                mix_conc_indices%cols(ind_non_can_vec(i))%col_1(2:mix_conc_indices%cols(ind_non_can_vec(i))%dim-2), num_mix_loc, &
                mixing_ratios_conc%cols(ind_non_can_vec(i))%col_1, c_mix)
            call cpu_time(cpu_prof_t1)
            cpu_mix_total = cpu_mix_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Reaction rate computation
            call cpu_time(cpu_prof_t0)
            call p_r_tilde(this%waters(&
                mix_react_indices%cols(ind_non_can_vec(i))%col_1(&
                2:mix_react_indices%cols(ind_non_can_vec(i))%dim-2)),&
                this%waters(idx)%solid_chemistry%reactive_zone%ind_var_act_species,&
                this%waters(idx)%solid_chemistry%reactive_zone%speciation_alg%comp_mat,&
                mixing_ratios_R%cols(ind_non_can_vec(i))%col_1(2:),&
                mix_react_indices%cols(ind_non_can_vec(i))%col_1(&
                mix_react_indices%cols(ind_non_can_vec(i))%dim-1),&
                mix_react_indices%cols(ind_non_can_vec(i))%col_1(&
                mix_react_indices%cols(ind_non_can_vec(i))%dim),&
                theta_r,Delta_t,r_tilde)
            call cpu_time(cpu_prof_t1)
            cpu_r_tilde_total = cpu_r_tilde_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Lumping of reaction mixing ratios for stability
            call cpu_time(cpu_prof_t0)
            call this%waters(idx)%modify_mix_ratios_reacts(&
                mixing_ratios_R%cols(ind_non_can_vec(i))%col_1,c_mix(1:n_v),&
                Delta_t,r_tilde,&
                mixing_ratios_R_modif%cols(ind_non_can_vec(i))%col_1,num_lump)
            this%num_lump=this%num_lump+num_lump
            c_hat=c_mix(1:n_v)+Delta_t*r_tilde
            call cpu_time(cpu_prof_t1)
            cpu_modify_mix_total = cpu_modify_mix_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Getters (get_conc_nc)
            call cpu_time(cpu_prof_t0)
            conc_nc=this%waters(idx)%get_conc_nc()
            call cpu_time(cpu_prof_t1)
            cpu_getters_total = cpu_getters_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Reactive mixing solver
            call cpu_time(cpu_prof_t0)
            call p_solver(this%waters(idx),this%waters(idx)%get_c1_old_old(),c_hat,&
                mixing_ratios_R_modif%cols(ind_non_can_vec(i))%col_1(1),&
                lambdas_R(ind_non_can_vec(i)),Delta_t,theta_r,conc_nc,conc_comp(1:n_p))
            call cpu_time(cpu_prof_t1)
            cpu_solver_total = cpu_solver_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Mass action law check for equilibrium reactions (ideal solution)
            if (n_eq > 0) then
                mal_residual = matmul(this%waters(idx)%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star, &
                    log10(conc_nc(1:n_p))) + &
                    this%waters(idx)%solid_chemistry%reactive_zone%speciation_alg%logK_star - &
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
            call this%waters(idx)%compute_salinity()
            call this%waters(idx)%compute_ionic_strength()
            call cpu_time(cpu_prof_t1)
            cpu_derived_total = cpu_derived_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Equilibrium reaction amounts
            call cpu_time(cpu_prof_t0)
            if (n_eq>0 .and. this%waters(idx)%indices_rk%num_cols>0) then
                call this%waters(idx)%compute_Re_kin(c_hat(n_p+1:n_v),Delta_t,&
                    lambdas_R(ind_non_can_vec(i)))
            else if (n_eq>0) then
                call this%waters(idx)%compute_Re(c_hat(1:n_v),Delta_t,&
                    lambdas_R(ind_non_can_vec(i)))
            end if
            call cpu_time(cpu_prof_t1)
            cpu_eq_rates_total = cpu_eq_rates_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Solid/gas chemistry
            call cpu_time(cpu_prof_t0)
            if (has_minerals_flag(i)) then
                call this%waters(idx)%solid_chemistry%compute_mass_bal_mins(Delta_t)
                call this%waters(idx)%solid_chemistry%compute_conc_minerals_iter(Delta_t)
            end if
            if (has_gases_flag(i)) then
                call this%waters(idx)%gas_chemistry%compute_conc_gases_iter(Delta_t,&
                    this%waters(idx)%volume,&
                    [this%waters(idx)%re_mean,this%waters(idx)%rk_mean])
                call this%waters(idx)%gas_chemistry%compute_vol_gas_species_conc()
                call this%waters(idx)%gas_chemistry%compute_log_act_coeffs_gases()
            end if
            call cpu_time(cpu_prof_t1)
            cpu_solid_gas_total = cpu_solid_gas_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Accumulate reaction amounts
            call cpu_time(cpu_prof_t0)
            this%waters(idx)%Rk_accum=this%waters(idx)%Rk_accum+this%waters(idx)%Rk
            this%waters(idx)%solid_chemistry%Rk_accum=&
                this%waters(idx)%solid_chemistry%Rk_accum+this%waters(idx)%solid_chemistry%Rk
            call cpu_time(cpu_prof_t1)
            cpu_Rk_accum_total = cpu_Rk_accum_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Output writing
            call cpu_time(cpu_prof_t0)
            if ((k==this%chem_out_options%time_steps(kk)) .and. (idx==this%chem_out_options%ind_waters(ii))) then
                write(unit,"(I10)", advance='no') this%chem_out_options%ind_waters(ii)
                do j=1,n_v
                    if (this%chem_syst%species(this%waters(idx)%ind_var_act_species(j))%cst_act_flag .eqv. .false.) then
                        write(unit,"(ES15.5)", advance='no') conc_nc(j)
                    end if
                end do
                write(unit,*)
                if (ii<this%chem_out_options%num_waters) then
                    ii=ii+1
                else if (kk<this%chem_out_options%num_time_steps) then
                    kk=kk+1
                    ii=1
                else
                    call cpu_time(cpu_prof_t1)
                    cpu_output_total = cpu_output_total + (cpu_prof_t1 - cpu_prof_t0)
                    deallocate(c_hat, r_tilde, c_mix, conc_nc)
                    exit
                end if
            end if
            call cpu_time(cpu_prof_t1)
            cpu_output_total = cpu_output_total + (cpu_prof_t1 - cpu_prof_t0)
            deallocate(c_hat, r_tilde, c_mix, conc_nc)
        end do
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
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Cache refresh:          ', cpu_cache_total,      ' s  (', cpu_cache_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Getters (conc_nc):      ', cpu_getters_total,    ' s  (', cpu_getters_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Modify mix + c_hat:     ', cpu_modify_mix_total, ' s  (', cpu_modify_mix_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Rk accumulation:        ', cpu_Rk_accum_total,   ' s  (', cpu_Rk_accum_total/60d0, ' min)'
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
        write(iunit_cpu,'(A,F8.2,A)') '    Cache refresh:           ', 100d0*cpu_cache_total/cpu_total_elapsed,      ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Getters (conc_nc):       ', 100d0*cpu_getters_total/cpu_total_elapsed,    ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Modify mix + c_hat:      ', 100d0*cpu_modify_mix_total/cpu_total_elapsed, ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Rk accumulation:         ', 100d0*cpu_Rk_accum_total/cpu_total_elapsed,   ' %'
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
    write(iunit_cpu,'(A,F12.6,A)') '  Avg cache refresh per time step:             ', cpu_cache_total / dble(num_time_loc), ' s'
    write(iunit_cpu,'(A,F12.6,A)') '  Avg getters per time step:                   ', cpu_getters_total / dble(num_time_loc), ' s'
    write(iunit_cpu,'(A,F12.6,A)') '  Avg modify_mix per time step:                ', cpu_modify_mix_total / dble(num_time_loc), ' s'
    write(iunit_cpu,'(A,F12.6,A)') '  Avg Rk accum per time step:                  ', cpu_Rk_accum_total / dble(num_time_loc), ' s'
    close(iunit_cpu)
    print *, 'CPU profiling written to: ', dir//root//'_cpu_profile.out'

!> Close output file and clean up
    close(unit)
    deallocate(conc_comp)
    deallocate(lambdas_R)
    deallocate(all_conc_old, all_ind_aq_sp)
    deallocate(n_p_cache, n_v_cache, n_eq_cache)
    deallocate(has_minerals_flag, has_gases_flag)
    write(*,*) "Number of lumpings: ", this%num_lump
end subroutine solve_reactive_mixing_ideal_cons
