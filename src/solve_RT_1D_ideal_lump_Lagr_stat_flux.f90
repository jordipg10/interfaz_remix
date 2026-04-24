!> \file solve_RT_1D_ideal_lump_Lagr_stat_flux.f90
!> \brief Solves 1D reactive transport with ideal speciation, lumping, and Lagrangian particle tracking
!> \details
!>   This subroutine solves the 1D reactive transport problem using:
!>   - **Ideal speciation**: Assumes ideal activity coefficients (γ = 1)
!>   - **Lumping**: Applies lumping technique to mixing ratios of reaction amounts for efficiency
!>   - **Lagrangian approach**: Tracks fluid parcels (particles) as they move through domain
!>   - **Stationary flux**: Flux field does not change with time
!>
!>   **Purpose:**
!>   Couple advective-dispersive transport with geochemical reactions by operator splitting:
!>   1. **Transport step**: Move particles and mix using mixing ratios
!>   2. **Reaction step**: Solve chemical equilibrium/kinetics at each location
!>
!>   **Computed Quantities:**
!>   - Aqueous species concentrations
!>   - Activities and activity coefficients
!>   - Equilibrium reaction rates
!>   - Kinetic reaction rates
!>   - Mineral volumetric fractions (if minerals present)
!>   - Gas concentrations and volumes (if gases present)
!>
!>   **Key Features:**
!>   - **Lumping**: Groups similar reactive components to reduce computational cost
!>   - **Lagrangian tracking**: Particles carry chemical composition through domain
!>   - **Ideal conditions**: Simplifies activity calculations (aᵢ = cᵢ)
!>   - **Operator splitting**: Decouples transport and reaction for efficiency
!>
!>   **Time Integration Methods:**
!>   - Euler Explicit (EE): θ = 0, explicit reactions
!>   - Euler Implicit (EI): θ = 1, implicit reactions
!>   - Crank-Nicolson (CN): θ = 0.5, semi-implicit
!>
!>   **Chemical System Types:**
!>   1. Equilibrium only: Fast reactions, algebraic equations
!>   2. Kinetic only: Slow reactions, ODEs
!>   3. Mixed: Both equilibrium and kinetic reactions
!>
!> \param[inout] this Reactive transport object (transient or stationary)
!> \param[in]    dir  Directory path for output files
!> \param[in]    root Root name for output files (base filename without extension)
!> \author Jordi Petchamé-Guerrero
!> \date November 2025
!> \ingroup reactive transport
subroutine solve_RT_ideal_lump_Lagr_stat_flux_1D(this,dir,root)
    use RT_m, only: RT_1D_transient_c, RT_c, move_particles_stat_flux_1D, &
        move_particles_stat_flux_EC_1D !< Import 1D reactive transport classes
    use aqueous_chemistry_m, only: aqueous_chemistry_c, & !< Import aqueous chemistry class
        reactive_mixing_iter_EE_eq_kin_ideal, reactive_mixing_iter_EE_kin, & !< Euler explicit solvers with lumping
        reactive_mixing_iter_EI_kin_anal_ideal_opt2, mixing_iter_comp_ideal, & !< Implicit solvers for kinetic-only and equilibrium-only systems
        compute_r_tilde_impl_opt1, compute_r_tilde_impl_opt2, compute_c_mix, compute_c_mix_global, &
        compute_r_tilde_impl_opt3, compute_r_tilde_impl_opt4, reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2, &
        reactive_mixing_iter_EI_eq_anal_ideal !< Utility functions for mixing calculations
    implicit none !< Enforce explicit variable declarations for type safety
!> \name Arguments
!> \{
    class(RT_1D_transient_c) :: this !< Reactive transport 1D object (polymorphic: transient or stationary) [-]
    character(len=*), intent(in) :: dir !< Directory for output file [string]
    character(len=*), intent(in) :: root !< Root name for output file (base filename without extension) [string]
!> \}
!> \name Loop Counters and Indices
!> \{
    integer(kind=4) :: i !< Counter for target waters (domain cells) [-]
    integer(kind=4) :: m !< Counter for target waters (domain cells) [-]
    integer(kind=4) :: j !< Counter for target solids or species indices [-]
    integer(kind=4) :: l !< Counter for reactive zones (spatial geochemical regions) [-]
    integer(kind=4) :: k !< Main time step counter for transport loop (1 to Num_time) [-]
    integer(kind=4) :: kk !< Output time step counter for chemistry output options [-]
    integer(kind=4) :: ii !< Output target water counter for chemistry output options [-]
!> \}
!> \name Problem Dimensions
!> \{
    integer(kind=4) :: num_tar_wat !< Total number of target waters (domain cells + boundary cells) [-]
    integer(kind=4) :: num_tar_sol !< Total number of target solids (mineral assemblages) [-]
    integer(kind=4) :: n_p !< Number of primary species (basis species for speciation) [-]
    integer(kind=4) :: n_v !< Total number of variable activity species (aqueous + solid) [-]
    integer(kind=4) :: n_v_aq !< Number of aqueous variable activity species (ions, complexes) [-]
    integer(kind=4) :: n_v_aq_2 !< Number of aqueous secondary variable activity species (from equilibrium reactions) [-]
    integer(kind=4) :: n_eq !< Number of equilibrium reactions in the domain [-]
    integer(kind=4) :: mix_ind !< Starting index for mixing waters (1 for explicit, 2 for implicit methods) [-]
    integer(kind=4) :: num_can_vec !< Number of canonical vectors (basis vectors in mixing ratios matrix) [-]
    integer(kind=4) :: num_non_can_vec !< Number of non-canonical vectors (dependent vectors requiring reactive transport) [-]
    integer(kind=4) :: unit !< Fortran I/O unit number for output file [-]
    integer(kind=4) :: k_flux !< Counter for flux field changes (if flux is time-dependent) [-]
    integer(kind=4) :: idx !< Cached index ind_non_can_vec(i) for current target water
    integer(kind=4) :: kk_snap !< Counter for chemistry snapshot output time steps
!> \}
!> \name Index Arrays
!> \{
    integer(kind=4), allocatable :: tar_gas_indices(:) !< Indices of target gases in each reactive zone (maps gases to zones) [-]
    integer(kind=4), allocatable :: tar_sol_indices(:) !< Indices of target solids in each reactive zone (maps solids to zones) [-]
    integer(kind=4), allocatable :: tar_wat_indices(:) !< Indices of target waters in each reactive zone (maps waters to zones) [-]
    integer(kind=4), allocatable :: perm(:) !< Permutation vector for reordering aqueous concentrations [-]
    integer(kind=4), allocatable :: ind_can_vec(:) !< Indices of canonical vectors in mixing ratios matrix (basis vectors) [-]
    integer(kind=4), allocatable :: ind_non_can_vec(:) !< Indices of non-canonical vectors requiring reactive transport solve [-]
!> \}
!> \name Time and Mixing Parameters
!> \{
    REAL(KIND=8) :: time !< Current simulation time [s] (accumulated from t=0)
    REAL(KIND=8) :: Delta_t !< Time step size [s] (can vary between steps)
    REAL(KIND=8) :: theta !< Time weighting factor [-]: 0=explicit, 1=implicit, 0.5=Crank-Nicolson
    REAL(KIND=8) :: cpu_t0, cpu_t1 !< CPU time markers for time step profiling [s]
    REAL(KIND=8) :: cpu_mix_t0, cpu_mix_t1, cpu_mix_total !< CPU time markers for compute_c_mix profiling [s]
    REAL(KIND=8) :: cpu_total_t0, cpu_total_t1 !< CPU time markers for total subroutine profiling [s]
    REAL(KIND=8) :: cpu_prof_t0, cpu_prof_t1 !< Reusable CPU time markers for phase profiling [s]
    REAL(KIND=8) :: cpu_update_old_total !< Accumulated CPU time for update_old_attributes [s]
    REAL(KIND=8) :: cpu_displace_total !< Accumulated CPU time for particle displacement [s]
    REAL(KIND=8) :: cpu_solver_total !< Accumulated CPU time for reactive mixing solver [s]
    REAL(KIND=8) :: cpu_derived_total !< Accumulated CPU time for derived quantities (pH, salinity, I) [s]
    REAL(KIND=8) :: cpu_eq_rates_total !< Accumulated CPU time for equilibrium reaction rates [s]
    REAL(KIND=8) :: cpu_solid_gas_total !< Accumulated CPU time for solid/gas chemistry [s]
    REAL(KIND=8) :: cpu_output_total !< Accumulated CPU time for output writing [s]
    REAL(KIND=8) :: cpu_total_elapsed !< Total subroutine elapsed CPU time [s]
    REAL(KIND=8) :: cpu_data_prep_total !< Accumulated CPU time for data preparation (dims, realloc, copy) [s]
    REAL(KIND=8) :: cpu_getters_total !< Accumulated CPU time for get_conc_nc + get_c1_old_old [s]
    REAL(KIND=8) :: cpu_snapshot_total !< Accumulated CPU time for write_chemistry snapshots [s]
    integer(kind=4) :: iunit_cpu !< Fortran I/O unit number for CPU profiling output file [-]
    REAL(KIND=8) :: y !< Sum of upstream mixing ratios [-] (for mass balance validation)
    real(kind=8) :: lambda_R !< Modified mixing ratio for reaction rates [-]
    real(kind=8) :: lambda_upstream !< Lumped lambda of upstream water (1.0 for boundary) [-]
!> \}
!> \name Concentration and Reaction Rate Arrays
!> \{
    REAL(KIND=8), allocatable :: conc_comp(:) !< Variable activity mobile species concentrations after transport mixing [mol/L]
    real(kind=8), allocatable :: lambdas_R(:) !< Mixed kinetic reaction rates for current cell [mol/L/s]
    REAL(KIND=8), allocatable :: c_hat(:) !< Concentrations from previous time step (before mixing) [mol/L]
    REAL(KIND=8), allocatable :: c_mix(:) !< Concentrations from previous time step (before mixing) [mol/L]
    REAL(KIND=8), allocatable :: conc_nc(:) !< Variable activity species concentrations (solution vector) [mol/L]
    REAL(KIND=8), allocatable :: r_tilde(:) !< Variable activity species concentrations (solution vector) [mol/L]
    REAL(KIND=8), allocatable :: lumped_lambdas(:) !< Variable activity species concentrations (solution vector) [mol/L]
!> \}
!> \name Chemistry State Objects
!> \{
    type(aqueous_chemistry_c), allocatable :: waters_new(:) !< Target waters at time step k+1 (after reaction step) [-]
    type(aqueous_chemistry_c), allocatable :: waters_old(:) !< Target waters at time step k (current time) [-]
    type(aqueous_chemistry_c), allocatable :: waters_old_old(:) !< Target waters at time step k-1 (for multistep methods) [-]
    type(aqueous_chemistry_c), allocatable :: mix_waters(:) !< Array of mixing waters (upstream contributors) [-]
    integer(kind=4) :: num_lump !< Number of lumping operations performed (for efficiency tracking) [-]
    integer(kind=4) :: j_mix !< loop index for mixing waters
    integer(kind=4) :: num_mix_loc !< number of mixing waters
    real(kind=8), allocatable :: conc_old_mix(:,:) !< conc_old from each mixing water
    integer(kind=4), allocatable :: ind_aq_sp_mix(:,:) !< indices_aq_species from each mixing water
!> \}
!> \name Pre-computed caches (invariant across time steps)
!> \{
    real(kind=8), allocatable :: all_conc_old(:,:) !< Cached conc_old for ALL waters, refreshed each time step
    integer(kind=4), allocatable :: all_ind_aq_sp(:,:) !< Cached indices_aq_species for ALL waters: all_ind_aq_sp(k,w) = position in water w's conc_old of aq_phase species k
    integer(kind=4) :: max_n_conc, max_n_idx, max_num_mix_1d !< Max sizes for pre-allocation
    integer(kind=4) :: max_n_p, max_n_v !< Max primary species and var act species across all target waters
    integer(kind=4) :: num_waters_1d, w !< Cached num_waters and loop index
    real(kind=8), allocatable :: c1_old_old_buf(:) !< Pre-allocated buffer for c1_old_old (avoids per-iteration allocation)
    real(kind=8), allocatable :: c1_downstream_buf(:) !< Pre-allocated buffer for downstream initial guess (current water's primary species)
    logical :: is_output_step !< Flag: is current time step an output step?
    logical :: has_cat_exch !< Flag: does current target water have cation exchange?
    logical :: has_gas !< Flag: does current target water have gas chemistry?

!> \}
!> \name Procedure Pointers (Runtime Polymorphism)
!> \{
    !> Procedure pointer for reactive mixing solver (selected at runtime based on chemical system)
    procedure(reactive_mixing_iter_EI_kin_anal_ideal_opt2), pointer :: p_solver=>null() !< Pointer to reactive mixing iteration solver
    !> Procedure pointer for kinetic rate mixing (selected based on integration method)
    procedure(compute_r_tilde_impl_opt1), pointer :: p_r_tilde=>null() !< Pointer to kinetic rate mixing function
    !> Procedure pointer for particle displacement under stationary flux (selected based on scheme)
    procedure(move_particles_stat_flux_1D), pointer :: p_displace_particles=>null() !< Pointer to particle displacement function (selected based on flux type and transport method)
!> \}
!> \section init_section Initialization
!> Initialize target waters for time stepping
    !allocate(waters_old(this%chemistry%num_waters)) !< Allocate array for old target waters
    !allocate(waters_old_old(this%chemistry%num_waters)) !< Allocate array for old old target waters
    call cpu_time(cpu_total_t0) !< Start total subroutine CPU timer
    n_eq=0 !< Initialize total number of equilibrium reactions
    do i=1,this%chemistry%num_target_waters !< Loop over all target waters
        n_eq=n_eq+this%target_waters(&
            i)%aq_chem%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions !< Count total equilibrium reactions
        !call waters_old(i)%copy_aq_chem(this%chemistry%waters(i)) !< Initialize time stepping variables
        !call waters_old_old(i)%copy_aq_chem(waters_old(i)) !< Initialize time stepping variables
    end do
    !waters_old=waters_old !< Copy current state to "old" (time level k)
    !waters_old_old=waters_old !< Copy to "old old" (time level k-1) for multistep schemes

!> \section type_select Type Selection and Setup
!> Select transient reactive transport type (polymorphic dispatch)
!select type (this)
!type is (RT_1D_transient_c) !< Transient (time-dependent) reactive transport
    !call this%get_transport_obj() !< Initialize transport object (e.g., set up time discretization)
    time=0d0 !< Initialize simulation time to zero [s]
    k_flux=1 !< Initialize flux change counter (for time-varying velocity fields) [-]
    kk=2 !< Initialize output time step counter (start at 2, as 1 is initial condition) [-]
    ii=1 !< Initialize output target water counter (start at first output location) [-]
    !< Set final output time step to total number of time steps (workaround for output timing)
    this%chemistry%chem_out_options%time_steps(this%chemistry%chem_out_options%num_time_steps)=this%transport%time_discr%Num_time
    !< Validate that user-specified output time steps do not exceed total number of time steps
    do i=2,this%chemistry%chem_out_options%num_time_steps-1
        if (this%chemistry%chem_out_options%time_steps(i) > this%transport%time_discr%Num_time) then
            print *, "Output time step index", this%chemistry%chem_out_options%time_steps(i), &
                "exceeds total number of time steps", this%transport%time_discr%Num_time
            error stop "Mismatch between out_opts.dat and discr_temp.dat: output time step exceeds total time steps"
        end if
    end do
    kk_snap=2 !< Initialize snapshot output counter (skip sentinel at index 1)
    unit=7 !< Assign Fortran I/O unit number for output file (arbitrary choice) [-]
    !< Open output file for formatted sequential writing
    open(unit,file=dir//root//'_output.out',form="formatted",access="sequential",status="unknown")
!> \section solver_select Solver Selection
!> Select reactive mixing subroutine based on:
!> - Chemical system type (equilibrium/kinetic/both)
!> - Time integration method (explicit/implicit)
!> - Jacobian calculation method (analytical/numerical)
!> \subsection eq_kin_case Equilibrium + Kinetic Reactions
    if (this%chemistry%chem_syst%num_kin_reacts>0 .and. n_eq>0) then !< Both equilibrium and kinetic reactions present
            if (this%chemistry%rk_avg_opt==1) then !< Option 1: Simple averaging of kinetic rates
                error stop "rk average option 1 not implemented yet" !< Not implemented
            else if (this%chemistry%rk_avg_opt==2) then !< Option 2: Improved kinetic rate averaging with lumping
                p_solver=>reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2 !< Assign implicit solver with lumping
            else
                error stop "rk average option not implemented yet" !< Unsupported option
            end if
!> \subsection eq_only_case Equilibrium-Only Reactions
    else if (n_eq>0) then !< Only equilibrium reactions (no kinetics)
        p_solver=>reactive_mixing_iter_EI_eq_anal_ideal !< Assign implicit solver (only implemented option for eq-only)
!> \subsection kin_only_case Kinetic-Only Reactions
    else !< Only kinetic reactions (no equilibrium)
            if (this%chemistry%rk_avg_opt==1) then !< Option 1: Concentration averaging
                error stop "rk average option 1 not implemented yet" !< Not implemented
            else if (this%chemistry%rk_avg_opt==2) then !< Option 2: Reaction rate averaging
                p_solver=>reactive_mixing_iter_EI_kin_anal_ideal_opt2 !< Assign implicit kinetic solver with lumping
            else
                error stop "rk average option not implemented yet" !< Unsupported option
            end if
    end if
    if (this%transport%spatial_discr%targets_flag==0) then !< Stationary flux scheme 1: Standard particle tracking
        p_displace_particles=>move_particles_stat_flux_1D !< Assign standard particle displacement function
    else if (this%transport%spatial_discr%targets_flag==1) then !< Stationary flux scheme 2: Particle tracking with Eulerian correction
        p_displace_particles=>move_particles_stat_flux_EC_1D !< Assign particle displacement function with Eulerian correction
    else
        error stop "Stationary flux Lagrangian scheme not implemented yet" !< Unsupported flux scheme
    end if
    allocate(lumped_lambdas(this%transport%mixing_ratios_R%num_cols)) !< Allocate array for lumped kinetic rates
    do i=1,this%transport%mixing_ratios_R%num_cols !< Loop over mixing ratios columns
        lumped_lambdas(i)=sum(this%transport%mixing_ratios_R%cols(i)%col_1) !< Sum mixing ratios for mass balance check
        !if (abs(lumped_lambdas(i)-1d0)>this%chemistry%CV_params%abs_tol) then !< Validate sum to unity within tolerance
        !    print *, "Mixing ratios do not sum to 1 in column ", i, " sum = ", lumped_lambdas(i)
        !    error stop "Mixing ratios validation failed"
        !end if
    end do
    !lumped_lambdas=>null() !< Initialize lumped kinetic rates pointer to null
    !> \section canonical_vectors Canonical Vector Extraction
    !> Extract canonical vectors (basis vectors) from mixing ratios matrix
    !> Canonical vectors represent independent mixing patterns that don't require reactive transport
    call this%transport%mixing_ratios_conc%get_can_vec(this%chemistry%CV_params%abs_tol,num_can_vec,ind_can_vec,num_non_can_vec,&
        ind_non_can_vec) !< Get indices of canonical and non-canonical vectors using tolerance
    !> \section init_conc Initialize Concentration History
    !> Set old concentrations in all target waters for time stepping
    do i=1,this%chemistry%num_waters !< Loop over all waters (domain + boundary)
        !print *, this%chemistry%waters(i)%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
        !call this%chemistry%waters_init(i)%set_conc_old() !< Copy current to old (time level k to k-1)
        call this%chemistry%waters(i)%set_conc_old() !< Copy current to old (time level k to k-1)
        !call this%chemistry%waters_init(i)%set_conc_old_old() !< Copy old to old_old (k-1 to k-2)
        call this%chemistry%waters(i)%set_conc_old_old() !< Copy old to old_old (k-1 to k-2)
        !call this%chemistry%waters_init(i)%solid_chemistry%set_conc_old() !< Copy solid concentrations to old
        call this%chemistry%waters(i)%solid_chemistry%set_conc_old() !< Copy solid concentrations to old
        !call this%chemistry%waters_init(i)%solid_chemistry%set_conc_old_old() !< Copy solid concentrations to old_old
        call this%chemistry%waters(i)%set_solid_chemistry_old() !< Copy current kinetic rates to old
    end do
    !allocate(mix_waters(this%chemistry%num_waters))
    !> \section output_header Write Output File Header
    !> Write species names as column headers in output file
    write(unit,"(A10,*(A15))") 'Water', &
        (this%chemistry%chem_syst%aq_phase%aq_species(&
        this%chemistry%chem_out_options%ind_aq_species(j))%name, & !< Species name for column j
        j=1,this%chemistry%chem_out_options%num_aq_species) !< Loop over species to output
    !> \section precompute_1d Pre-compute global caches for 1D solver
    num_waters_1d = this%chemistry%num_waters
    max_n_conc = 0; max_n_idx = 0; max_num_mix_1d = 0; max_n_p = 0; max_n_v = 0
    do i = 1, num_non_can_vec
        idx = ind_non_can_vec(i)
        associate(mix_col => this%transport%mix_conc_indices%cols(idx)%col_1, &
                  spec_alg => this%target_waters(idx)%aq_chem%solid_chemistry%reactive_zone%speciation_alg)
        j = size(this%chemistry%waters(mix_col(1))%conc_old)
        if (j > max_n_conc) max_n_conc = j
        j = this%transport%mix_conc_indices%cols(idx)%dim - 3
        if (j > max_num_mix_1d) max_num_mix_1d = j
        if (spec_alg%num_prim_species > max_n_p) max_n_p = spec_alg%num_prim_species
        if (spec_alg%num_var_act_species > max_n_v) max_n_v = spec_alg%num_var_act_species
        end associate
    end do
    allocate(all_conc_old(max_n_conc, num_waters_1d))
    !> Cache indices_aq_species for ALL waters in contiguous 2D array
    !> all_ind_aq_sp(k, w) = indices_aq_species(k) of water w = position in water w's conc_old of aq_phase species k
    max_n_idx = 0
    do w = 1, num_waters_1d
        j = this%chemistry%waters(w)%aq_phase%num_species
        if (j > max_n_idx) max_n_idx = j
    end do
    allocate(all_ind_aq_sp(max_n_idx, num_waters_1d))
    do w = 1, num_waters_1d
        all_ind_aq_sp(1:size(this%chemistry%waters(w)%indices_aq_species), w) = &
            this%chemistry%waters(w)%indices_aq_species
    end do
    !> Pre-allocate reusable buffers (avoids alloc/dealloc per inner iteration)
    allocate(conc_comp(max_n_p))
    allocate(c1_old_old_buf(max_n_p))
    allocate(c1_downstream_buf(max_n_p))
    cpu_mix_total = 0d0
    cpu_update_old_total = 0d0
    cpu_displace_total = 0d0
    cpu_solver_total = 0d0
    cpu_derived_total = 0d0
    cpu_eq_rates_total = 0d0
    cpu_solid_gas_total = 0d0
    cpu_output_total = 0d0
    cpu_data_prep_total = 0d0
    cpu_getters_total = 0d0
    cpu_snapshot_total = 0d0
    open(newunit=iunit_cpu, file=dir//trim(root)//'_cpu_profile.out', status='replace', action='write')
        !> \section time_loop Main Time Loop
        !> Loop over all time steps to solve reactive transport
        do k=1,this%transport%time_discr%Num_time !< Loop from time step 1 to final time step
            call cpu_time(cpu_t0) !< Start CPU timer for this time step
            Delta_t=this%transport%time_discr%get_Delta_t(k) !< Get time step size for step k
            time=time+Delta_t !< Advance simulation time

            !> Write time to output file if at output time step
            is_output_step = (k==this%chemistry%chem_out_options%time_steps(kk))
            if (is_output_step) then
               write(unit,"(/,2x,'t = ',*(ES15.5),/)") time !< Write current time to file
            end if
            !> \subsection update_old Update Previous Time Step Values
            !> Shift time levels: current → old, old → old_old for multistep methods
            call cpu_time(cpu_prof_t0)
            do i=1,this%chemistry%num_waters !< Loop over domain target waters only (not boundaries)
                !call this%target_waters(i)%update_old_attributes() !< Shift all old attributes
                call this%chemistry%waters(i)%update_old_attributes() !< Shift all old attributes
                ! call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%update_conc_old() !< Shift aqueous concentrations
                ! call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%solid_chemistry%update_conc_old() !< Shift solid concentrations
                ! call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%update_rk_old() !< Shift aqueous kinetic rates
                ! call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%solid_chemistry%update_rk_old() !< Shift solid kinetic rates
                ! call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%update_solid_chemistry_old() !< Shift aqueous concentrations (old_old)
                ! call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%update_indices_rk_old() !< Shift solid concentrations (old_old)
            end do
            call cpu_time(cpu_prof_t1)
            cpu_update_old_total = cpu_update_old_total + (cpu_prof_t1 - cpu_prof_t0)
            !call this%move_particles_stat_flux(k) !< Displace particles based on stationary flux field
            call cpu_time(cpu_prof_t0)
            call p_displace_particles(this,k) !< Displace particles based on stationary flux field
            call this%introduce_particle(k) !< Add new particles entering domain at time step k
            call cpu_time(cpu_prof_t1)
            cpu_displace_total = cpu_displace_total + (cpu_prof_t1 - cpu_prof_t0)
            !> Refresh global conc_old cache AFTER displacement (post-displacement conc_old)
            do w = 1, num_waters_1d
                all_conc_old(1:size(this%chemistry%waters(w)%conc_old), w) = &
                    this%chemistry%waters(w)%conc_old
            end do
            !> \subsection target_loop Target Waters Loop (Spatial Domain)
            !> Loop over all domain target waters to solve reactive transport at each location
            do i=1,num_non_can_vec !< Loop over target waters (spatial cells)
                idx = ind_non_can_vec(i) !< Cache non-canonical vector index
                associate(tw => this%target_waters(idx)%aq_chem, &
                          spec_alg => this%target_waters(idx)%aq_chem%solid_chemistry%reactive_zone%speciation_alg)
                !> \subsubsection extract_dims Extract Problem Dimensions
                call cpu_time(cpu_prof_t0)
                n_p=spec_alg%num_prim_species !< Number of primary species [-]
                n_v=spec_alg%num_var_act_species !< Total variable activity species [-]
                n_v_aq=spec_alg%num_aq_var_act_species !< Aqueous variable activity species [-]
                n_v_aq_2=spec_alg%num_aq_sec_var_act_species !< Secondary aqueous variable activity species [-]
                !> Validate target water ordering
                associate(mix_col => this%transport%mix_conc_indices%cols(idx)%col_1)
                if (mix_col(1) /= this%chemistry%tar_wat_indices(idx)) then
                        print *, mix_col(1), this%chemistry%tar_wat_indices(idx)
                        print *, "Problematic target water index: ", i
                        error stop "Waters not in the right order in mixing indices"
                    end if
                !> \subsubsection transport_mix Transport Mixing Step
                num_mix_loc=this%transport%mix_conc_indices%cols(idx)%dim-3
                call cpu_time(cpu_prof_t1)
                cpu_data_prep_total = cpu_data_prep_total + (cpu_prof_t1 - cpu_prof_t0)
                call cpu_time(cpu_mix_t0)
                call compute_c_mix_global(this%chemistry%waters(mix_col(1)),&
                    all_conc_old, all_ind_aq_sp, mix_col(2:), num_mix_loc,&
                    this%transport%mixing_ratios_conc%cols(idx)%col_1,&
                    c_hat)
                call cpu_time(cpu_mix_t1)
                cpu_mix_total = cpu_mix_total + (cpu_mix_t1 - cpu_mix_t0)
                end associate !< End mix_col association
                !> \subsubsection fill_c1_old_old Fill c1_old_old buffer inline (avoids allocating temporary)
                call cpu_time(cpu_prof_t0)
                conc_nc=tw%get_conc_nc() !< Get variable activity species concentrations
                c1_old_old_buf(1:n_p) = 1d0 !< Initialize (matches get_c1_old_old default)
                do j=1,spec_alg%num_aq_prim_species
                    c1_old_old_buf(j) = tw%conc_old_old(tw%ind_prim_species(j))
                end do
                has_cat_exch = (tw%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl > 0)
                has_gas = associated(tw%gas_chemistry)
                if (has_cat_exch) then
                    c1_old_old_buf(spec_alg%num_aq_prim_species+1) = &
                        tw%solid_chemistry%conc_old_old(tw%solid_chemistry%mineral_zone%num_minerals+1)
                end if
                if (has_gas) then
                    c1_old_old_buf(spec_alg%num_aq_prim_species+1:n_p) = &
                        tw%gas_chemistry%conc_old_old( &
                        tw%gas_chemistry%reactive_zone%gas_phase%num_gases_eq+1: &
                        tw%gas_chemistry%reactive_zone%gas_phase%num_species) / tw%volume
                end if
                call cpu_time(cpu_prof_t1)
                cpu_getters_total = cpu_getters_total + (cpu_prof_t1 - cpu_prof_t0)
                !> Build downstream initial guess using CURRENT water's primary species indices
                block
                    integer :: j_ds, ds_wat_idx
                    ds_wat_idx = this%chemistry%downstream_water_indices(idx)
                    c1_downstream_buf(1:n_p) = 1d0
                    do j_ds = 1, spec_alg%num_aq_prim_species
                        c1_downstream_buf(j_ds) = this%chemistry%waters(ds_wat_idx)%concentrations(tw%ind_prim_species(j_ds))
                    end do
                    if (has_cat_exch) then
                        c1_downstream_buf(spec_alg%num_aq_prim_species+1) = &
                            this%chemistry%waters(ds_wat_idx)%solid_chemistry%concentrations( &
                            tw%solid_chemistry%mineral_zone%num_minerals+1)
                    end if
                    if (has_gas) then
                        c1_downstream_buf(spec_alg%num_aq_prim_species+1:n_p) = &
                            this%chemistry%waters(ds_wat_idx)%gas_chemistry%concentrations( &
                            tw%gas_chemistry%reactive_zone%gas_phase%num_gases_eq+1: &
                            tw%gas_chemistry%reactive_zone%gas_phase%num_species) / tw%volume
                    end if
                end block
                !> \subsubsection reactive_mix Reactive Mixing Iteration
                call cpu_time(cpu_prof_t0)
                if (this%chemistry%upstream_water_indices(idx)==0) then
                    lambda_upstream = 1d0 !< Boundary upstream: pure source
                else
                    lambda_upstream = lumped_lambdas(this%chemistry%upstream_water_indices(idx))
                end if
                call p_solver(tw,& !< Current target water
                    c1_old_old_buf(1:n_p),& !< Primary species from k-1 [mol/L] (pre-allocated buffer)
                    c_hat(1:n_v),& !< Mixed var act species from transport (trim cst act tail)
                    lambda_upstream,& !< Lumped lambda from upstream water [-]
                    lumped_lambdas(idx),& !< Reaction mixing ratio for current target water [-]
                    Delta_t,& !< Time step size
                    this%transport%time_discr%theta_r,& !< Time weighting factor [-]
                    conc_nc,&
                    conc_comp(1:n_p),& !< Output: new concentrations (uses slice of pre-allocated buffer)
                    c1_downstream_buf(1:n_p)) !< Downstream water primary concs (current water's species) as initial guess
                call cpu_time(cpu_prof_t1)
                cpu_solver_total = cpu_solver_total + (cpu_prof_t1 - cpu_prof_t0)
                
                !> \subsubsection update_derived Update Derived Quantities
                call cpu_time(cpu_prof_t0)
                call tw%compute_pH()
                call tw%compute_salinity()
                call tw%compute_ionic_strength()
                call cpu_time(cpu_prof_t1)
                cpu_derived_total = cpu_derived_total + (cpu_prof_t1 - cpu_prof_t0)
                
                !> \subsubsection accum_rk Accumulate Reaction Amounts
                tw%Rk_accum=tw%Rk_accum+tw%Rk !< Accumulate aqueous reaction amounts [mol]
                tw%solid_chemistry%Rk_accum=tw%solid_chemistry%Rk_accum+tw%solid_chemistry%Rk !< Accumulate solid reaction amounts [mol]
                
                !> \subsubsection eq_rates Equilibrium Reaction Rates
                call cpu_time(cpu_prof_t0)
                if (spec_alg%num_eq_reactions>0 .and. &
                    tw%indices_rk%num_cols>0) then !< Both equilibrium and kinetic reactions present
                    call tw%compute_Re_kin(c_hat(n_p+1:n_v), Delta_t, lumped_lambdas(idx))
                else if (spec_alg%num_eq_reactions>0) then !< Only equilibrium reactions present
                    call tw%compute_Re(c_hat(1:n_v), Delta_t, lumped_lambdas(idx))
                end if
                call cpu_time(cpu_prof_t1)
                cpu_eq_rates_total = cpu_eq_rates_total + (cpu_prof_t1 - cpu_prof_t0)
                
                !> \subsubsection solid_chem Solid Chemistry State Variables
                call cpu_time(cpu_prof_t0)
                if (associated(tw%solid_chemistry%mineral_zone)) then
                    call tw%solid_chemistry%compute_mass_bal_mins(Delta_t)
                    call tw%solid_chemistry%compute_conc_minerals_iter(Delta_t)
                end if
                
                !> \subsubsection gas_species_chem Gas Chemistry State Variables
                if (associated(tw%gas_chemistry)) then
                    call tw%gas_chemistry%compute_conc_gases_iter(&
                        Delta_t, tw%volume, [tw%re_mean, tw%rk_mean])
                    call tw%gas_chemistry%compute_vol_gas_species_conc()
                    call tw%gas_chemistry%compute_log_act_coeffs_gases()
                end if
                call cpu_time(cpu_prof_t1)
                cpu_solid_gas_total = cpu_solid_gas_total + (cpu_prof_t1 - cpu_prof_t0)
                !> \subsubsection write_output Write Chemistry Output
                if (is_output_step .and. (idx==this%chemistry%chem_out_options%ind_waters(ii))) then
                    call cpu_time(cpu_prof_t0)
                    write(unit,"(I10,*(ES15.5))") this%chemistry%chem_out_options%ind_waters(ii), &
                        (tw%concentrations(tw%indices_aq_species(&
                        this%chemistry%chem_out_options%ind_aq_species(j))), &
                        j=1,this%chemistry%chem_out_options%num_aq_species)
                    if (ii<this%chemistry%chem_out_options%num_waters) then
                        ii=ii+1
                    else if (kk<this%chemistry%chem_out_options%num_time_steps) then
                        kk=kk+1
                        ii=1
                    else
                        continue
                    end if
                    call cpu_time(cpu_prof_t1)
                    cpu_output_total = cpu_output_total + (cpu_prof_t1 - cpu_prof_t0)
                end if
                !> \subsubsection cleanup_dealloc Cleanup and Deallocation
                deallocate(c_hat)
                end associate !< End tw/spec_alg association
        end do !< End target waters spatial loop
        !> Close c_hat output file after first time step
        !> \subsection write_chem_snapshot Write chemistry snapshot at output time steps
        !> (Snapshots disabled — write_chemistry is called once at the end via write_RT)
        call cpu_time(cpu_t1) !< End CPU timer for this time step
        !write(*,'(A,I10,A,F14.4,A)') '[CPU] time step ', k, ' took ', cpu_t1 - cpu_t0, ' s'
        !> Write aqueous concentrations after the first displacement
        
        !this%target_waters(this%chemistry%num_target_waters)%aq_chem%concentrations=&
        !    this%chemistry%wat_types(this%chemistry%num_wat_types)%concentrations !< Update downstream boundary concentrations
        !this%target_waters(this%chemistry%num_target_waters)%aq_chem%activities=&
        !    this%chemistry%wat_types(this%chemistry%num_wat_types)%activities !< Update downstream boundary concentrations
        end do !< End main time loop
    
    call cpu_time(cpu_total_t1) !< Stop total subroutine CPU timer
    cpu_total_elapsed = cpu_total_t1 - cpu_total_t0
    write(iunit_cpu,*) '====================================================================================='
    write(iunit_cpu,*) '  CPU TIME PROFILING: solve_RT_ideal_lump_Lagr_stat_flux_1D'
    write(iunit_cpu,*) '====================================================================================='
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Update old attributes:  ', cpu_update_old_total, ' s  (', cpu_update_old_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Particle displacement:  ', cpu_displace_total,   ' s  (', cpu_displace_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Data preparation:       ', cpu_data_prep_total,  ' s  (', cpu_data_prep_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Conservative mixing (compute_c_mix):', cpu_mix_total,        ' s  (', cpu_mix_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Getters (conc_nc+c1_oo):', cpu_getters_total,    ' s  (', cpu_getters_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Reactive mixing (p_solver):', cpu_solver_total,   ' s  (', cpu_solver_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Derived quantities:     ', cpu_derived_total,    ' s  (', cpu_derived_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Equilibrium rates:      ', cpu_eq_rates_total,   ' s  (', cpu_eq_rates_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Solid/gas chemistry:    ', cpu_solid_gas_total,  ' s  (', cpu_solid_gas_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Output writing:         ', cpu_output_total,     ' s  (', cpu_output_total/60d0, ' min)'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  Snapshot (write_chem):  ', cpu_snapshot_total,   ' s  (', cpu_snapshot_total/60d0, ' min)'
    write(iunit_cpu,*) '-------------------------------------------------------------------------------------'
    write(iunit_cpu,'(A,F14.4,A,F10.2,A)') '  TOTAL SUBROUTINE TIME:  ', cpu_total_elapsed,    ' s  (', cpu_total_elapsed/60d0, ' min)'
    write(iunit_cpu,*) '====================================================================================='
    if (cpu_total_elapsed > 0.0d0) then
        write(iunit_cpu,*) ''
        write(iunit_cpu,*) '  Breakdown (% of total):'
        write(iunit_cpu,'(A,F8.2,A)') '    Update old attributes:   ', 100d0*cpu_update_old_total/cpu_total_elapsed, ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Particle displacement:   ', 100d0*cpu_displace_total/cpu_total_elapsed,   ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Data preparation:        ', 100d0*cpu_data_prep_total/cpu_total_elapsed,  ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Conservative mixing (compute_c_mix):', 100d0*cpu_mix_total/cpu_total_elapsed,        ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Getters (conc_nc+c1_oo): ', 100d0*cpu_getters_total/cpu_total_elapsed,    ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Reactive mixing (p_solver):         ', 100d0*cpu_solver_total/cpu_total_elapsed,     ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Derived quantities:      ', 100d0*cpu_derived_total/cpu_total_elapsed,    ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Equilibrium rates:       ', 100d0*cpu_eq_rates_total/cpu_total_elapsed,   ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Solid/gas chemistry:     ', 100d0*cpu_solid_gas_total/cpu_total_elapsed,  ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Output writing:          ', 100d0*cpu_output_total/cpu_total_elapsed,     ' %'
        write(iunit_cpu,'(A,F8.2,A)') '    Snapshot (write_chem):   ', 100d0*cpu_snapshot_total/cpu_total_elapsed,   ' %'
        write(iunit_cpu,*) '====================================================================================='
    end if
    write(iunit_cpu,*) ''
    write(iunit_cpu,'(A,F12.6,A)') '  Avg conservative mixing (compute_c_mix) per time step:        ', cpu_mix_total / dble(this%transport%time_discr%Num_time), ' s'
    write(iunit_cpu,'(A,F12.6,A)') '  Avg conservative mixing (compute_c_mix) per time step/target water:  ', cpu_mix_total / dble(this%transport%time_discr%Num_time * num_non_can_vec), ' s'
    write(iunit_cpu,'(A,F12.6,A)') '  Avg reactive mixing (p_solver) per time step:      ', cpu_solver_total / dble(this%transport%time_discr%Num_time), ' s'
    write(iunit_cpu,'(A,F12.6,A)') '  Avg reactive mixing (p_solver) per time step/target water:', cpu_solver_total / dble(this%transport%time_discr%Num_time * num_non_can_vec), ' s'
    close(iunit_cpu) !< Close CPU profiling output file
    print *, 'CPU profiling written to: ', dir//root//'_cpu_profile.out'
    !> Close output file
    close(unit)
    !> Deallocate pre-allocated buffers
    deallocate(conc_comp)
    deallocate(c1_old_old_buf)
    deallocate(c1_downstream_buf)
    if (allocated(conc_nc)) deallocate(conc_nc)
    !end select !< End type selection
end subroutine !< End solve_RT_ideal_lump_Lagr_stat_flux_1D