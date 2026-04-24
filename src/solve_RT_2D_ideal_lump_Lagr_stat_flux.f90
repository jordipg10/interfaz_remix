!> \file solve_RT_2D_ideal_lump_Lagr_stat_flux.f90
!> \brief Solves 2D reactive transport with ideal speciation, lumping, stationary flux, and Lagrangian particle tracking
!> \details
!>   This subroutine solves the 2D reactive transport problem using:
!>   - **Ideal speciation**: Assumes ideal activity coefficients (γ = 1)
!>   - **Lumping**: Applies lumping technique to mixing ratios of reaction amounts for efficiency
!>   - **Lagrangian approach**: Tracks water parcels (particles) as they move through domain
!>   - **Stationary flux**: Flux field does not change with time
!>
!>   **Purpose:**
!>   Decouple transport and reaction steps using WMA:
!>   1. **Transport step**: Move particles and mix using mixing ratios
!>   2. **Reaction step**: Solve chemical equilibrium/kinetics at each location
!>
!>   **Computed Quantities:**
!>   - Species and component concentrations
!>   - Activities and activity coefficients
!>   - Equilibrium reaction rates and amounts
!>   - Kinetic reaction rates and amounts
!>   - Mineral volumetric fractions (if minerals present)
!>   - Gas concentrations and volumes (if gases present)
!>
!>   **Key Features:**
!>   - **Lumping**: Groups similar reactive components to reduce computational cost
!>   - **Lagrangian tracking**: Particles carry chemical composition through domain
!>   - **Ideal conditions**: Simplifies activity calculations (aᵢ = cᵢ)
!>   - **WMA**: Decouples transport and reactions for efficiency
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
subroutine solve_RT_ideal_lump_Lagr_stat_flux_2D(this,dir,root)
    use RT_m, only: RT_2D_transient_c, move_particles_stat_flux_2D
    use aqueous_chemistry_m, only: aqueous_chemistry_c, &
        reactive_mixing_iter_EI_kin_anal_ideal_opt2, &
        compute_c_mix_global, &
        reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2, &
        reactive_mixing_iter_EI_eq_anal_ideal
    implicit none
!> Arguments
    class(RT_2D_transient_c) :: this
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: root
!> Loop counters and indices
    integer(kind=4) :: i
    integer(kind=4) :: j
    integer(kind=4) :: k                  !< Main time step counter
    integer(kind=4) :: kk                 !< Output time step counter
    integer(kind=4) :: ii                 !< Output target water counter
!> Problem dimensions
    integer(kind=4) :: n_p                !< Number of primary species
    integer(kind=4) :: n_v                !< Total number of variable activity species
    integer(kind=4) :: n_eq               !< Number of equilibrium reactions in the domain
    integer(kind=4) :: num_can_vec        !< Number of canonical vectors
    integer(kind=4) :: num_non_can_vec    !< Number of non-canonical vectors
    integer(kind=4) :: unit               !< Fortran I/O unit number for output file
    integer(kind=4) :: idx                !< Cached index ind_non_can_vec(i)
    integer(kind=4) :: kk_snap            !< Counter for chemistry snapshot output time steps
    integer(kind=4) :: num_eq_loc         !< Cached equilibrium reaction count for current target water
    integer(kind=4) :: prev_n_v           !< Previous n_v for conditional conc_nc reallocation
    integer(kind=4) :: num_waters_loc     !< Cached this%chemistry%num_waters
    integer(kind=4) :: num_time_loc       !< Cached total number of time steps
    integer(kind=4) :: num_out_ts         !< Cached chem_out_options%num_time_steps
    integer(kind=4) :: num_out_waters     !< Cached chem_out_options%num_waters
    integer(kind=4) :: num_out_aq_sp      !< Cached chem_out_options%num_aq_species
    integer(kind=4) :: num_kin_rk_cols    !< Cached tw%indices_rk%num_cols
    integer(kind=4) :: num_mix_loc        !< Number of mixing waters for current target
    logical :: is_output_step             !< True if current time step is an output step
    logical :: do_profile_inner           !< True if inner-loop profiling is active this time step
    integer(kind=4) :: num_profiled_steps !< Count of time steps with inner-loop profiling active
    real(kind=8) :: lumped_lambda_loc     !< Cached lumped_lambdas(idx) for current target water
    real(kind=8) :: theta_r_loc           !< Cached this%chemistry%theta_r
!> Index arrays
    integer(kind=4), allocatable :: ind_can_vec(:)     !< Indices of canonical vectors
    integer(kind=4), allocatable :: ind_non_can_vec(:) !< Indices of non-canonical vectors
    logical, allocatable :: has_minerals_flag(:)        !< Pre-computed mineral presence flag per target water
    logical, allocatable :: has_gases_flag(:)           !< Pre-computed gas phase presence flag per target water
!> Time and profiling variables
    real(kind=8) :: time                   !< Current simulation time [s]
    real(kind=8) :: Delta_t                !< Time step size [s]
    real(kind=8) :: cpu_t0, cpu_t1         !< CPU time markers for time step profiling [s]
    real(kind=8) :: cpu_total_t0, cpu_total_t1 !< CPU time markers for total subroutine profiling [s]
    real(kind=8) :: cpu_prof_t0, cpu_prof_t1   !< Reusable CPU time markers for phase profiling [s]
    real(kind=8) :: cpu_mix_total          !< Accumulated CPU time for compute_c_mix [s]
    real(kind=8) :: cpu_update_old_total   !< Accumulated CPU time for update_old_attributes [s]
    real(kind=8) :: cpu_displace_total     !< Accumulated CPU time for particle displacement [s]
    real(kind=8) :: cpu_solver_total       !< Accumulated CPU time for reactive mixing solver [s]
    real(kind=8) :: cpu_derived_total      !< Accumulated CPU time for derived quantities [s]
    real(kind=8) :: cpu_eq_rates_total     !< Accumulated CPU time for equilibrium reaction rates [s]
    real(kind=8) :: cpu_solid_gas_total    !< Accumulated CPU time for solid/gas chemistry [s]
    real(kind=8) :: cpu_output_total       !< Accumulated CPU time for output writing [s]
    real(kind=8) :: cpu_total_elapsed      !< Total subroutine elapsed CPU time [s]
    real(kind=8) :: cpu_data_prep_total    !< Accumulated CPU time for data preparation [s]
    real(kind=8) :: cpu_getters_total      !< Accumulated CPU time for get_conc_nc + get_c1_old_old [s]
    real(kind=8) :: cpu_snapshot_total     !< Accumulated CPU time for write_chemistry snapshots [s]
    integer(kind=4) :: iunit_cpu           !< Fortran I/O unit number for CPU profiling output file
!> Concentration and reaction rate arrays
    real(kind=8), allocatable :: conc_comp(:)       !< Component concentrations after mixing [mol/L]
    real(kind=8), allocatable :: c_hat(:)            !< Mixed concentrations from transport [mol/L]
    real(kind=8), allocatable :: conc_nc(:)          !< Variable activity species concentrations [mol/L]
    real(kind=8), allocatable :: lumped_lambdas(:)   !< Lumped mixing ratios per target water [-]
    real(kind=8), allocatable :: c1_oo(:)            !< Persistent buffer for get_c1_old_old result [mol/L]
!> Pre-computed caches (invariant across time steps)
    integer(kind=4), allocatable :: n_p_cache(:)       !< Cached num_prim_species per non-canonical target water
    integer(kind=4), allocatable :: n_v_cache(:)       !< Cached num_var_act_species per non-canonical target water
    integer(kind=4), allocatable :: num_eq_cache(:)    !< Cached num_eq_reactions per non-canonical target water
    integer(kind=4), allocatable :: n_conc_cache(:)    !< Cached conc_old size per non-canonical target water
    integer(kind=4), allocatable :: n_idx_cache(:)     !< Cached indices_aq_species size per non-canonical target water
    integer(kind=4), allocatable :: all_ind_aq_sp(:,:) !< Cached indices_aq_species for ALL waters: all_ind_aq_sp(k,w) = position in water w's conc_old of aq_phase species k
    real(kind=8), allocatable :: all_conc_old(:,:)     !< Cached conc_old for ALL waters, refreshed each time step
    integer(kind=4) :: max_n_p     !< Max n_p across target waters (for pre-allocation)
    integer(kind=4) :: max_num_mix !< Max num_mix across target waters (for pre-allocation)
    integer(kind=4) :: max_n_conc  !< Max conc_old size across target waters (for pre-allocation)
    integer(kind=4) :: max_n_idx   !< Max indices_aq_species size across target waters (for pre-allocation)
    integer(kind=4) :: w           !< Loop index for waters
!> Procedure pointers
    procedure(reactive_mixing_iter_EI_kin_anal_ideal_opt2), pointer :: p_solver=>null()
    procedure(move_particles_stat_flux_2D), pointer :: p_displace_particles=>null()


!> -----------------------------------------------------------------------
!> Initialization
    !> Start total subroutine CPU timer for overall performance measurement
    call cpu_time(cpu_total_t0)
    !> Initialize total equilibrium reaction counter to zero
    n_eq=0
    !> Loop over all target waters to count total equilibrium reactions in the domain
    do i=1,this%chemistry%num_target_waters
        !> Accumulate equilibrium reactions from this target water's reactive zone
        n_eq=n_eq+this%target_waters(i)%aq_chem%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    !> End equilibrium reaction counting loop
    end do

!> -----------------------------------------------------------------------
!> Setup
    !> Initialize simulation time to zero [s]
    time=0d0
    !> Initialize output time step counter (starts at 2; index 1 is the initial condition)
    kk=2
    !> Initialize output target water counter to first water
    ii=1
    !> Cache total number of time steps locally for loop bounds
    num_time_loc = this%transport%time_discr%Num_time
    !> Cache number of output time steps from chemistry options
    num_out_ts = this%chemistry%chem_out_options%num_time_steps
    !> Cache number of output waters from chemistry options
    num_out_waters = this%chemistry%chem_out_options%num_waters
    !> Cache number of output aqueous species from chemistry options
    num_out_aq_sp = this%chemistry%chem_out_options%num_aq_species
    !> Force last output time step to coincide with the final simulation step
    this%chemistry%chem_out_options%time_steps(num_out_ts)=num_time_loc
    !> Validate that user-specified output time steps do not exceed total number of time steps
    do i=2,num_out_ts-1
        if (this%chemistry%chem_out_options%time_steps(i) > num_time_loc) then
            print *, "Output time step index", this%chemistry%chem_out_options%time_steps(i), &
                "exceeds total number of time steps", num_time_loc
            error stop "Mismatch between out_opts.dat and discr_temp.dat: output time step exceeds Num_time"
        end if
    end do
    !> Initialize chemistry snapshot counter (starts at 2; index 1 is initial condition)
    kk_snap=2
    !> Assign Fortran I/O unit number for main output file
    unit=7
    !> Open main output file for formatted sequential writing of species concentrations
    open(unit,file=dir//trim(root)//'_output.out',form="formatted",access="sequential",status="unknown")

!> -----------------------------------------------------------------------
!> Solver Selection
    !> Branch on chemical system type: mixed (eq+kin), equilibrium-only, or kinetic-only
    if (this%chemistry%chem_syst%num_kin_reacts>0 .and. n_eq>0) then
        !> Equilibrium + Kinetic Reactions
        !> Check rk averaging option for mixed equilibrium+kinetic system
        if (this%chemistry%rk_avg_opt==1) then
            !> Halt: rk averaging option 1 not yet implemented for mixed systems
            error stop "rk average option 1 not implemented yet"
        !> Use rk averaging option 2 for mixed equilibrium+kinetic system
        else if (this%chemistry%rk_avg_opt==2) then
            !> Assign mixed eq+kin analytical ideal solver with option 2
            p_solver=>reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2
        !> Unrecognized rk averaging option for mixed system
        else
            !> Halt: unknown rk averaging option
            error stop "rk average option not implemented yet"
        !> End rk averaging option selection for mixed system
        end if
    !> Check for equilibrium-only system (no kinetic reactions)
    else if (n_eq>0) then
        !> Equilibrium-Only Reactions
        !> Check if Euler Implicit (method 2) is selected for equilibrium reactions
        if (this%int_method_chem_reacts==2) then
            !> Assign equilibrium-only analytical ideal Euler Implicit solver
            p_solver=>reactive_mixing_iter_EI_eq_anal_ideal
        !> End integration method check for equilibrium-only system
        end if
    !> Kinetic-only system (no equilibrium reactions in domain)
    else
        !> Kinetic-Only Reactions
        !> Check rk averaging option for kinetic-only system
        if (this%chemistry%rk_avg_opt==1) then
            !> Halt: rk averaging option 1 not yet implemented for kinetic-only systems
            error stop "rk average option 1 not implemented yet"
        !> Use rk averaging option 2 for kinetic-only system
        else if (this%chemistry%rk_avg_opt==2) then
            !> Assign kinetic-only analytical ideal solver with option 2
            p_solver=>reactive_mixing_iter_EI_kin_anal_ideal_opt2
        !> Unrecognized rk averaging option for kinetic-only system
        else
            !> Halt: unknown rk averaging option
            error stop "rk average option not implemented yet"
        !> End rk averaging option selection for kinetic-only system
        end if
    !> End chemical system type branching
    end if

    !> Select particle displacement method based on spatial discretization target flag
    if (this%transport%spatial_discr%targets_flag==0) then
        !> Assign standard 2D stationary flux particle displacement routine
        p_displace_particles=>move_particles_stat_flux_2D
    !> Non-standard target flag: not yet supported
    else
        !> Halt: alternative Lagrangian schemes not implemented for stationary flux
        error stop "Stationary flux Lagrangian scheme not implemented yet"
    !> End particle displacement method selection
    end if

    !> Allocate array for lumped mixing ratios (one entry per mixing ratio column)
    allocate(lumped_lambdas(this%transport%mixing_ratios_R%num_cols))
    !> Loop over all mixing ratio columns to compute lumped (summed) values
    do i=1,this%transport%mixing_ratios_R%num_cols
        !> Sum all entries in mixing ratio column i to get the lumped mixing ratio
        lumped_lambdas(i)=sum(this%transport%mixing_ratios_R%cols(i)%col_1)
    !> End lumped mixing ratio computation loop
    end do

!> -----------------------------------------------------------------------
!> Canonical Vector Extraction
    !> Identify canonical and non-canonical vectors from mixing ratios using absolute tolerance
    call this%transport%mixing_ratios_conc%get_can_vec(this%chemistry%CV_params%abs_tol,num_can_vec,ind_can_vec,num_non_can_vec,&
        ind_non_can_vec)

!> -----------------------------------------------------------------------
!> Initialize Concentration History
    !> Loop over all waters to set initial previous-time-step concentration values
    do i=1,this%chemistry%num_waters
        !> Store current aqueous concentrations as "old" for time stepping
        call this%chemistry%waters(i)%set_conc_old()
        !> Store current aqueous concentrations as "old_old" (two steps back) for multi-step methods
        call this%chemistry%waters(i)%set_conc_old_old()
        !> Store current solid chemistry concentrations as "old" for time stepping
        call this%chemistry%waters(i)%solid_chemistry%set_conc_old()
        !> Store entire solid chemistry state as "old" for time stepping
        call this%chemistry%waters(i)%set_solid_chemistry_old()
    !> End concentration history initialization loop
    end do

!> -----------------------------------------------------------------------
!> Write Output File Header
    !> Write column headers with selected aqueous species names to the output file
    write(unit,"(2x,'Species :',5x,*(A15))") (this%chemistry%chem_syst%aq_phase%aq_species(&
        this%chemistry%chem_out_options%ind_aq_species(j))%name, &
        j=1,num_out_aq_sp)

!> -----------------------------------------------------------------------
!> Validate Water Ordering (one-time check)
    !> Loop over non-canonical vectors to verify mixing index consistency with target waters
    do i=1,num_non_can_vec
        !> Cache the non-canonical vector index for this iteration
        idx = ind_non_can_vec(i)
        !> Check if the first mixing index matches the expected target water index
        if (this%transport%mix_conc_indices%cols(idx)%col_1(1) /= this%chemistry%tar_wat_indices(idx)) then
            !> Print the mismatched mixing index and target water index for debugging
            print *, this%transport%mix_conc_indices%cols(idx)%col_1(1), this%chemistry%tar_wat_indices(idx)
            !> Print the problematic target water loop index
            print *, "Problematic target water index: ", i
            !> Halt: water ordering inconsistency detected between transport and chemistry
            error stop "Waters not in the right order in mixing indices"
        !> End water ordering mismatch check
        end if
    !> End water ordering validation loop
    end do

!> -----------------------------------------------------------------------
!> Pre-compute per-water metadata and caches
    !> Cache total number of waters locally for loop bounds
    num_waters_loc = this%chemistry%num_waters
    !> Allocate cache array for number of primary species per non-canonical target water
    allocate(n_p_cache(num_non_can_vec))
    !> Allocate cache array for number of variable activity species per non-canonical target water
    allocate(n_v_cache(num_non_can_vec))
    !> Allocate cache array for number of equilibrium reactions per non-canonical target water
    allocate(num_eq_cache(num_non_can_vec))
    !> Allocate cache array for conc_old size per non-canonical target water
    allocate(n_conc_cache(num_non_can_vec))
    !> Allocate cache array for indices_aq_species size per non-canonical target water
    allocate(n_idx_cache(num_non_can_vec))
    !> Initialize maximum dimension trackers to zero for pre-allocation sizing
    max_n_p = 0; max_num_mix = 0; max_n_conc = 0; max_n_idx = 0
    !> Loop over non-canonical vectors to populate per-water caches and find max dimensions
    do i = 1, num_non_can_vec
        !> Retrieve the global target water index for this non-canonical vector
        idx = ind_non_can_vec(i)
        !> Associate local aliases for speciation algebra and mixing column of this target water
        associate(spec_alg => this%target_waters(idx)%aq_chem%solid_chemistry%reactive_zone%speciation_alg, &
                  mix_col => this%transport%mix_conc_indices%cols(idx)%col_1)
        !> Cache number of primary species for this target water
        n_p_cache(i) = spec_alg%num_prim_species
        !> Cache number of variable activity species for this target water
        n_v_cache(i) = spec_alg%num_var_act_species
        !> Cache number of equilibrium reactions for this target water
        num_eq_cache(i) = spec_alg%num_eq_reactions
        !> Cache size of conc_old array for the first mixing water
        n_conc_cache(i) = size(this%chemistry%waters(mix_col(1))%conc_old)
        !> Cache size of indices_aq_species array for the first mixing water
        n_idx_cache(i) = size(this%chemistry%waters(mix_col(1))%indices_aq_species)
        !> Compute number of mixing waters (column dim minus 3 metadata entries)
        num_mix_loc = this%transport%mix_conc_indices%cols(idx)%dim - 3
        !> Update max primary species count if this water exceeds current max
        if (n_p_cache(i) > max_n_p) max_n_p = n_p_cache(i)
        !> Update max mixing waters count if this water exceeds current max
        if (num_mix_loc > max_num_mix) max_num_mix = num_mix_loc
        !> Update max conc_old size if this water exceeds current max
        if (n_conc_cache(i) > max_n_conc) max_n_conc = n_conc_cache(i)
        !> Update max indices_aq_species size if this water exceeds current max
        if (n_idx_cache(i) > max_n_idx) max_n_idx = n_idx_cache(i)
        !> End associate block for speciation algebra and mixing column aliases
        end associate
    !> End per-water metadata caching loop
    end do
    !> Cache indices_aq_species for ALL waters in contiguous 2D array (never changes)
    !> all_ind_aq_sp(k, w) = indices_aq_species(k) of water w = position in water w's conc_old of aq_phase species k
    allocate(all_ind_aq_sp(max_n_idx, num_waters_loc))
    do w = 1, num_waters_loc
        all_ind_aq_sp(1:size(this%chemistry%waters(w)%indices_aq_species), w) = &
            this%chemistry%waters(w)%indices_aq_species
    end do
    !> Pre-allocate working arrays to maximum sizes
    !> Allocate component concentration working array to maximum primary species count
    allocate(conc_comp(max_n_p))
    !> Allocate global conc_old 2D cache array (max conc size x number of waters)
    allocate(all_conc_old(max_n_conc, num_waters_loc))
    !> Initialise conditional-reallocation trackers
    !> Set previous n_v to zero to force initial allocation of conc_nc and c1_oo
    prev_n_v = 0
    !> Cache reaction time-weighting parameter theta_r locally for inner loop access
    theta_r_loc = this%transport%time_discr%theta_r
    !> Initialize all CPU profiling accumulators to zero [s]
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
    !> Pre-compute mineral/gas presence flags (invariant across time steps)
    !> Allocate boolean flag array for mineral zone presence per non-canonical vector
    allocate(has_minerals_flag(num_non_can_vec))
    !> Allocate boolean flag array for gas chemistry presence per non-canonical vector
    allocate(has_gases_flag(num_non_can_vec))
    !> Loop over non-canonical vectors to check mineral/gas associations per target water
    do i = 1, num_non_can_vec
        !> Retrieve the global target water index for this non-canonical vector
        idx = ind_non_can_vec(i)
        !> Check and store whether this target water has an associated mineral zone
        has_minerals_flag(i) = associated(this%target_waters(idx)%aq_chem%solid_chemistry%mineral_zone)
        !> Check and store whether this target water has an associated gas chemistry object
        has_gases_flag(i) = associated(this%target_waters(idx)%aq_chem%gas_chemistry)
    !> End mineral/gas flag pre-computation loop
    end do
    !> Initialize profiled time step counter to zero
    num_profiled_steps = 0
    !> Open CPU profiling output file (replaces any existing file)
    open(newunit=iunit_cpu, file=dir//trim(root)//'_cpu_profile.out', status='replace', action='write')

!> -----------------------------------------------------------------------
!> Main Time Loop
    !> Begin main time stepping loop from step 1 to total number of time steps
    do k=1,num_time_loc
        !> Record CPU time at start of this time step for per-step profiling
        call cpu_time(cpu_t0)
        !> Compute time step size for step k from the time discretization
        Delta_t=this%transport%time_discr%get_Delta_t(k)
        !> Advance simulation time by the current time step size
        time=time+Delta_t
        !> Print progress message every ~5% of total steps and at the first step

        !> Check if current time step matches the next designated output time step
        is_output_step = (k == this%chemistry%chem_out_options%time_steps(kk))
        !> If this is an output step, write timing header to the output file
        if (is_output_step) then
           !> Write formatted simulation time header to the output file
           write(unit,"(/,2x,'t = ',*(ES15.5),/)") time
        !> End output step time header writing
        end if
        !> Update Previous Time Step Values
        !> Profile: start CPU timer for update_old_attributes phase
        call cpu_time(cpu_prof_t0)
        !> Loop over all waters to update their previous-time-step attributes
        do i=1,num_waters_loc
            !> Shift current concentrations and state to "old" storage for this water
            call this%chemistry%waters(i)%update_old_attributes()
        !> End update_old_attributes loop over waters
        end do
        !> Profile: record CPU time after update_old_attributes phase
        call cpu_time(cpu_prof_t1)
        !> Profile: accumulate elapsed time for update_old_attributes phase
        cpu_update_old_total = cpu_update_old_total + (cpu_prof_t1 - cpu_prof_t0)
        !> Profile: start CPU timer for particle displacement phase
        call cpu_time(cpu_prof_t0)
        !> Move particles according to stationary flux field for time step k
        call p_displace_particles(this,k)
        !> Introduce new particles (e.g. boundary inflow) at time step k
        call this%introduce_particles(k)
        !> Profile: record CPU time after particle displacement phase
        call cpu_time(cpu_prof_t1)
        !> Profile: accumulate elapsed time for particle displacement phase
        cpu_displace_total = cpu_displace_total + (cpu_prof_t1 - cpu_prof_t0)
        !> Refresh global conc_old cache AFTER displacement (post-displacement conc_old)
        do w = 1, num_waters_loc
            all_conc_old(1:size(this%chemistry%waters(w)%conc_old), w) = &
                this%chemistry%waters(w)%conc_old
        end do
        !> Determine whether to profile the inner loop this time step (sample ~5% of steps)
        do_profile_inner = (mod(k, max(num_time_loc/20, 1)) == 0) .or. (k <= 2) .or. (k == num_time_loc)
        !> If profiling is active this step, increment the profiled step counter
        if (do_profile_inner) num_profiled_steps = num_profiled_steps + 1

        !> Target Waters Loop (Spatial Domain)
        do i=1,num_non_can_vec
            !> Retrieve the global target water index for this non-canonical vector
            idx = ind_non_can_vec(i)
            !> Associate tw as a shorthand alias for this target water's aqueous chemistry
            associate(tw => this%target_waters(idx)%aq_chem)
            !> Profile: start data preparation timer if profiling is active this step
            if (do_profile_inner) call cpu_time(cpu_prof_t0)
            !> Load cached number of primary species for this target water
            n_p = n_p_cache(i)
            !> Load cached number of variable activity species for this target water
            n_v = n_v_cache(i)
            !> Load cached number of equilibrium reactions for this target water
            num_eq_loc = num_eq_cache(i)
            !> Load cached lumped mixing ratio for this target water
            lumped_lambda_loc = lumped_lambdas(idx)
            !> Transport Mixing Step
            associate(mix_col => this%transport%mix_conc_indices%cols(idx)%col_1)
            !> Compute number of mixing waters (column dimension minus 3 metadata entries)
            num_mix_loc = this%transport%mix_conc_indices%cols(idx)%dim - 3
            !> Profile: end data prep timer and accumulate elapsed time
            if (do_profile_inner) then
                !> Profile: record CPU time after data preparation phase
                call cpu_time(cpu_prof_t1)
                !> Profile: accumulate elapsed data preparation CPU time
                cpu_data_prep_total = cpu_data_prep_total + (cpu_prof_t1 - cpu_prof_t0)
                !> Profile: reset start timer for next phase (conservative mixing)
                cpu_prof_t0 = cpu_prof_t1
            !> End data preparation profiling guard
            end if
            !> Compute conservative mixing of concentrations using mixing ratios and cached conc_old
            call compute_c_mix_global(this%chemistry%waters(mix_col(1)), &
                all_conc_old, all_ind_aq_sp, mix_col(2:), num_mix_loc, &
                this%transport%mixing_ratios_conc%cols(idx)%col_1, c_hat)
            !> Profile: end conservative mixing timer and accumulate elapsed time
            if (do_profile_inner) then
                !> Profile: record CPU time after conservative mixing phase
                call cpu_time(cpu_prof_t1)
                !> Profile: accumulate elapsed conservative mixing CPU time
                cpu_mix_total = cpu_mix_total + (cpu_prof_t1 - cpu_prof_t0)
                !> Profile: reset start timer for next phase (getters)
                cpu_prof_t0 = cpu_prof_t1
            !> End conservative mixing profiling guard
            end if
            !> End associate block for mixing column alias
            end associate
            !> Conditionally reallocate conc_nc only when n_v changes
            if (n_v /= prev_n_v) then
                !> Deallocate conc_nc if currently allocated (size mismatch)
                if (allocated(conc_nc)) deallocate(conc_nc)
                !> Deallocate c1_oo if currently allocated (size mismatch)
                if (allocated(c1_oo)) deallocate(c1_oo)
                !> Update tracked n_v to current value to avoid unnecessary future reallocations
                prev_n_v = n_v
            !> End conditional reallocation check
            end if
            !> Retrieve current variable activity species concentrations from target water
            conc_nc=tw%get_conc_nc()
            !> Retrieve old-old (two steps back) primary species concentrations from target water
            c1_oo=tw%get_c1_old_old()
            !> Profile: end getters timer and accumulate elapsed time
            if (do_profile_inner) then
                !> Profile: record CPU time after getter calls
                call cpu_time(cpu_prof_t1)
                !> Profile: accumulate elapsed getter CPU time
                cpu_getters_total = cpu_getters_total + (cpu_prof_t1 - cpu_prof_t0)
                !> Profile: reset start timer for next phase (reactive mixing solver)
                cpu_prof_t0 = cpu_prof_t1
            !> End getter profiling guard
            end if
            !> Reactive Mixing Iteration
            call p_solver(tw, c1_oo, c_hat(1:n_v), &
                lumped_lambda_loc, lumped_lambda_loc, &
                Delta_t, theta_r_loc, conc_nc, conc_comp)
            !> Profile: end reactive mixing solver timer and accumulate elapsed time
            if (do_profile_inner) then
                !> Profile: record CPU time after reactive mixing solver
                call cpu_time(cpu_prof_t1)
                !> Profile: accumulate elapsed reactive mixing solver CPU time
                cpu_solver_total = cpu_solver_total + (cpu_prof_t1 - cpu_prof_t0)
                !> Profile: reset start timer for next phase (derived quantities)
                cpu_prof_t0 = cpu_prof_t1
            !> End reactive mixing solver profiling guard
            end if
            !> Update Derived Quantities
            !> Compute pH from current aqueous concentrations
            call tw%compute_pH()
            !> Compute salinity from current aqueous concentrations
            call tw%compute_salinity()
            !> Compute ionic strength from current aqueous concentrations
            call tw%compute_ionic_strength()
            !> Profile: end derived quantities timer and accumulate elapsed time
            if (do_profile_inner) then
                !> Profile: record CPU time after derived quantity calculations
                call cpu_time(cpu_prof_t1)
                !> Profile: accumulate elapsed derived quantities CPU time
                cpu_derived_total = cpu_derived_total + (cpu_prof_t1 - cpu_prof_t0)
                !> Profile: reset start timer for next phase (equilibrium rates)
                cpu_prof_t0 = cpu_prof_t1
            !> End derived quantities profiling guard
            end if
            !> Accumulate Reaction Amounts
            !> Add current kinetic reaction amounts to running accumulator for aqueous chemistry
            tw%Rk_accum=tw%Rk_accum+tw%Rk
            !> Add current kinetic reaction amounts to running accumulator for solid chemistry
            tw%solid_chemistry%Rk_accum=tw%solid_chemistry%Rk_accum+tw%solid_chemistry%Rk
            !> Equilibrium Reaction Rates
            !> Cache number of kinetic reaction rate columns for this target water
            num_kin_rk_cols = tw%indices_rk%num_cols
            !> If both equilibrium reactions and kinetic rate columns exist, compute mixed eq rates
            if (num_eq_loc>0 .and. num_kin_rk_cols>0) then
                !> Compute equilibrium reaction rates accounting for kinetic coupling
                call tw%compute_Re_kin(c_hat(n_p+1:n_v), Delta_t, lumped_lambda_loc)
            !> If only equilibrium reactions exist (no kinetic rate columns), compute pure eq rates
            else if (num_eq_loc>0) then
                !> Compute equilibrium reaction rates from mixed secondary concentrations
                call tw%compute_Re(c_hat(1:n_v), Delta_t, lumped_lambda_loc)
            !> End equilibrium reaction rate computation branching
            end if
            !> Profile: end equilibrium rates timer and accumulate elapsed time
            if (do_profile_inner) then
                !> Profile: record CPU time after equilibrium rate calculations
                call cpu_time(cpu_prof_t1)
                !> Profile: accumulate elapsed equilibrium rates CPU time
                cpu_eq_rates_total = cpu_eq_rates_total + (cpu_prof_t1 - cpu_prof_t0)
                !> Profile: reset start timer for next phase (solid/gas chemistry)
                cpu_prof_t0 = cpu_prof_t1
            !> End equilibrium rates profiling guard
            end if
            !> Solid Chemistry State Variables
            if (has_minerals_flag(i)) then
                !> Compute mineral mass balances from kinetic reaction amounts and Delta_t
                call tw%solid_chemistry%compute_mass_bal_mins(Delta_t)
                !> Update mineral concentrations iteratively based on mass balance results
                call tw%solid_chemistry%compute_conc_minerals_iter(Delta_t)
            !> End mineral chemistry update conditional
            end if
            !> Gas Chemistry State Variables
            if (has_gases_flag(i)) then
                !> Compute gas concentrations iteratively from reaction rates and cell volume
                call tw%gas_chemistry%compute_conc_gases_iter( &
                    Delta_t, tw%volume, [tw%re_mean, tw%rk_mean])
                !> Compute volumetric gas species concentrations from updated gas amounts
                call tw%gas_chemistry%compute_vol_gas_species_conc()
                !> Compute logarithmic activity coefficients for gas species
                call tw%gas_chemistry%compute_log_act_coeffs_gases()
            !> End gas chemistry update conditional
            end if
            !> Profile: end solid/gas chemistry timer and accumulate elapsed time
            if (do_profile_inner) then
                !> Profile: record CPU time after solid/gas chemistry calculations
                call cpu_time(cpu_prof_t1)
                !> Profile: accumulate elapsed solid/gas chemistry CPU time
                cpu_solid_gas_total = cpu_solid_gas_total + (cpu_prof_t1 - cpu_prof_t0)
                !> Profile: reset start timer for next phase (output writing)
                cpu_prof_t0 = cpu_prof_t1
            !> End solid/gas chemistry profiling guard
            end if
            !> Write Chemistry Output
            if (is_output_step .and. &
                (idx==this%chemistry%chem_out_options%ind_waters(ii))) then
                !> Write water index and selected species concentrations to the output file
                write(unit,"(I10,*(ES15.5))") this%chemistry%chem_out_options%ind_waters(ii), &
                    (tw%concentrations(tw%indices_aq_species( &
                    this%chemistry%chem_out_options%ind_aq_species(j))), &
                    j=1,num_out_aq_sp)
                !> If more waters remain to output at this time step, advance water counter
                if (ii<num_out_waters) then
                    !> Advance to the next output water
                    ii=ii+1
                !> If all waters done but more output time steps remain, advance time step counter
                else if (kk<num_out_ts) then
                    !> Advance to the next output time step
                    kk=kk+1
                    !> Reset output water counter to first water for the new output time step
                    ii=1
                !> End output counter advancement branching
                end if
            !> End chemistry output writing conditional
            end if
            !> Profile: end output writing timer and accumulate elapsed time
            if (do_profile_inner) then
                !> Profile: record CPU time after output writing
                call cpu_time(cpu_prof_t1)
                !> Profile: accumulate elapsed output writing CPU time
                cpu_output_total = cpu_output_total + (cpu_prof_t1 - cpu_prof_t0)
            !> End output writing profiling guard
            end if
            !> Deallocate temporary mixed concentration array (re-allocated next iteration)
            deallocate(c_hat)
            !> End associate block for target water aqueous chemistry alias
            end associate
        !> End target waters inner loop
        end do

        !> Write chemistry snapshot at output time steps
        !> (Snapshots disabled — write_chemistry is called once at the end via write_RT)
        !> Record CPU time at end of this time step for per-step profiling
        call cpu_time(cpu_t1)
        !> Print elapsed CPU time for this individual time step
        write(*,'(A,I10,A,F14.4,A)') '[CPU] time step ', k, ' took ', cpu_t1 - cpu_t0, ' s'
    !> End main time stepping loop
    end do

!> -----------------------------------------------------------------------
!> Final Cleanup
    !> Safely deallocate component concentration working array if allocated
    if (allocated(conc_comp)) deallocate(conc_comp)
    !> Safely deallocate variable-activity species concentration array if allocated
    if (allocated(conc_nc)) deallocate(conc_nc)
    !> Safely deallocate old-old primary species concentration buffer if allocated
    if (allocated(c1_oo)) deallocate(c1_oo)
    !> Safely deallocate mineral presence flag array if allocated
    if (allocated(has_minerals_flag)) deallocate(has_minerals_flag)
    !> Safely deallocate gas presence flag array if allocated
    if (allocated(has_gases_flag)) deallocate(has_gases_flag)
    !> Record CPU time at the end of the subroutine for total elapsed calculation
    call cpu_time(cpu_total_t1)
    !> Compute total elapsed CPU time for the entire subroutine
    cpu_total_elapsed = cpu_total_t1 - cpu_total_t0
    !> Write CPU profiling report to file and close unit
    call write_cpu_profile_report( &
        iunit_cpu, &
        num_profiled_steps, num_time_loc, num_non_can_vec, &
        cpu_total_elapsed, &
        cpu_update_old_total, cpu_displace_total, &
        cpu_data_prep_total, cpu_mix_total, cpu_getters_total, &
        cpu_solver_total, cpu_derived_total, cpu_eq_rates_total, &
        cpu_solid_gas_total, cpu_output_total, cpu_snapshot_total, &
        dir, trim(root))
    !> Close the main species concentration output file
    close(unit)
!> End of subroutine
end subroutine
