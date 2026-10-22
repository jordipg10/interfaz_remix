!> \file solve_RT_1D_ideal_cons_Euler_stat_flux.f90
!> \brief Solve 1D reactive transport with ideal speciation, lumping, and Lagrangian particle tracking
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
!> \param[inout] this Reactive transport 1D object (transient or stationary)
!> \param[in] root Root path for output files
!> \author Jordi Petchamé-Guerrero
!> \date November 2025
subroutine solve_RT_ideal_cons_Euler_stat_flux_1D(this,dir,root)
    use RT_m, only: RT_1D_transient_c, RT_1D_stat_c, RT_c !< Import 1D reactive transport classes
    use chemistry_m, only: write_out_step_vars !< Shared per-output-step writer (conc/react_rate/vol_fract)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, write_celldbg_cell, & !< Import aqueous chemistry class
        reactive_mixing_iter_EE_eq_kin_ideal, reactive_mixing_iter_EE_kin, & !< Euler explicit solvers with lumping
        reactive_mixing_iter_EI_kin_anal_ideal_opt2, mixing_iter_comp_ideal, & !< Implicit solvers for kinetic-only and equilibrium-only systems
        compute_r_tilde_impl_opt1_idx, compute_r_tilde_impl_opt2_idx, compute_c_mix_by_zones_cons, compute_c_mix_global_cons, &
        compute_r_tilde_impl_opt3_idx, compute_r_tilde_impl_opt4_idx, compute_c_hat_patankar !< Utility functions for mixing calculations
    use arrays_m, only: int_array_c
    implicit none !< Enforce explicit variable declarations for type safety
!> @name Arguments
!> @{
    class(RT_1D_transient_c) :: this !< Reactive transport 1D object (polymorphic: transient or stationary) [-]
    character(len=*), intent(in) :: dir !< Directory for output file [string]
    character(len=*), intent(in) :: root !< Root name for output file (base filename without extension) [string]
!> @}
!> @name Loop Counters and Indices
!> @{
    integer(kind=4) :: i !< Counter for target waters (domain cells) [-]
    integer(kind=4) :: m !< Counter for target waters (domain cells) [-]
    integer(kind=4) :: j !< Counter for target solids or species indices [-]
    integer(kind=4) :: l !< Counter for reactive zones (spatial geochemical regions) [-]
    integer(kind=4) :: k !< Main time step counter for transport loop (1 to Num_time) [-]
    integer(kind=4) :: kk !< Output time step counter for chemistry output options [-]
    integer(kind=4) :: ii !< Output target water counter for chemistry output options [-]
!> @}
!> @name Problem Dimensions
!> @{
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
!> @}
!> @name Index Arrays
!> @{
    integer(kind=4), allocatable :: tar_gas_indices(:) !< Indices of target gases in each reactive zone (maps gases to zones) [-]
    integer(kind=4), allocatable :: tar_sol_indices(:) !< Indices of target solids in each reactive zone (maps solids to zones) [-]
    integer(kind=4), allocatable :: tar_wat_indices(:) !< Indices of target waters in each reactive zone (maps waters to zones) [-]
    integer(kind=4), allocatable :: perm(:) !< Permutation vector for reordering aqueous concentrations [-]
    integer(kind=4), allocatable :: ind_can_vec(:) !< Indices of canonical vectors in mixing ratios matrix (basis vectors) [-]
    integer(kind=4), allocatable :: ind_non_can_vec(:) !< Indices of non-canonical vectors requiring reactive transport solve [-]
!> @}
!> @name Time and Mixing Parameters
!> @{
    REAL(KIND=8) :: time !< Current simulation time [s] (accumulated from t=0)
    REAL(KIND=8) :: Delta_t !< Time step size [s] (can vary between steps)
    REAL(KIND=8) :: theta !< Time weighting factor [-]: 0=explicit, 1=implicit, 0.5=Crank-Nicolson
    REAL(KIND=8) :: y !< Sum of upstream mixing ratios [-] (for mass balance validation)
    real(kind=8) :: lambda_R !< Modified mixing ratio for reaction rates [-]
!> @}
!> @name Concentration and Reaction Rate Arrays
!> @{
    REAL(KIND=8), allocatable :: conc_comp(:) !< Variable activity mobile species concentrations after transport mixing [mol/L]
    real(kind=8), allocatable :: lambdas_R(:) !< Mixed kinetic reaction rates for current cell [mol/L/s]
    REAL(KIND=8), allocatable :: c_hat(:) !< Concentrations from previous time step (before mixing) [mol/L]
    REAL(KIND=8), allocatable :: c_mix(:) !< Concentrations from previous time step (before mixing) [mol/L]
    REAL(KIND=8), allocatable :: c_tilde(:) !< Conservative part of c_mix (no reaction rates) [mol/L]
    REAL(KIND=8), allocatable :: conc_nc(:) !< Variable activity species concentrations (solution vector) [mol/L]
    REAL(KIND=8), allocatable :: r_tilde(:) !< Variable activity species concentrations (solution vector) [mol/L]
    !> Domain-integrated reactive-mixing mass balance (per variable-activity species)
    integer(kind=4) :: iunit_celldbg              !< per-cell [CELL-DBG] trace file unit
    integer(kind=4) :: idx_ref                    !< reference target-water index (species names/dimension)
    integer(kind=4) :: n_v_ref                    !< reference number of variable-activity species
    real(kind=8), allocatable :: mb_lhs(:), mb_transp(:), mb_react(:) !< domain-summed terms
    real(kind=8), allocatable :: mb_lhs_c(:), mb_transp_c(:), mb_react_c(:) !< per-cell terms
!> @}
!> @name Chemistry State Objects
!> @{
    type(aqueous_chemistry_c), allocatable :: waters_new(:) !< Target waters at time step k+1 (after reaction step) [-]
    type(aqueous_chemistry_c), allocatable :: waters_old(:) !< Target waters at time step k (current time) [-]
    type(aqueous_chemistry_c), allocatable :: waters_old_old(:) !< Target waters at time step k-1 (for multistep methods) [-]
    type(aqueous_chemistry_c), allocatable :: mix_waters(:) !< Array of mixing waters (upstream contributors) [-]
    integer(kind=4) :: num_lump !< Number of lumping operations performed (for efficiency tracking) [-]
    integer(kind=4) :: j_mix !< loop index for mixing waters
    integer(kind=4) :: num_mix_loc !< number of mixing waters
    real(kind=8), allocatable :: conc_old_mix(:,:) !< conc_old from each mixing water
    integer(kind=4), allocatable :: ind_aq_sp_mix(:,:) !< indices_aq_species from each mixing water
!> @}
!> @name Procedure Pointers (Runtime Polymorphism)
!> @{
    !> Procedure pointer for reactive mixing solver (selected at runtime based on chemical system)
    procedure(reactive_mixing_iter_EI_kin_anal_ideal_opt2), pointer :: p_solver=>null() !< Pointer to reactive mixing iteration solver
    !> Procedure pointer for kinetic rate mixing (selected based on integration method)
    procedure(compute_r_tilde_impl_opt4_idx), pointer :: p_r_tilde=>null() !< Pointer to index-based kinetic rate mixing function
    procedure(compute_c_mix_by_zones_cons), pointer :: p_mixer
!> @}
!> @name By-reactive-zone previous-step map (rebuilt each time step)
!> @{
    type(int_array_c) :: wat_indices_zones !< cols(jz)%col_1 = global target-water indices in old reaction zone jz
    integer(kind=4) :: num_rz !< number of old reaction zones (set by link_target_waters_old_reaction_zones)
!> @}
!> @section init_section Initialization
!> Initialize target waters for time stepping
    if (this%chemistry%mix_by_zones_flag) then
        p_mixer => compute_c_mix_by_zones_cons
    else
        p_mixer => compute_c_mix_global_cons
    end if
    allocate(waters_old(this%chemistry%num_waters)) !< Allocate array for old target waters
    !allocate(waters_old_old(this%chemistry%num_waters)) !< Allocate array for old old target waters
    n_eq=0 !< Initialize total number of equilibrium reactions
    do i=1,this%chemistry%num_waters !< Loop over all target waters (domain + external)
        n_eq=n_eq+this%chemistry%waters(i)%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions !< Count total equilibrium reactions
        call waters_old(i)%copy_aq_chem(this%chemistry%waters(i)) !< Initialize time stepping variables
        !call waters_old_old(i)%copy_aq_chem(waters_old(i)) !< Initialize time stepping variables
    end do
    !waters_old=waters_old !< Copy current state to "old" (time level k)
    !waters_old_old=waters_old !< Copy to "old old" (time level k-1) for multistep schemes

!> @section type_select Type Selection and Setup
!> Select transient reactive transport type (polymorphic dispatch)
select type (this)
type is (RT_1D_transient_c) !< Transient (time-dependent) reactive transport
    time=0d0 !< Initialize simulation time to zero [s]
    k_flux=1 !< Initialize flux change counter (for time-varying velocity fields) [-]
    kk=2 !< Initialize output time step counter (start at 2, as 1 is initial condition) [-]
    ii=1 !< Initialize output target water counter (start at first output location) [-]
    unit=7 !< Assign Fortran I/O unit number for output file (arbitrary choice) [-]
    !< Open output file for formatted sequential writing
    block
        character(len=:), allocatable :: rtag
        if (this%chemistry%lump_flag) then
            rtag = '_Euler_lumped'
        else
            rtag = '_Euler_consistent'
        end if
        open(unit,file=dir//trim(root)//rtag//'_output.out',form="formatted",access="sequential",status="unknown")
    end block
!> @section solver_select Solver Selection
!> Select reactive mixing subroutine based on:
!> - Chemical system type (equilibrium/kinetic/both)
!> - Time integration method (explicit/implicit)
!> - Jacobian calculation method (analytical/numerical)
!> @subsection eq_kin_case Equilibrium + Kinetic Reactions
    if (this%chemistry%chem_syst%num_kin_reacts>0 .and. n_eq>0) then !< Both equilibrium and kinetic reactions present
        !if (this%int_method_chem_reacts==1) then !< Integration method 1: Euler Explicit (EE)
            !theta=0d0 !< Set θ=0 for fully explicit time integration [-]
            !p_r_tilde=>p_r_tilde_impl_opt1 !< Assign explicit rate mixing function
            !mix_ind=1 !< Mix all waters including downstream (explicit allows this) [-]
        !else if (this%chemistry%Jac_opt==1) then !< Jacobian option 1: Analytical Jacobian (exact derivatives)
            !mix_ind=2 !< Mix all waters except downstream (implicit requires this to avoid circular dependency) [-]
            ! if (this%int_method_chem_reacts==2) then !< Integration method 2: Euler Implicit (EI)
            !     theta=1d0 !< Set θ=1 for fully implicit time integration [-]
            ! else if (this%int_method_chem_reacts==3) then !< Integration method 3: Crank-Nicolson (CN)
            !     theta=5d-1 !< Set θ=0.5 for Crank-Nicolson (second-order accurate) [-]
            ! else
            !     error stop "Integration method for chemical reactions not implemented yet" !< Unsupported method
            ! end if
            !< Select kinetic rate averaging option
            if (this%chemistry%rk_avg_opt==1) then !< Option 1: Simple averaging of kinetic rates
                error stop "rk average option 1 not implemented yet" !< Not implemented
            else if (this%chemistry%rk_avg_opt==2) then !< Option 2: Improved kinetic rate averaging with lumping
                !p_solver=>reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2 !< Assign implicit solver with lumping
            else
                error stop "rk average option not implemented yet" !< Unsupported option
            end if
        ! else if (this%int_method_chem_reacts==3 .and. this%chemistry%Jac_opt==1) then !< Crank-Nicolson not implemented
        !     error stop "Crank-Nicolson for equilibrium+kinetic not implemented yet"
        ! else
        !     error stop "Integration method for chemical reactions not implemented yet"
        ! end if
!> @subsection eq_only_case Equilibrium-Only Reactions
    else if (n_eq>0) then !< Only equilibrium reactions (no kinetics)
        !p_solver=>mixing_iter_comp_ideal !< Assign equilibrium-only solver with component lumping
!> @subsection kin_only_case Kinetic-Only Reactions
    else !< Only kinetic reactions (no equilibrium)
        ! if (this%int_method_chem_reacts==1) then !< Integration method 1: Euler Explicit (EE)
        !     theta=0d0 !< Set θ=0 for fully explicit time integration [-]
        !     p_solver=>reactive_mixing_iter_EE_kin !< Assign explicit kinetic solver with lumping
        !     mix_ind=1 !< Mix all waters including downstream (explicit allows this) [-]
        ! else if (this%chemistry%Jac_opt==1) then !< Jacobian option 1: Analytical Jacobian
        !     if (this%int_method_chem_reacts==2) then !< Integration method 2: Euler Implicit (EI)
        !         theta=1d0 !< Set θ=1 for fully implicit time integration [-]
        !     else if (this%int_method_chem_reacts==3) then !< Integration method 3: Crank-Nicolson (CN)
        !         theta=5d-1 !< Set θ=0.5 for Crank-Nicolson (second-order accurate) [-]
        !     else
        !         error stop "Integration method for chemical reactions not implemented yet" !< Unsupported method
        !     end if
            !< Select kinetic rate averaging option for kinetic-only system
            if (this%chemistry%rk_avg_opt==1) then !< Option 1: Concentration averaging
                error stop "rk average option 1 not implemented yet" !< Not implemented
            else if (this%chemistry%rk_avg_opt==2) then !< Option 2: Reaction rate averaging
                p_solver=>reactive_mixing_iter_EI_kin_anal_ideal_opt2 !< Assign implicit kinetic solver with lumping
            else
                error stop "rk average option not implemented yet" !< Unsupported option
            end if
            !< Select estimation option for kinetic-only system
            if (this%chemistry%r_down_opt==1) then !< Option 1: Use previous time step values
                !mix_ind=2 !< Mix all waters except current
                p_r_tilde=>compute_r_tilde_impl_opt1_idx !< Assign index-based rate mixing (closest-upstream downstream estimate)
                !p_solver=>reactive_mixing_iter_EI_kin_anal_ideal_opt2 !< Assign implicit kinetic solver with lumping
                !mix_ind=1
                !compute_r_tilde=>compute_r_tilde_impl_opt1_bis !< Assign implicit rate mixing function with previous step values
                !p_solver=>reactive_mixing_iter_EI_kin_anal_ideal_opt2_bis !< Assign implicit kinetic solver with lumping
            else if (this%chemistry%r_down_opt==2) then !< Option 2: Reaction rate averaging
                p_r_tilde=>compute_r_tilde_impl_opt2_idx !< Assign index-based rate mixing (first-derivative downstream estimate)
            else if (this%chemistry%r_down_opt==3) then !< Option 3: No estimation (zero downstream rate)
                p_r_tilde=>compute_r_tilde_impl_opt3_idx !< Assign index-based rate mixing (second-derivative downstream estimate)
            else if (this%chemistry%r_down_opt==4) then !< Option 4: Extrapolation based on upstream rates
                p_r_tilde=>compute_r_tilde_impl_opt4_idx !< Assign index-based rate mixing (Eulerian consistent, current-cell solid/gas)
                !p_solver=>reactive_mixing_iter_EI_kin_anal_ideal_opt2 !< Assign implicit kinetic solver with lumping
            else
                error stop "Downstream reaction rate estimation option not implemented yet" !< Unsupported option
            end if
            !mix_ind=2 !< Mix all waters except downstream (implicit constraint) [-]
        ! else if (this%int_method_chem_reacts==3 .and. this%chemistry%Jac_opt==1) then !< Crank-Nicolson not implemented
        !     error stop "Crank-Nicolson for kinetic-only not implemented yet"
        ! else
        !     error stop "Integration method for chemical reactions not implemented yet"
        ! end if
    end if
    
    !allocate(lumped_lambdas(this%transport%mixing_ratios_R%num_cols)) !< Allocate array for lumped kinetic rates
    !do i=1,this%transport%mixing_ratios_R%num_cols !< Loop over mixing ratios columns
    !    lumped_lambdas(i)=sum(this%transport%mixing_ratios_R%cols(i)%col_1) !< Sum mixing ratios for mass balance check
    !    if (abs(lumped_lambdas(i)-1d0)>this%chemistry%CV_params%abs_tol) then !< Validate sum to unity within tolerance
    !        print *, "Mixing ratios do not sum to 1 in column ", i, " sum = ", lumped_lambdas(i)
    !        error stop "Mixing ratios validation failed"
    !    end if
    !end do
    !lumped_lambdas=>null() !< Initialize lumped kinetic rates pointer to null
    !> @section canonical_vectors Canonical Vector Extraction
    !> Extract canonical vectors (basis vectors) from mixing ratios matrix
    !> Canonical vectors represent independent mixing patterns that don't require reactive transport
    call this%transport%mixing_ratios_conc%get_can_vec(this%chemistry%CV_params%abs_tol,num_can_vec,ind_can_vec,num_non_can_vec,&
        ind_non_can_vec) !< Get indices of canonical and non-canonical vectors using tolerance
    !> @section init_conc Initialize Concentration History
    !> Set old concentrations in all target waters for time stepping
    do i=1,this%chemistry%num_waters_init !< Loop over all target waters (domain + boundary)
        !print *, this%chemistry%waters(i)%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
        call this%chemistry%waters_init(i)%set_conc_old() !< Copy current to old (time level k to k-1)
        call this%chemistry%waters(i)%set_conc_old() !< Copy current to old (time level k to k-1)
        call this%chemistry%waters_init(i)%set_conc_old_old() !< Copy old to old_old (k-1 to k-2)
        call this%chemistry%waters(i)%set_conc_old_old() !< Copy old to old_old (k-1 to k-2)
        call this%chemistry%waters_init(i)%solid_chemistry%set_conc_old() !< Copy solid concentrations to old
        call this%chemistry%waters(i)%solid_chemistry%set_conc_old() !< Copy solid concentrations to old
        call this%chemistry%waters_init(i)%solid_chemistry%set_conc_old_old() !< Copy solid concentrations to old_old
        call this%chemistry%waters(i)%set_solid_chemistry_old() !< Copy current kinetic rates to old
    end do
    !allocate(mix_waters(this%chemistry%num_waters))
    !> @section output_header Output headers/blocks are written per step by write_out_step_vars (honours chem_out_options%vars)
        !> Domain-integrated mass balance: reference cell (homogeneous reactive zone) + output file
        idx_ref=this%chemistry%tar_wat_indices(1)
        n_v_ref=this%chemistry%waters(idx_ref)%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
        allocate(mb_lhs(n_v_ref), mb_transp(n_v_ref), mb_react(n_v_ref))
        !> Per-cell [CELL-DBG] trace file (bulk calcite cells 3/8 + interface cells 15/16/17)
        open(newunit=iunit_celldbg, file=dir//root//'_Euler_consistent_cell_dbg.out', status='replace', action='write')
        !> @section time_loop Main Time Loop
        !> Loop over all time steps to solve reactive transport
        !> Write initial state (t=0) before the first time step, honouring chem_out_options%vars
        if (this%chemistry%chem_out_options%time_steps(1)==0) then
            call write_out_step_vars(this%chemistry, unit, 0d0)
        end if
        do k=1,this%transport%time_discr%Num_time !< Loop from time step 1 to final time step
            Delta_t=this%transport%time_discr%get_Delta_t(k) !< Get time step size for step k [s]
            time=time+Delta_t !< Advance simulation time [s]
            !> Reset the domain-integrated mass-balance accumulators for this time step
            mb_lhs=0d0
            mb_transp=0d0
            mb_react=0d0
            !> @subsection update_old Update Previous Time Step Values
            !> Shift time levels: current → old, old → old_old for multistep methods
            do i=1,this%chemistry%num_target_waters !< Loop over domain target waters only (not boundaries)
                call this%target_waters(i)%update_attributes() !< Shift all old attributes
                ! call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%update_conc_old() !< Shift aqueous concentrations
                ! call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%solid_chemistry%update_conc_old() !< Shift solid concentrations
                ! call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%update_rk_old() !< Shift aqueous kinetic rates
                ! call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%solid_chemistry%update_rk_old() !< Shift solid kinetic rates
                ! call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%update_solid_chemistry_old() !< Shift aqueous concentrations (old_old)
                ! call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%update_indices_rk_old() !< Shift solid concentrations (old_old)
            end do

            ! Build the target-water-indices-per-old-reaction-zone map used by the by-zone
            ! previous-time-step contribution to r_tilde (rebuilt each step from solid_chemistry_old).
            call this%chemistry%link_target_waters_old_reaction_zones(num_rz, wat_indices_zones)

            !call this%move_particles_stat_flux(k) !< Displace particles based on stationary flux field
            !call this%introduce_particle(k) !< Add new particles entering domain at time step k
            !> @subsection target_loop Target Waters Loop (Spatial Domain)
            !> Loop over all domain target waters to solve reactive transport at each location
            do i=1,this%chemistry%num_target_waters !< Loop over domain target waters (spatial cells)
                !> @subsubsection extract_dims Extract Problem Dimensions
                !> Get dimensions for current reactive zone
                n_p=this%chemistry%waters(&
                    this%chemistry%tar_wat_indices(i))%solid_chemistry%reactive_zone%speciation_alg%num_prim_species !< Number of primary species [-]
                n_v=this%chemistry%waters(&
                    this%chemistry%tar_wat_indices(i))%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species !< Total variable activity species [-]
                n_v_aq=this%chemistry%waters(&
                    this%chemistry%tar_wat_indices(i))%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species !< Aqueous variable activity species [-]
                n_v_aq_2=this%chemistry%waters(&
                    this%chemistry%tar_wat_indices(i))%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species !< Secondary aqueous variable activity species [-]
                
                !> @subsubsection check_exch Check for Cation Exchange Sites
                !> If exchange sites present, solver already assigned
                if (this%chemistry%waters(&
                    this%chemistry%tar_wat_indices(i))%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl &
                    >0) then !< Variable activity species include aqueous and solid (exchange sites)
                    continue !< Exchange sites present, no additional solver assignment needed
                end if
                
                !> @subsubsection allocate_conc Allocate Concentration Arrays
                allocate(c_hat(n_v)) !< Allocate solution vector for variable activity species [mol/L]
                !allocate(mix_waters(this%chemistry%num_target_waters)) !< Allocate solution vector for mixing waters [mol/L]
                allocate(conc_comp(n_p)) !< Allocate solution vector for components [mol/L]
                allocate(r_tilde(n_v)) !< Allocate solution vector for variable activity species [mol/L]
                !> Validate target water ordering (debugging check)
                if (this%transport%mix_conc_indices%cols(i)%col_1(1)/= &
                    this%chemistry%tar_wat_indices(i)) then
                        print *, this%transport%mix_conc_indices%cols(i)%col_1,&
                            this%chemistry%tar_wat_indices(i)
                        print *, i
                        error stop "Target waters not in the right order"
                    end if
                
                !> @subsubsection transport_mix Transport Mixing Step
                !> Compute mixed concentrations after advective-dispersive transport
                ! Page-4 two-stage per-reactive-zone consistent mix (within-zone speciation + cross-zone
                ! aqueous mixing); externals folded in stage 3. Uses conc_old (previous-step Eulerian state).
                ! The reaction mixing ratios share the concentration stencil (0 for the external waters).
                call p_mixer(this%chemistry%waters, this%chemistry%zone_waters(i,:), &
                    this%transport%mix_conc_indices%cols(i)%col_1(1:this%transport%mix_conc_indices%cols(i)%dim-2), &
                    this%transport%mixing_ratios_conc%cols(i)%col_1(1:this%transport%mix_conc_indices%cols(i)%dim-2), &
                    this%transport%mixing_ratios_R%cols(i)%col_1(1:this%transport%mix_conc_indices%cols(i)%dim-2), &
                    this%chemistry%wat_indices_zones, this%chemistry%ext_waters_indices, Delta_t, &
                    this%transport%time_discr%theta_r, c_mix, c_tilde, .true.)
                call p_r_tilde(this%chemistry%waters, &
                    this%transport%mix_conc_indices%cols(i)%col_1(&
                    2:this%transport%mix_conc_indices%cols(i)%dim-2),& !< whole waters array + mixing-water index list (no vector-subscript gather)
                    this%chemistry%waters(this%chemistry%tar_wat_indices(i)&
                    )%solid_chemistry%reactive_zone%ind_var_act_species,&
                    this%transport%mixing_ratios_R%cols(i)%col_1(&
                    2:),&
                    this%transport%mix_conc_indices%cols(i)%col_1(&
                    this%transport%mix_conc_indices%cols(i)%dim-1),&
                    this%transport%mix_conc_indices%cols(i)%col_1(&
                    this%transport%mix_conc_indices%cols(i)%dim),&
                    this%transport%time_discr%theta_r,r_tilde, & !< Time weighting factor [-]
                    curr_water_idx=this%chemistry%tar_wat_indices(i), & !< current (self) cell index
                    lambda_self=this%transport%mixing_ratios_R%cols(i)%col_1(1), & !< self reaction ratio (Eulerian by-reactive-zone previous-step term)
                    wat_indices_zones=wat_indices_zones) !< target-water indices grouped by old reactive zone
                !> Positivity-preserving reactive update (classical Patankar): keeps c_hat >= 0
                !> WITHOUT modifying the reaction mixing ratios (dispersion unchanged).
                call compute_c_hat_patankar(c_mix,Delta_t,r_tilde,c_hat) !< positivity-preserving mixed+reacted concentrations
                call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%set_conc_hat(c_hat) !< aqueous blocks of the Newton speciation input, for the conc_hat attribute
                conc_nc=this%chemistry%waters(this%chemistry%tar_wat_indices(i))%get_conc_var_act() !< Reorder mixed concentrations to original ordering
                !> @subsubsection reactive_mix Reactive Mixing Iteration
                !> Solve coupled transport-reaction system using selected solver
                !print *, this%chemistry%tar_wat_indices(i)
                !> Autentica chapuza
                if (this%chemistry%r_down_opt==4) then
                    lambda_R=this%transport%compute_mix_ratio_R_opt4(i) !< Precompute modified reaction lambdas
                else
                    !allocate(lambdas_R(this%transport%mixing_ratios_R%num_cols)) !< Allocate lambdas_R array
                    !do m=1,this%transport%mixing_ratios_R%num_cols
                        !lambdas_R(m)=this%transport%mixing_ratios_R%cols(m)%col_1(1) !< Get lambdas from initial reaction mixing ratios
                    !end do
                    lambda_R=this%transport%mixing_ratios_R%cols(i)%col_1(1) !< Get lambda for current target water
                    !lambdas_R=lambdas_R*theta !< Scale by time weighting factor
                end if
                call p_solver(this%chemistry%waters(this%chemistry%tar_wat_indices(i)),& !< Current target water
                    this%chemistry%waters(this%chemistry%tar_wat_indices(i))%get_conc_old_old_prim_new_rz(),& !< Primary species from k-1 [mol/L]
                    c_hat,& !< Mixed concentrations from transport
                    lambda_R,& !< Reaction mixing ratio for current target water [-]
                    Delta_t,& !< Time step size
                    this%transport%time_discr%theta_r,& !< Time weighting factor [-]
                    conc_nc,&
                    conc_comp) !< Output: new concentrations
                
                !> @subsubsection update_derived Update Derived Quantities
                !> Compute pH, salinity, and ionic strength from updated concentrations
                call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%compute_pH()
                call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%compute_saturation_mins() !< saturation indices of all minerals at the end-of-step state
                call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%compute_salinity()
                call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%compute_ionic_strength()
                
                !> @subsubsection accum_rk Accumulate Reaction Amounts
                !> Update cumulative aqueous reaction amounts for mass balance tracking
                this%chemistry%waters(this%chemistry%tar_wat_indices(i))%cum_rxn_amnt_kin=&
                    this%chemistry%waters(this%chemistry%tar_wat_indices(i))%cum_rxn_amnt_kin+&
                    this%chemistry%waters(this%chemistry%tar_wat_indices(i))%rxn_amnt_kin !< Accumulate aqueous reaction amounts [mol]
                !> Update cumulative solid reaction amounts for mass balance tracking
                this%chemistry%waters(this%chemistry%tar_wat_indices(i))%solid_chemistry%cum_rxn_amnt_kin=&
                    this%chemistry%waters(this%chemistry%tar_wat_indices(i))%solid_chemistry%cum_rxn_amnt_kin+&
                    this%chemistry%waters(this%chemistry%tar_wat_indices(i))%solid_chemistry%rxn_amnt_kin !< Accumulate solid reaction amounts [mol]
                
                !> @subsubsection eq_rates Equilibrium Reaction Rates
                !> Compute equilibrium reaction rates from mass balance equation
                if (this%chemistry%waters(&
                    this%chemistry%tar_wat_indices(i))%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions>0) &
                        then !< Equilibrium reactions present (kinetic term handled inside, zero if none)
                    call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%compute_Re(&
                        c_hat(1:n_v),& !< All variable activity species concentrations
                        Delta_t,& !< Time step size
                        lambda_R,& !< Reaction mixing ratio [-]
                        this%transport%time_discr%theta_r) !< Time weighting factor [-]
                end if

                !> [CELL-DBG] intermediate results for bulk calcite cells 3/8 and interface cells 15/16/17
                if (k<=10 .and. (this%chemistry%num_target_waters<=10 .or. i==3 .or. i==8 .or. i==15 .or. i==16 .or. i==17)) then
                    call write_celldbg_cell(this%chemistry%waters(this%chemistry%tar_wat_indices(i)), iunit_celldbg, k, i, n_v, c_mix, c_tilde, c_hat)
                end if

                !> Domain-integrated mass-balance contribution of this cell
                call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%get_react_mix_mass_bal_terms( &
                    c_mix(1:n_v), c_hat(1:n_v), lambda_R, Delta_t, this%transport%time_discr%theta_r, &
                    mb_lhs_c, mb_transp_c, mb_react_c)
                if (size(mb_lhs_c) == n_v_ref) then
                    mb_lhs = mb_lhs + mb_lhs_c
                    mb_transp = mb_transp + mb_transp_c
                    mb_react = mb_react + mb_react_c
                end if
                deallocate(mb_lhs_c, mb_transp_c, mb_react_c)

                !> @subsubsection solid_chem Solid Chemistry State Variablestate Variables
                !> Update mineral state if minerals present in reactive zone
                if (associated(this%chemistry%waters(&
                    this%chemistry%tar_wat_indices(i))%solid_chemistry%mineral_zone)) then
                    !> Compute mass volumetric fractions of minerals from mass balance equation
                    call this%chemistry%waters(&
                        this%chemistry%tar_wat_indices(i))%solid_chemistry%compute_mass_bal_mins(&
                        Delta_t) !< Time step size [s]
                    !> Compute concentrations of minerals from volumetric fractions
                    call this%chemistry%waters(&
                        this%chemistry%tar_wat_indices(i))%solid_chemistry%compute_conc_minerals_iter(&
                        Delta_t,lambda_R) !< Time step size [s], reaction mixing ratio
                end if
                
                !> @subsubsection gas_species_chem Gas Chemistry State Variables
                !> Update gas state if gas phase present
                if (associated(this%chemistry%waters(&
                    this%chemistry%tar_wat_indices(i))%gas_chemistry)) then
                    !> Compute concentrations of gases from reaction rates
                    call this%chemistry%waters(&
                        this%chemistry%tar_wat_indices(i))%gas_chemistry%compute_conc_gases_iter(&
                        Delta_t,& !< Time step size [s]
                        this%chemistry%waters(this%chemistry%tar_wat_indices(i))%volume,& !< Water volume [L]
                        [this%chemistry%waters(this%chemistry%tar_wat_indices(i))%rxn_rate_eq,& !< Mean eq. rates [mol/L/s]
                        this%chemistry%waters(this%chemistry%tar_wat_indices(i))%rxn_rate_kin_avg]) !< Mean kin. rates [mol/L/s]
                    !> Compute volume of gas phase from concentrations
                    call this%chemistry%waters(&
                        this%chemistry%tar_wat_indices(i))%gas_chemistry%compute_vol_gas_species_conc()
                    !> Compute activity coefficients of gases (fugacity coefficients)
                    call this%chemistry%waters(&
                        this%chemistry%tar_wat_indices(i))%gas_chemistry%compute_log_act_coeffs_gases()
                end if
                !> @subsubsection cleanup_dealloc Cleanup and Deallocation
                deallocate(c_hat) !< Free transport mixing array
                deallocate(c_mix) !< Free mixed concentrations array
                deallocate(c_tilde) !< Free conservative mixed concentrations array
                deallocate(conc_comp) !< Free concentration array
                deallocate(conc_nc) !< Free upstream kinetic rates array
                deallocate(r_tilde) !< Free reaction rates array
                !deallocate(lambdas_R) !< Free mixing waters array
                !deallocate(mix_waters) !< Free mixing waters array
            end do !< End target waters spatial loop
            !> @subsection particle_tracking Lagrangian Particle Tracking
            ! if (k<this%transport%time_discr%Num_time) then !< Skip particle movement on final time step
            !     !> Move particles according to velocity field (Lagrangian step)
            !     call this%move_particles_stat_flux(k) !< Displace particles based on stationary flux field
            !     !> Introduce new particles at inlet boundary
            !     call this%introduce_particle(k) !< Add new particles entering domain at time step k
            !     !> chapuza
            !     ! this%chemistry%waters(this%chemistry%bd_waters_indices(1)&
            !     !     )%concentrations=2d0*this%chemistry%wat_types(1)%concentrations - &
            !     !     this%chemistry%waters(this%chemistry%tar_wat_indices(1))%concentrations !< Update inlet boundary concentrations
            !     print *, "Inlet concentrations updated at time step ", k
            !     print *, this%chemistry%waters(this%chemistry%bd_waters_indices(1))%concentrations(2:)
            ! end if
            ! print *, this%chemistry%waters(1)%concentrations(2:) !< Debugging print statement
            ! print *, this%chemistry%waters(1)%pos !< Debugging print statement
            ! print *, this%chemistry%waters(1)%solid_chemistry%tar%id !< Debugging print statement
            ! print *, this%chemistry%waters(1)%indices_rk%cols(1)%col_1 !< Debugging print statement
            ! print *, this%chemistry%waters(1)%get_Sk_nc() !< Debugging print statement
            ! print *, this%chemistry%waters(2)%concentrations(2:) !< Debugging print statement
            ! print *, this%chemistry%waters(2)%pos !< Debugging print statement
            ! print *, this%chemistry%waters(2)%solid_chemistry%tar%id !< Debugging print statement
            ! print *, this%chemistry%waters(2)%indices_rk%cols(1)%col_1 !< Debugging print statement
            ! print *, this%chemistry%waters(2)%get_Sk_nc() !< Debugging print statement
            ! print *, this%chemistry%waters(3)%concentrations(2:) !< Debugging print statement
            ! print *, this%chemistry%waters(3)%pos !< Debugging print statement
            ! print *, this%chemistry%waters(3)%solid_chemistry%tar%id !< Debugging print statement
            ! print *, this%chemistry%waters(3)%indices_rk%cols(1)%col_1 !< Debugging print statement
            ! print *, this%chemistry%waters(3)%get_Sk_nc() !< Debugging print statement
            ! print *, this%chemistry%waters(16)%concentrations(2:) !< Debugging print statement
            ! print *, this%chemistry%waters(16)%pos !< Debugging print statement
            ! print *, this%chemistry%waters(16)%solid_chemistry%tar%id !< Debugging print statement
            ! print *, this%chemistry%waters(16)%indices_rk%cols(1)%col_1 !< Debugging print statement
            ! print *, this%chemistry%waters(16)%get_Sk_nc() !< Debugging print statement
            ! print *, this%chemistry%waters(17)%concentrations(2:) !< Debugging print statement
            ! print *, this%chemistry%waters(17)%pos !< Debugging print statement
            ! print *, this%chemistry%waters(17)%solid_chemistry%tar%id !< Debugging print statement
            ! print *, this%chemistry%waters(17)%indices_rk%cols(1)%col_1 !< Debugging print statement
            ! print *, this%chemistry%waters(17)%indices_rk%cols(2)%col_1 !< Debugging print statement
            ! print *, this%chemistry%waters(17)%get_Sk_nc() !< Debugging print statement
            ! print *, this%chemistry%waters(18)%concentrations(2:) !< Debugging print statement
            ! print *, this%chemistry%waters(18)%pos !< Debugging print statement
            ! print *, this%chemistry%waters(18)%solid_chemistry%tar%id !< Debugging print statement
            ! print *, this%chemistry%waters(18)%indices_rk%cols(1)%col_1 !< Debugging print statement
            ! print *, this%chemistry%waters(18)%get_Sk_nc() !< Debugging print statement
            ! print *, this%chemistry%waters(this%chemistry%num_waters-1)%concentrations(2:) !< Debugging print statement
            ! print *, this%chemistry%waters(this%chemistry%num_waters-1)%pos !< Debugging print statement
            ! print *, this%chemistry%waters(this%chemistry%num_waters-1)%solid_chemistry%tar%id !< Debugging print statement
            ! print *, this%chemistry%waters(this%chemistry%num_waters-1)%indices_rk%cols(1)%col_1 !< Debugging print statement
            ! print *, this%chemistry%waters(this%chemistry%num_waters-1)%get_Sk_nc() !< Debugging print statement
            ! print *, this%chemistry%waters(this%chemistry%num_waters)%concentrations(2:) !< Debugging print statement
            ! print *, this%chemistry%waters(this%chemistry%num_waters)%pos !< Debugging print statement
            ! print *, this%chemistry%waters(this%chemistry%num_waters)%solid_chemistry%tar%id !< Debugging print statement
            ! print *, this%chemistry%waters(this%chemistry%num_waters)%indices_rk%cols(1)%col_1 !< Debugging print statement
            ! print *, this%chemistry%waters(this%chemistry%num_waters)%get_Sk_nc() !< Debugging print statement
            !> Debugging: Validate mass balance after time step
            !> @subsection time_update Update Time Step History
            !> Shift time levels for next time step: k+1 → k, k → k-1
            !waters_old_old=waters_old !< Store k for k-1 (previous previous)
            waters_old=this%chemistry%waters !< Store k+1 for k (previous)
            !> Domain-integrated mass balance for this time step (per variable-activity species)
            call this%chemistry%write_react_mix_mass_bal_dom(idx_ref,time,n_v_ref,mb_lhs,mb_transp,mb_react)
            !> Write chemistry output honouring chem_out_options%vars (conc / react_rate / vol_fract)
            if (k==this%chemistry%chem_out_options%time_steps(kk)) then
                call write_out_step_vars(this%chemistry, unit, time)
                if (kk<this%chemistry%chem_out_options%num_time_steps) kk=kk+1
            end if
        end do !< End main time loop
    
    !> Close output file
    close(unit)
    close(iunit_celldbg)
    deallocate(mb_lhs, mb_transp, mb_react)
    end select !< End type selection
    !deallocate(lambdas_R) !< Free lumped kinetic rates array
    !deallocate(mix_waters) !< Free mixing waters array
    deallocate(waters_old) !< Free old target waters array
end subroutine !< End solve_RT_ideal_cons_Euler_stat_flux_1D