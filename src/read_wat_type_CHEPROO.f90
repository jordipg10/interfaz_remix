!> @file read_wat_type_CHEPROO.f90
!> @brief Reads a water type from CHEPROO data input format
!>
!> @details
!> This subroutine reads and initializes a single water type from CHEPROO-formatted input data.
!> CHEPROO is a geochemical modeling framework that describes aqueous chemistry with various
!> constraint types (icons) for primary species.
!>
!> Key functionality:
!> - Parses primary aqueous species with 4 icon types (constraint methods)
!> - Handles equilibrium constraints with minerals, gases, and surface species
!> - Computes initial concentrations through speciation calculations
!> - Supports multiple activity coefficient models and solving methods
!> - Integrates solid chemistry and gas chemistry
!>
!> Icon types (constraint methods):
!> - icon=1: Total concentration constraint (ctot specified)
!> - icon=2: Component constraint (for mass balance)
!> - icon=3: Direct concentration/activity constraint (e.g., pH)
!> - icon=4: Equilibrium with mineral, gas, or surface species
!>
!> Input file format:
!> ```
!> 'guess'                                    ! Header line
!> <species_name> <icon> <guess> <ctot> <constraint_name>
!> <species_name> <icon> <guess> <ctot> <constraint_name>
!> ...
!> '*'                                        ! Terminator
!> ```
!>
!> Concentration units assumption: **All input concentrations are in molalities**
!>
!> @param[in,out] this Aqueous chemistry object to initialize
!> @param[in] n_p_aq Number of primary aqueous species
!> @param[in] num_cstr Number of constraints (icon=4 entries)
!> @param[in] num_gas_zones Number of gas zones in the system
!> @param[in] model Activity coefficients model (0=ideal, 1+=non-ideal)
!> @param[in] Jac_opt Jacobian option: 0=incremental coefficients, 1=analytical
!> @param[in] unit File unit number for reading
!> @param[out] niter Number of iterations performed in speciation calculation
!> @param[out] CV_flag Convergence flag (TRUE if converged, FALSE otherwise)
!> @param[in] gas_species_chem Gas chemistry object (optional) - chapuza (temporary fix)
!>
!> @note Assumes concentration units in input file are molalities
!> @note Currently uses several "chapuza" (temporary fixes) for specific cases
!> @warning Surface chemistry and cation exchange implementation incomplete
!> @warning Multiple gas zones not yet implemented (num_gas_zones must be 0 or 1)
!>
!> @see initialise_conc_anal_ideal
!> @see initialise_conc_anal
!> @see initialise_conc_anal_exch
!> @see compute_c2_from_c1_Picard
!> @see compute_conc_surf_ideal
!>
!> Algorithm overview:
!> 1. Validate input parameters and allocate arrays
!> 2. Read species data line-by-line (icon, guess, ctot, constraint)
!> 3. Process icon types and update counters
!> 4. Handle equilibrium constraints (minerals, gases, surface species)
!> 5. Set up speciation algebra and reactive zone
!> 6. Compute initial concentrations based on icon distribution and model
!> 7. Handle species index swapping if needed (flag_Se)
!>
!> @author jordi Prat
!> @date 2025
!>
subroutine read_wat_type_CHEPROO(this,n_p_aq,num_cstr,num_gas_zones,model,Jac_opt,unit,niter,CV_flag,gas_chem)
    !> ================================================================
    !> Module imports - bring in required chemistry types
    !> ================================================================
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use gas_chemistry_m, only: gas_chemistry_c
    use solid_chemistry_m, only: solid_chemistry_c
    use aq_species_m, only: aq_species_c
    use species_m, only: species_c
    implicit none !> enforce explicit variable declarations
    
    !> ================================================================
    !> Input/output parameters
    !> ================================================================
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object to initialize (modified)
    integer(kind=4), intent(in) :: n_p_aq !> number of primary aqueous species
    integer(kind=4), intent(in) :: num_cstr !> number of constraints (icon=4 entries)
    integer(kind=4), intent(in) :: num_gas_zones !> number of gas zones in system
    integer(kind=4), intent(in) :: model !> activity coefficients model (0=ideal, 1+=non-ideal)
    integer(kind=4), intent(in) :: Jac_opt !> Jacobian option: 0=incremental coefficients, 1=analytical
    integer(kind=4), intent(in) :: unit !> file unit number for reading
    integer(kind=4), intent(out) :: niter !> number of iterations in speciation calculation (output)
    logical, intent(out) :: CV_flag !> convergence flag: TRUE if converged, FALSE otherwise (output)
    type(gas_chemistry_c), intent(in), optional :: gas_chem !> gas chemistry object (optional) - chapuza (temporary fix)
    !type(solid_chemistry_c), intent(in), optional :: surf_chem !> (COMMENTED) surface chemistry for exchange reactions
    !real(kind=8), intent(in), optional :: c1_surf !> (COMMENTED) surface species concentration - chapuza
    !real(kind=8), intent(in), optional :: CEC !> (COMMENTED) cation exchange capacity - chapuza
    !real(kind=8), intent(out), optional :: conc_exch(:) !> (COMMENTED) exchange cation concentrations - chapuza
    
    !> ================================================================
    !> Local variables - loop counters and indices
    !> ================================================================
    integer(kind=4) :: i,j,k,l,ll,m,niwtype,nbwtype,nrwtype,gas_ind,min_ind,n_comp_aq,ind_cstr,n_gas_species_constr,n_aq_comp,icon,&
        aq_sp_ind,ind_sp,gas_flag,n !> integer variables
        !> i: general loop counter
        !> j: mineral equilibrium counter
        !> k: primary aqueous species counter
        !> l: constraints counter (total)
        !> ll: heterogeneous constraints counter
        !> m: gas constraints counter
        !> n: equilibrium reactions counter
        !> niwtype, public, nbwtype, public, nrwtype: water type counters (not currently used)
        !> gas_ind: gas index in gas phase
        !> min_ind: mineral index (not currently used)
        !> n_comp_aq: number of aqueous components (not currently used)
        !> ind_cstr: index of constraint in chemical system
        !> n_gas_species_constr: number of gas constraints (not currently used)
        !> n_aq_comp: number of aqueous components (not currently used)
        !> icon: constraint type for current species (1-4)
        !> aq_sp_ind: aqueous species index in aqueous phase
        !> ind_sp: species index in chemical system
        !> gas_flag: flag indicating if gas phase is present (0=no, 1=yes)
    
    !> Allocatable integer arrays
    integer(kind=4), allocatable :: icons(:) !> icon types for each primary species (1-4)
    integer(kind=4), allocatable :: indices_constrains(:,:) !> constraint indices: (l,1)=species index, (l,2)=eq reaction index
    integer(kind=4), allocatable :: gas_indices(:) !> gas indices (not currently used)
    integer(kind=4), allocatable :: n_icons(:) !> count of each icon type (dimension 4)
    integer(kind=4), allocatable :: prim_indices(:) !> primary species indices (not currently used)
    integer(kind=4), allocatable :: swap(:) !> indices of species to swap for speciation algebra (dimension 2)
    integer(kind=4), allocatable :: aux_ind_cstr(:,:) !> auxiliary constraint indices for swapping
    integer(kind=4), allocatable :: ind_mins(:) !> mineral indices in chemical system
    integer(kind=4), allocatable :: ind_swap(:) !> indices in constraint array to swap (dimension 2)
    integer(kind=4), allocatable :: ind_non_flow_species(:) !> indices for non-flowing species
    integer(kind=4), allocatable :: ind_eq_reacts(:) !> indices for heterogeneous equilibrium reactions
    
    !> Allocatable real arrays
    real(kind=8), allocatable :: ctots(:) !> total concentrations for each primary species
    real(kind=8), allocatable :: c2_init(:) !> initial secondary species concentrations
    real(kind=8), allocatable :: c2_ig(:) !> initial guess for secondary species concentrations
    real(kind=8), allocatable :: a1(:) !> activities of primary species
    real(kind=8), allocatable :: c1(:) !> concentrations of primary species (molarities)
    real(kind=8), allocatable :: c1_aq(:) !> aqueous concentrations of primary species (molarities)
    real(kind=8), allocatable :: act_ads_cats_ig(:) !> initial guess for adsorbed cation activities
    real(kind=8), allocatable :: log_act_coeffs(:) !> log10 activity coefficients of secondary species
    
    !> Character variables for parsing
    character(len=256) :: prim_sp_name !> primary species name (not currently used)
    character(len=256) :: label !> section label from file
    character(len=256) :: aq_sp_name !> aqueous species name (not currently used)
    character(len=256), allocatable :: constrains(:) !> constraint names (not currently used)
    
    !> Real scalar variables
    real(kind=8) :: temp !> temperature (not currently used)
    real(kind=8) :: conc !> concentration variable (not currently used)
    real(kind=8) :: ionic_strength !> ionic strength (not currently used)
    real(kind=8) :: guess !> initial guess for species concentration
    real(kind=8) :: ctot !> total concentration for constraint
    
    !> Logical flags
    logical :: flag_gas !> TRUE if species is a gas (not currently used)
    logical :: flag_min !> TRUE if species is a mineral (not currently used)
    logical :: flag !> general flag for presence checks
    logical :: flag_comp !> TRUE to use component matrix with constant activity species
    logical :: flag_surf !> TRUE if surface chemistry (cation exchange) is present
    logical :: flag_Se !> TRUE if species swap is needed for speciation algebra
    logical :: flag_sp !> TRUE if species found in chemical system
    
    !> ================================================================
    !> Temporary chemistry objects
    !> ================================================================
    !type(reactive_zone_c) :: react_zone !> temporary reactive zone object (default)
    !type(gas_chemistry_c) :: gas_species_chem !> (COMMENTED) gas chemistry object - now optional parameter
    type(solid_chemistry_c) :: solid_chem !> temporary solid chemistry object
    !type(gas_species_c) :: gas !> temporary gas object (not currently used)
    !type(mineral_c) :: mineral !> temporary mineral object (not currently used)
    !type(mineral_zone_c) :: min_zone !> default mineral zone object
    type(aq_species_c) :: aq_species !> aqueous species object for parsing
    type(species_c) :: constrain !> constraint species object for parsing
    
    !> ================================================================
    !> SECTION 1: Surface chemistry initialization (reactive zone)
    !> ================================================================
    !print *, associated(this%gas_chemistry) !> (COMMENTED DEBUG) check if gas chemistry is associated

    !> Check if this water type has surface chemistry (cation exchange)
    if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then !> if surface complexes exist
        !> (COMMENTED) Original approach - set surface chemistry from reactive zones
        !call surf_chem%set_reactive_zone(this%reactive_zones(1)) !> chapuza, we assume that the first reactive zone is the one with surface chemistry
        !call this%set_solid_chemistry(surf_chem)
    else !> no surface complexes - use default setup
        !> (COMMENTED) Alternative approach - set defaults for solid chemistry
        !call solid_chem%set_reactive_zone(react_zone) !> we set reactive zone by default
        !call solid_chem%set_mineral_zone(min_zone) !> we set mineral zone by default
        !call this%set_solid_chemistry(solid_chem) !> we set solid chemistry by default
        !> Allocate array for non-flowing species (minerals, gases, surface species)
        !call this%solid_chemistry%reactive_zone%allocate_ind_non_flow_species(num_cstr) !> allocate with size = number of constraints
         allocate(ind_non_flow_species(num_cstr)) !> allocate with size = number of constraints
    end if
    
    !> (COMMENTED) Alternative mineral zone setup
    !if (this%solid_chemistry%reactive_zone%chem_syst%num_minerals>0) then !> if minerals exist in chemical system
        ! call min_zone%set_chem_syst_min_zone(this%solid_chemistry%reactive_zone%chem_syst) !> link chem system to mineral zone
        ! call this%solid_chemistry%set_mineral_zone(min_zone) !> link mineral zone to solid chemistry
    !end if
    
    allocate(ind_mins(num_cstr)) !> allocate mineral indices array - chapuza (temporary array, size = constraints)
    allocate(ind_eq_reacts(this%solid_chemistry%reactive_zone%chem_syst%speciation_alg%num_eq_reactions)) !> allocate indices equilibrium reactions array
    forall (i=1:this%solid_chemistry%reactive_zone%chem_syst%speciation_alg%num_eq_reactions) ind_eq_reacts(i)=i !> initialize to default indices (chapuza)
    !> ================================================================
    !> SECTION 2: Validation and array allocation
    !> ================================================================
    !> (COMMENTED) Alternative file reading structure
    !    read(unit,*) label !> read section label
    !    if (label=='INITIAL AND BOUNDARY WATER TYPES') then !> if reading water types section
    !        i=0 !> counter initial water types
    !        read(unit,*) model !> read activity coefficients model
    !        read(unit,*) this%num_init_wat_types, this%num_bd_wat_types, this%num_rech_wat_types !> read water type counts
        
        !> Validate number of primary aqueous species
        if (n_p_aq<0 .or. n_p_aq>this%aq_phase%num_species) then !> if n_p_aq out of valid range
            error stop "Number of primary aqueous species not valid" !> stop with error message
        !> Validate number of constraints  
        else if (num_cstr<0 .or. num_cstr>n_p_aq) then !> if num_cstr out of valid range (0 to n_p_aq)
            error stop "Number of constrains not valid" !> stop with error message
        else !> parameters valid - proceed with allocation
            !> Allocate arrays for constraint data
            allocate(prim_indices(n_p_aq),icons(n_p_aq),ctots(n_p_aq),indices_constrains(num_cstr,2))
                !> prim_indices: primary species indices array
                !> icons: icon types for each primary species (1-4)
                !> ctots: total concentrations for each primary species
                !> indices_constrains: (num_cstr, 2) - (species index, eq reaction index)
            
            !> (COMMENTED) Alternative gas chemistry initialization
            !if (this%solid_chemistry%reactive_zone%chem_syst%gas_phase%num_species>0) then !> if gas species exist
            !    call react_zone%set_gas_phase(this%solid_chemistry%reactive_zone%chem_syst%gas_phase) !> link gas phase
            !    call gas_species_chem%set_reactive_zone(react_zone) !> link reactive zone to gas chemistry
            !    call gas_species_chem%allocate_partial_pressures() !> allocate partial pressures array
            !    call gas_species_chem%allocate_conc_gases() !> allocate gas concentrations array
            !    call gas_species_chem%allocate_log_act_coeffs_gases() !> allocate gas activity coefficients array
            !    call gas_species_chem%set_temp(this%temp) !> set temperature for gas chemistry
            !    call gas_species_chem%set_volume(1d0) !> set arbitrary volume (1.0)
            !    call this%set_gas_chemistry(gas_species_chem) !> link gas chemistry to aqueous chemistry
            !end if
        end if
        
        !> (COMMENTED) Alternative concentration and activity arrays allocation
        !call this%allocate_conc_comp(n_p_aq) !> allocate component concentrations array
        !call this%allocate_log_act_coeffs() !> allocate log activity coefficients array
        !call this%set_log_act_coeffs() !> initialize log activity coefficients
        !call this%allocate_activities_aq_species() !> allocate aqueous species activities array
        
        !> ================================================================
        !> SECTION 3: Read icon-based constraints from file
        !> ================================================================
        read(unit,*) label !> read section label from file
        if (index(label,'guess')/=0) then !> check if label contains 'guess' substring !> file format: 'icon, guess, ctot, constrain'
            !> Initialize counters for different constraint types
            k=1 !> counter for primary aqueous species
            l=0 !> counter for total constraints
            ll=0 !> counter for heterogeneous (non-aqueous) constraints
            m=1 !> counter for gas constraints
            j=1 !> counter for mineral constraints
            n=this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts !> counter for equilibrium reactions in reactive zone
            allocate(n_icons(4)) !> allocate icon count array (one element per icon type 1-4)
            n_icons=0 !> initialize all icon counts to zero !> number of each icon option
            gas_flag=0 !> initialize gas flag - chapuza
            
            !> ================================================================
            !> Main loop: read icon-based constraints for each primary species
            !> ================================================================
            do !> infinite loop - exit when '*' found
                !> Read constraint data for current primary species
                !> Read constraint data for current primary species
                read(unit,*) aq_species%name, icon, guess, ctot, constrain%name
                    !> aq_species%name: aqueous species name
                    !> icon: constraint type (1-4)
                    !> guess: initial guess for concentration/activity
                    !> ctot: total concentration or constraint value
                    !> constrain%name: constraint species name (mineral, gas, or surface species)
                
                !> Check exit condition: '*' marks end of species list
                if (aq_species%name=='*') then !> if sentinel character found
                    exit !> exit reading loop
                else !> continue processing species
                    !> (COMMENTED) Alternative aqueous phase check
                    !call this%aq_phase%is_species_in_aq_phase(aq_species,flag,aq_sp_ind) !> check if species in aqueous phase
                    !if (flag.eqv..true.) then !> if species found
                    !prim_indices(k)=this%indices_aq_species(k) !> store species index
                    
                    !> (COMMENTED) Chapuza - special handling for water
                    !if (aq_species%name=='h2o(p)') then !> if species is water
                    !    guess=1d0/18d-3 !> set guess = 1/(18 g/mol) = molarity of pure water
                    !    ctot=guess !> impose water concentration constraint
                    !end if
                    
                        !> ================================================================
                        !> Process constraints based on icon type (1-4)
                        !> ================================================================
                        if (icon==1) then !> ICON 1: Total concentration constraint
                            n_icons(1)=n_icons(1)+1 !> increment icon 1 counter
                            !> (COMMENTED) Special water handling for icon 1
                            !if (aq_species%name=='h2o(p)') then !> if species is water
                            !    guess=1d0/18d-3 !> set guess for water molarity
                            !    ctot=guess !> impose water concentration constraint
                            !end if
                            
                        else if (icon==2) then !> ICON 2: Component constraint
                            n_icons(2)=n_icons(2)+1 !> increment icon 2 counter
                            n_aq_comp=n_aq_comp+1 !> increment aqueous component counter
                            
                        else if (icon==3) then !> ICON 3: Direct concentration/activity constraint
                            n_icons(3)=n_icons(3)+1 !> increment icon 3 counter
                            if (aq_species%name=='h+') then !> if species is H+ (proton)
                                !> (COMMENTED) Alternative proton index setting
                                !call this%aq_phase%set_ind_prot(this%indices_aq_species(k)) !> set proton index in aqueous phase
                                call this%set_pH(-log10(ctot)) !> set pH from concentration (pH = -log10[H+])
                            ! else if (aq_species%name=='h2o(p)') then !> (COMMENTED) if species is water
                            !     call this%compute_conc_ideal_water() !> compute ideal water concentration
                            end if
                            
                        else if (icon==4) then !> ICON 4: Equilibrium constraint (mineral, gas, or surface)
                            n_icons(4)=n_icons(4)+1 !> increment icon 4 counter
                            !> Check if constraint is an equilibrium reaction in chemical system
                            call this%solid_chemistry%reactive_zone%chem_syst%is_eq_reaction_in_chem_syst(constrain%name,flag,&
                                ind_cstr) !> search for equilibrium reaction, get flag and index
                            if (flag.eqv..true.) then !> if constraint found in chemical system
                                l=l+1 !> increment total constraint counter
                                !n=n+1 !> increment equilibrium reactions counter
                                indices_constrains(l,1)=k !> store species index in water type
                                indices_constrains(l,2)=ind_cstr !> store equilibrium reaction index in chem system
                                !ind_eq_reacts(n)=ind_cstr !> store index of heterogeneous equilibrium reaction
                                !> (COMMENTED) Alternative index adjustment
                                !ind_cstr=ind_cstr-this%solid_chemistry%reactive_zone%chem_syst%num_redox_eq_reacts-&
                                !this%solid_chemistry%reactive_zone%chem_syst%aq_phase%num_aq_complexes !> chapuza
                                
                                !> Get constraint species index in chemical system
                                call this%solid_chemistry%reactive_zone%chem_syst%is_species_in_chem_syst(constrain,flag_sp,ind_sp)
                                if (flag_sp.eqv..true.) then !> if constraint species found in chemical system
                                    !> (COMMENTED) Alternative heterogeneous constraint counter
                                    !ll=ll+1 !> increment counter of non-aqueous constraints
                                    !ind_cstr=ind_cstr-this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts !> adjust index to heterog eq reaction
                                    
                                    !> Check if this is a heterogeneous (non-aqueous) equilibrium reaction
                                    if (ind_cstr>this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts) then !> if heterogeneous
                                        ll=ll+1 !> increment heterogeneous equilibrium reactions counter
                                        !> Assign species index to non-flowing species array
                                        ind_non_flow_species(ll)=ind_sp
                                        n=n+1 !> increment equilibrium reactions counter
                                        ind_eq_reacts(n)=ind_cstr !> store index of heterogeneous equilibrium reaction
                                        !> Adjust index to heterogeneous equilibrium reaction (subtract aqueous reactions)
                                        ind_cstr=ind_cstr-this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts
                                        
                                        !> --------------------------------------------------------
                                        !> Determine type of heterogeneous constraint
                                        !> --------------------------------------------------------
                                        if (ind_cstr<=this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq_var_act) then
                                            !> Mineral with variable activity (dissolving/precipitating)
                                            ind_mins(j)=ind_cstr !> store mineral index - chapuza
                                            j=j+1 !> increment mineral counter
                                            !> Update mineral counters in reactive zone
                                            this%solid_chemistry%reactive_zone%num_minerals_var_act=&
                                                this%solid_chemistry%reactive_zone%num_minerals_var_act+1
                                            this%solid_chemistry%reactive_zone%num_minerals=&
                                                this%solid_chemistry%reactive_zone%num_minerals+1
                                            !> Update mineral counters in mineral zone
                                            this%solid_chemistry%mineral_zone%num_minerals_eq_var_act=&
                                                this%solid_chemistry%mineral_zone%num_minerals_eq_var_act+1
                                            this%solid_chemistry%mineral_zone%num_minerals_eq=&
                                                this%solid_chemistry%mineral_zone%num_minerals_eq+1
                                            this%solid_chemistry%mineral_zone%num_minerals=&
                                                this%solid_chemistry%mineral_zone%num_minerals+1
                                                
                                        else if (ind_cstr<=this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq) then
                                            !> Mineral with constant activity (pure phase at saturation)
                                            ind_mins(j)=ind_cstr !> store mineral index - chapuza repetida ademas
                                            j=j+1 !> increment mineral counter
                                            !> Update constant activity mineral counters in reactive zone
                                            this%solid_chemistry%reactive_zone%num_minerals_cst_act=&
                                                this%solid_chemistry%reactive_zone%num_minerals_cst_act+1
                                            this%solid_chemistry%reactive_zone%num_minerals=&
                                                this%solid_chemistry%reactive_zone%num_minerals+1
                                            !> Update constant activity mineral counters in mineral zone
                                            this%solid_chemistry%mineral_zone%num_minerals_eq_cst_act=&
                                                this%solid_chemistry%mineral_zone%num_minerals_eq_cst_act+1
                                            this%solid_chemistry%mineral_zone%num_minerals_eq=&
                                                this%solid_chemistry%mineral_zone%num_minerals_eq+1
                                            this%solid_chemistry%mineral_zone%num_minerals=&
                                                this%solid_chemistry%mineral_zone%num_minerals+1
                                                
                                        else if (ind_cstr<=this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq+&
                                            this%solid_chemistry%reactive_zone%chem_syst%cat_exch_zone%num_exch_cats) then
                                            !> Cation exchange site (surface chemistry)
                                            this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats=&
                                                this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats+1
                                                
                                        else if (ind_cstr<=this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq+&
                                            this%solid_chemistry%reactive_zone%chem_syst%cat_exch_zone%num_exch_cats+&
                                            this%solid_chemistry%reactive_zone%chem_syst%gas_phase%num_gases_eq_var_act) then
                                            !> Gas with variable activity (non-ideal gas phase)
                                            this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act=&
                                                this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act+1
                                            this%solid_chemistry%reactive_zone%gas_phase%num_var_act_species=&
                                                this%solid_chemistry%reactive_zone%gas_phase%num_var_act_species+1
                                            this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq=&
                                                this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq+1
                                            gas_flag=1 !> set flag indicating gas phase is present
                                            
                                        else
                                            !> Gas with constant activity (ideal gas phase)
                                            this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act=&
                                                this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1
                                            this%solid_chemistry%reactive_zone%gas_phase%num_cst_act_species=&
                                                this%solid_chemistry%reactive_zone%gas_phase%num_cst_act_species+1
                                            this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq=&
                                                this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq+1
                                            gas_flag=1 !> set flag indicating gas phase is present
                                        end if
                                    end if !> end heterogeneous check
                                else !> constraint species not found in chemical system
                                    error stop "Constrain not found in chemical system"
                                end if
                                        
                                !call constrain%is_gas(flag_gas)
                                !if (flag_gas==.true.) then
                                !    call gas%set_name(constrain%name)
                                !    call this%gas_chemistry%reactive_zone%gas_phase%is_gas_in_gas_phase(gas,flag,gas_ind)
                                !    this%gas_chemistry%activities(gas_ind)=ctot
                                !else if (ind_cstr<=this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq) then
                                !    call THIS%solid_chemistry%reactive_zone%ind_non_flow_species(ind_cstr)%copy_species(this%solid_chemistry%reactive_zone%chem_syst%minerals(this%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin+ind_cstr)%mineral)
                                !else
                                !    call THIS%solid_chemistry%reactive_zone%ind_non_flow_species(ind_cstr-this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq-this%solid_chemistry%reactive_zone%chem_syst%gas_phase%num_gases_eq_cst_act)%copy_species(this%solid_chemistry%reactive_zone%chem_syst%cat_exch%surf_compl(this%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin+ind_cstr)%mineral)
                                !end if
                                !if (ind_cstr>this%solid_chemistry%reactive_zone%chem_syst%num_minerals_cst_act .AND. ind_cstr<this%solid_chemistry%reactive_zone%chem_syst%num_redox_eq_reacts) then
                                !else if (ind_cstr>this%solid_chemistry%reactive_zone%chem_syst%speciation_alg%num_eq_reactions .AND. ind_cstr<this%solid_chemistry%reactive_zone%chem_syst%speciation_alg%num_eq_reactions) then
                                !end if
                                !call this%solid_chemistry%reactive_zone%chem_syst%gas_phase%is_gas_in_gas_phase(gas,flag,gas_ind)
                            else !> constraint not found in chemical system
                                error stop "Constrain not found in chemical system"
                            end if
                            
                        !> (COMMENTED) Icon 5 - not yet implemented
                        !else if (icon==5) then
                        !    n_icons(6)=n_icons(6)+1
                        
                        else !> invalid icon type
                            error stop "icon option not implemented yet"
                        end if
                        
                        !> ================================================================
                        !> Store concentration and icon data for current species
                        !> ================================================================
                        !this%concentrations(this%indices_aq_species(k))=guess !> (COMMENTED) alternative indexing
                        icons(k)=icon !> store icon type for this primary species
                        !ctots(this%indices_aq_species(k))=ctot !> (COMMENTED) alternative indexing
                        
                        !> Set concentration based on species type
                        if (icon==3 .and. aq_species%name=='h2o(p)') then !> if icon=3 and species is water
                            call this%compute_conc_ideal_water() !> compute ideal water concentration (55.5 M)
                        else !> for all other species
                            this%concentrations(k)=guess !> store initial guess as concentration
                        end if
                        !icons(k)=icon !> (COMMENTED) duplicate icon storage
                        ctots(k)=ctot !> store total concentration
                        k=k+1 !> increment primary species counter !> aqui hay que verificar dimension
                    !else !> (COMMENTED) species not in aqueous phase
                    !    error stop !> stop with error
                    !end if
                end if !> end of species processing
            end do !> end of icon reading loop
            
            !> ================================================================
            !> SECTION 4: Update chemistry counters and setup zones
            !> ================================================================
            call this%solid_chemistry%reactive_zone%cat_exch_zone%compute_num_surf_compl() !> compute number of surface complexes
            call this%solid_chemistry%reactive_zone%gas_phase%compute_num_species_phase() !> compute number of gas species
            !call this%solid_chemistry%reactive_zone%rearrange_ind_non_flow_species() !> rearrange non-flowing species array
            call this%solid_chemistry%reactive_zone%set_num_solids() !> set total number of solid species
            call this%solid_chemistry%reactive_zone%allocate_ind_non_flow_species(ll) !> allocate non-flowing species indices
            call this%solid_chemistry%reactive_zone%set_ind_non_flow_species(ind_non_flow_species(1:ll)) !> set non-flowing species indices
            call this%solid_chemistry%reactive_zone%allocate_ind_mins() !> allocate mineral indices array
            !> Set mineral indices in chemical system
            call this%solid_chemistry%reactive_zone%set_ind_mins_chem_syst(&
                ind_mins(1:this%solid_chemistry%reactive_zone%num_minerals)) !> pass only valid mineral indices
            !call this%solid_chemistry%reactive_zone%set_ind_eq_reacts() !> set indices of heterogeneous equilibrium reactions
            call this%solid_chemistry%compute_num_solids_solid_chem() !> compute total number of solids in solid chemistry
            
            !> Setup mineral zone if associated
            if (associated(this%solid_chemistry%mineral_zone)) then
                call this%solid_chemistry%mineral_zone%set_num_mins_kin_min_zone() !> set number of kinetic minerals
            end if
            
            !> ================================================================
            !> Handle gas chemistry based on gas_flag and num_gas_zones
            !> ================================================================
            if (gas_flag==1 .and. num_gas_zones==1) then !> if gas phase found and exactly 1 gas zone
                call this%set_gas_chemistry(gas_chem) !> link gas chemistry to aqueous chemistry - chapuza
            else if (gas_flag==1 .and. num_gas_zones==0) then !> if gas found but no zones allocated
                error stop "Gas zones not allocated but gas phase found"
            else if (gas_flag==1 .and. num_gas_zones>1) then !> if gas found and multiple zones
                error stop "Multiple gas zones not implemented yet"
            end if
        end if !> end of 'guess' label check
        
    !> ================================================================
    !> SECTION 5: Compute activity coefficient constants
    !> ================================================================
    do i=1,this%aq_phase%num_species !> loop over all aqueous species
        !> Compute activity coefficient constants for each species
        call this%aq_phase%aq_species(i)%params_act_coeff%compute_csts(this%aq_phase%aq_species(i)%valence,this%params_aq_sol,model)
    end do
    
    !> ================================================================
    !> SECTION 6: Set speciation algebra for reactive zone
    !> ================================================================
!> We set speciation algebra attribute in reactive zone
    allocate(swap(2)) !> allocate swap array for potential species swapping (dimension 2)
    flag_comp=.false. !> initialize to FALSE - we use the component matrix with constant activity species
    call this%set_spec_alg_aq_chem(flag_comp,flag_surf,flag_Se,swap,ind_eq_reacts(1:n)) !> setup speciation algebra, returns flag_Se and swap if needed
    
    !> (COMMENTED) Alternative speciation algebra setup approach
    ! call this%solid_chemistry%reactive_zone%speciation_alg%set_flag_comp(.false.) !> set component flag
    ! if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then !> if surface complexes exist
    !     flag_surf=.true. !> enable surface chemistry flag
    ! else !> no surface chemistry
    !     flag_surf=.false. !> disable flag
    !     call this%solid_chemistry%allocate_conc_solids() !> allocate solid concentrations
    !     call this%solid_chemistry%allocate_activities() !> allocate activities
    !     call this%solid_chemistry%allocate_log_act_coeffs_solid_chem() !> allocate log activity coefficients
    ! end if
    ! call this%solid_chemistry%reactive_zone%speciation_alg%set_flag_cat_exch(flag_surf) !> set exchange flag
    ! call this%solid_chemistry%reactive_zone%set_speciation_alg_dimensions() !> set algebra dimensions
    ! call this%set_ind_species() !> set species indices
    ! call this%solid_chemistry%reactive_zone%set_ind_eq_reacts() !> set equilibrium reaction indices
    ! call this%solid_chemistry%reactive_zone%set_stoich_mat_react_zone() !> set stoichiometric matrix
    ! call this%solid_chemistry%reactive_zone%set_ind_mins_stoich_mat() !> set mineral indices in stoich matrix
    ! call this%solid_chemistry%reactive_zone%set_ind_gases_stoich_mat() !> set gas indices in stoich matrix
    ! allocate(swap(2)) !> allocate swap array
    ! call this%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(flag_Se,swap) !> compute algebra arrays
    
    !> ================================================================
    !> Handle species swapping if needed (flag_Se = TRUE)
    !> ================================================================
    if (flag_Se .eqv. .true.) then !> if species swap is required for numerical stability
        !> We swap indices constrains (chapuza)
        allocate(ind_swap(2)) !> allocate swap indices array (dimension 2)
        aux_ind_cstr=indices_constrains !> backup original constraint indices
        !indices_constrains(:,2)=indices_constrains(:,2)-this%solid_chemistry%reactive_zone%chem_syst%num_redox_eq_reacts !> (COMMENTED) adjust indices
        
        !> Find indices in constraint array corresponding to swap species
        do i=1,num_cstr !> loop over all constraints
            if (indices_constrains(i,2)==swap(1)) then !> if constraint matches first swap species
                ind_swap(1)=i !> store constraint index
            else if (indices_constrains(i,2)==swap(2)) then !> if constraint matches second swap species
                ind_swap(2)=i !> store constraint index
            else
                continue !> skip to next constraint
            end if
        end do
        
        !> Perform the swap: exchange constraint data at the two indices
        indices_constrains(ind_swap(1),:)=aux_ind_cstr(ind_swap(2),:) !> swap constraint 1 with backup of constraint 2
        indices_constrains(ind_swap(2),:)=aux_ind_cstr(ind_swap(1),:) !> swap constraint 2 with backup of constraint 1
    end if
    
!> ================================================================
!> SECTION 7: Compute initial concentrations (in molarities)
!> ================================================================
!> Determine appropriate initialization method based on icon types, model, and Jacobian option
!> We compute initial concentrations (in molarities)
    !> aqui podrias usar polimorfismo (could use polymorphism here)
    !print *, n_icons([2,4]) !> (COMMENTED DEBUG) print counts of icon types 2 and 4
    
    !> ----------------------------------------------------------------
    !> Branch 1: If component constraints (icon=2) or equilibrium constraints (icon=4) exist
    !> ----------------------------------------------------------------
    if (sum(n_icons([2,4]))>0) then !> if there are icon 2 or icon 4 constraints
        
        if (flag_surf .eqv. .true.) then !> if surface chemistry (cation exchange) is present
            !> Surface chemistry branch
            if (Jac_opt==1) then !> if analytical Jacobian
                !> Initialize concentrations with analytical Jacobian and cation exchange
                call this%initialise_conc_anal_exch(icons,n_icons,indices_constrains,ctots,niter,CV_flag)
            else !> incremental Jacobian with surface chemistry not implemented
                error stop "Initialisation subroutine not implemented yet"
            end if
            
        else if (model==0) then !> if ideal activity coefficients model (no surface chemistry)
            !> Ideal model with analytical initialization
            call this%initialise_conc_anal_ideal(icons,n_icons,indices_constrains,ctots,niter,CV_flag)
            
        else !> non-ideal model (model > 0), no surface chemistry
            !> Non-ideal model branches
            if (Jac_opt==0) then !> if incremental coefficients
                !> Initialize with incremental Jacobian (numerical derivatives)
                call this%initialise_conc_incr_coeff(icons,n_icons,indices_constrains,ctots,niter,CV_flag)
            else if (Jac_opt==1) then !> if analytical Jacobian
                !> Initialize with analytical Jacobian
                call this%initialise_conc_anal(icons,n_icons,indices_constrains,ctots,niter,CV_flag)
            else !> invalid Jac_opt value
                error stop "Jac_opt not valid"
            end if
        end if
        
    !> ----------------------------------------------------------------
    !> Branch 2: If ideal model and only direct constraints (icon=1,3)
    !> ----------------------------------------------------------------
    else if (model==0) then !> if ideal model and no component/equilibrium constraints
        
        if (flag_surf .eqv. .true.) then !> if surface chemistry (cation exchange) is present
            !> Ideal model with surface chemistry (cation exchange)
            !> We compute activities of aqueous species
            call this%compute_activities_aq() !> compute activities from concentrations (ideal: activity = concentration)
            !> We get concentrations of primary species
            c1=this%get_c1() !> get array of primary species concentrations
            !> We set initial guess for concentrations of surface complexes
            allocate(act_ads_cats_ig(this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats))
                !> allocate array for adsorbed cation activities initial guess
            act_ads_cats_ig(1)=1d-6 !> first cation: very small activity - chapuza (temporary value)
            act_ads_cats_ig(2)=1d0-act_ads_cats_ig(1) !> second cation: remainder to sum to 1 - chapuza
            !call this%set_conc_sec_species(c2_ig) !> (COMMENTED) set secondary species concentrations
            
            !> Solve for surface species concentrations
            if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats==2) then !> if exactly 2 exchange cations
                !> Analytical solution for 2-cation exchange system
                call this%solid_chemistry%compute_conc_surf_ideal_bin(c1(&
                this%solid_chemistry%reactive_zone%cat_exch_zone%exch_cat_indices)) !> solve analytically for two exchange cations
            else !> if more than 2 exchange cations
                !> Numerical solution using Newton-Raphson method
                call this%solid_chemistry%compute_conc_surf_ideal_Newton(c1(&
                this%solid_chemistry%reactive_zone%cat_exch_zone%exch_cat_indices),act_ads_cats_ig,niter,CV_flag) !> iterate using Newton-Raphson
                !call this%solid_chemistry%compute_conc_surf_ideal_Picard(c1_aq(this%solid_chemistry%reactive_zone%cat_exch_zone%exch_cat_indices),act_ads_cats_ig,niter,CV_flag) !> (COMMENTED) alternative: Picard iteration
            end if
            
        else if (this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions>0) then !> if equilibrium reactions exist
            !> Ideal model with equilibrium reactions (no surface chemistry)
            call this%set_act_aq_species() !> set aqueous species activities from concentrations
            call this%compute_log_act_coeff_wat() !> compute log activity coefficient for water
            c1=this%get_c1() !> get activities of primary species
            log_act_coeffs=this%get_log_gamma() !> get log activity coefficients of secondary species
            allocate(c2_init(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
                !> allocate array for secondary species concentrations
            !> Compute secondary species concentrations from primary activities
            call this%compute_c2_from_c1_ideal(c1,log_act_coeffs,c2_init) !> ideal model: direct computation
            call this%compute_activities_aq() !> compute activities after secondary species are set

        else !> no equilibrium reactions, no surface chemistry
            !> Simple ideal model - only direct constraints
            call this%set_act_aq_species() !> set aqueous species activities
            call this%compute_log_act_coeff_wat() !> compute water log activity coefficient
            !call this%compute_molarities() !> (COMMENTED) compute molarities
        end if
        
    !> ----------------------------------------------------------------
    !> Branch 3: Non-ideal model with only direct constraints (icon=1,3)
    !> ----------------------------------------------------------------
    else !> non-ideal model (model > 0) with no component/equilibrium constraints
        !> Non-ideal model with direct constraints - use Picard iteration
        c1=this%get_c1() !> get concentrations of primary species
        allocate(c2_ig(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
            !> allocate array for secondary species initial guess
        c2_ig=1d-16 !> set very small initial guess - chapuza (temporary value)
        call this%set_conc_sec_species(c2_ig) !> set secondary species concentrations to initial guess
        !> Solve for secondary species using Picard iteration (activity coefficients updated each iteration)
        call this%compute_c2_from_c1_Picard(c1,c2_ig,c2_init,niter,CV_flag) !> iterate until convergence
    end if
    call this%compute_pH()
    call this%compute_salinity()

    !> Debug: print computed concentrations
    !print *, '[DEBUG read_wat_type_CHEPROO] Aqueous species concentrations:'
    do i = 1, this%aq_phase%num_species
        !print *, this%concentrations(i)
    end do
    !print *, '[DEBUG read_wat_type_CHEPROO] Aqueous species activities:'
    do i = 1, this%aq_phase%num_species
        !print *, this%activities(i)
    end do
    ! print *, '[DEBUG read_wat_type_CHEPROO] Secondary species concentrations:'
    ! do i = 1, this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
    !     print '(A,I4,A,I4,A,ES14.6)', '  sec  ', i, ' -> ind=', this%ind_sec_species(i), &
    !         '  c = ', this%concentrations(this%ind_sec_species(i))
    ! end do
    
end subroutine !> end of read_wat_type_CHEPROO