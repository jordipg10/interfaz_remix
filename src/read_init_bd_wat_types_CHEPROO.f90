!> @file read_init_bd_wat_types_CHEPROO.f90
!> @brief Reads initial and boundary water types from CHEPROO format input file
!> 
!> @details
!> This subroutine parses initial and boundary water type definitions from a CHEPROO-formatted
!> geochemical input file. It performs a two-pass reading strategy to handle complex chemical
!> compositions with multiple constraints and species.
!>
!> **Algorithm Overview:**
!>
!> The subroutine operates in several phases:
!>
!> **Phase 1: Pre-processing (First Pass)**
!> - Read activity coefficients model
!> - Read number of water types
!> - Allocate arrays for water types, solid types, reactive zones, gas zones
!> - For each water type:
!>   * Read index and temperature
!>   * Initialize aqueous phase and species indices
!>   * Read water type name
!>   * Read species with icon-based constraints:
!>     - icon=1: Fixed concentration
!>     - icon=2: Fixed total concentration
!>     - icon=3: Fixed activity (e.g., pH)
!>     - icon=4: Constrained by equilibrium reaction
!>   * Count primary aqueous species and constraints
!>
!> **Phase 2: Main Processing (Second Pass)**
!> - Rewind file and search for 'INITIAL AND BOUNDARY WATER TYPES' section
!> - For each water type:
!>   * Set up reactive zone with chemical system and CV parameters
!>   * Set up mineral zone
!>   * Handle cation exchange zones if present
!>   * Read detailed water type composition (delegates to read_wat_type_CHEPROO)
!>   * Set up solid chemistry associations
!>
!> **Phase 3: Post-processing**
!> - Eliminate constant activity species from chemical system component matrix
!> - Eliminate constant activity species from water type component matrices
!> - Set aqueous phase indices for kinetic mineral reactions
!> - Set aqueous phase indices for linear kinetic reactions
!>
!> @param[in,out] this Chemistry object containing all chemical system data
!> @param[in] unit File unit number for reading CHEPROO input file
!>
!> @note The term "chapuza" (Spanish for "workaround/hack") appears frequently in comments,
!>       indicating areas where the implementation uses temporary or non-ideal solutions
!>
!> @warning This subroutine modifies the chemical system's speciation algebra by eliminating
!>          constant activity species, which affects all subsequent calculations
!>
!> @warning The icon constraint system must use values 1-4 only; other values will cause
!>          program termination
!>
!> @see read_wat_type_CHEPROO For detailed water type composition reading
!> @see chemistry_c For the main chemistry object structure
!> @see reactive_zone_c For reactive zone setup
!>
!> @author jordi
!> @date November 2025
!>
subroutine read_init_bd_wat_types_CHEPROO(this,unit)
!> ================================================================
!> MODULE IMPORTS
!> ================================================================
    use chemistry_m, only: chemistry_c !> main chemistry class
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use reactive_zone_m, only: reactive_zone_c, compare_react_zones
    use mineral_zone_m, only: mineral_zone_c
    use solid_chemistry_m, only: solid_chemistry_c
    use gas_chemistry_m, only: gas_chemistry_c
    use aq_species_m, only: aq_species_c
    use gas_species_m, only: gas_species_c !> gas species class (constraint membership test)
    use mineral_m, only: mineral_c !> mineral class (constraint membership test)
    
    implicit none !> no implicit variable typing
    
!> ================================================================
!> SUBROUTINE PARAMETERS
!> ================================================================
    class(chemistry_c) :: this !> chemistry object (contains all chemical system data)
    integer(kind=4), intent(in) :: unit !> file unit number for reading input file
    !type(solid_chemistry_c), intent(inout), allocatable :: init_cat_exch_zones(:) !> (COMMENTED) initial cation exchange zones
    
!> ================================================================
!> LOCAL VARIABLES - Counters and Indices
!> ================================================================
    integer(kind=4) :: i !> loop counter for water types
    integer(kind=4) :: j !> water type index (may differ from loop counter i)
    integer(kind=4) :: k !> counter for primary aqueous species
    integer(kind=4) :: l !> counter for dissolved solids
    integer(kind=4) :: nwtype !> total number of water types
    integer(kind=4) :: icon !> constraint icon (1=activity, 2=concentration, 3=pH, 4=equilibrium)
    integer(kind=4) :: n_p_aq !> number of primary aqueous species (unused)
    integer(kind=4) :: gas_ind !> gas species index (unused)
    integer(kind=4) :: min_ind !> mineral index (unused)
    integer(kind=4) :: model !> activity coefficients model identifier
    integer(kind=4) :: niter !> number of iterations (output from read_wat_type_CHEPROO)
    integer(kind=4) :: sp_ind !> species index in aqueous phase
    integer(kind=4) :: num_gas_zones_loc !> local number of gas zones passed to read_wat_type_CHEPROO (1 when a dummy gas zone is provided)
    integer(kind=4) :: n_geq !> number of equilibrium gases in the chemical system (dummy gas zone size)
    integer(kind=4) :: n_geq_cst !> number of constant-activity equilibrium gases in the chemical system
    integer(kind=4) :: n_geq_var !> number of variable-activity equilibrium gases in the chemical system
    integer(kind=4) :: jg !> loop counter over equilibrium gases (dummy gas zone build)
    integer(kind=4) :: kk !> counter for constant-activity equilibrium gas indices (dummy gas zone build)
    logical :: build_dummy_gas !> TRUE when a dummy gas zone must be provided (gas used only as a speciation constraint)
    
!> ================================================================
!> LOCAL VARIABLES - Dynamic Arrays
!> ================================================================
    integer(kind=4), allocatable :: cols(:) !> column indices (size 2, purpose unclear)
    integer(kind=4), allocatable :: num_aq_prim_array(:) !> number of primary aqueous species per water type
    integer(kind=4), allocatable :: num_cstr_array(:) !> number of constraints (icon=4) per water type
    logical, allocatable :: wat_has_gas(:) !> TRUE if the water type has a gas as an icon=4 constraint (needs a gas zone)
    logical, allocatable :: wat_has_min(:) !> TRUE if the water type has a mineral as an icon=4 constraint (needs a mineral zone)
    real(kind=8), allocatable :: gas_pp(:,:) !> per-water-type partial pressures of the equilibrium gases (ctot of each icon=4 gas constraint)
    real(kind=8), allocatable :: min_act(:,:) !> per-water-type activities of the equilibrium minerals (ctot of each icon=4 mineral constraint; 1 for a pure phase)
    
!> ================================================================
!> LOCAL VARIABLES - Strings and Names
!> ================================================================
    character(len=256) :: prim_sp_name !> primary species name (unused)
    character(len=256) :: constrain !> constraint name (unused)
    character(len=256) :: label !> section label or keyword read from file
    character(len=256) :: name !> water type name
    
!> ================================================================
!> LOCAL VARIABLES - Real Numbers
!> ================================================================
    real(kind=8) :: guess !> initial guess for concentration (unused)
    real(kind=8) :: c_tot !> total concentration (unused)
    real(kind=8) :: temp !> temperature in Celsius (converted to Kelvin)
    real(kind=8) :: conc !> concentration value (unused)
    real(kind=8), allocatable :: eq_csts(:) !> equilibrium constants array (unused)
    
!> ================================================================
!> LOCAL VARIABLES - Logical Flags
!> ================================================================
    logical :: CV_flag !> control volume flag (output from read_wat_type_CHEPROO)
    logical :: flag !> general purpose flag for search operations
    logical :: flag_gas !> TRUE if an icon=4 constraint species is a gas in the chemical system
    logical :: flag_min !> TRUE if an icon=4 constraint species is a mineral in the chemical system
    logical :: flag_surf !> surface flag (unused)
    logical :: flag_comp !> component flag (unused)
    logical :: flag_Se !> selenium flag (unused)
    
!> ================================================================
!> LOCAL VARIABLES - Chemistry Objects
!> ================================================================
    type(reactive_zone_c) :: react_zone !> default reactive zone object (unused)
    type(reactive_zone_c), allocatable :: react_zones(:) !> array of reactive zones (one per water type)
    type(reactive_zone_c), allocatable :: gas_rz_tmp(:) !> local temp reactive zones backing the per-water-type gas zones during reading (their gas phase is copied into react_zones_wat_types after dedup)
    type(solid_chemistry_c), allocatable :: sol_types_loc(:) !> per-water-type solid types (local temp; distinct set copied into sol_types after dedup)
    type(gas_chemistry_c), allocatable :: gas_zones_loc(:) !> per-water-type gas zones (local temp; distinct set copied into gas_zones_wat_types after dedup)
    type(mineral_zone_c), allocatable :: min_zones(:) !> array of mineral zones (one per water type)
    type(solid_chemistry_c) :: solid_chem !> default solid chemistry object (unused)
    !type(solid_chemistry_c), allocatable :: solid_chems(:) !> (COMMENTED) array of solid chemistry objects
    type(aq_species_c) :: aq_species !> temporary aqueous species object for searching
    type(gas_species_c) :: gas_constr !> temporary gas species object for constraint membership test
    type(mineral_c) :: min_constr !> temporary mineral object for constraint membership test
    !type(mineral_c) :: mineral !> temporary mineral object (unused)
    !type(aq_phase_c) :: old_aq_phase !> old aqueous phase object (unused, for rearrangement)
    

!> ================================================================
!> SECTION 1: Read activity coefficients model
!> ================================================================
    !> Read activity coefficients model identifier
    read(unit,*) this%act_coeffs_model
        !> model identifier (e.g., 1=ideal, 2=Debye-Huckel, 3=Davies, etc.)
    
!> ================================================================
!> SECTION 2: Read number of water types and allocate structures
!> ================================================================
    !> Read total number of water types
    read(unit,*) nwtype !> number of initial and boundary water type definitions
    
!> ================================================================
!> SECTION 3: Cation exchange (surface) zones are read independently
!> ----------------------------------------------------------------
!> Surface (cation exchange) zones are read AFTER the water types, by
!> read_init_cat_exch_zones_CHEPROO, so they are independent of the water
!> types. The water type <-> exchanger coupling is resolved later through the
!> deferred solid-chemistry linking path (init_cat_exch_zones -> materials ->
!> target_solids -> waters); nothing about surfaces is done here.
!> ================================================================
    
!> ================================================================
!> SECTION 4: Allocate chemistry arrays for water types
!> ================================================================
    !> Allocate water types array (1:1 with the water types read from the file)
    call this%allocate_wat_types(nwtype) !> allocate water types array (size = nwtype)
    !> NOTE: sol_types, react_zones_wat_types and gas_zones_wat_types are NOT allocated
    !> here. Each water type's solid type / reactive zone / gas zone is built into a
    !> LOCAL temporary during reading; the persistent chemistry arrays are allocated at
    !> the END (SECTION 8c) with the DISTINCT counts, because several water types may
    !> share one solid type, gas zone and/or reactive zone.
    
    !> ================================================================
    !> SECTION 4b: Build the per-water-type gas zones
    !> ----------------------------------------------------------------
    !> The gas zones built here are associated ONLY to the water types
    !> (gas_zones_wat_types); this section is independent of how many global gas
    !> zones (this%num_gas_zones) the chemistry object has. Whenever the chemical
    !> system defines equilibrium gases (e.g. co2(g)), each water type is given
    !> its own gas zone holding all those equilibrium gases, backed by a local
    !> temporary reactive zone (gas_rz_tmp) during reading; after deduplication
    !> (SECTION 8c) the gas zone shares the water type's own react_zones_wat_types.
    !> This lets a water type use a gas as an equilibrium constraint (icon=4, e.g.
    !> hco3- constrained by co2(g)), which sets gas_flag=1 in read_wat_type_CHEPROO
    !> and requires a gas_chemistry object. Partial pressures are nominal (1 atm):
    !> the actual PCO2 of each water is set by its constraint total (ctot) during
    !> speciation.
    !> ================================================================
    build_dummy_gas=(this%chem_syst%gas_phase%num_gases_eq>0)
    num_gas_zones_loc=0
    if (build_dummy_gas) then
        num_gas_zones_loc=1 !> the water type's own gas zone is provided to read_wat_type_CHEPROO
        n_geq=this%chem_syst%gas_phase%num_gases_eq
        n_geq_cst=this%chem_syst%gas_phase%num_gases_eq_cst_act
        n_geq_var=this%chem_syst%gas_phase%num_gases_eq_var_act
        allocate(gas_rz_tmp(nwtype),gas_zones_loc(nwtype)) !> local per-water-type gas zones and their backing reactive zones (distinct set copied into gas_zones_wat_types/react_zones_wat_types after dedup)
        allocate(gas_pp(n_geq,nwtype)) !> per-water-type gas partial pressures
        gas_pp=1d0 !> default nominal 1 atm; overwritten by each water's icon=4 gas-constraint ctot
        do j=1,nwtype
            !> Reactive zone of the dummy gas zone (gas phase only)
            call gas_rz_tmp(j)%set_chem_syst_react_zone(this%chem_syst)
            call gas_rz_tmp(j)%set_CV_params(this%CV_params)
            call gas_rz_tmp(j)%gas_phase%allocate_gases(n_geq)
            call gas_rz_tmp(j)%gas_phase%set_num_var_act_species_phase(n_geq_var)
            call gas_rz_tmp(j)%gas_phase%set_num_cst_act_species_phase(n_geq_cst)
            call gas_rz_tmp(j)%gas_phase%compute_num_species_phase()
            call gas_rz_tmp(j)%gas_phase%set_num_gases_eq(n_geq)
            call gas_rz_tmp(j)%gas_phase%set_num_gases_eq_cst_act(n_geq_cst)
            call gas_rz_tmp(j)%gas_phase%set_num_gases_eq_var_act(n_geq_var)
            call gas_rz_tmp(j)%gas_phase%set_num_gases_kin(0)
            call gas_rz_tmp(j)%set_num_solids()
            call gas_rz_tmp(j)%set_num_non_flow_species()
            call gas_rz_tmp(j)%set_speciation_alg_dimensions(.true.)
            !> Dummy gas zone (gas_chemistry) linked to that reactive zone
            call gas_zones_loc(j)%set_reactive_zone(gas_rz_tmp(j))
            call gas_zones_loc(j)%set_temp(298.15d0)
            call gas_zones_loc(j)%set_volume(1d0)
            call gas_zones_loc(j)%allocate_partial_pressures()
            call gas_zones_loc(j)%allocate_conc_gases()
            call gas_zones_loc(j)%allocate_log_act_coeffs_gases()
            call gas_zones_loc(j)%allocate_var_act_species_indices(n_geq_var)
            call gas_zones_loc(j)%allocate_cst_act_species_indices(n_geq_cst)
            call gas_zones_loc(j)%allocate_ind_gases_eq_cst_act()
            !> Copy the equilibrium gases from the chemical system, set nominal partial pressures
            kk=0
            do jg=1,n_geq
                call gas_rz_tmp(j)%gas_phase%gases(jg)%copy_species(&
                    this%chem_syst%gas_phase%gases(jg))
                gas_zones_loc(j)%activities(jg)=1d0 !> nominal partial pressure (atm)
                if (jg<=this%chem_syst%gas_phase%num_gases_eq_cst_act) then
                    kk=kk+1
                    gas_zones_loc(j)%ind_gases_eq_cst_act(kk)=jg
                end if
            end do
            call gas_zones_loc(j)%compute_conc_gases_ideal()
            call gas_zones_loc(j)%compute_pressure()
            call gas_zones_loc(j)%compute_log_act_coeffs_gases()
            call gas_zones_loc(j)%set_indices_gases()
            call gas_rz_tmp(j)%set_ind_non_flow_species()
            call gas_rz_tmp(j)%set_ind_eq_reacts()
            call gas_rz_tmp(j)%allocate_ind_var_act_species()
            call gas_rz_tmp(j)%set_stoich_mat_react_zone()
            call gas_rz_tmp(j)%set_ind_gases_stoich_mat()
        end do
    end if

    !> Allocate local temporary arrays
    allocate(react_zones(nwtype),min_zones(nwtype),sol_types_loc(nwtype)) !> per-water-type local reactive/mineral zones and solid types
    allocate(cols(2)) !> allocate column indices array (size 2, purpose unclear)
    allocate(num_aq_prim_array(nwtype),num_cstr_array(nwtype)) !> allocate counters for species and constraints
    allocate(wat_has_gas(nwtype),wat_has_min(nwtype)) !> per-water-type gas/mineral constraint flags
    allocate(min_act(this%chem_syst%num_minerals,nwtype)) !> per-water-type mineral activities
    min_act=1d0 !> default activity 1 (pure phase); overwritten by each water's icon=4 mineral-constraint ctot
    
!> ================================================================
!> SECTION 5: Initialize counter arrays
!> ================================================================
    !> Initialize species and constraint counters to zero
    num_aq_prim_array=0 !> zero all elements (no primary species counted yet)
    num_cstr_array=0 !> zero all elements (no constraints counted yet)
    wat_has_gas=.false. !> no gas constraints detected yet
    wat_has_min=.false. !> no mineral constraints detected yet
    
!> ================================================================
!> SECTION 6: First pass - read water type structure and count species
!> ================================================================
     do i=1,nwtype !> loop over all water types (first pass)
        !> ----------------------------------------------------------------
        !> Read water type index and temperature
        !> ----------------------------------------------------------------
        read(unit,*) j, temp !> read index (j) and temperature in Celsius
        
        !> Validate water type index
        if (j<1 .or. j>nwtype) error stop !> index must be in range [1, nwtype]
        
        !ind_wat_type(i)=j !> (COMMENTED) store index mapping
        
        !> ----------------------------------------------------------------
        !> Initialize aqueous phase for this water type
        !> ----------------------------------------------------------------
        !> Set aqueous phase from chemical system
        call this%wat_types(j)%set_aq_phase(this%chem_syst%aq_phase)
            !> copy aqueous phase definition from chemical system
        
        !> Set default species indices (identity mapping initially)
        call this%wat_types(j)%set_indices_aq_species_aq_chem() !> indices_aq_species(i) = i
        
        !call this%wat_types(j)%set_ind_diss_solids_aq_chem() !> (COMMENTED) set default dissolved solids indices
        
        !> ----------------------------------------------------------------
        !> Set temperature and density
        !> ----------------------------------------------------------------
        !> Convert temperature from Celsius to Kelvin
        call this%wat_types(j)%set_temp(temp+273.15d0) !> T(K) = T(C) + 273.15
        
        !> Calculate water density at this temperature (Kell 1975 correlation, from the temperature just set)
        call this%wat_types(j)%compute_density_water() !> density = f(temperature) [kg/L]
        
        !> Debye-Hückel parameters A and B of this water type at its temperature and density
        !> (PHREEQC convention); read before the species constraints so that the speciation below
        !> (read_wat_type_CHEPROO) evaluates the activity coefficients at the water's temperature
        call this%wat_types(j)%compute_params_aq_sol() !> A, B = f(temperature, density)
        
        !call this%wat_types(j)%set_solid_chemistry(solid_chem) !> (COMMENTED) set solid chemistry
        
        !> ----------------------------------------------------------------
        !> Allocate concentration and activity arrays
        !> ----------------------------------------------------------------
        !> Allocate concentration array for aqueous species
        call this%wat_types(j)%allocate_conc_aq_species() !> allocate conc array (size = num aqueous species)
        
        !> Allocate log activity coefficients array
        call this%wat_types(j)%allocate_log_act_coeffs_aq_chem() !> allocate log(gamma) array
        
        !> Allocate activities array
        call this%wat_types(j)%allocate_activities_aq_species() !> allocate activities array (a = gamma * c)
        
        !> ----------------------------------------------------------------
        !> Read water type name
        !> ----------------------------------------------------------------
        read(unit,*) name !> read water type name/label
        call this%wat_types(j)%set_name(trim(name)) !> store trimmed name
        
        !> ----------------------------------------------------------------
        !> Read and process species constraints
        !> ----------------------------------------------------------------
        !> Read label to identify constraint format
        read(unit,*) label !> expect label containing 'icon' keyword
        
        if (index(label,'icon')/=0) then !> if label contains 'icon' substring
            !> Icon-based constraint format: 'icon, guess, ctot, constrain'
            
            !> Initialize species counters for this water type
            k=0 !> counter for primary aqueous species
            l=0 !> counter for dissolved solids (unused)
            
            do !> infinite loop - exit when sentinel '*' found
                !> Read species name, icon and the full constraint line (guess, ctot, constrain)
                read(unit,*) aq_species%name, icon, guess, c_tot, constrain !> species, icon, guess, ctot, constrain
                
                if (aq_species%name=='*') then !> if sentinel character found
                    exit !> exit species reading loop
                else !> process this species
                    !> ------------------------------------------------------------
                    !> Search for species in aqueous phase
                    !> ------------------------------------------------------------
                    call this%wat_types(j)%aq_phase%is_species_in_aq_phase(aq_species,flag,sp_ind)
                        !> Returns: flag (TRUE if found), sp_ind (index in aqueous phase)
                    
                    if (flag .eqv. .true.) then !> if species found in aqueous phase
                        !> Increment primary species counter
                        k=k+1 !> increment counter
                        
                        !> --------------------------------------------------------
                        !> Handle non-default species ordering
                        !> --------------------------------------------------------
                        if (sp_ind/=k) then !> if species index doesn't match counter (non-default order)
                            !> Store bidirectional index mapping
                            this%wat_types(j)%indices_aq_species(sp_ind)=k !> map aqueous phase index to water type index
                            this%wat_types(j)%indices_aq_phase(k)=sp_ind !> map water type index to aqueous phase index
                        end if
                        
                        !> --------------------------------------------------------
                        !> Handle special species (COMMENTED OUT)
                        !> --------------------------------------------------------
                        ! if (aq_species%name=='h2o') then !> if water species
                        !     !call this%wat_types(j)%set_ind_wat_aq_chem(k) !> set water index
                        ! else !> other species (dissolved solids)
                        !     l=l+1 !> increment dissolved solids counter
                        !     !this%wat_types(j)%ind_diss_solids(l)=k !> store dissolved solid index
                        !     ! if (aq_species%name=='h+') then !> if proton
                        !     !     call this%wat_types(j)%set_ind_prot_aq_chem(k) !> set proton index
                        !     ! !> aqui faltan especies (other special species missing)
                        !     ! end if
                        ! end if
                        
                        !> --------------------------------------------------------
                        !> Count primary aqueous species for this water type
                        !> --------------------------------------------------------
                        num_aq_prim_array(j)=num_aq_prim_array(j)+1
                            !> increment primary species counter for water type j
                        
                        !> --------------------------------------------------------
                        !> Process constraint icon and count constraints
                        !> --------------------------------------------------------
                        if (icon==4) then !> if icon=4 (equilibrium constraint)
                            !> Increment constraint counter
                            num_cstr_array(j)=num_cstr_array(j)+1
                                !> increment constraint counter for water type j
                            
                            !> Detect the TYPE of this heterogeneous constraint from the
                            !> constrain column: a gas means this water type needs a gas zone,
                            !> a mineral means it needs a mineral zone. (Surface complexes are
                            !> handled separately.)
                            gas_constr%name=trim(constrain) !> candidate gas name
                            call this%chem_syst%gas_phase%is_gas_in_gas_phase(gas_constr,flag_gas,gas_ind)
                            if (flag_gas) then
                                wat_has_gas(j)=.true. !> water type has a gas constraint
                                !> The gas's partial pressure is the ctot value on this
                                !> icon=4 line (whose constrain column names the gas).
                                if (build_dummy_gas .and. gas_ind<=n_geq) gas_pp(gas_ind,j)=c_tot
                            end if
                            min_constr%name=trim(constrain) !> candidate mineral name
                            call this%chem_syst%is_mineral_in_chem_syst(min_constr,flag_min,min_ind)
                            if (flag_min) then
                                wat_has_min(j)=.true. !> water type has a mineral constraint
                                !> The mineral's activity is the ctot value on this
                                !> icon=4 line (whose constrain column names the mineral).
                                min_act(min_ind,j)=c_tot
                            end if
                        
                        !> --------------------------------------------------------
                        !> Alternative icon handling (COMMENTED OUT)
                        !> --------------------------------------------------------
                        !if (icon==1) then !> icon=1: fixed activity
                        !    n_icons(1)=n_icons(1)+1 !> count icon=1 constraints
                        !else if (icon==2) then !> icon=2: fixed concentration
                        !    n_icons(2)=n_icons(2)+1 !> count icon=2 constraints
                        !    !n_aq_comp=n_aq_comp+1 !> increment aqueous components
                        !else if (icon==3) then !> icon=3: fixed pH (for H+)
                        !    n_icons(3)=n_icons(3)+1 !> count icon=3 constraints
                        !    if (aq_species%name=='h+') then !> if proton species
                        !        call this%aq_phase%set_ind_proton(sp_ind) !> set proton index
                        !        call this%set_pH(-log10(ctot)) !> set pH value
                        !    end if
                        !else if (icon==4) then !> icon=4: equilibrium constraint
                        !    n_icons(4)=n_icons(4)+1 !> count icon=4 constraints
                        !    !num_cstr=num_cstr+1 !> increment constraint counter
                        !    !call this%chem_syst%is_eq_reaction_in_chem_syst(constrain,flag,ind_cstr)
                        !    !if (flag==.true.) then !> if constraint reaction found
                        !    !    indices_constrains(l)=ind_cstr !> store constraint index
                        !    !    l=l+1 !> increment constraint counter
                        !    !else !> constraint reaction not found
                        !    !    error stop !> stop with error
                        !    !end if
                        !!else if (icon==5) then !> icon=5: (not implemented)
                        !!    n_icons(6)=n_icons(6)+1 !> count icon=5 constraints
                        else if (icon<1 .or. icon>4) then !> if icon out of valid range [1,4]
                            error stop "icon option not implemented yet" !> stop with error message
                        end if
                        
                        !> --------------------------------------------------------
                        !> Store constraint data (COMMENTED OUT)
                        !> --------------------------------------------------------
                        !this%concentrations(sp_ind)=guess !> store initial guess
                        !icons(sp_ind)=icon !> store icon type
                        !ctots(sp_ind)=ctot !> store total concentration
                        !constrains(sp_ind)=constrain !> store constraint name
                        !k=k+1 !> (COMMENTED) increment counter (redundant with k=k+1 above)
                        
                    else !> species not found in aqueous phase
                        !> --------------------------------------------------------
                        !> Error: Species not in aqueous phase
                        !> --------------------------------------------------------
                        error stop "Aqueous species not present in aqueous phase"
                            !> terminate program with error message
                    end if
                end if !> end if aq_species%name /= '*'
            end do !> end species reading loop
        else !> label does not contain 'icon'
            !> --------------------------------------------------------------------
            !> Error: Invalid water type format
            !> --------------------------------------------------------------------
            error stop "Error reading water type" !> terminate - expected 'icon' label
        end if
    end do !> end first pass water types loop
    
!> ================================================================
!> SECTION 6b: Apply the per-water-type gas partial pressures
!> ----------------------------------------------------------------
!> Each icon=4 gas constraint above recorded its gas partial pressure (the ctot
!> value on the line whose constrain column names the gas) into gas_pp. Set those
!> partial pressures on each water type's gas zone (replacing the nominal 1 atm
!> used while the zone structure was built) and recompute the ideal-gas
!> concentrations, total pressure and activity coefficients accordingly.
!> ================================================================
    if (build_dummy_gas) then
        do j=1,nwtype
            do jg=1,n_geq
                gas_zones_loc(j)%activities(jg)=gas_pp(jg,j)
            end do
            call gas_zones_loc(j)%compute_conc_gases_ideal()
            call gas_zones_loc(j)%compute_pressure()
            call gas_zones_loc(j)%compute_log_act_coeffs_gases()
        end do
    end if
    
!> ================================================================
!> SECTION 7: Second pass - search for detailed water type data
!> ================================================================
    !> Rewind file to beginning for second pass
    rewind(unit) !> reset file pointer to start
    
    do !> infinite loop - search for section header
        !> Read line looking for section header
        read(unit,*) label !> read label from current line
        
        if (label=='INITIAL AND BOUNDARY WATER TYPES') then !> if section header found
            !> ----------------------------------------------------------------
            !> Found section header - read detailed water type data
            !> ----------------------------------------------------------------
            
            !> Read activity coefficients model (again)
            read(unit,*) model !> model identifier (should match first pass value)
            
            !> Read number of water types (again)
            read(unit,*) nwtype !> number of water types (should match first pass value)
            
!> ================================================================
!> SECTION 8: Process each water type with full chemistry setup
!> ================================================================
            do i=1,nwtype !> loop over all water types (second pass)
                !> ------------------------------------------------------------
                !> Read water type index, temperature, and name
                !> ------------------------------------------------------------
                read(unit,*) j, temp !> read index and temperature in Celsius
                read(unit,*) name !> read water type name
                
                !> ------------------------------------------------------------
                !> Set up reactive zone for this water type
                !> ------------------------------------------------------------
                !> Set control volume parameters in reactive zone
                call react_zones(j)%set_CV_params(this%CV_params)
                    !> copy CV parameters to reactive zone
                
                !> Set chemical system in reactive zone
                call react_zones(j)%set_chem_syst_react_zone(this%chem_syst)
                    !> copy chemical system definition to reactive zone
                !> Attach the chemical system's cation-exchange zone so the water type
                !> equilibrates with the exchanger: its exchange reactions become part of the
                !> reactive zone's equilibrium reactions (num_eq_reactions), required by the
                !> speciation algebra to build the component matrix and by the exchange solve.
                if (this%chem_syst%cat_exch_zone%num_surf_compl>0) then
                    call react_zones(j)%set_cat_exch_zone()
                end if
                
                !> ------------------------------------------------------------
                !> Set up mineral zone for this water type
                !> ------------------------------------------------------------
                !> Set chemical system in mineral zone
                call min_zones(j)%set_chem_syst_min_zone(this%chem_syst)
                    !> copy chemical system definition to mineral zone
                
                !> ------------------------------------------------------------
                !> Set up solid chemistry associations (into LOCAL temporaries)
                !> ------------------------------------------------------------
                !> Link reactive zone to solid type
                call sol_types_loc(j)%set_reactive_zone(react_zones(j))
                    !> set default reactive zone for this solid type
                
                !> Link solid chemistry to water type
                call this%wat_types(j)%set_solid_chemistry(sol_types_loc(j)) !> chapuza (set solid type)
                
                !> ------------------------------------------------------------
                !> Read water type composition (independent of surface zones)
                !> ------------------------------------------------------------
                !> Set the (local) mineral zone for this solid type. The actual
                !> mineral/surface content is linked later via the deferred
                !> solid-chemistry path, so no cation exchange zone is used here.
                call sol_types_loc(j)%set_mineral_zone(min_zones(j))
                    !> set mineral zone for this solid type
                
                !> Read water type composition. When the chemical system has gases
                !> used only as constraints, the water type's own gas zone
                !> (gas_zones_loc) is supplied to the reader.
                if (build_dummy_gas) then
                    call this%wat_types(j)%read_wat_type_CHEPROO(num_aq_prim_array(j),num_cstr_array(j),num_gas_zones_loc,&
                        this%act_coeffs_model,this%Jac_opt,unit,niter,CV_flag,gas_zones_loc(j))
                else
                    call this%wat_types(j)%read_wat_type_CHEPROO(num_aq_prim_array(j),num_cstr_array(j),num_gas_zones_loc,&
                        this%act_coeffs_model,this%Jac_opt,unit,niter,CV_flag)
                end if
                
                !> ------------------------------------------------------------
                !> Populate the mineral zone from the reactive zone's minerals
                !> ------------------------------------------------------------
                !> read_wat_type_CHEPROO fills only the reactive zone (minerals + eq
                !> gases + surface complexes). The icon=4 mineral constraints are all
                !> EQUILIBRIUM minerals (no kinetic ones), so the mineral zone's
                !> equilibrium counts equal the reactive zone's mineral counts and the
                !> chemical-system mineral indices are copied directly (same order:
                !> variable-activity minerals first, then constant-activity).
                call min_zones(j)%set_num_mins_eq_min_zone(react_zones(j)%num_minerals)
                call min_zones(j)%set_num_mins_eq_cst_act_min_zone(react_zones(j)%num_minerals_cst_act)
                call min_zones(j)%set_num_mins_eq_var_act_min_zone(react_zones(j)%num_minerals_var_act)
                call min_zones(j)%set_num_mins_kin_cst_act_min_zone(0)
                call min_zones(j)%set_num_mins_kin_var_act_min_zone(0)
                call min_zones(j)%allocate_ind_chem_syst_min_zone(react_zones(j)%num_minerals)
                if (react_zones(j)%num_minerals>0) then
                    min_zones(j)%ind_min_chem_syst=react_zones(j)%ind_mins_chem_syst
                end if
                call min_zones(j)%allocate_ind_min_Sk(0)
            end do !> end second pass water types loop
            
            !> Exit search loop after processing section
            exit !> exit file search loop
            
         else !> label is not 'INITIAL AND BOUNDARY WATER TYPES'
            !> Continue searching for section header
            continue !> read next line
         end if
    end do !> end file search loop
    
!> ================================================================
!> SECTION 8c: Deduplicate the water-type solid / reactive / gas zones
!> ----------------------------------------------------------------
!> The per-water-type objects above were built one slot per water type. Two
!> water types with the SAME icon=4 constraint set share ONE reactive zone
!> (compared via compare_react_zones), and hence one solid type and (when they
!> have gas constraints) one gas zone. We now collapse the per-water-type slots
!> into the DISTINCT set, reallocate the chemistry arrays to the distinct counts
!> (they are NOT 1:1 with the water types) and re-link every water type to its
!> shared objects. All these links are pointers, so re-linking is mandatory
!> after the arrays are reallocated.
!> ================================================================
    block
        integer(kind=4) :: a !> water-type loop index
        integer(kind=4) :: b !> representative water type of a distinct mineral zone
        integer(kind=4) :: d !> distinct-index loop index
        integer(kind=4) :: n_rz !> number of distinct reactive zones
        integer(kind=4) :: n_mz !> number of distinct mineral zones
        integer(kind=4) :: n_st !> number of distinct solid types
        integer(kind=4) :: rz_map(nwtype) !> water type -> distinct reactive-zone index
        integer(kind=4) :: rz_rep(nwtype) !> distinct reactive-zone index -> representative water type
        integer(kind=4) :: mz_map(nwtype) !> water type -> distinct mineral-zone index
        integer(kind=4) :: mz_rep(nwtype) !> distinct mineral-zone index -> representative water type
        integer(kind=4) :: st_map(nwtype) !> water type -> distinct solid-type index
        integer(kind=4) :: st_rep(nwtype) !> distinct solid-type index -> representative water type
        integer(kind=4) :: n_gas !> number of gas zones (one per gas-bearing water type)
        integer(kind=4) :: gz_map(nwtype) !> water type -> gas-zone index (only set when wat_has_gas)
        integer(kind=4) :: p !> mineral position within a solid type's solids array
        logical :: matched !> TRUE if a water type matches an existing distinct object
        !> ------------------------------------------------------------
        !> 1a) Distinct REACTIVE zones (finest partition), compared by their
        !>     non-flowing species set (minerals + surface complexes + eq gases).
        !> ------------------------------------------------------------
        n_rz=0
        do a=1,nwtype
            matched=.false.
            do d=1,n_rz
                call compare_react_zones(react_zones(a),react_zones(rz_rep(d)),matched)
                if (matched) then
                    rz_map(a)=d
                    exit
                end if
            end do
            if (.not. matched) then
                n_rz=n_rz+1
                rz_rep(n_rz)=a
                rz_map(a)=n_rz
            end if
        end do
        !> ------------------------------------------------------------
        !> 1b) Distinct MINERAL zones (coarser: minerals only), compared by their
        !>     chemical-system mineral index set. Water types with the same minerals
        !>     but different gas/surface constraints share ONE mineral zone.
        !> ------------------------------------------------------------
        n_mz=0
        do a=1,nwtype
            matched=.false.
            do d=1,n_mz
                b=mz_rep(d)
                if (min_zones(a)%num_minerals==min_zones(b)%num_minerals) then
                    if (min_zones(a)%num_minerals==0) then
                        matched=.true.
                    else if (all(min_zones(a)%ind_min_chem_syst==min_zones(b)%ind_min_chem_syst)) then
                        matched=.true.
                    end if
                end if
                if (matched) then
                    mz_map(a)=d
                    exit
                end if
            end do
            if (.not. matched) then
                n_mz=n_mz+1
                mz_rep(n_mz)=a
                mz_map(a)=n_mz
            end if
        end do
        !> ------------------------------------------------------------
        !> 1c) Distinct SOLID types: a solid type is identified by the (reactive
        !>     zone, mineral zone) pair it points to. Two water types share a solid
        !>     type iff they share BOTH zones.
        !> ------------------------------------------------------------
        n_st=0
        do a=1,nwtype
            matched=.false.
            do d=1,n_st
                if (rz_map(a)==rz_map(st_rep(d)) .and. mz_map(a)==mz_map(st_rep(d))) then
                    st_map(a)=d
                    matched=.true.
                    exit
                end if
            end do
            if (.not. matched) then
                n_st=n_st+1
                st_rep(n_st)=a
                st_map(a)=n_st
            end if
        end do
        !> ------------------------------------------------------------
        !> 1d) Gas zones are NOT deduplicated: each gas-bearing water type keeps its
        !>     OWN gas zone, because its partial pressures (e.g. PCO2) are set
        !>     per-water later during speciation and cannot be shared. Consequently
        !>     two gas zones of two water types that share a reactive zone will point
        !>     to the SAME reactive zone.
        !> ------------------------------------------------------------
        n_gas=0
        if (build_dummy_gas) then
            do a=1,nwtype
                if (wat_has_gas(a)) then
                    n_gas=n_gas+1
                    gz_map(a)=n_gas
                end if
            end do
        end if
        !> ------------------------------------------------------------
        !> 2) Allocate the persistent chemistry arrays with their (independent)
        !>    DISTINCT counts. These counts need NOT be equal: several solid types
        !>    may point to the same shared reactive zone or mineral zone, and several
        !>    gas zones may point to the same shared reactive zone.
        !> ------------------------------------------------------------
        allocate(this%react_zones_wat_types(n_rz))
        allocate(this%min_zones_wat_types(n_mz))
        allocate(this%sol_types(n_st))
        if (build_dummy_gas) allocate(this%gas_zones_wat_types(n_gas)) !> one gas zone per gas-bearing water type
        !> ------------------------------------------------------------
        !> 3a) Copy the distinct reactive zones. For gas-bearing reactive zones,
        !>     transfer the equilibrium gas phase (structure) onto the shared zone;
        !>     the per-water gas zones are built and linked separately (step 3d).
        !> ------------------------------------------------------------
        do d=1,n_rz
            call this%react_zones_wat_types(d)%copy_react_zone(react_zones(rz_rep(d)))
            if (build_dummy_gas .and. wat_has_gas(rz_rep(d))) then
                this%react_zones_wat_types(d)%gas_phase=gas_rz_tmp(rz_rep(d))%gas_phase
            end if
        end do
        !> ------------------------------------------------------------
        !> 3b) Copy the distinct mineral zones.
        !> ------------------------------------------------------------
        do d=1,n_mz
            call this%min_zones_wat_types(d)%copy_mineral_zone(min_zones(mz_rep(d)))
        end do
        !> ------------------------------------------------------------
        !> 3c) Copy the distinct solid types and link each to its SHARED reactive
        !>     zone and mineral zone (many solid types may point to the same zone).
        !> ------------------------------------------------------------
        do d=1,n_st
            call this%sol_types(d)%copy_solid_chemistry(sol_types_loc(st_rep(d)))
            call this%sol_types(d)%set_reactive_zone(this%react_zones_wat_types(rz_map(st_rep(d))))
            call this%sol_types(d)%set_mineral_zone(this%min_zones_wat_types(mz_map(st_rep(d))))
            !> Set the solid activities: each equilibrium mineral's activity is the ctot
            !> on its icon=4 constraint line (captured per water type in min_act; 1 for a
            !> pure phase). Minerals occupy the first num_minerals positions of the solids
            !> array (any surface complexes follow, left at the default 1).
            call this%sol_types(d)%allocate_activities()
            this%sol_types(d)%activities=1d0
            do p=1,this%sol_types(d)%reactive_zone%num_minerals
                this%sol_types(d)%activities(p)=&
                    min_act(this%sol_types(d)%reactive_zone%ind_mins_chem_syst(p),st_rep(d))
            end do
        end do
        !> ------------------------------------------------------------
        !> 3d) Build one gas zone per gas-bearing water type and link it to that
        !>     water's SHARED reactive zone (the SAME one its solid type points to).
        !>     Two gas zones may therefore point to the same reactive zone.
        !> ------------------------------------------------------------
        if (build_dummy_gas) then
            do a=1,nwtype
                if (wat_has_gas(a)) then
                    call this%gas_zones_wat_types(gz_map(a))%copy_gas_chemistry(gas_zones_loc(a))
                    call this%gas_zones_wat_types(gz_map(a))%set_reactive_zone(&
                        this%react_zones_wat_types(rz_map(a)))
                end if
            end do
        end if
        !> ------------------------------------------------------------
        !> 4) Re-link every water type to its shared (deduplicated) objects. A water
        !>    type reaches its reactive zone through its solid type; a gas zone is
        !>    linked only to water types that have a gas constraint.
        !> ------------------------------------------------------------
        do a=1,nwtype
            call this%wat_types(a)%set_solid_chemistry(this%sol_types(st_map(a)))
            if (build_dummy_gas .and. wat_has_gas(a)) then
                call this%wat_types(a)%set_gas_chemistry(this%gas_zones_wat_types(gz_map(a)))
            end if
        end do
    end block
    
!> ================================================================
!> SECTION 9: Post-processing - eliminate constant activity species
!> ================================================================
    !> Eliminate constant activity species from chemical system component matrix
    call this%chem_syst%change_spec_alg_chem_syst(this%CV_params%zero)
        !> remove species with constant activity (activity coefficients ~= 1) from speciation algebra
        !> modifies component matrix and speciation algebra
        !> zero threshold from CV_params (control volume parameters)

    !> Rearrange mineral zone indices after speciation algebra change
    do i=1,size(this%min_zones_wat_types) !> loop over the distinct (deduplicated) mineral zones
        call this%min_zones_wat_types(i)%set_ind_min_Sk()
    end do

    !> Eliminate constant activity species from water type component matrices
    do i=1,nwtype !> loop over all water types
        call this%wat_types(i)%change_spec_alg_aq_chem()
            !> remove constant activity species from this water type's speciation algebra
            !> updates component matrix for aqueous chemistry
    end do
    
!> ================================================================
!> SECTION 10: Rearrange indices (COMMENTED OUT)
!> ================================================================
!> Chapuza (workaround for index rearrangement - not currently used)
    ! do i=1,nwtype !> loop over water types
    !     !> rearrange indices after speciation algebra changes
    !     ! print *, this%wat_types(i)%indices_aq_species !> debug print
    !     ! print *, this%wat_types(i)%ind_prim_species !> debug print
    !     ! print *, this%wat_types(i)%ind_diss_solids !> debug print
    !     ! call this%wat_types(i)%set_ind_prim_sec_species() !> set primary/secondary species indices
    !     !call this%wat_types(i)%set_indices_aq_species_aq_chem() !> reset species indices
    !     !call this%wat_types(i)%set_ind_diss_solids_aq_chem() !> reset dissolved solids indices
    ! end do
    !do i=1,this%num_bd_this%wat_types !> loop over boundary water types
    !    call this%bd_this%wat_types(i)%rearrange_state_vars(old_aq_phase)
    !        !> rearrange state variables based on old aqueous phase definition
    !end do
    
!> ================================================================
!> SECTION 11: Set aqueous phase indices for kinetic mineral reactions
!> ================================================================
!> Chapuza (workaround for setting kinetic mineral reaction indices)
    do i=1,this%chem_syst%num_minerals_kin !> loop over kinetic mineral reactions
        !> Set aqueous phase indices for reactants in mineral kinetic reactions
        call this%chem_syst%min_kin_reacts(i)%set_indices_aq_phase_min(this%chem_syst%aq_phase,this%chem_syst%species)
            !> map mineral reaction species to aqueous phase indices
            !> needed for kinetic rate calculations
    end do
    
!> ================================================================
!> SECTION 12: Set aqueous phase indices for linear kinetic reactions
!> ================================================================
!> Chapuza (workaround for setting linear kinetic reaction indices)
    do i=1,this%chem_syst%num_lin_kin_reacts !> loop over linear kinetic reactions
        !> Set aqueous phase indices for reactants in linear kinetic reactions
        call this%chem_syst%lin_kin_reacts(i)%set_index_aq_phase_lin(this%chem_syst%aq_phase,this%chem_syst%species)

    end do
    
    !print *, this%sol_types(1)%mineral_zone%num_minerals_eq, this%sol_types(1)%mineral_zone%num_minerals_kin
        !> (COMMENTED) debug print: number of equilibrium and kinetic minerals in first solid type
        
end subroutine