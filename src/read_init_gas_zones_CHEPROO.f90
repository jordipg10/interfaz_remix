!> @file read_gas_zones_CHEPROO.f90
!> @brief Reads gas zones from CHEPROO format input file
!> 
!> @details
!> This subroutine parses initial gas zone definitions from a CHEPROO-formatted
!> geochemical input file. Gas zones represent gas phase compartments containing
!> gas species at specified partial pressures, temperatures, and volumes.
!>
!> **Algorithm Overview:**
!>
!> The subroutine operates in several phases:
!>
!> **Phase 1: Pre-processing and First Pass**
!> - Read number of gas zones
!> - Handle special case: if ngtype=0, return immediately (no gas zones)
!> - Allocate arrays for gas zones and indices
!> - For each gas zone (first pass):
!>   * Read zone index, temperature (Celsius), and volume
!>   * Read zone name
!>   * Read gas species and partial pressures
!>   * Count gases by type:
!>     - Equilibrium gases with constant activity
!>     - Equilibrium gases with variable activity
!>     - Kinetic gases with constant activity
!>     - Kinetic gases with variable activity
!>   * Set up reactive zone structure
!>   * Configure speciation algebra dimensions
!>
!> **Phase 2: Second Pass - Detailed Gas Zone Reading**
!> - Rewind file and search for 'INITIAL AND BOUNDARY GAS ZONES' section
!> - For each gas zone:
!>   * Re-read zone index, temperature, volume, name
!>   * Read gas species names and partial pressures
!>   * Assign gas species from chemical system
!>   * Set activities (activities = partial pressures for gases)
!>   * Compute gas concentrations (ideal gas law)
!>   * Compute total pressure
!>   * Compute log activity coefficients
!>   * Set gas indices
!>   * Check if gas phase equals previous zone (avoid duplication)
!>   * If unique, increment reactive zones counter
!>
!> **Phase 3: Integration with Reactive Zones**
!> - If surface/cation exchange reactive zones exist:
!>   * Create combined reactive zones (gas + surface + gas+surface)
!>   * Total zones = num_surf + num_gas + (num_gas × num_surf)
!> - If no surface zones:
!>   * Create simple gas reactive zones
!>
!> **Gas Classification:**
!> - **Equilibrium gases (cst act)**: Participate in equilibrium reactions, constant partial pressure
!> - **Equilibrium gases (var act)**: Participate in equilibrium reactions, variable partial pressure
!> - **Kinetic gases**: Produced/consumed by kinetic reactions
!>
!> @param[in,out] this Chemistry object containing all chemical system data
!> @param[in] unit File unit number for reading CHEPROO input file
!> @param[out] ngrz Number of gas reactive zones (output)
!>
!> @note The term "chapuza" (Spanish for "workaround/hack") appears in comments,
!>       indicating areas where the implementation uses temporary or non-ideal solutions
!>
!> @note Gas activities are equal to partial pressures (a_i = P_i)
!>
!> @note Ideal gas behavior is assumed for concentration calculations: n_i = P_i*V/(RT), where n_i is in mol
!>
!> @warning If ngtype=0, the subroutine returns immediately without processing
!>
!> @see chemistry_c For the main chemistry object structure
!> @see gas_chemistry_c For gas chemistry zone structure
!> @see gas_phase_c For gas phase structure
!>
!> @author jordi
!> @date November 2025
!>
subroutine read_gas_zones_CHEPROO(this,unit,ngrz)
!> ================================================================
!> MODULE IMPORTS
!> ================================================================
    use chemistry_m, only: chemistry_c
    use gas_chemistry_m, only: gas_chemistry_c
    use reactive_zone_m, only: reactive_zone_c
    use gas_species_m, only: gas_species_c
    use gas_phase_m, only: gas_phase_c, are_gas_phases_equal
        !> chemistry classes: main chemistry, gas chemistry, reactive zone, gas species, gas phase
    
    implicit none !> no implicit variable typing
    
!> ================================================================
!> SUBROUTINE PARAMETERS
!> ================================================================
    class(chemistry_c) :: this !> chemistry object (contains all chemical system data)
    integer(kind=4), intent(in) :: unit !> file unit number for reading input file
    !type(gas_chemistry_c), intent(out), allocatable :: this%gas_zones(:)
        !> (COMMENTED) gas zones array (now part of chemistry object)
    integer(kind=4), intent(out) :: ngrz !> number of gas reactive zones (output)
    !type(reactive_zone_c), intent(inout), allocatable, optional :: reactive_zones(:)
        !> (COMMENTED) optional reactive zones array
    
!> ================================================================
!> LOCAL VARIABLES - Counters and Indices
!> ================================================================
    integer(kind=4) :: num_rz !> total number of reactive zones after combining
    integer(kind=4) :: i !> loop counter (general purpose)
    integer(kind=4) :: j !> loop counter (general purpose, nested loops)
    integer(kind=4) :: k !> counter for equilibrium gases with constant activity
    integer(kind=4) :: ngtype !> total number of gas zone types
    integer(kind=4) :: igtype !> current gas zone type index
    integer(kind=4) :: nrwtype !> number of reactive water types (unused)
    integer(kind=4) :: gas_ind !> index of gas in gas phase of chemical system
    integer(kind=4) :: num_gases_loc !> number of gases in current zone
    integer(kind=4) :: num_gases_var !> number of gases with variable activity in current zone
    integer(kind=4) :: num_gases_cst !> number of gases with constant activity in current zone
    integer(kind=4) :: n_gas_eq !> number of equilibrium gases in current zone
    integer(kind=4) :: n_gas_kin !> number of kinetic gases in current zone
    integer(kind=4) :: num_surf_rz !> number of existing surface reactive zones
    integer(kind=4) :: n_gas_eq_cst_act !> number of equilibrium gases with constant activity
    integer(kind=4) :: n_gas_eq_var_act !> number of equilibrium gases with variable activity
    integer(kind=4) :: num_gas_zones !> counter for processed gas zones
    
!> ================================================================
!> LOCAL VARIABLES - Dynamic Arrays
!> ================================================================
    integer(kind=4), allocatable :: ind_gases(:) !> indices of gases in gas phase of chemical system
    integer(kind=4), allocatable :: ind_gas_zones(:) !> indices of processed gas zones
    
!> ================================================================
!> LOCAL VARIABLES - Strings
!> ================================================================
    character(len=256) :: str !> temporary string for reading labels/headings
    character(len=256) :: constrain !> constraint name (unused)
    character(len=256) :: label !> section label for file navigation
    character(len=256) :: name !> gas zone name
    
!> ================================================================
!> LOCAL VARIABLES - Real Numbers (Gas Properties)
!> ================================================================
    real(kind=8) :: guess !> initial guess value (unused)
    real(kind=8) :: conc !> concentration value (unused)
    real(kind=8) :: temp !> temperature in Celsius (converted to Kelvin)
    real(kind=8) :: part_press !> partial pressure of gas species
    real(kind=8) :: vol !> volume of gas zone
    
!> ================================================================
!> LOCAL VARIABLES - Chemistry Objects
!> ================================================================
    type(gas_species_c) :: gas !> temporary gas species object for searching
    type(gas_phase_c) :: gas_phase !> temporary gas phase object (unused)
    type(reactive_zone_c) :: react_zone !> temporary reactive zone object (unused)
    type(reactive_zone_c), allocatable :: aux_react_zones(:) !> auxiliary array for copying existing reactive zones
    type(reactive_zone_c), allocatable :: react_zones(:) !> array of reactive zones (one per gas zone)
    
!> ================================================================
!> LOCAL VARIABLES - Logical Flags
!> ================================================================
    logical :: flag_comp !> flag to compute speciation algebra dimensions
    logical :: flag !> general purpose flag for search/comparison operations

!> ================================================================
!> SECTION 1: Read number of gas zones and validate
!> ================================================================
    !> Read total number of gas zone types
    read(unit,*) ngtype !> number of gas zones to be defined
    
    !> Validate number of gas zones
    if (ngtype<0) then !> if negative number
        error stop "Number of gas zones cannot be negative" !> stop with error message
    else if (ngtype==0) then !> if zero gas zones
        !> No gas zones to process - return immediately
        ngrz=0 !> set output parameter to zero
        return !> exit subroutine early
    end if
    
!> ================================================================
!> SECTION 2: Allocate arrays for gas zones
!> ================================================================
    !allocate(this%gas_zones(ngtype)) !> (COMMENTED) allocate gas zones array
    
    !> Allocate gas zones array in chemistry object
    call this%allocate_gas_zones(ngtype) !> allocate array (size = ngtype)
    
    !> Allocate index tracking arrays
    allocate(ind_gas_zones(ngtype)) !> allocate array for storing processed zone indices
    allocate(ind_gases(this%chem_syst%gas_phase%num_species)) !> allocate array for gas indices
    allocate(react_zones(ngtype)) !> allocate reactive zones array (one per gas zone)

!> ================================================================
!> SECTION 3: Initialize counters and flags
!> ================================================================
    !> Initialize gas zones counter
    num_gas_zones=0 !> counter for processed gas zones
    
    !> Initialize gas reactive zones counter
    ngrz=1 !> counter for gas reactive zones (starts at 1, not 0)

    !> Initialize speciation algebra computation flag
    flag_comp=.true. !> flag to compute speciation algebra dimensions (by default: TRUE)
    
!> ================================================================
!> SECTION 4: First pass - read gas zone structure and count species
!> ================================================================
    do !> infinite loop - exit when all zones processed or sentinel found
        !> ----------------------------------------------------------------
        !> Read gas zone header information
        !> ----------------------------------------------------------------
        read(unit,*) igtype, temp, vol !> read zone index, temperature (Celsius), and volume
        
        !> Validate gas zone index
        if (igtype<1 .or. igtype>ngtype) error stop "Index of gas zone out of range"
            !> index must be in range [1, ngtype]
        
        !> Read gas zone name
        read(unit,*) name !> read zone name/label
        call this%gas_zones(igtype)%set_name(trim(name)) !> store trimmed name in gas zone
        
        !> Read headings line
        read(unit,*) str !> read headings or sentinel character
        
        !> ----------------------------------------------------------------
        !> Check for sentinel or gas data
        !> ----------------------------------------------------------------
        if (str=='*') then !> if sentinel character found
            !> No gas species in this zone - check if more zones remain
            if (num_gas_zones<ngtype) then !> if more zones expected
                continue !> continue to next zone
            else !> all zones processed
                exit !> exit first pass loop
            end if
        else if (index(str,'gas')/=0) then !> if heading contains 'gas' keyword
            !> Valid gas zone - proceed with reading gas species list
            num_gas_zones=num_gas_zones+1 !> increment gas zones counter
            ind_gas_zones(num_gas_zones)=igtype !> store gas zone index in array
            
            !> ----------------------------------------------------------------
            !> Initialize gas species counters for this zone
            !> ----------------------------------------------------------------
            num_gases_loc=0 !> counter: total gases in this zone
            num_gases_var=0 !> counter: variable activity gases
            num_gases_cst=0 !> counter: constant activity gases
            n_gas_eq=0 !> counter: equilibrium gases
            n_gas_eq_cst_act=0 !> counter: equilibrium gases with constant activity
            n_gas_eq_var_act=0 !> counter: equilibrium gases with variable activity
            n_gas_kin=0 !> counter: kinetic gases
            
            !> ----------------------------------------------------------------
            !> Read gas species list for this zone
            !> ----------------------------------------------------------------
            do !> infinite loop - exit when gas list ends
                read(unit,*) gas%name, part_press !> read gas name and partial pressure (atm)
                
                !> Check for end of gas list
                if (gas%name=='*') then !> if sentinel character found
                    !> ============================================================
                    !> End of gas list - configure reactive zone
                    !> ============================================================
                    
                    !> Set chemical system reference for this reactive zone
                    call react_zones(igtype)%set_chem_syst_react_zone(this%chem_syst)
                        !> links to parent chemistry system
                    
                    !> Set control volume parameters
                    call react_zones(igtype)%set_CV_params(this%CV_params)
                        !> sets CV for reactive transport calculations
                    
                    !> Allocate gas phase arrays
                    call react_zones(igtype)%gas_phase%allocate_gases(num_gases_loc)
                        !> allocate array for all gases in this zone
                    
                    !> Set gas phase activity categories
                    call react_zones(igtype)%gas_phase%set_num_var_act_species_phase(num_gases_var)
                        !> set number of variable activity gas species
                    call react_zones(igtype)%gas_phase%set_num_cst_act_species_phase(num_gases_cst)
                        !> set number of constant activity gas species
                    call react_zones(igtype)%gas_phase%compute_num_species_phase()
                        !> compute total species count (should equal num_gases_loc)
                    
                    !> Set gas reaction type categories
                    call react_zones(igtype)%gas_phase%set_num_gases_eq(n_gas_eq)
                        !> set number of equilibrium gases
                    call react_zones(igtype)%gas_phase%set_num_gases_eq_cst_act(n_gas_eq_cst_act)
                        !> set number of constant activity equilibrium gases
                    call react_zones(igtype)%gas_phase%set_num_gases_eq_var_act(n_gas_eq_var_act)
                        !> set number of variable activity equilibrium gases
                    call react_zones(igtype)%gas_phase%set_num_gases_kin(n_gas_kin)
                        !> set number of kinetic gases
                    
                    !> Set non-flowing species counts
                    call react_zones(igtype)%set_num_solids()
                        !> sets number of solid species in reactive zone
                    call react_zones(igtype)%set_num_non_flow_species()
                        !> sets total non-flowing species (solids + constant activity gases)
                    
                    !> Configure speciation algebra dimensions
                    call react_zones(igtype)%set_speciation_alg_dimensions(flag_comp)
                        !> sets dimensions for speciation calculations
                    
                    !> (Commented out code - duplicate gas phase checking)
                    !> Note: The following code would check if this gas phase matches
                    !> a previous one to avoid duplication, but is currently disabled
                    !if (num_gas_zones>1) then
                    !    call are_gas_phases_equal(react_zones(igtype)%gas_phase,react_zones(ind_this%gas_zones(num_gas_zones-1))%gas_phase,flag)
                    !    if (flag==.true.) then
                    !        call this%gas_zones(igtype)%set_reactive_zone(react_zones(ind_this%gas_zones(num_gas_zones-1)))
                    !    else
                    !        call this%gas_zones(igtype)%set_reactive_zone(react_zones(igtype))
                    !    end if
                    !end if
                    
                    !> Set reactive zone for this gas zone
                    call this%gas_zones(igtype)%set_reactive_zone(react_zones(igtype))
                        !> assigns configured reactive zone to gas zone (note: workaround implementation)
                    
                    !> Set thermodynamic properties
                    call this%gas_zones(igtype)%set_temp(temp+273.15)
                        !> convert temperature from Celsius to Kelvin
                    call this%gas_zones(igtype)%set_volume(vol)
                        !> set zone volume (m³)
                    
                    !> Allocate gas-specific arrays in gas zone
                    call this%gas_zones(igtype)%allocate_partial_pressures()
                        !> allocate array for partial pressures (atm)
                    call this%gas_zones(igtype)%allocate_conc_gases()
                        !> allocate array for gas concentrations (mol)
                    call this%gas_zones(igtype)%allocate_log_act_coeffs_gases()
                        !> allocate array for logarithms of activity coefficients
                    call this%gas_zones(igtype)%allocate_var_act_species_indices(num_gases_var)
                        !> allocate array for variable activity species indices
                    call this%gas_zones(igtype)%allocate_cst_act_species_indices(num_gases_cst)
                        !> allocate array for constant activity species indices
                    call this%gas_zones(igtype)%allocate_ind_gases_eq_cst_act()
                        !> allocate array for constant activity equilibrium gas indices
                    
                    exit !> exit gas species reading loop - zone complete
                    
                !> ------------------------------------------------------------
                !> Process individual gas species
                !> ------------------------------------------------------------
                else !> gas%name is not sentinel - process this gas
                    !> Check if gas exists in chemical system gas phase
                    call this%chem_syst%gas_phase%is_gas_in_gas_phase(gas,flag,gas_ind)
                        !> searches for gas by name, returns flag and index
                    
                    if (flag .eqv. .true.) then !> if gas found in chemical system
                        !> ----------------------------------------------------
                        !> Gas found - categorize and count
                        !> ----------------------------------------------------
                        num_gases_loc=num_gases_loc+1 !> increment total gas counter
                        ind_gases(num_gases_loc)=gas_ind
                            !> store index of gas in chemical system gas phase
                        
                        !> Check if constant activity equilibrium gas
                        if (this%chem_syst%gas_phase%gases(gas_ind)%cst_act_flag.eqv. .true. .and. gas_ind<=&
                        this%chem_syst%gas_phase%num_gases_eq) then
                            !> Gas has constant activity AND is equilibrium gas
                            num_gases_cst=num_gases_cst+1 !> increment constant activity counter
                            n_gas_eq=n_gas_eq+1 !> increment equilibrium gas counter
                            n_gas_eq_cst_act=n_gas_eq_cst_act+1 !> increment constant activity equilibrium gas counter
                        
                        !> (Commented out code - modify equilibrium constant)
                        !> Note: This would adjust equilibrium constant for constant activity gas
                        !> but is currently disabled
                            !this%chem_syst%eq_reacts(this%chem_syst%aq_phase%num_aq_complexes+&
                            !this%chem_syst%num_minerals_eq+gas_ind)%eq_cst=this%chem_syst%eq_reacts(&
                            !this%chem_syst%aq_phase%num_aq_complexes+this%chem_syst%num_minerals_eq+gas_ind)%eq_cst/part_press
                        
                        !> Check if variable activity equilibrium gas
                        else if (this%chem_syst%gas_phase%gases(gas_ind)%cst_act_flag.eqv. .false. .and. gas_ind<=&
                            this%chem_syst%gas_phase%num_gases_eq) then
                            !> Gas has variable activity AND is equilibrium gas
                            num_gases_var=num_gases_var+1 !> increment variable activity counter
                            n_gas_eq=n_gas_eq+1 !> increment equilibrium gas counter
                            n_gas_eq_var_act=n_gas_eq_var_act+1 !> increment variable activity equilibrium gas counter
                        
                        !> Check if constant activity kinetic gas
                        else if (this%chem_syst%gas_phase%gases(gas_ind)%cst_act_flag.eqv. .true. .and. gas_ind>&
                            this%chem_syst%gas_phase%num_gases_eq) then
                            !> Gas has constant activity AND is kinetic gas (index > num_gases_eq)
                            num_gases_cst=num_gases_cst+1 !> increment constant activity counter
                            n_gas_kin=n_gas_kin+1 !> increment kinetic gas counter
                        else
                            !> Variable activity kinetic gas (default case)
                            num_gases_var=num_gases_var+1 !> increment variable activity counter
                            n_gas_kin=n_gas_kin+1 !> increment kinetic gas counter
                        end if
                        
                        !> (Commented out code - alternative categorization)
                        !if (gas_ind<=this%chem_syst%gas_phase%num_gases_eq) then
                        !    n_gas_eq=n_gas_eq+1
                        !else
                        !    n_gas_kin=n_gas_kin+1
                        !end if
                    
                    else !> gas not found
                        error stop "Gas not found in gas phase"
                            !> gas name does not exist in chemical system
                    end if !> end gas found check
                end if !> end sentinel check
            end do !> end gas species list loop
            
            if (num_gas_zones==ngtype) exit !> all zones read - exit first pass
        else !> str does not contain 'gas'
            exit !> no gas data - exit first pass
        end if !> end heading check
    end do !> end first pass outer loop

!> ================================================================
!> SECTION 5: Second pass - read full gas composition data
!> ================================================================
    
    rewind(unit) !> rewind file to beginning for second pass
    num_gas_zones=0 !> reset counter
    !num_gases_glob=1 !> (unused/commented)
    
    do !> search for gas zones section
        read(unit,*) label !> read section label
        if (label=='INITIAL AND BOUNDARY GAS ZONES') then !> found section
            read(unit,*) ngtype !> number of gas zones (redundant read)
            
            !> ============================================================
            !> Process each gas zone - read full data
            !> ============================================================
            do i=1,ngtype !> loop through all gas zones
                !> Read zone header (redundant but required by file format)
                read(unit,*) igtype, temp, vol !> index, temperature, volume
                read(unit,*) name !> zone name
                read(unit,*) str !> headings
                num_gas_zones=num_gas_zones+1 !> increment counter
                k=0 !> counter: constant activity equilibrium gases
                
                !> --------------------------------------------------------
                !> Read partial pressures for all gases in zone
                !> --------------------------------------------------------
                do j=1,this%gas_zones(igtype)%reactive_zone%gas_phase%num_species
                    read(unit,*) gas%name, part_press !> gas name and partial pressure (atm)
                    call this%gas_zones(igtype)%reactive_zone%gas_phase%gases(j)%copy_species(this%chem_syst%gas_phase%gases(&
                        ind_gases(j))) !> assign gas from chemical system
                    this%gas_zones(igtype)%activities(j)=part_press !> store partial pressure as activity (atm)
                    
                    if (ind_gases(j) <= this%chem_syst%gas_phase%num_gases_eq_cst_act) then
                        k=k+1 !> increment constant activity equilibrium counter
                        this%gas_zones(igtype)%ind_gases_eq_cst_act(k)=ind_gases(j) !> store index
                    end if
                end do
                
                !> --------------------------------------------------------
                !> Compute derived properties
                !> --------------------------------------------------------
                call this%gas_zones(igtype)%compute_conc_gases_ideal() !> compute concentrations (ideal gas law)
                call this%gas_zones(igtype)%compute_pressure() !> compute total pressure
                call this%gas_zones(igtype)%compute_log_act_coeffs_gases() !> compute log10(activity coefficients)
                call this%gas_zones(igtype)%set_indices_gases() !> set gas indices
                
                !> --------------------------------------------------------
                !> Check for duplicate gas phases
                !> --------------------------------------------------------
                if (num_gas_zones>1) then !> not first zone
                    call are_gas_phases_equal(this%gas_zones(igtype)%reactive_zone%gas_phase,this%gas_zones(&
                        ind_gas_zones(num_gas_zones-1))%reactive_zone%gas_phase,flag)
                        !> compare with previous zone gas phase
                    
                    if (flag .eqv. .true.) then !> identical gas phases
                        call this%gas_zones(igtype)%set_reactive_zone(this%gas_zones(ind_gas_zones(&
                            num_gas_zones-1))%reactive_zone)
                            !> reuse previous reactive zone
                    else !> different gas phases
                        !call this%gas_zones(igtype)%set_reactive_zone(react_zones(igtype))
                        ngrz=ngrz+1 !> increment gas reactive zones counter
                        call this%gas_zones(igtype)%reactive_zone%set_ind_non_flow_species()
                            !> set non-flowing species array
                        call this%gas_zones(igtype)%reactive_zone%set_ind_eq_reacts()
                        call this%gas_zones(igtype)%reactive_zone%allocate_ind_var_act_species()
                            !> set equilibrium reaction indices
                        call this%gas_zones(igtype)%reactive_zone%set_stoich_mat_react_zone()
                            !> set stoichiometric matrix
                        call this%gas_zones(igtype)%reactive_zone%set_ind_gases_stoich_mat()
                            !> set gas indices in stoichiometric matrix
                    end if
                else !> first zone
                    call this%gas_zones(igtype)%reactive_zone%set_ind_non_flow_species()
                        !> set non-flowing species array
                    call this%gas_zones(igtype)%reactive_zone%set_ind_eq_reacts()
                    call this%gas_zones(igtype)%reactive_zone%allocate_ind_var_act_species()
                        !> set equilibrium reaction indices
                    call this%gas_zones(igtype)%reactive_zone%set_stoich_mat_react_zone()
                        !> set stoichiometric matrix
                    call this%gas_zones(igtype)%reactive_zone%set_ind_gases_stoich_mat()
                        !> set gas indices in stoichiometric matrix
                end if
                
                read(unit,*) str !> read next line (may be sentinel or next zone)
            end do !> end gas zones loop
            
        else if (label=='end') then !> end of file marker
            backspace(unit) !> backspace to leave 'end' for next reader
            exit !> exit search loop
        else !> other section
            continue !> continue searching
        end if
    end do !> end second pass search loop

!> ================================================================
!> SECTION 6: Integrate gas zones with existing reactive zones
!> ================================================================
    
    !if (present(reactive_zones)) then !> (commented - always executes)
    
        !> ============================================================
        !> CASE 1: Existing reactive zones present (surface/cation exchange)
        !> ============================================================
        if (allocated(this%reactive_zones)) then !> reactive zones already exist
            !> Strategy: Create combinatorial reactive zones
            !> - First ngrz zones: gas-only reactive zones
            !> - Next num_surf_rz zones: existing surface/exchange zones (no gas)
            !> - Remaining zones: all combinations of gas zones × surface zones
            !> Total zones = ngrz + num_surf_rz + (ngrz × num_surf_rz)
            
            !> ------------------------------------------------------------
            !> Save existing surface/exchange reactive zones
            !> ------------------------------------------------------------
            num_surf_rz=size(this%reactive_zones) !> number of existing zones
            allocate(aux_react_zones(num_surf_rz)) !> temporary storage
            do i=1,num_surf_rz
                call aux_react_zones(i)%copy_react_zone(this%reactive_zones(i))
                    !> copy existing zones to temporary array
            end do
            
            !> ------------------------------------------------------------
            !> Compute total number of reactive zones
            !> ------------------------------------------------------------
            num_rz=num_surf_rz+ngtype*(1+num_surf_rz)
                !> num_surf_rz: surface-only zones
                !> ngtype: gas-only zones  
                !> ngtype*num_surf_rz: combined gas+surface zones
            
            !> ------------------------------------------------------------
            !> Reallocate reactive zones array with new size
            !> ------------------------------------------------------------
            !deallocate(this%reactive_zones) !> (commented - handled in allocate_reactive_zones)
            !allocate(this%reactive_zones(num_rz)) !> (commented - handled below)
            call this%allocate_reactive_zones(num_rz)
                !> allocate expanded reactive zones array
            
            !> ------------------------------------------------------------
            !> Restore surface/exchange zones at end of array
            !> ------------------------------------------------------------
            do i=1,num_surf_rz
                call this%reactive_zones(ngrz+i)%copy_react_zone(aux_react_zones(i))
                    !> place surface zones after gas zones
                    !> positions: ngrz+1 to ngrz+num_surf_rz
            end do
            
            !> ------------------------------------------------------------
            !> Set up gas-only reactive zones
            !> ------------------------------------------------------------
            do i=1,ngrz
                call this%reactive_zones(i)%set_chem_syst_react_zone(this%gas_zones(i)%reactive_zone%chem_syst)
                    !> link to chemical system
                call this%reactive_zones(i)%set_gas_phase(this%gas_zones(i)%reactive_zone%gas_phase)
                    !> set gas phase (no surface/exchange)
                !call this%gas_zones(i)%set_reactive_zone(this%reactive_zones(i))
                    !> (commented - circular reference issue)
            end do
            
            !> ------------------------------------------------------------
            !> Set up combined gas+surface reactive zones
            !> ------------------------------------------------------------
            do i=1,ngrz !> loop through gas zones
                do j=1,num_surf_rz !> loop through surface zones
                    !> Position in array: ngrz + num_surf_rz + (i-1)*num_surf_rz + j
                    !> Simplified: ngrz + i*num_surf_rz + j
                    call this%reactive_zones(ngrz+i*num_surf_rz+j)%set_chem_syst_react_zone(&
                        this%gas_zones(I)%reactive_zone%chem_syst)
                        !> link to chemical system (Note: capital I - possible bug?)
                    call this%reactive_zones(ngrz+i*num_surf_rz+j)%set_gas_phase(this%gas_zones(i)%reactive_zone%gas_phase)
                        !> set gas phase from gas zone i
                    call this%reactive_zones(ngrz+i*num_surf_rz+j)%set_cat_exch_zone(aux_react_zones(i)%cat_exch_zone)
                        !> set cation exchange zone from surface zone j
                end do
            end do
        
        !> ============================================================
        !> CASE 2: No existing reactive zones (gas-only simulation)
        !> ============================================================
        else !> no reactive zones allocated yet
            !> Simple case: create one reactive zone per gas zone
            !allocate(this%reactive_zones(ngrz)) !> (commented - handled below)
            call this%allocate_reactive_zones(ngrz)
                !> allocate reactive zones array
            do i=1,ngrz
                call this%reactive_zones(i)%copy_react_zone(this%gas_zones(i)%reactive_zone)
                    !> copy gas zone reactive zone to chemistry reactive zones array
                call this%gas_zones(i)%set_reactive_zone(this%reactive_zones(i))
                    !> link gas zone back to chemistry reactive zone (note: workaround)
            end do
        end if !> end reactive zones allocated check
    !end if !> end present(reactive_zones) check (commented - always executes)

end subroutine !> end read_gas_zones_CHEPROO