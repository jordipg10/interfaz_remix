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
!>     - icon=1: Fixed activity
!>     - icon=2: Fixed concentration
!>     - icon=3: Fixed pH (for H+)
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
!> @param[in] gas_species_chem Optional gas chemistry object (chapuza/workaround)
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
subroutine read_init_bd_wat_types_CHEPROO(this,unit,gas_chem)
!> ================================================================
!> MODULE IMPORTS
!> ================================================================
    use chemistry_m, only: chemistry_c !> main chemistry class
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use gas_chemistry_m, only: gas_chemistry_c
    use reactive_zone_m, only: reactive_zone_c
    use mineral_zone_m, only: mineral_zone_c
    use solid_chemistry_m, only: solid_chemistry_c
    use aq_species_m, only: aq_species_c
    
    implicit none !> no implicit variable typing
    
!> ================================================================
!> SUBROUTINE PARAMETERS
!> ================================================================
    class(chemistry_c) :: this !> chemistry object (contains all chemical system data)
    integer(kind=4), intent(in) :: unit !> file unit number for reading input file
    !type(solid_chemistry_c), intent(inout), allocatable :: init_cat_exch_zones(:) !> (COMMENTED) initial cation exchange zones
    type(gas_chemistry_c), intent(in), optional, target :: gas_chem !> optional gas chemistry object (chapuza/workaround)
    
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
    
!> ================================================================
!> LOCAL VARIABLES - Dynamic Arrays
!> ================================================================
    integer(kind=4), allocatable :: cols(:) !> column indices (size 2, purpose unclear)
    integer(kind=4), allocatable :: num_aq_prim_array(:) !> number of primary aqueous species per water type
    integer(kind=4), allocatable :: num_cstr_array(:) !> number of constraints (icon=4) per water type
    
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
    logical :: flag_surf !> surface flag (unused)
    logical :: flag_comp !> component flag (unused)
    logical :: flag_Se !> selenium flag (unused)
    
!> ================================================================
!> LOCAL VARIABLES - Chemistry Objects
!> ================================================================
    type(reactive_zone_c) :: react_zone !> default reactive zone object (unused)
    type(reactive_zone_c), allocatable :: react_zones(:) !> array of reactive zones (one per water type)
    !type(mineral_zone_c), allocatable :: min_zones(:) !> moved to this%min_zones_wat_types for persistence
    type(solid_chemistry_c) :: solid_chem !> default solid chemistry object (unused)
    !type(solid_chemistry_c), allocatable :: solid_chems(:) !> (COMMENTED) array of solid chemistry objects
    type(solid_chemistry_c), allocatable :: cat_exch_zones(:) !> cation exchange zones array
    type(aq_species_c) :: aq_species !> temporary aqueous species object for searching
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
!> SECTION 3: Handle cation exchange zones initialization
!> ================================================================
    !> Check if a single initial cation exchange zone exists
    if (this%num_init_cat_exch_zones==1) then !> if exactly one initial cat exch zone defined
        !> Replicate single cation exchange zone for all water types
        allocate(cat_exch_zones(nwtype)) !> allocate cat exch zones array (one per water type)
        
        do i=1,nwtype !> loop over all water types
            !> Copy initial cation exchange zone to each water type
            call cat_exch_zones(i)%copy_solid_chemistry(this%init_cat_exch_zones(1)) !> chapuza (copy zone 1)
        end do
        
        !deallocate(init_cat_exch_zones) !> (COMMENTED) deallocate original array
        
        !> Reallocate initial cat exch zones to match number of water types
        call this%allocate_init_cat_exch_zones(nwtype) !> allocate with new size = nwtype
    end if
    
!> ================================================================
!> SECTION 4: Allocate chemistry arrays for water types
!> ================================================================
    !> Allocate water types array
    call this%allocate_wat_types(nwtype) !> allocate water types array (size = nwtype)
    
    !> Allocate solid types associated with water types
    call this%allocate_wat_type_solids() !> allocate per-water-type solid chemistry templates (size = num_wat_types)
    
    !> Allocate reactive zones associated with water types
    call this%allocate_react_zones_wat_types() !> allocate reactive zones (chapuza - one per water type)
    
    !> Allocate gas zones associated with water types
    call this%allocate_gas_zones_wat_types() !> allocate gas zones (chapuza - one per water type)

    !> Allocate persistent mineral zones for water types
    call this%allocate_min_zones_wat_types() !> allocate mineral zones (one per water type, persistent)
    
    !> Allocate local temporary arrays
    allocate(react_zones(nwtype)) !> allocate reactive zones (size = nwtype)
    allocate(cols(2)) !> allocate column indices array (size 2, purpose unclear)
    allocate(num_aq_prim_array(nwtype),num_cstr_array(nwtype)) !> allocate counters for species and constraints
    
!> ================================================================
!> SECTION 5: Initialize counter arrays
!> ================================================================
    !> Initialize species and constraint counters to zero
    num_aq_prim_array=0 !> zero all elements (no primary species counted yet)
    num_cstr_array=0 !> zero all elements (no constraints counted yet)
    
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
        call this%wat_types(j)%set_temp(temp+273.18) !> T(K) = T(C) + 273.15 (approx 273.18)
        
        !> Calculate water density at this temperature
        call this%wat_types(j)%set_density() !> density = f(temperature)
        
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
                !> Read species name and constraint icon
                read(unit,*) aq_species%name, icon !> species name and icon value
                
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
                call this%react_zones_wat_types(j)%set_CV_params(this%CV_params)
                    !> copy CV parameters to reactive zone
                
                !> Set chemical system in reactive zone
                call this%react_zones_wat_types(j)%set_chem_syst_react_zone(this%chem_syst)
                    !> copy chemical system definition to reactive zone
                
                !> ------------------------------------------------------------
                !> Set up mineral zone for this water type
                !> ------------------------------------------------------------
                !> Set chemical system in mineral zone
                call this%min_zones_wat_types(j)%set_chem_syst_min_zone(this%chem_syst)
                    !> copy chemical system definition to mineral zone
                
                !> ------------------------------------------------------------
                !> Set up solid chemistry associations
                !> ------------------------------------------------------------
                !> Link reactive zone to solid type
                call this%wat_type_solids(j)%set_reactive_zone(this%react_zones_wat_types(j))
                    !> set default reactive zone for this solid type
                
                !call this%wat_type_solids(j)%set_mineral_zone(min_zones(j)) !> (COMMENTED) set default mineral zone
                
                !> ------------------------------------------------------------
                !> Handle gas chemistry (COMMENTED OUT)
                !> ------------------------------------------------------------
                ! if (present(gas_species_chem)) then !> if optional gas chemistry provided
                !    call this%wat_types(j)%set_gas_chemistry(gas_species_chem) !> chapuza (set gas chemistry)
                ! end if
                
                !> Link solid chemistry to water type
                call this%wat_types(j)%set_solid_chemistry(this%wat_type_solids(j)) !> chapuza (set solid type)
                
                !> ------------------------------------------------------------
                !> Handle cation exchange zones if present
                !> ------------------------------------------------------------
                if (this%num_init_cat_exch_zones>0) then !> if cation exchange zones exist
                    !> Assign cation exchange zone to solid type
                    call this%wat_type_solids(j)%copy_solid_chemistry(cat_exch_zones(j))
                        !> chapuza (copy cat exch zone to solid type)
                    
                    !> --------------------------------------------------------
                    !> Set mineral zone (conditional on existence)
                    !> --------------------------------------------------------
                    if (this%num_mineral_zones==0) then !> if no mineral zones defined
                        !> Use dummy mineral zone
                        call this%wat_type_solids(j)%set_mineral_zone(this%min_zone_dummy)
                            !> chapuza (use empty/dummy mineral zone)
                    else !> mineral zones exist
                        !> Use actual mineral zone
                        call this%wat_type_solids(j)%set_mineral_zone(this%min_zones_wat_types(j))
                            !> chapuza (set mineral zone for this water type)
                    end if
                    
                    !> --------------------------------------------------------
                    !> Read detailed water type composition with cat exch
                    !> --------------------------------------------------------
                    call this%wat_types(j)%read_wat_type_CHEPROO(num_aq_prim_array(j),num_cstr_array(j),this%num_gas_zones,&
                        this%act_coeffs_model,&
                        this%Jac_opt,unit,niter,CV_flag)
                        !> read full water type composition (species conc, constraints, etc.)
                        !> Parameters: num primary species, num constraints, num gas zones, act coeff model, Jacobian option
                        !> Returns: niter (iterations), CV_flag (control volume flag)
                    
                    !> Assign updated solid chemistry back to initial cat exch zones
                    call this%init_cat_exch_zones(j)%copy_solid_chemistry(this%wat_type_solids(j))
                        !> chapuza (copy solid type back to init cat exch zone)
                
                else !> no cation exchange zones
                    !> --------------------------------------------------------
                    !> Set mineral zone without cation exchange
                    !> --------------------------------------------------------
                    call this%wat_type_solids(j)%set_mineral_zone(this%min_zones_wat_types(j))
                        !> set mineral zone for this solid type
                    
                    !> --------------------------------------------------------
                    !> Read water type composition (conditional on gas chemistry)
                    !> --------------------------------------------------------
                    if (present(gas_chem)) then !> if optional gas chemistry provided
                        !> Read water type with gas chemistry
                        call this%wat_types(j)%read_wat_type_CHEPROO(num_aq_prim_array(j),num_cstr_array(j),this%num_gas_zones,&
                            this%act_coeffs_model,&
                            this%Jac_opt,unit,niter,CV_flag,gas_chem)
                            !> read with optional gas_species_chem argument
                    else !> no gas chemistry
                        !call this%wat_type_solids(j)%set_mineral_zone(min_zones(j)) !> (COMMENTED) redundant mineral zone set
                        
                        !> Read water type without gas chemistry
                        call this%wat_types(j)%read_wat_type_CHEPROO(num_aq_prim_array(j),num_cstr_array(j),this%num_gas_zones,&
                            this%act_coeffs_model,this%Jac_opt,unit,niter,CV_flag)
                            !> read without optional gas_species_chem argument
                    end if
                end if !> end if num_init_cat_exch_zones > 0
            end do !> end second pass water types loop
            
            !> Exit search loop after processing section
            exit !> exit file search loop
            
         else !> label is not 'INITIAL AND BOUNDARY WATER TYPES'
            !> Continue searching for section header
            continue !> read next line
         end if
    end do !> end file search loop
    
!> ================================================================
!> SECTION 9: Post-processing - eliminate constant activity species
!> ================================================================
    !> Eliminate constant activity species from chemical system component matrix
    call this%chem_syst%change_spec_alg_chem_syst(this%CV_params%zero)
        !> remove species with constant activity (activity coefficients ~= 1) from speciation algebra
        !> modifies component matrix and speciation algebra
        !> zero threshold from CV_params (control volume parameters)
    
    !> Rearrange mineral zone indices after speciation algebra change
    do i=1,nwtype
        call this%wat_type_solids(i)%mineral_zone%set_ind_min_Sk()
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