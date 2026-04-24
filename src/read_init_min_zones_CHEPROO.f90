!> \file read_init_min_zones_CHEPROO.f90
!> \brief Read initial mineral zones from CHEPROO format input file
!> \details
!>   Parses CHEPROO-formatted input to initialize mineral zones for reactive transport.
!>   
!>   **Key assumptions:**
!>   - All minerals are pure phases (activity = 1)
!>   - Each mineral zone has a unique set of minerals
!>   - Temperature specified in Celsius, converted to Kelvin
!>   
!>   **Input file format:**
!>   ```
!>   INITIAL MINERAL ZONES
!>   <nmtype>           ! Number of mineral zones
!>   <zone_id> <temp>   ! Zone index, temperature (°C)
!>   <name>             ! Zone name
!>   mineral            ! Header line
!>   <min_name> <vol_frac> <react_surf>  ! Mineral name, volume fraction [-], reactive surface [m²/m³]
!>   ...
!>   *                  ! End of zone marker
!>   ```
!>   
!>   **Two-pass algorithm:**
!>   1. First pass: Count minerals, classify as kinetic/equilibrium and constant/variable activity
!>   2. Second pass: Read properties (volume fractions, reactive surfaces)
!>   
!>   **Reactive zone creation:**
!>   - If gas zones exist: creates combined solid-gas reactive zones
!>   - Otherwise: creates pure solid reactive zones
!>   
!> \author jordi Petchamé-Guerrero
!> \date October 2025

!> \brief Read initial mineral zones from CHEPROO format file
!> \details
!>   Two-pass reading algorithm:
!>   
!>   **First pass:**
!>   - Counts minerals per zone
!>   - Classifies minerals (kinetic vs equilibrium, constant vs variable activity)
!>   - Allocates arrays and initializes zone structures
!>   
!>   **Second pass:**
!>   - Reads mineral properties (volume fractions, reactive surfaces)
!>   - Sets up stoichiometric matrices
!>   - Creates combined reactive zones if gas zones exist
!>   
!>   **Mineral classification:**
!>   - Kinetic minerals: index <= num_minerals_kin in chemical system
!>   - Equilibrium minerals: index > num_minerals_kin
!>   - Constant activity: cst_act_flag = TRUE
!>   - Variable activity: cst_act_flag = FALSE
!>   
!> \param[in,out] this Chemistry object to populate with mineral zones
!> \param[in] unit File unit number (must be open for reading)
!> \param[out] nmrz Number of mineral reactive zones created
!> \param[in] surf_chem Optional: surface chemistry for cation exchange coupling
subroutine read_init_min_zones_CHEPROO(this,unit,nmrz,surf_chem)
!> We assume all minerals are pure phases
    use chemistry_m, only: chemistry_c
    use solid_chemistry_m, only: solid_chemistry_c
    use mineral_zone_m, only: mineral_zone_c
    use mineral_m, only: mineral_c
    use reactive_zone_m, only: reactive_zone_c
    implicit none
    class(chemistry_c) :: this !> chemistry object
    integer(kind=4), intent(in) :: unit !> file unit number for reading CHEPROO input
    !type(solid_chemistry_c), intent(out), allocatable :: materials(:) !> initial mineral zones
    integer(kind=4), intent(out) :: nmrz !> number of mineral reactive zones ( we assume each mineral zone has different minerals )
    type(solid_chemistry_c), intent(in), optional :: surf_chem  !< Optional surface chemistry for coupling with cation exchange
    
    !< \name Local variables for mineral classification and counting
    !< @{
    integer(kind=4) :: num_mins_cst_kin   !< Number of constant activity kinetic minerals in current zone
    integer(kind=4) :: num_mins_var_kin   !< Number of variable activity kinetic minerals in current zone
    integer(kind=4) :: num_mins_loc_kin   !< Total kinetic minerals in current zone (local counter)
    integer(kind=4) :: num_surf_rz        !< Number of surface reactive zones
    integer(kind=4) :: num_gas_rz         !< Number of gas-only reactive zones
    integer(kind=4) :: num_rz_old         !< Previous number of reactive zones (before adding mineral zones)
    integer(kind=4) :: i,j,k              !< Loop indices for iterations and array traversals
    integer(kind=4) :: nmtype             !< Number of mineral zone types read from file
    integer(kind=4) :: nrwtype            !< Number of reactive water types (unused)
    integer(kind=4) :: icon               !< Icon/flag variable (unused)
    integer(kind=4) :: num_min_zones      !< Counter for mineral zones processed
    integer(kind=4) :: num_mins_rz        !< Number of minerals in reactive zone (unused)
    integer(kind=4) :: num_mins_glob      !< Global mineral count (unused)
    integer(kind=4) :: num_mins_loc       !< Number of minerals in current local zone
    integer(kind=4) :: num_mins_loc_eq    !< Number of equilibrium minerals in current local zone
    integer(kind=4) :: ind_rz             !< Reactive zone index (unused)
    integer(kind=4) :: imtype             !< Index of current mineral zone type being processed
    integer(kind=4) :: min_ind            !< Index of mineral in chemical system
    integer(kind=4) :: num_mins_vare_mean !< Number of variable activity equilibrium minerals
    integer(kind=4) :: num_mins_cst_eq    !< Number of constant activity equilibrium minerals
    integer(kind=4) :: num_rz             !< Total number of reactive zones (including combined zones)
    !< @}
    
    !< \name Arrays for mineral indexing
    !< @{
    integer(kind=4), allocatable :: init_min_indices(:)  !< Indices of initial minerals (unused)
    integer(kind=4), allocatable :: min_indices(:,:)     !< Mineral indices: (mineral_num, zone_num) - maps local to global indices
    !< @}
    
    !< \name String variables for parsing
    !< @{
    character(len=256) :: str        !< General string for reading headers and markers
    character(len=256) :: constrain  !< Constraint string (unused)
    character(len=256) :: label      !< Label for section identification
    character(len=256) :: min_name   !< Mineral name (unused, using mineral%name instead)
    character(len=256) :: name       !< Zone name read from file
    !< @}
    
    !< \name Real variables for properties
    !< @{
    real(kind=8) :: guess       !< Initial guess for solver (unused)
    real(kind=8) :: c_tot       !< Total concentration (unused)
    real(kind=8) :: temp        !< Temperature in Celsius (converted to Kelvin)
    real(kind=8) :: vol_frac    !< Volume fraction of mineral [-]
    real(kind=8) :: react_surf  !< Reactive surface area [m² mineral / m³ rock]
    real(kind=8) :: conc        !< Mineral concentration [M] (unused, set to 1.0 for pure phases)
    !< @}
    
    !< \name Logical flags
    !< @{
    logical :: min_flag   !< TRUE if mineral found in chemical system, FALSE otherwise
    logical :: flag_comp  !< TRUE to compute speciation algebra dimensions (component matrix flag)
    !< @}
    
    !< \name Type variables for zone construction
    !< @{
    type(mineral_c) :: mineral                           !< Temporary mineral object for reading
    type(reactive_zone_c) :: aux_react_zone              !< Auxiliary reactive zone (unused)
    type(mineral_zone_c), allocatable :: min_zones(:)    !< Array of mineral zones (unused)
    type(reactive_zone_c), allocatable :: react_zones(:) !< Array of reactive zones for this set
    type(reactive_zone_c), allocatable :: aux_react_zones(:) !< Auxiliary array for zone combination
    !< @}
    
    flag_comp=.true. !> flag to compute speciation algebra dimensions (by default: TRUE = no cst act species in component matrix)

    !< ========================================================================
    !< \subsection first_pass_header First Pass: Read number of mineral zones and validate
    !< ========================================================================
    read(unit,*) nmtype  !< Read number of mineral zone types from file
    
    if (nmtype<0) then  !< Validate positive count
        error stop "Number of mineral zones must be positive"
    else if (nmtype==0) then  !< No mineral zones to process
        nmrz=0 !> number of mineral reactive zones (set to zero)
        return  !< Exit subroutine early
    else  !< Process mineral zones
        !< ====================================================================
        !< \subsection first_pass_allocation Allocate arrays for zones
        !< ====================================================================
        !nmrz=1 !> counter mineral reactive zones
        allocate(react_zones(nmtype)) !> chapuza: allocate reactive zones array
        call this%allocate_mineral_zones(nmtype)  !< Allocate mineral zones in chemistry object
        allocate(min_indices(this%chem_syst%num_minerals,nmtype))  !< Allocate indexing array [mineral, zone]
        nmrz=nmtype !> number of mineral reactive zones equals zone types
        call this%set_num_materials(nmtype)  !< Set total number of materials
        call this%allocate_materials()  !< Allocate materials array in chemistry object
        
            !allocate(this%materials(nmtype))
            !< ================================================================
            !< \subsection first_pass_loop First pass loop: Count and classify minerals
            !< ================================================================
            num_min_zones=0 !> counter mineral zones processed
            do  !< Outer loop: iterate until all zones processed
                read(unit,*) imtype, temp !> index of mineral zone (1-based), temperature in Celsius
                if (imtype<1 .or. imtype>nmtype) error stop "Mineral zone index out of range"  !< Validate zone index
                read(unit,*) name !> name of mineral zone (for documentation)
                read(unit,*) str !> headings line (expecting "mineral" keyword or "*" end marker)
                    if (str=='*') then  !< End marker "*" encountered before mineral data
                        if (num_min_zones<nmtype) then  !< Not all zones processed yet
                            continue  !< Skip to next iteration
                        else  !< All zones processed
                            exit  !< Exit outer loop
                        end if
                    else if (index(str,'mineral')/=0) then  !< "mineral" keyword found - start of mineral list
                        !< Initialize counters for this zone
                        num_min_zones=num_min_zones+1  !< Increment zone counter
                        num_mins_loc=0 !> counter minerals in this zone (all types)
                        num_mins_loc_eq=0 !> counter equilibrium minerals in this zone
                        num_mins_loc_kin=0 !> counter kinetic minerals in this zone
                        num_mins_cst_eq=0  !< Constant activity equilibrium minerals
                        num_mins_cst_kin=0  !< Constant activity kinetic minerals
                        num_mins_vare_mean=0  !< Variable activity equilibrium minerals
                        num_mins_var_kin=0  !< Variable activity kinetic minerals
                        do  !< Inner loop: read mineral names and classify
                            read(unit,*) mineral%name  !< Read mineral name from file
                            if (mineral%name=='*') then  !< End marker "*" - finish this zone
                                call this%mineral_zones(imtype)%set_chem_syst_min_zone(this%chem_syst)
                                call this%mineral_zones(imtype)%set_id(imtype)
                                call this%materials(imtype)%set_mineral_zone(this%mineral_zones(imtype)) !> chapuza
                                call react_zones(imtype)%set_chem_syst_react_zone(this%chem_syst)
                                call this%materials(imtype)%mineral_zone%set_chem_syst_min_zone(this%chem_syst)
                                call react_zones(imtype)%set_CV_params(this%CV_params)
                                call react_zones(imtype)%allocate_ind_mins(num_mins_loc_eq)
                                call this%materials(imtype)%mineral_zone%allocate_ind_chem_syst_min_zone(num_mins_loc)
                                call this%materials(imtype)%mineral_zone%allocate_ind_min_Sk(num_mins_loc_kin)
                                call this%materials(imtype)%mineral_zone%set_num_mins_eq_min_zone(num_mins_loc_eq)
                                call this%materials(imtype)%mineral_zone%set_num_mins_eq_cst_act_min_zone(num_mins_cst_eq)
                                call this%materials(imtype)%mineral_zone%set_num_mins_kin_cst_act_min_zone(num_mins_cst_kin)
                                call this%materials(imtype)%mineral_zone%set_num_mins_eq_var_act_min_zone(num_mins_vare_mean)
                                !call react_zones(imtype)%allocate_ind_non_flow_species()
                                if (present(surf_chem)) then
                                    call react_zones(imtype)%set_cat_exch_zone(surf_chem%reactive_zone%cat_exch_zone)
                                    call this%materials(imtype)%set_CEC(surf_chem%CEC)
                                else
                                    !call react_zones(imtype)%set_surf_chem(this%surf_chem)
                                    !call this%materials(imtype)%set_surf_chem(this%surf_chem)
                                end if
                                !call react_zones(imtype)%set_num_mins(num_mins_loc)
                                call react_zones(imtype)%set_num_solids()
                                call react_zones(imtype)%set_num_mins_cst_act(num_mins_cst_eq)
                                !call this%materials(imtype)%mineral_zone%set_num_mins_cst_act_min_zone(num_mins_cst_eq+num_mins_cst_kin)
                                call react_zones(imtype)%set_num_mins_var_act(num_mins_vare_mean)
                                !call this%materials(imtype)%mineral_zone%set_num_mins_var_act_min_zone(num_mins_vare_mean+num_mins_var_kin)
                                call react_zones(imtype)%set_num_non_flow_species() !> sets number of non-flowing species in reactive zone
                                call react_zones(imtype)%set_speciation_alg_dimensions(flag_comp)
                                !> aqui habria que comprobar si se estan repitiendo zonas reactivas
                                call this%materials(imtype)%set_reactive_zone(react_zones(imtype))
                                call this%materials(imtype)%compute_num_solids_solid_chem()
                                call this%materials(imtype)%allocate_vol_fracts_mins()
                                call this%materials(imtype)%allocate_react_surfaces()
                                call this%materials(imtype)%allocate_conc_solids()
                                call this%materials(imtype)%allocate_log_act_coeffs_solid_chem()
                                call this%materials(imtype)%allocate_activities()
                                call this%materials(imtype)%set_temp(temp+273.18) !> Kelvin
                                call this%materials(imtype)%allocate_var_act_species_indices(num_mins_vare_mean+num_mins_var_kin)
                                call this%materials(imtype)%allocate_cst_act_species_indices(num_mins_cst_eq+num_mins_cst_kin)
                                call this%materials(imtype)%set_indices_solids()
                                exit
                            else
                                !> We check if mineral exists in the chemical system
                                call this%chem_syst%is_mineral_in_chem_syst(mineral,min_flag,min_ind)
                                if (min_flag .eqv. .true.) then
                                    num_mins_loc=num_mins_loc+1
                                    min_indices(num_mins_loc,num_min_zones)=min_ind !> we save index of mineral in chemical system
                                    if (this%chem_syst%minerals(min_ind)%mineral%cst_act_flag.eqv..true. .and. &
                                    min_ind>this%chem_syst%num_minerals_kin) then
                                        num_mins_cst_eq=num_mins_cst_eq+1
                                        num_mins_loc_eq=num_mins_loc_eq+1
                                    else if (this%chem_syst%minerals(min_ind)%mineral%cst_act_flag.eqv..false. .and. &
                                        min_ind>this%chem_syst%num_minerals_kin) then
                                        num_mins_vare_mean=num_mins_vare_mean+1
                                        num_mins_loc_eq=num_mins_loc_eq+1
                                    else if (this%chem_syst%minerals(min_ind)%mineral%cst_act_flag.eqv..true. .and. &
                                        min_ind<=this%chem_syst%num_minerals_kin) then
                                        num_mins_cst_kin=num_mins_cst_kin+1
                                        num_mins_loc_kin=num_mins_loc_kin+1
                                    else
                                        num_mins_var_kin=num_mins_var_kin+1
                                        num_mins_loc_kin=num_mins_loc_kin+1
                                    end if
                                else
                                    error stop "Mineral not found in chemical system"
                                end if
                            end if
                        end do
                    else
                        exit
                    end if
                !end do
                if (num_min_zones==nmtype) exit  !< All zones processed, exit outer loop
            end do  !< End first pass loop
            rewind(unit)  !< Reset file to beginning for second pass
            
        !< ====================================================================
        !< \subsection second_pass Second Pass: Read mineral properties
        !< ====================================================================
            num_min_zones=0 !> counter mineral zones (reset for second pass)
            do  !< Loop to find "INITIAL MINERAL ZONES" section
                read(unit,*) label  !< Read label from file
                if (label=='INITIAL MINERAL ZONES') then  !< Found correct section
                    read(unit,*) nmtype !> number of mineral zones (re-read for confirmation)
                    do  !< Loop through zones in second pass
                        read(unit,*) imtype, temp !> index of mineral zone, temperature (C) - re-read
                        if (imtype<1 .or. imtype>nmtype) error stop "Mineral zone index out of range"  !< Validate
                        read(unit,*) name !> name of mineral zone (re-read)
                        read(unit,*) str !> headings (re-read)
                        if (index(str,'mineral')/=0) then  !< "mineral" keyword found
                            num_min_zones=num_min_zones+1  !< Increment zone counter
                            num_mins_loc=0 !> counter minerals in this zone (reset)
                            num_mins_loc_kin=0 !> counter kinetic minerals in this zone (reset)
                            num_mins_loc_eq=0 !> counter equilibrium minerals in this zone (reset)
                            do  !< Loop through minerals in this zone
                                read(unit,*) mineral%name, vol_frac, react_surf !> Mineral name, volume fraction [-], reactive surface [m² mineral/m³ rock]
                                if (mineral%name=='*') then  !< End marker for this zone
                                    exit  !< Exit mineral loop
                                else  !< Valid mineral data
                                    num_mins_loc=num_mins_loc+1  !< Increment local mineral counter
                                    if (min_indices(num_mins_loc,num_min_zones)>this%chem_syst%num_minerals_kin) then !> Equilibrium mineral (index > num_minerals_kin)
                                        num_mins_loc_eq=num_mins_loc_eq+1 !> Increment equilibrium mineral counter
                                        !< Set index in mineral zone (equilibrium minerals placed after kinetic)
                                        this%materials(imtype)%mineral_zone%ind_min_chem_syst(&
                                            this%materials(imtype)%mineral_zone%num_minerals_kin+num_mins_loc_eq)=&
                                            min_indices(num_mins_loc,num_min_zones) !> Global index from chemical system
                                        !< Set index in reactive zone
                                        this%materials(imtype)%reactive_zone%ind_mins_chem_syst(num_mins_loc_eq)=&
                                            min_indices(num_mins_loc,num_min_zones) !> Same global index
                                        !< Set volume fraction [-] (dimensionless)
                                        this%materials(imtype)%vol_fracts_mins(&
                                            this%materials(imtype)%mineral_zone%num_minerals_kin+num_mins_loc_eq)=vol_frac
                                        !< Set reactive surface [m² mineral/m³ rock]
                                        this%materials(imtype)%react_surfaces(&
                                            this%materials(imtype)%mineral_zone%num_minerals_kin+num_mins_loc_eq)=react_surf
                                        !< Set concentration [M] (pure phase assumption: c = 1.0 M)
                                        this%materials(imtype)%concentrations(&
                                            this%materials(imtype)%mineral_zone%num_minerals_kin+num_mins_loc_eq)=1d0
                                        !< Set activity [-] (pure phase assumption: a = 1.0)
                                        this%materials(imtype)%activities(&
                                        this%materials(imtype)%mineral_zone%num_minerals_kin+num_mins_loc_eq)=1d0
                                    else !> Kinetic mineral (index <= num_minerals_kin)
                                        num_mins_loc_kin=num_mins_loc_kin+1 !> Increment kinetic mineral counter
                                        !< Set index in mineral zone (kinetic minerals come first)
                                        this%materials(imtype)%mineral_zone%ind_min_chem_syst(num_mins_loc_kin)=&
                                            min_indices(num_mins_loc,num_min_zones) !> Global index from chemical system
                                        !< Set volume fraction [-]
                                        this%materials(imtype)%vol_fracts_mins(num_mins_loc_kin)=vol_frac
                                        !< Set reactive surface [m² mineral/m³ rock]
                                        this%materials(imtype)%react_surfaces(num_mins_loc_kin)=react_surf
                                        !< Set concentration [M] (pure phase assumption)
                                        this%materials(imtype)%concentrations(num_mins_loc_kin)=1d0
                                        !< Set activity [-] (pure phase assumption)
                                        this%materials(imtype)%activities(num_mins_loc_kin)=1d0
                                    end if
                                end if
                            end do  !< End mineral loop for this zone
                            
                            !! \subsection setup_stoichiometry Set up stoichiometric matrices and reaction indices
                            !< Set non-flowing species indices
                            call this%materials(imtype)%reactive_zone%set_ind_non_flow_species()
                            !< Set equilibrium reaction indices
                            call this%materials(imtype)%reactive_zone%set_ind_eq_reacts()
                            call this%materials(imtype)%reactive_zone%allocate_ind_var_act_species()
                            !< Set stoichiometric matrix for reactive zone
                            call this%materials(imtype)%reactive_zone%set_stoich_mat_react_zone()
                            !< Set mineral indices in stoichiometric matrix
                            call this%materials(imtype)%reactive_zone%set_ind_mins_stoich_mat()
                            !< Set gas indices in stoichiometric matrix
                            call this%materials(imtype)%reactive_zone%set_ind_gases_stoich_mat()
                            !< Set mineral indices in mineral zone stoichiometric matrix
                            call this%materials(imtype)%mineral_zone%set_ind_min_Sk()
                            
                            if (num_min_zones==nmtype) then  !< Last zone processed
                                exit  !< Exit zone loop
                            end if
                        else  !< No more zones to process
                            exit  !< Exit loop
                        end if
                    end do  !< End zone processing loop
                    exit  !< Exit outer loop
                else  !< Not at "INITIAL MINERAL ZONES" label
                    continue  !< Keep searching
                end if
            end do  !< End search loop
            
                !! \subsection combine_zones Combine mineral zones with existing reactive zones
                if (allocated(this%reactive_zones)) then  !< Existing reactive zones present
                    num_rz_old=size(this%reactive_zones)  !< Save old number of reactive zones
                    allocate(aux_react_zones(num_rz_old))  !< Temporary storage for existing zones
                    num_gas_rz=0  !< Counter for gas-only reactive zones
                    do i=1,num_rz_old  !< Loop through existing reactive zones
                        !< Copy existing reactive zone to temporary storage
                        call aux_react_zones(i)%copy_react_zone(this%reactive_zones(i))
                        !< Check if zone has gases but no solids (pure gas zone)
                        if (this%reactive_zones(i)%gas_phase%num_gases_eq>0 .and. this%reactive_zones(i)%num_solids==0) then
                            num_gas_rz=num_gas_rz+1  !< Increment gas zone counter
                        end if
                    end do  !< End loop through existing zones
                    !num_surf_rz=(num_rz_old-num_gas_rz)/(1+num_gas_rz) !> Commented out calculation
                    !< Calculate total number of reactive zones
                    num_rz=num_gas_rz+nmtype*(1+num_gas_rz) !> Pure gas zones + mineral zones × (1 + gas zones)
                    !< Allocate new reactive zones array
                    call this%allocate_reactive_zones(num_rz)
                    
                    !< Copy pure gas zones from auxiliary storage
                    do i=1,num_gas_rz
                        call this%reactive_zones(i)%copy_react_zone(aux_react_zones(i))
                        call this%gas_zones(i)%set_reactive_zone(this%reactive_zones(i)) !> Set gas zone (temporary fix)
                    end do
                    
                    !< Copy mineral-only reactive zones
                    do i=1,nmtype
                        call this%reactive_zones(num_gas_rz+i)%copy_react_zone(this%materials(i)%reactive_zone)
                        call this%materials(i)%set_reactive_zone(this%reactive_zones(num_gas_rz+i)) !> Set material zone (temporary fix)
                    end do
                    
                    !< Create combined gas + mineral reactive zones
                    do i=1,num_gas_rz  !< Loop through gas zones
                        !do j=1,num_surf_rz !> Old surface zone logic (commented out)
                        !    call this%reactive_zones(num_gas_rz+nmtype+i*num_surf_rz+j)%copy_react_zone(aux_react_zones(&
                        !    num_gas_rz+num_surf_rz+(i-1)*num_surf_rz+j))
                        !end do
                        do j=1,nmtype  !< Loop through mineral zones
                            !< Assign mineral reactive zone to combined zone
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%copy_react_zone(&
                                this%materials(j)%reactive_zone)
                            !< Add gas phase to combined zone
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_gas_phase(&
                                this%reactive_zones(i)%gas_phase)
                            !call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_gas_phase(& !> Duplicate (commented)
                            !    this%reactive_zones(i)%gas_phase)
                            !< Update non-flowing species for combined zone
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_ind_non_flow_species()
                            !< Update speciation algebra dimensions
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_speciation_alg_dimensions()
                            !< Update equilibrium reaction indices
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_ind_eq_reacts()
                            !> Update variable activity species indices
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%allocate_ind_var_act_species()
                            !< Update stoichiometric matrix
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_stoich_mat_react_zone()
                            !< Update mineral indices in stoichiometric matrix
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_ind_mins_stoich_mat()
                            !< Update gas indices in stoichiometric matrix
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_ind_gases_stoich_mat()
                        end do  !< End mineral zone loop
                    end do  !< End gas zone loop
                else  !< No existing reactive zones
                    !< Create new reactive zones array with only mineral zones
                    call this%allocate_reactive_zones(nmtype)
                    do i=1,nmtype  !< Loop through mineral zones
                        !< Assign mineral reactive zone to reactive zones array
                        call this%reactive_zones(i)%copy_react_zone(this%materials(i)%reactive_zone)
                        !< Set material reactive zone pointer (temporary fix)
                        call this%materials(i)%set_reactive_zone(this%reactive_zones(i))
                    end do  !< End mineral zone loop
                end if  !< End reactive zone combination logic
    end if  !< End file unit validity check
end subroutine read_init_min_zones_CHEPROO  !< End subroutine