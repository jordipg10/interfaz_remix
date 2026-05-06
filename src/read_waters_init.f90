!> @file read_waters_init.f90
!> @brief Reads initial target waters and their associated multiphase chemistry from file
!> @details This subroutine reads target water initialization data from a formatted file,
!>          establishing the initial chemical state for reactive transport simulations.
!>          It handles complex scenarios including:
!>          - Domain (internal) waters, boundary waters, and recharge waters
!>          - Individual target waters or ranges of target waters
!>          - Linking waters to solid chemistry (minerals) and gas chemistry
!>          - Reactive zone detection and assignment
!>          - Species index management for numerical stability
!>          - Stoichiometric matrix computation for kinetic reactions
!>
!>          File format expects labels and data in specific order:
!>          - "WATERS" label followed by counts and associations
!>          - Each line: water_index water_type solid_index gas_index flag
!>          - Supports range notation (e.g., "1-10" for waters 1 through 10)
!>          - Flag: 0=boundary, 1=domain, 2=recharge
!>
!>          The routine builds arrays of waters, target_solids, and target_gases,
!>          linking them appropriately and computing necessary speciation algorithms.
!>
!>          ALGORITHM OVERVIEW:
!>          1. Initialize water type counts and gas zone counts
!>          2. For systems without materials, create default reactive/mineral zones
!>          3. Open input file and read "WATERS" section
!>          4. Allocate target waters, boundary/recharge indices arrays
!>          5. For each line in file:
!>             a. Parse water/solid/gas indices (single or range)
!>             b. Validate indices and water type flag
!>             c. Assign indices to boundary/domain/recharge arrays
!>             d. For each water in range:
!>                - Call read_tar_wat_line to initialize water with multiphase chemistry
!>                - Detect reactive zone changes for optimization
!>                - Compute speciation algorithms if new zone
!>                - Handle species swapping for numerical stability
!>          6. Close file and set external waters indices
!>
!> @param[in,out] this Chemistry object to populate with target waters
!> @param[in] root Root filename (path without extension) for input file
!> @param[in] nsrz Number of solid reactive zones in the system
!> @param[in] ngrz Number of gas reactive zones in the system
!>
!> @note File is opened and closed within this subroutine (unit 57)
!> @warning Assumes water types, solid types, and gas types are already initialized
!> @warning File format is rigid - incorrect format will cause read errors or crashes
!> @see read_tar_wat_line for detailed water initialization logic
!>
subroutine read_waters_init(this,root,nsrz,ngrz,num_tar)
    use chemistry_m, only: chemistry_c
    use solid_chemistry_m, only: solid_chemistry_c
    use reactive_zone_m, only: reactive_zone_c
    implicit none
    class(chemistry_c) :: this !> chemistry object to populate
    character(len=*), intent(in) :: root !> root filename (without extension) for input files
    integer(kind=4), intent(in) :: nsrz !> number of solid reactive zones
    integer(kind=4), intent(in) :: ngrz !> number of gas reactive zones
    integer(kind=4), intent(in) :: num_tar !> expected number of target (domain) waters from the mesh
    
    !> Loop counters and index variables
    integer(kind=4) :: ind_rech      !> counter for recharge waters processed
    integer(kind=4) :: ind_dom       !> counter for domain (internal) waters processed
    integer(kind=4) :: counter_swap  !> counter for species index swaps (unused)
    integer(kind=4) :: ind           !> general purpose index variable
    integer(kind=4) :: i,j,k,m       !> loop indices
    
    !> Water type and count variables
    integer(kind=4) :: nwtype        !> number of water types in system
    integer(kind=4) :: num_wat   !> total number of target waters to read
    integer(kind=4) :: tar_wat_ind   !> current target water index being processed
    integer(kind=4) :: wtype         !> water type index for current target water
    
    !> Solid chemistry variables
    integer(kind=4) :: istype        !> solid type index (deprecated/unused)
    integer(kind=4) :: nstype        !> number of solid types (deprecated/unused)
    integer(kind=4) :: tar_sol_ind   !> target solid index for current water
    integer(kind=4) :: iszn          !> solid zone index for current water
    integer(kind=4) :: num_sol_types !> number of solid types (unused)
    integer(kind=4) :: first_sol     !> first solid index in a range
    integer(kind=4) :: last_sol      !> last solid index in a range
    integer(kind=4) :: int_sol_size  !> size of solid index range
    integer(kind=4) :: ind_bar_sol   !> position of '-' character in solid range string
    
    !> Gas chemistry variables
    integer(kind=4) :: ngzns         !> number of gas zones
    integer(kind=4) :: tar_gas_ind   !> target gas index for current water
    integer(kind=4) :: igzn          !> gas zone index for current water
    integer(kind=4) :: num_gas_types !> number of gas types (unused)
    integer(kind=4) :: first_gas     !> first gas index in a range
    integer(kind=4) :: last_gas      !> last gas index in a range
    integer(kind=4) :: int_gas_size  !> size of gas index range
    integer(kind=4) :: ind_bar_gas   !> position of '-' character in gas range string
    
    !> Boundary and recharge water variables
    integer(kind=4) :: nbwtype       !> number of boundary water types (unused)
    integer(kind=4) :: bwtype        !> boundary water type index (unused)
    integer(kind=4) :: num_wat_rech !> number of recharge target waters
    integer(kind=4) :: num_wat_bd   !> number of boundary target waters
    integer(kind=4) :: ind_bd        !> counter for boundary waters processed
    integer(kind=4) :: ind_bar_wat   !> position of '-' character in water range string
    
    !> Water range parsing variables
    integer(kind=4) :: first_wat     !> first water index in a range
    integer(kind=4) :: last_wat      !> last water index in a range
    integer(kind=4) :: int_wat_size  !> size of water index range
    integer(kind=4) :: flag_wat_type !> flag: 0=boundary, 1=domain, 2=recharge
    
    !> Mixing and auxiliary variables
    integer(kind=4) :: mix_wat_ind   !> mixing water index (unused)
    integer(kind=4) :: aux_col       !> auxiliary column index (unused)
    integer(kind=4) :: aux_istype    !> auxiliary solid type from previous iteration
    integer(kind=4) :: aux_igzn      !> auxiliary gas zone from previous iteration
    integer(kind=4) :: unit          !> file unit number for I/O
    
    !> Allocatable index arrays
    integer(kind=4), allocatable :: ind_tar_solids(:) !> indices of target solids for range
    integer(kind=4), allocatable :: ind_tar_gases(:)  !> indices of target gases for range
    integer(kind=4), allocatable :: ind_sol_zones(:)  !> solid zone IDs for range
    integer(kind=4), allocatable :: ind_gas_zones(:)  !> gas zone IDs for range
    integer(kind=4), allocatable :: swap(:)           !> species indices to swap (size 2)
    integer(kind=4), allocatable :: aux_swap(:)       !> auxiliary swap array (size 2)
    
    !> Concentration and composition arrays (unused in current implementation)
    real(kind=8), allocatable :: c_nc(:)        !> non-constant concentration (unused)
    real(kind=8), allocatable :: u_init(:,:)    !> initial component concentrations (unused)
    real(kind=8), allocatable :: c1_init(:)     !> initial primary species concentrations (unused)
    real(kind=8), allocatable :: c2_init(:)     !> initial secondary species concentrations (unused)
    real(kind=8), allocatable :: c2_ig(:)       !> initial gas concentrations (unused)
    real(kind=8), allocatable :: gamma_2aq(:)   !> aqueous activity coefficients (unused)
    
    !> String variables for parsing
    character(len=256) :: label      !> label read from file (e.g., "TARGET WATERS")
    character(len=256) :: str        !> general purpose string variable (unused)
    character(len=256) :: int_wat    !> water index or range string from file
    character(len=256) :: int_sol    !> solid index or range string from file
    character(len=256) :: int_gas    !> gas index or range string from file
    character(len=:), allocatable :: int_wat_trim !> trimmed water index string
    character(len=:), allocatable :: int_sol_trim !> trimmed solid index string
    character(len=:), allocatable :: int_gas_trim !> trimmed gas index string
    character(len=:), allocatable :: first_str    !> first index in range (as string)
    character(len=:), allocatable :: last_str     !> last index in range (as string)
    
    !> Logical flags
    logical :: flag_comp      !> flag for component-based speciation (unused)
    logical :: flag_surf      !> flag for surface complexation (unused)
    logical :: flag_aq_phase  !> flag for aqueous phase (unused)
    logical :: flag           !> general purpose flag for domain water detection
    logical :: flag_Se        !> flag indicating species swap needed for stability
    
    !> Default chemistry objects
    type(solid_chemistry_c), target :: solid_chem !> default solid chemistry for waters without solids
    !type(reactive_zone_c), target :: react_zone   !> default reactive zone object
    type(reactive_zone_c), allocatable :: aux_react_zones(:) !> temporary storage for reactive zones (TEMPORARY FIX)
    !type(mineral_zone_c), target :: min_zone      !> default mineral zone object
    
    !> ============================================================================
    !> SECTION 1: INITIALIZATION - Set up basic counts and validate configuration
    !> ============================================================================
    
    nwtype=this%num_wat_types !> get number of predefined water types from chemistry object
    
    !> (COMMENTED BLOCK) Old logic for determining number of solid types
    !> This block handled cases with/without solid reactive zones differently
    ! if (nsrz>0) then
    !     nstype=this%num_mineral_zones !> NOTE: Comment indicates this should use materials instead
    ! else
    !     nstype=0 !> TEMPORARY FIX: no solid reactive zones
    !     !allocate(this%target_solids(1)) !> (commented) allocate single default target solid
    !     !allocate(this%mineral_zones(1)) !> (commented) allocate single default mineral zone
    !     !call this%target_solids(1)%set_reactive_zone(this%reactive_zones(1)) !> (commented) link to first reactive zone
    !     !call this%mineral_zones(1)%set_chem_syst_min_zone(this%chem_syst) !> (commented) set chemical system
    !     !call this%target_solids(1)%set_mineral_zone(this%mineral_zones(1)) !> (commented) link mineral zone
    ! end if
    
    !> Determine number of gas zones based on gas reactive zones count
    if (ngrz>0) then
        ngzns=this%num_gas_zones !> get number of predefined gas zones from chemistry object
    else
        ngzns=0 !> TEMPORARY FIX: no gas zones present in system
    end if
    
    !flag_comp=.true. !> (commented) component-based speciation flag - by default true

    
    !> ============================================================================
    !> SECTION 2: SPECIAL CASE - Systems without solid materials
    !> ============================================================================
    !> For aqueous-only or gas-only systems (no minerals), create a default
    !> reactive zone and mineral zone to maintain code consistency
    
    if (this%num_materials==0) then !> TEMPORARY FIX for systems without solid materials
        !> Step 2.1: Save existing reactive zones to temporary array
        allocate(aux_react_zones(this%num_reactive_zones)) !> allocate temporary storage
        do i=1,this%num_reactive_zones
            call aux_react_zones(i)%copy_react_zone(this%reactive_zones(i)) !> copy existing zones
        end do
        
        !> Step 2.2: Expand reactive zones array to include one additional default zone
        if (allocated(this%reactive_zones)) deallocate(this%reactive_zones) !> deallocate old array
        call this%allocate_reactive_zones(this%num_reactive_zones+1) !> allocate with +1 size
        do i=1,size(aux_react_zones)
            call this%reactive_zones(i)%copy_react_zone(aux_react_zones(i)) !> restore original zones
        end do
        
        !> Step 2.3: Initialize the new default reactive zone (last position)
        call this%reactive_zones(this%num_reactive_zones)%set_chem_syst_react_zone(this%chem_syst) !> link chemical system
        call this%reactive_zones(this%num_reactive_zones)%set_CV_params(this%CV_params) !> set control volume parameters
        call this%reactive_zones(this%num_reactive_zones)%set_speciation_alg_dimensions(.true.) !> set speciation algebra dimensions
        call this%reactive_zones(this%num_reactive_zones)%set_ind_eq_reacts() !> TEMPORARY FIX: set equilibrium reaction indices
        call this%reactive_zones(this%num_reactive_zones)%allocate_ind_var_act_species() !> TEMPORARY FIX: set equilibrium reaction indices
        call this%reactive_zones(this%num_reactive_zones)%set_stoich_mat_react_zone() !> TEMPORARY FIX: set stoichiometric matrix
        call this%reactive_zones(this%num_reactive_zones)%set_ind_gases_stoich_mat() !> TEMPORARY FIX: set gas indices in stoich matrix
        call this%reactive_zones(this%num_reactive_zones)%set_ind_mins_stoich_mat() !> TEMPORARY FIX: set mineral indices in stoich matrix
        allocate(swap(2)) !> TEMPORARY FIX: allocate swap array for species index swapping
        call this%reactive_zones(this%num_reactive_zones)%compute_speciation_alg_arrays(flag_Se,swap) !> compute speciation arrays
        deallocate(aux_react_zones) !> clean up temporary storage
        
        !> Step 2.2b: Re-associate gas zone and target gas reactive_zone pointers
        !> The reallocation above invalidated all pointers into the old reactive_zones array
        if (ngrz > 0) then
            do i=1,this%num_gas_zones
                call this%gas_zones(i)%set_reactive_zone(this%reactive_zones(i))
            end do
            do i=1,this%num_target_gases
                call this%target_gases(i)%set_reactive_zone(this%reactive_zones(this%target_gases(i)%id))
                call this%target_gases_init(i)%set_reactive_zone(this%reactive_zones(this%target_gases_init(i)%id))
            end do
        end if
        
        !> Step 2.4: Create default mineral zone
        call this%allocate_mineral_zones(1) !> TEMPORARY FIX: allocate single mineral zone
        call this%mineral_zones(1)%set_chem_syst_min_zone(this%chem_syst) !> link chemical system to mineral zone
        
        !> Step 2.5: Create default solid chemistry object
        call solid_chem%set_reactive_zone(this%reactive_zones(this%num_reactive_zones)) !> link default reactive zone
        call solid_chem%set_mineral_zone(this%mineral_zones(1)) !> link default mineral zone
        
        !> (COMMENTED BLOCK) Alternative approach: allocate all target solids with default chemistry
        ! call this%allocate_target_solids() !> allocate target solids array
        ! do i=1,this%num_target_solids
        !     this%target_solids(i)=solid_chem !> assign default solid chemistry to all
        ! end do
    end if
    
    !> ============================================================================
    !> SECTION 3: Initialize counters and open input file
    !> ============================================================================
    
    ind_rech=0 !> initialize counter for recharge waters processed
    ind_bd=0   !> initialize counter for boundary waters processed
    ind_dom=0  !> initialize counter for domain (internal) waters processed

    unit=57 !> set file unit number (arbitrary choice, must not conflict with other units)
    open(unit, file=root//'_tar_wat.dat', status='old', action='read') !> open target waters file for reading
    
    !> ============================================================================
    !> SECTION 4: Main file reading loop - Process labels and data
    !> ============================================================================
    do !> infinite loop - exit when "end" label encountered
        read(unit,*) label !> read next label from file
        if (label=='end') then !> end of file marker
            exit !> exit file reading loop
        else if (label=='WATERS') then !> waters section found
            !> Step 4.1: Read counts and allocate arrays
            read(unit,*) num_wat !> read total number of waters from file
            call this%allocate_waters(num_wat) !> allocate array for current target waters
            call this%allocate_waters_init(num_wat) !> allocate array for initial target waters (t=0)
            read(unit,*) num_wat_rech !> read number of recharge (infiltration) waters
            call this%allocate_rech_waters_indices(num_wat_rech) !> allocate array for recharge water indices
            read(unit,*) num_wat_bd !> read number of boundary waters
            call this%allocate_bd_waters_indices(num_wat_bd) !> allocate array for boundary water indices
            !> Validate that the number of target waters declared in the file matches
            !> the number of targets in the mesh provided by the caller.
            if (num_wat-num_wat_rech-num_wat_bd /= num_tar) then
                write(*,*) 'Error: el numero de target waters en ',root//'_tar_wat.dat',&
                    ' (',num_wat-num_wat_rech-num_wat_bd,') no coincide con el numero de targets de la malla (',num_tar,').'
                error stop "El numero de target waters no coincide con el numero de targets de la malla"
            end if
            call this%allocate_tar_wat_indices(num_wat-num_wat_rech-num_wat_bd) !> allocate array for target water indices
            !call this%allocate_ext_waters_indices(num_wat_bd+num_wat_rech) !> allocate external waters (boundary + recharge)
            !call this%allocate_waters(this%num_target_waters+this%num_rech_waters+this%num_bd_waters) !> allocate domain waters (total - external)
            !call this%allocate_waters_init(this%num_target_waters_init+this%num_rech_waters+this%num_bd_waters) !> allocate initial domain waters
            !> Step 4.2: Special allocation for systems without materials
            if (this%num_materials==0) then !> no solid materials in system (aqueous/gas only)
                call this%allocate_target_solids(this%num_target_waters) !> TEMPORARY FIX: assume 1:1 with domain waters
                do i=1,this%num_target_solids !> loop over all allocated target solids
                    call this%target_solids(i)%copy_solid_chemistry(solid_chem) !> assign default solid chemistry to each
                end do  
            end if
            
            !> Step 4.3: Special allocation for systems without gas zones
            if (ngzns==0) then !> no gas zones defined in system
                call this%allocate_target_gases(this%num_target_waters) !> TEMPORARY FIX: assume 1:1 with domain waters
            end if
            
            !> Step 4.4: Initialize auxiliary indices for optimization
            aux_istype=0 !> initialize auxiliary solid type index (0 = first iteration, no previous zone)
            aux_igzn=0 !> initialize auxiliary gas zone index (0 = first iteration, no previous zone)
            
            !if (num_wat_rech>0 .or. num_tar_wat_bd>0) then !> (commented) old check for external waters
                !> Step 4.5: Main loop - read and process each target water line
                do !> loop until all target waters processed (exits when last_wat or tar_wat_ind == num_tar_wat)
                    !> Read one data line from file: water_id water_type solid_id gas_id flag
                    read(unit,*) int_wat, wtype, int_sol, int_gas, flag_wat_type
                    !print *, int_wat, wtype, public, int_sol, int_gas, flag_wat_type !> (commented) debug output
                    
                    flag=.true. !> initialize flag to true (will be set false for domain waters later)
                    
                    !> Trim whitespace from input strings for proper parsing
                    int_wat_trim=trim(int_wat) !> remove trailing spaces from water index string
                    int_sol_trim=trim(int_sol) !> remove trailing spaces from solid index string  
                    int_gas_trim=trim(int_gas) !> remove trailing spaces from gas index string
                    
                    !> Find hyphen position to detect range notation (e.g., "5-10" means indices 5 through 10)
                    ind_bar_wat=index(int_wat_trim,'-') !> find '-' in water string (0 if single index)
                    ind_bar_sol=index(int_sol_trim,'-') !> find '-' in solid string (0 if single index)
                    ind_bar_gas=index(int_gas_trim,'-') !> find '-' in gas string (0 if single index)
                    
                    !> Validate water type index
                    if (wtype<1 .or. wtype>nwtype) then !> water type must be in valid range
                        error stop "Water type index out of bounds" !> abort if invalid
                    !> (COMMENTED BLOCK) Additional validation options - not currently enforced
                    ! else if (int_sol<0 .or. int_sol>this%num_materials) then
                    !     error stop "Solid zone index out of bounds"
                    ! else if (int_gas<0 .or. int_gas>this%num_gas_zones) then
                    !     error stop "Gas zone index out of bounds"
                    !> (COMMENTED BLOCK) Old water type flag handling - replaced by logic below
                    !else if (flag_wat_type==0) then !> boundary water
                    !    !ind_bd=ind_bd+1 !> counter boundary waters
                    !    !this%bd_waters_indices(ind_bd)=tar_wat_ind
                    !else if (flag_wat_type==2) then !> external water
                    !    !ind_ext=ind_ext+1 !> counter external waters
                    !    !this%ext_waters_indices(ind_ext)=tar_wat_ind
                    !else if (flag_wat_type==1) then !> domain water
                    !    !ind_dom=ind_dom+1 !> counter domain waters
                        !this%tar_wat_indices(ind_dom)=tar_wat_ind
                        !flag=.false.
                    else if (flag_wat_type<0 .or. flag_wat_type>2) then !> validate flag value
                        error stop "Water type flag not implemented" !> flag must be 0, 1, or 2
                    end if
                    
                    !> ================================================================
                    !> BRANCH A: Range of target waters (e.g., "5-10")
                    !> ================================================================
                    if (ind_bar_wat>0) then !> hyphen found in water string = range notation
                        !> Parse range: extract first and last indices
                        first_str=int_wat_trim(1:ind_bar_wat-1) !> substring before '-' (first index)
                        last_str=int_wat_trim(ind_bar_wat+1:) !> substring after '-' (last index)
                        read(first_str,*) first_wat !> convert first index string to integer
                        read(last_str,*) last_wat !> convert last index string to integer
                        int_wat_size=last_wat-first_wat+1 !> calculate size of range (inclusive)
                        
                        !> Validate range bounds
                        if (first_wat<1 .or. first_wat>last_wat .or. last_wat>this%num_waters) then
                            error stop "water index out of bounds" !> range must be valid and within total
                        !> (COMMENTED BLOCK) Additional validation options
                        !else if (tar_wat_ind<1 .or. tar_wat_ind>this%num_tar_waters) then
                        !    error stop "Target water index out of bounds"
                        !else if (wtype<1 .or. wtype>nwtype) then
                        !    error stop "Water type index out of bounds"
                        !else if (ind_tar_sol<0 .or. ind_tar_sol>this%num_target_solids) then
                        !    error stop "Target solid index out of bounds"
                        !else if (int_gas<0 .or. int_gas>ngzns) then
                        !    error stop "Gas type index out of bounds"
                        
                        !> Assign water indices to appropriate category based on flag
                        else if (flag_wat_type==0) then !> BOUNDARY WATER
                            do i=1,int_wat_size !> loop over all waters in range
                                this%bd_waters_indices(ind_bd+i)=first_wat+i-1 !> store boundary water index
                            end do
                            !this%bd_waters_indices_init=this%bd_waters_indices !> copy current to initial array
                            ind_bd=ind_bd+int_wat_size !> increment boundary water counter by range size
                        else if (flag_wat_type==2) then !> RECHARGE WATER (infiltration/precipitation)
                            !this%ext_waters_indices(ind_ext)=tar_wat_ind !> (commented) old single index approach
                            do i=1,int_wat_size !> loop over all waters in range
                                this%rech_waters_indices(ind_rech+i)=first_wat+i-1 !> store recharge water index
                            end do
                            ind_rech=ind_rech+int_wat_size !> increment recharge water counter by range size
                        else if (flag_wat_type==1) then !> TARGET WATER (internal)
                            !this%tar_wat_indices(ind_dom)=tar_wat_ind !> (commented) old single index approach
                            do i=1,int_wat_size !> loop over all waters in range
                                this%tar_wat_indices(ind_dom+i)=first_wat+i-1 !> store target water index
                            end do
                            !this%tar_wat_indices_init=this%tar_wat_indices !> copy current to initial array
                            ind_dom=ind_dom+int_wat_size !> increment domain water counter by range size
                            flag=.false. !> set flag to false for domain waters
                        !else
                        !    error stop "Water type flag out of bounds"
                        end if
                        
                        !> ============================================================
                        !> Parse target solids (mineral phase) - range or single index
                        !> ============================================================
                        if (ind_bar_sol>0) then !> hyphen found in solid string = range notation
                            !> Parse solid range: extract first and last indices
                            first_str=int_sol_trim(1:ind_bar_sol-1) !> substring before '-' (first solid)
                            last_str=int_sol_trim(ind_bar_sol+1:) !> substring after '-' (last solid)
                            read(first_str,*) first_sol !> convert first solid index to integer
                            read(last_str,*) last_sol !> convert last solid index to integer
                            int_sol_size=last_sol-first_sol+1 !> calculate solid range size (inclusive)
                            
                            !> Validate solid range bounds
                            if (first_sol<0 .or. first_sol>last_sol .or. last_sol>this%num_target_solids) then
                                error stop "Target solid index out of bounds" !> range must be valid
                            end if
                            
                            !> Verify consistency: solid range size must match water range size
                            if (int_sol_size/=int_wat_size) then
                                error stop "Dimension error: number of target solids not consistent with target waters interval"
                            end if
                            
                            !> Build arrays of solid indices and corresponding zone IDs
                            allocate(ind_tar_solids(int_sol_size), ind_sol_zones(int_sol_size)) !> allocate for range
                            ind_tar_solids = [(i, i=first_sol, last_sol)] !> create sequential array [first_sol, first_sol+1, ..., last_sol]
                            do i=1,int_sol_size !> loop over solid range
                                ind_sol_zones(i)=this%target_solids(ind_tar_solids(i))%id !> get solid zone ID for each target solid
                            end do
                            ! print *, "ind_tar_solids", ind_tar_solids !> (commented) debug output
                            ! print *, "ind_sol_zones", ind_sol_zones !> (commented) debug output
                        else !> single target solid index (no range)
                            read(int_sol_trim,*) tar_sol_ind !> convert solid index string to integer
                            
                            !> Validate single solid index
                            if (tar_sol_ind<0 .or. tar_sol_ind>this%num_target_solids) then
                                error stop "Target solid index out of bounds" !> solid index must be valid
                            else
                                !> All waters in range point to same single solid
                                allocate(ind_tar_solids(int_wat_size),ind_sol_zones(int_wat_size)) !> allocate arrays sized for water range
                                ind_tar_solids=tar_sol_ind !> set all elements to same solid index
                                if (tar_sol_ind>0) then !> solid index is valid (>0 means solid present)
                                    ind_sol_zones=this%target_solids(tar_sol_ind)%id !> get solid zone ID for this solid
                                else !> tar_sol_ind=0 means no solid phase
                                    ind_sol_zones=0 !> set zone ID to 0 (no solid chemistry)
                                end if
                            end if
                        end if
                        
                        !> ============================================================
                        !> Parse target gases (gas phase) - range or single index
                        !> ============================================================
                        if (ind_bar_gas>0) then !> hyphen found in gas string = range notation
                            !> Parse gas range: extract first and last indices
                            first_str=int_gas_trim(1:ind_bar_gas-1) !> substring before '-' (first gas)
                            last_str=int_gas_trim(ind_bar_gas+1:) !> substring after '-' (last gas)
                            read(first_str,*) first_gas !> convert first gas index to integer
                            read(last_str,*) last_gas !> convert last gas index to integer
                            int_gas_size=last_gas-first_gas+1 !> calculate gas range size (inclusive)
                            
                            !> Validate gas range bounds
                            if (first_gas<0 .or. first_gas>last_gas .or. last_gas>this%num_target_gases) then
                                error stop "Target gas index out of bounds" !> range must be valid
                            end if
                            
                            !> Verify consistency: gas range size must match water range size
                            if (int_gas_size/=int_wat_size) then
                                error stop "Dimension error: number of target gases not consistent with target waters interval"
                            end if
                            
                            !> Build arrays of gas indices and corresponding zone IDs
                            allocate(ind_tar_gases(int_gas_size), ind_gas_zones(int_gas_size)) !> allocate for range
                            ind_tar_gases = [(i, i=first_gas, last_gas)] !> create sequential array [first_gas, first_gas+1, ..., last_gas]
                            do i=1,int_gas_size !> loop over gas range
                                ind_gas_zones(i)=this%target_gases(ind_tar_gases(i))%id !> get gas zone ID for each target gas
                            end do
                            ! print *, "ind_tar_gases", ind_tar_gases !> (commented) debug output
                            ! print *, "ind_gas_zones", ind_gas_zones !> (commented) debug output
                        else !> single target gas index (no range)
                            read(int_gas_trim,*) tar_gas_ind !> convert gas index string to integer
                            
                            !> Validate single gas index
                            if (tar_gas_ind<0 .or. tar_gas_ind>this%num_target_gases) then
                                error stop "Target gas index out of bounds" !> gas index must be valid
                            else
                                !> All waters in range point to same single gas
                                allocate(ind_tar_gases(int_wat_size),ind_gas_zones(int_wat_size)) !> allocate arrays sized for water range
                                ind_tar_gases=tar_gas_ind !> set all elements to same gas index
                                if (tar_gas_ind>0) then !> gas index is valid (>0 means gas present)
                                    ind_gas_zones=this%target_gases(tar_gas_ind)%id !> get gas zone ID for this gas
                                else !> tar_gas_ind=0 means no gas phase
                                    ind_gas_zones=0 !> set zone ID to 0 (no gas chemistry)
                                end if
                            end if
                        end if
                        !> (COMMENTED BLOCK) Alternative approach using get_tar_sol_ind method - not currently used
                        ! if (int_sol>0) then
                        !     ind_tar_sol=this%get_tar_sol_ind(int_sol) !> retrieve target solids for given solid zone
                        !     if (size(ind_tar_sol)/=int_size) then
                        !         error stop "Dimension error: number of target solids associated to solid zone not consistent with &
                        !             int_wat size"
                        !     end if
                        ! else
                        !     allocate(ind_tar_sol(int_size))
                        !     ind_tar_sol=0 !> set to zero if no solids
                        ! end if
                        ! if (int_gas>0) then
                        !     ind_tar_gases=this%get_tar_gas_ind(int_gas) !> retrieve target gases for given gas zone
                        !     if (size(ind_tar_gases)/=int_size) then
                        !         error stop "Dimension error: number of target gases associated to gas zone not consistent with & 
                        !             int_wat size"
                        !     end if
                        ! else
                        !     allocate(ind_tar_gases(int_size))
                        !     ind_tar_gases=0 !> set to zero if no gases
                        ! end if
                        
                        !> ============================================================
                        !> Initialize each target water in range
                        !> ============================================================
                        do i=1,int_wat_size !> loop over all waters in range [first_wat, last_wat]
                            !> Call read_tar_wat_line to initialize water chemistry for this water
                            !> @param flag: flag external water (false=domain, true=boundary/recharge)
                            !> @param ind_sol_zones(i): solid zone ID for water index (first_wat+i-1)
                            !> @param ind_gas_zones(i): gas zone ID for water index (first_wat+i-1)
                            !> @param first_wat+i-1: absolute target water index
                            !> @param wtype: water type index (1..nwtype)
                            !> @param ind_tar_solids(i): target solid index for this water
                            !> @param ind_tar_gases(i): target gas index for this water
                            !> @param aux_istype: auxiliary solid type index (output)
                            !> @param aux_igzn: auxiliary gas zone index (output)
                            call this%read_tar_wat_line(flag,ind_sol_zones(i),ind_gas_zones(i),first_wat+i-1,wtype,&
                                ind_tar_solids(i),ind_tar_gases(i),aux_istype,aux_igzn)                        
                        end do
                        ! print *, this%waters(2)%solid_chemistry%rk_new !> (commented) debug output
                        ! print *, this%target_solids(1)%rk_new !> (commented) debug output
                        
                        !> Check if we've processed all target waters
                        if (last_wat==this%num_waters) then !> last water in range matches total count
                            exit !> exit main reading loop - all waters initialized
                        end if
                        
                        !> Clean up temporary arrays for this range iteration
                        deallocate(ind_tar_solids,ind_tar_gases,ind_sol_zones,ind_gas_zones) !> free memory before next iteration
                        
                    !> ================================================================
                    !> BRANCH B: Single target water (no range notation)
                    !> ================================================================
                    else !> no hyphen found in water string = single water index
                        !> Parse single indices for water, solid, and gas
                        read(int_wat_trim,*) tar_wat_ind !> convert water index string to integer
                        read(int_sol_trim,*) tar_sol_ind !> convert solid index string to integer
                        read(int_gas_trim,*) tar_gas_ind !> convert gas index string to integer
                        
                        !> Validate indices and dimension consistency
                        if (tar_wat_ind<1 .or. tar_wat_ind>this%num_waters) then
                            error stop "Target water index out of bounds" !> water index must be valid
                        else if (ind_bar_sol>0) then !> solid has range notation but water doesn't
                            error stop "Dimension error: a single target water cannot point to a range of target solids"
                        else if (ind_bar_gas>0) then !> gas has range notation but water doesn't
                            error stop "Dimension error: a single target water cannot point to a range of target gases"
                        else if (tar_sol_ind<0 .or. tar_sol_ind>this%num_target_solids) then
                            error stop "Target solid index out of bounds" !> solid index must be valid
                        else if (tar_gas_ind<0 .or. tar_gas_ind>this%num_target_gases) then
                            error stop "Target gas index out of bounds" !> gas index must be valid
                        !> Assign single water to appropriate category based on flag
                        else if (flag_wat_type==0) then !> BOUNDARY WATER (Dirichlet BC)
                            ind_bd=ind_bd+1 !> increment boundary water counter
                            this%bd_waters_indices(ind_bd)=tar_wat_ind !> store boundary water index
                            !this%bd_waters_indices_init(ind_bd)=this%bd_waters_indices(ind_bd) !> copy to initial array
                        else if (flag_wat_type==2) then !> RECHARGE WATER (infiltration/precipitation)
                            ind_rech=ind_rech+1 !> increment recharge water counter
                            this%rech_waters_indices(ind_rech)=tar_wat_ind !> store recharge water index
                        else if (flag_wat_type==1) then !> DOMAIN WATER (internal/mobile)
                            ind_dom=ind_dom+1 !> increment domain water counter
                            this%tar_wat_indices(ind_dom)=tar_wat_ind !> store domain water index
                            !this%tar_wat_indices_init(ind_dom)=tar_wat_ind !> store domain water index
                            flag=.false. !> set flag to false for domain waters
                            !> (COMMENTED) Old approach using loop_read_tar_wat_init - replaced by read_tar_wat_line
                            !call this%loop_read_tar_wat_init(flag,this%wat_types,init_sol_types,init_gas_types,nsrz,ngrz,tar_wat_ind,wtype, public,&
                            !    istype, public,int_gas,aux_istype, public,aux_int_gas,solid_chem)
                        end if
                        
                        !> Extract solid zone index from target solid
                        if (tar_sol_ind>0) then !> solid index is valid (>0 means solid present)
                            iszn=this%target_solids(tar_sol_ind)%id !> get solid zone ID from target solid
                            !> (COMMENTED) Old validation approach
                            ! if (size(ind_sol_zn)/=1) then
                            !     error stop "Dimension error: number of target solids associated to solid zone must be 1"
                            ! end if
                        else !> tar_sol_ind=0 means no solid phase
                            !allocate(ind_tar_solids(1)) !> (commented) old allocation approach
                            iszn=0 !> set zone ID to 0 (no solid chemistry)
                        end if
                        
                        !> Extract gas zone index from target gas
                        if (tar_gas_ind>0) then !> gas index is valid (>0 means gas present)
                            igzn=this%target_gases(tar_gas_ind)%id !> get gas zone ID from target gas
                            !> (COMMENTED) Old validation approach
                            ! if (size(ind_tar_gases)/=1) then
                            !     error stop "Dimension error: number of target gases associated to gas zone must be 1"
                            ! end if
                        else !> tar_gas_ind=0 means no gas phase
                            !allocate(ind_tar_gases(1)) !> (commented) old allocation approach
                            igzn=0 !> set zone ID to 0 (no gas chemistry)
                        end if
                        
                        !> Initialize single target water chemistry
                        !> Call read_tar_wat_line with single indices
                        !> @param flag: flag external water (false=domain, true=boundary/recharge)
                        !> @param iszn: solid zone ID for this water
                        !> @param igzn: gas zone ID for this water
                        !> @param tar_wat_ind: absolute target water index
                        !> @param wtype: water type index (1..nwtype)
                        !> @param tar_sol_ind: target solid index
                        !> @param tar_gas_ind: target gas index
                        !> @param aux_istype: auxiliary solid type index (output)
                        !> @param aux_igzn: auxiliary gas zone index (output)
                        call this%read_tar_wat_line(flag,iszn,igzn,tar_wat_ind,wtype,&
                                tar_sol_ind,tar_gas_ind,aux_istype,aux_igzn)
                            !> (COMMENTED) Old approach setting initial solid types
                            !this%waters_init(ind_tar_solids)=init_sol_types(sol_zone) !> set initial solid type
                            !call this%target_solids_init(ind_dom)%set_target(ind_dom) !> set target index
                        !print *, this%waters(2)%solid_chemistry%rk_new !> (commented) debug output
                        !print *, this%target_solids(1)%rk_new !> (commented) debug output
                        
                        !> Check if we've processed all target waters
                        if (tar_wat_ind==this%num_waters) then !> single water index matches total count
                            exit !> exit main reading loop - all waters initialized
                        end if
                        !deallocate(ind_tar_solids,ind_tar_gases) !> (commented) no deallocation needed for single water
                    end if !> end of single water branch
                    
                    !> (COMMENTED) Old reaction rate computation approach - not currently used
                    !allocate(rk(this%waters(tar_wat_ind)%indices_rk%num_cols)) !> chapuza (temporary fix)
                    !call this%waters(tar_wat_ind)%compute_rk(rk) !> compute reaction rates
                    
                    !> ================================================================
                    !> (COMMENTED BLOCK) Old target water initialization approach
                    !> ================================================================
                    !> This large commented section shows an older implementation strategy
                    !> that set target waters from water types and reactive zones directly
                    !> rather than using read_tar_wat_line. Preserved for reference.
                !    !aux_istype=istype !> store solid type for optimization
                !    this%waters(tar_wat_ind)=water_types(wtype) !> assign water from type library
                !    ! if (counter_swap==0) then
                !    !     call this%waters(tar_wat_ind)%set_aq_phase(this%chem_syst%aq_phase)
                !    !     !call this%waters(tar_wat_ind)%set_indices_aq_species_aq_chem()
                !    ! end if
                !    if (istype>0) then !> solid chemistry is present
                !        this%target_solids(tar_wat_ind)=init_sol_types(istype) !> assign solid from type library
                !        !> chapuza intercambio (temporary fix for cation exchange)
                !        if (this%reactive_zones(ngzns+nstype*int_gas+1)%cat_exch_zone%num_surf_compl>0) then
                !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+1))
                !        else
                !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+istype))
                !        end if
                !        if (int_gas>0) then !> both solid and gas present
                !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+istype))
                !            this%target_gases(tar_wat_ind)=init_gas_types(int_gas) !> assign gas from type library
                !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+istype))
                !            call this%waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
                !        else !> solid only, no gas
                !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+istype))
                !        end if
                !        call this%waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(tar_wat_ind))
                !    else !> no solid chemistry
                !        if (int_gas>0) then !> gas only, no solid
                !            this%target_gases(tar_wat_ind)=init_gas_types(int_gas)
                !            ! call this%reactive_zones(int_gas)%set_ind_eq_reacts() !> chapuza
                !            ! call this%reactive_zones(int_gas)%set_stoich_mat_react_zone() !> chapuza
                !            ! call this%reactive_zones(int_gas)%set_ind_gases_stoich_mat() !> chapuza
                !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(int_gas))
                !            call this%waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
                !            call solid_chem%set_reactive_zone(this%reactive_zones(int_gas))
                !        else !> neither solid nor gas - aqueous only
                !            ! call react_zone%set_ind_eq_reacts() !> chapuza
                !            ! call react_zone%set_stoich_mat_react_zone() !> chapuza
                !            ! call react_zone%set_ind_gases_stoich_mat() !> chapuza
                !            ! call react_zone%set_ind_mins_stoich_mat() !> chapuza
                !            call solid_chem%set_reactive_zone(react_zone)
                !        end if
                !        call solid_chem%set_mineral_zone(min_zone)
                !        call this%waters(tar_wat_ind)%set_solid_chemistry(solid_chem)
                !    end if
                !    !> Detect and handle new reactive zones for optimization
                !    if (aux_istype==0 .or. aux_istype/=istype .or. aux_int_gas/=int_gas) then !> assume target waters grouped by reactive zones
                !        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_speciation_alg_dimensions(flag_comp)
                !        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_ind_eq_reacts() !> chapuza
                !        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_stoich_mat_react_zone() !> chapuza
                !        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_ind_gases_stoich_mat() !> chapuza
                !        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_ind_mins_stoich_mat() !> chapuza
                !        call this%waters(tar_wat_ind)%set_ind_species() !> set species indices
                !        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
                !            flag_Se,swap) !> compute speciation algorithm arrays
                !        if (flag_Se.eqv..true.) then !> swap species indices if needed
                !            aux_swap(1)=this%waters(tar_wat_ind)%ind_var_act_species(swap(1))
                !            aux_swap(2)=this%waters(tar_wat_ind)%ind_var_act_species(swap(2))
                !            !     this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
                !            !print *, this%waters(tar_wat_ind)%ind_var_act_species
                !            !print *, this%waters(tar_wat_ind)%ind_sec_species
                !            this%waters(tar_wat_ind)%ind_var_act_species(swap(1))=aux_swap(2)
                !            this%waters(tar_wat_ind)%ind_var_act_species(swap(2))=aux_swap(1)
                !            ! this%waters(tar_wat_ind)%ind_prim_species(swap(1))=aux_swap(2) !> index of primary species
                !            ! this%waters(tar_wat_ind)%ind_sec_species(swap(2)-&
                !            !     this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=&
                !            !     aux_swap(1) !> index of secondary species
                !        end if
                !    else if (aux_istype>0 .or. aux_int_gas>0) then !> reactive zone unchanged - reuse previous indices
                !        this%waters(tar_wat_ind)%ind_var_act_species=this%waters(tar_wat_ind-1)%ind_var_act_species
                !        !this%waters(tar_wat_ind)%ind_sec_species=this%waters(tar_wat_ind-1)%ind_sec_species
                !    end if
                !    print *, this%waters(tar_wat_ind)%ind_var_act_species
                !!> Chapuza (temporary fix for mineral kinetics)
                !    !if (associated(this%waters(tar_wat_ind)%solid_chemistry%mineral_zone)) then
                !        if (this%waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin<&
                !            this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin) then
                !            call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod(&
                !                this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                !                this%waters(tar_wat_ind)%solid_chemistry%mineral_zone%ind_min_chem_syst(1:&
                !                this%waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin))
                !        else
                !            call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod()
                !        end if
                !    !end if
                !    call this%waters(tar_wat_ind)%allocate_reaction_rates() !> allocate reaction rate arrays
                !    call this%waters(tar_wat_ind)%set_indices_rk() !> set reaction rate indices
                !    aux_istype=istype !> store current solid type for next iteration
                !    aux_int_gas=int_gas !> store current gas zone for next iteration
                
                end do !> end of main reading loop
            
            !> ================================================================
            !> (COMMENTED BLOCK) Alternative logic for systems without recharge/boundary waters
            !> ================================================================
            !> This commented section shows an alternative implementation for cases
            !> where all waters are domain waters (no external forcing).
            !else !> (COMMENTED) branch for systems without recharge or boundary waters
            !    do! i=1,this%num_tar_waters !> (COMMENTED) loop over all target waters
            !        read(unit,*) int_wat, wtype, public, int_sol, int_gas !> (COMMENTED) read water data
            !        int_wat_trim=trim(int_wat)
            !        ind_bar=index(int_wat_trim,'-')
            !        if (wtype<1 .or. wtype>nwtype) then
            !            error stop "Water type index out of bounds"
            !        else if (ind_tar_solids<0 .or. ind_tar_solids>this%num_target_solids) then
            !            error stop "Target solid index out of bounds"
            !        else if (ind_tar_gases<0 .or. ind_tar_gases>this%num_target_gases) then
            !            error stop "Target gas index out of bounds"
            !        end if
            !        if (ind_bar>0) then !> (COMMENTED) range notation detected
            !            first_str=int_wat_trim(1:ind_bar-1) !> first target of int_wat
            !            last_str=int_wat_trim(ind_bar+1:) !> last target of int_wat
            !            read(first_str,*) first
            !            read(last_str,*) last
            !            if (first<1 .or. first>last .or. last>num_tar_wat) then
            !                error stop "Target water index out of bounds"
            !            !else if (tar_wat_ind<1 .or. tar_wat_ind>this%num_tar_waters) then
            !            !    error stop "Target water index out of bounds"
            !            !else if (wtype<1 .or. wtype>nwtype) then
            !            !    error stop "Water type index out of bounds"
            !            !else if (ind_tar_solids<0 .or. ind_tar_solids>this%num_target_solids) then
            !            !    error stop "Target solid index out of bounds"
            !            !else if (int_gas<0 .or. int_gas>ngzns) then
            !            !    error stop "Gas type index out of bounds"
            !            !else if (flag_wat_type==0) then !> boundary water
            !            !    ind_bd=ind_bd+1 !> counter boundary waters
            !            !    this%bd_waters_indices(ind_bd)=tar_wat_ind
            !            !else if (flag_wat_type==2) then !> external water
            !            !    ind_ext=ind_ext+1 !> counter external waters
            !            !    this%ext_waters_indices(ind_ext)=tar_wat_ind
            !            !else if (flag_wat_type==1) then !> domain water
            !            !    ind_dom=ind_dom+1 !> counter domain waters
            !            !    this%tar_wat_indices(ind_dom)=tar_wat_ind
            !            !else
            !            !    error stop "Water type flag out of bounds"
            !            end if
            !            do i=first,last !> (COMMENTED) loop over range
            !                this%tar_wat_indices(i)=i !> all waters are domain waters
            !                call this%read_tar_wat_line(.false.,nsrz,ngrz,i,wtype, public,&
            !                    int_sol,int_gas,aux_istype, public,aux_int_gas) !> initialize each water
            !            end do
            !            if (last==num_tar_wat) then !> all waters processed
            !                exit
            !            end if
            !        else !> (COMMENTED) single target water
            !            read(int_wat_trim,*) tar_wat_ind !> target water index
            !            if (tar_wat_ind<1 .or. tar_wat_ind>num_tar_wat) then
            !                error stop "Target water index out of bounds"
            !            else
            !                this%tar_wat_indices(tar_wat_ind)=tar_wat_ind !> domain water
            !                call this%read_tar_wat_line(.false.,nsrz,ngrz,tar_wat_ind,wtype, public,&
            !                    int_sol,int_gas,aux_istype, public,aux_int_gas) !> initialize water
            !                !call this%loop_read_tar_wat_init(flag,this%wat_types,init_sol_types,init_gas_types,nsrz,ngrz,tar_wat_ind,wtype, public,&
            !                !    istype, public,int_gas,aux_istype, public,aux_int_gas,solid_chem)
            !                !this%waters_init(ind_tar_solids)=init_sol_types(sol_zone) !> we set the initial solid type
            !                !call this%target_solids_init(ind_dom)%set_target(ind_dom) !> we set the target index
            !            end if
            !            if (tar_wat_ind==num_tar_wat) then !> all waters processed
            !                exit
            !            end if
            !        end if
            !    end do
            !    !do i=1,this%num_tar_waters
            !    !    read(unit,*) tar_wat_ind, wtype, public, istype, public, int_gas
            !    !    if (tar_wat_ind<1 .or. tar_wat_ind>this%num_tar_waters) then
            !    !        error stop "Target water index out of bounds"
            !    !    else if (wtype<1 .or. wtype>nwtype) then
            !    !        error stop "Water type index out of bounds"
            !    !    else if (istype<0 .or. istype>nstype) then
            !    !        error stop "Solid type index out of bounds"
            !    !    else if (int_gas<0 .or. int_gas>ngzns) then
            !    !        error stop "Gas type index out of bounds"
            !    !    else
            !    !        ind_dom=ind_dom+1
            !    !        this%tar_wat_indices(ind_dom)=tar_wat_ind
            !    !    end if
            !    !    call this%loop_read_tar_wat_init(.false.,this%wat_types,init_sol_types,init_gas_types,nsrz,ngrz,tar_wat_ind,&
            !    !        wtype, public,istype, public,int_gas,aux_istype, public,aux_int_gas,solid_chem)
            !    !!    !aux_istype=istype
            !    !!    this%waters(tar_wat_ind)=water_types(wtype)
            !    !!    ! if (counter_swap==0) then
            !    !!    !     call this%waters(tar_wat_ind)%set_aq_phase(this%chem_syst%aq_phase)
            !    !!    !     !call this%waters(tar_wat_ind)%set_indices_aq_species_aq_chem()
            !    !!    ! end if
            !    !!    if (istype>0) then
            !    !!        this%target_solids(tar_wat_ind)=init_sol_types(istype)
            !    !!        !> chapuza intercambio
            !    !!        if (this%reactive_zones(ngzns+nstype*int_gas+1)%cat_exch_zone%num_surf_compl>0) then
            !    !!            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+1))
            !    !!        else
            !    !!            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+istype))
            !    !!        end if
            !    !!        if (int_gas>0) then
            !    !!            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+istype))
            !    !!            this%target_gases(tar_wat_ind)=init_gas_types(int_gas)
            !    !!            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+istype))
            !    !!            call this%waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
            !    !!        else
            !    !!            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+istype))
            !    !!        end if
            !    !!        call this%waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(tar_wat_ind))
            !    !!    else if (int_gas>0) then
            !    !!        this%target_gases(tar_wat_ind)=init_gas_types(int_gas)
            !    !!        call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(int_gas))
            !    !!        call this%waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
            !    !!        call solid_chem%set_reactive_zone(this%reactive_zones(int_gas))
            !    !!        call this%waters(tar_wat_ind)%set_solid_chemistry(solid_chem)
            !    !!    else
            !    !!        call solid_chem%set_reactive_zone(react_zone)
            !    !!        call this%waters(tar_wat_ind)%set_solid_chemistry(solid_chem)
            !    !!    end if
            !    !!    !> we check if there is a new reactive zone
            !    !!    if (aux_istype==0 .or. aux_istype/=istype .or. aux_int_gas/=int_gas) then !> we assume target waters are grouped by their reactive zones
            !    !!        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_speciation_alg_dimensions(flag_comp)
            !    !!        call this%waters(tar_wat_ind)%set_ind_species()
            !    !!        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
            !    !!            flag_Se,swap)
            !    !!        if (flag_Se.eqv..true.) then !> we swap indices of species
            !    !!            aux_swap(1)=this%waters(tar_wat_ind)%ind_var_act_species(swap(1))
            !    !!            aux_swap(2)=this%waters(tar_wat_ind)%ind_var_act_species(swap(2))
            !    !!            !     this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
            !    !!            !print *, this%waters(tar_wat_ind)%ind_var_act_species
            !    !!            !print *, this%waters(tar_wat_ind)%ind_sec_species
            !    !!            this%waters(tar_wat_ind)%ind_var_act_species(swap(1))=aux_swap(2)
            !    !!            this%waters(tar_wat_ind)%ind_var_act_species(swap(2))=aux_swap(1)
            !    !!            ! this%waters(tar_wat_ind)%ind_prim_species(swap(1))=aux_swap(2) !> index of primary species
            !    !!            ! this%waters(tar_wat_ind)%ind_sec_species(swap(2)-&
            !    !!            !     this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=&
            !    !!            !     aux_swap(1) !> index of secondary species
            !    !!        end if
            !    !!    else if (aux_istype>0 .or. aux_int_gas>0) then !> indices remain the same because reactive zone is the same
            !    !!        this%waters(tar_wat_ind)%ind_var_act_species=this%waters(tar_wat_ind-1)%ind_var_act_species
            !    !!        !this%waters(tar_wat_ind)%ind_sec_species=this%waters(tar_wat_ind-1)%ind_sec_species
            !    !!    end if
            !    !!    ! if (flag_Se.eqv..true.) then
            !    !!    !     counter_swap=counter_swap+1
            !    !!    !     aux_swap=swap
            !    !!    !     !print *, this%waters(tar_wat_ind)%indices_aq_species
            !    !!    !     do j=tar_wat_ind+1,this%num_tar_waters
            !    !!    !         call this%waters(j)%set_aq_phase(this%chem_syst%aq_phase)
            !    !!    !         call this%waters(j)%solid_chemistry%reactive_zone%set_speciation_alg_dimensions(&
            !    !!    !             flag_comp)
            !    !!    !         call this%waters(j)%set_ind_prim_sec_species()
            !    !!    !         this%waters(j)%ind_prim_species(swap(1))=swap(2)
            !    !!    !         this%waters(j)%ind_sec_species(swap(2)-&
            !    !!    !             this%waters(j)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=swap(1) !> chapuza
            !    !!    !         !call this%waters(j)%set_indices_aq_species_aq_chem()
            !    !!    !         ! this%waters(j)%indices_aq_species(swap(2))=swap(1)
            !    !!    !         ! this%waters(j)%indices_aq_species(swap(1))=swap(2)
            !    !!    !     end do
            !    !!    ! end if
            !    !!    !print *, this%waters(tar_wat_ind)%ind_var_act_species
            !    !!    !print *, this%waters(tar_wat_ind)%ind_sec_species
            !    !!    ! if (aux_swap(1)>0 .AND. aux_swap(2)>0) then
            !    !!    !     !call this%waters(tar_wat_ind)%set_indices_aq_species_aq_chem()
            !    !!    !     this%waters(tar_wat_ind)%ind_prim_species(aux_swap(1))=aux_swap(2)
            !    !!    !     this%waters(tar_wat_ind)%ind_sec_species(swap(2)-&
            !    !!    !         this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=swap(1) !> chapuza
            !    !!    ! end if
            !    !!    if (this%waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin<&
            !    !!        this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin) then
            !    !!        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod(&
            !    !!            this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
            !    !!            this%waters(tar_wat_ind)%solid_chemistry%mineral_zone%ind_min_chem_syst(1:&
            !    !!            this%waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin))
            !    !!    else
            !    !!        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod()
            !    !!    end if
            !    !!    call this%waters(tar_wat_ind)%allocate_reaction_rates()
            !    !!    call this%waters(tar_wat_ind)%set_indices_rk()
            !    !!    aux_istype=istype
            !    !!    aux_int_gas=int_gas
            !    !end do
            !end if
        
        !> ================================================================
        !> Section 5: Non-TARGET WATERS labels (skipped)
        !> ================================================================
        else !> label is not "TARGET WATERS"
            continue !> skip to next label - read next iteration
        end if
        
    end do !> end of infinite label reading loop
    
    !> ================================================================
    !> Section 6: File closing and finalization
    !> ================================================================
    close(unit) !> close target waters file (*_tar_wat.dat) - initialization complete
    
    !> (COMMENTED) Old approach storing initial state copies
    !this%waters_init=this%waters !> copy current to initial state
    !print *, this%waters_init(3)%id !> debug output
    !do i=1,this%num_tar_waters !> loop over all waters
    !    call this%waters_init(i)%set_solid_chemistry() !> initialize solid chemistry for initial state
    !this%target_solids_init=this%target_solids !> copy solids to initial state
    !if (allocated(this%target_gases)) then !> if gases are present
    !    this%target_gases_init=this%target_gases !> copy gases to initial state
    !end if
    !call this%write_aq_comps_init(root) !> (commented) write aqueous components of initial target waters
    
    call this%allocate_ext_waters_indices(num_wat_bd+num_wat_rech) !> allocate external waters (boundary + recharge)
    call this%set_ext_waters_indices() !> compute external waters array from boundary and recharge arrays
    this%tar_wat_indices_init=this%tar_wat_indices !> store initial target water indices
end subroutine !> end of read_waters_init subroutine