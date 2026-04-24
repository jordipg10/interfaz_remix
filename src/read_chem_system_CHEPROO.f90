!> \file read_chem_system_CHEPROO.f90
!! \brief Reads a complete chemical system from CHEPROO-formatted input file
!!
!! \details This subroutine parses a CHEPROO-format chemical system database file using a
!! **two-pass reading algorithm**:
!!
!! **First pass**: Iterates through all sections to count species and reactions of each type.
!! This determines allocation sizes for all arrays.
!!
!! **Second pass**: Iterates through sections again to read actual data (names, properties, flags)
!! and populate the chemical system object.
!!
!! After reading, the routine calls specialized database readers (master25, kinetics, Monod) to
!! complete the chemical system definition with stoichiometric coefficients, equilibrium constants,
!! and kinetic parameters.
!!
!! **CHEPROO format structure**:
!! The input file is organized into labeled sections, each starting with a label line and ending
!! with a line containing '*'. The sections are:
!! - PRIMARY AQUEOUS SPECIES: Primary aqueous species (basis species)
!! - AQUEOUS COMPLEXES: Secondary aqueous species (derived from primary via equilibrium reactions)
!! - MINERALS: Mineral species with equilibrium/kinetic and constant/variable activity flags
!! - GASES: Gas species with equilibrium/kinetic and constant/variable activity flags
!! - SURFACE COMPLEXES: Surface complexation and cation exchange species
!! - REDOX REACTIONS: Redox reactions (equilibrium or kinetic)
!! - LINEAR REACTIONS: Linear kinetic reactions (radioactive decay, etc.)
!!
!! Special handling:
!! - h2o(p) or h2o is treated as constant activity (water activity ≈ 1)
!! - h+ is identified as the proton species
!! - Exchangeable cations are parsed from surface complex names (format: X-cat or X2-cat)
!!
!! \param[in,out] this Chemical system object to be populated
!! \param[in] path_DB Path to the directory containing CHEPROO database files
!! \param[in] unit File unit number for reading the CHEPROO system file (already opened)
!!
!! \author Jordi Petchamé-Guerrero
!! \date 2024
!!
!! \see read_master25, read_kinetics_DB, read_Monod_DB
!!
!> Suponemos que el archivo ya ha sido abierto
subroutine read_chem_system_CHEPROO(this,path_DB,unit)
    use chem_system_m, only: chem_system_c
    use species_m, only: species_c
    use aq_species_m, only: aq_species_c
    use mineral_m, only: mineral_c
    use surf_compl_m, only: surface_c
    use gas_species_m, only: gas_species_c
    implicit none
    
!> \name Subroutine arguments
!! @{
    class(chem_system_c) :: this                !< [in,out] Chemical system object to be populated with species and reactions
    character(len=*), intent(in) :: path_DB     !< [in] Path to directory containing database files (master25, kinetics, Monod)
    integer(kind=4), intent(in) :: unit         !< [in] File unit number for reading CHEPROO input file (must be already open)
!! @}
    
!> \name Stoichiometric and thermodynamic arrays
!! @{
    real(kind=8), allocatable :: Sk(:,:)        !< Stoichiometric matrix for kinetic reactions [n_species × n_k] [-]
    real(kind=8), allocatable :: logK(:)        !< Log10 of equilibrium constants for equilibrium reactions [n_eq] [-]
    real(kind=8), allocatable :: gamma_1(:)     !< First Debye-Hückel activity coefficient parameter [n_species] [-]
    real(kind=8), allocatable :: gamma_2(:)     !< Second Debye-Hückel activity coefficient parameter [n_species] [-]
!! @}

!> \name Loop and index counters
!! @{
    integer(kind=4) :: i                        !< General loop counter (reused for different sections) [-]
    integer(kind=4) :: j                        !< Secondary loop counter (reused for different sections) [-]
    integer(kind=4) :: k                        !< Tertiary loop counter (reused for different sections) [-]
    integer(kind=4) :: l                        !< Quaternary loop counter (reused for different sections) [-]
    integer(kind=4) :: ind_var_act_sp           !< Running index for variable activity species during second pass [-]
    integer(kind=4) :: ind_diss_solids          !< Running index for dissolved solids (aqueous species) during second pass [-]
    integer(kind=4) :: exch_cat_ind             !< Index of exchangeable cation in aqueous species array [-]
    integer(kind=4) :: tmp_index                !< Temporary index variable (reused) [-]
!! @}

!> \name Species counters (total counts across all phases)
!! @{
    integer(kind=4) :: num_sp                   !< Total number of species in chemical system (all phases) [-]
    integer(kind=4) :: num_aq_sp                !< Number of aqueous species (primary + complexes) [-]
    integer(kind=4) :: num_sec_aq_sp            !< Number of secondary aqueous species (complexes) [-]
    integer(kind=4) :: num_aq_compl             !< Number of aqueous complexes (same as num_sec_aq_sp) [-]
    integer(kind=4) :: num_var_act_sp           !< Number of species with variable activity [-]
    integer(kind=4) :: num_cst_act_sp           !< Number of species with constant activity (e.g., h2o, pure minerals) [-]
    integer(kind=4) :: num_mins                 !< Total number of mineral species [-]
    integer(kind=4) :: num_mins_eq              !< Number of minerals at equilibrium [-]
    integer(kind=4) :: num_eq_cst_act_mins      !< Number of equilibrium minerals with constant activity [-]
    integer(kind=4) :: num_eq_var_act_mins      !< Number of equilibrium minerals with variable activity [-]
    integer(kind=4) :: num_kin_cst_act_mins     !< Number of kinetic minerals with constant activity [-]
    integer(kind=4) :: num_kin_var_act_mins     !< Number of kinetic minerals with variable activity [-]
    integer(kind=4) :: num_var_act_mins         !< Total number of minerals with variable activity [-]
    integer(kind=4) :: num_gases                !< Total number of gas species [-]
    integer(kind=4) :: num_var_act_gases        !< Number of gases with variable activity (partial pressure) [-]
    integer(kind=4) :: num_cst_act_gases        !< Number of gases with constant activity (atmospheric) [-]
    integer(kind=4) :: num_surf_compl           !< Number of surface complexation species [-]
    integer(kind=4) :: num_exch_cats            !< Number of exchangeable cations [-]
!! @}

!> \name Reaction counters
!! @{
    integer(kind=4) :: n_eq                     !< Total number of equilibrium reactions [-]
    integer(kind=4) :: n_k                      !< Total number of kinetic reactions [-]
    integer(kind=4) :: n_eq_homog               !< Number of homogeneous equilibrium reactions [-]
    integer(kind=4) :: n_min_kin                !< Number of kinetic mineral reactions (dissolution/precipitation) [-]
    integer(kind=4) :: n_gas_eq                 !< Total number of gas equilibrium reactions [-]
    integer(kind=4) :: n_gas_eq_cst_act         !< Number of gas equilibrium reactions with constant activity [-]
    integer(kind=4) :: n_gas_eq_var_act         !< Number of gas equilibrium reactions with variable activity [-]
    integer(kind=4) :: n_gas_kin                !< Total number of gas kinetic reactions [-]
    integer(kind=4) :: n_gas_kin_cst_act        !< Number of gas kinetic reactions with constant activity [-]
    integer(kind=4) :: n_gas_kin_var_act        !< Number of gas kinetic reactions with variable activity [-]
    integer(kind=4) :: n_redox                  !< Total number of redox reactions (eq + kin) [-]
    integer(kind=4) :: n_redox_eq               !< Number of redox equilibrium reactions [-]
    integer(kind=4) :: n_redox_kin              !< Number of redox kinetic reactions (Monod) [-]
    integer(kind=4) :: n_lin_eq                 !< Number of linear equilibrium reactions [-]
    integer(kind=4) :: n_lin_kin                !< Number of linear kinetic reactions (e.g., radioactive decay) [-]
    integer(kind=4) :: n_r                      !< Total number of reactions (temporary variable) [-]
!! @}

!> \name Index arrays
!! @{
    integer(kind=4), allocatable :: n_tar(:)                !< Number of target species for reactions [n_reactions] [-]
    integer(kind=4), allocatable :: mins_eq_indices(:)      !< Indices of equilibrium minerals [num_mins_eq] [-]
    integer(kind=4), allocatable :: gases_eq_indices(:)     !< Indices of equilibrium gases [n_gas_eq] [-]
    integer(kind=4), allocatable :: indices_lin_reacts(:,:) !< Indices of species in linear reactions [n_lin_kin+n_lin_eq, 2] [-]
    integer(kind=4) :: num_mins_eq_indices                  !< Counter for equilibrium mineral indices [-]
!! @}

!> \name Temporary reaction and species properties
!! @{
    integer(kind=4) :: exch_cat_valence         !< Valence of exchangeable cation [-]
    integer(kind=4) :: kin_react_type           !< Type of kinetic reaction (1=mineral, 4=redox, 5=linear) [-]
    real(kind=8) :: aux                         !< Auxiliary variable for temporary calculations (reused) [-]
    real(kind=8) :: conc                        !< Concentration (temporary, for reading) [M]
    real(kind=8) :: temp                        !< Temperature (temporary, for reading) [K]
    real(kind=8) :: SI                          !< Saturation index (temporary, for reading) [-]
    real(kind=8) :: lambda                      !< Decay constant for linear kinetic reactions [1/T]
    real(kind=8) :: yield                       !< Yield coefficient for redox reactions [-]
!! @}

!> \name String variables for reading and parsing
!! @{
    character(len=256) :: str                   !< General string for reading lines from file
    character(len=256) :: str1                  !< First string in multi-column reads (species name, reaction name)
    character(len=256) :: str2                  !< Second string in multi-column reads (species name, product name)
    character(len=256) :: str3                  !< Third string in multi-column reads
    character(len=256) :: str4                  !< Fourth string in multi-column reads
    character(len=256) :: str5                  !< Fifth string in multi-column reads
    character(len=256) :: Monod_name            !< Name of Monod kinetic reaction (redox)
    character(len=256) :: file_kin_params       !< Filename for kinetic parameters database
    character(len=256) :: label                 !< Section label in CHEPROO file (e.g., "PRIMARY AQUEOUS SPECIES")
    character(len=:), allocatable :: str_block_trim    !< Trimmed string for block reading
    character(len=:), allocatable :: str_trim          !< Trimmed string for species/mineral/gas names
    character(len=:), allocatable :: valence_str       !< String representation of valence
    character(len=:), allocatable :: exch_cat_val      !< String representation of exchangeable cation valence
    character(len=:), allocatable :: exch_cat_name     !< Name of exchangeable cation (parsed from surface complex)
    character(len=:), allocatable :: str_trim_1        !< First trimmed string (for linear reactions)
    character(len=:), allocatable :: str_trim_2        !< Second trimmed string (for linear reactions)
!! @}

!> \name String arrays for species names
!! @{
    character(len=256), allocatable :: aq_species_str(:)     !< Array of aqueous species names [num_aq_sp]
    character(len=256), allocatable :: prim_species_str(:)   !< Array of primary species names [num_primary]
    character(len=256), allocatable :: cst_act_species_str(:)!< Array of constant activity species names [num_cst_act_sp]
    character(len=256), allocatable :: minerals_str(:)       !< Array of mineral names [num_mins]
    character(len=256), allocatable :: solid_species_str(:)  !< Array of solid species names [num_solids]
    character(len=256), allocatable :: kin_react_names(:)    !< Array of kinetic reaction names [n_k]
!! @}

!> \name Logical flags for reading and classification
!! @{
    logical :: flag                             !< General flag for existence checks
    logical :: eq_label                         !< Flag indicating if reaction is at equilibrium (.true.) or kinetic (.false.)
    logical :: exch_cat_flag                    !< Flag indicating if species is an exchangeable cation
    logical :: cst_act_label                    !< Flag indicating if species has constant activity
    logical :: flag_1                           !< First flag for species search in linear reactions
    logical :: flag_2                           !< Second flag for species search in linear reactions
    logical :: flag_surf                        !< Flag indicating if surface complexation is present
!! @}

!> \name Derived type objects and arrays for chemical species
!! @{
    type(species_c) :: species                  !< Temporary species object for reading
    type(species_c), allocatable :: aux_Species(:)              !< Auxiliary species object for rearranging indices
    type(species_c), allocatable :: surf_compl(:) !< Array of surface complexation species [num_surf_compl]
    type(aq_species_c) :: exch_cat              !< Temporary exchangeable cation object
    type(aq_species_c) :: aq_sp_1               !< First aqueous species in linear reactions
    type(aq_species_c) :: aq_sp_2               !< Second aqueous species in linear reactions
    type(aq_species_c), allocatable :: aq_species(:)  !< Array of aqueous species objects [num_aq_sp]
    type(aq_species_c), allocatable :: exch_cats(:)   !< Array of exchangeable cation objects [num_exch_cats]
    type(aq_species_c), allocatable :: prim_species(:)!< Array of primary species objects [num_primary]
    type(mineral_c) :: mineral                  !< Temporary mineral object for reading
    type(mineral_c), allocatable :: mins(:)     !< Array of mineral objects [num_mins]
    type(surface_c) :: cat_exch_obj             !< Temporary cation exchange object
    type(gas_species_c) :: gas                          !< Temporary gas object for reading
    type(gas_species_c), allocatable :: gases(:)        !< Array of gas objects [num_gases]
    !class(kin_params_c), pointer :: p_kin_params=>null()     !< Pointer to kinetic parameters object
    !class(kin_reaction_c), pointer :: p_kin_react=>null()    !< Pointer to kinetic reaction object
    !type(kin_reaction_poly_c) :: kin_react_ptr                !< Kinetic reaction pointer object (polymorphic)
    !class(kin_reaction_poly_c), allocatable :: kin_reacts(:) !< Array of kinetic reaction objects [n_k]
    !type(eq_reaction_c) :: eq_react                           !< Temporary equilibrium reaction object
    !type(eq_reaction_c), allocatable :: eq_reacts(:)          !< Array of equilibrium reaction objects [n_eq]
!! @}
    
!> \name File unit numbers for database reads
!! @{
    integer(kind=4) :: unit_master_25           !< File unit for master25 database (equilibrium reactions)
    integer(kind=4) :: unit_kinetics            !< File unit for kinetics database (mineral reactions)
    integer(kind=4) :: unit_redox               !< File unit for Monod database (redox reactions)
!! @}
    
    logical :: flag_comp                        !< Flag for component speciation algebra (set to .false. initially)
    !logical :: flag_surf                        !< Flag for surface complexation (set based on num_surf_compl>0)
    integer(kind=4) :: wat_flag                 !< Flag indicating water presence (0=no water, 1=water present)
    
!> \subsection init_counters Initialization of all counters to zero
!! All species and reaction counters are initialized to zero before first pass through CHEPROO file.
!! This ensures proper counting during the first iteration.
    num_sp=0                                    !< Initialize total species counter
    num_aq_sp=0                                 !< Initialize aqueous species counter
    num_aq_compl=0                              !< Initialize aqueous complexes counter
    num_cst_act_sp=0                            !< Initialize constant activity species counter
    num_cst_act_gases=0                         !< Initialize constant activity gases counter
    num_var_act_gases=0                         !< Initialize variable activity gases counter
    num_eq_cst_act_mins=0                       !< Initialize equilibrium constant activity minerals counter
    num_eq_var_act_mins=0                       !< Initialize equilibrium variable activity minerals counter
    num_kin_cst_act_mins=0                      !< Initialize kinetic constant activity minerals counter
    num_kin_var_act_mins=0                      !< Initialize kinetic variable activity minerals counter
    num_var_act_sp=0                            !< Initialize variable activity species counter
    num_surf_compl=0                            !< Initialize surface complexation species counter
    num_exch_cats=0                             !< Initialize exchangeable cations counter
    num_gases=0                                 !< Initialize gas species counter
    num_mins=0                                  !< Initialize mineral species counter
    num_mins_eq=0                               !< Initialize equilibrium minerals counter
    n_eq=0                                      !< Initialize total equilibrium reactions counter
    n_k=0                                       !< Initialize total kinetic reactions counter
    n_min_kin=0                                 !< Initialize kinetic mineral reactions counter
    n_gas_eq=0                                  !< Initialize gas equilibrium reactions counter
    n_gas_eq_cst_act=0                          !< Initialize gas equilibrium constant activity counter
    n_gas_eq_var_act=0                          !< Initialize gas equilibrium variable activity counter
    n_gas_kin_cst_act=0                         !< Initialize gas kinetic constant activity counter
    n_gas_kin_var_act=0                         !< Initialize gas kinetic variable activity counter
    n_gas_kin=0                                 !< Initialize gas kinetic reactions counter
    n_lin_kin=0                                 !< Initialize linear kinetic reactions counter
    n_lin_eq=0                                  !< Initialize linear equilibrium reactions counter
    n_redox_eq=0                                !< Initialize redox equilibrium reactions counter
    n_redox_kin=0                               !< Initialize redox kinetic reactions counter
    
!> \subsection first_pass First iteration through CHEPROO file (counting phase)
!! This loop reads through all sections of the CHEPROO file to count the number of species and
!! reactions of each type. This information is used to allocate arrays in the chemical system
!! object before the second pass fills in the actual data.
!!
!! The loop continues until it encounters the 'end' label. For each section label, it reads
!! the species/reactions within that section (until encountering '*') and increments the
!! appropriate counters based on flags (eq_label, cst_act_label).
    do                                          !< Start first iteration loop
        read(unit,*) label                      !< Read section label
        if (label=='end') then                  !< Check for end of file marker
            rewind(unit)                        !< Rewind file to beginning for second pass
            exit                                !< Exit first iteration loop
        
!> \name PRIMARY AQUEOUS SPECIES section (first pass)
!! Reads primary aqueous species and counts them. Special handling for water (h2o or h2o(p))
!! which is treated as constant activity species (activity ≈ 1).
!! All other primary species are variable activity and also counted as total aqueous species.
        else if (label=='PRIMARY AQUEOUS SPECIES') then
            do                                  !< Loop through primary aqueous species
                read(unit,*) str                !< Read species name
                if (str=='*') then              !< Check for section end marker
                    exit                        !< Exit species loop
                else if (str=='h2o(p)') then    !< Special case: water in primary species
                !> We assume water has constant activity (activity ≈ 1 for dilute solutions)
                    num_cst_act_sp=num_cst_act_sp+1        !< Increment constant activity species counter
                    this%aq_phase%wat_flag=1                !< Set water flag (1 = present)
                    this%aq_phase%ind_wat=num_aq_sp+1      !< Store water index (will be next species)
                else                            !< All other primary species
                    num_var_act_sp=num_var_act_sp+1        !< Increment variable activity species counter
                end if
                num_aq_sp=num_aq_sp+1           !< Increment total aqueous species counter
                num_sp=num_sp+1                 !< Increment total species counter
            end do
        
!> \name AQUEOUS COMPLEXES section (first pass)
!! Reads aqueous complexes (secondary species formed from primary via equilibrium reactions).
!! Special handling for h2o (alternative water notation) as constant activity.
!! Each complex corresponds to one equilibrium reaction.
        else if (label=='AQUEOUS COMPLEXES') then
            do                                  !< Loop through aqueous complexes
                read(unit,*) str                !< Read complex name
                if (str=='*') then              !< Check for section end marker
                    exit                        !< Exit complexes loop
                else if (str=='h2o') then       !< Alternative water notation
                !> We assume water has constant activity (activity ≈ 1)
                    num_cst_act_sp=num_cst_act_sp+1        !< Increment constant activity species counter
                    this%aq_phase%wat_flag=1                !< Set water flag (1 = present)
                    this%aq_phase%ind_wat=num_aq_sp+1      !< Store water index
                else                            !< All other aqueous complexes
                    num_var_act_sp=num_var_act_sp+1        !< Increment variable activity species counter
                end if
                num_aq_sp=num_aq_sp+1           !< Increment total aqueous species counter
                num_aq_compl=num_aq_compl+1     !< Increment aqueous complexes counter
                num_sp=num_sp+1                 !< Increment total species counter
                n_eq=n_eq+1                     !< Increment equilibrium reactions counter (one per complex)
            end do
            
!> \name MINERALS section (first pass)
!! Reads mineral species with two flags:
!! - eq_label: .true. = equilibrium, .false. = kinetic
!! - cst_act_label: .true. = constant activity (pure solid), .false. = variable activity
!! This gives 4 combinations that are counted separately for proper array allocation.
        else if (label=='MINERALS') then
            do                                  !< Loop through minerals
                read(unit,*) str, eq_label, cst_act_label  !< Read mineral name and flags
                if (str=='*') exit              !< Check for section end marker
                num_mins=num_mins+1             !< Increment total minerals counter
                num_sp=num_sp+1                 !< Increment total species counter
                
                !> Case 1: Equilibrium mineral with constant activity (pure solid at equilibrium)
                if (cst_act_label.eqv..true. .AND. eq_label.eqv..true.) then
                    n_eq=n_eq+1                         !< Increment equilibrium reactions counter
                    num_mins_eq=num_mins_eq+1           !< Increment equilibrium minerals counter
                    num_eq_cst_act_mins=num_eq_cst_act_mins+1  !< Increment eq const act minerals
                    num_cst_act_sp=num_cst_act_sp+1    !< Increment constant activity species counter
                
                !> Case 2: Equilibrium mineral with variable activity (solid solution at equilibrium)
                else if (cst_act_label.eqv..false. .AND. eq_label.eqv..true.) then
                    n_eq=n_eq+1                         !< Increment equilibrium reactions counter
                    num_mins_eq=num_mins_eq+1           !< Increment equilibrium minerals counter
                    num_eq_var_act_mins=num_eq_var_act_mins+1  !< Increment eq var act minerals
                    num_var_act_sp=num_var_act_sp+1    !< Increment variable activity species counter
                
                !> Case 3: Kinetic mineral with constant activity (pure solid with kinetic dissolution/precipitation)
                else if (cst_act_label.eqv..true. .AND. eq_label.eqv..false.) then
                    num_cst_act_sp=num_cst_act_sp+1    !< Increment constant activity species counter
                    num_kin_cst_act_mins=num_kin_cst_act_mins+1  !< Increment kin const act minerals
                    n_k=n_k+1                           !< Increment kinetic reactions counter
                    n_min_kin=n_min_kin+1               !< Increment kinetic minerals counter
                
                !> Case 4: Kinetic mineral with variable activity (solid solution with kinetic reactions)
                else
                    num_kin_var_act_mins=num_kin_var_act_mins+1  !< Increment kin var act minerals
                    num_var_act_sp=num_var_act_sp+1    !< Increment variable activity species counter
                    n_k=n_k+1                           !< Increment kinetic reactions counter
                    n_min_kin=n_min_kin+1               !< Increment kinetic minerals counter
                end if
            end do
            
!> \name GASES section (first pass)
!! Reads gas species with two flags (same structure as MINERALS):
!! - eq_label: .true. = equilibrium (Henry's law), .false. = kinetic
!! - cst_act_label: .true. = constant activity (atmospheric), .false. = variable activity
!! Four combinations are counted separately for proper array allocation.
        else if (label=='GASES') then
            do                                  !< Loop through gases
                read(unit,*) str, eq_label, cst_act_label  !< Read gas name and flags
                if (str=='*') exit              !< Check for section end marker
                num_gases=num_gases+1           !< Increment total gases counter
                num_sp=num_sp+1                 !< Increment total species counter
                
                !> Case 1: Equilibrium gas with constant activity (atmospheric gas at equilibrium)
                if (cst_act_label.eqv..true. .AND. eq_label.eqv..true.) then
                    n_eq=n_eq+1                         !< Increment equilibrium reactions counter
                    n_gas_eq=n_gas_eq+1                 !< Increment gas equilibrium reactions counter
                    n_gas_eq_cst_act=n_gas_eq_cst_act+1  !< Increment gas eq const act counter
                    num_cst_act_gases=num_cst_act_gases+1
                    num_cst_act_sp=num_cst_act_sp+1
                else if (cst_act_label.eqv..false. .AND. eq_label.eqv..true.) then
                    n_eq=n_eq+1
                    n_gas_eq=n_gas_eq+1
                    n_gas_eq_var_act=n_gas_eq_var_act+1
                    num_var_act_gases=num_var_act_gases+1
                    num_var_act_sp=num_var_act_sp+1
                else if (cst_act_label.eqv..true. .AND. eq_label.eqv..false.) then 
                    n_gas_kin=n_gas_kin+1
                    n_gas_kin_cst_act=n_gas_kin_cst_act+1
                    num_cst_act_gases=num_cst_act_gases+1
                    num_cst_act_sp=num_cst_act_sp+1
                else
                    n_gas_kin=n_gas_kin+1
                    n_gas_kin_var_act=n_gas_kin_var_act+1
                    num_var_act_gases=num_var_act_gases+1
                    num_var_act_sp=num_var_act_sp+1
                end if
            end do 
        else if (label=='SURFACE COMPLEXES') then
            do
                read(unit,*) str
                if (str=='*') exit
                str_trim=trim(str)
                num_surf_compl=num_surf_compl+1
                num_var_act_sp=num_var_act_sp+1
                num_sp=num_sp+1
                exch_cat_ind=index(str_trim,'-')
                if (exch_cat_ind>0 .and. exch_cat_ind<len(str_trim)) then
                    num_exch_cats=num_exch_cats+1
                    n_eq=n_eq+1
                end if
            end do
        else if (label=='REDOX REACTIONS') then
            do
                read(unit,*) str, eq_label
                if (str=='*') then
                    exit
                !else if (str=='aerobic degradation DOC (review)') then
                !    n_redox_eq=n_redox_eq+1
                !    n_eq=n_eq+1
                !else
                !    n_redox_kin=n_redox_kin+1
                !    n_k=n_k+1
                end if
                if (eq_label.eqv..true.) then
                    n_redox_eq=n_redox_eq+1
                    n_eq=n_eq+1
                else
                    n_redox_kin=n_redox_kin+1
                    n_k=n_k+1
                end if
            end do
        else if (label=='LINEAR REACTIONS') then
            do
                read(unit,*) str1, str2, lambda, eq_label
                if (str1=='*') exit
                if (eq_label.eqv..true.) then
                    n_lin_eq=n_lin_eq+1
                    n_eq=n_eq+1
                else
                    n_lin_kin=n_lin_kin+1
                    n_k=n_k+1
                end if
            end do
        else 
            continue
        end if
    end do
!> We set & allocate attributes
    !> Aqueous phase
    call this%aq_phase%allocate_aq_species(num_aq_sp)
    call this%aq_phase%allocate_ind_diss_solids()
    call this%aq_phase%set_num_aq_complexes(num_aq_compl)
    call this%aq_phase%set_num_var_act_species_phase(num_aq_sp-this%aq_phase%wat_flag)
    call this%aq_phase%set_num_cst_act_species_phase(this%aq_phase%wat_flag)
    !> Gas phase
    call this%gas_phase%allocate_gases(num_gases)
    call this%gas_phase%set_num_gases_eq(n_gas_eq) !> 
    call this%gas_phase%set_num_gases_eq_cst_act(n_gas_eq_cst_act) !> 
    call this%gas_phase%set_num_gases_eq_var_act(n_gas_eq_var_act) !> 
    call this%gas_phase%set_num_gases_kin_cst_act(n_gas_kin_cst_act) !> 
    call this%gas_phase%set_num_gases_kin_var_act(n_gas_kin_var_act) !> 
    call this%gas_phase%set_num_gases_kin(n_gas_kin) !> 
    call this%gas_phase%set_num_var_act_species_phase(num_var_act_gases)
    call this%gas_phase%set_num_cst_act_species_phase(num_cst_act_gases)
    !> Cation exchange
    call this%cat_exch_zone%allocate_surf_compl(num_surf_compl)
    call this%cat_exch_zone%allocate_exch_cat_indices(num_exch_cats)
    call this%cat_exch_zone%set_num_exch_cats(num_exch_cats)
    !> Chemical system
    call this%allocate_cst_act_sp_indices(num_cst_act_sp)
    call this%allocate_var_act_sp_indices(num_var_act_sp)
    call this%allocate_species()
    call this%allocate_minerals(num_mins)
    call this%set_num_minerals_eq(num_mins_eq)
    call this%set_num_minerals_eq_cst_act(num_eq_cst_act_mins)
    call this%set_num_minerals_eq_var_act(num_eq_var_act_mins)
    call this%allocate_reacts(n_eq,n_k)
    call this%set_num_redox_eq_reacts(n_redox_eq)
    call this%set_num_aq_eq_reacts()
    call this%allocate_min_kin_reacts(n_min_kin)
    call this%set_num_minerals_kin_cst_act(num_kin_cst_act_mins)
    call this%set_num_minerals_kin_var_act(num_kin_var_act_mins)
    call this%allocate_redox_kin_reacts(n_redox_kin)
    call this%allocate_lin_kin_reacts(n_lin_kin)
    call this%set_num_aq_kin_reacts()
    call this%compute_num_reacts()
    call this%compute_num_solids_chem_syst()
    !> Speciatrion algebra    
    call this%speciation_alg%set_flag_comp(.false.) !> chapuza
    if (num_surf_compl>0) then
        flag_surf=.true.
    else
        flag_surf=.false.
    end if
    call this%speciation_alg%set_flag_cat_exch(flag_surf)
    call this%speciation_alg%set_dimensions(this%speciation_alg%num_species,this%speciation_alg%num_eq_reactions,this%num_cst_act_species,this%aq_phase%num_species,&
        this%aq_phase%num_species-this%aq_phase%wat_flag,this%num_minerals_kin,num_gases-n_gas_eq)
!> Second iteration of chemical system file
    ind_var_act_sp=0 !< counter variable activity species
    ind_diss_solids=0 !> counter dissolved solids
    i=0 !> counter aqueous species
    do
        read(unit,*) label
        if (label=='end') then
            exit
        else if (label=='PRIMARY AQUEOUS SPECIES' .OR. label=='AQUEOUS COMPLEXES') then !> suponemos ordenadas en primarias y secundarias
            do
                read(unit,*) str
                if (str=='*') exit
                i=i+1
                str_trim=trim(str)
                call this%aq_phase%aq_species(i)%set_name(str_trim)
                if (str=='h2o(p)' .or. str=='h2o') then
                    call this%aq_phase%aq_species(i)%set_cst_act_flag(.true.)
                    this%cst_act_sp_indices(1)=i
                else
                    if (str=='h+') then
                        call this%aq_phase%set_ind_prot(i)
                    end if
                    ind_var_act_sp=ind_var_act_sp+1
                    ind_diss_solids=ind_diss_solids+1
                    call this%aq_phase%aq_species(i)%set_cst_act_flag(.false.)
                    this%var_act_sp_indices(ind_var_act_sp)=i
                    this%aq_phase%ind_diss_solids(ind_diss_solids)=i
                end if
                call this%species(i)%copy_species(this%aq_phase%aq_species(i))
            end do
        else if (label=='MINERALS') then
            i=0 !> counter minerals equilibrium var activity
            j=0 !> counter minerals kinetic variable activity
            k=0 !> counter minerals equilibrium constant activity
            l=0 !> counter minerals kinetic constant activity
            do
                read(unit,*) str, eq_label, cst_act_label
                if (str=='*') exit
                str_trim=trim(str)
                if ((eq_label .eqv. .true.) .and. (cst_act_label .eqv. .false.)) then !> equilibrium variable activity
                    i=i+1
                    call this%minerals(this%num_minerals_kin+i)%set_phase_name(str_trim)
                    call this%minerals(this%num_minerals_kin+i)%mineral%set_name(str_trim)
                    call this%minerals(this%num_minerals_kin+i)%mineral%set_cst_act_flag(cst_act_label)
                    call this%minerals(this%num_minerals_kin+i)%mineral%set_valence(0)
                else if (eq_label .eqv. .true. .and. cst_act_label .eqv. .true.) then !> equilibrium constant activity
                    k=k+1
                    call this%minerals(this%num_minerals-this%num_minerals_eq_cst_act+k)%set_phase_name(str_trim)
                    call this%minerals(this%num_minerals-this%num_minerals_eq_cst_act+k)%mineral%set_name(str_trim)
                    call this%minerals(this%num_minerals-this%num_minerals_eq_cst_act+k)%mineral%set_cst_act_flag(cst_act_label)
                    call this%minerals(this%num_minerals-this%num_minerals_eq_cst_act+k)%mineral%set_valence(0)
                else if (eq_label .eqv. .false. .and. cst_act_label .eqv. .false.) then !> kinetic variable activity
                    j=j+1
                    call this%minerals(j)%set_phase_name(str_trim)
                    call this%minerals(j)%mineral%set_name(str_trim)
                    call this%minerals(j)%mineral%set_cst_act_flag(cst_act_label)
                    call this%minerals(j)%mineral%set_valence(0)
                else !> kinetic constant activity
                    l=l+1
                    call this%minerals(this%num_minerals_kin_var_act+l)%set_phase_name(str_trim)
                    call this%minerals(this%num_minerals_kin_var_act+l)%mineral%set_name(str_trim)
                    call this%minerals(this%num_minerals_kin_var_act+l)%mineral%set_cst_act_flag(cst_act_label)
                    call this%minerals(this%num_minerals_kin_var_act+l)%mineral%set_valence(0)
                end if
            end do
        else if (label=='GASES') then
            i=0 !> counter gases equilibrium with constant activity
            j=0 !> counter gases equilibrium with variable activity
            k=0 !> counter gases kinetic with constant activity
            l=0 !> counter gases equilibrium with variable activity
            !> first gases in equilibrium with constant activity, then equilibrium with variable activity, & 
            !! then kinetic with constant activity, then kinetic with variable activity
            do
                read(unit,*) str, eq_label, cst_act_label
                if (str=='*') exit
                str_trim=trim(str)
                if ((eq_label .eqv. .true.) .and. (cst_act_label .eqv. .true.)) then
                    i=i+1
                    call this%gas_phase%gases(i)%set_name(str_trim)
                    call this%gas_phase%gases(i)%set_cst_act_flag(cst_act_label)
                    call this%gas_phase%gases(i)%set_valence(0)
                else if ((eq_label .eqv. .true.) .and. (cst_act_label .eqv. .false.)) then
                    j=j+1
                    call this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+j)%set_name(str_trim)
                    call this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+j)%set_cst_act_flag(cst_act_label)
                    call this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+j)%set_valence(0)
                else if ((eq_label .eqv. .false.) .and. (cst_act_label .eqv. .true.)) then
                    k=k+1
                    call this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+this%gas_phase%num_gases_eq_var_act+k)%set_name(&
                        str_trim)
                    call this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+&
                        this%gas_phase%num_gases_eq_var_act+k)%set_cst_act_flag(cst_act_label)
                    call this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+this%gas_phase%num_gases_eq_var_act+&
                        k)%set_valence(0)
                else
                    l=l+1
                    call this%gas_phase%gases(this%gas_phase%num_species-this%gas_phase%num_gases_kin_var_act+l)%set_name(str_trim)
                    call this%gas_phase%gases(this%gas_phase%num_species-this%gas_phase%num_gases_kin_var_act+l)%set_cst_act_flag(&
                        cst_act_label)
                    call this%gas_phase%gases(this%gas_phase%num_species-this%gas_phase%num_gases_kin_var_act+l)%set_valence(0)
                end if
            end do
        else if (label=='SURFACE COMPLEXES') then
            i=0 !> counter surface complexes
            j=0 !> counter exchangeable cations
            do
                read(unit,*) str
                if (str=='*') exit
                i=i+1
                str_trim=trim(str)
                call this%cat_exch_zone%surf_compl(i)%set_name(str_trim)
                call this%cat_exch_zone%surf_compl(i)%set_cst_act_flag(.false.)
                exch_cat_ind=index(str_trim,'-')
                if (exch_cat_ind>0 .and. exch_cat_ind<len(str_trim)) then !> exchangeable cation
                    j=j+1
                    if (exch_cat_ind>1) then
                        if (exch_cat_ind==2) then
                            exch_cat_valence=1
                            allocate(character(len=(len(str_trim)-exch_cat_ind+1)) :: exch_cat_name)
                            write(exch_cat_name(1:len(str_trim)-exch_cat_ind),"(A)") str_trim(exch_cat_ind+1:len(str_trim))
                            write(exch_cat_name(len(str_trim)-exch_cat_ind+1:len(str_trim)-exch_cat_ind+1),"(A1)") '+'
                        else
                            exch_cat_val=trim(str(2:exch_cat_ind-1))
                            read(exch_cat_val,*) exch_cat_valence
                            allocate(character(len=(len(str_trim)-exch_cat_ind+1+len(exch_cat_val))) :: exch_cat_name)
                            write(exch_cat_name(1:len(str_trim)-exch_cat_ind),"(A)") str_trim(exch_cat_ind+1:len(str_trim))
                            write(exch_cat_name(len(str_trim)-exch_cat_ind+1:len(str_trim)-exch_cat_ind+1),"(A1)") '+'
                            write(exch_cat_name(len(exch_cat_name)-len(exch_cat_val)+1:len(exch_cat_name)),"(A)") exch_cat_val
                        end if
                    end if
                    call exch_cat%set_name(exch_cat_name)
                    call exch_cat%set_valence(exch_cat_valence)
                    call this%aq_phase%is_species_in_aq_phase(exch_cat,exch_cat_flag,exch_cat_ind)
                    if (exch_cat_flag.eqv..true.) then
                        this%cat_exch_zone%exch_cat_indices(j)=exch_cat_ind
                    else
                        error stop "Exchangeable cation is not in the chemical system"
                    end if
                    deallocate(exch_cat_name)
                else if (exch_cat_ind>0) then !> free surface complex
                    ! call this%cat_exch%exch_cats(j)%set_name(str_trim)
                    ! call this%cat_exch%exch_cats(j)%set_cst_act_flag(.false.)
                    call this%cat_exch_zone%surf_compl(i)%set_valence(-1)
                else
                    error stop "Surface complex not well formatted"
                end if
                
            end do
        else if (label=='REDOX REACTIONS') then
            i=0 !> counter redox kinetic reactions
            j=0 !> counter redox equilibrium reactions
            do
                read(unit,*) str, eq_label!, yield !> reads name of reaction, equilibrium flag, and yield
                if (str=='*') then
                    exit
                end if
                str_trim=trim(str)
                if (eq_label.eqv..false.) then !> redox kinetic reactions
                    i=i+1
                    call this%kin_reacts(this%num_lin_kin_reacts+i)%set_kin_reaction(this%redox_kin_reacts(i))
                    call this%redox_kin_reacts(i)%set_react_name(str_trim)
                    call this%redox_kin_reacts(i)%set_react_type(4)
                    if (str=='aerobic degradation DOC (review)') then !> chapuza
                        call this%redox_kin_reacts(i)%set_yield(yield)
                    end if
                else !> redox eq reacts (they are the first in equilibrium reactions attribute of chemical system)
                    j=j+1
                    call this%eq_reacts(j)%set_react_name(str_trim)
                    call this%eq_reacts(j)%set_react_type(4)
                end if
            end do
        else if (label=='LINEAR REACTIONS') then
            i=0 !> counter linear kinetic reactions
            j=0 !> counter linear equilibrium reactions
            allocate(indices_lin_reacts(n_lin_kin+n_lin_eq,2)) !> first kinetic, then equilibrium
            do
                read(unit,*) str1, str2, lambda, eq_label
                if (str1=='*') exit
                str_trim_1=trim(str1)
                str_trim_2=trim(str2)
                call aq_sp_1%set_name(str_trim_1)
                call aq_sp_2%set_name(str_trim_2)
                if (eq_label.eqv..false.) then
                    i=i+1
                    call this%aq_phase%is_species_in_aq_phase(aq_sp_1,flag_1,indices_lin_reacts(i,1))
                    call this%aq_phase%is_species_in_aq_phase(aq_sp_2,flag_2,indices_lin_reacts(i,2))
                    if (flag_1.eqv..false. .or. flag_2.eqv..false.) then
                        error stop "Species in linear kinetic reaction not found in chemical system"
                    else
                        call this%kin_reacts(i)%set_kin_reaction(this%lin_kin_reacts(i))
                        call this%lin_kin_reacts(i)%set_react_name(str_trim_1//'-->'//str_trim_2)
                        call this%lin_kin_reacts(i)%set_react_type(5)
                        call this%lin_kin_reacts(i)%allocate_reaction(2)
                        call this%lin_kin_reacts(i)%set_lambda(lambda)
                    end if
                else
                    !j=j+1
                    !call this%eq_reacts(j)%set_react_name(str_trim_1//str_trim_2)
                    !call this%eq_reacts(j)%set_react_type(5)
                end if
            end do
        else 
            continue
        end if
    end do
!> We read databases
    !> Master25 (equilibrium reactions)
    unit_master_25=2
    call this%read_master25(path_DB,unit_master_25)
    !> Kinetics (minerals)
    unit_kinetics=3
    if (this%num_minerals_kin>0) then
        call this%read_kinetics_DB(path_DB,unit_kinetics)
    end if
    !> Monod reactions
    unit_redox=4
    if (this%num_redox_kin_reacts>0 .or. this%num_redox_eq_reacts>0) then
        call this%read_Monod_DB(path_DB,unit_redox)
    end if
    !> Linear kinetcir eactions
    if (this%num_lin_kin_reacts>0) then
        do i=1,this%num_lin_kin_reacts
            call this%is_species_in_chem_syst(this%aq_phase%aq_species(indices_lin_reacts(i,1)),flag_1,indices_lin_reacts(i,1))
            call this%is_species_in_chem_syst(this%aq_phase%aq_species(indices_lin_reacts(i,2)),flag_2,indices_lin_reacts(i,2))
            call this%lin_kin_reacts(i)%set_all_species(indices_lin_reacts(i,:))
            call this%lin_kin_reacts(i)%set_stoichiometry([-1d0,1d0])
        end do
    end if
!> We rearrange species and equilibrium reactions
    aux_Species=this%species
    call this%rearrange_species()
    call this%compute_z2()
    !if (this%speciation_alg%flag_comp .eqv. .true.) then
        call this%rearrange_eq_reacts()
    !end if
        
        do i=1,this%speciation_alg%num_eq_reactions
            call this%eq_reacts(i)%rearrange_species_indices(aux_Species,this%species)
        end do
        do i=1,this%num_minerals_kin
            call this%min_kin_reacts(i)%rearrange_species_indices(aux_Species,this%species)
        end do
        do i=1,this%num_redox_kin_reacts
            call this%redox_kin_reacts(i)%rearrange_species_indices(aux_Species,this%species)
        end do
        do i=1,this%num_lin_kin_reacts
            call this%lin_kin_reacts(i)%rearrange_species_indices(aux_Species,this%species)
        end do
!> We set stoichiometric matrices
    call this%set_stoich_mat()
    call this%set_stoich_mat_gas()
    call this%set_stoich_mat_sol()
end subroutine