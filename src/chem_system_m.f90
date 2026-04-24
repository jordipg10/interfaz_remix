!> @file chem_system_m.f90
!> @brief Chemical system module containing the chemical theory (like a "chemistry book")
!> @details This module defines the chemical system class which contains species, phases, and reactions.
!> It provides the fundamental chemical framework including aqueous species, minerals, gases, and reactions.
module chem_system_m
    use eq_reaction_m, only: eq_reaction_c !< Import equilibrium reaction classes
    use kin_reaction_m, only: kin_reaction_c, kin_reaction_poly_c !< Import kinetic reaction classes
    use lin_kin_reaction_m, only: lin_kin_reaction_c !< Import linear kinetic reaction class
    use redox_kin_reaction_m, only: redox_kin_c !< Import redox kinetic reaction and Monod parameters
    use kin_mineral_m, only: kin_mineral_c !< Import kinetic mineral and mineral classes
    use mineral_m, only: mineral_c !< Import mineral class
    use aq_phase_m, only: aq_phase_c !< Import aqueous phase class
    use aq_species_m, only: aq_species_c !< Import aqueous phase and species classes
    use species_m, only: species_c !< Import species class and comparison function
    use gas_phase_m, only: gas_phase_c !< Import gas phase and gas classes
    use surf_compl_m, only: cat_exch_zone_c !< Import surface complexation classes
    use speciation_algebra_m, only: speciation_algebra_s !< Import speciation algebra class
    use params_spec_vol_m, only: params_spec_vol_Redlich_c !< Import specific volume parameters
    implicit none !< Declare implicit none to avoid implicit typing
    save !< Save all module variables between calls
    private !< Default accessibility is private
    !> @class chem_system_c
    !> @brief Chemical system class
    !> @details This type contains all chemical information including species, reactions, phases, and stoichiometric matrices
    type, public :: chem_system_c
        integer(kind=4) :: num_reacts !< Total number of reactions (equilibrium + kinetic)
        type(species_c), allocatable :: species(:)  !< Array of all species in the system
        real(kind=8), allocatable :: z2(:) !< Squared charges (valences) of aqueous species for ionic strength calculations
        type(aq_phase_c) :: aq_phase !< Aqueous phase object (assumption: only 1 aqueous phase)
        !integer(kind=4) :: num_gas_phases=0 !< Number of gas phases
        type(gas_phase_c) :: gas_phase !< Gas phase object (assumption: only 1 gas phase)
        integer(kind=4) :: num_var_act_species !< Number of variable activity species
        integer(kind=4), allocatable :: var_act_sp_indices(:) !< Indices of variable activity species in "species" array
        integer(kind=4) :: num_minerals=0 !< Total number of minerals
        integer(kind=4) :: num_minerals_eq=0 !< Number of minerals in equilibrium
        integer(kind=4) :: num_minerals_eq_cst_act=0 !< Number of equilibrium minerals with constant activity
        integer(kind=4) :: num_minerals_eq_var_act=0 !< Number of equilibrium minerals with variable activity
        integer(kind=4) :: num_minerals_kin=0 !< Number of minerals NOT in equilibrium (kinetic)
        integer(kind=4) :: num_minerals_kin_cst_act=0 !< Number of kinetic minerals with constant activity
        integer(kind=4) :: num_minerals_kin_var_act=0 !< Number of kinetic minerals with variable activity
        type(mineral_c), allocatable :: minerals(:) !< Array of minerals (order: kin var act, kin cst act, eq var act, eq cst act)
        integer(kind=4) :: num_cst_act_species !< Number of constant activity species
        integer(kind=4), allocatable :: cst_act_sp_indices(:) !< Indices of constant activity species in "species" array
        integer(kind=4) :: num_solids=0 !< Total number of solid phases (minerals + surface complexes)
        !integer(kind=4) :: num_surf_compl=0 !< Number of surface complexes
        type(cat_exch_zone_c) :: cat_exch_zone !< Surface complex objects for surface complexation
        real(kind=8), allocatable :: stoich_mat(:,:) !< Global stoichiometric matrix (S)
        real(kind=8), allocatable :: stoich_mat_sol(:,:) !< Solid stoichiometric matrix (S_s)
        real(kind=8), allocatable :: stoich_mat_gas(:,:) !< Gas stoichiometric matrix (S_g)
        real(kind=8), allocatable :: Se(:,:) !< Equilibrium stoichiometric matrix
        real(kind=8), allocatable :: Sk(:,:) !< Kinetic stoichiometric matrix
        integer(kind=4) :: num_eq_reacts_init=0 !< Total number of equilibrium reactions at initialization
        integer(kind=4) :: num_aq_eq_reacts=0 !< Number of aqueous equilibrium reactions
        integer(kind=4) :: num_aq_kin_reacts=0 !< Number of aqueous kinetic reactions
        integer(kind=4) :: num_redox_eq_reacts=0 !< Number of redox aqueous equilibrium reactions
        type(eq_reaction_c), allocatable :: eq_reacts(:) !< Array of equilibrium reactions
        type(speciation_algebra_s) :: speciation_alg !< Speciation algebra object for solving chemical equilibria
        integer(kind=4) :: num_kin_reacts=0 !< Total number of kinetic reactions
        type(kin_reaction_poly_c), allocatable :: kin_reacts(:) !< Array of kinetic reactions (polymorphic pointer array)
        integer(kind=4) :: num_lin_kin_reacts=0 !< Number of linear kinetic reactions
        type(lin_kin_reaction_c), allocatable :: lin_kin_reacts(:) !< Array of linear kinetic reactions (assumed aqueous)
        type(kin_mineral_c), allocatable :: min_kin_reacts(:) !< Array of mineral kinetic reactions
        integer(kind=4) :: num_redox_kin_reacts=0 !< Number of Monod (redox) kinetic reactions
        type(redox_kin_c), allocatable :: redox_kin_reacts(:) !< Array of redox kinetic reactions
        integer(kind=4) :: num_gas_kin_reacts=0 !< Number of gas kinetic reactions
    contains
    !> @name Set Methods
    !> @brief Setter procedures for chemical system attributes
    !> @{
        !procedure :: set_num_species !< Set the number of species
        procedure :: set_num_minerals !< Set the number of minerals
        procedure :: set_num_minerals_eq !< Set the number of equilibrium minerals
        procedure :: set_num_minerals_eq_cst_act !< Set the number of constant activity equilibrium minerals
        procedure :: set_num_minerals_eq_var_act !< Set the number of variable activity equilibrium minerals
        procedure :: set_num_minerals_kin_cst_act !< Set the number of constant activity kinetic minerals
        procedure :: set_num_minerals_kin_var_act !< Set the number of variable activity kinetic minerals
        procedure :: set_num_kin_reacts !< Set the number of kinetic reactions
        procedure :: set_num_gas_kin_reacts !< Set the number of gas kinetic reactions
        procedure :: set_num_lin_kin_reacts !< Set the number of linear kinetic reactions
        procedure :: set_num_minerals_kin !< Set the number of kinetic minerals
        procedure :: set_num_redox_kin_reacts !< Set the number of redox kinetic reactions
        procedure :: set_num_redox_eq_reacts !< Set the number of redox equilibrium reactions
        procedure :: set_num_aq_eq_reacts !< Set the number of aqueous equilibrium reactions
        procedure :: set_num_aq_kin_reacts !< Set the number of aqueous kinetic reactions
        procedure :: set_num_cst_act_species !< Set the number of constant activity species
        procedure :: set_num_var_act_species !< Set the number of variable activity species
        procedure :: set_species !< Set the species array
        procedure :: set_cat_exch_zone !< Set the cation exchange zones
        procedure :: set_eq_reacts !< Set the equilibrium reactions array
        procedure :: set_kin_reacts !< Set the kinetic reactions array
        procedure :: set_stoich_mat !< Set the global stoichiometric matrix
        procedure :: set_stoich_mat_gas !< Set the gas stoichiometric matrix
        procedure :: set_stoich_mat_sol !< Set the solid stoichiometric matrix
    !> @}
    !> @name Allocate Methods
    !> @brief Allocation procedures for chemical system arrays
    !> @{
        procedure :: allocate_species !< Allocate the species array
        procedure :: allocate_cst_act_sp_indices !< Allocate the constant activity species indices array
        procedure :: allocate_var_act_sp_indices !< Allocate the variable activity species indices array
        procedure :: allocate_reacts !< Allocate both equilibrium and kinetic reaction arrays
        procedure :: allocate_eq_reacts !< Allocate the equilibrium reactions array
        procedure :: allocate_kin_reacts !< Allocate the kinetic reactions array
        procedure :: allocate_redox_kin_reacts !< Allocate the redox kinetic reactions array
        procedure :: allocate_minerals !< Allocate the minerals array
        procedure :: allocate_min_kin_reacts !< Allocate the mineral kinetic reactions array
        procedure :: allocate_lin_kin_reacts !< Allocate the linear kinetic reactions array
        !procedure :: allocate_cat_exch_zones !< Allocate the linear kinetic reactions array
    !>@}
    !> @name Compute Methods
    !> @brief Computation procedures for derived quantities
    !> @{
        procedure :: compute_num_kin_reacts !< Compute total number of kinetic reactions
        procedure :: compute_num_reacts !< Compute total number of reactions
        procedure :: compute_num_species !< Compute total number of species
        procedure :: compute_z2 !< Compute squared charges of species
        procedure :: compute_num_solids_chem_syst !< Compute total number of solids
    !> @}
    !> @name Read Methods
    !> @brief Input/reading procedures for chemical system data
    !> @{
        procedure :: read_chem_system_CHEPROO !< Read chemical system from CHEPROO format
        procedure :: read_chem_system_PFLOTRAN !< Read chemical system from PFLOTRAN format
        procedure :: read_master25 !< Read chemical system from master25 database
        procedure :: read_kinetics_DB !< Read kinetics database
        procedure :: read_Monod_DB !< Read Monod (redox) reactions
        procedure :: read_PHREEQC_DB_opc1 !< Read PHREEQC database option 1
        procedure :: read_PHREEQC_DB_opc2 !< Read PHREEQC database option 2
    !> @}
    !> @name Query Methods
    !> @brief Query procedures to check system contents
    !> @{
        procedure :: is_mineral_in_chem_syst !< Check if a mineral exists in the system
        procedure :: is_reaction_in_chem_syst !< Check if a reaction exists in the system
        procedure :: is_eq_reaction_in_chem_syst !< Check if an equilibrium reaction exists in the system
        procedure :: is_kin_reaction_in_chem_syst !< Check if a kinetic reaction exists in the system
        procedure :: is_species_in_chem_syst !< Check if a species exists in the system
    !> @}
    !> @name Get Methods
    !> @brief Getter procedures for chemical system attributes
    !> @{
        procedure :: get_eq_csts !< Get equilibrium constants of all equilibrium reactions
        procedure :: get_mineral_index_by_name !< Get mineral index by name
    !> @}
    !> @name Rearrange Methods
    !> @brief Procedures to rearrange arrays in specific orders
    !> @{
        procedure :: rearrange_eq_reacts !< Rearrange equilibrium reactions in canonical order
        procedure :: rearrange_species !< Rearrange species based on speciation algebra
    !> @}
    !> @name Utility Methods
    !> @brief Utility and helper procedures
    !> @{
        procedure :: compute_eq_csts_gases_cst_act !< Compute equilibrium constants for constant activity gases
        procedure :: change_spec_alg_chem_syst !< Change speciation algebra settings
    !> @}
    end type
!> @brief Interface declarations for external subroutines
!> @details This interface block declares the signatures of subroutines implemented in separate files
    interface
        !> @brief Read generic chemical system from file
        !> @param[in,out] this Chemical system object
        !> @param[in] path Directory path  
        !> @param[in] filename Name of the file
        subroutine read_chem_system(this,path,filename)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this
            character(len=*), intent(in) :: path
            character(len=*), intent(in) :: filename
        end subroutine
        
        !> @brief Read chemical system from CHEPROO format database
        !> @param[in,out] this Chemical system object
        !> @param[in] path_DB Path to CHEPROO database directory
        !> @param[in] unit File unit number for reading
        subroutine read_chem_system_CHEPROO(this,path_DB,unit)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            character(len=*), intent(in) :: path_DB !< Path to database
            integer(kind=4), intent(in) :: unit !< File unit number
        end subroutine
        
        !> @brief Read chemical system from PFLOTRAN format database
        !> @param[in,out] this Chemical system object
        !> @param[in] path Path to PFLOTRAN database directory
        !> @param[in] unit File unit number for reading
        subroutine read_chem_system_PFLOTRAN(this,path,unit)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            character(len=*), intent(in) :: path !< Path to database
            integer(kind=4), intent(in) :: unit !< File unit number
        end subroutine

        !> @brief Read chemical system from master25 thermodynamic database
        !> @param[in,out] this Chemical system object
        !> @param[in] path Path to master25 database file
        !> @param[in] unit File unit number for reading
        subroutine read_master25(this,path,unit)
            import chem_system_c
            !import reaction_c
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            character(len=*), intent(in) :: path !< Path to database file
            integer(kind=4), intent(in) :: unit !< File unit number
        end subroutine
        
        !> @brief Read kinetics database containing kinetic reaction data
        !> @param[in,out] this Chemical system object
        !> @param[in] path Path to kinetics database file
        !> @param[in] unit File unit number for reading
        subroutine read_kinetics_DB(this,path,unit)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            character(len=*), intent(in) :: path !< Path to database file
            integer(kind=4), intent(in) :: unit !< File unit number
        end subroutine
        
        !> @brief Read Monod (microbial redox) reactions from file
        !> @param[in,out] this Chemical system object
        !> @param[in] path Path to Monod reactions file
        !> @param[in] unit File unit number for reading
        subroutine read_Monod_DB(this,path,unit)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            character(len=*), intent(in) :: path !< Path to Monod file
            integer(kind=4), intent(in) :: unit !< File unit number
        end subroutine
        
        !> @brief Read PHREEQC database using option 1 format
        !> @param[in,out] this Chemical system object
        !> @param[in] filename Name of PHREEQC database file
        subroutine read_PHREEQC_DB_opc1(this,filename)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            character(len=*), intent(in) :: filename !< Database filename
        end subroutine
        
        !> @brief Read PHREEQC database using option 2 format
        !> @param[in,out] this Chemical system object
        !> @param[in] filename Name of PHREEQC database file
        subroutine read_PHREEQC_DB_opc2(this,filename)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            character(len=*), intent(in) :: filename !< Database filename
        end subroutine
        
        !> @brief Read PFLOTRAN database file
        !> @param[in,out] this Chemical system object
        !> @param[in] unit File unit number for reading
        !> @param[in] filename Name of PFLOTRAN database file
        subroutine read_PFLOTRAN_DB(this,unit,filename)
        import chem_system_c
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in) :: unit !< File unit number
            character(len=*), intent(in) :: filename !< Database filename
        end subroutine
        

        
      
       
     
       
  
    
        !> @brief Set global stoichiometric matrix for all reactions
        !> @param[in,out] this Chemical system object
        subroutine set_stoich_mat(this)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this !< Chemical system object
        end subroutine
        
        !> @brief Set gas phase stoichiometric matrix
        !> @param[in,out] this Chemical system object
        subroutine set_stoich_mat_gas(this)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this !< Chemical system object
        end subroutine
        
        !> @brief Set solid phase stoichiometric matrix
        !> @param[in,out] this Chemical system object
        subroutine set_stoich_mat_sol(this)
            import chem_system_c
            implicit none
            class(chem_system_c) :: this !< Chemical system object
        end subroutine
!> @brief Marker for end of interface declarations
    end interface
!> @brief Separator line for implementation section  
!***************************************************************************************************************************************************!
    contains
!> @brief Separator line for SET methods section
!*********************** SET ***********************************************************************************************************************!
       
        
        !> @brief Set the number of reactions in the chemical system
        !> @details If num_reacts is provided, validates and sets it. Otherwise computes from equilibrium and kinetic reactions.
        !> @param[in,out] this Chemical system object
        !> @param[in] num_reacts Optional number of reactions
        subroutine set_num_reacts(this,num_reacts)
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in), optional :: num_reacts !< Optional number of reactions to set
            if (present(num_reacts)) then !> Check if argument is provided
                if (num_reacts>this%speciation_alg%num_species) then !> Validate that reactions don't exceed species count
                    error stop "Number of reactions cannot be greater than number of species" !> Error if validation fails
                else if (num_reacts<0) then !> Validate that number is non-negative
                    error stop "Number of reactions cannot be negative" !> Error if negative
                else
                    this%num_reacts=num_reacts !> Set the attribute value
                end if
            else !> If argument not provided
                this%num_reacts=this%speciation_alg%num_eq_reactions+this%num_kin_reacts !> Compute as sum of equilibrium and kinetic reactions
            end if
        end subroutine
        
      

        !> @brief Set the number of kinetic reactions
        !> @details Validates that kinetic reactions don't exceed total reactions
        !> @param[in,out] this Chemical system object
        !> @param[in] num_kin_reacts Number of kinetic reactions
        subroutine set_num_kin_reacts(this,num_kin_reacts)
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in) :: num_kin_reacts !< Number of kinetic reactions
            if (num_kin_reacts>this%num_reacts) error stop "Number of kinetic reactions cannot be greater than number of reactions" !> Validate and error if exceeded
            this%num_kin_reacts=num_kin_reacts !> Assign value to attribute
        end subroutine

        !> @brief Set the number of gas kinetic reactions
        !> @details Validates that gas kinetic reactions don't exceed total kinetic reactions
        !> @param[in,out] this Chemical system object
        !> @param[in] num_gas_kin_reacts Number of gas kinetic reactions
        subroutine set_num_gas_kin_reacts(this,num_gas_kin_reacts)
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in) :: num_gas_kin_reacts !< Number of gas kinetic reactions
            if (num_gas_kin_reacts>this%num_kin_reacts) error stop "Number of gas kinetic reactions cannot be greater than number & 
                of kinetic reactions" !> Validate against total kinetic reactions
            this%num_gas_kin_reacts=num_gas_kin_reacts !> Assign value to attribute
        end subroutine
        
        !> @brief Set the species array and update num_species
        !> @details Validates array size against current num_species if already allocated
        !> @param[in,out] this Chemical system object
        !> @param[in] species Array of species to assign
        subroutine set_species(this,species)
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            class(species_c), intent(in) :: species(:) !< Array of species to set
                        
            if (allocated(this%species) .and. size(species)/=this%speciation_alg%num_species) then !> Check size consistency
                error stop "Wrong number of species" !> Error if sizes don't match
            else
                this%species=species !> Assign species array
                this%speciation_alg%num_species=size(species) !> Update number of species
            end if
        end subroutine
        
        !> @brief Set the number of constant activity species
        !> @details If argument provided, uses it. Otherwise computes from indices array size.
        !> @param[in,out] this Chemical system object
        !> @param[in] num_cst_act_species Optional number of constant activity species
        subroutine set_num_cst_act_species(this,num_cst_act_species)
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in), optional :: num_cst_act_species !< Optional number to set
            if (present(num_cst_act_species)) then !> If argument provided
                this%num_cst_act_species=num_cst_act_species !> Use provided value
            else !> If not provided
                this%num_cst_act_species=size(this%cst_act_sp_indices) !> Compute from array size
            end if
        end subroutine
        
        
        !> @brief Set the cation exchange object
        !> @param[in,out] this Chemical system object
        !> @param[in] cat_exch Cation exchange object to assign
        subroutine set_cat_exch_zone(this,cat_exch_zone)
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            class(cat_exch_zone_c), intent(in) :: cat_exch_zone !< Cation exchange object
            this%cat_exch_zone=cat_exch_zone !> Assign cation exchange object
        end subroutine
                
        !> @brief Set the number of variable activity species
        !> @details Computes from optional argument, indices array, or difference between total and constant activity species
        !> @param[in,out] this Chemical system object
        !> @param[in] num_var_act_species Optional number of variable activity species
        subroutine set_num_var_act_species(this,num_var_act_species)
        !> This subroutine sets the "num_var_act_species" attribute 
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in), optional :: num_var_act_species !< Optional number to set
            if (present(num_var_act_species)) then !> If argument provided
                this%num_var_act_species=num_var_act_species !> Use provided value
            else if (allocated(this%var_act_sp_indices)) then !> If indices array allocated
                this%num_var_act_species=size(this%var_act_sp_indices) !> Compute from array size
            else if (allocated(this%species) .and. allocated(this%cst_act_sp_indices)) then !> If species and cst act indices allocated
                this%num_var_act_species=this%speciation_alg%num_species-this%num_cst_act_species !> Compute as difference
            else !> If unable to determine
                error stop "Unable to compute the number of variable activity species" !> Error message
            end if
        end subroutine
        
       
        
        !> @brief Set the total number of minerals
        !> @details Validates that the number is non-negative
        !> @param[in,out] this Chemical system object
        !> @param[in] num_minerals Number of minerals to set
        subroutine set_num_minerals(this,num_minerals)
        !> This subroutine sets the "num_minerals"
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in) :: num_minerals !< Number of minerals
            if (num_minerals<0) error stop "Number of minerals cannot be negative" !> Validate non-negative
            this%num_minerals=num_minerals !> Assign value
        end subroutine
        
        !> @brief Set the number of equilibrium minerals
        !> @details Validates against total number of minerals
        !> @param[in,out] this Chemical system object
        !> @param[in] num_minerals_eq Number of equilibrium minerals
        subroutine set_num_minerals_eq(this,num_minerals_eq)
        !> This subroutine sets the "num_minerals_eq" attribute
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in) :: num_minerals_eq !< Number of equilibrium minerals
            if (num_minerals_eq>this%num_minerals .AND. allocated(this%minerals)) error stop "Number of minerals in equilibrium &
                cannot be greater than number of minerals" !> Validate against total minerals
            this%num_minerals_eq=num_minerals_eq !> Assign value
        end subroutine
        
        !> @brief Set the number of constant activity equilibrium minerals
        !> @details Validates against total equilibrium minerals
        !> @param[in,out] this Chemical system object
        !> @param[in] num_minerals_eq_cst_act Number of constant activity equilibrium minerals
        subroutine set_num_minerals_eq_cst_act(this,num_minerals_eq_cst_act)
        !> This subroutine sets the "num_minerals_eq_cst_act" attribute
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in) :: num_minerals_eq_cst_act !< Number of cst act eq minerals
            if (num_minerals_eq_cst_act>this%num_minerals_eq .AND. allocated(this%minerals)) then !> Check bound
                error stop "Number of minerals in equilibrium with constant activity cannot be greater than number of minerals in & 
                    equilibrium" !> Error if exceeds eq minerals
            end if
            this%num_minerals_eq_cst_act=num_minerals_eq_cst_act !> Assign value
        end subroutine
        
        !> @brief Set the number of variable activity equilibrium minerals
        !> @details Validates against total equilibrium minerals
        !> @param[in,out] this Chemical system object
        !> @param[in] num_minerals_eq_var_act Number of variable activity equilibrium minerals
        subroutine set_num_minerals_eq_var_act(this,num_minerals_eq_var_act)
        !> This subroutine sets the "num_minerals_eq_var_act" attribute
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in) :: num_minerals_eq_var_act !< Number of var act eq minerals
            if (num_minerals_eq_var_act>this%num_minerals_eq .AND. allocated(this%minerals)) then !> Check bound
                error stop "Number of minerals in equilibrium with constant activity cannot be greater than number of minerals in & 
                    equilibrium" !> Error if exceeds eq minerals
            end if
            this%num_minerals_eq_var_act=num_minerals_eq_var_act !> Assign value
        end subroutine
        
        !> @brief Set the number of constant activity kinetic minerals
        !> @details Validates against total kinetic minerals
        !> @param[in,out] this Chemical system object
        !> @param[in] num_minerals_kin_cst_act Number of constant activity kinetic minerals
        subroutine set_num_minerals_kin_cst_act(this,num_minerals_kin_cst_act)
        !> This subroutine sets the "num_minerals_kin_cst_act" attribute
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in) :: num_minerals_kin_cst_act !< Number of cst act kin minerals
            if (num_minerals_kin_cst_act>this%num_minerals_kin .AND. allocated(this%minerals)) then !> Check bound
                error stop "Number of minerals NOT in equilibrium with constant activity cannot be greater than number of minerals &
                    NOT in equilibrium" !> Error if exceeds kin minerals
            end if
            this%num_minerals_kin_cst_act=num_minerals_kin_cst_act !> Assign value
        end subroutine
        
        !> @brief Set the number of variable activity kinetic minerals
        !> @details Validates against total kinetic minerals
        !> @param[in,out] this Chemical system object
        !> @param[in] num_minerals_kin_var_act Number of variable activity kinetic minerals
        subroutine set_num_minerals_kin_var_act(this,num_minerals_kin_var_act)
        !> This subroutine sets the "num_minerals_kin_var_act" attribute
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in) :: num_minerals_kin_var_act !< Number of var act kin minerals
            if (num_minerals_kin_var_act>this%num_minerals_kin .AND. allocated(this%minerals)) then !> Check bound
                error stop "Number of minerals NOT in equilibrium with constant activity cannot be greater than number of minerals &
                    NOT in equilibrium" !> Error if exceeds kin minerals
            end if
            this%num_minerals_kin_var_act=num_minerals_kin_var_act !> Assign value
        end subroutine
        
        !> @brief Set the equilibrium reactions array
        !> @details Deallocates existing array if present, assigns new array, updates count
        !> @param[in,out] this Chemical system object
        !> @param[in] eq_reacts Array of equilibrium reactions
        subroutine set_eq_reacts(this,eq_reacts)
        !> This subroutine sets the "eq_reacts" and "num_eq_reacts" attributes 
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            class(eq_reaction_c), intent(in) :: eq_reacts(:) !< Array of equilibrium reactions
            if (allocated(this%eq_reacts)) then !> If array already allocated
                deallocate(this%eq_reacts) !> Deallocate existing array
            end if
            this%eq_reacts=eq_reacts !> Assign new array
            this%speciation_alg%num_eq_reactions=size(this%eq_reacts) !> Update count from array size
        end subroutine
        
        !> @brief Set the kinetic reactions array
        !> @details Deallocates existing array if present, assigns new array, updates count
        !> @param[in,out] this Chemical system object
        !> @param[in] kin_reacts Array of kinetic reactions (polymorphic)
        subroutine set_kin_reacts(this,kin_reacts)
        !> This subroutine sets the "kin_reacts" and "num_kin_reacts" attributes 
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            type(kin_reaction_poly_c), intent(in) :: kin_reacts(:) !< Array of kinetic reactions
            if (allocated(this%kin_reacts)) then !> If array already allocated
                deallocate(this%kin_reacts) !> Deallocate existing array
            end if
            this%kin_reacts=kin_reacts !> Assign new array
            this%num_kin_reacts=size(this%kin_reacts) !> Update count from array size
        end subroutine
        
       !> @brief Set the number of linear kinetic reactions
       !> @param[in,out] this Chemical system object
       !> @param[in] num_lin_kin_reacts Number of linear kinetic reactions
       subroutine set_num_lin_kin_reacts(this,num_lin_kin_reacts)
        !> This subroutine sets the "num_lin_kin_reacts" attribute 
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in) :: num_lin_kin_reacts !< Number of linear kinetic reactions
            this%num_lin_kin_reacts=num_lin_kin_reacts !> Assign value
       end subroutine
       
       !> @brief Set the number of kinetic minerals (NOT in equilibrium)
       !> @param[in,out] this Chemical system object
       !> @param[in] num_minerals_kin Number of kinetic minerals
       subroutine set_num_minerals_kin(this,num_minerals_kin)
        !> This subroutine sets the "num_minerals_kin" attribute
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in) :: num_minerals_kin !< Number of kinetic minerals
            this%num_minerals_kin=num_minerals_kin !> Assign value
       end subroutine
       
       !> @brief Set the number of redox kinetic reactions (Monod reactions)
       !> @param[in,out] this Chemical system object
       !> @param[in] num_redox_kin_reacts Number of redox kinetic reactions
       subroutine set_num_redox_kin_reacts(this,num_redox_kin_reacts)
        !> This subroutine sets the "num_redox_kin_reacts" attribute
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in) :: num_redox_kin_reacts !< Number of redox kinetic reactions
            this%num_redox_kin_reacts=num_redox_kin_reacts !> Assign value
       end subroutine
       
       !> @brief Set the number of redox equilibrium reactions
       !> @param[in,out] this Chemical system object
       !> @param[in] num_redox_eq_reacts Number of redox equilibrium reactions
       subroutine set_num_redox_eq_reacts(this,num_redox_eq_reacts)
        !> This subroutine sets the "num_redox_eq_reacts" attribute
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in) :: num_redox_eq_reacts !< Number of redox equilibrium reactions
            this%num_redox_eq_reacts=num_redox_eq_reacts !> Assign value
       end subroutine
       
       !> @brief Set the number of aqueous equilibrium reactions
       !> @details If not provided, computes as redox eq reactions plus aqueous complexes
       !> @param[in,out] this Chemical system object
       !> @param[in] num_aq_eq_reacts Optional number of aqueous equilibrium reactions
       subroutine set_num_aq_eq_reacts(this,num_aq_eq_reacts)
       !> This subroutine sets the "num_aq_eq_reacts" attribute
           class(chem_system_c) :: this !< Chemical system object
           integer(kind=4), intent(in), optional :: num_aq_eq_reacts !< Optional number to set
           if (present(num_aq_eq_reacts)) then !> If argument provided
               this%num_aq_eq_reacts=num_aq_eq_reacts !> Use provided value
           else !> If not provided
               this%num_aq_eq_reacts=this%num_redox_eq_reacts+this%aq_phase%num_aq_complexes !> Compute as sum
           end if
       end subroutine
       
       !> @brief Set the number of aqueous kinetic reactions
       !> @details If not provided, computes as redox kin reactions plus linear kin reactions
       !> @param[in,out] this Chemical system object
       !> @param[in] num_aq_kin_reacts Optional number of aqueous kinetic reactions
       subroutine set_num_aq_kin_reacts(this,num_aq_kin_reacts)
       !> This subroutine sets the "num_aq_kin_reacts" attribute
           class(chem_system_c) :: this !< Chemical system object
           integer(kind=4), intent(in), optional :: num_aq_kin_reacts !< Optional number to set
           if (present(num_aq_kin_reacts)) then !> If argument provided
               this%num_aq_kin_reacts=num_aq_kin_reacts !> Use provided value
           else !> If not provided
               this%num_aq_kin_reacts=this%num_redox_kin_reacts+this%num_lin_kin_reacts !> Compute as sum
           end if
       end subroutine
       
!> @brief Separator line for ALLOCATE methods section
!*********************** ALLOCATE ******************************************************************************************************************!
        !> @brief Allocate the species array
        !> @details Sets num_species if provided or computes it, deallocates existing array, allocates new array
        !> @param[in,out] this Chemical system object
        !> @param[in] num_species Optional number of species
        subroutine allocate_species(this,num_species)
        !< This subroutine allocates the attribute "species"
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in), optional :: num_species !< Optional number to allocate
            if (present(num_species)) then !> If argument provided
                if (num_species<0) then !> Validate non-negative
                    error stop "Number of species cannot be negative" !> Error if negative
                else
                    this%speciation_alg%num_species=num_species !> Set attribute
                end if
            else !> If not provided
                call this%compute_num_species() !> Compute from other attributes
            end if
            if (allocated(this%species)) then !> If already allocated
                deallocate(this%species) !> Deallocate existing array
            end if
            allocate(this%species(this%speciation_alg%num_species)) !> Allocate array with correct size
        end subroutine
        
        !> @brief Allocate the constant activity species indices array
        !> @details Sets num_cst_act_species if provided, deallocates existing, allocates new
        !> @param[in,out] this Chemical system object
        !> @param[in] num_cst_act_species Optional number of constant activity species
        subroutine allocate_cst_act_sp_indices(this,num_cst_act_species)
        !< This subroutine allocates the attribute "cst_act_sp_indices"
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in), optional :: num_cst_act_species !< Optional number to allocate
            if (present(num_cst_act_species)) then !> If argument provided
                this%num_cst_act_species=num_cst_act_species !> Set attribute
            end if
            if (allocated(this%cst_act_sp_indices)) then !> If already allocated
                deallocate(this%cst_act_sp_indices) !> Deallocate existing array
            end if
            allocate(this%cst_act_sp_indices(this%num_cst_act_species)) !> Allocate with correct size
        end subroutine
        
        !> @brief Allocate the variable activity species indices array
        !> @details Sets num_var_act_species if provided, allocates array
        !> @param[in,out] this Chemical system object
        !> @param[in] num_var_act_species Optional number of variable activity species
        subroutine allocate_var_act_sp_indices(this,num_var_act_species)
        !< This subroutine allocates the attribute "var_act_sp_indices"
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in), optional :: num_var_act_species !< Optional number to allocate
            if (present(num_var_act_species)) then !> If argument provided
                this%num_var_act_species=num_var_act_species !> Set attribute
            end if
            allocate(this%var_act_sp_indices(this%num_var_act_species)) !> Allocate with correct size
        end subroutine
        
        !> @brief Allocate both equilibrium and kinetic reaction arrays
        !> @details Calls individual allocate methods for each reaction type
        !> @param[in,out] this Chemical system object
        !> @param[in] num_eq_reacts Number of equilibrium reactions
        !> @param[in] num_kin_reacts Number of kinetic reactions
        subroutine allocate_reacts(this,num_eq_reacts,num_kin_reacts)
        !< This subroutine allocates the attributes "eq_reacts" & "kin_reacts"
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in) :: num_eq_reacts,num_kin_reacts !< Reaction counts
            call this%allocate_eq_reacts(num_eq_reacts) !> Allocate equilibrium reactions
            call this%allocate_kin_reacts(num_kin_reacts) !> Allocate kinetic reactions
            !call this%compute_num_reacts() !> (Commented) Compute total reactions
        end subroutine
        
        !> @brief Allocate the equilibrium reactions array
        !> @details Sets num_eq_reacts if provided, allocates array
        !> @param[in,out] this Chemical system object
        !> @param[in] num_eq_reacts Optional number of equilibrium reactions
        subroutine allocate_eq_reacts(this,num_eq_reacts)
        !< This subroutine allocates the attribute "eq_reacts"
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in), optional :: num_eq_reacts !< Optional number to allocate
            if (present(num_eq_reacts)) then !> If argument provided
                this%speciation_alg%num_eq_reactions=num_eq_reacts !> Set attribute
            end if
            allocate(this%eq_reacts(this%speciation_alg%num_eq_reactions)) !> Allocate with correct size
        end subroutine
        
        !> @brief Allocate the kinetic reactions array
        !> @details Sets num_kin_reacts if provided, allocates array
        !> @param[in,out] this Chemical system object
        !> @param[in] num_kin_reacts Optional number of kinetic reactions
        subroutine allocate_kin_reacts(this,num_kin_reacts)
        !< This subroutine allocates the attribute "kin_reacts"
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in), optional :: num_kin_reacts !< Optional number to allocate
            if (present(num_kin_reacts)) then !> If argument provided
                this%num_kin_reacts=num_kin_reacts !> Set attribute
            end if
            allocate(this%kin_reacts(this%num_kin_reacts)) !> Allocate with correct size
        end subroutine
        
        !> @brief Allocate the redox kinetic reactions array
        !> @details Sets num_redox_kin_reacts if provided, allocates array
        !> @param[in,out] this Chemical system object
        !> @param[in] num_redox_kin_reacts Optional number of redox kinetic reactions
        subroutine allocate_redox_kin_reacts(this,num_redox_kin_reacts)
        !< This subroutine allocates the attribute "redox_kin_reacts"
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in), optional :: num_redox_kin_reacts !< Optional number to allocate
            if (present(num_redox_kin_reacts)) then !> If argument provided
                this%num_redox_kin_reacts=num_redox_kin_reacts !> Set attribute
            end if
            allocate(this%redox_kin_reacts(this%num_redox_kin_reacts)) !> Allocate with correct size
        end subroutine
        
        !> @brief Allocate the linear kinetic reactions array
        !> @details Sets num_lin_kin_reacts if provided, allocates array
        !> @param[in,out] this Chemical system object
        !> @param[in] num_lin_kin_reacts Optional number of linear kinetic reactions
        subroutine allocate_lin_kin_reacts(this,num_lin_kin_reacts)
        !< This subroutine allocates the attribute "lin_kin_reacts"
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in), optional :: num_lin_kin_reacts !< Optional number to allocate
            if (present(num_lin_kin_reacts)) then !> If argument provided
                this%num_lin_kin_reacts=num_lin_kin_reacts !> Set attribute
            end if
            allocate(this%lin_kin_reacts(this%num_lin_kin_reacts)) !> Allocate with correct size
        end subroutine
        
        !> @brief Allocate the minerals array
        !> @details Sets num_minerals if provided, allocates array
        !> @param[in,out] this Chemical system object
        !> @param[in] num_minerals Optional number of minerals
        subroutine allocate_minerals(this,num_minerals)
        !< This subroutine allocates the attribute "minerals"
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in), optional :: num_minerals !< Optional number to allocate
            if (present(num_minerals)) then !> If argument provided
                call this%set_num_minerals(num_minerals) !> Set using setter method
            end if
            allocate(this%minerals(this%num_minerals)) !> Allocate with correct size
        end subroutine
        
        !> @brief Allocate the mineral kinetic reactions array
        !> @details Requires minerals array to be allocated first, sets num_minerals_kin if provided
        !> @param[in,out] this Chemical system object
        !> @param[in] num_min_kin Optional number of mineral kinetic reactions
        subroutine allocate_min_kin_reacts(this,num_min_kin)
        !< This subroutine allocates the attribute "min_kin_reacts"
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            integer(kind=4), intent(in), optional :: num_min_kin !< Optional number to allocate
            if (.not. allocated(this%minerals)) error stop !> Error if minerals not allocated
            if (present(num_min_kin)) then !> If argument provided
                this%num_minerals_kin=num_min_kin !> Set attribute
            end if
            allocate(this%min_kin_reacts(this%num_minerals_kin)) !> Allocate with correct size
        end subroutine        

!> @brief Separator line for COMPUTE methods section
!*********************** COMPUTE ***********************************************************************************************************************!
        !> @brief Compute total number of species
        !> @details Computes as sum of variable activity and constant activity species
        !> @param[in,out] this Chemical system object
        subroutine compute_num_species(this)
        !< This subroutine computes the attribute "num_species"
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            this%speciation_alg%num_species=this%num_var_act_species+this%num_cst_act_species !> Sum var act and cst act species
        end subroutine
        
        !> @brief Compute total number of kinetic reactions
        !> @details Computes as sum of linear kin, mineral kin, and redox kin reactions
        !> @param[in,out] this Chemical system object
        subroutine compute_num_kin_reacts(this)
        !< This subroutine computes the attribute "num_kin_reacts"
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            this%num_kin_reacts=this%num_lin_kin_reacts+this%num_minerals_kin+this%num_redox_kin_reacts !> Sum all kinetic reaction types
        end subroutine
        
        
        
        
       
        
        
!> @brief Separator line for IS (query) methods section
!********************************** IS *************************************************************************************************************!
        !> @brief Check if a mineral belongs to the chemical system
        !> @details Searches minerals array by name, returns flag and optional index
        !> @param[in] this Chemical system object
        !> @param[in] mineral Mineral to search for
        !> @param[out] flag TRUE if mineral found, FALSE otherwise
        !> @param[out] mineral_ind Optional index in minerals array (0 if not found)
        subroutine is_mineral_in_chem_syst(this,mineral,flag,mineral_ind)
        !> This subroutine checks if a mineral belongs to the chemical system
            implicit none
            class(chem_system_c), intent(in) :: this !< chemical system
            class(mineral_c), intent(in) :: mineral !< mineral
            logical, intent(out) :: flag !< TRUE if mineral belongs to chemical system, FALSE otherwise
            integer(kind=4), intent(out), optional :: mineral_ind !> index of mineral in "minerals" attribute (if not belongs: 0)
            
            integer(kind=4) :: i !> Loop counter
            
            flag=.false. !> Initialize to not found
            if (present(mineral_ind)) then !> If index requested
                mineral_ind=0 !> Initialize to zero (not found)
            end if
            do i=1,this%num_minerals !> Loop through all minerals
                if (mineral%name==this%minerals(i)%name) then !> Compare names
                    flag=.true. !> Set found flag
                    if (present(mineral_ind)) then !> If index requested
                        mineral_ind=i !> Return index
                    end if
                    exit !> Exit loop early (found)
                end if
            end do
        end subroutine
        
        !> @brief Check if a species belongs to the chemical system
        !> @details Searches species array by name, returns flag and optional index
        !> @param[in] this Chemical system object
        !> @param[in] species Species to search for
        !> @param[out] flag TRUE if species found, FALSE otherwise
        !> @param[out] species_ind Optional index in species array (0 if not found)
        subroutine is_species_in_chem_syst(this,species,flag,species_ind)
        !> This subroutine checks if a species belongs to the chemical system
            implicit none
            class(chem_system_c), intent(in) :: this !< chemical system
            class(species_c), intent(in) :: species !< species
            logical, intent(out) :: flag !> TRUE if species belongs to chemical system, FALSE otherwise
            integer(kind=4), intent(out), optional :: species_ind !> species index in "species" attribute (if not belongs: 0)
            
            integer(kind=4) :: i !> Loop counter
            
            flag=.false. !> Initialize to not found
            if (present(species_ind)) then !> If index requested
                species_ind=0 !> Initialize to zero (not found)
            end if
            do i=1,this%speciation_alg%num_species !> Loop through all species
                if (species%name==this%species(i)%name) then !> Compare names
                    flag=.true. !> Set found flag
                    if (present(species_ind)) then !> If index requested
                        species_ind=i !> Return index
                    end if
                    exit !> Exit loop early (found)
                end if
            end do
        end subroutine
        
        
        
        !> @brief Check if an equilibrium reaction belongs to the chemical system
        !> @details Searches eq_reacts array by reaction name, returns flag and optional index
        !> @param[in] this Chemical system object
        !> @param[in] react_name Reaction name to search for
        !> @param[out] flag TRUE if reaction found, FALSE otherwise
        !> @param[out] react_ind Optional index in eq_reacts array (0 if not found)
        subroutine is_eq_reaction_in_chem_syst(this,react_name,flag,react_ind)
        !> This subroutine checks if an equilibrium reaction belongs to the chemical system
            implicit none
            class(chem_system_c), intent(in) :: this !< chemical system
            character(len=*), intent(in) :: react_name !< reaction name
            logical, intent(out) :: flag !> TRUE if reaction belongs to chemical system, FALSE otherwise
            integer(kind=4), intent(out), optional :: react_ind !< index in attribute "eq_reacts" (0 if not present)
            
            integer(kind=4) :: i,sp_ind !> Loop counter and species index
            integer(kind=4), allocatable :: sp_indices(:) !> Species indices array
            logical :: sp_flag !> Species flag
            
            flag=.false. !> Initialize to not found
            if (present(react_ind)) then !> If index requested
                react_ind=0 !> Initialize to zero (not found)
            end if
            do i=1,this%speciation_alg%num_eq_reactions !> Loop through all equilibrium reactions
                if (this%eq_reacts(i)%name==react_name) then !> Compare names
                    flag=.true. !> Set found flag
                    if (present(react_ind)) then !> If index requested
                        react_ind=i !> Return index
                    end if
                    exit !> Exit loop early (found)
                end if
            end do
        end subroutine
        
        !> @brief Check if a kinetic reaction belongs to the chemical system
        !> @details Searches kin_reacts array by reaction name, returns flag and optional index
        !> @param[in] this Chemical system object
        !> @param[in] react_name Reaction name to search for
        !> @param[out] flag TRUE if reaction found, FALSE otherwise
        !> @param[out] react_ind Optional index in kin_reacts array (0 if not found)
         subroutine is_kin_reaction_in_chem_syst(this,react_name,flag,react_ind)
         !> This subroutine checks if a kinetic reaction belongs to the chemical system
             implicit none
             class(chem_system_c), intent(in) :: this !< chemical system
             character(len=*), intent(in) :: react_name !< reaction name
             logical, intent(out) :: flag !> TRUE if reaction belongs to chemical system, FALSE otherwise
             integer(kind=4), intent(out), optional :: react_ind !< index in attribute "kin_reacts" (0 if not present)
            
             integer(kind=4) :: i !> Loop counter
            
             flag=.false. !> Initialize to not found
             if (present(react_ind)) then !> If index requested
                 react_ind=0 !> Initialize to zero (not found)
             end if
             !> Search linear kinetic reactions (indices 1..num_lin_kin_reacts in kin_reacts)
             do i=1,this%num_lin_kin_reacts
                 if (this%lin_kin_reacts(i)%name==react_name) then
                     flag=.true.
                     if (present(react_ind)) react_ind=i
                     return
                 end if
             end do
             !> Search redox kinetic reactions (indices num_lin_kin+1..num_lin_kin+num_redox_kin in kin_reacts)
             do i=1,this%num_redox_kin_reacts
                 if (this%redox_kin_reacts(i)%name==react_name) then
                     flag=.true.
                     if (present(react_ind)) react_ind=this%num_lin_kin_reacts+i
                     return
                 end if
             end do
             !> Search mineral kinetic reactions (indices after redox in kin_reacts)
             do i=1,this%num_minerals_kin
                 if (this%min_kin_reacts(i)%name==react_name) then
                     flag=.true.
                     if (present(react_ind)) react_ind=this%num_lin_kin_reacts+this%num_redox_kin_reacts+i
                     return
                 end if
             end do
         end subroutine
!> @brief Separator line for GET methods section
!*********************** GET ***********************************************************************************************************************!
        !> @brief Get equilibrium constants of all equilibrium reactions
        !> @details Allocates and fills array with eq_cst values from all eq_reacts
        !> @param[in] this Chemical system object
        !> @return K Array of equilibrium constants
        function get_eq_csts(this) result(K)
        !> This function returns equilibrium constants of equilibrium reactions
            implicit none
            class(chem_system_c), intent(in) :: this !< Chemical system object
            real(kind=8), allocatable :: K(:) !< Array of equilibrium constants
            
            integer(kind=4) :: i !> Loop counter

            allocate(K(this%speciation_alg%num_eq_reactions)) !> Allocate array for all eq reactions
            do i=1,this%speciation_alg%num_eq_reactions !> Loop through all equilibrium reactions
                K(i)=this%eq_reacts(i)%eq_cst !> Extract equilibrium constant
            end do
        end function
!> @brief Separator line for REARRANGE methods section
!*********************** REARRANGE *****************************************************************************************************************!
        !> @brief Rearrange equilibrium reactions in canonical order
        !> @details Reorders eq_reacts array in specific sequence based on reaction types and activity flags.
        !> If flag_comp=TRUE, the canonical order is:
        !>   1. Constant activity minerals in equilibrium
        !>   2. Constant activity gases in equilibrium  
        !>   3. Redox equilibrium reactions
        !>   4. Aqueous complexes
        !>   5. Variable activity minerals in equilibrium
        !>   6. Cation exchange
        !>   7. Variable activity gases in equilibrium
        !> If flag_comp=FALSE, the canonical order is:
        !>   1. Redox equilibrium reactions
        !>   2. Aqueous complexes  
        !>   3. Variable activity minerals in equilibrium
        !>   4. Constant activity minerals in equilibrium
        !>   5. Cation exchange
        !>   6. Variable activity gases in equilibrium
        !>   7. Constant activity gases in equilibrium
        !> @param[in,out] this Chemical system object
        subroutine rearrange_eq_reacts(this)
            implicit none
            class(chem_system_c) :: this !> chemical system
            
            integer(kind=4) :: i,ind_min_cst_act,ind_aq,ind_gas_var_act,ind_surf,ind_gas_cst_act,ind_min_var_act,ind_redox !> Loop counter and index pointers
            type(eq_reaction_c), allocatable :: aux_eq_reacts(:) !> Temporary copy of reactions
            
            aux_eq_reacts=this%eq_reacts !> Save original reactions
            deallocate(this%eq_reacts) !> Deallocate original array
            allocate(this%eq_reacts(this%speciation_alg%num_eq_reactions)) !> Reallocate with same size
            if (this%speciation_alg%flag_comp .eqv. .true.) then !> If component matrix formulation
                !> we initialise counters
                ind_min_cst_act=1 !> constant activity minerals in equilibrium
                ind_gas_cst_act=ind_min_cst_act+this%num_minerals_eq_cst_act !> constant activity gases in equilibrium
                ind_redox=ind_gas_cst_act+THIS%gas_phase%num_gases_eq_cst_act !> redox equilibrium reactions
                ind_aq=ind_redox+this%num_redox_eq_reacts !> aqueous complexes
                ind_min_var_act=ind_aq+this%aq_phase%num_aq_complexes !> variable activity minerals
                ind_surf=ind_min_var_act+THIS%num_minerals_eq_var_act !> cation exchange
                ind_gas_var_act=ind_surf+this%cat_exch_zone%num_exch_cats !> variable activity gases
                do i=1,this%speciation_alg%num_eq_reactions !> Loop through all reactions
                    if (aux_eq_reacts(i)%react_type==2) then !> mineral dissolution/precipitation (type 2)
                        if (this%species(aux_eq_reacts(i)%species_ind(aux_eq_reacts(i)%num_species))%cst_act_flag.eqv..true.) then !> Check last species for constant activity
                            this%eq_reacts(ind_min_cst_act)=aux_eq_reacts(i) !> Place in cst act mineral section
                            ind_min_cst_act=ind_min_cst_act+1 !> Increment index
                        else !> Variable activity mineral
                            this%eq_reacts(ind_min_var_act)=aux_eq_reacts(i) !> Place in var act mineral section
                            ind_min_var_act=ind_min_var_act+1 !> Increment index
                        end if
                    else if (aux_eq_reacts(i)%react_type==1) then !> aqueous complex (type 1)
                        this%eq_reacts(ind_aq)=aux_eq_reacts(i) !> Place in aqueous complex section
                        ind_aq=ind_aq+1 !> Increment index
                    else if (aux_eq_reacts(i)%react_type==6) then !> gas (type 6)
                        if (this%species(aux_eq_reacts(i)%species_ind(aux_eq_reacts(i)%num_species))%cst_act_flag.eqv..true.) then !> Check last species for constant activity
                            this%eq_reacts(ind_gas_cst_act)=aux_eq_reacts(i) !> Place in cst act gas section
                            ind_gas_cst_act=ind_gas_cst_act+1 !> Increment index
                        else !> Variable activity gas
                            this%eq_reacts(ind_gas_var_act)=aux_eq_reacts(i) !> Place in var act gas section
                            ind_gas_var_act=ind_gas_var_act+1 !> Increment index
                        end if
                    else if (aux_eq_reacts(i)%react_type==3) then !> cation exchange (type 3)
                        this%eq_reacts(ind_surf)=aux_eq_reacts(i) !> Place in surface section
                        ind_surf=ind_surf+1 !> Increment index
                    else if (aux_eq_reacts(i)%react_type==4) then !> redox (type 4)
                        this%eq_reacts(ind_redox)=aux_eq_reacts(i) !> Place in redox section
                        ind_redox=ind_redox+1 !> Increment index
                    end if
                end do
            else
                !> If not component matrix formulation
                !> we initialise counters
                ind_redox=1 !> redox equilibrium reactions
                ind_aq=ind_redox+this%num_redox_eq_reacts !> aqueous complexes
                ind_min_var_act=ind_aq+this%aq_phase%num_aq_complexes !> variable activity minerals
                ind_min_cst_act=ind_min_var_act+THIS%num_minerals_eq_var_act !> constant activity minerals
                ind_surf=ind_min_cst_act+THIS%num_minerals_eq_cst_act !> cation exchange
                ind_gas_var_act=ind_surf+this%cat_exch_zone%num_exch_cats !> variable activity gases
                ind_gas_cst_act=ind_gas_var_act+THIS%gas_phase%num_gases_eq_var_act !> constant activity gases
                do i=1,this%speciation_alg%num_eq_reactions !> Loop through all reactions
                    if (aux_eq_reacts(i)%react_type==4) then !> redox (type 4)
                        this%eq_reacts(ind_redox)=aux_eq_reacts(i) !> Place in redox section
                        ind_redox=ind_redox+1 !> Increment index
                    else if (aux_eq_reacts(i)%react_type==1) then !> aqueous complex (type 1)
                        this%eq_reacts(ind_aq)=aux_eq_reacts(i) !> Place in aqueous complex section
                        ind_aq=ind_aq+1 !> Increment index
                    else if (aux_eq_reacts(i)%react_type==2) then !> mineral dissolution/precipitation (type 2)
                        if (this%species(aux_eq_reacts(i)%species_ind(aux_eq_reacts(i)%num_species))%cst_act_flag.eqv..true.) then !> Check last species for constant activity
                            this%eq_reacts(ind_min_cst_act)=aux_eq_reacts(i) !> Place in cst act mineral section
                            ind_min_cst_act=ind_min_cst_act+1 !> Increment index
                        else !> Variable activity mineral
                            this%eq_reacts(ind_min_var_act)=aux_eq_reacts(i) !> Place in var act mineral section
                            ind_min_var_act=ind_min_var_act+1 !> Increment index
                        end if
                    else if (aux_eq_reacts(i)%react_type==3) then !> cation exchange (type 3)
                        this%eq_reacts(ind_surf)=aux_eq_reacts(i) !> Place in surface section
                        ind_surf=ind_surf+1 !> Increment index
                    else if (aux_eq_reacts(i)%react_type==6) then !> gas (type 6)
                        if (this%species(aux_eq_reacts(i)%species_ind(aux_eq_reacts(i)%num_species))%cst_act_flag.eqv..true.) then !> Check last species for constant activity
                            this%eq_reacts(ind_gas_cst_act)=aux_eq_reacts(i) !> Place in cst act gas section
                            ind_gas_cst_act=ind_gas_cst_act+1 !> Increment index
                        else !> Variable activity gas
                            this%eq_reacts(ind_gas_var_act)=aux_eq_reacts(i) !> Place in var act gas section
                            ind_gas_var_act=ind_gas_var_act+1 !> Increment index
                        end if
                    end if
                end do
            end if
        end subroutine
        
        !> @brief Rearrange species array based on speciation algebra and surface complexation
        !> @details This complex subroutine reorders species array in different sequences depending on:
        !>   - Whether component matrix formulation is used (flag_comp)
        !>   - Whether cation exchange is present (flag_cat_exch)
        !> Four distinct ordering schemes are implemented for different combinations.
        !> Also populates var_act_sp_indices and cst_act_sp_indices arrays.
        !> @param[in,out] this Chemical system object
        subroutine rearrange_species(this)
        !< This subroutine rearranges the "species" attribute depending on the definition of the component matrix and the presence &
        !!    of surface complexes
            implicit none
            class(chem_system_c) :: this !> chemical system
            
            integer(kind=4) :: i,num_sp,num_aq_sec,num_var_act_sp,num_cst_act_sp !> Counters for species positioning
            num_sp=0 !> counter number of species
            num_var_act_sp=0 !> counter number of variable activity species
            num_cst_act_sp=this%aq_phase%wat_flag !> counter number of constant activity species (starts with water if present)
            if ((this%speciation_alg%flag_comp .eqv. .false.) .and. (this%speciation_alg%flag_cat_exch .eqv. .true.)) then
            !<      primary aqueous species
            !<      free surface
            !<      var act minerals NOT in equilibrium
            !<      cst act minerals NOT in equilibrium
            !<      var act gases NOT in equilibrium
            !<      cst act gases NOT in equilibrium
            !<      aqueous complexes
            !<      var act minerals in equilibrium
            !<      cst act minerals in equilibrium
            !<      surface complexes
            !<      var act gases in equilibrium
            !<      cst act gases in equilibrium
                do i=1,this%speciation_alg%num_aq_prim_species !> Loop through primary aqueous species
                    call this%species(i)%copy_species(this%aq_phase%aq_species(i)) !> Assign primary aqueous species to beginning of array
                end do
                call this%species(this%speciation_alg%num_aq_prim_species+1)%copy_species(this%cat_exch_zone%surf_compl(1)) !> free surface (immediately after primary aq)
                num_sp=num_sp+this%speciation_alg%num_aq_prim_species+1 !> Update species counter
                do i=1,this%num_minerals_kin !> Loop through kinetic minerals
                    call this%species(num_sp+i)%copy_species(this%minerals(i)%mineral) !> Assign kinetic minerals
                end do
                num_sp=num_sp+this%num_minerals_kin !> Update counter
                do i=1,this%gas_phase%num_gases_kin_var_act !> Loop through kinetic variable activity gases
                    call this%species(num_sp+i)%copy_species(this%gas_phase%gases(this%gas_phase%num_species-&
                        this%gas_phase%num_gases_kin_var_act+i)) !> Assign kinetic gases (skip equilibrium gases)
                end do
                num_sp=num_sp+this%gas_phase%num_gases_kin_var_act !> Update counter
                do i=1,this%gas_phase%num_gases_kin_cst_act !> Loop through kinetic constant activity gases
                    call this%species(num_sp+i)%copy_species(this%gas_phase%gases(this%gas_phase%num_gases_eq+i)) !> Assign kinetic gases (skip equilibrium gases)
                end do
                num_sp=num_sp+this%gas_phase%num_gases_kin_cst_act !> Update counter
                do i=1,this%aq_phase%num_aq_complexes !> Loop through aqueous complexes
                    call this%species(num_sp+i)%copy_species(this%aq_phase%aq_species(this%speciation_alg%num_aq_prim_species+i)) !> Assign secondary aqueous species
                end do
                num_sp=num_sp+this%aq_phase%num_aq_complexes !> Update counter
                do i=1,this%num_minerals_eq !> Loop through equilibrium minerals
                    call this%species(num_sp+i)%copy_species(this%minerals(this%num_minerals_kin+i)%mineral) !> Assign eq minerals (skip kinetic minerals)
                end do
                num_sp=num_sp+this%num_minerals_eq !> Update counter  
                do i=1,this%cat_exch_zone%num_exch_cats !> Loop through exchangeable cations
                    call this%species(num_sp+i)%copy_species(this%cat_exch_zone%surf_compl(1+i)) !> Assign surface complexes (skip free surface)
                end do
                num_sp=num_sp+this%cat_exch_zone%num_exch_cats !> Update counter
                do i=1,this%gas_phase%num_gases_eq_var_act !> Loop through equilibrium variable activity gases
                    call this%species(num_sp+i)%copy_species(this%gas_phase%gases(&
                        this%gas_phase%num_gases_eq_cst_act+&
                        i)) !> Assign gases equilibrium variable activity (skip kinetic gases)
                end do
                num_sp=num_sp+this%gas_phase%num_gases_eq_var_act !> Update counter
                do i=1,this%gas_phase%num_gases_eq_cst_act !> Loop through equilibrium constant activity gases
                    call this%species(num_sp+i)%copy_species(this%gas_phase%gases(&
                        i)) !> Assign gases equilibrium constant activity (skip kinetic gases)
                end do
                num_sp=num_sp+this%gas_phase%num_gases_eq_cst_act !> Update counter
            else if ((this%speciation_alg%flag_comp.eqv..true.) .and. (this%speciation_alg%flag_cat_exch.eqv..true.)) then
                !<      primary aqueous species
                !<      free surface
                !<      kinetic variable activity minerals
                !>      kinetic variable activity gases
                !<      secondary variable activity aqueous species
                !<      equilibrium variable activity minerals
                !<      surface complexes
                !<      variable activity gases in equilibrium
                !<      ideal water
                !<      kinetic constant activity minerals
                !<      equilibrium constant activity minerals
                !<      kinetic constant activity gases
                !<      constant activity gases in equilibrium
                call this%species(this%num_var_act_species+1)%copy_species(this%aq_phase%aq_species(this%aq_phase%ind_wat))
                do i=1,this%speciation_alg%num_aq_prim_species
                    call this%species(i)%copy_species(this%aq_phase%aq_species(this%aq_phase%ind_diss_solids(i))) !> chapuza porque asumes que los iones ya están ordenados en primarios y secundarios
                end do
                call this%species(this%speciation_alg%num_aq_prim_species+1)%copy_species(this%cat_exch_zone%surf_compl(1)) !> free surface
                num_var_act_sp=this%speciation_alg%num_aq_prim_species+1 !> update counter
                do i=1,this%num_minerals_kin_var_act
                    call this%species(num_var_act_sp+i)%copy_species(this%minerals(i)%mineral)
                end do
                num_var_act_sp=num_var_act_sp+this%num_minerals_kin_var_act
                do i=1,this%gas_phase%num_gases_kin_var_act
                    call this%species(num_var_act_sp+i)%copy_species(this%gas_phase%gases(&
                        this%gas_phase%num_species-this%gas_phase%num_gases_kin_var_act+i))
                end do
                num_var_act_sp=num_var_act_sp+this%gas_phase%num_gases_kin_var_act
                do i=1,this%speciation_alg%num_aq_sec_var_act_species
                    call this%species(this%speciation_alg%num_prim_species+i)%copy_species(&
                        this%aq_phase%aq_species(this%aq_phase%ind_diss_solids(this%speciation_alg%num_aq_prim_species+i))) !> secondary variable activity aqueous species (chapuza porque asumes que los iones ya están ordenados en primarios y secundarios)
                end do
                num_var_act_sp=num_var_act_sp+this%speciation_alg%num_aq_var_act_species
                do i=1,this%num_minerals_eq_var_act
                    call this%species(num_var_act_sp+i)%copy_species(this%minerals(this%num_minerals_kin+i)%mineral)
                end do
                num_var_act_sp=num_var_act_sp+this%num_minerals_eq_var_act
                ! do i=1,this%num_minerals
                !     if (this%minerals(i)%mineral%cst_act_flag.eqv..false.) then
                !         call this%species(num_var_act_sp+1)%copy_species(this%minerals(i)%mineral)
                !         num_var_act_sp=num_var_act_sp+1
                !     else
                !         call this%species(this%num_var_act_species+num_cst_act_sp+1)%copy_species(this%minerals(i)%mineral)
                !         num_cst_act_sp=num_cst_act_sp+1
                !     end if
                ! end do
                do i=1,this%cat_exch_zone%num_exch_cats
                    if (this%cat_exch_zone%surf_compl(1+i)%cst_act_flag.eqv..false.) then
                        call this%species(num_var_act_sp+1)%copy_species(this%cat_exch_zone%surf_compl(1+i))
                        num_var_act_sp=num_var_act_sp+1
                    else
                        call this%species(this%num_var_act_species+num_cst_act_sp+1)%copy_species(this%cat_exch_zone%surf_compl(1+i))
                        num_cst_act_sp=num_cst_act_sp+1
                    end if
                end do
                do i=1,this%gas_phase%num_gases_eq_var_act
                    call this%species(num_var_act_sp+i)%copy_species(this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+i))
                end do
                num_var_act_sp=num_var_act_sp+this%gas_phase%num_gases_eq_var_act
                ! do i=1,this%gas_phase%num_species
                !     if (this%gas_phase%gases(i)%cst_act_flag.eqv..false.) then
                !         call this%species(num_var_act_sp+1)%copy_species(this%gas_phase%gases(i))
                !         num_var_act_sp=num_var_act_sp+1
                !     else
                !         call this%species(this%num_var_act_species+num_cst_act_sp+1)%copy_species(this%gas_phase%gases(i))
                !         num_cst_act_sp=num_cst_act_sp+1
                !     end if
                ! end do
                do i=1,this%num_minerals_kin_cst_act
                    call this%species(num_cst_act_sp+i)%copy_species(this%minerals(this%num_minerals_kin_var_act+i)%mineral)
                end do
                num_cst_act_sp=num_cst_act_sp+this%num_minerals_kin_cst_act
                do i=1,this%num_minerals_eq_cst_act
                    call this%species(num_cst_act_sp+i)%copy_species(this%minerals(this%num_minerals_kin+&
                        this%num_minerals_eq_var_act+i)%mineral)
                end do
                num_cst_act_sp=num_cst_act_sp+this%num_minerals_eq_cst_act
                do i=1,this%gas_phase%num_gases_kin_cst_act
                    call this%species(num_cst_act_sp+i)%copy_species(this%gas_phase%gases(&
                        this%gas_phase%num_gases_eq+i))
                end do
                num_cst_act_sp=num_cst_act_sp+this%gas_phase%num_gases_kin_cst_act
                do i=1,this%gas_phase%num_gases_eq_cst_act
                    call this%species(num_cst_act_sp+i)%copy_species(this%gas_phase%gases(i))
                end do
                num_cst_act_sp=num_cst_act_sp+this%gas_phase%num_gases_eq_cst_act
            else if (this%speciation_alg%flag_comp .eqv. .true. .and. this%speciation_alg%flag_cat_exch .eqv. .false.) then
                !<      primary aqueous species
                !<      kinetic variable activity minerals
                !>      kinetic variable activity gases
                !<      secondary variable activity aqueous species
                !<      equilibrium variable activity minerals
                !<      variable activity gases in equilibrium
                !<      ideal water
                !<      kinetic constant activity minerals
                !<      equilibrium constant activity minerals
                !<      kinetic constant activity gases
                !<      constant activity gases in equilibrium
                call this%species(this%num_var_act_species+1)%copy_species(this%aq_phase%aq_species(this%aq_phase%ind_wat))
                do i=1,this%speciation_alg%num_aq_prim_species
                    call this%species(i)%copy_species(this%aq_phase%aq_species(this%aq_phase%ind_diss_solids(i))) !> chapuza porque asumes que los iones ya están ordenados en primarios y secundarios
                end do
                !call this%species(this%speciation_alg%num_aq_prim_species+1)%copy_species(this%cat_exch_zone%surf_compl(1)) !> free surface
                num_var_act_sp=this%speciation_alg%num_aq_prim_species !> update counter
                do i=1,this%num_minerals_kin_var_act
                    call this%species(num_var_act_sp+i)%copy_species(this%minerals(i)%mineral)
                end do
                num_var_act_sp=num_var_act_sp+this%num_minerals_kin_var_act
                do i=1,this%gas_phase%num_gases_kin_var_act
                    call this%species(num_var_act_sp+i)%copy_species(this%gas_phase%gases(&
                        this%gas_phase%num_species-this%gas_phase%num_gases_kin_var_act+i))
                end do
                num_var_act_sp=num_var_act_sp+this%gas_phase%num_gases_kin_var_act
                do i=1,this%speciation_alg%num_aq_sec_var_act_species
                    call this%species(this%speciation_alg%num_prim_species+i)%copy_species(&
                        this%aq_phase%aq_species(this%aq_phase%ind_diss_solids(this%speciation_alg%num_aq_prim_species+i))) !> secondary variable activity aqueous species (chapuza porque asumes que los iones ya están ordenados en primarios y secundarios)
                end do
                num_var_act_sp=num_var_act_sp+this%speciation_alg%num_aq_sec_var_act_species
                do i=1,this%num_minerals_eq_var_act
                    call this%species(num_var_act_sp+i)%copy_species(this%minerals(this%num_minerals_kin+i)%mineral)
                end do
                num_var_act_sp=num_var_act_sp+this%num_minerals_eq_var_act
                ! do i=1,this%num_minerals
                !     if (this%minerals(i)%mineral%cst_act_flag.eqv..false.) then
                !         call this%species(num_var_act_sp+1)%copy_species(this%minerals(i)%mineral)
                !         num_var_act_sp=num_var_act_sp+1
                !     else
                !         call this%species(this%num_var_act_species+num_cst_act_sp+1)%copy_species(this%minerals(i)%mineral)
                !         num_cst_act_sp=num_cst_act_sp+1
                !     end if
                ! end do
                !do i=1,this%cat_exch_zone%num_exch_cats
                !    if (this%cat_exch_zone%surf_compl(1+i)%cst_act_flag.eqv..false.) then
                !        call this%species(num_var_act_sp+1)%copy_species(this%cat_exch_zone%surf_compl(1+i))
                !        num_var_act_sp=num_var_act_sp+1
                !    else
                !        call this%species(this%num_var_act_species+num_cst_act_sp+1)%copy_species(this%cat_exch_zone%surf_compl(1+i))
                !        num_cst_act_sp=num_cst_act_sp+1
                !    end if
                !end do
                do i=1,this%gas_phase%num_gases_eq_var_act
                    call this%species(num_var_act_sp+i)%copy_species(this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+i))
                end do
                num_var_act_sp=num_var_act_sp+this%gas_phase%num_gases_eq_var_act
                num_sp=num_var_act_sp+num_cst_act_sp
                ! do i=1,this%gas_phase%num_species
                !     if (this%gas_phase%gases(i)%cst_act_flag.eqv..false.) then
                !         call this%species(num_var_act_sp+1)%copy_species(this%gas_phase%gases(i))
                !         num_var_act_sp=num_var_act_sp+1
                !     else
                !         call this%species(this%num_var_act_species+num_cst_act_sp+1)%copy_species(this%gas_phase%gases(i))
                !         num_cst_act_sp=num_cst_act_sp+1
                !     end if
                ! end do
                do i=1,this%num_minerals_kin_cst_act
                    call this%species(num_sp+i)%copy_species(this%minerals(this%num_minerals_kin_var_act+i)%mineral)
                end do
                num_sp=num_sp+this%num_minerals_kin_cst_act
                do i=1,this%num_minerals_eq_cst_act
                    call this%species(num_sp+i)%copy_species(this%minerals(this%num_minerals_kin+&
                        this%num_minerals_eq_var_act+i)%mineral)
                end do
                num_sp=num_sp+this%num_minerals_eq_cst_act
                do i=1,this%gas_phase%num_gases_kin_cst_act
                    call this%species(num_sp+i)%copy_species(this%gas_phase%gases(&
                        this%gas_phase%num_gases_eq+i))
                end do
                num_sp=num_sp+this%gas_phase%num_gases_kin_cst_act
                do i=1,this%gas_phase%num_gases_eq_cst_act
                    call this%species(num_sp+i)%copy_species(this%gas_phase%gases(i))
                end do
                num_sp=num_sp+this%gas_phase%num_gases_eq_cst_act
            else
            !<      primary aqueous species
            !<      var act minerals NOT in equilibrium
            !<      cst act minerals NOT in equilibrium
            !<      gases var act NOT in equilibrium
            !<      gases cst act NOT in equilibrium
            !<      aqueous complexes
            !<      var act minerals in equilibrium
            !<      cst act minerals in equilibrium
            !<      var act gases in equilibrium
            !<      cst act gases in equilibrium
                do i=1,this%speciation_alg%num_aq_prim_species
                    call this%species(i)%copy_species(this%aq_phase%aq_species(i))
                    if (this%aq_phase%aq_species(i)%cst_act_flag .eqv. .true.) then
                        this%cst_act_sp_indices(1)=i
                    else
                        num_var_act_sp=num_var_act_sp+1
                        this%var_act_sp_indices(num_var_act_sp)=i
                    end if
                end do
                num_sp=num_sp+this%speciation_alg%num_aq_prim_species
                do i=1,this%num_minerals_kin_var_act
                    call this%species(num_sp+i)%copy_species(this%minerals(i)%mineral)
                    this%var_act_sp_indices(num_var_act_sp+i)=num_sp+i
                end do
                num_var_act_sp=num_var_act_sp+this%num_minerals_kin_var_act
                num_sp=num_sp+this%num_minerals_kin_var_act
                do i=1,this%num_minerals_kin_cst_act
                    call this%species(num_sp+i)%copy_species(this%minerals(this%num_minerals_kin_var_act+i)%mineral)
                    this%cst_act_sp_indices(num_cst_act_sp+i)=num_sp+i
                end do
                num_sp=num_sp+this%num_minerals_kin_cst_act
                num_cst_act_sp=num_cst_act_sp+this%num_minerals_kin_cst_act
                do i=1,this%gas_phase%num_gases_kin_var_act
                    call this%species(num_sp+i)%copy_species(this%gas_phase%gases(this%gas_phase%num_gases_eq+&
                        this%gas_phase%num_gases_kin_cst_act+i))
                    this%var_act_sp_indices(num_var_act_sp+i)=num_sp+i
                end do
                num_var_act_sp=num_var_act_sp+this%gas_phase%num_gases_kin_var_act
                num_sp=num_sp+this%gas_phase%num_gases_kin_var_act
                do i=1,this%gas_phase%num_gases_kin_cst_act
                    call this%species(num_sp+i)%copy_species(this%gas_phase%gases(this%gas_phase%num_gases_eq+i))
                    this%cst_act_sp_indices(num_cst_act_sp+i)=num_sp+i
                end do
                num_sp=num_sp+this%gas_phase%num_gases_kin_cst_act
                num_cst_act_sp=num_cst_act_sp+this%gas_phase%num_gases_kin_cst_act
                do i=1,this%speciation_alg%num_sec_aq_species
                    call this%species(num_sp+i)%copy_species(this%aq_phase%aq_species(this%speciation_alg%num_aq_prim_species+i))
                    if (this%aq_phase%aq_species(this%speciation_alg%num_aq_prim_species+i)%cst_act_flag .eqv. .true.) then
                        this%cst_act_sp_indices(1)=num_sp+i
                    else
                        num_var_act_sp=num_var_act_sp+1
                        this%var_act_sp_indices(num_var_act_sp)=num_sp+i
                    end if
                end do
                num_sp=num_sp+this%speciation_alg%num_sec_aq_species
                do i=1,this%num_minerals_eq_var_act
                    call this%species(num_sp+i)%copy_species(this%minerals(this%num_minerals_kin+i)%mineral)
                    this%var_act_sp_indices(num_var_act_sp+i)=num_sp+i
                end do
                num_sp=num_sp+this%num_minerals_eq_var_act
                num_var_act_sp=num_var_act_sp+this%num_minerals_eq_var_act
                do i=1,this%num_minerals_eq_cst_act
                    call this%species(num_sp+i)%copy_species(this%minerals(this%num_minerals_kin+&
                        this%num_minerals_eq_var_act+i)%mineral)
                    this%cst_act_sp_indices(num_cst_act_sp+i)=num_sp+i
                end do
                num_sp=num_sp+this%num_minerals_eq_cst_act
                num_cst_act_sp=num_cst_act_sp+this%num_minerals_eq_cst_act
                do i=1,this%gas_phase%num_gases_eq_var_act
                    call this%species(num_sp+i)%copy_species(this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+i))
                    this%var_act_sp_indices(num_var_act_sp+i)=num_sp+i
                end do
                num_var_act_sp=num_var_act_sp+this%gas_phase%num_gases_eq_var_act
                num_sp=num_sp+this%gas_phase%num_gases_eq_var_act
                do i=1,this%gas_phase%num_gases_eq_cst_act
                    call this%species(num_sp+i)%copy_species(this%gas_phase%gases(i))
                    this%cst_act_sp_indices(num_cst_act_sp+i)=num_sp+i
                end do
                num_sp=num_sp+this%gas_phase%num_gases_eq_cst_act
                num_cst_act_sp=num_cst_act_sp+this%gas_phase%num_gases_eq_cst_act
            end if
        end subroutine
        
        !> @brief Compute squared charges (valences) of all species
        !> @details Allocates z2 array and fills with square of each species' valence for ionic strength calculations
        !> @param[in,out] this Chemical system object
        subroutine compute_z2(this)
        !> This subroutine computes attribute "z2"
            implicit none
            class(chem_system_c) :: this !< chemical system
            
            integer(kind=4) :: i !> Loop counter
            if (.not. allocated(this%z2)) then !> If z2 not yet allocated
                allocate(this%z2(this%speciation_alg%num_species)) !> Allocate for all species
            end if
            do i=1,this%speciation_alg%num_species !> Loop through all species
                this%z2(i)=this%species(i)%valence**2 !> Compute squared valence
            end do
        end subroutine
        
        !> @brief Compute total number of solid phases
        !> @details Computes as sum of minerals and surface complexes
        !> @param[in,out] this Chemical system object
        subroutine compute_num_solids_chem_syst(this)
        !> This subroutine computes the "num_solids" attribute 
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            this%num_solids=this%num_minerals+this%cat_exch_zone%num_surf_compl !> Sum minerals and surface complexes
        end subroutine
        
        !> @brief Compute total number of reactions
        !> @details Computes as sum of equilibrium and kinetic reactions
        !> @param[in,out] this Chemical system object
        subroutine compute_num_reacts(this)
        !> This subroutine computes the "num_reacts" attribute 
            implicit none
            class(chem_system_c) :: this !< Chemical system object
            this%num_reacts=this%speciation_alg%num_eq_reactions+this%num_kin_reacts !> Sum equilibrium and kinetic reactions
        end subroutine
        
        !> @brief Compute equilibrium constants for constant activity gases
        !> @details Modifies eq_cst values by dividing by partial pressures for gases with constant activity
        !> @param[in,out] this Chemical system object
        !> @param[in] part_pressures Partial pressures of constant activity equilibrium gases
        !> @param[in] ind_gases_eq_cst_act Indices of constant activity equilibrium gases in gas phase
        subroutine compute_eq_csts_gases_cst_act(this,part_pressures,ind_gases_eq_cst_act)
            class(chem_system_c) :: this !< Chemical system object
            real(kind=8), intent(in) :: part_pressures(:) !> partial pressures of gases in equilibrium with constant activity
            integer(kind=4), intent(in) :: ind_gases_eq_cst_act(:) !> indices of gases in equilibrium with constant activity in gas phase object
            
            integer(kind=4) :: i !> Loop counter
            do i=1,size(ind_gases_eq_cst_act) !> Loop through constant activity gases
                this%eq_reacts(this%num_aq_eq_reacts+&
                    this%num_minerals_eq+ind_gases_eq_cst_act(i))%eq_cst = this%eq_reacts(&
                    this%num_aq_eq_reacts+this%num_minerals_eq + & 
                    ind_gases_eq_cst_act(i))%eq_cst/part_pressures(ind_gases_eq_cst_act(i)) !> Divide eq_cst by partial pressure
            end do
        end subroutine
        
   

    !> @brief Change speciation algebra to exclude constant activity species from component matrix
    !> @details Modifies the chemical system to use a different speciation algebra formulation.
    !> Eliminates constant activity species from component matrix, rearranges species and reactions,
    !> recomputes stoichiometric matrices, and updates speciation algebra arrays.
    !> @param[in,out] this Chemical system object
    !> @param[in] tol Tolerance for speciation algorithm convergence
    subroutine change_spec_alg_chem_syst(this,tol)
    !< This subroutine changes the speciation algebra object of the chemical system from considering constant activity species in &
        !! the component matrix to not considering them
    implicit none
    class(chem_system_c) :: this !> chemical system
    real(kind=8), intent(in) :: tol !> tolerance for speciation algorithm

    logical :: flag_Se !> flag for computing Se matrix
    integer(kind=4), allocatable :: cols(:) !> columns of Se matrix to be swapped
    real(kind=8), allocatable :: eq_csts(:) !> equilibrium constants of equilibrium reactions
    type(species_c), allocatable :: aux_Species(:) !> old species array before rearrangement
    integer(kind=4) :: i !> loop counter

    allocate(cols(2)) !> Allocate columns array for swapping

    call this%speciation_alg%elim_cst_act_species(this%cat_exch_zone%num_surf_compl,&
        this%aq_phase%num_species,this%num_minerals_kin,&
        this%gas_phase%num_species-this%gas_phase%num_gases_eq) !> Eliminate constant activity species from speciation algebra
    aux_Species=this%species !> Save old species ordering before rearrangement
    call this%rearrange_species() !> Rearrange species array based on new formulation
    call this%compute_z2() !> Recompute squared charges
    call this%rearrange_eq_reacts() !> Rearrange equilibrium reactions
    !> Rearrange species indices in all reactions
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
    call this%set_stoich_mat() !> Rebuild global stoichiometric matrix
    call this%set_stoich_mat_gas() !> Rebuild gas stoichiometric matrix
    call this%set_stoich_mat_sol() !> Rebuild solid stoichiometric matrix
    eq_csts=this%get_eq_csts() !> Get equilibrium constants
    call this%speciation_alg%compute_arrays(this%Se,eq_csts,tol,flag_Se,cols) !> Recompute speciation algebra arrays with new formulation
    end subroutine

    subroutine is_redox_kin_reaction_in_chem_syst(this,react_name,flag,react_ind)
        !> @brief Check if redox kinetic reaction is in chemical system
        !> @details Searches for redox kinetic reaction by name in chemical system's redox kinetic reactions
        !> @param[in] this Chemical system object
        !> @param[in] react_name Name of redox kinetic reaction to search for
        !> @param[out] flag TRUE if reaction found, FALSE otherwise
        !> @param[out] react_ind (optional) Index of found reaction in redox kinetic reactions array
        implicit none
        class(chem_system_c), intent(in) :: this !> Chemical system object
        character(len=*), intent(in) :: react_name !> Name of redox kinetic reaction to search for
        logical, intent(out) :: flag !> TRUE if reaction found, FALSE otherwise
        integer(kind=4), intent(out), optional :: react_ind !> Index of found reaction in redox kinetic reactions array

        integer(kind=4) :: i !> Loop counter

        flag=.false. !> Initialize found flag to FALSE
        if (present(react_ind)) then !> If index requested
            react_ind=0 !> Initialize to zero (not found)
        end if
        do i=1,this%num_redox_kin_reacts !> Loop through all redox kinetic reactions
            if (this%redox_kin_reacts(i)%name==react_name) then !> Compare names
                flag=.true. !> Set found flag
                if (present(react_ind)) then !> If index requested
                    react_ind=i !> Return index
                end if
                exit !> Exit loop early (found)
            end if
        end do
    end subroutine is_redox_kin_reaction_in_chem_syst
    
    subroutine get_mineral_index_by_name(this,min_name,min_index)
        !> @brief Get index of mineral by name in chemical system
        !> @details Searches for mineral by name in chemical system's minerals array
        !> @param[in] this Chemical system object
        !> @param[in] min_name Name of mineral to search for
        !> @param[out] min_index Index of found mineral in minerals array
        implicit none
        class(chem_system_c), intent(in) :: this !> Chemical system object
        character(len=*), intent(in) :: min_name !> Name of mineral to search for
        integer(kind=4), intent(out) :: min_index !> Index of found mineral in minerals array
        
        integer(kind=4) :: i !> Loop counter
        !found=.false. !> Initialize found flag to FALSE
        min_index=0 !> Initialize index to 0 (not found)
        do i=1,this%num_minerals !> Loop through all minerals
            if (this%minerals(i)%mineral%name==min_name) then !> Compare names
                !found=.true. !> Set found flag
                min_index=i !> Return index
                exit !> Exit loop early (found)
            end if
        end do
        if (min_index==0) then
            write(*,*) 'Error: Mineral ',trim(min_name),' not found in chemical system.'
            error stop
        end if
    end subroutine get_mineral_index_by_name
    
    subroutine is_reaction_in_chem_syst(this,react_name,flag,react_ind)
        !> @brief Check if equilibrium reaction is in chemical system
        !> @details Searches for equilibrium reaction by name in chemical system's equilibrium reactions
        !> @param[in] this Chemical system object
        !> @param[in] react_name Name of equilibrium reaction to search for
        !> @param[out] flag 0 if not found, 1 if equilibrium, 2 if kinetic
        !> @param[out] react_ind: Index of found reaction in equilibrium reactions array
        implicit none
        class(chem_system_c), intent(in) :: this !> Chemical system object
        character(len=*), intent(in) :: react_name !> Name of equilibrium reaction to search for
        integer(kind=4), intent(out) :: flag !> 0 if not found, 1 if equilibrium, 2 if kinetic
        integer(kind=4), intent(out) :: react_ind !> Index of found reaction in equilibrium or kinetic reactions array (0 if not found)
        integer(kind=4) :: i !> Loop counter
        flag=0 !> Initialize found flag to 0 (not found)
        react_ind=0 !> Initialize to zero (not found)
        do i=1,this%speciation_alg%num_eq_reactions !> Loop through all equilibrium reactions
            if (this%eq_reacts(i)%name==react_name) then !> Compare names
                flag=1 !> Set found flag
                react_ind=i !> Return index
                exit !> Exit loop early (found)
            end if
        end do
        if (flag==0) then !> If not found in equilibrium reactions
            do i=1,this%num_kin_reacts !> Loop through all kinetic reactions
                if (.not. associated(this%kin_reacts(i)%kin_reaction)) cycle
                if (this%kin_reacts(i)%kin_reaction%name==react_name) then !> Compare names
                    flag=2 !> Set found flag
                    react_ind=i !> Return index
                    exit !> Exit loop early (found)
                end if
            end do
        end if
        if (flag==0) then
            write(*,*) 'Error: Reaction ',trim(react_name),' not found in chemical system.'
            error stop
        end if
    end subroutine is_reaction_in_chem_syst
    
    !subroutine allocate_cat_exch_zones(this,ncez)
    !class(chem_system_c), intent(inout) :: this !> Chemical system object
    !integer(kind=4), intent(in), optional :: ncez !> Number of gas phases to allocate
    !if (present(ncez)) then
    !    if (ncez < 0) then
    !        error stop 'Error: Number of cation exchange zones must be non-negative.'
    !    end if
    !    this%num_cat_exch_zones=ncez
    !end if
    !if (.not. allocated(this%cat_exch_zones)) then
    !    allocate(this%cat_exch_zones(this%num_cat_exch_zones))
    !end if
    !end subroutine allocate_cat_exch_zones
!> @brief End of chemical system module
end module