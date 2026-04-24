!> \file reaction_m.f90
!> \brief Abstract reaction superclass module
!> \details
!>   This module defines the abstract base class for all chemical reactions.
!>   
!>   Reaction is the fundamental superclass for:
!>   - Equilibrium reactions (aqueous complexation, mineral dissolution, cation exchange)
!>   - Kinetic reactions (linear, redox, Monod, mineral kinetics)
!>   - Gas phase reactions
!>   
!>   Key features:
!>   - Polymorphic species array (reactants and products)
!>   - Stoichiometric coefficients (positive for products, negative for reactants)
!>   - Equilibrium constant (K_eq) with temperature dependence
!>   - Reaction classification system via react_type
!>   
!>   Reaction types:
!>   - 1: Aqueous complexation (secondary species formation)
!>   - 2: Mineral dissolution/precipitation
!>   - 3: Cation exchange
!>   - 4: Redox reaction
!>   - 5: Linear kinetic reaction
!>   - 6: Gas phase reaction
!>   - 7: Nonlinear kinetic reaction
!>   
!>   Temperature dependence of K_eq:
!>   \f[ \log_{10}(K) = a_1 + a_2 T + \frac{a_3}{T} + a_4 \log_{10}(T) + \frac{a_5}{T^2} + a_6 T^2 \f]
!>   where coeffs_logK_T = [a_1, a_2, a_3, a_4, a_5, a_6]
!>   
!> \author Jordi Petchamé-Guerrero
!> \date October 2025
module reaction_m
    use species_m, only: species_c !< Import species base class for polymorphic species array
    implicit none !< Enforce explicit variable declarations for type safety
    private !< Make all entities private by default
    public :: write_reaction_sup !< Public procedure for writing reaction information (available to all derived types)
    !> \brief Abstract reaction superclass
    !> \details
    !>   Base class for all chemical reactions in the system.
    !>   
    !>   This abstract type defines the core attributes and operations common to all reactions:
    !>   - Species participation (reactants and products)
    !>   - Stoichiometric coefficients
    !>   - Equilibrium constants
    !>   - Temperature dependence
    !>   
    !>   Stoichiometric convention:
    !>   - Negative coefficients: Reactants (consumed)
    !>   - Positive coefficients: Products (formed)
    !>   
    !>   Example (calcite dissolution):
    !>   CaCO3(s) = Ca²⁺ + CO3²⁻
    !>   - species = [CaCO3, Ca²⁺, CO3²⁻]
    !>   - stoichiometry = [-1.0, 1.0, 1.0]
    !>   
    !>   Derived classes:
    !>   - eq_reaction_c: Equilibrium reactions (instantaneous)
    !>   - kin_reaction_c: Kinetic reactions (rate-based)
    !>   - lin_kin_reaction_c: Linear kinetics
    !>   - redox_kin_reaction_c: Redox kinetics
    !>   - kin_mineral_c: Mineral kinetics
    type, public, abstract :: reaction_c !> reaction superclass
        integer(kind=4), allocatable :: species_ind(:)              !< Indices in chemical system of species involved (reactants + products)
        real(kind=8), allocatable :: stoichiometry(:)            !< Stoichiometric coefficients (same order as species)
        character(len=256) :: name                                !< Reaction name (e.g., "calcite_dissolution")
        integer(kind=4) :: num_species                            !< Number of species involved
        integer(kind=4) :: react_type                             !< Reaction type code (1-7, see module header)
        real(kind=8) :: eq_cst                                    !< Equilibrium constant K_eq at reference T
        real(kind=8), allocatable :: coeffs_logK_T(:)             !< 6 coefficients for log(K) temperature dependence
        real(kind=8) :: delta_h                                   !< Enthalpy of reaction [kJ/mol]
    contains
    !> Set procedures
        procedure :: set_stoichiometry                    !< Set stoichiometric coefficients
        procedure :: set_all_species                      !< Set all species at once
        procedure :: set_species_indices_from_names                  !< Set species indices from array
        procedure :: rearrange_species_indices                  !< Set species indices from array
        procedure :: set_single_species                   !< Set individual species by index
        procedure :: set_eq_cst                           !< Set equilibrium constant
        procedure :: set_delta_h                          !< Set enthalpy
        procedure :: set_react_name                       !< Set reaction name
        procedure :: set_react_type                       !< Set reaction type code
    !> Allocate/deallocate procedures
        procedure :: allocate_reaction                    !< Allocate species and stoichiometry arrays
        procedure :: deallocate_reaction                  !< Deallocate all allocatable components
    !> Query procedures
        procedure :: is_species_in_react                  !< Check if single species participates
        procedure :: are_species_in_react                 !< Check if multiple species participate
    !> Computation procedures
        procedure :: compute_logK_dep_T                   !< Compute K_eq at given temperature
    !> Write procedures
        procedure :: write_reaction=>write_reaction_sup   !< Write reaction to file
    !> Utility procedures
        procedure :: copy_attributes                      !< Copy all attributes from another reaction
        procedure :: change_sign_stoichiometry            !< Reverse reaction direction
    end type
    
    abstract interface
       
    end interface
    
contains
    
       !> @brief Allocate reaction arrays
       !> @details Allocates memory for species, stoichiometry, and temperature coefficients.
       !> If num_species is provided, updates this%num_species before allocation.
       !> Always allocates 6 coefficients for temperature-dependent K_eq calculation.
       !> @param[in,out] this Reaction object to allocate
       !> @param[in] num_species Optional number of species (updates this%num_species if present) [-]
       subroutine allocate_reaction(this,num_species)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(inout) :: this !< Reaction object to allocate [-]
            integer(kind=4), intent(in), optional :: num_species !< Optional number of species to set [-]
            
            !> Update number of species if provided as argument
            if (present(num_species)) then
                this%num_species=num_species !< Set number of species from optional argument [-]
            end if
            !> Allocate all dynamic arrays: stoichiometry coefficients, species array, and 6 temperature coefficients
            allocate(this%stoichiometry(this%num_species),this%species_ind(this%num_species),this%coeffs_logK_T(6)) !< Allocate arrays [-]
       end subroutine
       
       !> @brief Deallocate reaction arrays
       !> @details Frees memory for stoichiometry, species, and temperature coefficient arrays.
       !> Should be called before destroying reaction object to prevent memory leaks.
       !> @param[in,out] this Reaction object to deallocate
       subroutine deallocate_reaction(this)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(inout) :: this !< Reaction object to deallocate [-]
            
            !> Deallocate all dynamic arrays in order
            deallocate(this%stoichiometry,this%species_ind,this%coeffs_logK_T) !< Free memory for all allocatable components [-]
       end subroutine
       
       !> @brief Set stoichiometric coefficients
       !> @details Assigns stoichiometric coefficients for all species in reaction.
       !> Validates that array size matches num_species if species array is allocated.
       !> Convention: negative for reactants, positive for products.
       !> @param[in,out] this Reaction object to configure
       !> @param[in] stoichiometry Array of stoichiometric coefficients (same order as species) [-]
       subroutine set_stoichiometry(this,stoichiometry)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(inout) :: this !< Reaction object to configure [-]
            real(kind=8), intent(in) :: stoichiometry(:) !< Stoichiometric coefficients array (negative=reactants, positive=products) [-]
            
            !> Validate dimension if species array already allocated
            if (allocated(this%species_ind) .and. size(stoichiometry)/=this%num_species) then
                error stop "Dimension error in 'set_stoichiometry'" !< Halt if size mismatch detected
            end if
            this%stoichiometry=stoichiometry !< Assign stoichiometric coefficients [-]
       end subroutine
       
       !> @brief Set all species at once
       !> @details Assigns entire species array via polymorphic assignment.
       !> Species array must be pre-allocated with correct size.
       !> @param[in,out] this Reaction object to configure
       !> @param[in] species Array of species objects (polymorphic) [-]
       subroutine set_all_species(this,species_ind)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(inout) :: this !< Reaction object to configure [-]
            integer(kind=4), intent(in) :: species_ind(:) !< Array of species (reactants and products) [-]
            
            this%species_ind=species_ind !< Assign entire species array via polymorphic copy [-]
       end subroutine
       
       !> @brief Set individual species by index
       !> @details Assigns a single species at specified position in species array.
       !> Validates that index is within valid range [1, num_species].
       !> Uses polymorphic copy_species method for proper type handling.
       !> @param[in,out] this Reaction object to configure
       !> @param[in] species Species object to assign [-]
       !> @param[in] index Position in species array (1-based indexing) [-]
       subroutine set_single_species(this,react_ind,chem_syst_ind)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(inout) :: this !< Reaction object to configure [-]
            integer(kind=4), intent(in) :: react_ind !< Index of reactant species [-]
            integer(kind=4), intent(in) :: chem_syst_ind !< Index of species in chemical system [-]
            
            !> Validate that index is within valid bounds
            if (react_ind<1 .or. react_ind>this%num_species) error stop "Index out of bounds in set_single_species" !< Halt if invalid index
            this%species_ind(react_ind)=chem_syst_ind
            !call this%species(react_ind)%copy_species(chem_syst_species(chem_syst_ind)) !< Assign species using polymorphic method [-]
       end subroutine

       !> @brief Set species names from array
       !> @details Assigns names to all species from string array.
       !> Validates that array size matches num_species.
       !> Iterates through species calling set_name for each.
       !> @param[in,out] this Reaction object to configure
       !> @param[in] species_indices Array of species indices (integers) [-]
       subroutine set_species_indices_from_names(this,species_names,chem_syst_species)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(inout) :: this !< Reaction object to configure [-]
            character(len=*), intent(in) :: species_names(:) !< Array of species names (must match num_species) [-]
            class(species_c), intent(in) :: chem_syst_species(:) !< Array of species from the chemical system [-]
            integer(kind=4) :: i,j !< Loop counters [-]
            logical :: found !< Flag indicating if species was found [-]
            
            !> Validate dimension match
            if (size(species_names).ne.this%num_species) error stop "Dimension error in 'set_species_indices'" !< Halt if size mismatch
            !> Allocate species_ind if not already allocated
            if (.not. allocated(this%species_ind)) then
                allocate(this%species_ind(this%num_species))
            end if
            !> Iterate over all species names and look up in the chemical system
            do i=1,this%num_species
                found=.false.

                do j=1,size(chem_syst_species)
                    if (species_names(i)==chem_syst_species(j)%name) then
                        this%species_ind(i)=j
                        found=.true.
                        exit
                    end if
                end do
                if (.not. found) then
                    error stop "Species '"//trim(species_names(i))//"' not found in chemical system"
                end if
            end do
       end subroutine

       !> @brief Set equilibrium constant
       !> @details Assigns equilibrium constant K_eq at reference temperature.
       !> For equilibrium reactions: K_eq = ∏[products]^ν / ∏[reactants]^ν
       !> @param[in,out] this Reaction object to configure
       !> @param[in] eq_cst Equilibrium constant K_eq (dimensionless) [-]
       subroutine set_eq_cst(this,eq_cst)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(inout) :: this !< Reaction object to configure [-]
            real(kind=8), intent(in) :: eq_cst !< Equilibrium constant value at reference temperature [-]
            
            this%eq_cst=eq_cst !< Assign equilibrium constant [-]
       end subroutine
       
       !> @brief Set enthalpy of reaction
       !> @details Assigns standard enthalpy change ΔH for the reaction.
       !> Used for temperature-dependent K_eq calculations via van't Hoff equation.
       !> @param[in,out] this Reaction object to configure
       !> @param[in] delta_h Enthalpy of reaction [kJ/mol]
       subroutine set_delta_h(this,delta_h)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(inout) :: this !< Reaction object to configure [-]
            real(kind=8), intent(in) :: delta_h !< Standard enthalpy change of reaction [kJ/mol]
            
            this%delta_h=delta_h !< Assign enthalpy of reaction [kJ/mol]
       end subroutine
       
       !> @brief Set reaction name
       !> @details Assigns descriptive name to reaction for identification.
       !> Examples: "calcite_dissolution", "pyrite_oxidation", etc.
       !> @param[in,out] this Reaction object to configure
       !> @param[in] name Reaction name string (up to 256 characters) [-]
       subroutine set_react_name(this,name)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(inout) :: this !< Reaction object to configure [-]
            character(len=*), intent(in) :: name !< Reaction name for identification [-]
            
            this%name=name !< Assign reaction name [-]
       end subroutine
       
       !> @brief Set reaction type code
       !> @details Assigns reaction classification: 1=aqueous, 2=mineral, 3=exchange,
       !> 4=redox, 5=linear kinetic, 6=gas, 7=nonlinear kinetic.
       !> See module header for complete type descriptions.
       !> @param[in,out] this Reaction object to configure
       !> @param[in] react_type Reaction type code (1-7) [-]
       subroutine set_react_type(this,react_type)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(inout) :: this !< Reaction object to configure [-]
            integer(kind=4), intent(in) :: react_type !< Reaction type classification code (1-7) [-]
            
            this%react_type=react_type !< Assign reaction type code [-]
       end subroutine
       
       !> @brief Check if single species participates in reaction
       !> @details Searches for species by name in reaction's species list.
       !> Returns flag=true if found, along with optional index position.
       !> Uses name comparison for species matching.
       !> @param[in] this Reaction object to search
       !> @param[in] species Species to search for [-]
       !> @param[out] flag True if species found, false otherwise [-]
       !> @param[out] species_ind Optional: index of species in reaction (0 if not found) [-]
       subroutine is_species_in_react(this,species_ind_chem_syst,flag,species_ind_react)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(in) :: this !< Reaction object to search [-]
            integer(kind=4), intent(in) :: species_ind_chem_syst !< Species to search for in reaction [-]
            logical, intent(out) :: flag !< True if species found in reaction, false otherwise [-]
            integer(kind=4), intent(out), optional :: species_ind_react !< Optional: index of species in reaction (0 if not found) [-]
            
            integer(kind=4) :: i !< Loop counter for iterating over reaction species [-]
            
            flag=.false. !< Initialize flag to false (not found) [-]
            !> Initialize optional species index to zero if present
            if (present(species_ind_react)) then
                species_ind_react=0 !< Set to zero indicating species not yet found [-]
            end if
            !> Search for species by name in reaction's species list
            do i=1,this%num_species !< Iterate through all species in reaction
                if (species_ind_chem_syst==this%species_ind(i)) then !< Compare species names for match
                    flag=.true. !< Set flag to true indicating species found [-]
                    if (present(species_ind_react)) then
                        species_ind_react=i !< Store index where species was found [-]
                    end if
                    exit !< Exit loop once species found (no need to continue searching)
                end if
            end do
       end subroutine
    !******************************************************************************************************************************
       !> \brief Check if multiple species participate in the reaction
       !> \details
       !>   Searches for a set of species in the reaction.
       !>   Returns true only if ALL species are found.
       !>   
       !>   Algorithm:
       !>   - Searches for each query species in reaction's species list
       !>   - Returns indices in order of query species array
       !>   - Sets flag=false if any species is not found
       !>   
       !>   species_ind must be pre-allocated with size matching species array.
       !> \param[in] this The reaction object (const)
       !> \param[in] species Array of species to search for
       !> \param[out] flag True if ALL species participate, false otherwise
       !> \param[out] species_ind Indices in reaction (0 if species not found)
       !> @brief Check if multiple species participate in reaction
       !> @details Searches for set of species in reaction's species list.
       !> Returns true only if ALL query species are found.
       !> species_ind array must be pre-allocated with size matching species array.
       !> @param[in] this Reaction object to query
       !> @param[in] species Array of species to search for [-]
       !> @param[out] flag True if ALL species found, false if any missing [-]
       !> @param[out] species_ind Indices in reaction (0 if any species not found) [-]
       subroutine are_species_in_react(this,species_ind_chem_syst,flag,species_ind_react)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(in) :: this !< Reaction object to search [-]
            integer(kind=4), intent(in) :: species_ind_chem_syst(:) !< Array of species indices to search for in reaction [-]
            logical, intent(out) :: flag !< True if all species found, false if any missing [-]
            integer(kind=4), intent(out) :: species_ind_react(:) !< Indices of species in reaction (pre-allocated) [-]
            
            integer(kind=4) :: i,j !< Loop counters: i for reaction species, j for query species [-]
            
            !> Validate that species and species_ind arrays have matching dimensions
            if (size(species_ind_chem_syst)/=size(species_ind_react)) then
                error stop "Dimension error in subroutine 'are_species_in_react'" !< Halt if array size mismatch
            end if
            flag=.true. !< Initialize flag to true (assume all found until proven otherwise) [-]
            i=1 !< Initialize reaction species index to 1 [-]
            j=1 !< Initialize query species index to 1 [-]
            !> Nested search loop: for each query species, search entire reaction species list
            do
                if (species_ind_chem_syst(j)==this%species_ind(i)) then !< Check if current query species matches current reaction species
                    species_ind_react(j)=i !< Store index where query species j was found in reaction [-]
                    if (j<size(species_ind_chem_syst)) then !< Check if more query species remain
                        j=j+1 !< Move to next query species [-]
                        i=1 !< Reset reaction species index to start new search [-]
                    else
                        exit !< All query species found, exit successfully
                    end if
                else if (i<this%num_species) then !< Current query species not matched, more reaction species to check
                    i=i+1 !< Move to next reaction species [-]
                else !< Reached end of reaction species without finding match
                    flag=.false. !< Set flag to false indicating not all species found [-]
                    exit !< Exit loop early since at least one query species missing
                end if
            end do
       end subroutine
    !******************************************************************************************************************************
        !> \brief Copy all attributes from another reaction
        !> \details
        !>   Deep copy of all reaction attributes:
        !>   - species array (polymorphic copy)
        !>   - stoichiometry coefficients
        !>   - name and type
        !>   - thermodynamic properties (eq_cst, delta_h, coeffs_logK_T)
        !>   
        !>   Both reactions must have arrays allocated before calling.
        !> \param[in,out] this The reaction object (destination)
        !> \param[in] reaction Source reaction to copy from
        !> @brief Copy all attributes from another reaction
        !> @details Deep copy of all reaction properties including polymorphic species array.
        !> Both reactions must have arrays allocated before calling.
        !> Copies: species, stoichiometry, name, type, public, thermodynamic properties.
        !> @param[in,out] this Destination reaction object
        !> @param[in] reaction Source reaction to copy from [-]
        subroutine copy_attributes(this,reaction)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(inout) :: this !< Destination reaction (modified) [-]
            class(reaction_c), intent(in) :: reaction !< Source reaction to copy from (const) [-]
            
            this%species_ind=reaction%species_ind !< Copy polymorphic species array [-]
            this%stoichiometry=reaction%stoichiometry !< Copy stoichiometric coefficients [-]
            this%name=reaction%name !< Copy reaction name [-]
            this%num_species=reaction%num_species !< Copy number of species [-]
            this%react_type=reaction%react_type !< Copy reaction type code [-]
            this%delta_h=reaction%delta_h !< Copy enthalpy of reaction [kJ/mol]
            this%eq_cst=reaction%eq_cst !< Copy equilibrium constant [-]
            this%coeffs_logK_T=reaction%coeffs_logK_T !< Copy temperature dependence coefficients [-]
        end subroutine
    !******************************************************************************************************************************
        !> \brief Compute temperature-dependent equilibrium constant
        !> \details
        !>   Computes K_eq at specified temperature using PHREEQC analytical expression:
        !>   \f[ \log_{10}(K) = a_1 + a_2 T + \frac{a_3}{T} + a_4 \log_{10}(T) + \frac{a_5}{T^2} + a_6 T^2 \f]
        !>   
        !>   where T is temperature [K] and coeffs_logK_T = [a_1, ..., a_6]
        !>   
        !>   Updates this%eq_cst with K = 10^(logK)
        !> \param[in,out] this The reaction object (modified: eq_cst is updated)
        !> \param[in] temp Temperature [K]
        !> @brief Compute temperature-dependent equilibrium constant
        !> @details Calculates K_eq at specified temperature using PHREEQC analytical expression:
        !> log₁₀(K) = a₁ + a₂T + a₃/T + a₄log₁₀(T) + a₅/T² + a₆T²
        !> Updates this%eq_cst with computed K = 10^(logK).
        !> @param[in,out] this Reaction object (eq_cst is updated)
        !> @param[in] temp Temperature of solution [K]
        subroutine compute_logK_dep_T(this,temp)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(inout) :: this !< Reaction object (eq_cst modified) [-]
            real(kind=8), intent(in) :: temp !< Temperature for K_eq calculation [K]
            
            real(kind=8) :: logK !< Logarithm base 10 of equilibrium constant [-]
            
            !> Validate that exactly 6 temperature coefficients are available
            if (size(this%coeffs_logK_T)/=6) error stop "There must be 6 coefficients to compute logK(T)" !< Halt if wrong array size
            !> Compute log₁₀(K) using PHREEQC analytical expression with 6 terms
            logK=this%coeffs_logK_T(1)+this%coeffs_logK_T(2)*temp+this%coeffs_logK_T(3)/temp+this%coeffs_logK_T(4)*log10(temp)+&
                this%coeffs_logK_T(5)/(temp**2)+this%coeffs_logK_T(6)*(temp**2) !< Calculate logK from polynomial [-]
            this%eq_cst=10**logK !< Convert from log₁₀(K) to K via antilog [-]
        end subroutine
    !******************************************************************************************************************************
        !> \brief Allocate temperature coefficient array
        !> \details
        !>   Allocates the 6-element array for temperature dependence coefficients.
        !>   Based on PHREEQC database format.
        !> \param[in,out] this The reaction object
        !> @brief Allocate temperature coefficient array
        !> @details Allocates 6-element array for PHREEQC temperature dependence coefficients.
        !> Based on PHREEQC database format for analytical log(K) vs T expressions.
        !> @param[in,out] this Reaction object to allocate
        subroutine allocate_coeffs_logK_T(this)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(inout) :: this !< Reaction object to allocate coefficients array [-]
            
            allocate(this%coeffs_logK_T(6)) !< Allocate 6-element array for temperature dependence coefficients [-]
        end subroutine
    !******************************************************************************************************************************
        !> \brief Reverse reaction direction
        !> \details
        !>   Changes the sign of all stoichiometric coefficients.
        !>   Effectively reverses the reaction: products become reactants and vice versa.
        !>   
        !>   Example:
        !>   Forward:  A + B → C  (stoich = [-1, -1, 1])
        !>   Reverse:  C → A + B  (stoich = [1, 1, -1])
        !> \param[in,out] this The reaction object
        !> @brief Reverse reaction direction
        !> @details Changes sign of all stoichiometric coefficients to reverse reaction.
        !> Products become reactants and vice versa.
        !> Example: A + B → C ([-1,-1,1]) becomes C → A + B ([1,1,-1])
        !> @param[in,out] this Reaction object to reverse
        subroutine change_sign_stoichiometry(this)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(inout) :: this !< Reaction object to reverse direction [-]
            integer(kind=4) :: i !< Loop counter for iterating over species [-]
            
            !> Iterate through all species and negate their stoichiometric coefficients
            do i=1,this%num_species !< Loop through all species in reaction
                this%stoichiometry(i)=-this%stoichiometry(i) !< Negate coefficient to reverse direction [-]
            end do
        end subroutine
    !******************************************************************************************************************************
        !> \brief Write reaction to file
        !> \details
        !>   Writes reaction information to specified file unit:
        !>   - Reaction name
        !>   - Reaction type code
        !>   - Number of species
        !>   - Equilibrium constant
        !>   
        !>   Format suitable for reading back or debugging.
        !> \param[in,out] this The reaction object
        !> \param[in] unit File unit number (must be open for writing)
        !> @brief Write reaction to file
        !> @details Writes reaction information to specified file unit in formatted output.
        !> Includes: name (30 chars), type code, number of species, equilibrium constant.
        !> File unit must be open for writing before calling.
        !> @param[in] this Reaction object to write
        !> @param[in] unit File unit number (must be open for writing) [-]
        subroutine write_reaction_sup(this,unit)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(in) :: this !< Reaction object to write to file [-]
            integer(kind=4), intent(in) :: unit !< File unit number for output (must be open) [-]
            
            write(unit,"(10x,A30/)") this%name !< Write reaction name (30 character field) [-]
            ! write(unit,"(10x,I5/)") this%react_type !< Write reaction type code (integer) [-]
            ! write(unit,"(10x,I5/)") this%num_species !< Write number of species (integer) [-]
            ! write(unit,"(10x,ES15.5/)") this%eq_cst !< Write equilibrium constant (scientific notation) [-]
        end subroutine
        
        subroutine rearrange_Species_indices(this,old_Species,new_Species)
            implicit none !< Enforce explicit variable declarations
            class(reaction_c), intent(inout) :: this !< Reaction object to configure [-]
            class(species_c), intent(in) :: old_Species(:) !< Array of old species in chemical system [-]
            class(species_c), intent(in) :: new_Species(:) !< Array of new species in chemical system [-]
            
            integer(kind=4) :: i,j !< Loop counters [-]
            integer(kind=4), allocatable :: old_species_ind(:) !< Temporary copy of old species indices [-]
            logical :: found !< Flag indicating if species was found [-]
            
            if (size(old_Species)/=size(new_species)) then
                error stop "Dimension error in 'rearrange_Species_indices'" !< Halt if size mismatch detected
            end if
            !> Save old indices before remapping
            old_species_ind=this%species_ind
            !> For each species in the reaction, find its new index in the rearranged species array
            do i=1,this%num_species
                found=.false.
                do j=1,size(new_Species)
                    if (old_Species(old_species_ind(i))%name==new_Species(j)%name) then
                        this%species_ind(i)=j !< Update to new index [-]
                        found=.true.
                        exit
                    end if
                end do
                if (.not. found) then
                    error stop "Species not found in new species array in 'rearrange_Species_indices'"
                end if
            end do
        end subroutine

end module