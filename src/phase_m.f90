!> @file phase_m.f90
!> @brief Abstract phase base class module for multiphase reactive transport
!> @details This module defines the abstract base class phase_c that represents a thermodynamic
!> phase (aqueous, gas, solid) containing multiple chemical species. It provides the fundamental
!> properties and operations common to all phases.
!>
!> @par Phase Concept:
!> A phase is a physically distinct and chemically homogeneous region of matter:
!> - **Aqueous phase**: Dissolved species in water (ions, complexes)
!> - **Gas phase**: Gaseous species (CO₂, O₂, CH₄)
!> - **Solid phase**: Mineral species (calcite, quartz, clays)
!>
!> @par Activity Classification:
!> Species in a phase are classified by activity behavior:
!> - **Variable activity**: Activity depends on concentration (most ions, dissolved gases)
!> - **Constant activity**: Activity fixed at 1 (pure minerals, water in dilute solutions)
!>
!> @par Species Count Relationship:
!> num_species = num_var_act_species + num_cst_act_species
!>
!> @par Usage:
!> This is an abstract class - concrete implementations include:
!> - aq_phase_c (aqueous phase)
!> - gas_phase_c (gas phase)
!> - mineral_c (mineral phase)
!> - biofilm_c (biofilm phase)
!> - surface_c (surface phase)
!>
!> @see aq_phase_m For aqueous phase implementation
!> @see gas_phase_m For gas phase implementation
!> @author Jordi Petchamé-Guerrero
!> @date November 2025

module phase_m
    implicit none !< Enforce explicit variable declarations for type safety
    save !< Preserve module variables between procedure calls
    private !< Default visibility for module entities is private
    !> @brief Abstract base class for thermodynamic phases
    !> @details This type serves as the foundation for all phase implementations.
    !> It tracks the number and types of species, allowing phases to organize species
    !> by activity behavior for efficient speciation calculations.
    !>
    !> @par Design Pattern:
    !> Abstract base class (cannot be instantiated directly)
    !> Concrete subclasses must implement phase-specific behavior
    !>
    !> @par Species Organization:
    !> Species are typically arranged:
    !> 1. Variable activity species first (require activity coefficient updates)
    !> 2. Constant activity species last (activity = 1, no updates needed)
    !>
    !> @par Invariant:
    !> num_species = num_var_act_species + num_cst_act_species
    !> This relationship must always hold true.
    type, public, abstract :: phase_c !< Abstract phase base class (cannot be instantiated)
        integer(kind=4) :: num_species=0 !< Total number of species in phase [-] (default: 0)
        integer(kind=4) :: num_var_act_species=0 !< Number of variable activity species [-] (default: 0, activity depends on concentration)
        integer(kind=4) :: num_cst_act_species=0 !< Number of constant activity species [-] (default: 0, activity = 1)
        character(len=256)  :: name !< Phase name identifier (e.g., 'aqueous', 'gas', 'calcite') [string]
    contains
    !> @name Setter Methods
    !> @{
        procedure :: set_num_species_phase !< Set total number of species with validation
        procedure :: set_num_var_act_species_phase !< Set number of variable activity species with validation
        procedure :: set_num_cst_act_species_phase !< Set number of constant activity species with validation
        procedure :: set_phase_name !< Set phase name identifier
    !> @}
    !> @name Computation Methods
    !> @{
        procedure :: compute_num_species_phase !< Compute total species from variable + constant activity counts
    !> @}
    end type
    
    contains
        !> @brief Set phase name identifier
        !> @details Assigns a descriptive name to the phase for identification and output.
        !> The name is typically the phase type (e.g., 'aqueous', 'gas') or mineral name.
        !>
        !> @par Example Names:
        !> - Aqueous phase: 'aqueous', 'solution', 'water'
        !> - Gas phase: 'gas', 'vapor'
        !> - Solid phases: 'calcite', 'quartz', 'kaolinite'
        !>
        !> @param[inout] this Phase object
        !> @param[in] name Phase name (up to 256 characters) [string]
        subroutine set_phase_name(this,name)
            implicit none !< Enforce explicit variable declarations
            class(phase_c), intent(inout) :: this !< Phase object [-]
            character(len=*), intent(in) :: name !< Phase name identifier [string]
            this%name=name !< Assign phase name [string]
        end subroutine 
        
        !> @brief Set total number of species in phase
        !> @details Assigns the total count of species present in the phase.
        !> Performs validation to ensure the count is non-negative.
        !>
        !> @par Validation:
        !> - num_species >= 0 (terminates if negative)
        !>
        !> @par Relationship:
        !> After setting, should satisfy: num_species = num_var_act_species + num_cst_act_species
        !>
        !> @param[inout] this Phase object
        !> @param[in] num_species Total number of species [-]
        subroutine set_num_species_phase(this,num_species)
            implicit none !< Enforce explicit variable declarations
            class(phase_c), intent(inout) :: this !< Phase object [-]
            integer(kind=4), intent(in) :: num_species !< Total number of species [-]
            
            !> Validate that species count is non-negative
            if (num_species<0) error stop "Number of species cannot be negative" !< Terminate if negative
            this%num_species=num_species !< Assign total species count [-]
        end subroutine 
        
        !> @brief Compute total species count from activity categories
        !> @details Calculates total number of species by summing variable and constant
        !> activity species counts. This enforces the fundamental relationship:
        !> num_species = num_var_act_species + num_cst_act_species
        !>
        !> @par Use Case:
        !> Called after setting num_var_act_species and num_cst_act_species to automatically
        !> update the total count, ensuring consistency.
        !>
        !> @par Example:
        !> If phase has 8 variable activity species and 2 constant activity species,
        !> this will set num_species = 10.
        !>
        !> @param[inout] this Phase object
        subroutine compute_num_species_phase(this)
            implicit none !< Enforce explicit variable declarations
            class(phase_c), intent(inout) :: this !< Phase object [-]
            
            !> Compute total species as sum of variable and constant activity species
            this%num_species=this%num_var_act_species+this%num_cst_act_species !< Total = variable + constant [-]
        end subroutine 
        
        !> @brief Set number of variable activity species
        !> @details Assigns the count of species whose activity depends on concentration.
        !> These species require activity coefficient calculations in speciation.
        !>
        !> @par Variable Activity Species Examples:
        !> - Aqueous phase: Most ions (Na⁺, Cl⁻, Ca²⁺), dissolved gases
        !> - Gas phase: All gas species at variable pressures
        !> - Solid phase: Solid solutions with variable composition
        !>
        !> @par Validation:
        !> - num_var_act_species >= 0 (non-negative)
        !> - num_var_act_species <= num_species (cannot exceed total)
        !>
        !> @param[inout] this Phase object
        !> @param[in] num_var_act_species Number of variable activity species [-]
        subroutine set_num_var_act_species_phase(this,num_var_act_species)
            implicit none !< Enforce explicit variable declarations
            class(phase_c), intent(inout) :: this !< Phase object [-]
            integer(kind=4), intent(in) :: num_var_act_species !< Number of variable activity species [-]
            
            !> Validate non-negative count
            if (num_var_act_species<0) then
                error stop "Number of variable activity species cannot be negative" !< Terminate if negative
            !> Validate does not exceed total species
            else if (num_var_act_species>this%num_species) then
                error stop "Number of variable activity species cannot be greater than number of species" !< Terminate if too large
            !> All validation passed, assign value
            else
                this%num_var_act_species=num_var_act_species !< Assign variable activity species count [-]
            end if
        end subroutine 
        
        !> @brief Set number of constant activity species
        !> @details Assigns the count of species with fixed activity (typically activity = 1).
        !> These species do not require activity coefficient updates, improving computational efficiency.
        !>
        !> @par Constant Activity Species Examples:
        !> - Aqueous phase: Water in dilute solutions (aH₂O = 1)
        !> - Gas phase: (rare, usually all gases have variable activity)
        !> - Solid phase: Pure minerals (acalcite = 1, aquartz = 1)
        !>
        !> @par Thermodynamic Justification:
        !> For pure phases or components in large excess, activity is conventionally set to 1.
        !> This simplifies equilibrium calculations without significant loss of accuracy.
        !>
        !> @par Validation:
        !> - num_cst_act_species >= 0 (non-negative)
        !> - num_cst_act_species <= num_species (cannot exceed total)
        !>
        !> @param[inout] this Phase object
        !> @param[in] num_cst_act_species Number of constant activity species [-]
        subroutine set_num_cst_act_species_phase(this,num_cst_act_species)
            implicit none !< Enforce explicit variable declarations
            class(phase_c), intent(inout) :: this !< Phase object [-]
            integer(kind=4), intent(in) :: num_cst_act_species !< Number of constant activity species [-]
            
            !> Validate non-negative count
            if (num_cst_act_species<0) then
                error stop "Number of species cannot be negative" !< Terminate if negative (error message could be more specific)
            !> Validate does not exceed total species
            else if (num_cst_act_species>this%num_species) then
                error stop "Number of constant activity species cannot be greater than number of species" !< Terminate if too large
            !> All validation passed, assign value
            else
                this%num_cst_act_species=num_cst_act_species !< Assign constant activity species count [-]
            end if
        end subroutine 
end module