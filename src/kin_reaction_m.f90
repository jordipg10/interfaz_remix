!> \file kin_reaction_m.f90
!> \brief Defines abstract and concrete types for kinetic reactions in reactive transport
!>
!> \details
!> This module provides the foundational infrastructure for kinetic reactions in
!> geochemical reactive transport simulations. It defines abstract base classes and
!> container types that are extended by specific kinetic reaction implementations
!> (mineral dissolution/precipitation, Monod reactions, redox reactions, etc.).
!>
!> **Key Components:**
!> 1. **kin_reaction_c (abstract):** Base class for all kinetic reactions
!>    - Stores aqueous species that control reaction rates (e.g., H⁺, O₂, nutrients)
!>    - Provides interface for parameter I/O
!>    - Extends the base reaction_c class
!>
!> 2. **kin_reaction_poly_c:** Polymorphic container for kinetic reaction pointers
!>    - Workaround for Fortran's limited pointer array support
!>    - Allows arrays of heterogeneous kinetic reaction types
!>
!> **Kinetic Reaction Rate Dependence:**
!> Kinetic reactions typically depend on concentrations/activities of aqueous species:
!> \f[
!> r = f(c_1, c_2, \ldots, c_n)
!> \f]
!> where \f$ c_i \f$ are aqueous species concentrations that control the rate.
!>
!> Examples:
!> - Mineral dissolution: Rate depends on pH (H⁺ concentration)
!> - Monod reactions: Rate depends on electron acceptor and donor concentrations
!> - Redox reactions: Rate depends on oxidant and reductant species
!>
!> **Inheritance Hierarchy:**
!> ```
!> reaction_c (base)
!>   └─ kin_reaction_c (abstract - this module)
!>       ├─ lin_kin_reaction_c (linear kinetics)
!>       ├─ kin_mineral_c (mineral reactions)
!>       └─ redox_kin_reaction_c (redox reactions)
!> ```
!>
!> \note This is an abstract module - concrete implementations must define write_params
!> \note The indices_aq_phase array maps to positions in the aqueous phase species array
!>
!> \see reaction_m For base reaction class
!> \see kin_params_m For kinetic parameter classes
!> \see kin_mineral_m For mineral kinetic reactions
!> \see redox_kin_reaction_m For redox kinetic reactions
!>
!> \author Generated documentation
!> \date 2024

module kin_reaction_m
    use reaction_m, only: reaction_c !> Import base reaction class and write procedure
        !> reaction_c: Base class for all reactions (equilibrium and kinetic)
        !> write_reaction_sup: Procedure to write reaction stoichiometry and equilibrium constant
    use kin_params_m, only: kin_params_c !> Import base kinetic parameters class
        !> kin_params_c: Contains kinetic rate constants and reaction-specific parameters
    implicit none !> Enforce explicit variable declaration (no implicit typing)
    save !> Save module variables between procedure calls (persistent state)
    private !> Private module scope: internal details hidden from outside modules
!> ================================================================================
!> TYPE DEFINITION: kin_reaction_c (Abstract Kinetic Reaction)
!> ================================================================================
    
    !> \brief Abstract base class for kinetic reactions
    !>
    !> \details
    !> This abstract class extends reaction_c to add kinetic-specific functionality.
    !> All kinetic reactions (minerals, redox, biodegradation) inherit from this class.
    !>
    !> **Key Features:**
    !> - Tracks aqueous species that control reaction rate (e.g., catalysts, inhibitors)
    !> - Provides allocation procedures for species indices
    !> - Defines abstract interface for parameter I/O
    !>
    !> **Rate-Controlling Species:**
    !> Many kinetic reactions depend on specific aqueous species:
    !> - Mineral dissolution: H⁺, OH⁻, organic ligands
    !> - Biodegradation: O₂, NO₃⁻, SO₄²⁻ (electron acceptors)
    !> - Nitrification: NH₄⁺, O₂
    !>
    !> The `indices_aq_phase` array stores the indices of these species in the
    !> aqueous phase species array, allowing efficient rate calculations.
    !>
    !> \extends reaction_c
    
    type, public, abstract, extends(reaction_c) :: kin_reaction_c !> Kinetic reaction abstract subclass
        !> Extends reaction_c (which contains stoichiometry, species, equilibrium constant)
        !> with kinetic-specific members
        
        !> ----------------------------------------------------------------
        !> Rate-controlling aqueous species
        !> ----------------------------------------------------------------
        integer(kind=4) :: num_aq_rk=0 !> Number of aqueous species relevant for kinetic reaction rate [-]
            !> Number of aqueous species whose concentrations/activities affect the reaction rate
            !> Default: 0 (no rate-controlling species defined)
            !> Examples:
            !>   - Acid-catalyzed dissolution: num_aq_rk = 1 (H⁺)
            !>   - Monod reaction: num_aq_rk = 2 (e.g., glucose + O₂)
            !>   - Complex catalysis: num_aq_rk > 2
        
        integer(kind=4), allocatable :: indices_aq_phase(:) !> Indices of rate-controlling species in aqueous phase [-]
            !> Array dimension: (num_aq_rk)
            !> Maps rate-controlling species to their position in aqueous phase species array
            !> indices_aq_phase(i) = index of i-th rate-controlling species in aq_phase%aq_species(:)
            !> Example: If H⁺ is 3rd species in aqueous phase, indices_aq_phase(1) = 3
        integer(kind=4), allocatable :: indices_react_species(:) !> Reaction-species indices for each aqueous species [-]
            !> Array dimension: (num_aq_rk)
            !> indices_react_species(i) = position in stoichiometry/species_ind for the i-th aqueous species
            !> Needed when not all species in the reaction are aqueous (e.g., H2O with constant activity)
    
    contains
        !> ----------------------------------------------------------------
        !> Public procedures
        !> ----------------------------------------------------------------
        procedure :: set_num_aq_rk !> Set number of aqueous species relevant for kinetic rate
        procedure :: allocate_indices_aq_phase_kin_react !> Allocate indices array for aqueous phase species
        procedure(write_params), public, deferred :: write_params !> Abstract deferred procedure to write kinetic parameters
            !> Must be implemented by concrete subclasses (e.g., kin_mineral_c, redox_kin_reaction_c)
            !> Writes kinetic parameters (rate constants, activation energy, etc.) to file
    end type kin_reaction_c
    
!> ================================================================================
!> TYPE DEFINITION: kin_reaction_poly_c (Polymorphic Container)
!> ================================================================================
    
    !> \brief Polymorphic container for kinetic reaction pointers
    !>
    !> \details
    !> This type provides a workaround for Fortran's limitation with arrays of
    !> polymorphic pointers. It allows creation of arrays containing pointers to
    !> different kinetic reaction subclasses.
    !>
    !> **Use Case:**
    !> ```fortran
    !> type(kin_reaction_poly_c), allocatable :: kin_reactions(:)
    !> allocate(kin_reactions(n))
    !> ! Can point to different subclasses:
    !> kin_reactions(1)%kin_reaction => my_mineral_reaction
    !> kin_reactions(2)%kin_reaction => my_redox_reaction
    !> ```
    !>
    !> **Why Needed:**
    !> Fortran does not directly support arrays like:
    !> ```fortran
    !> class(kin_reaction_c), pointer :: reactions(:)  ! Not allowed
    !> ```
    !> Instead, we use an array of containers, each holding a pointer.
    !>
    !> \note This is a common Fortran design pattern for polymorphic arrays
    
    type, public :: kin_reaction_poly_c !> Polymorphic container for kinetic reaction pointers
        !> Ad hoc class to create vector of kinetic reaction pointers (Fortran workaround)
        
        class(kin_reaction_c), pointer :: kin_reaction => null() !> Pointer to kinetic reaction object
            !> Can point to any concrete subclass of kin_reaction_c
            !> Examples: kin_mineral_c, redox_kin_reaction_c, lin_kin_reaction_c
    
    contains
        procedure :: set_kin_reaction !> Set the kinetic reaction pointer in the container
    end type kin_reaction_poly_c
    
!> ================================================================================
!> ABSTRACT INTERFACE: write_params
!> ================================================================================
    
    !> \brief Abstract interface for writing kinetic reaction parameters to file
    !>
    !> \details
    !> This abstract interface must be implemented by all concrete kinetic reaction
    !> subclasses. It defines the signature for writing kinetic parameters (rate
    !> constants, activation energies, experimental coefficients, etc.) to output files.
    !>
    !> **Implementation Requirements:**
    !> Each concrete subclass must provide a write_params implementation that writes:
    !> - Kinetic rate constants
    !> - Temperature dependence parameters (activation energy)
    !> - Catalytic coefficients
    !> - Any other reaction-specific parameters
    !>
    !> **Example Implementations:**
    !> - kin_mineral_c: Writes rate constants, activation energy, catalytic exponents
    !> - redox_kin_reaction_c: Writes Monod constants, inhibition parameters
    !> - lin_kin_reaction_c: Writes linear rate constant
    
    abstract interface
        !> \brief Write kinetic reaction parameters to file
        !>
        !> \details
        !> Abstract procedure that must be implemented by concrete subclasses.
        !> Writes kinetic parameters to the specified file unit in a format
        !> suitable for input/output and post-processing.
        !>
        !> \param[in] this Kinetic reaction object (intent: in)
        !> \param[in] unit File unit number for writing (Fortran I/O unit)
        !>
        !> \pre File must be opened for writing on the specified unit
        !> \post Kinetic parameters written to file in implementation-specific format
        
        subroutine write_params(this,unit)
            import kin_reaction_c !> Import kin_reaction_c type into interface scope
            implicit none !> Enforce explicit variable declaration
            class(kin_reaction_c), intent(in) :: this !> Kinetic reaction object (intent: in implied)
            integer(kind=4), intent(in) :: unit !> File unit number for I/O
        end subroutine
        
    end interface
    
!> ================================================================================
!> INTERFACE DECLARATIONS (Reserved)
!> ================================================================================
    
    interface 
        !> Reserved for future extensions, operator overloads, or generic interfaces
        !> Currently empty - may be used for:
        !> - Generic constructors
        !> - Operator overloading (e.g., ==, +)
        !> - Type-bound generic procedures
    end interface
    
!> ================================================================================
!> MODULE PROCEDURES
!> ================================================================================
    contains
        
!> --------------------------------------------------------------------------------
!> Subroutine: set_kin_reaction
!> --------------------------------------------------------------------------------
!> \brief Set the kinetic reaction pointer in the polymorphic container
!>
!> \details
!> Associates the container's kin_reaction pointer with a target kinetic reaction
!> object. This allows the container to hold a reference to any kinetic reaction
!> subclass.
!>
!> **Usage Pattern:**
!> ```fortran
!> type(kin_reaction_poly_c) :: container
!> class(kin_mineral_c), target :: my_mineral
!> call container%set_kin_reaction(my_mineral)
!> ! Now container%kin_reaction points to my_mineral
!> ```
!>
!> \param[in,out] this Polymorphic container object
!> \param[in] kin_reaction Kinetic reaction object to assign (must be target)
!>
!> \note The kin_reaction argument must have the TARGET attribute
!> \note The pointer association persists until reassigned or nullified
        
        subroutine set_kin_reaction(this,kin_reaction)
            class(kin_reaction_poly_c) :: this !> Polymorphic container (intent: inout)
            class(kin_reaction_c), intent(in), target :: kin_reaction !> Kinetic reaction to assign
                !> Must be declared with TARGET attribute to allow pointer association
            
            this%kin_reaction=>kin_reaction !> Associate pointer with target kinetic reaction
                !> Establishes pointer connection using Fortran pointer assignment (=>)
        end subroutine

!> --------------------------------------------------------------------------------
!> Subroutine: set_num_aq_rk
!> --------------------------------------------------------------------------------
!> \brief Set the number of aqueous species relevant for kinetic reaction rate
!>
!> \details
!> Sets the number of aqueous species whose concentrations/activities affect the
!> kinetic reaction rate. Validates that the number is non-negative.
!>
!> **Examples:**
!> - Mineral dissolution with H⁺ catalysis: num_aq_rk = 1
!> - Aerobic respiration (glucose + O₂): num_aq_rk = 2
!> - Multi-catalyzed reaction: num_aq_rk = 3+
!>
!> \param[in,out] this Kinetic reaction object
!> \param[in] num_aq_rk Number of aqueous species (must be ≥ 0)
!>
!> \pre num_aq_rk ≥ 0
!> \post this%num_aq_rk is set to the input value
!>
!> \throws error stop if num_aq_rk < 0
        
        subroutine set_num_aq_rk(this,num_aq_rk)
            class(kin_reaction_c) :: this !> Kinetic reaction object (intent: inout)
            integer(kind=4), intent(in) :: num_aq_rk !> Number of aqueous species
                !> Must be non-negative (0 allowed for rate-independent reactions)
            !if (prese
                if (num_aq_rk<0) then !> Validate input is non-negative
                    error stop "Number of aqueous species relevant for kinetic reaction rate must be positive"
                        !> Terminate program with error message if validation fails
                else !> Input is valid
                    this%num_aq_rk=num_aq_rk !> Assign value to object member
                        !> Store number of rate-controlling aqueous species
                end if
        end subroutine

!> --------------------------------------------------------------------------------
!> Subroutine: allocate_indices_aq_phase_kin_react
!> --------------------------------------------------------------------------------
!> \brief Allocate indices array for aqueous phase species controlling kinetic rates
!>
!> \details
!> Allocates the indices_aq_phase array which stores the positions of rate-controlling
!> aqueous species in the aqueous phase species array. Optionally sets num_aq_rk
!> before allocation.
!>
!> **Two Usage Patterns:**
!>
!> 1. With optional argument (sets num_aq_rk and allocates):
!> ```fortran
!> call kin_react%allocate_indices_aq_phase_kin_react(num_aq_rk=2)
!> ```
!>
!> 2. Without optional argument (num_aq_rk must be set beforehand):
!> ```fortran
!> call kin_react%set_num_aq_rk(2)
!> call kin_react%allocate_indices_aq_phase_kin_react()
!> ```
!>
!> \param[in,out] this Kinetic reaction object
!> \param[in] num_aq_rk [optional] Number of aqueous species (if provided, sets num_aq_rk)
!>
!> \pre If num_aq_rk not provided, this%num_aq_rk must be already set
!> \post indices_aq_phase array allocated with size this%num_aq_rk
!>
!> \note After allocation, indices must be filled with actual species positions
        
        subroutine allocate_indices_aq_phase_kin_react(this,num_aq_rk)
            class(kin_reaction_c) :: this !> Kinetic reaction object (intent: inout)
            integer(kind=4), intent(in), optional :: num_aq_rk !> [optional] Number of aqueous species
                !> If present, calls set_num_aq_rk before allocation
            
            if (present(num_aq_rk)) then !> Check if optional argument provided
                call this%set_num_aq_rk(num_aq_rk) !> Set num_aq_rk using validation procedure
                    !> Validates input is non-negative and assigns to this%num_aq_rk
            end if !> End optional argument check
            
            allocate(this%indices_aq_phase(this%num_aq_rk)) !> Allocate indices array
                !> Array dimension: (num_aq_rk)
                !> Will store integer indices mapping to aq_phase%aq_species(:)
                !> After allocation, array contains uninitialized values - must be filled
            allocate(this%indices_react_species(this%num_aq_rk)) !> Allocate parallel stoichiometry index array
        end subroutine
            
end module kin_reaction_m