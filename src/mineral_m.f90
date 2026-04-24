!> \file mineral_m.f90
!> \brief Mineral phase module
!> \details
!>   Defines mineral phase as extension of phase class.
!>   A mineral phase consists of a single solid (mineral) species.
!>   Inherits thermodynamic and chemical properties from phase_c.
!>   Used in reactive transport for mineral dissolution/precipitation reactions.
!>
!> \author jordi Petchamé-Guerrero
!> \date October 2025

module mineral_m
    use phase_m, only: phase_c                                                !< Base phase class
    use solid_m, only: solid_species_c                                                !< Solid species class
    implicit none
    save
    private
    !> \brief Mineral phase class
    !> \details
    !>   - Extends phase_c to represent pure mineral phases
    !>   - Contains single solid species (mineral) that defines the phase
    !>   - Used for mineral dissolution/precipitation reactions
    !>   - Thermodynamic properties inherited from phase_c
    !>   - Chemical composition defined by mineral species
    type, public, extends(phase_c) :: mineral_c
        type(solid_species_c) :: mineral                                              !< Mineral species defining this phase
    contains
        !> No type-bound procedures currently defined
    end type
    
!> \section pflotran_reference PFLOTRAN Reference
!> \details
!>   The following shows the PFLOTRAN mineral reaction type structure for reference.
!>   This is commented out as it's not currently used in this implementation.
!>   
!>   Mineral reaction structure would include:
!>   - id: Unique identifier for the mineral reaction
!>   - itype: Type of mineral reaction
!>   - name: Name of the mineral
!>   - molar_volume: Molar volume of the mineral (m³/mol)
!>   - molar_weight: Molar weight of the mineral (kg/mol)
!>   - print_me: Flag to control output
!>   - dbaserxn: Pointer to database reaction data
!>   - tstrxn: Pointer to transition state reaction data
!>   - next: Pointer to next mineral reaction (linked list)

!> PFLOTRAN:
    
      !type, public :: mineral_rxn_type
      !>  PetscInt :: id
      !>  PetscInt :: itype
      !>  character(len=MAXWORDLENGTH) :: name
      !>  PetscReal :: molar_volume
      !>  PetscReal :: molar_weight
      !>  PetscBool :: print_me
      !>  type(database_rxn_type), pointer :: dbaserxn
      !>  type(transition_state_rxn_type), pointer :: tstrxn
      !>  type(mineral_rxn_type), pointer :: next
      !end type mineral_rxn_type
  !>  

!*************************************************************************************************!>    
    !> \brief Interface declarations (currently empty)
    interface
      
    end interface
    
    
    contains
    
        !> \brief Set molar volume of mineral (COMMENTED OUT)
        !> \details
        !>   This subroutine would set the molar volume of a mineral phase.
        !>   Includes validation to ensure molar volume is non-negative.
        !>   Currently not implemented - commented out for future use.
        !>
        !> \param[inout] this Mineral object
        !> \param[in] mol_vol Molar volume to set (must be >= 0)
        !>
        !> Example implementation:
        !> \code{.f90}
        !> if (mol_vol < 0) error stop "Molar volume cannot be negative"
        !> this%mol_vol = mol_vol
        !> \endcode
        
        !subroutine set_mol_vol(this,mol_vol)
        !>    implicit none
        !>    class(mineral_c) :: this
        !>    real(kind=8), intent(in) :: mol_vol
        !>    if (mol_vol<0d0) error stop "Molar volume cannot be negative"
        !>    this%mol_vol=mol_vol
        !end subroutine
        
    !> \brief Compare two mineral sets for equality (COMMENTED OUT)
    !> \details
    !>   This subroutine would compare two mineral objects to determine if they are equal.
    !>   Equality criteria would include:
    !>   - Same number of mineral species
    !>   - Same number of equilibrium minerals
    !>   - Same constant activity species
    !>   - All minerals from first set present in second set
    !>   
    !>   Note: Current implementation has bugs (references gas_phase instead of minerals).
    !>   Needs correction before activation.
    !>
    !> \param[in] min_1 First mineral object to compare
    !> \param[in] min_2 Second mineral object to compare
    !> \param[out] flag TRUE if minerals are equal, FALSE otherwise
    !>
    !> Algorithm:
    !> 1. Check if number of species match
    !> 2. Check if number of equilibrium species match
    !> 3. Check if constant activity species match
    !> 4. Verify each mineral from min_1 exists in min_2
    
    !subroutine are_minerals_equal(min_1,min_2,flag)
    !    class(mineral_c), intent(in) :: min_1 !> first set of minerals
    !    class(mineral_c), intent(in) :: min_2 !> second set of minerals
    !    logical, intent(out) :: flag !> TRUE if minerals are equal, FALSE otherwise
    !    
    !    integer(kind=4) :: i
    !    logical :: min_flag !> mineral flag
    !
    !    flag=.true. !> set flag to true by default
    !    if (min_1%/=gas_phase_2%num_species .or. gas_phase_1%num_gases_eq/=gas_phase_2%num_gases_eq) then
    !        flag=.false.
    !    else if (gas_phase_1%num_cst_act_species/=gas_phase_2%num_cst_act_species .or. gas_phase_1%num_gases_eq_cst_act/=gas_phase_2%num_gases_eq_cst_act) then
    !        flag=.false.
    !    else
    !        do i=1,gas_phase_1%num_species
    !            call gas_phase_2%is_gas_in_gas_phase(gas_phase_1%gases(i),min_flag)
    !            if (min_flag .eqv. .false.) then
    !                flag=.false.
    !                exit
    !            end if
    !        end do
    !    end if
    !end subroutine
end module