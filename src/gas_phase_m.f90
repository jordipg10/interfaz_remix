!> \file gas_phase_m.f90
!> \brief Gas phase module for gaseous species management
!> \details
!>   This module defines the gas phase class for managing gases in equilibrium
!>   and kinetic reactions with constant or variable activity.
!>   
!>   **Key features:**
!>   - Equilibrium gases (constant/variable activity)
!>   - Kinetic gases (constant/variable activity)
!>   - Gas ordering conventions
!>   - Activity coefficient management
!>   
!>   **Gas ordering in gases array:**
!>   1. Equilibrium gases with constant activity
!>   2. Equilibrium gases with variable activity
!>   3. Kinetic gases with constant activity
!>   4. Kinetic gases with variable activity
!>   
!>   **Applications:**
!>   - Gas dissolution/exsolution reactions
!>   - CO2 degassing
!>   - Methane production/consumption
!>   - Oxygen transport
!>   
!> \author Jordi Petchamé-Guerrero
!> \date October 2025
module gas_phase_m
    use phase_m, only: phase_c !< Import generic phase class
    use gas_species_m, only: gas_species_c !< Import gas species class
    implicit none !< Enforce explicit variable declarations
    save !< Preserve module state between calls
    private !< Private module scope: internal details hidden from outside modules
    public :: are_gas_phases_equal !< Publicly expose equality check function for gas phases
    !> \brief Gas phase type for gaseous species management
    !> \details
    !>   Defines a gas phase containing gases in equilibrium and kinetic reactions.
    !>   Extends the generic phase_c class.
    !>   
    !>   **Gas ordering convention in gases array:**
    !>   1. Equilibrium gases with constant activity (buffered partial pressure)
    !>   2. Equilibrium gases with variable activity (free to evolve)
    !>   3. Kinetic gases with constant activity
    !>   4. Kinetic gases with variable activity
    !>   
    !>   **Constant activity gases:**
    !>   - Fixed partial pressure (e.g., atmospheric O2)
    !>   - Act as infinite reservoir/sink
    !>   
    !>   **Variable activity gases:**
    !>   - Partial pressure changes with reaction progress
    !>   - Finite amount in system
    !>   
    !>   **Usage example:**
    !>   ```fortran
    !>   type(gas_phase_c) :: gas_phase
    !>   call gas_phase%set_num_gases_eq(2)      ! CO2, O2
    !>   call gas_phase%set_num_gases_eq_cst_act(1)  ! O2 constant
    !>   call gas_phase%set_num_gases_eq_var_act(1)  ! CO2 variable
    !>   call gas_phase%allocate_gases(2)
    !>   ```
    type, public, extends(phase_c) :: gas_phase_c !< Gas phase (subclass of phase class)
        !> \brief Gas species array with specific ordering
        !> \details
        !>   Ordering: (1) equilibrium + constant activity, (2) equilibrium + variable activity, (3) kinetic + constant activity, (4) kinetic + variable activity
        type(gas_species_c), allocatable :: gases(:) !< Gases array with ordering convention [-]
        integer(kind=4) :: num_gases_eq=0 !< Number of gases in equilibrium reactions [-]
        integer(kind=4) :: num_gases_kin=0 !< Number of gases in kinetic reactions [-]
        integer(kind=4) :: num_gases_eq_cst_act=0 !< Number of equilibrium gases with constant activity [-]
        integer(kind=4) :: num_gases_eq_var_act=0 !< Number of equilibrium gases with variable activity [-]
        integer(kind=4) :: num_gases_kin_cst_act=0 !< Number of kinetic gases with constant activity [-]
        integer(kind=4) :: num_gases_kin_var_act=0 !< Number of kinetic gases with variable activity [-]
    contains
    !> Set procedures
        procedure :: set_num_gases_eq !< Set number of equilibrium gases
        procedure :: set_num_gases_eq_cst_act !< Set number of equilibrium gases with constant activity
        procedure :: set_num_gases_eq_var_act !< Set number of equilibrium gases with variable activity
        procedure :: set_num_gases_kin_cst_act !< Set number of kinetic gases with constant activity
        procedure :: set_num_gases_kin_var_act !< Set number of kinetic gases with variable activity
        procedure :: set_num_gases_kin !< Set number of kinetic gases
    !> Allocate procedures
        procedure :: allocate_gases !< Allocate gases array
    !> Query procedures
        procedure :: is_gas_in_gas_phase !< Check if gas exists in phase
    !> Compute procedures (currently disabled)
        !procedure :: compute_log_act_coeffs_gas_phase
        !procedure :: compute_log_Jacobian_act_coeffs_gas_phase
    end type
    
!>PFLOTRAN:
    
  !type, public :: gas_species_type
  !>  PetscInt :: id
  !>  character(len=MAXWORDLENGTH) :: name
  !>  PetscReal :: itype
  !>  PetscReal :: molar_volume
  !>  PetscReal :: molar_weight
  !>  PetscBool :: print_me
  !>  type(database_rxn_type), pointer :: dbaserxn
  !>  type(gas_species_type), pointer :: next
  !end type gas_species_type
  !
  !type, public :: gas_type
  !
  !>  PetscInt :: ngas
  !>  PetscInt :: nactive_gas
  !>  PetscInt :: npassive_gas
  !
  !>  type(gas_species_type), pointer :: list
  !
  !>  !> gas species names
  !>  character(len=MAXWORDLENGTH), pointer :: active_names(:)
  !>  character(len=MAXWORDLENGTH), pointer :: passive_names(:)
  !>  PetscBool :: print_all
  !>  PetscBool :: print_concentration
  !>  PetscBool :: print_partial_pressure
  !>  PetscBool, pointer :: active_print_me(:)
  !>  PetscBool, pointer :: passive_print_me(:)
  !
  !>  PetscInt, pointer :: acteqspecid(:,:)   !> (0:ncomp in rxn)
  !>  PetscReal, pointer :: acteqstoich(:,:)
  !>  PetscInt, pointer :: acteqh2oid(:)       !> id of water, if present
  !>  PetscReal, pointer :: acteqh2ostoich(:)  !> stoichiometry of water, if present
  !>  PetscReal, pointer :: acteqlogK(:)
  !>  PetscReal, pointer :: acteqlogKcoef(:,:)
  !
  !>  PetscReal, pointer :: actmolarwt(:)
  !>  PetscReal, pointer :: pasmolarwt(:)
  !
  !>  PetscInt, pointer :: paseqspecid(:,:)   !> (0:ncomp in rxn)
  !>  PetscReal, pointer :: paseqstoich(:,:)
  !>  PetscInt, pointer :: paseqh2oid(:)       !> id of water, if present
  !>  PetscReal, pointer :: paseqh2ostoich(:)  !> stoichiometry of water, if present
  !>  PetscReal, pointer :: paseqlogK(:)
  !>  PetscReal, pointer :: paseqlogKcoef(:,:)
  !
  !end type gas_type
    
    contains 
        !subroutine compute_log_act_coeffs_gas_phase(this,ionic_strength,log_act_coeffs)
        !    import gas_phase_c
        !    implicit none
        !    class(gas_phase_c) :: this
        !    real(kind=8), intent(in) :: ionic_strength
        !    real(kind=8), intent(out) :: log_act_coeffs(:) !> must be allocated
        !end subroutine
        !
        !subroutine compute_log_Jacobian_act_coeffs_gas_phase(this,ionic_strength,log_act_coeffs,conc,log_Jacobian_act_coeffs)
        !    import gas_phase_c
        !    implicit none
        !    class(gas_phase_c) :: this
        !    real(kind=8), intent(in) :: ionic_strength
        !    real(kind=8), intent(in) :: log_act_coeffs(:)
        !    real(kind=8), intent(in) :: conc(:) !> concentration of gas species in a given target
        !    real(kind=8), intent(out) :: log_Jacobian_act_coeffs(:,:) !> must be allocated
        !end subroutine
    
        !> \brief Allocate gases array
        !> \details
        !>   Allocates memory for the gases array.
        !>   
        !>   **Workflow:**
        !>   1. Set num_species (if provided)
        !>   2. Deallocate existing array (if allocated)
        !>   3. Allocate new array with current num_species
        !>   
        !>   **Safety:**
        !>   Prevents memory leaks by deallocating before reallocation.
        !> \param[in,out] this Gas phase object
        !> \param[in] num_species Optional: total number of gas species [-]
        subroutine allocate_gases(this,num_species)
            implicit none
            class(gas_phase_c) :: this !< Gas phase object
            integer(kind=4), intent(in), optional :: num_species !< Optional: total number of gas species [-]
            !< Set number of species if provided
            if (present(num_species)) then
                this%num_species=num_species !< Update species count
            end if
            !< Deallocate existing array if allocated
            if (allocated(this%gases)) then
                deallocate(this%gases) !< Free memory
            end if
            !< Allocate array with current size
            allocate(this%gases(this%num_species))
        end subroutine
        
        !> \brief Set number of equilibrium gases
        !> \details
        !>   Sets the total count of gases participating in equilibrium reactions.
        !>   
        !>   **Relationship:**
        !>   num_gases_eq = num_gases_eq_cst_act + num_gases_eq_var_act
        !>   
        !>   **Examples:**
        !>   - CO2(g) dissolution: num_gases_eq = 1
        !>   - CO2 + O2 system: num_gases_eq = 2
        !> \param[in,out] this Gas phase object
        !> \param[in] num_gases_eq Number of equilibrium gases [-]
        subroutine set_num_gases_eq(this,num_gases_eq)
            implicit none
            class(gas_phase_c) :: this !< Gas phase object
            integer(kind=4), intent(in) :: num_gases_eq !< Number of equilibrium gases [-]
            !< Assign equilibrium gas count
            this%num_gases_eq=num_gases_eq
        end subroutine
        
        !> \brief Set number of equilibrium gases with constant activity
        !> \details
        !>   Sets the count of equilibrium gases with fixed partial pressure.
        !>   
        !>   **Constant activity gases:**
        !>   - Fixed partial pressure (e.g., P_O2 = 0.21 atm)
        !>   - Act as infinite reservoir
        !>   - Do not consume/produce during simulation
        !>   
        !>   **Typical examples:**
        !>   - Atmospheric O2 in open system
        !>   - Fixed CO2 in controlled experiment
        !> \param[in,out] this Gas phase object
        !> \param[in] num_gases_eq_cst_act Number of equilibrium gases with constant activity [-]
        subroutine set_num_gases_eq_cst_act(this,num_gases_eq_cst_act)
            implicit none
            class(gas_phase_c) :: this !< Gas phase object
            integer(kind=4), intent(in) :: num_gases_eq_cst_act !< Number of equilibrium gases with constant activity [-]
            !< Assign constant activity equilibrium gas count
            this%num_gases_eq_cst_act=num_gases_eq_cst_act
        end subroutine
        
        !> \brief Set number of equilibrium gases with variable activity
        !> \details
        !>   Sets the count of equilibrium gases with evolving partial pressure.
        !>   
        !>   **Variable activity gases:**
        !>   - Partial pressure changes during simulation
        !>   - Finite amount in system
        !>   - Consume/produce according to reactions
        !>   
        !>   **Typical examples:**
        !>   - CO2 in closed system
        !>   - CH4 production in sediments
        !>   - H2S generation in anoxic conditions
        !> \param[in,out] this Gas phase object
        !> \param[in] num_gases_eq_var_act Number of equilibrium gases with variable activity [-]
        subroutine set_num_gases_eq_var_act(this,num_gases_eq_var_act)
            implicit none
            class(gas_phase_c) :: this !< Gas phase object
            integer(kind=4), intent(in) :: num_gases_eq_var_act !< Number of equilibrium gases with variable activity [-]
            !< Assign variable activity equilibrium gas count
            this%num_gases_eq_var_act=num_gases_eq_var_act
        end subroutine
        
        !> \brief Set number of kinetic gases with constant activity
        !> \details
        !>   Sets the count of kinetic gases with fixed partial pressure.
        !>   
        !>   **Kinetic gases with constant activity:**
        !>   - Participate in rate-limited reactions
        !>   - Fixed partial pressure (buffered)
        !>   - Reaction rate depends on P_gas
        !>   
        !>   **Usage example:**
        !>   Aerobic respiration with fixed O2 from atmosphere.
        !> \param[in,out] this Gas phase object
        !> \param[in] num_gases_kin_cst_act Number of kinetic gases with constant activity [-]
        subroutine set_num_gases_kin_cst_act(this,num_gases_kin_cst_act)
            implicit none
            class(gas_phase_c) :: this !< Gas phase object
            integer(kind=4), intent(in) :: num_gases_kin_cst_act !< Number of kinetic gases with constant activity [-]
            !< Assign constant activity kinetic gas count
            this%num_gases_kin_cst_act=num_gases_kin_cst_act
        end subroutine
        
        !> \brief Set number of kinetic gases with variable activity
        !> \details
        !>   Sets the count of kinetic gases with evolving partial pressure.
        !>   
        !>   **Kinetic gases with variable activity:**
        !>   - Participate in rate-limited reactions
        !>   - Partial pressure evolves during simulation
        !>   - Both P_gas and rate affect reaction progress
        !>   
        !>   **Usage example:**
        !>   Methane production in closed anaerobic system.
        !> \param[in,out] this Gas phase object
        !> \param[in] num_gases_kin_var_act Number of kinetic gases with variable activity [-]
        subroutine set_num_gases_kin_var_act(this,num_gases_kin_var_act)
            implicit none
            class(gas_phase_c) :: this !< Gas phase object
            integer(kind=4), intent(in) :: num_gases_kin_var_act !< Number of kinetic gases with variable activity [-]
            !< Assign variable activity kinetic gas count
            this%num_gases_kin_var_act=num_gases_kin_var_act
        end subroutine

        !> \brief Set total number of kinetic gases
        !> \details
        !>   Sets the total count of gases participating in kinetic reactions.
        !>   
        !>   **Relationship:**
        !>   num_gases_kin = num_gases_kin_cst_act + num_gases_kin_var_act
        !>   
        !>   **Kinetic reactions:**
        !>   - Rate-limited (not instantaneous equilibrium)
        !>   - Require rate law specification
        !>   - Time-dependent evolution
        !> \param[in,out] this Gas phase object
        !> \param[in] num_gases_kin Number of kinetic gases [-]
        subroutine set_num_gases_kin(this,num_gases_kin)
            implicit none
            class(gas_phase_c) :: this !< Gas phase object
            integer(kind=4), intent(in) :: num_gases_kin !< Number of kinetic gases [-]
            !< Assign kinetic gas count
            this%num_gases_kin=num_gases_kin
        end subroutine

        
        !> \brief Check if gas exists in gas phase
        !> \details
        !>   Searches for a gas by name in the gases array.
        !>   
        !>   **Search algorithm:**
        !>   1. Initialize flag to FALSE and index to 0
        !>   2. Loop through all gas species
        !>   3. Compare names (case-sensitive)
        !>   4. If match found: set flag=TRUE, store index, exit loop
        !>   
        !>   **Use cases:**
        !>   - Validate gas exists before operations
        !>   - Get gas index for array access
        !>   - Check reactive zone composition
        !>   - Verify input data consistency
        !> \param[in] this Gas phase object
        !> \param[in] gas Gas to search for
        !> \param[out] flag TRUE if gas found, FALSE otherwise
        !> \param[out] gas_ind Optional: index in gases array (0 if not found)
        subroutine is_gas_in_gas_phase(this,gas,flag,gas_ind) !< Checks if gas belongs to gas phase
            implicit none
            class(gas_phase_c), intent(in) :: this !< Gas phase object
            class(gas_species_c), intent(in) :: gas !< Gas to search for
            logical, intent(out) :: flag !< TRUE if gas belongs to gas phase, FALSE otherwise
            integer(kind=4), intent(out), optional :: gas_ind !< Optional: index of gas in gas phase (0 if not present)
            
            integer(kind=4) :: i !< Loop index for gas species
            
            !< Initialize search result to "not found"
            flag=.false.
            !< Initialize index to zero if requested
            if (present(gas_ind)) then
                gas_ind=0 !< Not found by default
            end if
            !< Loop through all gas species in phase
            do i=1,this%num_species
                !< Compare gas names (case-sensitive)
                if (gas%name==this%gases(i)%name) then
                    flag=.true. !< Mark as found
                    !< Store index if requested
                    if (present(gas_ind)) then
                        gas_ind=i !< Position in gases array
                    end if
                    exit !< Stop searching (found)
                end if
            end do
        end subroutine
        
    !> \brief Check if two gas phases are equal
    !> \details
    !>   Compares two gas phases for equality by checking:
    !>   1. Number of species
    !>   2. Number of equilibrium gases
    !>   3. Number of constant activity species
    !>   4. Number of equilibrium gases with constant activity
    !>   5. Presence of all gases from phase 1 in phase 2
    !>   
    !>   **Equality criteria:**
    !>   - Same total number of species
    !>   - Same number of equilibrium gases
    !>   - Same number of constant activity species
    !>   - Same number of equilibrium constant activity gases
    !>   - All gases from phase 1 exist in phase 2
    !>   
    !>   **Note:**
    !>   This is a unidirectional check (phase 1 ⊆ phase 2).
    !>   For bidirectional equality, also check phase 2 ⊆ phase 1.
    !>   
    !>   **Use cases:**
    !>   - Compare reactive zones
    !>   - Validate reactive transport setup
    !>   - Check for zone changes during simulation
    !> \param[in] gas_phase_1 First gas phase to compare
    !> \param[in] gas_phase_2 Second gas phase to compare
    !> \param[out] flag TRUE if gas phases equal, FALSE otherwise
    subroutine are_gas_phases_equal(gas_phase_1,gas_phase_2,flag)
    class(gas_phase_c), intent(in) :: gas_phase_1 !< First gas phase
    class(gas_phase_c), intent(in) :: gas_phase_2 !< Second gas phase
    logical, intent(out) :: flag !< TRUE if equal, FALSE otherwise

    integer(kind=4) :: i !< Loop index for gas species
    logical :: gas_flag !< Flag indicating if individual gas found

    !< Initialize result to "equal" (assume equal until proven otherwise)
    flag=.true.
    !< Check if total species count and equilibrium gas count match
    if (gas_phase_1%num_species/=gas_phase_2%num_species .or. gas_phase_1%num_gases_eq/=gas_phase_2%num_gases_eq) then
        flag=.false. !< Different species or equilibrium gas counts
    !< Check if constant activity species counts match
    else if (gas_phase_1%num_cst_act_species/=gas_phase_2%num_cst_act_species .or. gas_phase_1%num_gases_eq_cst_act/=&
            gas_phase_2%num_gases_eq_cst_act) then
        flag=.false. !< Different constant activity counts
    !< Check if all gases from phase 1 exist in phase 2
    else
        !< Loop through all gases in phase 1
        do i=1,gas_phase_1%num_species
            !< Check if gas exists in phase 2
            call gas_phase_2%is_gas_in_gas_phase(gas_phase_1%gases(i),gas_flag)
            !< If any gas not found in phase 2, phases are not equal
            if (gas_flag .eqv. .false.) then
                flag=.false. !< Gas from phase 1 missing in phase 2
                exit !< Stop checking (already not equal)
            end if
        end do
    end if
    end subroutine
end module