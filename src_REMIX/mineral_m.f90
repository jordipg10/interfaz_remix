!> Mineral phase subclass
module mineral_m
    use phase_m, only: phase_c
    use solid_m, only: solid_c
    implicit none
    save
    type, public, extends(phase_c) :: mineral_c !> mineral class
        type(solid_c) :: mineral !> mineral that defines mineral phase
    contains
    end type
    
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
    interface
      
    end interface
    
    
    contains
        !subroutine set_mol_vol(this,mol_vol)
        !>    implicit none
        !>    class(mineral_c) :: this
        !>    real(kind=8), intent(in) :: mol_vol
        !>    if (mol_vol<0d0) error stop "Molar volume cannot be negative"
        !>    this%mol_vol=mol_vol
        !end subroutine
        
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