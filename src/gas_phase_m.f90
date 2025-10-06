!> Gas phase subclass: contains gases in gas phase
module gas_phase_m
    use phase_m, only: phase_c
    use gas_m, only: gas_c
    implicit none
    save
    type, public, extends(phase_c) :: gas_phase_c
        !> first gases in equilibrium with constant activity, then equilibrium with variable activity, & 
        !! then kinetic with constant activity, then kinetic with variable activity
        type(gas_c), allocatable :: gases(:) !> gases (first equilibrium, then kinetic)
        integer(kind=4) :: num_gases_eq=0 !> number of gases in equilibrium
        integer(kind=4) :: num_gases_kin=0 !> number of gases in kinetic reactions
        integer(kind=4) :: num_gases_eq_cst_act=0 !> number of gases in equilibrium with constant activity
        integer(kind=4) :: num_gases_eq_var_act=0 !> number of gases in equilibrium with variable activity
        integer(kind=4) :: num_gases_kin_cst_act=0 !> number of gases in equilibrium with variable activity
        integer(kind=4) :: num_gases_kin_var_act=0 !> number of gases in equilibrium with variable activity
    contains
    !> Set
        procedure, public :: set_num_gases_eq
        procedure, public :: set_num_gases_eq_cst_act
        procedure, public :: set_num_gases_eq_var_act
        procedure, public :: set_num_gases_kin_cst_act
        procedure, public :: set_num_gases_kin_var_act
        procedure, public :: set_num_gases_kin
    !> Allocate
        procedure, public :: allocate_gases
    !> Is
        procedure, public :: is_gas_in_gas_phase
    !> Compute
        !procedure, public :: compute_log_act_coeffs_gas_phase
        !procedure, public :: compute_log_Jacobian_act_coeffs_gas_phase
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
        !subroutine compute_log_act_coeffs_gas_phase(this,ionic_act,log_act_coeffs)
        !    import gas_phase_c
        !    implicit none
        !    class(gas_phase_c) :: this
        !    real(kind=8), intent(in) :: ionic_act
        !    real(kind=8), intent(out) :: log_act_coeffs(:) !> must be allocated
        !end subroutine
        !
        !subroutine compute_log_Jacobian_act_coeffs_gas_phase(this,ionic_act,log_act_coeffs,conc,log_Jacobian_act_coeffs)
        !    import gas_phase_c
        !    implicit none
        !    class(gas_phase_c) :: this
        !    real(kind=8), intent(in) :: ionic_act
        !    real(kind=8), intent(in) :: log_act_coeffs(:)
        !    real(kind=8), intent(in) :: conc(:) !> concentration of gas species in a given target
        !    real(kind=8), intent(out) :: log_Jacobian_act_coeffs(:,:) !> must be allocated
        !end subroutine
    
        subroutine allocate_gases(this,num_species)
            implicit none
            class(gas_phase_c) :: this
            integer(kind=4), intent(in), optional :: num_species
            if (present(num_species)) then
                this%num_species=num_species
            end if
            if (allocated(this%gases)) then
                deallocate(this%gases)
            end if
            allocate(this%gases(this%num_species))
        end subroutine
        
        subroutine set_num_gases_eq(this,num_gases_eq)
            implicit none
            class(gas_phase_c) :: this
            integer(kind=4), intent(in) :: num_gases_eq
            this%num_gases_eq=num_gases_eq
        end subroutine
        
        subroutine set_num_gases_eq_cst_act(this,num_gases_eq_cst_act)
            implicit none
            class(gas_phase_c) :: this
            integer(kind=4), intent(in) :: num_gases_eq_cst_act
            this%num_gases_eq_cst_act=num_gases_eq_cst_act
        end subroutine
        
        subroutine set_num_gases_eq_var_act(this,num_gases_eq_var_act)
            implicit none
            class(gas_phase_c) :: this
            integer(kind=4), intent(in) :: num_gases_eq_var_act
            this%num_gases_eq_var_act=num_gases_eq_var_act
        end subroutine
        
        subroutine set_num_gases_kin_cst_act(this,num_gases_kin_cst_act)
            implicit none
            class(gas_phase_c) :: this
            integer(kind=4), intent(in) :: num_gases_kin_cst_act
            this%num_gases_kin_cst_act=num_gases_kin_cst_act
        end subroutine
        
        subroutine set_num_gases_kin_var_act(this,num_gases_kin_var_act)
            implicit none
            class(gas_phase_c) :: this
            integer(kind=4), intent(in) :: num_gases_kin_var_act
            this%num_gases_kin_var_act=num_gases_kin_var_act
        end subroutine

        subroutine set_num_gases_kin(this,num_gases_kin)
            implicit none
            class(gas_phase_c) :: this
            integer(kind=4), intent(in) :: num_gases_kin
            this%num_gases_kin=num_gases_kin
        end subroutine

        
        subroutine is_gas_in_gas_phase(this,gas,flag,gas_ind) !> checks if gas belongs to gas phase
            implicit none
            class(gas_phase_c), intent(in) :: this !> gas phase
            class(gas_c), intent(in) :: gas !> gas
            logical, intent(out) :: flag !> TRUE if it belongs to gas phase, FALSE otherwise
            integer(kind=4), intent(out), optional :: gas_ind !> index gas in gas phase (0 if not present)
            
            integer(kind=4) :: i
            
            flag=.false.
            if (present(gas_ind)) then
                gas_ind=0
            end if
            do i=1,this%num_species
                if (gas%name==this%gases(i)%name) then
                    flag=.true.
                    if (present(gas_ind)) then
                        gas_ind=i
                    end if
                    exit
                end if
            end do
        end subroutine
        
    subroutine are_gas_phases_equal(gas_phase_1,gas_phase_2,flag)
    class(gas_phase_c), intent(in) :: gas_phase_1 !> first gas phase
    class(gas_phase_c), intent(in) :: gas_phase_2 !> second gas phase
    logical, intent(out) :: flag !> TRUE if equal, FALSE otherwise

    integer(kind=4) :: i
    logical :: gas_flag !> gas flag

    flag=.true.
    if (gas_phase_1%num_species/=gas_phase_2%num_species .or. gas_phase_1%num_gases_eq/=gas_phase_2%num_gases_eq) then
        flag=.false.
    else if (gas_phase_1%num_cst_act_species/=gas_phase_2%num_cst_act_species .or. gas_phase_1%num_gases_eq_cst_act/=&
            gas_phase_2%num_gases_eq_cst_act) then
        flag=.false.
    else
        do i=1,gas_phase_1%num_species
            call gas_phase_2%is_gas_in_gas_phase(gas_phase_1%gases(i),gas_flag)
            if (gas_flag .eqv. .false.) then
                flag=.false.
                exit
            end if
        end do
    end if
    end subroutine
end module