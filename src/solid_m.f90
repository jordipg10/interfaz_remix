module solid_m
    use species_m
    implicit none
    save
    type, public, extends(species_c) :: solid_c
    end type
    
!> PFLOTRAN:
    
!>  type, public :: solid_solution_type
!>    character(len=MAXWORDLENGTH) :: name
!>    PetscInt :: num_stoich_solid
!>    character(len=MAXWORDLENGTH), pointer :: stoich_solid_names(:)
!>    PetscInt, pointer :: stoich_solid_ids(:)
!#if 0
!>    PetscInt :: num_end_member
!>    type(stoichiometric_solid_type), pointer :: stoich_solid
!#endif
!>    type(solid_solution_type), pointer :: next
!>  end type solid_solution_type
!
!#if 0
!>  type, public :: stoichiometric_solid_type
!>    type(mineral_rxn_type), pointer :: mineral !> stoichiometric solid
!>    type(mineral_rxn_type), pointer :: end_members
!>    type(stoichiometric_solid_type), pointer :: next
!>  end type stoichiometric_solid_type
!
!>  type, public :: solid_solution_rxn_type
!>    character(len=MAXSTRINGLENGTH) :: database_filename
!>    PetscInt :: num_dbase_temperatures
!>    PetscReal, pointer :: dbase_temperatures(:)
!>    type(solid_solution_type), pointer :: list
!>    type(mineral_type), pointer :: mineral
!>  end type solid_solution_rxn_type
end module