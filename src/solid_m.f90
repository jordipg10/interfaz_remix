!> \file solid_m.f90
!> \brief Solid species module for solid species (minerals, surface complexes)
!> \details
!>   This module defines the solid species class for representing solid species.
!>   
!>   **Key features:**
!>   - Extends generic species_c class
!>   - Represents minerals, surface complexes, precipitates
!>   - Used in equilibrium and kinetic reactions
!>   
!>   **Applications:**
!>   - Mineral dissolution/precipitation
!>   - Surface complexation
!>   - Cation exchange sites
!>   - Reactive transport with mineral reactions
!>   
!>   **Inheritance:**
!>   solid_species_c → species_c (inherits name, valence, activity flags)
!>   
!>   **Usage example:**
!>   ```fortran
!>   type(solid_species_c) :: calcite
!>   calcite%name = "CaCO3"
!>   calcite%cst_act_flag = .false.  ! Variable activity
!>   ```
!>   
!> \author Jordi Petchamé-Guerrero
!> \date October 2025
module solid_m
    use species_m, only: species_c !< Import generic species class
    implicit none !< Enforce explicit variable declarations
    save !< Preserve module state between calls
    private !< Default accessibility is private
    !> \brief Solid species type for minerals and surface complexes
    !> \details
    !>   Defines a solid species by extending the generic species_c class.
    !>   
    !>   **Inherits from species_c:**
    !>   - name: Solid name (e.g., "Calcite", "Quartz")
    !>   - valence: Charge (typically 0 for minerals)
    !>   - cst_act_flag: TRUE if constant activity (infinite reservoir)
    !>   - Other species attributes
    !>   
    !>   **Solid types:**
    !>   - Minerals: CaCO3, SiO2, FeOOH
    !>   - Surface complexes: X-Na, X2-Ca (cation exchange)
    !>   - Precipitates: Variable activity solids
    !>   
    !>   **Constant vs Variable activity:**
    !>   - Constant: Pure phase with infinite reservoir (activity = 1)
    !>   - Variable: Finite amount, activity can change
    !>   
    !>   **No additional attributes:**
    !>   All properties inherited from species_c.
    type, public, extends(species_c) :: solid_species_c !< Solid species (extends species_c)
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