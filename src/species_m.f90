!> \file species_m.f90
!> \brief Defines the species, element, and microorganism classes for chemical modeling.
!> \details 
!>   This module provides object-oriented definitions and routines for chemical species, 
!>   including elements, general species, and microorganisms. It supports activity coefficients, 
!>   specific volume calculations, and species classification (gas, surface complex, etc). 
!>   All routines and variables are documented inline for clarity and maintainability.
!>
!>   Class hierarchy:
!>   - element_c: Base class for chemical elements
!>   - species_c: Extends element_c, adds activity and thermodynamic parameters
!>   - microorganism_c: Extends species_c, adds biological growth/decay rates
!>
!>   Key capabilities:
!>   - Species name, valence, molecular weight, molar volume
!>   - Activity coefficient parameters (Debye-Hückel, Davies, etc.)
!>   - Specific volume parameters (Redlich equation)
!>   - Constant activity flag for pure phases and water
!>   - Gas species detection (via "(g)" suffix)
!>   - Surface complex detection (via "x-" or "x2-" prefix)
!>   - Array operations (append, compare)
!>
!> \author jordi Petchamé-Guerrero
!> \date 2025
module species_m
    use params_act_coeff_m      !< Module for activity coefficient parameters
    use params_spec_vol_m       !< Module for specific volume parameters
    implicit none
    save                        !< Preserve module variables between calls
    private                     !< Default accessibility is private
    !> \brief Base class for chemical elements
    !> \details 
    !>   Stores fundamental properties of a chemical element:
    !>   - Name (character string, max 256 characters)
    !>   - Valence (oxidation state, default 0)
    !>   - Molecular weight in kg/mol
    !>   - Molar volume in L/mol
    !>
    !>   Used as base class for species_c.
    type, public :: element_c
        character(len=256)  :: name              !< Element name (including oxidation state if applicable)
        integer(kind=4)     :: valence=0         !< Oxidation state (default 0 for neutral)
        real(kind=8)        :: molecular_weight=0d0  !< Molecular weight [kg/mol]
        real(kind=8)        :: mol_vol=0d0       !< Molar volume [L/mol]
    contains
        procedure :: set_name            !< Set element name
        procedure :: set_valence         !< Set oxidation state
        procedure :: set_molecular_weight !< Set molecular weight [kg/mol]
        procedure :: set_mol_vol         !< Set molar volume [L/mol]
    end type
    
    !> \brief Class for chemical species (extends element_c)
    !> \details 
    !>   Extends element_c to include:
    !>   - Pointer to base element (for PHREEQC compatibility)
    !>   - Constant activity flag (for pure phases, water)
    !>   - Activity coefficient parameters (Debye-Hückel, Davies, etc.)
    !>   - Specific volume parameters (Redlich equation)
    !>
    !>   Supports detection of:
    !>   - Gas species (name contains "(g)")
    !>   - Surface complexes (name contains "x-" or "x2-")
    type, public, extends(element_c) :: species_c
        class(element_c), pointer :: element     !< Pointer to base element (PHREEQC-style)
        logical :: cst_act_flag                  !< TRUE if species has constant activity (e.g., H2O, pure solids)
        type(params_act_coeff_c) :: params_act_coeff         !< Activity coefficient parameters
        type(params_spec_vol_Redlich_c) :: params_spec_vol_Redlich !< Redlich specific volume parameters
    contains
        procedure :: set_cst_act_flag    !< Set constant activity flag
        procedure :: set_element         !< Associate species with base element
        procedure :: set_params_act_coeff !< Set activity coefficient parameters
        procedure :: set_params_spec_vol_Redlich !< Set Redlich specific volume parameters
        procedure :: copy_species      !< Copy all properties from another species
        procedure :: read_species        !< Parse species from string
        procedure :: is_gas              !< Check if species is a gas (name contains "(g)")
        procedure :: is_surf_compl       !< Check if species is a surface complex (name contains "x-" or "x2-")
    end type
    
    !> \brief Class for microorganisms (extends species_c)
    !> \details 
    !>   Extends species_c to include biological parameters:
    !>   - Maximum growth rate [1/d]
    !>   - Decay rate [1/d]
    !>   - Decay flag (whether to consider decay)
    !>
    !>   Used for modeling microbial-mediated reactions in reactive transport.
    type, public, extends(species_c) :: microorganism_c
        real(kind=8) :: max_growth_rate          !< Maximum specific growth rate [1/d]
        real(kind=8) :: decay_rate               !< Decay rate [1/d]
        logical :: decay_flag                    !< TRUE if decay is modeled, FALSE otherwise
        real(kind=8) :: yield                   !< Biomass yield coefficient: rate of biomass produced to substrate consumed [M/M]
    end type
    
!> \brief PFLOTRAN species type (alternative implementation, currently not used)
!> \details 
!>   This commented-out type definition shows an alternative species structure
!>   used in PFLOTRAN reactive transport code. Included for reference.
!>
!>   Fields would include:
!>   - name: Species name (max word length)
!>   - id: Unique identifier
!>   - molar_weight: Molecular weight
!>   - mnrl_molar_density: Mineral molar density [mol/m³-mineral]
!>   - solubility_limit: Maximum solubility [mol/m³-liquid]
!>   - ele_kd: Distribution coefficient [m³-water/m³-bulk]
!>   - radioactive: Flag for radioactive species
!>   - print_me: Flag for output
!>   - next: Pointer to next species in linked list
  !type, public :: species_type
  !  character(len=MAXWORDLENGTH) :: name
  !  PetscInt :: id
  !  PetscReal :: molar_weight
  !  PetscReal :: mnrl_molar_density  !> [mol/m^3-mnrl]
  !  PetscReal :: solubility_limit    !> [mol/m^3-liq]
  !  PetscReal :: ele_kd              !> [m^3-water/m^3-bulk
  !  PetscBool :: radioactive
  !  PetscBool :: print_me
  !  type(species_type), pointer :: next
  !end type species_type
    
    !> \brief Interface for species comparison and reading operations
    !> \details 
    !>   Defines generic interfaces for:
    !>   - compare_species_arrays: Compare two species arrays element-wise
    !>   - read_species: Parse species from string (implementation elsewhere)
    interface
        !> \brief Compare two arrays of species for equality
        !> \param[in]  species_array_1 First array of species to compare
        !> \param[in]  species_array_2 Second array of species to compare
        !> \param[out] flag            TRUE if arrays are identical, FALSE otherwise
        !> \details 
        !>   Compares species arrays element-wise by name.
        !>   Arrays must have same size or error is raised.
        subroutine compare_species_arrays(species_array_1,species_array_2,flag)
            import species_c
            implicit none
            class(species_c), intent(in) :: species_array_1(:)  !< First species array
            class(species_c), intent(in) :: species_array_2(:)  !< Second species array
            logical, intent(out) :: flag                        !< TRUE if equal, FALSE otherwise
        end subroutine
        
        !> \brief Read species from string (interface)
        !> \param[in,out] this Species object to populate
        !> \param[in]     str  Input string to parse
        !> \details 
        !>   Parses species properties from formatted string.
        !>   Implementation must be provided elsewhere in codebase.
        subroutine read_species(this,str)
            import species_c
            implicit none
            class(species_c) :: this                            !< Species object to populate
            character(len=*), intent(in) :: str                 !< Input string containing species data
        end subroutine
    end interface
    
    contains
    
        !> \brief Set the name of an element
        !> \param[in,out] this Element object to modify
        !> \param[in]     name Element name to assign
        !> \details 
        !>   Simple assignment of element name string.
        !>   Name can include oxidation state information.
        subroutine set_name(this,name)
            implicit none
            class(element_c) :: this                            !< Element object
            character(len=*), intent(in) :: name                !< Name to assign
            this%name=name                                      !< Direct assignment
        end subroutine
        
        !> \brief Set the valence (oxidation state) of an element
        !> \param[in,out] this    Element object to modify
        !> \param[in]     valence Valence/oxidation state to assign
        !> \details 
        !>   Sets the oxidation state of the element.
        !>   Examples: Fe3+ would have valence=3, Cl- would have valence=-1
        subroutine set_valence(this,valence)
            implicit none
            class(element_c) :: this                            !< Element object
            integer(kind=4), intent(in) :: valence              !< Oxidation state
            this%valence=valence                                !< Direct assignment
        end subroutine
        
        !> \brief Set the element pointer for a species
        !> \param[in,out] this    Species object to modify
        !> \param[in]     element Element to associate with this species
        !> \details 
        !>   Creates a pointer association between species and its base element.
        !>   Used in PHREEQC-style formulations where species derive from master elements.
        subroutine set_element(this,element)
            implicit none
            class(species_c) :: this                            !< Species object
            class(element_c), intent(in), target :: element     !< Element to point to (must be target)
            this%element=>element                               !< Pointer assignment
        end subroutine
        
        !> \brief Set the molecular weight of an element
        !> \param[in,out] this             Element object to modify
        !> \param[in]     molecular_weight Molecular weight [kg/mol]
        !> \details 
        !>   Assigns molecular weight in kg/mol.
        !>   No validation - caller must ensure positive value.
        subroutine set_molecular_weight(this,molecular_weight)
            implicit none
            class(element_c) :: this                            !< Element object
            real(kind=8), intent(in) :: molecular_weight        !< Molecular weight [kg/mol]
            this%molecular_weight=molecular_weight              !< Direct assignment
        end subroutine
        
        !> \brief Set constant activity flag for a species
        !> \param[in,out] this          Species object to modify
        !> \param[in]     cst_act_flag  TRUE if species has constant activity
        !> \details 
        !>   Constant activity species include:
        !>   - Pure water (H2O)
        !>   - Pure mineral phases
        !>   - Species in buffered systems
        !>
        !>   These species don't participate in activity calculations.
        subroutine set_cst_act_flag(this,cst_act_flag)
            implicit none
            class(species_c) :: this                            !< Species object
            logical, intent(in) :: cst_act_flag                 !< Constant activity flag
            this%cst_act_flag=cst_act_flag                      !< Direct assignment
        end subroutine
        
        !> \brief Set activity coefficient parameters for a species
        !> \param[in,out] this              Species object to modify
        !> \param[in]     params_act_coeff  Activity coefficient parameters
        !> \details 
        !>   Activity coefficient parameters include:
        !>   - Type (Debye-Hückel, Davies, Pitzer, etc.)
        !>   - Ion size parameter
        !>   - Debye-Hückel a and b parameters
        !>   - Temperature dependence
        !>
        !>   Used to compute activity coefficients from ionic strength.
        subroutine set_params_act_coeff(this,params_act_coeff)
            implicit none
            class(species_c) :: this                            !< Species object
            type(params_act_coeff_c) :: params_act_coeff        !< Activity coefficient parameters
            this%params_act_coeff=params_act_coeff              !< Direct assignment (derived type)
        end subroutine
        
        !> \brief Set Redlich specific volume parameters for a species
        !> \param[in,out] this                      Species object to modify
        !> \param[in]     params_spec_vol_Redlich   Redlich equation parameters
        !> \details 
        !>   Redlich equation parameters for computing specific volume:
        !>   - Partial molar volume at infinite dilution
        !>   - Debye-Hückel limiting slope
        !>   - Ion size parameter
        !>
        !>   Used for pressure effects on thermodynamic properties.
        subroutine set_params_spec_vol_Redlich(this,params_spec_vol_Redlich)
            implicit none
            class(species_c) :: this                            !< Species object
            type(params_spec_vol_Redlich_c) :: params_spec_vol_Redlich !< Redlich parameters
            this%params_spec_vol_Redlich=params_spec_vol_Redlich       !< Direct assignment (derived type)
        end subroutine
        
        !> \brief Assign all properties from another species
        !> \param[in,out] this    Species object to assign to (destination)
        !> \param[in]     species Source species to copy from
        !> \details 
        !>   Performs deep copy of all species properties:
        !>   1. Name string
        !>   2. Constant activity flag
        !>   3. Valence (oxidation state)
        !>   4. Activity coefficient parameters
        !>   5. Redlich specific volume parameters
        !>   6. Element pointer
        !>   7. Molecular weight
        !>   8. Molar volume
        !>
        !>   Used for creating copies and array operations.
        subroutine copy_species(this,species)
            implicit none
            class(species_c) :: this                            !< Destination species
            class(species_c), intent(in) :: species             !< Source species to copy
            this%name=species%name                              !< Copy name
            this%cst_act_flag=species%cst_act_flag              !< Copy constant activity flag
            this%valence=species%valence                        !< Copy valence
            this%params_act_coeff=species%params_act_coeff      !< Copy activity coefficient params
            this%params_spec_vol_Redlich=species%params_spec_vol_Redlich !< Copy specific volume params
            this%element=>species%element                       !< Copy element pointer
            this%molecular_weight=species%molecular_weight      !< Copy molecular weight
            this%mol_vol=species%mol_vol                        !< Copy molar volume
        end subroutine
        
        !> \brief Append a species to a dynamically allocated array
        !> \param[in]     this          Species to append
        !> \param[in,out] species_array Array to append to (reallocated)
        !> \details 
        !>   Dynamically grows array by one element and appends species.
        !>   
        !>   Algorithm:
        !>   1. Copy existing array to temporary (aux_array)
        !>   2. Deallocate original array
        !>   3. If aux_array has elements:
        !>      - Allocate new array with size+1
        !>      - Copy all existing elements
        !>      - Append new species at end
        !>   4. If aux_array is empty:
        !>      - Allocate size 1
        !>      - Assign new species
        !>
        !>   Note: O(n) operation - inefficient for large arrays with many appends.
        subroutine append_species(this,species_array)
            implicit none
            class(species_c), intent(in) :: this                !< Species to append
            type(species_c), intent(inout), allocatable :: species_array(:) !< Array to grow
            integer(kind=4) :: i                                !< Loop index
            type(species_c), allocatable :: aux_array(:)        !< Temporary copy of existing array
            
            aux_array=species_array                             !< Copy existing array
            deallocate(species_array)                           !< Free original
            if (size(aux_array)>0) then                         !< If array was non-empty
                allocate(species_array(size(aux_array)+1))      !< Allocate size+1
                do i=1,size(species_array)-1                    !< Loop over existing elements
                    call species_array(i)%copy_species(aux_array(i)) !< Copy each element
                end do
                call species_array(size(species_array))%copy_species(this) !< Append new species
            else                                                !< Array was empty
                allocate(species_array(1))                      !< Allocate size 1
                call species_array(1)%copy_species(this)      !< Assign new species
            end if
        end subroutine
        
        !> \brief Check if species is a gas
        !> \param[in]  this Species object to check
        !> \param[out] flag TRUE if species is a gas, FALSE otherwise
        !> \details 
        !>   Detects gas species by searching for "(g)" suffix in name.
        !>   
        !>   Examples:
        !>   - "CO2(g)" → TRUE
        !>   - "O2(g)"  → TRUE
        !>   - "Ca+2"   → FALSE
        !>
        !>   Uses Fortran intrinsic index() function which returns:
        !>   - Position > 0 if substring found
        !>   - 0 if substring not found
        subroutine is_gas(this,flag)
            implicit none
            class(species_c), intent(in) :: this                !< Species to check
            LOGICAL, intent(out) :: flag                        !< TRUE if gas, FALSE otherwise
            integer(KIND=4) :: ind                              !< Index position of "(g)" in name
            
            ind=index(this%name,'(g)')                          !< Search for "(g)" substring
            if (ind>0) then                                     !< If "(g)" found
                flag=.true.                                     !< Mark as gas
            else                                                !< If "(g)" not found
                flag=.false.                                    !< Not a gas
            end if
        end subroutine
        
        !> \brief Check if species is a surface complex
        !> \param[in]  this Species object to check
        !> \param[out] flag TRUE if species is a surface complex, FALSE otherwise
        !> \details 
        !>   Detects surface complexes by searching for "x-" or "x2-" prefix in name.
        !>   
        !>   Surface complexes bind to mineral surfaces via exchange sites.
        !>   Naming convention:
        !>   - "x-" prefix: Single-site complex
        !>   - "x2-" prefix: Double-site complex
        !>
        !>   Examples:
        !>   - "x-SOH"   → TRUE (single site)
        !>   - "x2-SO4"  → TRUE (double site)
        !>   - "Ca+2"    → FALSE (aqueous ion)
        subroutine is_surf_compl(this,flag)
            implicit none
            class(species_c), intent(in) :: this                !< Species to check
            LOGICAL, intent(out) :: flag                        !< TRUE if surface complex, FALSE otherwise
            integer(KIND=4) :: ind                              !< Index position of "x-"
            integer(KIND=4) :: ind2                             !< Index position of "x2-"
            
            ind=index(this%name,'x-')                           !< Search for "x-" (single site)
            ind2=index(this%name,'x2-')                         !< Search for "x2-" (double site)
            if (ind>0 .or. ind2>0) then                         !< If either prefix found
                flag=.true.                                     !< Mark as surface complex
            else                                                !< If neither prefix found
                flag=.false.                                    !< Not a surface complex
            end if
        end subroutine
        
        !> \brief Set the molar volume of an element
        !> \param[in,out] this    Element object to modify
        !> \param[in]     mol_vol Molar volume [L/mol]
        !> \details 
        !>   Sets molar volume with validation.
        !>   Raises error if negative value provided.
        !>
        !>   Molar volume used for:
        !>   - Pressure effects on equilibrium
        !>   - Density calculations
        !>   - Volume change during reactions
        subroutine set_mol_vol(this,mol_vol)
            implicit none
            class(element_c) :: this                            !< Element object
            real(kind=8), intent(in) :: mol_vol                 !< Molar volume [L/mol]
            
            if (mol_vol<0d0) error stop "Molar volume cannot be negative" !< Validate non-negative
            this%mol_vol=mol_vol                                !< Assign if valid
        end subroutine
        
        !> \brief Check if two species arrays are equal
        !> \param[in]  species_array_1 First array of species to compare
        !> \param[in]  species_array_2 Second array of species to compare
        !> \param[out] flag            TRUE if arrays are identical, FALSE otherwise
        !> \details 
        !>   Compares two species arrays element-wise by name only.
        !>   
        !>   Requirements:
        !>   - Both arrays must have same size (enforced by error check)
        !>   - Comparison is case-sensitive
        !>   - Only names are compared (not other properties)
        !>   
        !>   Algorithm:
        !>   1. Check arrays have same size
        !>   2. If empty arrays, return TRUE
        !>   3. Loop through elements
        !>   4. Compare names
        !>   5. Exit early if mismatch found
        subroutine are_species_arrays_equal(species_array_1,species_array_2,flag)
            implicit none
            class(species_c), intent(in) :: species_array_1(:)  !< First species array
            class(species_c), intent(in) :: species_array_2(:)  !< Second species array
            logical, intent(out) :: flag                        !< TRUE if equal, FALSE otherwise
            integer(kind=4) :: i                                !< Loop index
            integer(kind=4) :: n                                !< Array size
            
            flag=.true.                                         !< Initialize as equal
            n=size(species_array_1)                             !< Get size of first array
            if(size(species_array_2)/=n) error stop "Dimension error in compare_species_arrays" !< Validate equal sizes
            if (n>0) then                                       !< If arrays are non-empty
                do i=1,n                                        !< Loop over all elements
                    if (species_array_1(i)%name/=species_array_2(i)%name) then !< Compare names
                        flag=.false.                            !< Mark as not equal
                        exit                                    !< Exit early (optimization)
                    else                                        !< Names match
                        continue                                !< Check next element
                    end if
                end do
            end if
        end subroutine
        
        !> \brief Check if two individual species are equal
        !> \param[in]  species_1 First species to compare
        !> \param[in]  species_2 Second species to compare
        !> \param[out] flag      TRUE if species are identical, FALSE otherwise
        !> \details 
        !>   Compares two species across multiple properties:
        !>   - Name (character string)
        !>   - Constant activity flag (logical)
        !>   - Valence (integer)
        !>   - Molecular weight (real)
        !>   - Molar volume (real)
        !>
        !>   Commented-out checks for:
        !>   - Element properties
        !>   - Activity coefficient parameters
        !>   - Specific volume parameters
        !>   
        !>   These could be enabled for stricter equality checking.
        subroutine are_species_equal(species_1,species_2,flag)
            implicit none
            class(species_c), intent(in) :: species_1           !< First species
            class(species_c), intent(in) :: species_2           !< Second species
            logical, intent(out) :: flag                        !< TRUE if equal, FALSE otherwise
            
            flag=.true.                                         !< Initialize as equal
            !> Check name and constant activity flag
            if (species_1%name/=species_2%name .or. species_1%cst_act_flag .neqv. species_2%cst_act_flag) then
                flag=.false.                                    !< Names or flags differ
            !> Check valence and molecular weight
            else if (species_1%valence/=species_2%valence .or. species_1%molecular_weight/=species_2%molecular_weight) then
                flag=.false.                                    !< Valence or molecular weight differs
            !> Check molar volume
            else if (species_1%mol_vol/=species_2%mol_vol) then
                flag=.false.                                    !< Molar volume differs
            !> Additional checks (currently commented out):
            !else if (species_1%element%name/=species_2%element%name .or. species_1%element%valence/=species_2%element%valence) then
            !    flag=.false.                                   !< Element properties differ
            !else if (species_1%params_act_coeff%name/=species_2%params_act_coeff%name .or. species_1%params_act_coeff%type/=species_2%params_act_coeff%type) then
            !    flag=.false.                                   !< Activity coefficient parameters differ
            !else if (species_1%params_spec_vol_Redlich%name/=species_2%params_spec_vol_Redlich%name .or. species_1%params_spec_vol_Redlich%type/=species_2%params_spec_vol_Redlich%type) then
            !    flag=.false.                                   !< Specific volume parameters differ
            end if
        end subroutine
        
end module