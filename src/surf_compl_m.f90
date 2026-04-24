!> \file surf_compl_m.f90
!> \brief Surface complexation module for surface sites and complexes
!> \details
!>   This module defines surface complexation and cation exchange properties.
!>   
!>   **Key components:**
!>   - surface_c: Surface phase with free sites and surface complexes
!>   - cat_exch_zone_c: Cation exchange (extends surface_c)
!>   
!>   **Applications:**
!>   - Ion exchange on clay minerals
!>   - Surface complexation on metal oxides
!>   - Adsorption reactions in reactive transport
!>   
!>   **Activity coefficient conventions:**
!>   - Gaines-Thomas
!>   - Gapon
!>   - Vanselow
!>   
!> \author Jordi Petchamé-Guerrero
!> \date October 2025
module surf_compl_m
    use solid_m, only: solid_species_c !< Import solid phase class
    use aq_species_m, only: aq_species_c !< Import phase and aqueous species classes
    use exch_sites_conv_m, only: exch_sites_conv_c !< Import exchange site conventions
    use phase_m, only: phase_c !< Import generic phase class
    implicit none !< Enforce explicit variable declarations
    save !< Preserve module state between calls
    private !< Default accessibility is private
    !> \brief Surface phase type for surface complexation
    !> \details
    !>   Defines a surface phase with free sites and bound surface complexes.
    !>   Extends the generic phase_c class.
    !>   
    !>   **Components:**
    !>   - free_site: Unoccupied surface sites
    !>   - surf_compl: Array of surface complexes (occupied sites)
    !>   
    !>   **Usage example:**
    !>   ```fortran
    !>   type(surface_c) :: surface
    !>   call surface%set_num_surf_compl(3)
    !>   call surface%allocate_surf_compl()
    !>   ```
    type, public, extends(phase_c) :: surface_c !< Surface phase (subclass of phase class)
        type(solid_species_c) :: free_site !< Free (unoccupied) surface site
        integer(kind=4) :: num_surf_compl=0 !< Number of surface complexes [-]
        type(solid_species_c), allocatable :: surf_compl(:) !< Array of surface complexes
    contains
        procedure :: set_num_surf_compl !< Set number of surface complexes
        procedure :: set_surf_compl !< Set surface complexes array
        procedure :: allocate_surf_compl !< Allocate surface complexes array
        procedure :: is_surf_compl_in !< Check if surface complex exists
    end type
    
    !> \brief Cation exchange type for ion exchange reactions
    !> \details
    !>   Defines cation exchange properties, extending surface_c.
    !>   Used for modeling ion exchange on clay minerals and zeolites.
    !>   
    !>   **Components:**
    !>   - num_exch_cats: Number of exchangeable cations
    !>   - exch_cat_indices: Indices linking to aqueous phase
    !>   - convention: Activity coefficient convention (Gaines-Thomas, Gapon, Vanselow)
    !>   
    !>   **Activity coefficient models:**
    !>   - Gaines-Thomas: Rational activity coefficients
    !>   - Gapon: Activity = mole fraction
    !>   - Vanselow: Activity = equivalent fraction
    !>   
    !>   **Usage example:**
    !>   ```fortran
    !>   type(cat_exch_zone_c) :: cation_exchange
    !>   call cation_exchange%set_num_exch_cats(2)  ! Na+, Ca2+
    !>   call cation_exchange%allocate_exch_cat_indices()
    !>   ```
    type, public, extends(surface_c) :: cat_exch_zone_c !< Cation exchange (subclass of surface class)
        integer(kind=4) :: num_exch_cats=0 !< Number of exchangeable cations [-]
        integer(kind=4), allocatable :: exch_cat_indices(:) !< Indices of exchangeable cations in aqueous phase object [-]
        class(exch_sites_conv_c), pointer :: convention !< Convention to compute activity coefficients
    contains
    !> Set procedures
        procedure :: set_convention !< Set activity coefficient convention
        procedure :: set_num_exch_cats !< Set number of exchangeable cations
    !> Allocate procedures
        procedure :: allocate_exch_cat_indices !< Allocate exchangeable cation indices array
    !> Compute procedures
        procedure :: compute_log_act_coeffs_ads_cats !< Compute log activity coefficients
        procedure :: compute_log_Jacobian_act_coeffs_ads_cats !< Compute Jacobian of log activity coefficients
        procedure :: compute_num_surf_compl !< Compute number of surface complexes
    end type
!>************************************************************************************************!
!> CHEPROO:
!type, public,::t_parentsurface
!> 
!character(len=100)                       :: name          !> Name of interphase
!> 
!type(t_pspecies), pointer, dimension(:)  :: pspecies      !> Pointer to chemical species
!> 
!type(t_phase), pointer                   :: paqph
!> 
!real*8, pointer, dimension(:,:)          :: propsite      !> Value of the properties of the sites
!> 
!character(len=100), pointer, dimension(:):: namepropsite  !> Name of sites properties
!> 
!integer                                  :: numsite       !> Number of site of the interphase
!
!integer                                  :: numsp         !> Total number of species of the interp
!
!integer                                  :: numpropsite   !> Number of properties of sites
!> 
!integer, pointer, dimension(:)           :: numspsite     !> Number of species for site [nsite]
!
!integer, pointer, dimension(:)           :: idxoh
!
!logical                                  :: locksp        !> .true. if species vector was allocate
!
!logical                                  :: lockpaqph     !> .true. if aqueous phase class was
!
!logical                                  :: lockpropsite
!
!logical                                  :: lockidxoh
!> 
!end type t_parentsurface
    
!type, public,::t_surface_cexch
!> 
!private                             ::
!> 
!type (t_parentsurface), pointer     :: pp        !> Poiter to parent surface
!> 
!real*8, dimension(:,:), pointer     :: cec       !> Cation exchange capacity 
!> 
!logical                             :: lockcec
!> 
!end type t_surface_cexch
!>************************************************************************************************!
!> PFLOTRAN:
  !type, public :: surface_complex_type
  !>  PetscInt :: id
  !>  character(len=MAXWORDLENGTH) :: name
  !>  character(len=MAXWORDLENGTH) :: free_site_name
  !>  PetscReal :: free_site_stoich
  !>  PetscReal :: Z
  !>  PetscReal :: forward_rate
  !>  PetscReal :: backward_rate
  !>  PetscBool :: print_me
  !>  !> pointer that can be used to index the master list
  !>  type(surface_complex_type), pointer :: ptr
  !>  type(database_rxn_type), pointer :: dbaserxn
  !>  type(surface_complex_type), pointer :: next
  !end type surface_complex_type
    
!> type, public :: surface_complexation_rxn_type
!>    PetscInt :: id
!>    PetscInt :: itype
!>    PetscInt :: free_site_id
!>    character(len=MAXWORDLENGTH) :: free_site_name
!>    PetscBool :: free_site_print_me
!>    PetscBool :: site_density_print_me
!>    PetscInt :: surface_itype
!>    PetscInt :: mineral_id
!>    character(len=MAXWORDLENGTH) :: surface_name
!>    PetscReal :: site_density !> site density in mol/m^3 bulk
!>    PetscReal, pointer :: rates(:)
!>    PetscReal, pointer :: site_fractions(:)
!>    PetscReal :: kinmr_scale_factor
!>    type(surface_complex_type), pointer :: complex_list
!>    type (surface_complexation_rxn_type), pointer :: next
!end type surface_complexation_rxn_type
    
!type, public :: ion_exchange_rxn_type
!>  PetscInt :: id
!>  character(len=MAXWORDLENGTH) :: mineral_name
!>  type(ion_exchange_cation_type), pointer :: cation_list
!>  PetscReal :: CEC
!>  type(ion_exchange_rxn_type), pointer :: next
!end type ion_exchange_rxn_type
!
!type, public :: ion_exchange_cation_type
!>  character(len=MAXWORDLENGTH) :: name
!>  PetscReal :: k
!>  type(ion_exchange_cation_type), pointer :: next
!end type ion_exchange_cation_type
    
    contains
        !> \brief Compute log activity coefficients for adsorbed cations
        !> \details
        !>   Computes log10 activity coefficients for all exchangeable cations
        !>   using the selected convention (Gaines-Thomas, Gapon, or Vanselow).
        !>   
        !>   **Algorithm:**
        !>   1. Loop through all exchangeable cations
        !>   2. For each cation, call convention-specific computation
        !>   3. Store result in log_act_coeffs array
        !>   
        !>   **Activity coefficient depends on:**
        !>   - Cation valence
        !>   - Cation exchange capacity (CEC)
        !>   - Selected convention
        !> \param[in] this Cation exchange object
        !> \param[in] valences Valences of exchangeable cations [num_exch_cats]
        !> \param[in] CEC Cation exchange capacity [meq/g] or [mol/kg]
        !> \param[out] log_act_coeffs Log10 activity coefficients [-] (must be pre-allocated, dim=num_exch_cats)
        subroutine compute_log_act_coeffs_ads_cats(this,valences,CEC,log_act_coeffs)
            implicit none
            class(cat_exch_zone_c), intent(in) :: this !< Cation exchange object
            integer(kind=4), intent(in) :: valences(:) !< Valences of exchangeable cations (dim=num_exch_cats)
            REAL(kind=8), intent(in) :: CEC !< Cation exchange capacity [meq/g] or [mol/kg]
            real(kind=8), intent(out) :: log_act_coeffs(:) !< Log10 activity coefficients [-] (must be allocated previously, dim=num_exch_cats)
    
            integer(kind=4) :: i !< Loop index for exchangeable cations
        
            !< Loop through all exchangeable cations
            do i=1,this%num_exch_cats
                !< Compute log activity coefficient using selected convention
                call this%convention%compute_log_act_coeff_ads_cat(valences(i),CEC,log_act_coeffs(i))
            end do
        end subroutine
        
        !> \brief Compute Jacobian of log activity coefficients for adsorbed cations
        !> \details
        !>   Computes the Jacobian matrix of log10 activity coefficients with respect
        !>   to species concentrations. Required for Newton-Raphson iterations.
        !>   
        !>   **Matrix dimensions:**
        !>   - Rows: num_exch_cats (exchangeable cations)
        !>   - Columns: num_species (all species in system)
        !>   
        !>   **Current implementation:**
        !>   Sets all elements to zero (constant activity coefficients assumed).
        !>   
        !>   **Future implementation:**
        !>   Should compute: d(log γ_i)/d(c_j) for concentration-dependent models.
        !> \param[in,out] this Cation exchange object
        !> \param[in] log_act_coeffs Log10 activity coefficients [-] (dim=num_exch_cats)
        !> \param[out] log_Jacobian_act_coeffs Jacobian matrix [-] (must be allocated, dim=num_exch_cats x num_species)
        subroutine compute_log_Jacobian_act_coeffs_ads_cats(this,log_act_coeffs,log_Jacobian_act_coeffs)
            implicit none
            class(cat_exch_zone_c) :: this !< Cation exchange object
            real(kind=8), intent(in) :: log_act_coeffs(:) !< Log10 activity coefficients [-] (dim=num_exch_cats)
            real(kind=8), intent(out) :: log_Jacobian_act_coeffs(:,:) !< Jacobian matrix [-] (must be allocated, dim=num_exch_cats x num_species)
            !< Initialize Jacobian to zero (constant activity coefficients)
            log_Jacobian_act_coeffs=0d0
        end subroutine
    
        !> \brief Compute total number of surface complexes
        !> \details
        !>   Computes num_surf_compl from num_exch_cats.
        !>   
        !>   **Formula:**
        !>   num_surf_compl = num_exch_cats + 1
        !>   
        !>   The +1 accounts for the free site (unoccupied exchange site).
        !>   
        !>   **Example:**
        !>   - If 2 exchangeable cations (Na+, Ca2+)
        !>   - Then 3 surface complexes: X-Na, X2-Ca, X (free site)
        !> \param[in,out] this Cation exchange object
        subroutine compute_num_surf_compl(this)
            implicit none
            class(cat_exch_zone_c) :: this !< Cation exchange object
            !< Check if exchangeable cations exist
            if (this%num_exch_cats>0) then
                !< Total surface complexes = exchangeable cations + free site
                this%num_surf_compl=this%num_exch_cats+1
            end if
        end subroutine
      
        !> \brief Set number of exchangeable cations
        !> \details
        !>   Sets the count of exchangeable cations in the cation exchange zone.
        !>   
        !>   **Typical values:**
        !>   - Clay minerals: 2-3 (e.g., Na+, K+, Ca2+)
        !>   - Zeolites: 1-4
        !> \param[in,out] this Cation exchange object
        !> \param[in] num_exch_cats Number of exchangeable cations [-]
        subroutine set_num_exch_cats(this,num_exch_cats)
            implicit none
            class(cat_exch_zone_c) :: this !< Cation exchange object
            integer(kind=4), intent(in), optional :: num_exch_cats !< Number of exchangeable cations [-]
            !< Assign number of exchangeable cations
            this%num_exch_cats=num_exch_cats
        end subroutine
                
        !> \brief Allocate array for exchangeable cation indices
        !> \details
        !>   Allocates the exch_cat_indices array that maps exchangeable cations
        !>   to their positions in the aqueous phase species array.
        !>   
        !>   **Purpose:**
        !>   Links adsorbed cations to their aqueous counterparts for:
        !>   - Mass balance calculations
        !>   - Activity coefficient computations
        !>   - Jacobian assembly
        !> \param[in,out] this Cation exchange object
        !> \param[in] num_exch_cats Optional: number of exchangeable cations [-]
        subroutine allocate_exch_cat_indices(this,num_exch_cats)
            implicit none
            class(cat_exch_zone_c) :: this !< Cation exchange object
            integer(kind=4), intent(in), optional :: num_exch_cats !< Optional: number of exchangeable cations [-]
            !< Set number of exchangeable cations if provided
            if (present(num_exch_cats)) then
                this%num_exch_cats=num_exch_cats !< Update count
            end if
            !< Allocate indices array with current num_exch_cats
            allocate(this%exch_cat_indices(this%num_exch_cats))
        end subroutine
        
        !> \brief Set number of surface complexes
        !> \details
        !>   Sets the count of surface complexes in the surface phase.
        !>   
        !>   **Includes:**
        !>   - Occupied surface sites (adsorbed species)
        !>   - Free sites (unoccupied)
        !> \param[in,out] this Surface object
        !> \param[in] num_surf_compl Number of surface complexes [-]
        subroutine set_num_surf_compl(this,num_surf_compl)
            implicit none
            class(surface_c) :: this !< Surface object
            integer(kind=4), intent(in) :: num_surf_compl !< Number of surface complexes [-]
            !< Assign number of surface complexes
            this%num_surf_compl=num_surf_compl
        end subroutine
        
        !> \brief Set surface complexes array
        !> \details
        !>   Assigns the surf_compl array with provided surface complex data.
        !>   
        !>   **Validation:**
        !>   - If array allocated: checks size doesn't exceed num_surf_compl
        !>   - Stops execution if size mismatch detected
        !>   
        !>   **Safety:**
        !>   Prevents buffer overflow by validating array bounds.
        !> \param[in,out] this Surface object
        !> \param[in] surf_compl Array of surface complexes (solid_species_c type)
        subroutine set_surf_compl(this,surf_compl)
            implicit none
            class(surface_c) :: this !< Surface object
            class(solid_species_c), intent(in) :: surf_compl(:) !< Array of surface complexes
            
            !< Validate array size if already allocated
            if (allocated(this%surf_compl) .and. size(surf_compl)>this%num_surf_compl) then
                error stop "Number of surface complexes is wrong" !< Size exceeds limit
            else
                !< Assign surface complexes array
                this%surf_compl=surf_compl
            end if
        end subroutine
        
        !> \brief Allocate surface complexes array
        !> \details
        !>   Allocates memory for surf_compl array.
        !>   
        !>   **Workflow:**
        !>   1. Validate num_surf_compl (if provided)
        !>   2. Deallocate existing array (if allocated)
        !>   3. Allocate new array with current num_surf_compl
        !>   
        !>   **Validation:**
        !>   - num_surf_compl must be non-negative
        !>   - Stops execution if negative value provided
        !> \param[in,out] this Surface object
        !> \param[in] num_surf_compl Optional: number of surface complexes [-]
        subroutine allocate_surf_compl(this,num_surf_compl)
            implicit none
            class(surface_c) :: this !< Surface object
            integer(kind=4), intent(in), optional :: num_surf_compl !< Optional: number of surface complexes [-]
            
            !< Validate and set num_surf_compl if provided
            if (present(num_surf_compl) .and. num_surf_compl>=0) then
                this%num_surf_compl=num_surf_compl !< Valid non-negative value
            else if (present(num_surf_compl)) then
                error stop "Number of surface complexes must be non-negative" !< Negative value error
            else
                continue !< Use existing num_surf_compl
            end if
            !< Deallocate existing array if allocated
            if (allocated(this%surf_compl)) then
                deallocate(this%surf_compl) !< Free memory
            end if
            !< Allocate array with current size
            allocate(this%surf_compl(this%num_surf_compl))
        end subroutine
                
        !> \brief Set activity coefficient convention for cation exchange
        !> \details
        !>   Associates a pointer to the activity coefficient convention object.
        !>   
        !>   **Available conventions:**
        !>   - Gaines-Thomas: Rational activity coefficients
        !>   - Gapon: Activity proportional to equivalent fraction
        !>   - Vanselow: Activity proportional to mole fraction
        !>   
        !>   **Note:**
        !>   Convention object must remain in scope for lifetime of cat_exch_zone_c.
        !> \param[in,out] this Cation exchange object
        !> \param[in] convention Activity coefficient convention (target, must persist)
        subroutine set_convention(this,convention)
            implicit none
            class(cat_exch_zone_c) :: this !< Cation exchange object
            class(exch_sites_conv_c), intent(in), target :: convention !< Activity coefficient convention (must be target)
            !< Associate pointer to convention
            this%convention=>convention
        end subroutine
        
        !> \brief Check if surface complex exists in surface phase
        !> \details
        !>   Searches for a surface complex by name in the surf_compl array.
        !>   
        !>   **Search algorithm:**
        !>   1. Initialize flag to FALSE and index to 0
        !>   2. Loop through all surface complexes
        !>   3. Compare names (case-sensitive)
        !>   4. If match found: set flag=TRUE, store index, exit loop
        !>   
        !>   **Use cases:**
        !>   - Validate surface complex exists before operations
        !>   - Get index for array access
        !>   - Check reactive zone composition
        !> \param[in] this Surface object
        !> \param[in] surf_compl Surface complex to search for
        !> \param[out] flag TRUE if surface complex found, FALSE otherwise
        !> \param[out] surf_compl_ind Optional: index in surf_compl array (0 if not found)
        subroutine is_surf_compl_in(this,surf_compl,flag,surf_compl_ind)
            implicit none
            class(surface_c), intent(in) :: this !< Surface object
            class(solid_species_c), intent(in) :: surf_compl !< Surface complex to search for
            logical, intent(out) :: flag !< TRUE if found, FALSE otherwise
            integer(kind=4), intent(out), optional :: surf_compl_ind !< Optional: index in array (0 if not found)
            
            integer(kind=4) :: i !< Loop index
            
            !< Initialize search result to "not found"
            flag=.false.
            !< Initialize index to zero if requested
            if (present(surf_compl_ind)) then
                surf_compl_ind=0 !< Not found by default
            end if
            !< Loop through all surface complexes
            do i=1,this%num_surf_compl
                !< Compare names (case-sensitive)
                if (surf_compl%name==this%surf_compl(i)%name) then
                    flag=.true. !< Mark as found
                    !< Store index if requested
                    if (present(surf_compl_ind)) then
                        surf_compl_ind=i !< Position in array
                    end if
                    exit !< Stop searching (found)
                end if
            end do
        end subroutine
end module