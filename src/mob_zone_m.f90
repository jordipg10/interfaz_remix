!> @file mob_zone_m.f90
!> @brief Mobile zone module for Multi-Rate Mass Transfer (MRMT) modeling
!> @details This module defines the mobile zone type and methods for managing mobile zone
!> properties in MRMT reactive transport simulations. The mobile zone represents the
!> advective-dispersive domain where solute transport occurs, as opposed to immobile zones
!> where mass transfer is diffusion-limited.
!>
!> @par Mobile Zone Characteristics:
!> - Advection and dispersion occur in mobile zone
!> - Exchange with multiple immobile zones via mass transfer
!> - Mobile porosity φ_m defines fraction of pore space in mobile domain
!> - Concentrations tracked spatially (discretization points) and for multiple species
!>
!> @see MRMT_m Multi-Rate Mass Transfer module
!> @see imm_zone_m Immobile zone module
!> @author Generated documentation
!> @date November 2025

!> @brief Mobile zone module for MRMT transport modeling
module mob_zone_m
    use imm_zone_m, only: imm_zone_c !< Import only immobile zone type for pointer array
    implicit none !< Enforce explicit variable declarations for type safety
    private !< Make all entities private by default
    
    !> @brief Wrapper type to hold pointer to immobile zone
    !> @details This wrapper type is necessary to create an array of pointers in Fortran.
    !> Each wrapper contains a pointer to an immobile zone object, allowing dynamic
    !> association of multiple immobile zones with a mobile zone.
    type, public :: imm_zone_ptr
        class(imm_zone_c), pointer :: ptr => null() !< Pointer to immobile zone object [-]
    end type
    
    !> @brief Mobile zone class for MRMT reactive transport
    !> @details Defines mobile zone properties including concentrations (spatial and species),
    !> mobile porosity, and methods for initialization and configuration.
    type, public :: mob_zone_c
        !class(spatial_discr_c), pointer :: spatial_discr !< Commented: Spatial discretization object for mobile zone (future feature) [-]
        !class(aq_phase_c), pointer :: aq_phase !< Commented: Aqueous phase object for mobile zone (future feature) [-]
        !integer(kind=8) :: num_mob_species !< Commented: Number of mobile species (redundant, use array size) [-]
        
        real(kind=8), allocatable :: conc_init(:) !< Initial concentration array: conc_init(species, spatial_point) [M/L³]
        real(kind=8), allocatable :: conc(:) !< Current concentration array: conc(species, spatial_point) [M/L³]
        real(kind=8) :: mob_por !< Mobile zone porosity φ_m (fraction of pore space in mobile domain, 0 < φ_m < 1) [-]
        type(imm_zone_ptr), allocatable :: imm_zones(:) !< Array of pointers to immobile zone objects associated with this mobile zone [-]
    contains
        procedure :: set_conc_init !< Setter method: initialize concentration arrays for mobile zone
        procedure :: set_mob_por !< Setter method: set mobile porosity with validation
        procedure :: set_imm_zones !< Setter method: set array of immobile zones with validation
    end type
    
contains
    
    !> @brief Set initial concentration in mobile zone
    !> @details Initializes the mobile zone concentration array with validation.
    !> Ensures that the provided concentration array has positive dimensions and
    !> non-negative concentration values (physical requirement). Also initializes
    !> the current concentration array to the same values.
    !>
    !> @param[in,out] this Mobile zone object being initialized
    !> @param[in] conc_init Initial concentration array to assign [M/L³]
    !>
    !> @par Array Dimensions:
    !> - Dimension 1: Number of species in mobile zone [-]
    !> - Dimension 2: Number of spatial discretization points [-]
    !>
    !> @par Physical Constraints:
    !> - Array dimensions must be positive (at least one species and one point)
    !> - All concentration values must be non-negative (physical requirement)
    !> - Current concentration initialized to match initial concentration
    !>
    !> @note This replaces the original incomplete validation logic
    subroutine set_conc_init(this,conc_init)
        implicit none !< Enforce explicit variable declarations for type safety
        class(mob_zone_c), intent(inout) :: this !< Mobile zone object to be initialized [-]
        real(kind=8), intent(in) :: conc_init(:) !< Input initial concentration array: conc_init(species, spatial_point) [M/L³]
        
        !> Validate array dimensions and concentration values
        if (size(conc_init)<=0) then
            error stop "Initial concentration array must have positive dimensions" !< Halt if dimensions invalid [-]
        else if (any(conc_init<0d0)) then
            error stop "Initial concentrations must be non-negative" !< Halt if negative concentrations found [-]
        end if
        
        this%conc_init=conc_init !< Assign validated initial concentration array [M/L³]
        this%conc=conc_init !< Initialize current concentration to initial values (direct assignment) [M/L³]
    end subroutine set_conc_init

    !> @brief Set mobile porosity
    !> @details Sets the mobile zone porosity with physical validation.
    !> Mobile porosity φ_m represents the fraction of total pore space that participates
    !> in advective-dispersive transport (mobile domain). Must be positive and typically
    !> less than total porosity.
    !>
    !> @param[in,out] this Mobile zone object being configured
    !> @param[in] mob_por Mobile porosity value to assign [-]
    !>
    !> @par Physical Constraints:
    !> - Mobile porosity must be positive: φ_m > 0
    !> - Typically: 0 < φ_m < φ_total (less than total porosity)
    !> - Represents fraction of pore space in mobile domain
    !>
    !> @par MRMT Context:
    !> In Multi-Rate Mass Transfer models, total porosity is partitioned:
    !> φ_total = φ_m + Σ φ_im,i
    !> where φ_im,i are immobile zone porosities
    elemental subroutine set_mob_por(this,mob_por)
        implicit none !< Enforce explicit variable declarations for type safety
        class(mob_zone_c), intent(inout) :: this !< Mobile zone object to be configured [-]
        real(kind=8), intent(in) :: mob_por !< Input mobile porosity value (must be positive) [-]
        
        !> Validate that mobile porosity is physically meaningful (positive value)
        if (mob_por<=0d0) error stop "Mobile porosity must be positive" !< Halt if non-positive porosity [-]
        this%mob_por=mob_por !< Assign validated mobile porosity value [-]
    end subroutine set_mob_por
    
    !> @brief Set array of immobile zones
    !> @details Allocates and initializes the array of immobile zone pointers with validation.
    !> This method sets up the association between a mobile zone and multiple immobile zones
    !> for Multi-Rate Mass Transfer (MRMT) modeling.
    !>
    !> @param[in,out] this Mobile zone object being configured
    !> @param[in] n_imm_zones Number of immobile zones to allocate [-]
    !>
    !> @par Physical Constraints:
    !> - Number of immobile zones must be positive: n_imm_zones > 0
    !> - Typical MRMT models use 1-10 immobile zones
    !> - Each immobile zone will have different exchange rates (multi-rate behavior)
    !>
    !> @par MRMT Context:
    !> Multiple immobile zones with different exchange rates α_i create the multi-rate
    !> mass transfer behavior. Common distributions:
    !> - Uniform spacing: α_i = α_min * (α_max/α_min)^((i-1)/(n-1))
    !> - Log-uniform spacing in exchange rate space
    !>
    !> @note This only allocates the pointer array. Individual zones must be allocated
    !> and initialized separately using allocate(this%imm_zones(i)%ptr)
    subroutine set_imm_zones(this,n_imm_zones)
        implicit none !< Enforce explicit variable declarations for type safety
        class(mob_zone_c), intent(inout) :: this !< Mobile zone object to be configured [-]
        integer(kind=4), intent(in) :: n_imm_zones !< Number of immobile zones (must be positive) [-]
        
        !> Validate that number of immobile zones is physically meaningful (positive)
        if (n_imm_zones<=0) error stop "Number of immobile zones must be positive" !< Halt if non-positive count [-]
        
        !> Deallocate and reallocate array (automatic deallocation on assignment in F2003+)
        if (allocated(this%imm_zones)) deallocate(this%imm_zones) !< Free previously allocated immobile zones array [-]
        allocate(this%imm_zones(n_imm_zones)) !< Allocate pointer wrapper array for immobile zones [-]
        
        !> Pointers initialized to null by default (=> null() in type definition)
        !> Individual zones must be allocated separately: allocate(this%imm_zones(i)%ptr)
    end subroutine set_imm_zones
    
end module