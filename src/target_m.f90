!> \file target_m.f90
!> \brief Target module for spatial discretization points
!> \details
!>   This module defines the target class representing spatial points
!>   (cells or nodes) in the computational mesh.
!>   
!>   **Key features:**
!>   - Physical and dimensionless coordinates
!>   - Boundary flag for boundary conditions
!>   - Geometric measures (length, area, volume)
!>   - Thickness for 2D/3D applications
!>   
!>   **Applications:**
!>   - Finite difference/volume discretization
!>   - Reactive transport mesh management
!>   - Boundary condition specification
!>   - Spatial integration
!>   
!>   **Target types:**
!>   - Interior targets: Within domain
!>   - Boundary targets: At domain edges (Dirichlet, Neumann, etc.)
!>   
!>   **TODO:**
!>   - Clarify subdomain integration
!>   
!> \author Jordi Petchamé-Guerrero
!> \date October 2025
module target_m
    !> Module for target class
    !> This module defines the target class which contains information about the target coordinates.
    !> The target can be a cell or a node in the mesh.
    !> FALTA ACLARAR LO DEL SUBDOMAIN
    !use subdomain_m, only: subdomain_c
    implicit none !< Enforce explicit variable declarations
    save !< Preserve module state between calls
    private !< Default accessibility is private
    !> \brief Target type for spatial discretization points
    !> \details
    !>   Defines a spatial point (cell or node) in the computational mesh.
    !>   
    !>   **Coordinate systems:**
    !>   - Physical coordinates: coord (m)
    !>   - Dimensionless coordinates: coord_D = coord/L_c (normalized by characteristic length)
    !>   
    !>   **Geometric properties:**
    !>   - measure: Length (1D), Area (2D), or Volume (3D)
    !>   - thickness: For 2D/3D cross-sectional areas
    !>   
    !>   **Dimensionality:**
    !>   - 1D: coord(1) = x, measure = Δx
    !>   - 2D: coord(1:2) = [x,y], measure = Δx·Δy
    !>   - 3D: coord(1:3) = [x,y,z], measure = Δx·Δy·Δz
    !>   
    !>   **Boundary targets:**
    !>   - is_boundary = TRUE for domain boundaries
    !>   - Used to apply Dirichlet, Neumann, or mixed BCs
    !>   
    !>   **Usage example:**
    !>   ```fortran
    !>   type(target_c) :: cell
    !>   call cell%set_id(10)
    !>   call cell%set_coordinates([0.5d0, 1.0d0])  ! 2D point
    !>   call cell%set_measure(0.01d0)              ! Area = 0.01 m²
    !>   call cell%compute_dimless_coords(100.0d0)  ! L_c = 100 m
    !>   ```
    type, public :: target_c !< Target for spatial discretization
        real(kind=8), allocatable :: coord(:) !< Physical space coordinates [m] (size = ndim)
        real(kind=8), allocatable :: coord_D(:) !< Dimensionless coordinates [-] (coord/L_c, size = ndim)
        !real(kind=8) :: y !> y coordinate (optional, for 2D)
        !real(kind=8) :: z !> z coordinate (optional, for 3D)
        integer(kind=4) :: id !< Target ID (must be non-negative) [-]
        logical :: is_boundary !< TRUE if boundary target, FALSE if interior [-]
        !class(subdomain_c), pointer :: subdomain !> pointer to the subdomain class
        real(kind=8) :: measure !< Measure: length [m] (1D), area [m²] (2D), volume [m³] (3D)
        real(kind=8) :: thickness !< Thickness [m] (optional, for 2D or 3D cross-sections)
    contains
        procedure :: set_coordinates !< Set physical coordinates
        procedure :: set_id !< Set target identifier
        procedure :: set_boundary_flag !< Set boundary condition flag
        !procedure :: set_subdomain !< Set subdomain pointer (future)
        procedure :: set_measure !< Set geometric measure
        procedure :: set_thickness !< Set target thickness
        procedure :: compute_dimless_coords !< Compute dimensionless coordinates
    end type target_c

    contains

    !> \brief Set physical space coordinates
    !> \details
    !>   Assigns physical coordinates to the target point.
    !>   
    !>   **Dimensionality:**
    !>   - 1D: coord = [x]
    !>   - 2D: coord = [x, y]
    !>   - 3D: coord = [x, y, z]
    !>   
    !>   **Note:**
    !>   coord array must be allocated before calling.
    !>   Typically size(coord) = 1, 2, or 3.
    !> \param[in,out] this Target object
    !> \param[in] coord Physical space coordinates [m] (size = ndim)
    subroutine set_coordinates(this,coord)
        implicit none
        class(target_c) :: this !< Target object to modify
        real(kind=8), intent(in) :: coord(:) !< Physical space coordinates [m]
        !real(kind=8), intent(in) :: x, y, z !> coordinates
        !< Assign coordinates array
        this%coord = coord
        !this%y = y
        !this%z = z
    end subroutine set_coordinates

    !> \brief Set target identifier
    !> \details
    !>   Assigns a unique integer ID to the target.
    !>   
    !>   **ID convention:**
    !>   - Must be non-negative (id >= 0)
    !>   - Should be unique within the mesh
    !>   - Used for indexing in arrays
    !>   
    !>   **Typical usage:**
    !>   - Sequential numbering: 1, 2, 3, ..., n_targets
    !>   - ID = 0 may be reserved for special purposes
    !> \param[in,out] this Target object
    !> \param[in] id Target identifier (non-negative) [-]
    subroutine set_id(this, id)
        implicit none
        class(target_c) :: this !< Target object to modify
        integer(kind=4), intent(in) :: id !< Target ID (non-negative) [-]
        !< Assign target identifier
        this%id = id
    end subroutine set_id

    !> \brief Set boundary condition flag
    !> \details
    !>   Marks the target as a boundary or interior point.
    !>   
    !>   **Flag values:**
    !>   - TRUE: Target is on domain boundary (apply BCs)
    !>   - FALSE: Target is in domain interior (solve PDEs)
    !>   
    !>   **Boundary types:**
    !>   - Dirichlet: Fixed concentration/pressure
    !>   - Neumann: Fixed flux
    !>   - Robin/Mixed: Combination of both
    !>   
    !>   **Usage:**
    !>   Essential for applying boundary conditions in solvers.
    !> \param[in,out] this Target object
    !> \param[in] is_boundary TRUE if boundary target, FALSE if interior [-]
    subroutine set_boundary_flag(this, is_boundary)
        implicit none
        class(target_c) :: this !< Target object to modify
        logical, intent(in) :: is_boundary !< TRUE if boundary, FALSE if interior [-]
        !< Assign boundary flag
        this%is_boundary = is_boundary
    end subroutine set_boundary_flag

    !subroutine set_subdomain(this, subdomain)
    !    implicit none
    !    class(target_c) :: this
    !    type(subdomain_c), intent(in), target :: subdomain !> pointer to the subdomain class
    !    this%subdomain => subdomain
    !end subroutine set_subdomain

    !> \brief Set geometric measure of target
    !> \details
    !>   Sets the measure (size) of the target depending on dimensionality.
    !>   
    !>   **Measure definition:**
    !>   - 1D: Length Δx [m]
    !>   - 2D: Area Δx·Δy [m²]
    !>   - 3D: Volume Δx·Δy·Δz [m³]
    !>   
    !>   **Usage in simulations:**
    !>   - Mass balance calculations
    !>   - Volume-averaged concentrations
    !>   - Integration over domain
    !>   
    !>   **Default:**
    !>   If not provided, measure = 0 (must be set explicitly for valid calculations).
    !> \param[in,out] this Target object
    !> \param[in] measure Optional: geometric measure [m], [m²], or [m³]
    subroutine set_measure(this, measure)
        implicit none
        class(target_c) :: this !< Target object to modify
        real(kind=8), intent(in), optional :: measure !< Optional: measure [m], [m²], or [m³]
        !< Check if measure provided
        if (present(measure)) then
            this%measure = measure !< Assign provided measure
        else
            this%measure = 0d0 !< Default measure is 0 (must be set explicitly)
        end if
    end subroutine set_measure
    
    !> \brief Set target thickness
    !> \details
    !>   Sets the thickness for 2D or 3D cross-sectional calculations.
    !>   
    !>   **Applications:**
    !>   - 2D models: Thickness perpendicular to plane (e.g., aquifer thickness)
    !>   - 3D models: Layer thickness in stratified systems
    !>   - Flux calculations: Cross-sectional area = length × thickness
    !>   
    !>   **Default:**
    !>   If not provided, thickness = 1.0 m (unit thickness).
    !>   
    !>   **Example:**
    !>   In 2D aquifer model, thickness represents aquifer depth.
    !> \param[in,out] this Target object
    !> \param[in] thickness Optional: target thickness [m]
    subroutine set_thickness(this,thickness)
        implicit none
        class(target_c) :: this !< Target object to modify
        real(kind=8), intent(in), optional :: thickness !< Optional: thickness [m]
        !< Check if thickness provided
        if (present(thickness)) then
            this%thickness = thickness !< Assign provided thickness
        else
            this%thickness = 1d0 !< Default thickness is 1 m (unit thickness)
        end if
    end subroutine set_thickness
    
    !> \brief Compute dimensionless coordinates
    !> \details
    !>   Computes dimensionless coordinates by normalizing with characteristic length.
    !>   
    !>   **Formula:**
    !>   coord_D = coord / L_c
    !>   
    !>   **Where:**
    !>   - coord: Physical coordinates [m]
    !>   - L_c: Characteristic length scale [m]
    !>   - coord_D: Dimensionless coordinates [-]
    !>   
    !>   **Characteristic length examples:**
    !>   - Domain length: L_c = L_domain
    !>   - Penetration depth: L_c = sqrt(D·t)
    !>   - Reaction length: L_c = v/k
    !>   
    !>   **Benefits of dimensionless coordinates:**
    !>   - Scale-independent analysis
    !>   - Improved numerical conditioning
    !>   - Universal solution curves
    !>   
    !>   **Validation:**
    !>   L_c must be positive; stops execution if L_c <= 0.
    !> \param[in,out] this Target object
    !> \param[in] L_c Characteristic length for normalization [m] (must be > 0)
    subroutine compute_dimless_coords(this, L_c)
        implicit none
        class(target_c) :: this !< Target object to modify
        real(kind=8), intent(in) :: L_c !< Characteristic length [m] (must be positive)
        !< Validate characteristic length is positive
        if (L_c <= 0d0) error stop "Characteristic length must be positive"
        !< Allocate dimensionless coordinates array with same size as physical coordinates
        allocate(this%coord_D(size(this%coord)))
        !< Compute dimensionless coordinates: coord_D = coord / L_c
        this%coord_D = this%coord / L_c
    end subroutine compute_dimless_coords
end module target_m