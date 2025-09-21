module target_m
    !> Module for target class
    !> This module defines the target class which contains information about the target coordinates.
    !> The target can be a cell or a node in the mesh.
    !> FALTA ACLARAR LO DEL SUBDOMAIN
    !use subdomain_m, only: subdomain_c
    implicit none
    save
    type, public :: target_c
        real(kind=8), allocatable :: coord(:) !> space coordinates
        real(kind=8), allocatable :: coord_D(:) !> dimensionless space coordinates
        !real(kind=8) :: y !> y coordinate (optional, for 2D)
        !real(kind=8) :: z !> z coordinate (optional, for 3D)
        integer(kind=4) :: id !> target ID (must be non-negative)
        logical :: is_boundary !> TRUE if target is a boundary target, FALSE otherwise
        !class(subdomain_c), pointer :: subdomain !> pointer to the subdomain class
        real(kind=8) :: measure !> measure of the target (length in 1D, area in 2D, volume in 3D)
        real(kind=8) :: thickness !> thickness of the target (optional, for 2D or 3D)
    contains
        procedure, public :: set_coordinates
        procedure, public :: set_id
        procedure, public :: set_boundary_flag
        !procedure, public :: set_subdomain
        procedure, public :: set_measure
        procedure, public :: set_thickness
        procedure, public :: compute_dimless_coords
    end type target_c

    contains

    subroutine set_coordinates(this,coord)
        implicit none
        class(target_c) :: this
        real(kind=8), intent(in) :: coord(:) !> space coordinates
        !real(kind=8), intent(in) :: x, y, z !> coordinates
        this%coord = coord
        !this%y = y
        !this%z = z
    end subroutine set_coordinates

    subroutine set_id(this, id)
        implicit none
        class(target_c) :: this
        integer(kind=4), intent(in) :: id !> target ID
        this%id = id
    end subroutine set_id

    subroutine set_boundary_flag(this, is_boundary)
        implicit none
        class(target_c) :: this
        logical, intent(in) :: is_boundary !> TRUE if target is a boundary target, FALSE otherwise
        this%is_boundary = is_boundary
    end subroutine set_boundary_flag

    !subroutine set_subdomain(this, subdomain)
    !    implicit none
    !    class(target_c) :: this
    !    type(subdomain_c), intent(in), target :: subdomain !> pointer to the subdomain class
    !    this%subdomain => subdomain
    !end subroutine set_subdomain

    subroutine set_measure(this, measure)
        implicit none
        class(target_c) :: this
        real(kind=8), intent(in), optional :: measure !> measure of the target
        if (present(measure)) then
            this%measure = measure
        else
            this%measure = 0d0 !> Default measure is 0
        end if
    end subroutine set_measure
    
    subroutine set_thickness(this,thickness)
        implicit none
        class(target_c) :: this
        real(kind=8), intent(in), optional :: thickness !> thickness of the target
        if (present(thickness)) then
            this%thickness = thickness
        else
            this%thickness = 1d0 !> Default thickness is 1
        end if
    end subroutine set_thickness
    
    subroutine compute_dimless_coords(this, L_c)
        implicit none
        class(target_c) :: this
        real(kind=8), intent(in) :: L_c !> characteristic length for dimensionless coordinates
        if (L_c <= 0d0) error stop "Characteristic length must be positive"
        allocate(this%coord_D(size(this%coord)))
        this%coord_D = this%coord / L_c !> Compute dimensionless coordinates
    end subroutine compute_dimless_coords
end module target_m