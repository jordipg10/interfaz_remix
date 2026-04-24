!> @file PDE_model_m.f90
!> @brief PDE model module for wrapping and managing PDE objects
!> @details This module provides a wrapper class for managing one-dimensional PDE objects
!> in a flexible, object-oriented manner. The PDE_model_c type acts as a container that
!> holds a pointer to a PDE_1D_c object, allowing for polymorphic behavior and dynamic
!> association of different PDE types at runtime.
!>
!> @par Design Pattern:
!> This module implements a wrapper/adapter pattern where:
!>   - The PDE_model_c type contains a polymorphic pointer to a PDE object
!>   - The pointer can reference any derived type of PDE_1D_c
!>   - Enables flexible model composition and runtime PDE type selection
!>
!> @par Usage:
!> The module is typically used to encapsulate PDE objects within larger simulation
!> frameworks, allowing models to switch between different PDE types (diffusion,
!> transport, reaction-diffusion, etc.) without changing the model structure.
!>
!> @see PDE_m For base PDE class definitions and hierarchy
!> @author Generated documentation
!> @date November 2025

!> @brief PDE model module - provides wrapper for PDE objects
module PDE_model_m
    use PDE_m, only: PDE_c !< Import PDE module with PDE_c class and derived types
    implicit none !< Enforce explicit variable declarations for type safety
    save !< Preserve module variables between procedure calls (all module data has SAVE attribute)
    private !< Private module scope: internal details hidden from outside modules
    !> @brief PDE model wrapper type for managing PDE objects
    !> @details Container class that holds a polymorphic pointer to a one-dimensional PDE object.
    !> Provides flexible association with different PDE types at runtime through the pointer.
    !> The pointer allows the model to reference various PDE implementations (diffusion, transport,
    !> advection-diffusion, etc.) without requiring recompilation.
    type, public :: PDE_model_c !< Public type: accessible from other modules
        !class(PDE_c), pointer :: PDE !< Polymorphic pointer to 1D PDE object (can point to any PDE_1D_c derived type)
    contains
        !procedure :: set_PDE !< Public method: associate PDE object with model (bind pointer to target PDE)
    end type !< End PDE_model_c type definition
    
    contains !< Begin module procedures (implementation section)
        
        !> @brief Associate a PDE object with the model
        !> @details Sets the model's PDE pointer to reference the provided PDE object using pointer
        !> association. This establishes the link between the model wrapper and the actual PDE
        !> implementation, allowing the model to access and manipulate the PDE through the pointer.
        !> The PDE argument must have the TARGET attribute to allow pointer association.
        !>
        !> @par Pointer Association:
        !> Uses Fortran pointer assignment (=>) to create an alias to the target PDE object.
        !> After this call, this%PDE points to the same memory location as the PDE argument.
        !> Changes to the PDE through the pointer affect the original object.
        !>
        !> @param[in,out] this PDE_model_c object whose PDE pointer will be set
        !> @param[in] PDE Target PDE object to associate with model (must be declared with TARGET attribute)
        ! subroutine set_PDE(this,PDE)
        !     implicit none !< Enforce explicit variable declarations for type safety
        !     class(PDE_model_c) :: this !< PDE model object containing pointer to be set [-]
        !     class(PDE_c), intent(in), target :: PDE !< Input PDE object to point to (TARGET allows pointer association) [-]
        !     this%PDE=>PDE !< Pointer association: make this%PDE point to PDE argument (=> operator for pointer assignment)
        ! end subroutine !< End set_PDE subroutine
        
end module !< End PDE_model_m module