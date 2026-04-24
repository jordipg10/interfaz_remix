!> \file conc_m.f90
!> \brief Concentration module for diffusion and transport problems
!> \details Defines the concentration class and related procedures for handling concentrations, boundary conditions, 
!> recharge, and external values in chemical transport simulations. This module provides a structured approach to
!> managing different types of concentration data needed in reactive transport modeling:
!> - Internal domain concentrations (conc)
!> - External/source concentrations (conc_ext)
!> - Boundary concentrations for Dirichlet conditions (conc_bd)
!> - Recharge concentrations and flags (conc_r, conc_r_flag)
!>
!> The module supports various boundary condition types and recharge scenarios commonly encountered in
!> subsurface transport modeling, including:
!> - Dirichlet boundary conditions (specified concentrations at boundaries)
!> - External sources/sinks (pumping wells, injection, discharge)
!> - Recharge zones (infiltration with associated solute concentrations)

module conc_m
    implicit none !> Enforce explicit variable declarations
    save !> Preserve module variables between procedure calls
    private !> Default accessibility is private
    !> \class conc_c
    !> \brief Concentration class for diffusion or transport problems
    !> \details Stores all concentration-related data for a chemical transport problem, including:
    !> - Domain concentrations (internal solution field)
    !> - External concentrations (sources/sinks outside domain)
    !> - Boundary concentrations (Dirichlet boundary conditions)
    !> - Recharge concentrations (infiltration/percolation sources)
    !> - Recharge flags (indicator for active recharge zones)
    !>
    !> This class provides setter methods to initialize and update each concentration array,
    !> enabling flexible management of different concentration types in reactive transport simulations.
    !> All concentration arrays are allocatable to accommodate variable problem sizes.
    
    type, public :: conc_c
        !> \var conc
        !> \brief Primary concentration array for the domain (c)
        !> \details Stores concentration values at spatial nodes within the computational domain.
        !> This is the main solution field for the transport equation. Units: [M/L³]
        real(kind=8), allocatable :: conc(:)
        
        !> \var conc_ext
        !> \brief External concentration array (c_e)
        !> \details Stores concentration values for external sources or sinks (e.g., pumping wells,
        !> injection points, discharge locations). These represent concentrations outside the primary
        !> computational domain but affecting the solution through source/sink terms. Units: [M/L³]
        real(kind=8), allocatable :: conc_ext(:)
        
        !> \var conc_bd
        !> \brief Boundary concentration array (c_b)
        !> \details Stores concentration values prescribed at Dirichlet (fixed concentration) boundaries.
        !> These values enforce known concentrations at domain boundaries. Units: [M/L³]
        real(kind=8), allocatable :: conc_bd(:)
        
        !> \var conc_r
        !> \brief Recharge concentration array (c_r)
        !> \details Stores concentration values associated with recharge (infiltration, percolation).
        !> Used when water enters the domain from above or other recharge zones with a specified
        !> concentration. Units: [M/L³]
        real(kind=8), allocatable :: conc_r(:)
        
        !> \var conc_r_flag
        !> \brief Recharge flag array
        !> \details Integer flags indicating active recharge zones: 1 if recharge rate r > 0 (active),
        !> 0 otherwise (inactive). Used to determine where recharge concentrations should be applied.
        !> Dimensionless [-]
        integer(kind=4), allocatable :: conc_r_flag(:)
    contains
        !> \brief Set concentration array
        !> \details Assigns values to the primary domain concentration array
        procedure :: set_conc
        
        !> \brief Set external concentration array
        !> \details Assigns values to external source/sink concentrations
        procedure :: set_conc_ext
        
        !> \brief Set boundary concentration array
        !> \details Assigns values to Dirichlet boundary concentrations
        procedure :: set_conc_bd
        
        !> \brief Set recharge concentration and flag arrays
        !> \details Assigns values to recharge concentrations and associated flags
        procedure :: set_conc_r
    end type conc_c
    
    contains
    
    !> \brief Set the concentration array
    !> \details Assigns values to the primary domain concentration array (this%conc).
    !> This subroutine copies the input concentration array to the object's internal storage.
    !> The array must be allocated and sized appropriately before calling this subroutine.
    !> Typically used to initialize concentrations or update them during time stepping.
    !>
    !> \param[in,out] this Concentration object to modify
    !> \param[in] conc Concentration array to set [M/L³]
    
    subroutine set_conc(this, conc)
        implicit none !> Enforce explicit variable declarations
        class(conc_c) :: this !> Concentration object (modified in place)
        real(kind=8), intent(in) :: conc(:) !> Input concentration array to copy [M/L³]
        this%conc = conc !> Assign input array to object's concentration array [M/L³]
    end subroutine set_conc
    
    !> \brief Set the external concentration array
    !> \details Assigns values to the external concentration array (this%conc_ext).
    !> External concentrations represent solute concentrations in sources or sinks outside
    !> the primary computational domain, such as:
    !> - Pumping well concentrations
    !> - Injection well concentrations
    !> - Discharge point concentrations
    !> - River or stream concentrations (for river-aquifer interaction)
    !>
    !> \param[in,out] this Concentration object to modify
    !> \param[in] conc_ext External concentration array to set [M/L³]
    
    subroutine set_conc_ext(this, conc_ext)
        implicit none !> Enforce explicit variable declarations
        class(conc_c) :: this !> Concentration object (modified in place)
        real(kind=8), intent(in) :: conc_ext(:) !> Input external concentration array to copy [M/L³]
        this%conc_ext = conc_ext !> Assign input array to object's external concentration array [M/L³]
    end subroutine set_conc_ext
    
    !> \brief Set the boundary concentration array
    !> \details Assigns values to the boundary concentration array (this%conc_bd).
    !> Boundary concentrations are used to enforce Dirichlet (first-type) boundary conditions,
    !> where concentrations are specified at domain boundaries. Common applications include:
    !> - Fixed concentration at inflow boundaries
    !> - Specified concentration at constant head boundaries
    !> - Known concentration at domain edges (e.g., contact with a reservoir)
    !>
    !> These values remain constant or can be updated for time-dependent boundary conditions.
    !>
    !> \param[in,out] this Concentration object to modify
    !> \param[in] conc_bd Boundary concentration array to set [M/L³]
    
    subroutine set_conc_bd(this, conc_bd)
        implicit none !> Enforce explicit variable declarations
        class(conc_c) :: this !> Concentration object (modified in place)
        real(kind=8), intent(in) :: conc_bd(:) !> Input boundary concentration array to copy [M/L³]
        this%conc_bd = conc_bd !> Assign input array to object's boundary concentration array [M/L³]
    end subroutine set_conc_bd
    
    !> \brief Set the recharge concentration and flag arrays
    !> \details Assigns values to both the recharge concentration array (this%conc_r) and
    !> the recharge flag array (this%conc_r_flag). Recharge represents infiltration of water
    !> (and associated solutes) into the domain, typically from:
    !> - Precipitation/rainfall infiltration
    !> - Irrigation return flow
    !> - Surface water percolation
    !> - Artificial recharge operations
    !>
    !> The flag array indicates active recharge zones:
    !> - conc_r_flag = 1: Active recharge (recharge rate r > 0), use conc_r
    !> - conc_r_flag = 0: No recharge (r = 0), ignore conc_r
    !>
    !> This allows for spatially variable recharge patterns where only certain zones
    !> receive recharge with associated solute input.
    !>
    !> \param[in,out] this Concentration object to modify
    !> \param[in] conc_r Recharge concentration array to set [M/L³]
    !> \param[in] conc_r_flag Recharge flag array to set: 1 if r>0, 0 otherwise [-]
    
    subroutine set_conc_r(this, conc_r, conc_r_flag)
        implicit none !> Enforce explicit variable declarations
        class(conc_c) :: this !> Concentration object (modified in place)
        real(kind=8), intent(in) :: conc_r(:) !> Input recharge concentration array to copy [M/L³]
        integer(kind=4), intent(in) :: conc_r_flag(:) !> Input recharge flag array to copy: 1 if recharge active, 0 otherwise [-]
        this%conc_r = conc_r !> Assign input array to object's recharge concentration array [M/L³]
        this%conc_r_flag = conc_r_flag !> Assign input flags to object's recharge flag array [-]
    end subroutine set_conc_r
    
end module conc_m !> End of concentration module