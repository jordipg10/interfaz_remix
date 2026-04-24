!> \file stability_parameters_m.f90
!> \brief Defines stability parameter classes and routines for time step control in numerical simulations
!> \details 
!>   This module provides an abstract superclass for stability parameters used to ensure
!>   numerical stability in time-dependent simulations (diffusion, advection, reaction-transport).
!>
!>   Key concepts:
!>   - Critical time step (Δt_crit): Maximum time step for stable explicit schemes
!>   - Dimensionless critical time step: Δt_crit / t_c (characteristic time scaling)
!>   - CFL (Courant-Friedrichs-Lewy) condition for advection
!>   - Diffusion stability criterion (Δt ≤ phi * Δx² / (2D) for explicit schemes)
!>   - von Neumann stability analysis for various numerical schemes
!>
!>   Abstract class structure:
!>   - stab_params_c: Base class for all stability parameter types
!>   - Subclasses implement compute_stab_params() for specific physics:
!>     * diff_stab_params_m: Diffusion stability
!>     * stab_params_flow_m: Advective flow stability
!>     * transport_stab_params_m: Combined advection-diffusion-reaction
!>
!>   Dependencies:
!>   - time_discr_m: Time discretization (time steps, schemes)
!>   - properties_m: Physical properties (diffusivity, velocity, etc.)
!>   - spatial_discr_m: Mesh/grid (via imports in interface)
!>
!> \author jordi Petchamé-Guerrero
!> \date 2025
module stability_parameters_m
    use time_discr_m, only: time_discr_c  !< Time discretization routines and types
    use properties_m, only: props_c  !< Physical property definitions
    use spatial_discr_m, only: spatial_discr_c  !< Spatial discretization (mesh/grid)
    implicit none
    save                  !< Preserve module variables between calls
    private               !< Private module scope
    !> \brief Abstract base class for stability parameters
    !> \details 
    !>   Defines common interface for computing stability constraints on time stepping.
    !>   
    !>   Member variables:
    !>   - Delta_t_crit: Critical time step for stability
    !>   - Delta_t_D_crit: Dimensionless critical time step (scaled by characteristic time)
    !>
    !>   Deferred procedures (must be implemented by subclasses):
    !>   - compute_stab_params(): Compute critical time step based on mesh, properties
    !>
    !>   Implemented procedures:
    !>   - compute_Delta_t_D_crit(): Compute dimensionless critical time step
    !>
    !>   Usage pattern:
    !>   1. Extend this class for specific physics (diffusion, advection, etc.)
    !>   2. Implement compute_stab_params() with stability criterion
    !>   3. Call compute_stab_params() each time mesh or properties change
    !>   4. Use Delta_t_crit to constrain time step selection
    type, public, abstract :: stab_params_c
        real(kind=8) :: Delta_t_crit     !< Critical time step for numerical stability
        real(kind=8) :: Delta_t_D_crit   !< Dimensionless critical time step = Δt_crit / t_c [-]
    contains
        !procedure(compute_stab_params), public, deferred :: compute_stab_params !< Compute stability parameters (deferred)
        procedure :: compute_Delta_t_D_crit                             !< Compute dimensionless critical time step
    end type
    
    !> \brief Abstract interface for computing stability parameters
    !> \details 
    !>   Defines the signature that all subclasses must implement.
    !>   
    !>   Inputs:
    !>   - props_obj: Physical properties (diffusivity, velocity, reaction rates)
    !>   - mesh: Spatial discretization (grid spacing, cell volumes)
    !>   - time_step: Current time step size [s]
    !>
    !>   Outputs:
    !>   - this%Delta_t_crit: Updated critical time step
    !>
    !>   Typical stability criteria:
    !>   - Diffusion (explicit): Δt ≤ Δx² / (2D) where D = diffusivity
    !>   - Advection (CFL): Δt ≤ Δx / |v| where v = velocity
    !>   - Reaction: Δt ≤ 1 / k where k = reaction rate
    !>   - Combined: Δt = min(Δt_diff, Δt_adv, Δt_react)
    abstract interface
        !> \brief Compute stability parameters (interface definition)
        !> \param[in,out] this      Stability parameters object to update
        !> \param[in]     props_obj Physical properties object
        !> \param[in]     mesh      Spatial discretization object
        !> \param[in]     time_step Current time step size [s]
        !> \details 
        !>   Each subclass implements specific stability criterion.
        !>   Must update this%Delta_t_crit based on provided inputs.
        subroutine compute_stab_params(this,props_obj,mesh,time_step)
            import stab_params_c    !< Import base class
            import props_c          !< Import properties class
            import spatial_discr_c  !< Import spatial discretization class
            import time_discr_c     !< Import time discretization class
            implicit none
            class(stab_params_c) :: this                        !< Stability parameters object (modified)
            class(props_c), intent(in) :: props_obj             !< Physical properties (read-only)
            class(spatial_discr_c), intent(in) :: mesh          !< Mesh/grid (read-only)
            real(kind=8), intent(in) :: time_step               !< Time step size [s]
        end subroutine
    end interface
    
    contains
    
        !> \brief Compute dimensionless critical time step
        !> \param[in,out] this Stability parameters object to update
        !> \param[in]     t_c  Characteristic time scale [s]
        !> \details 
        !>   Computes dimensionless critical time step:
        !>   \f[
        !>     \Delta t_D^{\text{crit}} = \frac{\Delta t^{\text{crit}}}{t_c}
        !>   \f]
        !>
        !>   Where:
        !>   - Δt_crit: Critical time step from stability criterion [s]
        !>   - t_c: Characteristic time scale of problem [s]
        !>
        !>   Characteristic time examples:
        !>   - Diffusion: t_c = L² / D (length² / diffusivity)
        !>   - Advection: t_c = L / v (length / velocity)
        !>   - Reaction: t_c = 1 / k (inverse reaction rate)
        !>
        !>   Dimensionless form useful for:
        !>   - Non-dimensional analysis
        !>   - Comparing different physical processes
        !>   - Scaling studies
        !>
        !>   Example:
        !>   \code{.f90}
        !>   real(kind=8) :: t_c
        !>   t_c = mesh%L**2 / props%diffusivity  ! Diffusive characteristic time
        !>   call stab_params%compute_Delta_t_D_crit(t_c)
        !>   print *, 'Dimensionless critical time step:', stab_params%Delta_t_D_crit
        !>   \endcode
        subroutine compute_Delta_t_D_crit(this,t_c)
            implicit none
            class(stab_params_c) :: this                        !< Stability parameters object
            real(kind=8), intent(in) :: t_c                     !< Characteristic time
            
            this%Delta_t_D_crit = this%Delta_t_crit / t_c       !< Compute dimensionless critical time step
        end subroutine
        
end module