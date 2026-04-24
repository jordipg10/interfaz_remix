!> \file diff_stab_params_m.f90
!> \brief Diffusion stability parameters module
!> \details Defines stability parameters for 1D diffusion equations to ensure numerical stability
!> in finite difference or finite volume schemes. The module computes the diffusion stability parameter
!> (beta) which controls the Courant-Friedrichs-Lewy (CFL) condition for explicit time integration.
!>
!> The stability parameter for diffusion is:
!> \f[
!> \beta = \frac{D \cdot \Delta t}{\phi \cdot (\Delta x)^2}
!> \f]
!> where:
!> - D is the dispersion/diffusion coefficient [L²/T]
!> - Δt is the time step [T]
!> - φ is the porosity [-]
!> - Δx is the spatial grid spacing [L]
!>
!> For explicit schemes, stability requires β ≤ 0.5 in 1D.
!> This module extends the base stability parameters class to handle diffusion-specific stability analysis.

module diff_stab_params_m
    use stability_parameters_m, only: stab_params_c !> Import base stability parameters class
    use diff_props_heterog_m, only: diff_props_heterog_1D_c, diff_props_heterog_2D_c !> Import heterogeneous diffusion properties class
    use spatial_discr_rad_m, only: spatial_discr_rad_c !> Import radial spatial discretization class
    use spatial_discr_m, only: spatial_discr_c !> Import spatial discretization class
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c !> Import 1D mesh classes
    use spatial_discr_2D_m, only: mesh_2D_Euler_homog_c, mesh_2D_Euler_heterog_c !> Import 2D mesh classes
    use properties_m, only: props_c !> Import vector operations module
    implicit none !> Enforce explicit variable declarations
    save !> Preserve module variables between procedure calls
    private !> Private module scope: internal details hidden from outside modules
    !> \class stab_params_diff_c
    !> \brief 1D diffusion equation stability parameters subclass
    !> \details Extends the base stability parameters class (stab_params_c) to store and compute
    !> stability parameters specific to diffusion/dispersion processes. The key parameter is beta,
    !> which quantifies the ratio of diffusive transport over a time step to the spatial grid scale.
    !>
    !> For explicit time integration schemes (e.g., forward Euler), the stability criterion is:
    !> β ≤ 0.5 for 1D problems
    !>
    !> Violating this condition leads to unbounded numerical oscillations and divergence.
    !> The module computes the maximum β across all grid cells to ensure global stability.
    
    type, public, extends(stab_params_c) :: stab_params_diff_c !> 1D diffusion equation stability parameters subclass extending base stability class
        !> \var beta
        !> \brief Dispersion stability parameter (diffusion number)
        !> \details The dimensionless dispersion stability parameter:
        !> \f$ \beta = \frac{D \cdot \Delta t}{\phi \cdot (\Delta x)^2} \f$
        !>
        !> Physical interpretation:
        !> - Ratio of diffusive distance (√(D·Δt)) to grid spacing Δx
        !> - Quantifies how far solute diffuses in one time step relative to cell size
        !> - β → 0: Very small time step or large grid, highly stable
        !> - β → 0.5: Approaching stability limit for explicit schemes
        !> - β > 0.5: Unstable for explicit schemes (oscillations, divergence)
        !>
        !> For heterogeneous properties, this stores the maximum β across all cells.
        !> Dimensionless [-]
        real(kind=8) :: beta !> Dispersion stability parameter: β = D·Δt/(φ·Δx²), must satisfy β ≤ 0.5 for stability [-]
    contains
        !> \brief Compute diffusion stability parameters
        !> \details Calculates the stability parameter for the current time step and grid
        procedure :: compute_stab_params_diff_1D !> Compute diffusion stability parameters for 1D problems
        procedure :: compute_stab_params_diff_2D !> Compute diffusion stability parameters for 2D problems
    end type
    
    contains
    
        !> \brief Compute diffusion stability parameters
        !> \details Calculates the diffusion stability parameter β for a given time step and spatial grid.
        !> The subroutine:
        !> 1. Extracts diffusion properties (dispersion, porosity) from the properties object
        !> 2. Computes β = D·Δt/(φ·Δx²) for each computational cell
        !> 3. Finds the maximum β across all cells (most restrictive constraint)
        !> 4. Checks if β ≤ 0.5 (stability criterion for explicit schemes)
        !> 5. Issues a warning if stability condition is violated
        !>
        !> For heterogeneous properties, each cell may have different D, φ, and Δx, so β varies spatially.
        !> The maximum β determines overall stability.
        !>
        !> \param[in,out] this Diffusion stability parameters object (beta is computed and stored)
        !> \param[in] props_obj Properties object containing diffusion coefficients and porosity
        !> \param[in] mesh Spatial discretization object containing grid spacing information
        !> \param[in] time_step Time step size Δt [T]
        
        subroutine compute_stab_params_diff_1D(this,props_obj,mesh,time_step)
            implicit none !> Enforce explicit variable declarations
            class(stab_params_diff_c) :: this !> Diffusion stability parameters object (beta will be set)
            class(diff_props_heterog_1D_c), intent(in) :: props_obj !> Properties object (base class, will be cast to diff_props_heterog_1D_c)
            class(spatial_discr_c), intent(in) :: mesh !> Spatial discretization object with grid information
            !class(time_discr_c), intent(in) :: time_discr !> Commented out: time discretization object (not currently used)
            real(kind=8), intent(in) :: time_step !> Time step size Δt for stability analysis [T]
            
            integer(kind=4) :: i !> Loop counter for iterating over grid cells [-]
            real(kind=8) :: beta_max !> Maximum stability parameter across all cells (most restrictive) [-]
            real(kind=8) :: beta !> Temporary storage for stability parameter at current cell [-]
                        
            select type (mesh) !> Type-selective block: cast props_obj to specific diffusion properties type
            type is (mesh_1D_Euler_homog_c) !> Check if properties object is heterogeneous diffusion properties
                !if (props_obj%homog_flag .eqv. .true.) then !> Check if properties are actually homogeneous (stored as arrays but uniform values)
                    beta_max=props_obj%diff_cent(1)*time_step/(props_obj%porosity(1)*mesh%Delta_x**2) !> Compute β for first cell: β₁ = D₁·Δt/(φ₁·Δx₁²) [-]
                    do i=2,mesh%Num_targets-mesh%targets_flag !> Loop over remaining computational cells (2 to N_cells)
                        beta=props_obj%diff_cent(i)*time_step/(props_obj%porosity(i)*mesh%Delta_x**2) !> Compute β for cell i: βᵢ = Dᵢ·Δt/(φᵢ·Δxᵢ²) [-]
                        if (beta>beta_max) then !> Check if current cell has larger β than previous maximum
                            beta_max=beta !> Update maximum β to current cell's value (most restrictive condition) [-]
                        end if !> End maximum comparison block
                    end do !> End loop over cells
                    !this%beta=beta_max !> Store maximum β in stability parameters object (global stability constraint) [-]
                !end if !> End homogeneous properties block
            type is (mesh_1D_Euler_heterog_c) !> Check if mesh is radial spatial discretization
                ! beta_max=(props%disp_cent(1)*time_step/props%porosity(1))*(&
                !     1d0/mesh%get_cell_size(1)**2+1d0/mesh%get_cell_size(1)**2) !> Compute β for first cell: β₁ = D₁·Δt/(φ₁·Δx₁²) [-]
                beta_max=props_obj%diff_cent(1)*time_step/(&
                    props_obj%porosity(1)*mesh%Delta_x(1)**2) !> Compute β for first cell: β₁ = D₁·Δt/(φ₁·Δx₁²) [-]
                do i=2,mesh%Num_targets-mesh%targets_flag !> Loop over remaining computational cells (2 to N_cells)
                    beta=props_obj%diff_cent(i)*time_step/(&
                        props_obj%porosity(i)*mesh%Delta_x(i)**2) !> Compute β for cell i: βᵢ = Dᵢ·Δt/(φᵢ·Δxᵢ²) [-]
                    if (beta>beta_max) then !> Check if current cell has larger β than previous maximum
                        beta_max=beta !> Update maximum β to current cell's value (most restrictive condition) [-]
                    end if !> End maximum comparison block
                end do !> End loop over cells
                !this%beta=beta_max !> Store maximum β in stability parameters object (global stability constraint) [-]
            end select !> End type-selective block
            this%beta=beta_max !> Store maximum β in stability parameters object (global stability constraint) [-]
            if (this%beta>5d-1) then !> Check if stability parameter exceeds critical value (β > 0.5 means unstable for explicit schemes)
                print *, "Unstable transport", this%beta !> Warning message: print stability parameter value to alert user
                !error stop "Dispersion condition violated" !> Commented out: would halt execution if stability violated (currently only warns)
            end if !> End stability check block
        end subroutine !> End of compute_stab_params_diff subroutine

        subroutine compute_stab_params_diff_2D(this,props_obj,mesh,time_step)
            implicit none !> Enforce explicit variable declarations
            class(stab_params_diff_c) :: this !> Diffusion stability parameters object (beta will be set)
            class(diff_props_heterog_2D_c), intent(in) :: props_obj !> Properties object (base class, will be cast to diff_props_heterog_2D_c)
            class(spatial_discr_c), intent(in) :: mesh !> Spatial discretization object with grid information
            real(kind=8), intent(in) :: time_step !> Time step size Δt for stability analysis [T]
            
            integer(kind=4) :: i !> Loop counter for iterating over grid cells [-]
            real(kind=8) :: beta_max !> Maximum stability parameter across all cells (most restrictive) [-]
            real(kind=8) :: beta !> Temporary storage for stability parameter at current cell [-]
                        
            select type (mesh) !> Type-selective block: cast props_obj to specific diffusion properties type
            type is (mesh_2D_Euler_homog_c) !> Check if properties object is heterogeneous diffusion properties
                beta_max=props_obj%diff_cent(1)*time_step*mesh%sq_hypot/(&
                    props_obj%porosity(1)*mesh%get_max_cell_size()**2) !> Compute β for first cell: β₁ = D₁·Δt/(φ₁·Δx₁²) [-]
                do i=2,mesh%Num_targets-mesh%targets_flag !> Loop over remaining computational cells (2 to N_cells)
                    beta=props_obj%diff_cent(i)*time_step*mesh%sq_hypot/(&
                        props_obj%porosity(i)*mesh%get_max_cell_size()**2) !> Compute β for cell i: βᵢ = Dᵢ·Δt/(φᵢ·Δxᵢ²) [-]
                    if (beta>beta_max) then !> Check if current cell has larger β than previous maximum
                        beta_max=beta !> Update maximum β to current cell's value (most restrictive condition) [-]
                    end if !> End maximum comparison block
                end do !> End loop over cells
            type is (mesh_2D_Euler_heterog_c) !> Check if mesh is radial spatial discretization
                beta_max=props_obj%diff_cent(1)*time_step*mesh%sq_hypot(1)/(&
                    props_obj%porosity(1)*mesh%get_cell_size(1)**2) !> Compute β for first cell: β₁ = D₁·Δt/(φ₁·Δx₁²) [-]
                do i=2,mesh%Num_targets-mesh%targets_flag !> Loop over remaining computational cells (2 to N_cells)
                    beta=props_obj%diff_cent(i)*time_step*mesh%sq_hypot(i)/(&
                        props_obj%porosity(i)*mesh%get_cell_size(i)**2) !> Compute β for cell i: βᵢ = Dᵢ·Δt/(φᵢ·Δxᵢ²) [-]
                    if (beta>beta_max) then !> Check if current cell has larger β than previous maximum
                        beta_max=beta !> Update maximum β to current cell's value (most restrictive condition
                    end if !> End maximum comparison block
                end do !> End loop over cells
            end select !> End type-selective block
            this%beta=beta_max !> Store maximum β in stability parameters object (global stability constraint) [-]
            if (this%beta>5d-1) then !> Check if stability parameter exceeds critical value (β > 0.5 means unstable for explicit schemes)
                print *, "Unstable transport", this%beta !> Warning message: print stability parameter value to alert user
            end if !> End stability check block
        end subroutine !> End of compute_stab_params_diff_2D subroutine
        
end module !> End of diff_stab_params_m module