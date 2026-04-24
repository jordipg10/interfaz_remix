!****************************************************************************************************************************************************
!> \file compute_trans_mat_flow_transient.f90
!> \brief Computes transmission matrix for transient radial flow problems
!>
!> \details
!> This subroutine constructs the transmission (or transmissibility) matrix for transient groundwater
!> flow in radial coordinates. The transmission matrix represents the geometric and spatial discretization
!> contributions to the discrete flow equation.
!>
!> **Physical Context:**
!>
!> For radial flow (e.g., pumping well, injection well), the flow equation in cylindrical coordinates is:
!> \f[
!>    S_s \frac{\partial h}{\partial t} = \frac{1}{r^{d-1}} \frac{\partial}{\partial r}\left(r^{d-1} K \frac{\partial h}{\partial r}\right) + Q
!> \f]
!>
!> Where:
!> - \f$ S_s \f$ = Specific storage [1/m]
!> - \f$ h \f$ = Hydraulic head [m]
!> - \f$ r \f$ = Radial coordinate [m]
!> - \f$ d \f$ = Spatial dimension (d=2 for 2D radial, d=3 for 3D radial)
!> - \f$ K \f$ = Hydraulic conductivity [m/s]
!> - \f$ Q \f$ = Source/sink term [1/s]
!>
!> **Dimensionless Formulation:**
!>
!> Using characteristic scales:
!> - \f$ r_D = r / L_c \f$ (dimensionless radius)
!> - \f$ h_D = h / H_c \f$ (dimensionless head)
!> - \f$ t_D = t K / (S_s L_c^2) \f$ (dimensionless time)
!>
!> The dimensionless equation becomes:
!> \f[
!>    \frac{\partial h_D}{\partial t_D} = \frac{1}{r_D^{d-1}} \frac{\partial}{\partial r_D}\left(r_D^{d-1} \frac{\partial h_D}{\partial r_D}\right)
!> \f]
!>
!> **Spatial Discretization:**
!>
!> Using finite differences on a radial grid with cell centers at r_i:
!> \f[
!>    \frac{\partial}{\partial r_D}\left(r_D^{d-1} \frac{\partial h_D}{\partial r_D}\right)\Bigg|_{r_i} \approx 
!>    \frac{r_{i+1/2}^{d-1}(h_{i+1}-h_i)/h_i - r_{i-1/2}^{d-1}(h_i-h_{i-1})/h_{i-1}}{\Delta r_i}
!> \f]
!>
!> **Transmission Matrix Structure:**
!>
!> The transmission matrix T is tridiagonal with:
!> \f[
!>    T_{i,i-1} = \frac{2r_i^{d-1}}{h_i} \quad \text{(sub-diagonal)}
!> \f]
!> \f[
!>    T_{i,i+1} = \frac{2r_i^{d-1}}{h_i} \quad \text{(super-diagonal)}
!> \f]
!> \f[
!>    T_{i,i} = -(T_{i,i-1} + T_{i,i+1}) \quad \text{(diagonal)}
!> \f]
!>
!> This structure ensures conservation of mass at each grid point.
!>
!> \param[in,out] this Flow transient object containing spatial discretization and transmission matrix
!>
!> \note The transmission matrix is reallocated if adaptive refinement is active
!> \warning Only implemented for radial spatial discretization (spatial_discr_rad_c)
!> \see flow_transient_c, spatial_discr_rad_c, allocate_trans_mat
!****************************************************************************************************************************************************
subroutine compute_trans_mat_flow_transient(this)
    !use PDE_m, only: PDE_1D_c
    !use transport_m, only: transpoRT_c, diffusion_1D_c
    use flow_transient_m, only: flow_transient_c  !< Import transient flow class
    use spatial_discr_rad_m, only: spatial_discr_rad_c
    use spatial_discr_m, only: spatial_discr_c
                                                    !< Import radial discretization classes
    implicit none
    
!****************************************************************************************************************************************************
!> \subsection arguments Argument Declarations
    
    class(flow_transient_c) :: this !< Flow transient object containing:
                                    !< - spatial_discr: Radial mesh discretization
                                    !< - trans_mat: Transmission matrix (tridiagonal)
                                    !< - dimless: Flag for dimensionless formulation
                                    !< - char_params_flow: Characteristic flow parameters
                                    
!****************************************************************************************************************************************************
!> \subsection local_variables Local Variable Declarations
    
    real(kind=8) :: r_i     !< Dimensionless radial coordinate at cell interface [dimensionless]
                            !< r_i = r_{i+1/2} (interface between cells i and i+1)
                            
    real(kind=8) :: h_i     !< Dimensionless radial step size [dimensionless]
                            !< h_i = r_{D,i+1} - r_{D,i}
                            
    real(kind=8) :: r_min_D !< Dimensionless minimum radius (inner boundary) [dimensionless]
                            !< r_min_D = r_min / L_c (not currently used)
                            
    integer(kind=4) :: i    !< Loop counter for grid points
    integer(kind=4) :: n    !< Total number of grid points (targets)
    integer(kind=4) :: opcion !< Option variable (declared but not used)
    
!****************************************************************************************************************************************************
!> \subsection preprocessing Preprocessing: Initialize dimensions and handle mesh adaptation
    
    !> Extract total number of grid points from spatial discretization
    n=this%spatial_discr%Num_targets  !< n = total number of targets (grid points)

    !> \subsection adaptive_mesh Check for adaptive mesh refinement
    !> If adaptive refinement is active, reallocate transmission matrix with new dimensions
    if (this%spatial_discr%adapt_ref.eq.1) then  !< Check if adaptive refinement flag is set
        !> Deallocate existing transmission matrix arrays
        deallocate(this%trans_mat%sub,this%trans_mat%diag,this%trans_mat%super)
                                            !< Free memory for old sub-diagonal, diagonal, super-diagonal
        
        !> Reallocate transmission matrix with updated mesh size
        call this%allocate_trans_mat()      !< Allocate new arrays based on current n
    end if
    
!****************************************************************************************************************************************************
!> \subsection main_computation Main Computation: Assemble transmission matrix elements
    
    !> \subsection polymorphic_dispatch Polymorphic type selection for radial discretization
    !> Use polymorphic dispatch to access radial-specific mesh properties
    select type (mesh=>this%spatial_discr)  !< Associate mesh with spatial_discr for type-specific access
    type is (spatial_discr_rad_c)           !< Case: Radial spatial discretization
        
        !> \subsection dimensionless_check Check for dimensionless formulation
        if (this%dimless.eqv..true.) then   !< If using dimensionless variables (r_D, h_D, t_D)
            
            !> Legacy code: dimensionless minimum radius (commented out)
            !r_min_D=mesh%r_min/this%char_params_flow%char_length
            !< (Commented) Would compute r_min_D = r_min / L_c for reference
            
            !> \subsection compute_subdiagonal Compute sub-diagonal elements of transmission matrix
            !> Loop over interior interfaces (i = 1 to n-1)
            do i=1,n-1  !< Loop over all cell interfaces
                
                !> Compute dimensionless radial coordinate at interface i+1/2
                r_i=mesh%r_min_D+sum(mesh%Delta_r_D(1:i)) !> $r_{D,i}$ radial coordinate
                                            !< r_D,i = r_min_D + Σ(Δr_D,j) for j=1..i
                                            !< This gives r at interface between cells i and i+1
                                            
                !> Compute dimensionless radial step from cell center i to i+1
                h_i=(mesh%targets(i+1)%coord_D(1)-mesh%targets(i)%coord_D(1)) !> $h_{D,i}$ radial step
                                            !< h_D,i = r_{D,i+1} - r_{D,i}
                                            !< Radial distance between adjacent cell centers
                                            
                !> Compute sub-diagonal transmission coefficient
                this%trans_mat%sub(i)=2d0*(r_i**(mesh%dim-1))/h_i !> sub-diagonal element
                                            !< T_{i,i-1} = 2 * r_i^(d-1) / h_i
                                            !< Factor r_i^(d-1) accounts for radial geometry
                                            !< - d=2: cylindrical (2D radial), factor = r_i
                                            !< - d=3: spherical (3D radial), factor = r_i²
                                            !< Factor 2 comes from averaging at interface
            end do  !< End loop over interfaces
        end if      !< End dimensionless check
    end select      !< End type selection
    
    !> \subsection symmetry_super Super-diagonal equals sub-diagonal (symmetric matrix)
    !> For radial flow with uniform properties, transmission matrix is symmetric
    this%trans_mat%super=this%trans_mat%sub !< T_{i,i+1} = T_{i,i-1} for all i
                                            !< Symmetry ensures flux continuity at interfaces
    
    !> \subsection compute_diagonal Compute diagonal elements from mass conservation
    !> Legacy initialization (commented): 
    !this%trans_mat%diag=0d0 !> initialise diagonal elements
    
    !> Diagonal elements ensure row sum equals zero (mass conservation)
    this%trans_mat%diag(2:n-1)=-this%trans_mat%sub(1:n-2)-this%trans_mat%super(2:n-1)
                                            !< T_{i,i} = -(T_{i,i-1} + T_{i,i+1}) for i=2..n-1
                                            !< Interior points: outflow = inflow at steady state
                                            !< Boundaries (i=1, i=n) handled separately by BCs
    
end subroutine !< End compute_trans_mat_flow_transient