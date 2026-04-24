!> \file transport_properties_heterog_m.f90
!> \brief Defines heterogeneous transport property classes for advection-dispersion problems
!> \details 
!>   This module extends diffusion properties to include advective transport, implementing
!>   heterogeneous (spatially variable) transport properties for 1D advection-dispersion-reaction.
!>
!>   Key Features:
!>   - Extends diff_props_heterog_1D_c to add flux/velocity fields
!>   - Supports spatially variable flux (Darcy velocity)
!>   - Computes dispersion from dispersivity and flux (D = α_L * |q|)
!>   - Handles constant and variable flux boundary conditions
!>   - Transient flux capabilities via time_fct_real_c
!>   - Linear and nonlinear flux profiles
!>   - Radial flow geometry support
!>   - Source/sink term computation from flux divergence
!>
!>   Transport Equation (1D):
!>   \f[
!>     \phi \frac{\partial c}{\partial t} = 
!>     \frac{\partial}{\partial x}\left(D \frac{\partial c}{\partial x}\right) - 
!>     \frac{\partial (qc)}{\partial x} + r
!>   \f]
!>   Where:
!>   - φ: Porosity
!>   - c: Concentration
!>   - D: Dispersion coefficient = α_L |q| (mechanical dispersion)
!>   - q: Darcy flux (volumetric flow rate per unit area)
!>   - α_L: Longitudinal dispersivity
!>   - r: Source/sink term
!>
!>   Flux Storage:
!>   - flux_int: Flux at cell interfaces (size: Num_cells + 1)
!>   - flux_cent: Flux at cell centers (size: Num_cells)
!>   - Staggered grid approach for numerical stability
!>
!>   File Format (_tpt_props.dat):
!>   ```
!>   source_flag [value]     ! Source term: logical [value if constant]
!>   porosity_flag [value]   ! Porosity: logical [value if constant]
!>   flux_flag value trans   ! Flux: logical value transient_flag
!>   alpha_L                 ! Longitudinal dispersivity
!>   ```
!>
!> \author jordi Petchamé-Guerrero
!> \date 2025
module transport_properties_heterog_m
    use spatial_discr_m, only: spatial_discr_c  !< Spatial discretization base class
    use diff_props_heterog_m, only: diff_props_heterog_1D_c, diff_props_heterog_2D_c  !< Diffusion properties base class
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c  !< 1D mesh classes
    use spatial_discr_rad_m, only: spatial_discr_rad_c  !< Radial mesh class
    use polynomials_m, only: real_poly_1D, der_real_poly_1D, sec_der_real_poly_1D  !< Polynomial evaluation
    use vectors_m, only: sum_squares  !< Vector utilities
    use time_fct_m, only: time_fct_real_c  !< Time-dependent function class
    implicit none  !< Require explicit variable declarations
    save           !< Preserve module variables between calls
    private      !< Private module scope by default
    public :: are_tpt_props_homog !< Public entities follow
    !> \brief Heterogeneous transport properties class for advection-dispersion
    !> \details 
    !>   Extends diff_props_heterog_1D_c to include advective flux and mechanical dispersion.
    !>   
    !>   Additional Member Variables:
    !>   - flux_int: Flux at cell interfaces [L/T] (Num_cells + 1)
    !>   - flux_cent: Flux at cell centers [L/T] (Num_cells)
    !>   - long_dispersivity: Longitudinal dispersivity α_L [L]
    !>   - cst_flux_flag: TRUE if flux is constant in space
    !>   - flux_trans: Time-dependent flux function
    !>
    !>   Inherited from diff_props_heterog_1D_c:
    !>   - porosity: Porosity φ [-]
    !>   - diffusion_int: Dispersion at interfaces [L²/T]
    !>   - diffusion_cent: Dispersion at centers [L²/T]
    !>   - source_term: Source/sink term [1/T]
    !>
    !>   Dispersion Relation:
    !>   \f[ D = \alpha_L |q| \f]
    !>   Mechanical dispersion scales with flux magnitude
    !>
    !>   Flux Computation:
    !>   - Linear: q(x) = q_0 + ∫r dx (mass balance)
    !>   - Nonlinear: q(x) = polynomial (user-specified)
    !>   - Radial: q(r) = Q / (2πrb) (cylindrical geometry)
    !>
    !>   Procedures:
    !>   - set_cst_flux_flag: Set constant flux indicator
    !>   - update_flux_int_trans: Update flux from transient series
    !>   - read_props: Read properties from file
    !>   - compute_dispersion: Calculate D from α_L and q
    !>   - compute_flux_lin: Linear flux from source term
    !>   - compute_flux_nonlin: Polynomial flux profile
    !>   - compute_flux_rad: Radial flux (cylindrical geometry)
    !>   - compute_source_term: Source from flux divergence
    !>   - are_props_homog: Check spatial homogeneity
    !>   - allocate_flux: Allocate flux arrays
    type, public, extends(diff_props_heterog_1D_c) :: tpt_props_heterog_1D_c
        real(kind=8), allocatable :: flux_int(:)     !< Flux at interfaces [L/T]
        real(kind=8), allocatable :: flux_cent(:)    !< Flux at cell centers [L/T]
        real(kind=8) :: long_dispersivity            !< Longitudinal dispersivity α_L [L]
        real(kind=8), allocatable :: disp_int(:)  !< hydrodynamic dispersion at interfaces [L²/T]
        real(kind=8), allocatable :: disp_cent(:) !< hydrodynamic dispersion at centers [L²/T]
        logical :: cst_flux_flag                     !< TRUE if flux constant in space
        type(time_fct_real_c) :: flux_trans          !< Transient flux (homogeneous)
    contains
        procedure :: set_cst_flux_flag           !< Set constant flux flag
        procedure :: update_flux_int_trans       !< Update flux from time series
        procedure :: read_props=>read_tpt_props_heterog_1D  !< Read from file
        procedure :: compute_dispersion_1D          !< Compute D = α_L |q|
        procedure :: compute_flux_lin            !< Linear flux profile
        procedure :: compute_flux_nonlin         !< Polynomial flux profile
        procedure :: compute_flux_rad            !< Radial flux
        procedure :: compute_source_term         !< Source from ∂q/∂x
        procedure :: allocate_tpt_props=>allocate_tpt_props_1D               !< Allocate flux arrays
    end type

    type, public, extends(tpt_props_heterog_1D_c) :: tpt_props_heterog_2D_c
        !> 2D heterogeneous transport properties
        real(kind=8) :: transv_dispersivity !< Transverse dispersivity α_T [L]
        real(kind=8) :: anisotropy_ratio   !< Anisotropy ratio α_T / α_L [-]
        real(kind=8), allocatable :: flux_int_y(:)  !< Flux in y-direction at interfaces [L/T]
        real(kind=8), allocatable :: flux_cent_y(:) !< Flux in y-direction at centers [L/T]
        real(kind=8), allocatable :: disp_tensor_int(:,:,:)  !< Hydrodynamic dispersion tensor at interfaces [L²/T]
        real(kind=8), allocatable :: disp_tensor_cent(:,:,:) !< Hydrodynamic dispersion tensor at centers [L²/T]
    contains
        !> Additional procedures for 2D transport properties can be defined here
        procedure :: allocate_tpt_props=>allocate_tpt_props_2D       !< Allocate 2D flux arrays
        procedure :: allocate_disp_tensor                      !< Allocate 2D flux arrays with separate Nx, Ny
        procedure :: compute_disp_tensor_2D       !< Compute 2D dispersion tensor
        procedure :: read_props=>read_tpt_props_heterog_2D  !< Read from file
    end type
    contains
        !> \brief Read source term from file (legacy/alternative format)
        !> \param[in,out] this     Transport properties object to populate
        !> \param[in]     filename Input file path
        !> \param[in]     mesh     Spatial discretization object
        !> \details 
        !>   Reads source/sink term from separate file.
        !>   
        !>   File format:
        !>   ```
        !>   cst_flag [value]  ! If TRUE, constant value follows
        !>   ```
        !>   or
        !>   ```
        !>   FALSE
        !>   r_1 r_2 ... r_N   ! Spatially variable source terms
        !>   ```
        !>
        !>   Validates array size against mesh targets.
        !>   Sets source_term_order=0 for constant source.
        !>
        !>   Note: This appears to be alternative to reading from _tpt_props.dat
        subroutine read_source_term_tpt(this,filename,mesh)
            implicit none
            class(tpt_props_heterog_1D_c) :: this              !< Transport properties object (modified)
            character(len=*), intent(in) :: filename        !< Input filename
            class(spatial_discr_c), intent(in) :: mesh      !< Spatial discretization object
            
            real(kind=8) :: r                               !< Constant source term value
            logical :: cst_source_term                      !< Flag: TRUE if constant source
            
            open(unit=1, file=filename, status='old', action='read')  !< Open input file
            read(1,*) cst_source_term                       !< Read constant source flag
            if (cst_source_term .eqv. .true.) then          !< Case: constant source term
                backspace(1)                                !< Rewind to re-read line with value
                read(1,*) cst_source_term, r                !< Read flag and constant value
                allocate(this%source_term(mesh%Num_targets - mesh%targets_flag))  !< Allocate array
                this%source_term = r                        !< Assign constant to all cells
                this%source_term_order = 0                  !< Order 0 = constant
            else if (allocated(this%flux_int) .and. allocated(this%flux_cent)) then  !< Flux already allocated
                continue                                    !< Skip (source computed from flux)
            else                                            !< Case: spatially variable source
                read(1,*) this%source_term                  !< Read source term array
                if (size(this%source_term) /= mesh%Num_targets - mesh%targets_flag) error stop "Dimension error in source term"  !< Validate size
            end if
            close(1)                                        !< Close file
        end subroutine

        !> \brief Set all transport properties programmatically (legacy/unused)
        !> \param[in,out] this            Transport properties object to modify
        !> \param[in]     porosity        Porosity at cell centers [-]
        !> \param[in]     diffusion_int  Dispersion at interfaces [L²/T]
        !> \param[in]     diffusion_cent Dispersion at centers [L²/T]
        !> \param[in]     flux_int        Flux at interfaces [L/T]
        !> \param[in]     flux_cent       Flux at centers [L/T]
        !> \details 
        !>   Direct setter for all transport properties (alternative to file I/O).
        !>   
        !>   Validates array dimensions:
        !>   - size(flux_cent) = size(porosity) = Num_cells
        !>   - size(flux_int) = size(porosity) + 1 = Num_cells + 1
        !>   - size(diffusion_cent) = size(porosity)
        !>   - size(diffusion_int) = size(flux_int)
        !>
        !>   Terminates program if dimension mismatch detected.
        !>
        !>   Note: Appears to be unused in current codebase.
        subroutine set_tpt_props_heterog(this,porosity,diffusion_int,diffusion_cent,flux_int,flux_cent)
            implicit none
            class(tpt_props_heterog_1D_c) :: this                  !< Transport properties object (modified)
            real(kind=8), intent(in) :: porosity(:)             !< Porosity array [-]
            real(kind=8), intent(in) :: diffusion_int(:)       !< Dispersion at interfaces [L²/T]
            real(kind=8), intent(in) :: diffusion_cent(:)      !< Dispersion at centers [L²/T]
            real(kind=8), intent(in) :: flux_int(:)             !< Flux at interfaces [L/T]
            real(kind=8), intent(in) :: flux_cent(:)            !< Flux at centers [L/T]
            
            if (size(flux_cent) /= size(porosity)) then         !< Validate flux_cent vs porosity
                error stop "Dimensions of porosity and flux at centres must be the same"
            else if (size(flux_int) /= size(flux_cent) + 1) then  !< Validate flux_int size
                error stop "Dimensions of fluxes are wrong"
            else if (size(porosity) /= size(diffusion_cent)) then  !< Validate diffusion_cent
                error stop "Dimensions of porosity and dispersion at cell centres must be the same"
            else if (size(diffusion_int) /= size(flux_int)) then   !< Validate diffusion_int
                error stop "Dimensions of dispersion at interfaces and flux at interfaces must be the same"
            else                                                !< All validations passed
                this%porosity = porosity                        !< Assign porosity
                this%diff_int = diffusion_int            !< Assign interface dispersion
                this%diff_cent = diffusion_cent          !< Assign center dispersion
                this%flux_int = flux_int                        !< Assign interface flux
                this%flux_cent = flux_cent                      !< Assign center flux
            end if
        end subroutine
        
        !> \brief Read transport properties from file (main I/O routine)
        !> \param[in,out] this           Transport properties object to populate
        !> \param[in]     root           Root filename (appends '_tpt_props.dat')
        !> \param[in]     spatial_discr  Spatial discretization object (optional)
        !> \details 
        !>   Primary file I/O routine for transport properties.
        !>   Reads from <root>_tpt_props.dat and <root>_flow_inf.dat (if transient).
        !>   
        !>   File format (_tpt_props.dat):
        !>   ```
        !>   source_flag [value]              ! Line 1: Source term
        !>   porosity_flag [value]            ! Line 2: Porosity
        !>   flux_flag [value] trans_flag     ! Line 3: Flux parameters
        !>   alpha_L                          ! Line 4: Longitudinal dispersivity
        !>   ```
        !>
        !>   Source Term (Line 1):
        !>   - TRUE value: Constant source = value, sets source_term_order=0
        !>   - FALSE: Heterogeneous (sets homog_flag=FALSE)
        !>
        !>   Porosity (Line 2):
        !>   - TRUE value: Constant porosity = value
        !>   - FALSE: Allocates array (values not read from this file)
        !>
        !>   Flux (Line 3):
        !>   - flux_flag=TRUE, value, trans_flag: Constant flux q=value
        !>   - flux_flag=FALSE, trans_flag: Non-constant spatial flux
        !>   - trans_flag=TRUE: Read transient series from _flow_inf.dat
        !>
        !>   Longitudinal Dispersivity (Line 4):
        !>   - α_L ≥ 0 (validation performed)
        !>   - Dispersion computed as D = α_L |q|
        !>
        !>   Initializes:
        !>   - homog_flag = TRUE (default, updated if heterogeneous)
        !>   - stat_flag = TRUE (default, updated if transient)
        !>
        !>   Example file:
        !>   ```
        !>   TRUE 0.0
        !>   TRUE 0.3
        !>   TRUE 1.0 FALSE
        !>   0.1
        !>   ```
        subroutine read_tpt_props_heterog_1D(this,root,spatial_discr)
            implicit none
            class(tpt_props_heterog_1D_c) :: this                  !< Transport properties object (modified)
            character(len=*), intent(in) :: root                !< Root filename
            class(spatial_discr_c), intent(in), optional :: spatial_discr  !< Spatial discretization
            
            integer(kind=4) :: n_flux, Num_cells               !< Flux array size, number of cells
            real(kind=8), parameter :: epsilon = 1d-12         !< Numerical tolerance
            real(kind=8) :: phi, D, q, r, alpha_L              !< Porosity, dispersion, flux, source, dispersivity
            logical :: flag, trans_flag                        !< General flag, transient flag
            character(len=100) :: label                         !< Input line buffer
            
            call this%set_homog_flag(.true.)                   !< Initialize as homogeneous (default)
            call this%set_stat_flag(.true.)                    !< Initialize as stationary (default)
            
            open(unit=1, file=root//'_tpt_props.dat', status='old', action='read')  !< Open transport properties file
            
            do 
                read(1,*) label  !< Read line
                if (trim(label) == 'TRANSPORT PROPERTIES') then  !< Skip empty lines
                    !> === Read source term (Line 1) ===
                    read(1,*) flag                                      !< Read source term flag
                    if (allocated(this%source_term)) deallocate(this%source_term)  !< Guard: pre-allocated by allocate_tpt_props
                    allocate(this%source_term(spatial_discr%Num_targets))
                    if (flag .eqv. .true.) then                         !< Case: constant source term
                        backspace(1)                                    !< Rewind to re-read with value
                        read(1,*) flag, r                               !< Read flag and constant value
                        this%source_term = r                            !< Assign constant to all cells
                        this%source_term_order = 0                      !< Order 0 = constant in space
                    else                                                !< Case: heterogeneous source term
                        call this%set_homog_flag(.false.)               !< Properties are heterogeneous
                    end if
                    
                    !> === Read porosity (Line 2) ===
                    read(1,*) flag                                      !< Read porosity flag
                    if (flag .eqv. .true.) then                         !< Case: constant porosity
                        backspace(1)                                    !< Rewind to re-read with value
                        read(1,*) flag, phi                             !< Read flag and constant value
                        if (allocated(this%porosity)) deallocate(this%porosity)  !< Guard
                        allocate(this%porosity(spatial_discr%Num_targets-spatial_discr%targets_flag))  !< Allocate porosity array
                        this%porosity = phi                             !< Assign constant to all cells
                    else                                                !< Case: heterogeneous porosity
                        if (.not. allocated(this%porosity)) &
                            allocate(this%porosity(spatial_discr%Num_targets-spatial_discr%targets_flag))
                    end if
                    
                    !> === Read flux (Line 3) ===
                    read(1,*) flag                                      !< Read constant flux flag
                    if ((flag .eqv. .true.) .and. (allocated(this%flux_int) .eqv. .true.)) then
                        !< Flux pre-allocated by allocate_tpt_props: just assign values
                        call this%set_cst_flux_flag(flag)
                        backspace(1)
                        read(1,*) flag, q, trans_flag
                        this%flux_int  = q
                        this%flux_cent = q
                    else
                        call this%set_cst_flux_flag(flag)
                        Num_cells = spatial_discr%Num_targets - spatial_discr%targets_flag
                        if (allocated(this%flux_cent)) deallocate(this%flux_cent)  !< Guard
                        if (allocated(this%flux_int))  deallocate(this%flux_int)   !< Guard
                        allocate(this%flux_cent(Num_cells), this%flux_int(Num_cells + 1))
                        
                        if (flag .eqv. .false.) then                    !< Case: non-constant spatial flux
                            call this%set_homog_flag(flag)
                            backspace(1)
                            read(1,*) flag, trans_flag                  !< Read flags (flux value not in file)
                        else                                            !< Case: constant spatial flux
                            backspace(1)                                !< Rewind line
                            read(1,*) flag, q, trans_flag               !< Read flag, flux value, transient flag
                            this%flux_int  = q                           !< Assign constant to interface flux
                            this%flux_cent = q                          !< Assign constant to center flux
                        end if
                        
                        if (trans_flag .eqv. .true.) then               !< Case: transient flux
                            call this%set_stat_flag(.false.)            !< Properties are time-dependent
                            open(unit=60, file=root//'_flow_inf.dat', status='old', action='read')  !< Open flow file
                            call this%flux_trans%read_time_series(60)   !< Read time series (assumes homogeneous)
                            close(60)                                   !< Close flow file
                        end if
                    end if
                    
                    !> === Read longitudinal dispersivity (Line 4) ===
                    read(1,*) alpha_L                                   !< Read dispersivity value
                    if (alpha_L < 0d0) then                             !< Validate non-negative
                        error stop "Error: longitudinal dispersivity must be non-negative"
                    else                                                !< Valid dispersivity
                        this%long_dispersivity = alpha_L                !< Store dispersivity
                    end if
                else if (trim(label) == 'end') then
                    exit                                              !< Skip empty lines
                else
                    continue
                    !error stop "Error reading transport properties file"
                end if
            end do
            close(1)                                            !< Close transport properties file
        end subroutine
        
        !> \brief Compute linear flux profile from source term
        !> \param[in,out] this             Transport properties object to update
        !> \param[in]     q_inf            Flux entering domain (boundary condition) [L/T]
        !> \param[in]     spatial_discr_obj Spatial discretization object
        !> \details 
        !>   Computes flux distribution assuming linear profile (steady-state mass balance).
        !>   
        !>   Mass Balance (1D):
        !>   \f[
        !>     \frac{dq}{dx} = r(x)
        !>   \f]
        !>   Integrated form:
        !>   \f[
        !>     q(x) = q_0 + \int_0^x r(\xi) d\xi
        !>   \f]
        !>
        !>   Discretization:
        !>   - flux_int(1) = q_inf (inflow boundary)
        !>   - flux_int(i+1) = flux_int(i) + r_i * Δx_i (accumulation)
        !>   - flux_cent(i) = (flux_int(i) + flux_int(i+1)) / 2 (average)
        !>
        !>   Allocates flux arrays (Num_cells + 1 interfaces, Num_cells centers).
        !>   Supports both uniform and non-uniform meshes via polymorphism.
        !>
        !>   Used for:
        !>   - Recharge/discharge scenarios
        !>   - Evapotranspiration
        !>   - Injection/extraction wells
        subroutine compute_flux_lin(this,q_inf,spatial_discr_obj)
            implicit none
            class(tpt_props_heterog_1D_c) :: this                  !< Transport properties object (modified)
            real(kind=8), intent(in) :: q_inf                   !< Inflow flux [L/T]
            class(spatial_discr_c), intent(in) :: spatial_discr_obj  !< Spatial discretization
            
            integer(kind=4) :: Num_cells, i                     !< Number of cells, loop index
            
            Num_cells = spatial_discr_obj%Num_targets - spatial_discr_obj%targets_flag  !< Compute cell count
            
            if (.not. allocated(this%flux_int)) then
                allocate(this%flux_int(Num_cells + 1))              !< Allocate flux arrays
            end if
            this%flux_int(1) = q_inf                            !< Set inflow boundary flux
            
            select type (spatial_discr_obj)                     !< Polymorphic dispatch on mesh type
            type is (mesh_1D_Euler_homog_c)                     !< Case: uniform mesh
                do i = 1, Num_cells                             !< Loop over cells
                    this%flux_int(i+1) = this%flux_int(i) + this%source_term(i) * spatial_discr_obj%Delta_x  !< Accumulate flux: q_{i+1} = q_i + r_i Δx
                    this%flux_cent(i) = 5d-1 * (this%flux_int(i) + this%flux_int(i+1))  !< Average to cell center: q_c = (q_i + q_{i+1})/2
                end do
            type is (mesh_1D_Euler_heterog_c)                   !< Case: non-uniform mesh
                do i = 1, Num_cells                             !< Loop over cells
                    this%flux_int(i+1) = this%flux_int(i) + this%source_term(i) * spatial_discr_obj%Delta_x(i)  !< Accumulate with variable Δx_i
                    this%flux_cent(i) = 5d-1 * (this%flux_int(i) + this%flux_int(i+1))  !< Average to cell center
                end do
            end select
        end subroutine
        
        !> \brief Compute nonlinear (polynomial) flux profile
        !> \param[in,out] this             Transport properties object to update
        !> \param[in]     flux_coeffs      Polynomial coefficients (decreasing degree order)
        !> \param[in]     spatial_discr_obj Spatial discretization object
        !> \details 
        !>   Evaluates user-specified polynomial flux at mesh points.
        !>   Assumes domain starts at x = x_0 (spatial_discr_obj%init_point(1)).
        !>   
        !>   Flux Polynomial:
        !>   \f[
        !>     q(x) = a_0 x^n + a_1 x^{n-1} + \cdots + a_n
        !>   \f]
        !>   where flux_coeffs = [a_0, a_1, ..., a_n] (length = deg + 1)
        !>
        !>   Evaluation Points:
        !>   - flux_int(1): At x = x_0 (inflow boundary)
        !>   - flux_int(i+1): At x = x_i (interface i)
        !>   - flux_cent(i): At x = x_{i-1/2} (cell center i)
        !>
        !>   Uniform Mesh:
        !>   - x_i = x_0 + i * Δx
        !>   - x_{i-1/2} = x_0 + (2i-1) * Δx / 2
        !>
        !>   Non-uniform Mesh:
        !>   - x_i = x_{i-1} + Δx_i (cumulative)
        !>   - x_{i-1/2} = x_{i-1} + Δx_i / 2
        !>
        !>   Used for analytical/verification problems with known flux profiles.
        subroutine compute_flux_nonlin(this,flux_coeffs,spatial_discr_obj)
            implicit none
            class(tpt_props_heterog_1D_c) :: this                  !< Transport properties object (modified)
            real(kind=8), intent(in) :: flux_coeffs(:)          !< Polynomial coefficients (decreasing order)
            class(spatial_discr_c), intent(in) :: spatial_discr_obj  !< Spatial discretization
            
            integer(kind=4) :: i, deg, Num_cells                !< Loop index, polynomial degree, cell count
            real(kind=8) :: x, dqn_dx, d2qn_dx2, Delta_x_n, x0  !< Position, derivatives, cell size, initial x
            
            deg = size(flux_coeffs) - 1                         !< Degree = # coefficients - 1
            Num_cells = spatial_discr_obj%Num_targets - spatial_discr_obj%targets_flag  !< Compute cell count
            
            allocate(this%flux_cent(Num_cells), this%flux_int(Num_cells + 1))  !< Allocate flux arrays
            x0 = spatial_discr_obj%init_point(1)                !< Get initial x-coordinate
            this%flux_int(1) = real_poly_1D(flux_coeffs, x0)    !< Evaluate polynomial at x_0
            
            select type (spatial_discr_obj)                     !< Polymorphic dispatch on mesh type
            type is (mesh_1D_Euler_homog_c)                     !< Case: uniform mesh
                do i = 1, Num_cells                             !< Loop over cells
                    this%flux_cent(i) = real_poly_1D(flux_coeffs, x0 + spatial_discr_obj%Delta_x * (2*i - 1) / 2d0)  !< Evaluate at cell center
                    this%flux_int(i+1) = real_poly_1D(flux_coeffs, x0 + spatial_discr_obj%Delta_x * i)  !< Evaluate at interface i+1
                end do
            type is (mesh_1D_Euler_heterog_c)                   !< Case: non-uniform mesh
                x = x0                                          !< Initialize position
                do i = 1, Num_cells                             !< Loop over cells
                    this%flux_cent(i) = real_poly_1D(flux_coeffs, x + spatial_discr_obj%Delta_x(i) / 2d0)  !< Evaluate at center
                    this%flux_int(i+1) = real_poly_1D(flux_coeffs, x + spatial_discr_obj%Delta_x(i))  !< Evaluate at interface
                    x = x + spatial_discr_obj%Delta_x(i)        !< Advance position
                end do
            end select
        end subroutine
        
        !> \brief Compute source term from flux polynomial divergence
        !> \param[in,out] this             Transport properties object to update
        !> \param[in]     spatial_discr_obj Spatial discretization object
        !> \param[in]     flux_coeffs      Flux polynomial coefficients (decreasing order)
        !> \details 
        !>   Computes source/sink term from flux divergence via mass balance.
        !>   
        !>   Mass Conservation:
        !>   \f[
        !>     r(x) = \frac{dq}{dx}
        !>   \f]
        !>
        !>   For polynomial flux q(x) = Σ a_i x^i:
        !>   \f[
        !>     \frac{dq}{dx} = \Sigma i \cdot a_i x^{i-1}
        !>   \f]
        !>
        !>   Implementation:
        !>   - Linear flux (deg=1): dq/dx = constant (a_0)
        !>   - Quadratic flux (deg=2): dq/dx evaluated at cell centers
        !>   - Higher degree: Not implemented (error)
        !>
        !>   Alternative (scheme 2):
        !>   Finite difference approximation:
        !>   \f[
        !>     r_i \approx \frac{q_{i+1} - q_i}{\Delta x_i}
        !>   \f]
        !>
        !>   Only for uniform mesh (mesh_1D_Euler_homog_c).
        subroutine compute_source_term(this,spatial_discr_obj,flux_coeffs)
            implicit none
            class(tpt_props_heterog_1D_c) :: this                  !< Transport properties object (modified)
            class(spatial_discr_c), intent(in) :: spatial_discr_obj  !< Spatial discretization
            real(kind=8), intent(in) :: flux_coeffs(:)          !< Flux polynomial coefficients (decreasing order)
            
            integer(kind=4) :: i                                !< Loop index
            real(kind=8) :: dq_dx                               !< Flux derivative dq/dx
            
            allocate(this%source_term(spatial_discr_obj%Num_targets - spatial_discr_obj%targets_flag))  !< Allocate source array
            
            select type (spatial_discr_obj)                     !< Polymorphic dispatch
            type is (mesh_1D_Euler_homog_c)                     !< Case: uniform mesh only
                if (spatial_discr_obj%scheme == 1 .and. spatial_discr_obj%targets_flag == 0) then  !< Traditional FD scheme
                    if (size(flux_coeffs) == 2) then            !< Case: linear flux q(x) = a*x + b
                        dq_dx = flux_coeffs(1)                  !< dq/dx = a (constant)
                        do i = 1, size(this%flux_cent)          !< Assign to all cells
                            this%source_term(i) = dq_dx         !< r_i = constant
                        end do
                    else if (size(flux_coeffs) == 3) then       !< Case: quadratic flux q(x) = a*x² + b*x + c
                        do i = 1, size(this%flux_cent)          !< Loop over cells
                            dq_dx = der_real_poly_1D(flux_coeffs, spatial_discr_obj%Delta_x * (2*i - 1) / 2d0)  !< Evaluate dq/dx at cell center
                            this%source_term(i) = dq_dx         !< r_i = dq/dx(x_{i-1/2})
                        end do
                    else                                        !< Case: cubic or higher
                        error stop "Subroutine 'compute_source_term' not implemented yet for cubic fluxes"
                    end if
                else if (spatial_discr_obj%scheme == 2 .and. spatial_discr_obj%targets_flag == 0) then  !< IFDS scheme
                    do i = 1, size(this%source_term)            !< Loop over cells
                        this%source_term(i) = (this%flux_int(i+1) - this%flux_int(i)) / spatial_discr_obj%Delta_x  !< Finite difference: (q_{i+1} - q_i) / Δx
                    end do
                else                                            !< Unsupported scheme
                    error stop "Subroutine 'compute_source_term' not implemented yet for this scheme"
                end if
            end select
        end subroutine
        
        !> \brief Check if transport properties are spatially homogeneous
        !> \param[in,out] this Transport properties object to check
        !> \details 
        !>   Determines if properties are uniform across spatial domain.
        !>   Inherits diffusion homogeneity check, then adds flux uniformity test.
        !>   
        !>   Homogeneity Criteria:
        !>   1. Diffusion properties homogeneous (porosity, diffusivity constant)
        !>   2. Flux spatially uniform: |q(x) - q_0| < ε for all x
        !>
        !>   Algorithm:
        !>   - Call parent are_diff_props_homog() to check diffusion
        !>   - If homogeneous so far, check flux variance
        !>   - Compare all flux_int values against first value
        !>   - Set homog_flag=.false. if any difference exceeds tolerance
        !>
        !>   Tolerance: ε = 10^-12 (floating-point precision threshold)
        !>
        !>   Used for:
        !>   - Numerical method selection (homogeneous allows simpler schemes)
        !>   - Performance optimization (avoid spatial loops if uniform)
        subroutine are_tpt_props_homog(this)
            implicit none
            class(tpt_props_heterog_1D_c) :: this                  !< Transport properties object (modified)
            
            integer(kind=4) :: i                                !< Loop index
            real(kind=8), parameter :: eps=1d-12                !< Tolerance for homogeneity test [-]
            
            !call this%diff_props_heterog_1D_c%are_props_homog()    !< Check diffusion properties first (sets homog_flag)
            if (this%homog_flag.eqv..true.) then                !< If diffusion is homogeneous
                do i=2,size(this%flux_int)                      !< Loop over all flux nodes (skip first)
                    if (abs(this%flux_int(1)-this%flux_int(i))>eps) then  !< Check if flux differs from reference
                        this%homog_flag=.false.                 !< Mark as heterogeneous
                        exit                                    !< Exit loop early (non-homogeneous detected)
                    end if
                end do
            end if
        end subroutine
        
        !subroutine allocate_flux_change(this,Num_chg)
        !implicit none
        !class(tpt_props_heterog_1D_c) :: this !> heterogeneous transport properties object
        !integer(kind=4), intent(in) :: Num_chg !> number of time steps where flux changes
        !if (allocated(this%flux_change)) then
        !    error stop "Flux change already allocated"
        !end if
        !this%num_flux_chg=Num_chg
        !allocate(this%flux_change(this%num_flux_chg))
        !this%flux_change=0 !> initialize to zero
        !end subroutine
        
        !subroutine read_flux_change(this,filename)
        !implicit none
        !class(tpt_props_heterog_1D_c) :: this !> heterogeneous transport properties object
        !character(len=*), intent(in) :: filename !> file name
        !
        !integer(kind=4) :: i,Num_chg
        !
        !open(unit=1,file=filename,status='old',action='read')
        !read(1,*) Num_chg
        !!if (allocated(this%flux_change)) then
        !!    error stop "Flux change already allocated"
        !!end if
        !call this%allocate_flux_change(Num_chg)
        !do i=1,this%num_flux_chg
        !    read(1,*) this%flux_change(i)
        !end do
        !close(1)
        !end subroutine
        
        !> \brief Update flux from transient time series
        !> \param[in,out] this Transport properties object to update
        !> \param[in]     k    Time step index
        !> \details 
        !>   Updates interface flux field from pre-loaded time series.
        !>   Used for time-dependent boundary conditions (flow variability).
        !>   
        !>   Time Series Structure:
        !>   - flux_trans%time_series(k): Flux value at time step k [L/T]
        !>   - Loaded from _flow_inf.dat file (see read_tpt_props_heterog)
        !>
        !>   Usage Pattern:
        !>   1. Load time series once during initialization
        !>   2. Call this routine at each time step to update flux
        !>   3. Flux remains constant until next call
        !>
        !>   Assumptions:
        !>   - Time series pre-allocated and populated
        !>   - Index k is valid: 1 ≤ k ≤ num_time_steps
        !>   - Flux assumed spatially homogeneous (single value applied to all interfaces)
        subroutine update_flux_int_trans(this,k)
        implicit none
        class(tpt_props_heterog_1D_c) :: this                      !< Transport properties object (modified)
        integer(kind=4), intent(in) :: k                        !< Time step index [-]

        this%flux_int=this%flux_trans%time_series(k)            !< Update flux array with time series value at step k
        end subroutine
        
        !> \brief Set constant flux indicator flag
        !> \param[in,out] this Transport properties object to modify
        !> \param[in]     flag Constant flux status (.true. = constant, .false. = variable)
        !> \details 
        !>   Controls whether flux is treated as constant in time.
        !>   Affects time integration strategy and performance optimization.
        !>   
        !>   Flag Meaning:
        !>   - .true.: Flux independent of time (steady-state flow)
        !>   - .false.: Flux varies with time (transient flow)
        !>
        !>   Impact on Simulation:
        !>   - Constant flux: Compute once, reuse throughout simulation
        !>   - Variable flux: Update from time series at each step (update_flux_int_trans)
        !>
        !>   Typically set during initialization based on:
        !>   - Input file specification (stat_flag line)
        !>   - Presence/absence of _flow_inf.dat file
        subroutine set_cst_flux_flag(this,flag)
        implicit none
        class(tpt_props_heterog_1D_c) :: this                      !< Transport properties object (modified)
        logical, intent(in) :: flag                             !< Constant flux flag: .true.=steady, .false.=transient
        this%cst_flux_flag=flag                                 !< Store flag value
        end subroutine
        
        

        !> \brief Compute radial flux distribution from volumetric flow rate
        !> \param[in,out] this Transport properties object to update
        !> \param[in]     mesh Radial mesh discretization object
        !> \param[in]     Q    Volumetric flow rate [L³/T]
        !> \details 
        !>   Computes radial Darcy flux distribution for cylindrical geometry.
        !>   Used in radial flow problems (wells, injection/extraction).
        !>   
        !>   Radial Darcy Flux:
        !>   \f[
        !>     q(r) = \frac{Q}{2\pi r b}
        !>   \f]
        !>   where:
        !>   - Q: Volumetric flow rate [L³/T]
        !>   - r: Radial coordinate [L]
        !>   - b: Aquifer thickness [L]
        !>
        !>   Evaluation Points:
        !>   - flux_int(1): At inner radius r_min
        !>   - flux_int(i+1): At r_{i} = r_{i-1} + Δr_i
        !>   - flux_cent(i): At cell center r_{i-1/2} = r_{i-1} + Δr_i/2
        !>
        !>   Assumptions:
        !>   - Axisymmetric flow (no θ or z variation)
        !>   - Steady-state radial flow
        !>   - Homogeneous aquifer properties
        subroutine compute_flux_rad(this,mesh,Q)
            implicit none
            class(tpt_props_heterog_1D_c) :: this                  !< Transport properties object (modified)
            class(spatial_discr_rad_c) :: mesh                  !< Radial mesh discretization
            real(kind=8), intent(in) :: Q                       !< Volumetric flow rate [L³/T]
            
            integer(kind=4) :: i                                !< Loop index
            real(kind=8) :: r_prev,r_i                          !< Previous and current radial coordinates [L]
            real(kind=8), parameter :: pi=4d0*atan(1d0)        !< π constant [-]

            r_i=mesh%r_min                                      !< Initialize at inner radius
            this%flux_int(1)=Q/(2d0*pi*mesh%r_min*mesh%targets(1)%thickness)  !< Flux at inner boundary: q(r_min)
            do i=1,mesh%Num_targets-mesh%targets_flag         !< Loop over radial cells
                this%flux_cent(i)=Q/(2d0*pi*(r_i+mesh%Delta_r(i)/2d0)*mesh%targets(i)%thickness)  !< Flux at cell center
                r_i=r_i+mesh%Delta_r(i)                         !< Update radial position: r_{i} = r_{i-1} + Δr_i
                this%flux_int(i+1)=Q/(2d0*pi*r_i*mesh%targets(i)%thickness)  !< Flux at outer interface: q(r_i)
            end do
        end subroutine

        !> \brief Compute mechanical dispersion coefficient from dispersivity and flux
        !> \param[in,out] this   Transport properties object to update
        !> \param[in]     scheme Spatial discretization scheme (1=traditional FD, 2=IFDS)
        !> \details 
        !>   Computes hydrodynamic dispersion from longitudinal dispersivity.
        !>   Mechanical dispersion is proportional to flow velocity (flux).
        !>   
        !>   Dispersion-Velocity Relationship:
        !>   \f[
        !>     D(x) = \alpha_L |q(x)|
        !>   \f]
        !>   where:
        !>   - D: Dispersion coefficient [L²/T]
        !>   - α_L: Longitudinal dispersivity [L]
        !>   - q: Darcy flux [L/T]
        !>
        !>   Evaluation Points (scheme-independent):
        !>   - diffusion_int: At cell interfaces (used with flux_int)
        !>   - diffusion_cent: At cell centers (used with flux_cent)
        !>
        !>   Physics:
        !>   - Mechanical dispersion arises from velocity heterogeneity at pore scale
        !>   - Larger flux → larger dispersion (mixing intensifies)
        !>   - Assumes linear relationship (valid for most groundwater flows)
        !>
        !>   Error Handling:
        !>   - Requires flux arrays allocated (compute_flux_* called first)
        subroutine compute_dispersion_1D(this,scheme)
            implicit none
            class(tpt_props_heterog_1D_c) :: this                  !< Transport properties object (modified)
            integer(kind=4), intent(in) :: scheme               !< Discretization scheme (1=FD, 2=IFDS) [-]
            
            integer(kind=4) :: i                                !< Loop index (unused, kept for compatibility)

            if (allocated(this%flux_int) .and. allocated(this%flux_cent)) then  !< Verify flux arrays exist
                this%disp_int=this%long_dispersivity*this%flux_int        !< D_int = α_L * q_int at interfaces
                this%disp_cent=this%long_dispersivity*this%flux_cent      !< D_cent = α_L * q_cent at centers
            else                                                !< Flux not computed yet
                error stop "Flux not allocated"
            end if
            ! if (scheme==1) then !> traditional FD scheme
            !     if (allocated(this%dispersivity) .and. allocated(this%flux_cent)) then
            !         allocate(this%dispersion(size(this%flux_cent)))
            !         this%dispersion=this%dispersivity*this%flux_cent
            !     else if (allocated(this%dispersivity)) then
            !         error stop "Flux not allocated"
            !     else
            !         error stop "Dispersivity not allocated"
            !     end if
            ! else if (scheme==2) then !> proposed scheme
            !     if (allocated(this%dispersivity) .and. allocated(this%flux_int)) then
            !         allocate(this%dispersion(size(this%flux_int)))
            !         this%dispersion=this%dispersivity*this%flux_int
            !     else if (allocated(this%dispersivity)) then
            !         error stop "Flux not allocated"
            !     else
            !         error stop "Dispersivity not allocated"
            !     end if
            ! else
            !     error stop "Subroutine 'compute_dispersion' not implemented yet for this scheme"
            ! end if
        end subroutine
        
        !> \brief Allocate memory for flux arrays
        !> \param[in,out] this      Transport properties object to modify
        !> \param[in]     Num_nodes Number of mesh nodes (cell interfaces)
        !> \details 
        !>   Allocates storage for flux fields on staggered grid.
        !>   Deallocates existing arrays to allow reallocation (e.g., mesh refinement).
        !>   
        !>   Array Dimensions:
        !>   - flux_int(Num_nodes): Flux at cell interfaces [L/T]
        !>     * Num_nodes = Num_cells + 1 (includes both boundaries)
        !>   - flux_cent(Num_nodes-1): Flux at cell centers [L/T]
        !>     * Num_nodes - 1 = Num_cells
        !>
        !>   Staggered Grid Structure:
        !>   \verbatim
        !>     Interface:  |----1----|----2----|--- ... ---|----N+1---|
        !>     Cell center:      1         2      ...           N
        !>                   ^         ^                      ^
        !>                flux_int   flux_int              flux_int
        !>                         flux_cent             flux_cent
        !>   \endverbatim
        !>
        !>   Safe Reallocation:
        !>   - Checks if arrays already allocated
        !>   - Deallocates before reallocating (prevents memory leak)
        !>   - Values not initialized (must be set by compute_flux_* routines)
        subroutine allocate_tpt_props_1D(this,Num_cells)
            class(tpt_props_heterog_1D_c) :: this                  !< Transport properties object (modified)
            integer(kind=4), intent(in) :: Num_cells            !< Number of mesh cells [-]
            if (allocated(this%flux_int)) then                  !< Check if interface flux already allocated
                deallocate(this%flux_int)                       !< Deallocate old array (allow resize)
            end if
            if (allocated(this%flux_cent)) then                 !< Check if center flux already allocated
                deallocate(this%flux_cent)                      !< Deallocate old array
            end if
            if (allocated(this%disp_cent)) then              !< Check if interface diffusion already allocated
                deallocate(this%disp_cent)                  !< Deallocate old array
            end if
            if (allocated(this%disp_int)) then               !< Check if center diffusion already allocated
                deallocate(this%disp_int)                   !< Deallocate old array
            end if
            call this%allocate_diff_props(Num_cells) !> Call base class allocation for common properties (source term, etc.)
            allocate(this%flux_int(Num_cells+1))                  !< Allocate interface flux (N+1 nodes)
            allocate(this%flux_cent(Num_cells))               !< Allocate center flux (N cells)
            allocate(this%disp_int(Num_cells+1))                  !< Allocate interface diffusion (N+1 nodes)
            allocate(this%disp_cent(Num_cells))               !< Allocate center diffusion (N cells
        end subroutine

        subroutine allocate_tpt_props_2D(this,Num_cells)
            class(tpt_props_heterog_2D_c) :: this
            integer(kind=4), intent(in) :: Num_cells
            call allocate_tpt_props_1D(this, Num_cells)
            if (allocated(this%flux_int_y))  deallocate(this%flux_int_y)
            if (allocated(this%flux_cent_y)) deallocate(this%flux_cent_y)
            allocate(this%flux_int_y(Num_cells+1))
            allocate(this%flux_cent_y(Num_cells))
            this%flux_int_y  = 0d0
            this%flux_cent_y = 0d0
        end subroutine

        !> \brief Allocate dispersion tensor arrays with correct rectangular shapes
        !> \details
        !>   disp_tensor_int  (Nx+1, Ny+1, 3)  — x/y-interface pairs
        !>           disp_tensor_cent (Nx,   Ny,   3)   — cell centres
        subroutine allocate_disp_tensor(this, Num_cells_x, Num_cells_y)
            class(tpt_props_heterog_2D_c) :: this
            integer(kind=4), intent(in) :: Num_cells_x  !< Number of cells in x-direction [-]
            integer(kind=4), intent(in) :: Num_cells_y  !< Number of cells in y-direction [-]

            if (allocated(this%disp_tensor_int))  deallocate(this%disp_tensor_int)
            if (allocated(this%disp_tensor_cent)) deallocate(this%disp_tensor_cent)

            allocate(this%disp_tensor_int (Num_cells_x + 1, Num_cells_y + 1, 3))  !< (Nx+1) x (Ny+1) x 3
            allocate(this%disp_tensor_cent(Num_cells_x,     Num_cells_y,     3))  !< Nx x Ny x 3

            this%disp_tensor_int  = 0d0
            this%disp_tensor_cent = 0d0
        end subroutine

        !> \brief Compute 2D mechanical dispersion tensor from dispersivities and flux
        !> \details  D_xx = alpha_L*qx^2/|q| + alpha_T*qy^2/|q|
        !>           D_yy = alpha_L*qy^2/|q| + alpha_T*qx^2/|q|
        !>           D_xy = (alpha_L - alpha_T)*qx*qy/|q|
        !> Stored in disp_tensor_int(i,j,1:3) = [D_xx, D_yy, D_xy] at each interface node (i,j)
        !> and disp_tensor_cent(i,j,1:3) at each cell centre (i,j).
        subroutine compute_disp_tensor_2D(this)
            implicit none
            class(tpt_props_heterog_2D_c) :: this

            integer(kind=4) :: i, j, Nx, Ny
            real(kind=8) :: qx, qy, q_mag

            if (.not. allocated(this%disp_tensor_int)) &
                error stop "disp_tensor_int not allocated: call allocate_disp_tensor first"
            if (.not. allocated(this%disp_tensor_cent)) &
                error stop "disp_tensor_cent not allocated: call allocate_disp_tensor first"

            !> --- interface nodes ---
            Nx = size(this%disp_tensor_int, 1) - 1  !< Num_cells_x
            Ny = size(this%disp_tensor_int, 2) - 1  !< Num_cells_y
            do j = 1, Ny + 1
                do i = 1, Nx + 1
                    qx    = this%flux_int(min(i, size(this%flux_int)))
                    qy    = this%flux_int_y(min(j, size(this%flux_int_y)))
                    q_mag = sqrt(qx**2 + qy**2)
                    if (q_mag > 0d0) then
                        this%disp_tensor_int(i,j,1) = this%long_dispersivity   * qx**2 / q_mag &
                                                    + this%transv_dispersivity * qy**2 / q_mag
                        this%disp_tensor_int(i,j,2) = this%long_dispersivity   * qy**2 / q_mag &
                                                    + this%transv_dispersivity * qx**2 / q_mag
                        this%disp_tensor_int(i,j,3) = (this%long_dispersivity - this%transv_dispersivity) &
                                                    * qx * qy / q_mag
                    else
                        this%disp_tensor_int(i,j,:) = 0d0
                    end if
                end do
            end do

            !> --- cell centres ---
            do j = 1, Ny
                do i = 1, Nx
                    qx    = this%flux_cent((j-1)*Nx + i)
                    qy    = this%flux_cent_y((j-1)*Nx + i)
                    q_mag = sqrt(qx**2 + qy**2)
                    if (q_mag > 0d0) then
                        this%disp_tensor_cent(i,j,1) = this%long_dispersivity   * qx**2 / q_mag &
                                                     + this%transv_dispersivity * qy**2 / q_mag
                        this%disp_tensor_cent(i,j,2) = this%long_dispersivity   * qy**2 / q_mag &
                                                     + this%transv_dispersivity * qx**2 / q_mag
                        this%disp_tensor_cent(i,j,3) = (this%long_dispersivity - this%transv_dispersivity) &
                                                     * qx * qy / q_mag
                    else
                        this%disp_tensor_cent(i,j,:) = 0d0
                    end if
                end do
            end do
        end subroutine

        !> \brief Read 2D transport properties from file
        subroutine read_tpt_props_heterog_2D(this, root, spatial_discr)
            implicit none
            class(tpt_props_heterog_2D_c) :: this
            character(len=*), intent(in) :: root
            class(spatial_discr_c), intent(in), optional :: spatial_discr

            integer(kind=4) :: Num_cells
            real(kind=8) :: phi, qx, qy, r, alpha_L, alpha_T
            logical :: flag, trans_flag
            character(len=100) :: label

            call this%set_homog_flag(.true.)
            call this%set_stat_flag(.true.)

            open(unit=1, file=root//'_tpt_props.dat', status='old', action='read')
            do
                read(1,*) label
                if (trim(label) == 'TRANSPORT PROPERTIES') then
                    !> === Source term ===
                    read(1,*) flag
                    if (allocated(this%source_term)) deallocate(this%source_term)  !< Guard
                    allocate(this%source_term(spatial_discr%Num_targets))
                    if (flag .eqv. .true.) then
                        backspace(1)
                        read(1,*) flag, r
                        this%source_term = r
                        this%source_term_order = 0
                    else
                        call this%set_homog_flag(.false.)
                    end if
                    !> === Porosity ===
                    read(1,*) flag
                    if (flag .eqv. .true.) then
                        backspace(1)
                        read(1,*) flag, phi
                        if (allocated(this%porosity)) deallocate(this%porosity)  !< Guard
                        allocate(this%porosity(spatial_discr%Num_targets-spatial_discr%targets_flag))
                        this%porosity = phi
                    else
                        if (.not. allocated(this%porosity)) &
                            allocate(this%porosity(spatial_discr%Num_targets-spatial_discr%targets_flag))
                    end if
                    !> === Flux ===
                    Num_cells = spatial_discr%Num_targets - spatial_discr%targets_flag
                    read(1,*) flag
                    call this%set_cst_flux_flag(flag)
                    if (flag .eqv. .false.) then
                        call this%set_homog_flag(flag)
                        backspace(1)
                        read(1,*) flag, trans_flag
                    else
                        backspace(1)
                        read(1,*) flag, qx, qy, trans_flag
                        !> flux arrays already allocated by allocate_tpt_props; just fill values
                        if (.not. allocated(this%flux_int)) then
                            if (allocated(this%flux_int))  deallocate(this%flux_int)
                            allocate(this%flux_int(Num_cells+1))
                        end if
                        if (.not. allocated(this%flux_cent)) then
                            if (allocated(this%flux_cent)) deallocate(this%flux_cent)
                            allocate(this%flux_cent(Num_cells))
                        end if
                        if (.not. allocated(this%flux_int_y)) then
                            if (allocated(this%flux_int_y))  deallocate(this%flux_int_y)
                            allocate(this%flux_int_y(Num_cells+1))
                        end if
                        if (.not. allocated(this%flux_cent_y)) then
                            if (allocated(this%flux_cent_y)) deallocate(this%flux_cent_y)
                            allocate(this%flux_cent_y(Num_cells))
                        end if
                        this%flux_int    = qx
                        this%flux_int_y  = qy
                        this%flux_cent   = qx
                        this%flux_cent_y = qy
                    end if
                    if (trans_flag .eqv. .true.) then
                        call this%set_stat_flag(.false.)
                        open(unit=60, file=root//'_flow_inf.dat', status='old', action='read')
                        call this%flux_trans%read_time_series(60)
                        close(60)
                    end if
                    !> === Dispersivities (alpha_L, alpha_T) ===
                    read(1,*) alpha_L, alpha_T
                    if (alpha_L < 0d0 .or. alpha_T < 0d0) &
                        error stop "Error: dispersivities must be non-negative"
                    this%long_dispersivity   = alpha_L
                    this%transv_dispersivity = alpha_T
                else if (trim(label) == 'end') then
                    exit
                else
                    continue
                end if
            end do
            close(1)
        end subroutine

end module