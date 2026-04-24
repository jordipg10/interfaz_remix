!> \file transport_m.f90
!> \brief 1D steady-state transport equation module
!> \details This module provides the transport_1D_c class for simulating 1D steady-state transport
!> with advection, dispersion, and reaction terms. The module extends the diffusion class to include
!> advective transport, making it suitable for solute transport in flowing groundwater or surface water.
!>
!> The governing equation is the steady-state advection-dispersion equation (ADE):
!> \f[
!> \frac{\partial}{\partial x}\left(D \frac{\partial c}{\partial x}\right) - 
!> \frac{\partial}{\partial x}(q c) + r c_r = 0
!> \f]
!> where:
!> - c is solute concentration
!> - c_r is external solute concentration
!> - D is dispersion coefficient (mechanical dispersion + molecular diffusion) [L²/T]
!> - q is Darcy flux [L/T]
!> - r is sink/source term [1/T]
!>
!> This module includes:
!> - Matrix assembly for finite difference/volume discretization
!> - Initialization from input files
!> - Output writing routines
!> - Mass balance error analysis for quality control
!> - Support for heterogeneous transport properties (spatially variable D, v)
!>
!> Boundary conditions supported:
!> - Dirichlet (fixed concentration)
!> - Neumann (fixed dispersive flux)
!> - Robin (fixed mass flux)
!>
!> The class extends diffusion_1D_c, inheriting diffusion/dispersion functionality and adding
!> advection-specific methods and properties.

module transport_m
    use diffusion_m, only: diffusion_1D_c !> Import base diffusion class and PDE base class
    use transport_properties_heterog_m, only: tpt_props_heterog_1D_c !> Import heterogeneous transport properties class
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c !> Import 1D spatial discretization classes
    implicit none !> Enforce explicit variable declarations
    save !> Preserve module variables between procedure calls
    private !> Private module scope: internal details hidden from outside modules
    !> \class transport_1D_c
    !> \brief 1D steady-state transport equation class
    !> \details Encapsulates the 1D steady-state advection-dispersion-reaction equation (ADE).
    !> Extends the diffusion_1D_c class to add advective transport capabilities.
    !>
    !> This class represents transport of solutes in 1D domains with:
    !> - Advection: Movement with bulk fluid flow (velocity field v)
    !> - Dispersion: Spreading due to molecular diffusion + mechanical dispersion (coefficient D)
    !> - Reaction: Source/sink terms (R) including decay, production, etc.
    !>
    !> Key features:
    !> - Heterogeneous properties: D(x), q(x), r(x) can vary spatially
    !> - Multiple boundary condition types (Dirichlet, Neumann, Robin)
    !> - Mass balance checking for solution verification
    !> - Integration with reactive transport models
    !>
    !> The class provides methods for:
    !> - Setting transport properties (velocity, dispersion, sources)
    !> - Assembling coefficient matrices for numerical solution
    !> - Computing mass balance errors for quality assurance
    !> - Reading input files and writing output files
    !>
    !> Inherits from diffusion_1D_c:
    !> - PDE framework (matrices, boundary conditions)
    !> - Spatial discretization (mesh, cell sizes)
    !> - Concentration fields
    !> - Time discretization (for transient extensions)
    
    type, public, extends(diffusion_1D_c) :: transport_1D_c !> 1D transport equation class, extends diffusion_1D_c to add advection
        !> \var tpt_props_heterog
        !> \brief Heterogeneous transport properties object
        !> \details Stores spatially variable transport properties including:
        !> - Flow field: q(x) at cell interfaces and centers [L/T]
        !> - Dispersion coefficient: D(x) (molecular diffusion + mechanical dispersion) [L²/T]
        !> - Source/sink terms: r(x) (reactions, decay, production) [M/(L³·T)]
        !> - Porosity: φ(x) (affects retardation and storage) [-]
        !>
        !> This object manages the spatial heterogeneity of transport parameters,
        !> allowing realistic representation of layered or variable media.
        type(tpt_props_heterog_1D_c) :: tpt_props_heterog !> Heterogeneous transport properties object containing velocity, dispersion, sources
    contains
    !> \brief Set heterogeneous transport properties object
    !> \details Assigns the transport properties to this transport equation instance
        procedure :: set_tpt_props_heterog_obj
        
    !> \brief Set concentration flag for reaction term
    !> \details Determines which cells have active source/sink terms based on R(x)
        procedure :: set_conc_r_flag=>set_conc_r_flag_tpt
        
    !> \brief Compute transport matrix for PDE
    !> \details Assembles coefficient matrix including advection and dispersion terms
        procedure :: compute_trans_mat_PDE=>compute_trans_mat_tpt
        
    !> \brief Compute mass balance error for ADE with Dirichlet discharge
    !> \details Checks conservation of mass by comparing fluxes and accumulation
        procedure :: mass_balance_error_ADE_stat_Dirichlet_discharge
        
    !> \brief Initialize PDE for transport
    !> \details Reads input files, sets up mesh, properties, and boundary conditions
        procedure :: initialise_PDE=>initialise_transport_1D
        
    !> \brief Write PDE solution to output file
    !> \details Outputs concentration profiles and diagnostics to files
        procedure :: write_PDE=>write_transport_1D
    end type
    
!**************************************************************************************************
!> \name External Procedure Interfaces
!> \brief Interface declarations for transport equation procedures defined in separate files
!> \details These interfaces declare the signatures of subroutines and functions that are
!> implemented in separate source files. This allows the compiler to check argument types
!> and counts when these procedures are called, improving code safety.
!**************************************************************************************************
    interface
        !> \brief Compute transport matrix for transport_1D_c object
        !> \details Assembles the coefficient matrix for the discretized transport equation,
        !> including both advection and dispersion terms. The matrix incorporates:
        !> - Advective fluxes: upwind or central differencing
        !> - Dispersive fluxes: central differencing
        !> - Boundary conditions: Dirichlet, Neumann, or Robin
        !> - Source/sink terms
        !>
        !> The assembled matrix is used in the linear system Ax = b for steady-state solutions
        !> or in time-stepping schemes for transient solutions.
        !>
        !> \param[in,out] this Transport equation object (matrix components are populated)
        
        subroutine compute_trans_mat_tpt(this)
            import transport_1D_c !> Import transport class definition from parent module
            implicit none !> Enforce explicit variable declarations
            class(transport_1D_c) :: this !> Transport equation object with matrix components to compute
        end subroutine
        
        !> \brief Initialize transport PDE from input files
        !> \details Reads all necessary input data for the transport problem:
        !> - Spatial discretization (mesh file: root_mesh.dat)
        !> - Transport properties (properties file: root_tpt_props.dat)
        !> - Boundary conditions (BC file: root_BCs.dat)
        !> - Initial conditions (IC file: root_IC.dat)
        !> - Source/sink terms (source file: root_source.dat)
        !>
        !> After reading, the subroutine:
        !> - Allocates arrays for concentrations and fluxes
        !> - Sets up the spatial mesh
        !> - Initializes property fields
        !> - Applies boundary and initial conditions
        !>
        !> \param[in,out] this Transport equation object to initialize
        !> \param[in] root Root name for input files (e.g., "problem" for problem_mesh.dat, etc.)
        !> \param[in] mesh_type Integer code for mesh type (e.g., 1 for uniform, 2 for non-uniform)
        
        subroutine initialise_transport_1D(this,path,root,mesh_type)
            import transport_1D_c !> Import transport class definition
            implicit none !> Enforce explicit variable declarations
            class(transport_1D_c) :: this !> Transport equation object to populate with input data
            character(len=*), intent(in) :: path !> Directory path for input files
            character(len=*), intent(in) :: root !> Root filename prefix for all input files
            integer(kind=4), intent(in) :: mesh_type !> Integer code for mesh type (e.g., 1 for uniform, 2 for non-uniform)
        end subroutine
        
        !> \brief Write transport PDE solution to output files
        !> \details Outputs the computed transport solution to formatted text files:
        !> - Concentration profiles: root_conc.dat (c vs. x for each time)
        !> - Time series: root_timeseries.dat (c vs. t at specific locations)
        !> - Mass balance: root_massbalance.dat (conservation diagnostics)
        !> - Diagnostics: root_diagnostics.dat (convergence, iterations, errors)
        !>
        !> Output format is typically column-based ASCII for easy plotting with Python, MATLAB, etc.
        !>
        !> \param[in] this Transport equation object containing solution
        !> \param[in] root Root name for output files
        !> \param[in] Time_out Array of output time snapshots [T]
        !> \param[in] output Solution array: output(spatial_point, time_snapshot) [M/L³]
        
        subroutine write_transport_1D(this)
            import transport_1D_c !> Import transport class definition
            implicit none !> Enforce explicit variable declarations
            class(transport_1D_c), intent(in) :: this !> Transport equation object with solution data
            ! character(len=*), intent(in) :: root !> Root filename prefix for output files
            ! real(kind=8), intent(in) :: Time_out(:) !> Array of output times: t₁, t₂, ..., tₙ [T]
            ! real(kind=8), intent(in) :: output(:,:) !> Solution array: c(x,t) with dimensions (N_space, N_time) [M/L³]
        end subroutine
        
        !> \brief Compute mass balance error for ADE with Dirichlet discharge boundary
        !> \details Calculates the relative mass balance error to verify conservation of mass.
        !> For steady-state transport with Dirichlet boundaries and discharge at the outlet:
        !>
        !> Mass balance: Inflow + Sources = Outflow + Sinks
        !>
        !> \f[
        !> \text{Error} = \frac{|\text{Inflow} + \text{Sources} - \text{Outflow} - \text{Sinks}|}
        !>                     {|\text{Inflow} + \text{Sources}|}
        !> \f]
        !>
        !> Components:
        !> - Inflow: Advective + dispersive flux at inlet boundary
        !> - Outflow: Advective + dispersive flux at outlet (discharge) boundary
        !> - Sources: Integrated source term over domain
        !> - Sinks: Integrated sink term over domain
        !>
        !> A small error (< 1e-6) indicates good numerical accuracy and correct implementation.
        !> Large errors suggest discretization problems, incorrect BCs, or coding errors.
        !>
        !> \param[in] this Transport equation object with solution
        !> \return mass_bal_err Relative mass balance error (dimensionless, typically 0-1) [-]
        
        function mass_balance_error_ADE_stat_Dirichlet_discharge(this) result(mass_bal_err)
            import transport_1D_c !> Import transport class definition
            implicit none !> Enforce explicit variable declarations
            class(transport_1D_c), intent(in) :: this !> Transport equation object with solution to check
            real(kind=8) :: mass_bal_err !> Computed mass balance error (0 = perfect conservation) [-]
        end function
    end interface
    
!*****************************************************************************************************************************
!> \name Module Procedures
!> \brief Implementation of transport equation methods
!*****************************************************************************************************************************
    contains
    
        !> \brief Set heterogeneous transport properties object
        !> \details Assigns a transport properties object to the transport equation instance.
        !> This copies the properties (velocity, dispersion, sources, porosity) from the input
        !> object to the internal storage.
        !>
        !> The properties object contains:
        !> - Velocity field: v(x) [L/T]
        !> - Dispersion coefficient: D(x) [L²/T]
        !> - Source terms: R(x) [M/(L³·T)]
        !> - Porosity: φ(x) [-]
        !>
        !> This must be called before matrix assembly or solution.
        !>
        !> \param[in,out] this Transport equation object (properties are set)
        !> \param[in] tpt_props_heterog Transport properties object to copy
        
        subroutine set_tpt_props_heterog_obj(this,tpt_props_heterog)
            implicit none !> Enforce explicit variable declarations
            class(transport_1D_c) :: this !> Transport equation object (modified in place)
            class(tpt_props_heterog_1D_c), intent(in)  :: tpt_props_heterog !> Input transport properties to assign
            this%tpt_props_heterog=tpt_props_heterog !> Copy transport properties to object's internal storage
        end subroutine
        
        

        !> \brief Set concentration flag for reaction term in transport equation
        !> \details Determines which computational cells have active source/sink terms.
        !> Creates a flag array where:
        !> - conc_r_flag(i) = 1 if source_term(i) > 0 (active source/sink at cell i)
        !> - conc_r_flag(i) = 0 if source_term(i) ≤ 0 (no source/sink at cell i)
        !>
        !> This flag is used to apply recharge concentrations or source concentrations
        !> only where sources are active, avoiding unnecessary computations in cells
        !> without sources.
        !>
        !> Common applications:
        !> - Recharge zones (infiltration with solute input)
        !> - Injection wells (point sources)
        !> - Distributed sources (e.g., fertilizer application)
        !>
        !> \param[in,out] this Transport equation object (conc_r_flag array is allocated and set)
        
        subroutine set_conc_r_flag_tpt(this)
            implicit none !> Enforce explicit variable declarations
            class(transport_1D_c) :: this !> Transport equation object (flag array will be created)
            integer(kind=4) :: i !> Loop counter for iterating over spatial cells [-]
            allocate(this%conc_r_flag(this%spatial_discr%Num_targets)) !> Allocate flag array with one element per computational cell
            this%conc_r_flag=0 !> Initialize all flags to 0 (no source/sink by default) [-]
            do i=1,this%spatial_discr%Num_targets !> Loop over all computational cells
                if (this%tpt_props_heterog%source_term(i)>0) then !> Check if cell i has positive source term (active source)
                    this%conc_r_flag(i)=1 !> Set flag to 1 (active) for this cell [-]
                end if !> End source term check
            end do !> End loop over cells
        end subroutine !> End of set_conc_r_flag_tpt subroutine 

!**************************************************************************************************
!> \name Commented-Out Analytical Solutions (Legacy/Reference)
!> \brief Verification benchmark functions currently disabled
!> \details These commented functions provide analytical solutions for specific steady-state
!> transport problems. They are preserved for:
!> - Code verification (compare numerical vs. analytical solutions)
!> - Regression testing
!> - Educational reference
!>
!> Available analytical solutions (commented out, not compiled):
!> 1. Linear flow: q(x) = -x, D = const, c(0) = 1, c(L) = 0
!>    - Uses error function (erf)
!>    - Solution: c(x) = -erf(x/√(2D))/erf(L/√(2D)) + 1
!> 2. Quadratic flow: q(x) = -x², D = const, c(0) = 1, c(L) = 0
!>    - Uses incomplete gamma function
!>    - More complex analytical form
!>
!> To use: uncomment, verify special function availability (erf, gamma), recompile.
!**************************************************************************************************

        ! function anal_sol_tpt_1D_stat_flujo_lin(this,x) result(conc)
        !     !> \brief Analytical solution for 1D steady transport with linear flow field
        !     !> \details Problem specification:
        !     !>   - D = constant (uniform dispersion coefficient) [L²/T]
        !     !>   - q(x) = -x (linear Darcy velocity field) [L/T]
        !     !>   - c(0) = 1 (Dirichlet BC at inlet) [M/L³]
        !     !>   - c(L) = 0 (Dirichlet BC at outlet) [M/L³]
        !     !> 
        !     !> Analytical solution:
        !     !> \f[
        !     !>   c(x) = 1 - \frac{\text{erf}(x/\sqrt{2D})}{\text{erf}(L/\sqrt{2D})}
        !     !> \f]
        !     !>
        !     !> This can be used to verify numerical solutions for problems with
        !     !> varying velocity fields.
        !     !>
        !     !> \param[in] this Transport equation object (provides D and L)
        !     !> \param[in] x Spatial position where to evaluate solution [L]
        !     !> \return conc Analytical concentration at position x
        !         class(transport_1D_c), intent(in) :: this !> Transport object with properties
        !         real(kind=8), intent(in) :: x !> Evaluation position [L]
        !         real(kind=8) :: conc !> Computed concentration
                
        !         real(kind=8) :: D !> Dispersion coefficient [L²/T]
        !         real(kind=8), parameter :: pi=4d0*atan(1d0) !> π constant (for reference, not used here) [-]
                
        !         D=this%tpt_props_heterog%dispersion(1) !> Extract dispersion from properties [L²/T]
    
        !         if (this%dimensionless.eqv..true.) then !> Check if problem is dimensionless
        !             conc=-erf(x/sqrt(2d0))/erf(1d0/sqrt(2d0)) + 1 !> Dimensionless form (D=1, L=1) [-]
        !         else !> Dimensional problem
        !             conc=-erf(x/sqrt(2d0*D))/erf(this%spatial_discr%measure/sqrt(2d0*D)) + 1d0 !> Full dimensional form
        !         end if !> End dimensionless check
        !     end function !> End analytical solution for linear flow
            
       
        !     function der_anal_sol_tpt_1D_stat_flujo_lin(this,x) result(der_conc)
        !     !> \brief Derivative of analytical solution for 1D steady transport with linear flow
        !     !> \details Problem specification:
        !     !>   - D = constant (uniform dispersion) [L²/T]
        !     !>   - q(x) = -x (linear velocity) [L/T]
        !     !>   - c(0) = 1, c(L) = 0 (Dirichlet BCs)
        !     !>
        !     !> Analytical derivative:
        !     !> \f[
        !     !>   \frac{dc}{dx} = -\frac{\sqrt{2}}{\text{erf}(L/\sqrt{2D})\sqrt{\pi D}}
        !     !>                    \exp\left(-\frac{x^2}{2D}\right)
        !     !> \f]
        !     !>
        !     !> Useful for:
        !     !> - Verifying flux computations
        !     !> - Checking gradient accuracy
        !     !> - Boundary condition implementation
        !     !>
        !     !> \param[in] this Transport equation object (provides D and L)
        !     !> \param[in] x Spatial position where to evaluate derivative [L]
        !     !> \return der_conc Analytical concentration gradient dc/dx at position x
        !         class(transport_1D_c), intent(in) :: this !> Transport object with properties
        !         real(kind=8), intent(in) :: x !> Evaluation position [L]
        !         real(kind=8) :: der_conc !> Computed concentration gradient
                
        !         real(kind=8) :: D,L !> Dispersion [L²/T] and domain length [L]
        !         real(kind=8), parameter :: pi=4d0*atan(1d0) !> π = 3.14159... [-]
                
        !         L=this%spatial_discr%measure !> Extract domain length [L]
        !         D=this%tpt_props_heterog%dispersion(1) !> Extract dispersion coefficient [L²/T]
        !         der_conc=(-sqrt(2d0)/(erf(this%spatial_discr%measure/sqrt(2d0*D))*sqrt(pi*D)))*exp(-(x**2)/(2d0*D)) !> Compute analytical derivative using erf and exponential
        !     end function !> End analytical derivative for linear flow
    
            
        
        !     function der_anal_sol_tpt_1D_stat_flujo_cuad(this,x) result(der_conc)
        !     !> \brief Derivative of analytical solution for 1D steady transport with quadratic flow
        !     !> \details Problem specification:
        !     !>   - D = constant (uniform dispersion) [L²/T]
        !     !>   - q(x) = -x² (quadratic Darcy velocity field) [L/T]
        !     !>   - c(0) = 1, c(L) = 0 (Dirichlet BCs)
        !     !>
        !     !> Analytical derivative (more complex than linear case):
        !     !> \f[
        !     !>   \frac{dc}{dx} = -\frac{3^{2/3} \exp(-x^3/(3D))}
        !     !>                         {\Gamma(1/3) - \text{incompl\_gamma\_term}}
        !     !> \f]
        !     !>
        !     !> where incompl_gamma_term = 0.327336564991358 is a precomputed value
        !     !> of the incomplete gamma function for specific problem parameters.
        !     !>
        !     !> This more challenging case tests numerical schemes under strongly
        !     !> varying velocity conditions.
        !     !>
        !     !> \param[in] this Transport equation object (provides D and L)
        !     !> \param[in] x Spatial position where to evaluate derivative [L]
        !     !> \return der_conc Analytical concentration gradient dc/dx at position x
        !         class(transport_1D_c), intent(in) :: this !> Transport object with properties
        !         real(kind=8), intent(in) :: x !> Evaluation position [L]
        !         real(kind=8) :: der_conc !> Computed concentration gradient
                
        !         real(kind=8) :: D,L !> Dispersion [L²/T] and domain length [L]
        !         real(kind=8), parameter :: pi=4d0*atan(1d0),incompl_gamma_term=0.327336564991358 !> π [-] and precomputed incomplete gamma function value [-]
                
        !         D=this%tpt_props_heterog%dispersion(1) !> Extract dispersion coefficient [L²/T]
        !         L=this%spatial_discr%measure !> Extract domain length [L]
    
        !         der_conc=-exp(-(x**3)/(3d0*D))*(3d0**(2d0/3d0))/(gamma(1d0/3d0)-incompl_gamma_term) !> Compute analytical derivative using gamma and exponential functions
        !     end function !> End analytical derivative for quadratic flow
        
end module !> End of transport_m module