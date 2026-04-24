!> \file main_PDE.f90
!> \brief Main driver subroutine for solving 1D partial differential equations (PDEs)
!> \details This subroutine serves as the top-level driver for solving 1D PDEs including:
!> - Diffusion equations
!> - Advection-dispersion equations (transport)
!> - Flow equations
!> - Reactive transport equations
!>
!> The subroutine follows a three-phase workflow:
!> 1. **Pre-processing**: Initialize PDE object from input files
!> 2. **Processing**: Solve the PDE and write output
!> 3. **Post-processing**: (implicit in solve_write_PDE)
!>
!> The PDE class handles:
!> - Spatial discretization (1D finite difference/volume methods)
!> - Time discretization (explicit/implicit schemes)
!> - Boundary conditions (Dirichlet, Neumann, Robin)
!> - Solution output and visualization
!>
!> \param[in,out] this PDE_1D_c object containing all PDE data (mesh, properties, solution)
!> \param[in] root Root directory path for input/output files
!>
!> \note This is a high-level interface. The actual numerical methods are implemented
!> in the solve_write_PDE method of the PDE_1D_c class.
!>
!> \warning Ensure input files are properly formatted before calling this subroutine.
!>
!> \see PDE_m::PDE_1D_c for the PDE class definition
!> \see solve_write_PDE for the solution and output method

subroutine main_PDE(this,path,root,mesh_type)
    use PDE_m, only: PDE_c !< Import 1D PDE class definition
    implicit none !< Enforce explicit variable declarations

!> \subsection variables Variable Declarations
    !> \param[in,out] this
    !> \brief PDE object to be initialized and solved
    !> \details This polymorphic object can represent various 1D PDEs:
    !> - Diffusion (heat/mass transfer)
    !> - Transport (advection-dispersion)
    !> - Flow (groundwater/surface water)
    !> - Reactive transport (chemistry + transport)
    !>
    !> The object encapsulates:
    !> - Spatial mesh (node coordinates, cell sizes)
    !> - Material properties (diffusivity, velocity, porosity)
    !> - Boundary conditions (types and values)
    !> - Time discretization parameters
    !> - Solution arrays (concentration, velocity, etc.)
    class(PDE_c) :: this !< PDE object (polymorphic - can be diffusion, transport, etc.)
    character(len=*), intent(in) :: path !< Directory path for input files
    character(len=*), intent(in) :: root !< Root directory path for I/O operations
    integer(kind=4), intent(in) :: mesh_type !< Integer code for mesh type (e.g., 1 for uniform, 2 for non-uniform)
    !> \param[in] root
    !> \brief Root directory path for input and output files
    !> \details This path is used to locate:
    !> - Input files: mesh.dat, properties.dat, BCs.dat, etc.
    !> - Output files: solution.dat, mass_balance.dat, etc.
    !>
    !> Example: root = './simulations/case1/'
    !character(len=*), intent(in) :: root !< Root directory path for I/O operations
    
    !> \var N_t
    !> \brief Number of time steps [-]
    !> \details Total number of temporal discretization points for transient simulations.
    !> For steady-state problems, this may be 1.
    integer(kind=4) :: N_t !< Number of time steps
    
    !> \var Time_out
    !> \brief Output time points array [T]
    !> \details Array containing specific times at which solution is written to files.
    !> Allows selective output instead of writing at every time step.
    real(kind=8), allocatable :: Time_out(:) !< Output time points [T]
    
    !> \var Delta_t
    !> \brief Time step sizes array [T]
    !> \details Array of time increments for each step. Can be:
    !> - Constant: Delta_t(i) = constant for all i
    !> - Variable: Adaptive time stepping based on stability/accuracy criteria
    real(kind=8), allocatable :: Delta_t(:) !< Time step sizes [T]
    
    !> \var Delta_r
    !> \brief Radial grid spacing array [L]
    !> \details For radial/cylindrical coordinate systems. Spacing between nodes in r-direction.
    !> Used in cylindrical/spherical diffusion problems.
    real(kind=8), allocatable :: Delta_r(:) !< Radial grid spacing [L]
    
    !> \var Delta_r_D
    !> \brief Radial grid spacing for diffusion [L]
    !> \details Radial spacing specifically for diffusion coefficient evaluation.
    !> May differ from Delta_r if staggered grid is used.
    real(kind=8), allocatable :: Delta_r_D(:) !< Radial grid spacing for diffusion [L]
    
    !> \var Final_time
    !> \brief Final simulation time [T]
    !> \details Total duration of the simulation. For steady-state, represents pseudo-time
    !> until convergence.
    real(kind=8) :: Final_time !< Final simulation time [T]
    
    !> \var measure
    !> \brief Generic measurement variable
    !> \details Purpose depends on context - may store norms, errors, or diagnostic values.
    real(kind=8) :: measure !< Generic measurement variable
    
    !> \var Delta_x
    !> \brief Spatial grid spacing [L]
    !> \details Uniform spacing between nodes in Cartesian x-direction.
    !> For uniform grids: Delta_x = length / (N_cells - 1)
    real(kind=8) :: Delta_x !< Spatial grid spacing [L]
    
    !> \var length
    !> \brief Domain length [L]
    !> \details Total length of the 1D computational domain from x=0 to x=L.
    real(kind=8) :: length !< Domain length [L]
    
    !> \var L2_norm_vi
    !> \brief L2 norm of velocity interface values [L/T]
    !> \details Euclidean norm of velocity at cell interfaces. Used for:
    !> - Stability analysis (CFL condition)
    !> - Convergence checking
    !> - Solution diagnostics
    real(kind=8) :: L2_norm_vi !< L2 norm of interface velocities [L/T]
    
    !> \var radius
    !> \brief Radial coordinate or domain radius [L]
    !> \details For radial/cylindrical problems, represents:
    !> - Current radial position
    !> - Outer boundary radius
    !> - Characteristic length scale
    real(kind=8) :: radius !< Radius [L]
    
    !> \var a
    !> \brief Generic parameter (problem-dependent)
    !> \details Multipurpose variable that may represent:
    !> - Diffusion coefficient [L²/T]
    !> - Reaction rate constant [1/T]
    !> - Geometric parameter [-]
    real(kind=8) :: a !< Generic parameter
    
    !> \var Delta_r_0
    !> \brief Initial/reference radial spacing [L]
    !> \details Reference value for radial grid generation, used as base spacing
    !> before stretching or refinement.
    real(kind=8) :: Delta_r_0 !< Initial radial spacing [L]
    
    !> \var scheme
    !> \brief Numerical scheme identifier [-]
    !> \details Integer flag specifying time integration method:
    !> - 1: Explicit Euler (forward difference)
    !> - 2: Implicit Euler (backward difference)
    !> - 3: Crank-Nicolson (centered)
    !> - 4: Runge-Kutta variants
    integer(kind=4) :: scheme !< Numerical scheme flag
    
    !> \var int_method
    !> \brief Integration method identifier [-]
    !> \details Flag for spatial integration/discretization:
    !> - 1: Finite differences
    !> - 2: Finite volumes
    !> - 3: Finite elements
    integer(kind=4) :: int_method !< Integration method flag
    
    !> \var Num_time
    !> \brief Number of time steps (synonym for N_t) [-]
    integer(kind=4) :: Num_time !< Number of time steps
    
    !> \var BCs
    !> \brief Boundary condition type identifier [-]
    !> \details Specifies BC types at domain boundaries:
    !> - 1: Dirichlet (fixed value)
    !> - 2: Neumann (fixed flux)
    !> - 3: Robin (mixed)
    integer(kind=4) :: BCs !< Boundary condition type flag
    
    !> \var opcion_BCs
    !> \brief Boundary condition option/variant [-]
    !> \details Sub-option within BC type (e.g., time-dependent vs constant)
    integer(kind=4) :: opcion_BCs !< BC option flag
    
    !> \var parameters_flag
    !> \brief Parameters input flag [-]
    !> \details Indicates how parameters are specified:
    !> - 1: Read from file
    !> - 2: Hardcoded defaults
    !> - 3: Computed from other inputs
    integer(kind=4) :: parameters_flag !< Parameters specification flag
    
    !> \var i
    !> \brief Loop index for general iterations [-]
    integer(kind=4) :: i !< Loop index
    
    !> \var j
    !> \brief Secondary loop index for nested iterations [-]
    integer(kind=4) :: j !< Secondary loop index
    
    !> \var Num_targets
    !> \brief Number of target points for output [-]
    !> \details Number of spatial locations where solution is monitored/output
    integer(kind=4) :: Num_targets !< Number of target output points
    
    !> \var eqn_flag
    !> \brief Equation type identifier [-]
    !> \details Specifies which PDE is being solved:
    !> - 1: Diffusion
    !> - 2: Advection-diffusion
    !> - 3: Transport with reaction
    integer(kind=4) :: eqn_flag !< Equation type flag
    
    !> \var dim
    !> \brief Spatial dimension [-]
    !> \details Always 1 for this 1D solver (included for generality)
    integer(kind=4) :: dim !< Spatial dimension (=1)
    
    !> \var niter_max
    !> \brief Maximum number of iterations for iterative solvers [-]
    !> \details Upper limit for:
    !> - Nonlinear solvers (Newton, Picard)
    !> - Iterative linear solvers (GMRES, BiCGSTAB)
    integer(kind=4) :: niter_max !< Maximum iterations
    
    !> \var niter
    !> \brief Current iteration counter [-]
    !> \details Tracks iterations in iterative solution methods
    integer(kind=4) :: niter !< Current iteration number
    
    !> \var Dirichlet_BC_location
    !> \brief Location flag for Dirichlet boundary condition [-]
    !> \details Specifies where Dirichlet BC is applied:
    !> - 1: Left boundary (x=0)
    !> - 2: Right boundary (x=L)
    !> - 3: Both boundaries
    integer(kind=4) :: Dirichlet_BC_location !< Dirichlet BC location flag
    
    !> \var targets_flag
    !> \brief Target points specification flag [-]
    !> \details Indicates how target output points are defined:
    !> - 1: Read from file
    !> - 2: Uniform distribution
    !> - 3: All nodes
    integer(kind=4) :: targets_flag !< Target points flag
    
!****************************************************************************************************************************************************
!> \subsection preprocessing Pre-processing Phase
!> Initialize the PDE object from input files located in the root directory.
!> This reads:
!> - Spatial mesh definition
!> - Material properties (diffusivity, velocity, porosity)
!> - Boundary conditions (types and values)
!> - Time discretization parameters
!> - Initial conditions
!> - Output settings
    call this%initialise_PDE(path,root,mesh_type) !< Initialize PDE from input files in root directory
    
!****************************************************************************************************************************************************
!> \subsection processing Processing Phase
!> Solve the PDE and write output to files. This method:
!> - Assembles coefficient matrices
!> - Applies boundary conditions
!> - Advances solution in time (explicit or implicit)
!> - Writes solution at specified output times
!> - Computes mass balance and other diagnostics
!>
!> \param[out] Time_out Array of times at which output was written [T]
    call this%solve_write_PDE(Time_out) !< Solve PDE and write output at specified times
    
!****************************************************************************************************************************************************
!> End of main_PDE subroutine
end subroutine