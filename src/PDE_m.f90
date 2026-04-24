!> \file PDE_m.f90
!> \brief One-dimensional partial differential equation (PDE) module
!> \details
!>   Abstract base class for solving 1D advection-dispersion-reaction PDEs.
!>   
!>   Governing PDE form:
!>   \f[
!>   \mathbf{T} c + g = 0
!>   \f]
!>   
!>   Where:
!>   - \f$\mathbf{T}\f$: Transition (transport) matrix (tridiagonal)
!>   - \f$c\f$: Concentration vector at discretization points
!>   - \f$g\f$: Source/sink term vector
!>   
!>   Source term composition:
!>   \f[
!>   g = \mathbf{R} c_{rech} + \mathbf{B} c_{BD}
!>   \f]
!>   - \f$\mathbf{R}\f$: Recharge matrix
!>   - \f$c_{rech}\f$: Recharge concentration
!>   - \f$\mathbf{B}\f$: Boundary matrix
!>   - \f$c_{BD}\f$: Boundary condition concentrations
!>   
!>   Applications:
!>   - Groundwater contaminant transport
!>   - Heat transfer
!>   - Reactive transport (base for coupled chemistry)
!>   - Diffusion processes
!>   
!>   Spatial discretization options:
!>   - Finite differences (Euler, Crank-Nicolson)
!>   - Cartesian or radial coordinates
!>   - Homogeneous or heterogeneous properties
!>   
!>   Solution methods:
!>   1. Direct numerical solution (linear system solver)
!>   2. Eigenvalue decomposition (analytical)
!>
!> \author jordi Petchamé-Guerrero
!> \date October 2025

module PDE_m
    use BCs_m, only: BCs_1D_c, BCs_2D_c                                                 !< Boundary conditions type
    use arrays_m, only: tridiag_matrix_c, diag_matrix_c, real_array_c, pentadiag_matrix_c  !< Linear system solvers and matrix types
    !use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c  !< 1D spatial discretization
    use spatial_discr_rad_m, only: spatial_discr_rad_c                        !< Radial spatial discretization
    use spatial_discr_m, only: spatial_discr_c                        !< Radial spatial discretization
    implicit none
    save
    private
    !> \brief Abstract 1D PDE superclass
    !> \details
    !>   Base class for one-dimensional partial differential equations.
    !>   Provides common infrastructure for transport, diffusion, and reaction PDEs.
    !>   
    !>   Matrix structure:
    !>   - Transition matrix T: tridiagonal (advection + dispersion)
    !>   - Recharge matrix R: diagonal (source/sink terms)
    !>   - External matrix E: general (additional terms)
    !>   - Boundary matrix B: vector (size 2 for 1D)
    !>   
    !>   Discretization:
    !>   - Polymorphic spatial_discr allows various mesh types
    !>   - Supports Cartesian and radial coordinates
    !>   - Handles homogeneous/heterogeneous material properties
    !>   
    !>   Solution approaches:
    !>   - Method 1: Numerical (direct linear system solve)
    !>   - Method 2: Eigendecomposition (analytical/semi-analytical)
    !>   
    !>   Dimensionless option:
    !>   - Flag to work with dimensionless variables
    !>   - Useful for stability analysis and benchmarking
    !>   
    !>   Deferred procedures (must be implemented by subclasses):
    !>   - initialise_PDE: Setup problem
    !>   - compute_trans_mat_PDE: Build transport matrix
    !>   - compute_source_term_PDE: Compute source/sink vector
    !>   - compute_rech_mat_PDE: Build recharge matrix
    !>   - solve_PDE_1D: Solve the PDE
    !>   - write_PDE: Output results
    type, public, abstract :: PDE_c
        class(spatial_discr_c), allocatable :: spatial_discr                      !< Spatial discretization (polymorphic)
        !type(diag_matrix_c) :: F_mat                                          !< Storage matrix (diagonal)
        type(diag_matrix_c) :: rech_mat                                       !< Recharge matrix R (diagonal)
        type(real_array_c) :: ext_mat                                         !< External matrix E
        real(kind=8), allocatable :: bd_mat(:)                                !< Boundary matrix B
        real(kind=8), allocatable :: source_term_PDE(:)                       !< Source term vector g
        logical :: dimless                                                    !< Dimensionless PDE flag (TRUE if dimensionless, FALSE otherwise)
        integer(kind=4) :: sol_method                                         !< Solution method: 1=Numerical, 2=Eigendecomposition
    contains
    !> Set procedures
        procedure :: set_spatial_discr                                !< Set spatial discretization
        !procedure(set_BCs), deferred :: set_BCs                                          !< Set boundary conditions
        procedure :: set_sol_method                                   !< Set solution method
        procedure :: set_dimless_flag                                 !< Set dimensionless flag
    !> Allocate procedures
        procedure(allocate_trans_mat), deferred :: allocate_trans_mat                               !< Allocate transition matrix
        procedure :: allocate_rech_mat                                !< Allocate recharge matrix
        procedure(allocate_bd_mat), deferred :: allocate_bd_mat                                  !< Allocate boundary matrix
        procedure :: allocate_source_term_PDE                         !< Allocate source term vector
        procedure(allocate_arrays_PDE), deferred :: allocate_arrays_PDE  !< Allocate all arrays
    !> Update procedures
        ! procedure :: update_trans_mat                                 !< Update transition matrix
    !> Initialise (deferred)
        procedure(initialise_PDE), deferred :: initialise_PDE         !< Initialize PDE (subclass-specific)
    !> Computation procedures (deferred)
        procedure(compute_trans_mat_PDE), deferred :: compute_trans_mat_PDE  !< Compute transport matrix
        procedure(write_PDE), deferred :: write_PDE             !< Write results to file
        procedure(compute_source_term_PDE), deferred :: compute_source_term_PDE  !< Compute source term
        procedure(compute_rech_mat_PDE), deferred :: compute_rech_mat_PDE  !< Compute recharge matrix
        procedure(solve_PDE), deferred :: solve_PDE             !< Solve PDE
        !procedure(compute_F_mat_PDE), public, deferred :: compute_F_mat_PDE   !< Compute storage matrix (deferred)
    !> Concrete procedures
        procedure :: solve_write_PDE                               !< Solve and write results
        procedure :: main_PDE                                         !< Main PDE driver routine
        ! procedure :: compute_flux_inf                                 !< Compute inflow flux
    end type

    type, public, abstract, extends(PDE_c) :: PDE_1D_c
        type(BCs_1D_c) :: BCs                                                    !< Boundary conditions
        type(tridiag_matrix_c) :: trans_mat                                   !< Transition matrix T (tridiagonal)
    contains
        procedure :: set_BCs_1D
        procedure :: solve_PDE=>solve_PDE_1D_stat                             !< Solve steady-state PDE
    !> Allocate procedures
        procedure :: allocate_trans_mat=>allocate_trans_mat_1D                              !< Allocate transition matrix
        procedure :: allocate_bd_mat=>allocate_bd_mat_1D                                  !< Allocate boundary matrix
        procedure :: update_trans_mat_1D                             !< Update transition matrix
        procedure :: allocate_arrays_PDE=>allocate_arrays_PDE_1D  !< Allocate all arrays
        ! procedure :: compute_flux_inf                                 !< Compute inflow flux
        !procedure :: compute_rech_mat_PDE=>compute_rech_mat_PDE_1D                                 !< Compute outflow flux
        !procedure :: compute_source_term_PDE=>compute_source_term_PDE_1D                                 !< Compute source/sink term
    end type

    type, public, abstract, extends(PDE_c) :: PDE_2D_c
        type(BCs_2D_c) :: BCs                                                    !< Boundary conditions
        type(pentadiag_matrix_c) :: trans_mat                                 !< Transition matrix T (pentadiagonal)
    contains
        procedure :: set_BCs_2D
        !procedure :: solve_PDE=>solve_PDE_2D_stat                             !< Solve steady-state PDE
    !> Allocate procedures
        procedure :: allocate_trans_mat=>allocate_trans_mat_2D                              !< Allocate transition matrix
        procedure :: allocate_bd_mat=>allocate_bd_mat_2D                                  !< Allocate boundary matrix
        !procedure :: compute_rech_mat_PDE=>compute_rech_mat_PDE_2D                                 !< Compute outflow flux
        !procedure :: compute_source_term_PDE=>compute_source_term_PDE_2D                                 !< Compute source/sink term
    end type
!****************************************************************************************************************************************************
    !> \brief Abstract interfaces for deferred procedures
    !> \details
    !>   These interfaces define the signatures of procedures that must be
    !>   implemented by concrete subclasses of PDE_c.
    abstract interface
        !> \brief Compute transition (transport) matrix
        !> \param[inout] this PDE object
        subroutine compute_trans_mat_PDE(this)
            import PDE_c
            class(PDE_c) :: this                                           !< PDE object
        end subroutine

        ! subroutine compute_F_mat_PDE(this)
        !     import PDE_c
        !     class(PDE_c) :: this                                           !< PDE object
        ! end subroutine
        
        !> \brief Initialize PDE from input files
        !> \param[inout] this PDE object
        !> \param[in] root Root name for input files
        subroutine initialise_PDE(this,path,root,mesh_type)
            import PDE_c
            class(PDE_c) :: this                                           !< PDE object
            character(len=*), intent(in) :: path                              !< Input file root name
            character(len=*), intent(in) :: root                              !< Input file root name
            integer(kind=4), intent(in) :: mesh_type
        end subroutine
        
        
        !> \brief Write PDE results to output file
        !> \param[in] this PDE object
        !> \param[in] root Root name for output files
        !> \param[in] Time_out Output time points
        !> \param[in] output Solution array (space × time)
        subroutine write_PDE(this)
            import PDE_c
            class(PDE_c), intent(in) :: this                               !< PDE object
            ! character(len=*), intent(in) :: root                              !< Output file root name
            ! real(kind=8), intent(in) :: Time_out(:)                           !< Time points for output
            ! real(kind=8), intent(in) :: output(:,:)                           !< Solution (Num_targets × Num_times)
        end subroutine

        subroutine compute_source_term_PDE(this)
            import PDE_c
            class(PDE_c) :: this                               !< PDE object
        end subroutine
        
        !> \brief Compute source/sink term vector
        !> \param[inout] this PDE object
        subroutine compute_source_term_PDE_1D(this)
            import PDE_1D_c
            class(PDE_1D_c) :: this                                           !< PDE object
            !integer(kind=4), intent(in), optional :: k                       !< [OPTIONAL] Time step index
        end subroutine

        subroutine compute_source_term_PDE_2D(this)
            import PDE_2D_c
            class(PDE_2D_c) :: this                                           !< PDE object
            !integer(kind=4), intent(in), optional :: k                       !< [OPTIONAL] Time step index
        end subroutine
        
        !> \brief Compute recharge matrix
        !> \param[inout] this PDE object
        subroutine compute_rech_mat_PDE(this)
            import PDE_c
            class(PDE_c) :: this                                           !< PDE object
            !integer(kind=4), intent(in), optional :: k                       !< [OPTIONAL] Time step index
        end subroutine

        subroutine compute_rech_mat_PDE_1D(this)
            import PDE_1D_c
            class(PDE_1D_c) :: this                                           !< PDE object
            !integer(kind=4), intent(in), optional :: k                       !< [OPTIONAL] Time step index
        end subroutine

        subroutine compute_rech_mat_PDE_2D(this)
            import PDE_2D_c
            class(PDE_2D_c) :: this                                           !< PDE object
        end subroutine
        
        !> \brief Solve the PDE system
        !> \param[inout] this PDE object
        !> \param[in] Time_out Output time points
        !> \param[out] output Solution array
        subroutine solve_PDE(this)
            import PDE_c
            class(PDE_c) :: this                                           !< PDE object
            !real(kind=8), intent(in) :: Time_out(:)                           !< Time points for solution
            !real(kind=8), intent(out) :: output(:,:)                          !< Solution array (Num_targets × Num_times)
        end subroutine

        subroutine allocate_trans_mat(this)
            import PDE_c
            class(PDE_c) :: this                                           !< PDE object
        end subroutine

        subroutine allocate_bd_mat(this)
            import PDE_c
            class(PDE_c) :: this                                           !< PDE object
        end subroutine

        subroutine allocate_arrays_PDE(this)
            import PDE_c
            class(PDE_c) :: this                                           !< PDE object
        end subroutine
    end interface
!****************************************************************************************************************************************************
    !> \brief Concrete procedure interfaces
    !> \details
    !>   External interfaces for non-deferred procedures.
    !>   These are implemented within the module.
    interface
        
        
        
        
        
        
        !> \brief Solve steady-state PDE (EXTERNAL)
        !> \param[inout] this PDE object
        subroutine solve_PDE_1D_stat(this)
            import PDE_1D_c
            class(PDE_1D_c) :: this                                           !< PDE object
        end subroutine
        
        !> \brief Solve PDE and write results (EXTERNAL)
        !> \param[inout] this PDE object
        !> \param[in] Time_out Output time points
        subroutine solve_write_PDE(this,Time_out)
            import PDE_c
            class(PDE_c) :: this                                           !< PDE object
            real(kind=8), intent(in) :: Time_out(:)                           !< Time points for output
        end subroutine
        
        !> \brief Main PDE driver routine (EXTERNAL)
        !> \param[inout] this PDE object
        !> \param[in] root Root name for input/output files
        !> \param[in] mesh_type Integer code for mesh type (e.g., 1 for uniform, 2 for non-uniform)
        subroutine main_PDE(this,path,root,mesh_type)
            import PDE_c
            class(PDE_c) :: this                                           !< PDE object
            character(len=*), intent(in) :: path                              !< File root name
            character(len=*), intent(in) :: root                              !< File root name
            integer(kind=4), intent(in) :: mesh_type
        end subroutine

        ! subroutine set_BCs(this,BCs_obj)
        !     import PDE_c
        !     class(PDE_c) :: this                                           !< PDE object
        !     class(BCs_c), intent(in) :: BCs_obj                               !< Boundary conditions
        ! end subroutine
    end interface
!****************************************************************************************************************************************************
    contains
    
        !> \brief Set spatial discretization
        !> \details
        !>   Associates the PDE with a spatial discretization object.
        !>   The discretization defines the mesh, nodes, and connectivity.
        !> \param[inout] this PDE object
        !> \param[in] spatial_discr_obj Spatial discretization (polymorphic)
        subroutine set_spatial_discr(this,spatial_discr_obj)
            implicit none
            class(PDE_c) :: this                                           !< PDE object
            class(spatial_discr_c), intent(in) :: spatial_discr_obj          !< Spatial discretization
            if (allocated(this%spatial_discr)) deallocate(this%spatial_discr)
            allocate(this%spatial_discr, source=spatial_discr_obj)           !< Deep copy into allocatable
        end subroutine 
        
        !> \brief Set boundary conditions
        !> \details
        !>   Assigns boundary condition object to PDE.
        !>   BC types: Dirichlet, Neumann, Robin, PMF, etc.
        !> \param[inout] this PDE object
        !> \param[in] BCs_obj Boundary conditions object
        subroutine set_BCs_1D(this,BCs_obj)
            implicit none
            class(PDE_1D_c) :: this                                           !< PDE object
            class(BCs_1D_c), intent(in) :: BCs_obj                               !< Boundary conditions
            this%BCs=BCs_obj                                                  !< Copy BC data
        end subroutine

        subroutine set_BCs_2D(this,BCs_obj)
            implicit none
            class(PDE_2D_c) :: this                                           !< PDE object
            class(BCs_2D_c), intent(in) :: BCs_obj                               !< Boundary conditions
            this%BCs=BCs_obj                                                  !< Copy BC data
        end subroutine
        
        !> \brief Allocate transition matrix
        !> \details
        !>   Allocates tridiagonal transition matrix T.
        !>   Size determined by number of discretization targets.
        !> \param[inout] this PDE object
        subroutine allocate_trans_mat_1D(this)
            implicit none
            class(PDE_1D_c) :: this                                           !< PDE object
            call this%trans_mat%allocate_array(this%spatial_discr%Num_targets)  !< Allocate n×n tridiagonal
        end subroutine

        subroutine allocate_trans_mat_2D(this)
            implicit none
            class(PDE_2D_c) :: this                                           !< PDE object
            call this%trans_mat%allocate_array(this%spatial_discr%Num_targets)  !< Allocate n×n tridiagonal
        end subroutine
        
        !> \brief Allocate recharge matrix
        !> \details
        !>   Allocates diagonal recharge matrix R.
        !>   Initialized to zero (no recharge initially).
        !> \param[inout] this PDE object
        subroutine allocate_rech_mat(this)
            implicit none
            class(PDE_c) :: this                                           !< PDE object
            call this%rech_mat%allocate_array(this%spatial_discr%Num_targets)  !< Allocate n×n diagonal
            this%rech_mat%diag=0d0                                            !< Initialize to zero
        end subroutine
        
        !> \brief Allocate boundary matrix
        !> \details
        !>   Allocates boundary matrix B for two boundaries (left and right).
        !>   Initialized to zero (no boundary influence initially).
        !> \param[inout] this PDE object
        subroutine allocate_bd_mat_1D(this)
            implicit none
            class(PDE_1D_c) :: this                                           !< PDE object
            allocate(this%bd_mat(2))                                          !< Two boundaries
            this%bd_mat=0d0                                                   !< Initialize to zero
        end subroutine

        subroutine allocate_bd_mat_2D(this)
            implicit none
            class(PDE_2D_c) :: this                                           !< PDE object
            allocate(this%bd_mat(4))                                          !< Four boundaries
            this%bd_mat=0d0                                                   !< Initialize to zero
        end subroutine
        
        !> \brief Update transition matrix
        !> \details
        !>   Replaces existing transition matrix with new one.
        !>   Used when transport properties change (e.g., time-varying flow).
        !> \param[inout] this PDE object
        !> \param[in] trans_mat New transition matrix
        subroutine update_trans_mat_1D(this,trans_mat)
            implicit none
            class(PDE_1D_c) :: this                                           !< PDE object
            class(tridiag_matrix_c), intent(in)  :: trans_mat                 !< New transition matrix
            this%trans_mat=trans_mat                                          !< Replace matrix
        end subroutine
        
        !> \brief Allocate source term vector
        !> \details
        !>   Allocates source/sink term vector g.
        !>   Initialized to zero (no sources/sinks initially).
        !> \param[inout] this PDE object
        subroutine allocate_source_term_PDE(this)
            implicit none
            class(PDE_c) :: this                                           !< PDE object
            allocate(this%source_term_PDE(this%spatial_discr%Num_targets))    !< Allocate vector
            this%source_term_PDE=0d0                                          !< Initialize to zero
        end subroutine
        
        !> \brief Set solution method
        !> \details
        !>   Selects algorithm for solving PDE:
        !>   - Method 1: Numerical (direct linear solve)
        !>   - Method 2: Eigendecomposition (analytical)
        !>   
        !>   Validates that method is implemented.
        !> \param[inout] this PDE object
        !> \param[in] method Solution method (1 or 2)
        subroutine set_sol_method(this,method)
            implicit none
            class(PDE_c) :: this                                           !< PDE object
            integer(kind=4), intent(in) :: method                             !< Solution method
            if (method<0 .or. method>2) error stop "Solution method not implemented"  !< Validate range
            this%sol_method=method                                            !< Set method
        end subroutine
        
        !> \brief Allocate all PDE arrays (stationary version)
        !> \details
        !>   Convenience routine to allocate all matrix/vector components:
        !>   - Transition matrix
        !>   - Recharge matrix
        !>   - Boundary matrix
        !>   - Source term vector
        !>   
        !>   Used for steady-state problems.
        !> \param[inout] this PDE object
        subroutine allocate_arrays_PDE_1D(this)
            implicit none
            class(PDE_1D_c) :: this                                           !< PDE object
            call this%allocate_trans_mat()                                    !< Allocate T
            call this%allocate_rech_mat()                                     !< Allocate R
            call this%allocate_bd_mat()                                       !< Allocate B
            call this%allocate_source_term_PDE()                              !< Allocate g
        end subroutine

        !> \brief Compute inflow flux from flow rate
        !> \details
        !>   Converts volumetric flow rate (caudal_inf) to flux based on geometry:
        !>   
        !>   For 2D radial coordinates:
        !>   \f[
        !>   \text{flux} = \frac{Q}{2\pi r_{min} b}
        !>   \f]
        !>   Where:
        !>   - Q: volumetric flow rate
        !>   - r_min: inner radius
        !>   - b: aquifer thickness
        !>   
        !>   For Cartesian 1D:
        !>   \f[
        !>   \text{flux} = Q
        !>   \f]
        !>   (assumes unit cross-sectional area or Q already as specific discharge)
        !>   
        !>   Polymorphic dispatch based on mesh type.
        !>
        !> \param[inout] this PDE object
        ! subroutine compute_flux_inf(this)
        ! implicit none
        ! class(PDE_1D_c) :: this                                               !< PDE object
        ! !real(kind=8), intent(in), optional :: phi                            !< Porosity (for specific discharge calculation)
        ! real(kind=8), parameter :: pi=4d0*atan(1d0)                           !< π constant
        ! select type (mesh=>this%spatial_discr)                                !< Polymorphic dispatch
        ! class is (spatial_discr_rad_c)                                        !< Radial coordinates
        !     if (mesh%dim==2) then
        !         !> 2D radial: flux = Q / (2πr*b)
        !         this%BCs%flux_inf=this%BCs%caudal_inf/(2d0*pi*mesh%r_min*mesh%targets(1)%thickness) !> Specific discharge [L/T]
        !     end if
        ! class default                                                         !< Cartesian or other
        !     this%BCs%flux_inf=this%BCs%caudal_inf                             !< flux = Q (assumes unit area)
        !     !print *, this%BCs%flux_inf
        ! !type is (mesh_1D_Euler_homog_c)
        ! !    this%BCs%flux_inf=this%BCs%caudal_inf!/this%spatial_discr%Delta_x
        ! !class is (mesh_1D_Euler_heterog_c)
        ! !    this%BCs%flux_inf=this%BCs%caudal_inf!/this%spatial_discr%Delta_x
        ! !
        ! end select
        ! end subroutine
        
        !> \brief Set dimensionless flag
        !> \details
        !>   Sets flag indicating whether PDE uses dimensionless variables.
        !>   
        !>   Dimensionless formulation:
        !>   - Useful for parameter studies and scaling analysis
        !>   - Reduces number of independent parameters
        !>   - Facilitates comparison with analytical solutions
        !>   - May improve numerical stability
        !>
        !> \param[inout] this PDE object
        !> \param[in] dimless_flag TRUE for dimensionless, FALSE for dimensional
        subroutine set_dimless_flag(this,dimless_flag)
            implicit none
            class(PDE_c) :: this                                           !< PDE object
            logical, intent(in) :: dimless_flag                               !< Dimensionless flag
            this%dimless=dimless_flag                                         !< Set flag
        end subroutine
        
end module PDE_m