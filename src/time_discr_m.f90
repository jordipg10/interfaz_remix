!> \file time_discr_m.f90
!> \brief Defines time discretization classes and routines for numerical time integration
!> \details 
!>   This module provides a comprehensive framework for time discretization in numerical
!>   simulations of transient PDEs (partial differential equations).
!>
!>   Key Features:
!>   - Abstract base class (time_discr_c) defining time discretization interface
!>   - Uniform time stepping (time_discr_homog_c) for constant Δt
!>   - Non-uniform time stepping (time_discr_heterog_c) for adaptive/variable Δt
!>   - Multiple integration methods: Euler explicit/implicit, Crank-Nicolson, RKF45
!>   - Dimensionless time step computation for non-dimensional analysis
!>   - File I/O for time discretization parameters
!>
!>   Integration Methods (int_method):
!>   - 1: Euler Explicit (forward Euler, 1st order, conditionally stable)
!>   - 2: Euler Implicit (backward Euler, 1st order, unconditionally stable)
!>   - 3: Crank-Nicolson (trapezoidal rule, 2nd order, unconditionally stable)
!>   - 4: RKF45 (Runge-Kutta-Fehlberg 4(5), adaptive, embedded error estimate)
!>
!>   Class Hierarchy:
!>   ```
!>   time_discr_c (abstract)
!>       ├── time_discr_homog_c (uniform Δt)
!>       └── time_discr_heterog_c (variable Δt)
!>   ```
!>
!>   Typical Usage:
!>   1. Instantiate concrete class (homog or heterog)
!>   2. Read parameters from file or set programmatically
!>   3. Compute derived quantities (dimensionless time, final time)
!>   4. Query time step values during time integration loop
!>
!>   Mathematical Formulas:
!>   - Final time: t_f = N * Δt (uniform) or t_f = Σ Δt_i (non-uniform)
!>   - Dimensionless time step: Δt_D = Δt / t_c (t_c = characteristic time)
!>   - Dimensionless final time: t_f_D = t_f / t_c
!>
!> \author jordi Petchamé-Guerrero
!> \date 2025
module time_discr_m
    implicit none         !< Require explicit variable declarations
    save                  !< Preserve module variables between calls
    private               !< Default accessibility is private
    !> \brief Abstract base class for time discretization
    !> \details 
    !>   Defines the common interface for all time discretization schemes.
    !>   
    !>   Member Variables:
    !>   - Final_time: Total simulation time
    !>   - Final_time_D: Dimensionless final time [-]
    !>   - Num_time: Number of time steps [-]
    !>   - int_method: Integration method ID (1-4)
    !>
    !>   Implemented Procedures:
    !>   - set_Final_time: Set final simulation time
    !>   - set_Num_time: Set number of time steps
    !>   - set_int_method: Set integration method
    !>   - compute_Final_time: Calculate final time from time steps
    !>   - compute_Final_time_D: Calculate dimensionless final time
    !>   - compute_Num_time: Calculate number of time steps
    !>
    !>   Deferred Procedures (must be implemented by subclasses):
    !>   - read_time_discr: Read time discretization from file
    !>   - get_Delta_t: Get time step value at index k
    !>   - get_Delta_t_D: Get dimensionless time step at index k
    !>   - set_Delta_t: Set time step values
    !>   - compute_dimless_Delta_t: Compute dimensionless time steps
    !>
    !>   Design Pattern:
    !>   Abstract Factory pattern - subclasses provide concrete implementations
    !>   for uniform vs. non-uniform time stepping strategies
    type, public, abstract :: time_discr_c
        real(kind=8) :: theta_t  !< Time integration weighting factor for transport (e.g., for theta-methods)
        real(kind=8) :: theta_r  !< Time integration weighting factor for reactions
        real(kind=8) :: Final_time      !< Final simulation time
        real(kind=8) :: Final_time_D    !< Dimensionless final time = t_f / t_c [-]
        integer(kind=4) :: Num_time     !< Number of time steps [-]
        integer(kind=4) :: int_method   !< Integration method: 1=EE, 2=EfI, 3=CN, 4=RKF45
    contains
        procedure :: set_theta_t                           !< Set time integration weighting factor for transport
        procedure :: set_theta_r                           !< Set time integration weighting factor for reactions
        procedure :: set_Final_time                          !< Set final time value
        procedure :: set_Num_time                            !< Set number of time steps
        procedure :: set_int_method                          !< Set integration method
        procedure :: compute_Final_time                      !< Compute final time from steps
        procedure :: compute_Final_time_D                    !< Compute dimensionless final time
        procedure :: compute_Num_time                        !< Compute number of steps
        procedure(read_time_discr), public, deferred :: read_time_discr           !< Read from file (deferred)
        procedure(get_Delta_t), public, deferred :: get_Delta_t                   !< Get time step (deferred)
        procedure(get_Delta_t_D), public, deferred :: get_Delta_t_D               !< Get dimless step (deferred)
        procedure(set_Delta_t), public, deferred :: set_Delta_t                   !< Set time steps (deferred)
        procedure(compute_dimless_Delta_t), public, deferred :: compute_dimless_Delta_t  !< Compute dimless (deferred)
    end type
!****************************************************************************************************************************************************
    !> \brief Abstract interfaces for deferred procedures
    !> \details 
    !>   Defines the signatures that subclasses must implement.
    !>   These interfaces ensure consistent API across uniform and non-uniform implementations.
    abstract interface
        !> \brief Read time discretization parameters from input file
        !> \param[in,out] this      Time discretization object to populate
        !> \param[in]     filename  Path to input file
        !> \details 
        !>   File format varies by implementation:
        !>   
        !>   Uniform (homog):
        !>   ```
        !>   int_method    ! 1-4
        !>   Delta_t       ! Time step [T]
        !>   Num_time      ! Number of steps
        !>   ```
        !>   
        !>   Non-uniform (heterog):
        !>   ```
        !>   size_t_vec           ! Number of time breakpoints
        !>   t_1 t_2 ... t_n      ! Time values [T]
        !>   n_1 n_2 ... n_{n-1}  ! Steps per interval
        !>   ```
        subroutine read_time_discr(this,filename)
            import time_discr_c                         !< Import abstract class
            implicit none
            class(time_discr_c) :: this                 !< Time discretization object (modified)
            character(len=*), intent(in) :: filename    !< Input filename
        end subroutine
        
        !> \brief Get time step value at specific index
        !> \param[in] this Time discretization object
        !> \param[in] k    Optional time step index (default: 1)
        !> \return Delta_t Time step value [T]
        !> \details 
        !>   For uniform stepping: returns same Δt regardless of k
        !>   For non-uniform: returns Δt(k) at index k
        !>   
        !>   If k not provided, returns first time step
        function get_Delta_t(this,k) result(Delta_t)
            import time_discr_c                         !< Import abstract class
            implicit none
            class(time_discr_c) :: this                 !< Time discretization object
            integer(kind=4), intent(in), optional :: k  !< Time step index (optional)
            real(kind=8) :: Delta_t                     !< Time step value [T]
        end function
        
        !> \brief Get dimensionless time step value at specific index
        !> \param[in] this Time discretization object
        !> \param[in] k    Optional time step index (default: 1)
        !> \return Delta_t_D Dimensionless time step [-]
        !> \details 
        !>   Returns Δt_D = Δt / t_c where t_c is characteristic time
        !>   Must call compute_dimless_Delta_t first to populate Delta_t_D
        function get_Delta_t_D(this,k) result(Delta_t_D)
            import time_discr_c                         !< Import abstract class
            implicit none
            class(time_discr_c) :: this                 !< Time discretization object
            integer(kind=4), intent(in), optional :: k  !< Time step index (optional)
            real(kind=8) :: Delta_t_D                   !< Dimensionless time step [-]
        end function
        
        !> \brief Set time step values programmatically
        !> \param[in,out] this    Time discretization object to modify
        !> \param[in]     Delta_t Array of time step values [T]
        !> \details 
        !>   For uniform: only first element used, must be positive
        !>   For non-uniform: all elements used, all must be positive
        !>   
        !>   Validates input (all Δt > 0) before assignment
        subroutine set_Delta_t(this,Delta_t)
            import time_discr_c                         !< Import abstract class
            implicit none
            class(time_discr_c) :: this                 !< Time discretization object (modified)
            real(kind=8), intent(in) :: Delta_t(:)      !< Time step values [T]
        end subroutine
        
        !> \brief Compute dimensionless time step values
        !> \param[in,out] this      Time discretization object to update
        !> \param[in]     char_time Characteristic time scale [T]
        !> \details 
        !>   Computes Δt_D = Δt / t_c for non-dimensional analysis
        !>   
        !>   Characteristic time examples:
        !>   - Diffusion: t_c = L² / D
        !>   - Advection: t_c = L / v
        !>   - Reaction: t_c = 1 / k
        subroutine compute_dimless_Delta_t(this,char_time)
        import time_discr_c                             !< Import abstract class
        implicit none
        class(time_discr_c) :: this                     !< Time discretization object (modified)
        real(kind=8), intent(in) :: char_time           !< Characteristic time [T]
        end subroutine
    end interface
!****************************************************************************************************************************************************
    !> \brief Concrete class for uniform (constant) time stepping
    !> \details 
    !>   Implements time_discr_c for uniform time steps.
    !>   
    !>   Additional Member Variables:
    !>   - Delta_t: Constant time step value [T]
    !>   - Delta_t_D: Constant dimensionless time step [-]
    !>
    !>   Characteristics:
    !>   - All time steps have same value: Δt_1 = Δt_2 = ... = Δt_N
    !>   - Final time: t_f = N * Δt
    !>   - Simple, efficient, common in many simulations
    !>   - Suitable when solution varies smoothly in time
    !>
    !>   Implementation Details:
    !>   - Scalar storage (not array) for efficiency
    !>   - get_Delta_t returns same value regardless of index k
    !>   Preferred for most problems unless adaptivity required
    type, public, extends(time_discr_c) :: time_discr_homog_c
        real(kind=8) :: Delta_t     !< Uniform time step [T]
        real(kind=8) :: Delta_t_D   !< Dimensionless uniform time step [-]
    contains
        procedure :: read_time_discr=>read_time_discr_homog                 !< Read from file (uniform format)
        procedure :: read_time_discr_WMA                                    !< Read WMA time discretisation from file
        procedure :: set_Delta_t=>set_Delta_t_homog                         !< Set uniform time step
        procedure :: get_Delta_t=>get_Delta_t_homog                         !< Get uniform time step
        procedure :: get_Delta_t_D=>get_Delta_t_D_homog                     !< Get dimless uniform step
        procedure :: compute_dimless_Delta_t=>compute_dimless_Delta_t_homog !< Compute dimless uniform
    end type
    
    !> \brief Concrete class for non-uniform (variable) time stepping
    !> \details 
    !>   Implements time_discr_c for non-uniform time steps.
    !>   
    !>   Additional Member Variables:
    !>   - Delta_t(:): Array of time step values [T]
    !>   - Delta_t_D(:): Array of dimensionless time steps [-]
    !>
    !>   Characteristics:
    !>   - Each time step can have different value: Δt_i ≠ Δt_j
    !>   - Final time: t_f = Σ Δt_i
    !>   - Number of steps: N = size(Delta_t)
    !>   - Useful for adaptive refinement or multi-scale problems
    !>
    !>   File Format:
    !>   Defines time steps via breakpoints and subdivision counts:
    !>   - Specify time values [t_0, t_1, ..., t_n]
    !>   - For each interval [t_i, t_{i+1}], specify number of uniform substeps
    !>   - Automatically constructs Delta_t array
    !>
    !>   Example:
    !>   Time values [0, 1, 10] with [10, 5] steps gives:
    !>   - 10 steps of Δt=0.1 in [0,1]
    !>   - 5 steps of Δt=1.8 in [1,10]
    type, public, extends(time_discr_c) :: time_discr_heterog_c
        real(kind=8), allocatable :: Delta_t(:)     !< Non-uniform time steps [T]
        real(kind=8), allocatable :: Delta_t_D(:)   !< Dimensionless non-uniform steps [-]
    contains
        procedure :: read_time_discr=>read_time_discr_heterog                 !< Read from file (heterog format)
        procedure :: set_Delta_t=>set_Delta_t_heterog                         !< Set non-uniform steps
        procedure :: get_Delta_t=>get_Delta_t_heterog                         !< Get step at index k
        procedure :: get_Delta_t_D=>get_Delta_t_D_heterog                     !< Get dimless step at k
        procedure :: compute_dimless_Delta_t=>compute_dimless_Delta_t_heterog !< Compute dimless array
    end type
!****************************************************************************************************************************************************
    contains
        !> \brief Set final simulation time value
        !> \param[in,out] this       Time discretization object to modify
        !> \param[in]     Final_time Final time value to set [T]
        !> \details 
        !>   Directly assigns the final simulation time.
        !>   No validation performed - assumes user provides valid positive value.
        !>   
        !>   Used when final time is known a priori (common in simulations).
        !>   Alternative: compute final time from Num_time and Delta_t.
        !>
        !>   Example:
        !>   \code{.f90}
        !>   call time_discr%set_Final_time(100.0d0)  ! 100 time units
        !>   \endcode
        subroutine set_Final_time(this,Final_time)
            implicit none
            class(time_discr_c) :: this                 !< Time discretization object (modified)
            real(kind=8), intent(in) :: Final_time      !< Final time [T]
            
            this%Final_time = Final_time                !< Assign final time to object member
        end subroutine
        
        !> \brief Set number of time steps
        !> \param[in,out] this     Time discretization object to modify
        !> \param[in]     Num_time Number of time steps
        !> \details 
        !>   Sets the total number of time steps for the simulation.
        !>   Validates that Num_time ≥ 1 (at least one step required).
        !>   
        !>   Program terminates with error if validation fails.
        !>   
        !>   Used when number of steps is known (e.g., from input file).
        !>   Alternative: compute from Final_time and Delta_t.
        !>
        !>   Example:
        !>   \code{.f90}
        !>   call time_discr%set_Num_time(1000)  ! 1000 time steps
        !>   \endcode
        subroutine set_Num_time(this,Num_time)
            implicit none
            class(time_discr_c) :: this                 !< Time discretization object (modified)
            integer(kind=4) :: Num_time                 !< Number of time steps
            
            if (Num_time < 1) error stop "Number of time steps must be positive"  !< Validate positive count
            this%Num_time = Num_time                    !< Assign to object member
        end subroutine

        subroutine set_theta_t(this,theta)
            implicit none
            class(time_discr_c) :: this                 !< Time discretization object (modified)
            real(kind=8), intent(in) :: theta           !< Time integration weighting factor
            
            if (theta < 0d0 .or. theta > 1d0) error stop "Theta must be in [0,1]"  !< Validate range [0,1]
            this%theta_t = theta                        !< Assign to object member
        end subroutine

        subroutine set_theta_r(this,theta)
            implicit none
            class(time_discr_c) :: this                 !< Time discretization object (modified)
            real(kind=8), intent(in) :: theta           !< Reaction time integration weighting factor
            
            if (theta < 0d0 .or. theta > 1d0) error stop "Theta must be in [0,1]"
            this%theta_r = theta
        end subroutine
        
        !> \brief Set time integration method
        !> \param[in,out] this       Time discretization object to modify
        !> \param[in]     int_method Integration method ID (1-4)
        !> \details 
        !>   Sets the numerical integration method for time stepping:
        !>   - 1: Euler Explicit (forward Euler)
        !>   - 2: Euler Implicit (backward Euler)
        !>   - 3: Crank-Nicolson (trapezoidal)
        !>   - 4: RKF45 (Runge-Kutta-Fehlberg adaptive)
        !>
        !>   Validates that int_method ∈ [1,4].
        !>   Program terminates if invalid method specified.
        !>
        !>   Example:
        !>   \code{.f90}
        !>   call time_discr%set_int_method(3)  ! Crank-Nicolson
        !>   \endcode
        subroutine set_int_method(this,int_method)
            implicit none
            class(time_discr_c) :: this                     !< Time discretization object (modified)
            integer(kind=4), intent(in) :: int_method       !< Integration method ID
            
            if (int_method > 4 .or. int_method < 1) error stop "Integration method not implemented yet"  !< Validate range [1,4]
            this%int_method = int_method                    !< Assign method to object member
        end subroutine
        
        !> \brief Compute final time from time steps
        !> \param[in,out] this Time discretization object to update
        !> \details 
        !>   Calculates final simulation time based on time step data:
        !>   
        !>   For uniform time stepping:
        !>   \f[ t_f = N \cdot \Delta t \f]
        !>   
        !>   For non-uniform time stepping:
        !>   \f[ t_f = \sum_{i=1}^{N} \Delta t_i \f]
        !>
        !>   Uses Fortran's SELECT TYPE for polymorphic dispatch.
        !>   Automatically determines which formula to apply based on concrete type.
        !>
        !>   Prerequisite: Delta_t and Num_time must be set first.
        !>
        !>   Example:
        !>   \code{.f90}
        !>   call time_discr%set_Delta_t([0.1d0])
        !>   call time_discr%set_Num_time(100)
        !>   call time_discr%compute_Final_time()  ! Final_time = 10.0
        !>   \endcode
        subroutine compute_Final_time(this)
            implicit none
            class(time_discr_c) :: this                 !< Time discretization object (modified)
            
            select type (this)                          !< Polymorphic type selection
            type is (time_discr_homog_c)                !< Case: uniform time stepping
                this%Final_time = this%Num_time * this%Delta_t  !< Multiply: t_f = N * Δt
            type is (time_discr_heterog_c)              !< Case: non-uniform time stepping
                this%Final_time = sum(this%Delta_t)     !< Sum all time steps: t_f = Σ Δt_i
            end select
        end subroutine
        
        !> \brief Compute number of time steps from final time and step size
        !> \param[in,out] this Time discretization object to update
        !> \details 
        !>   Calculates number of time steps based on time step data:
        !>   
        !>   For uniform time stepping:
        !>   \f[ N = \text{round}(t_f / \Delta t) \f]
        !>   Uses NINT (nearest integer) to round to whole steps.
        !>   
        !>   For non-uniform time stepping:
        !>   \f[ N = \text{size}(\Delta t) \f]
        !>   Simply counts array elements.
        !>
        !>   Prerequisite:
        !>   - Uniform: Final_time and Delta_t must be set
        !>   - Non-uniform: Delta_t array must be allocated
        !>
        !>   Example:
        !>   \code{.f90}
        !>   call time_discr%set_Final_time(10.0d0)
        !>   call time_discr%set_Delta_t([0.1d0])
        !>   call time_discr%compute_Num_time()  ! Num_time = 100
        !>   \endcode
        subroutine compute_Num_time(this)
            implicit none
            class(time_discr_c) :: this                 !< Time discretization object (modified)
            
            select type (this)                          !< Polymorphic type selection
            type is (time_discr_homog_c)                !< Case: uniform time stepping
                this%Num_time = nint(this%Final_time / this%Delta_t)  !< Round to nearest integer
            type is (time_discr_heterog_c)              !< Case: non-uniform time stepping
                this%Num_time = size(this%Delta_t)      !< Count array elements
            end select
        end subroutine
        
        !> \brief Set uniform time step value
        !> \param[in,out] this    Uniform time discretization object to modify
        !> \param[in]     Delta_t Array of time step values (only first element used) [T]
        !> \details 
        !>   For uniform stepping, only one time step value needed.
        !>   Takes first element of Delta_t array: Δt = Delta_t(1)
        !>   
        !>   Validates that Δt > 0 (strictly positive).
        !>   Program terminates if validation fails.
        !>
        !>   Array parameter used for interface consistency with heterogeneous case.
        !>
        !>   Example:
        !>   \code{.f90}
        !>   call time_discr%set_Delta_t([0.1d0])  ! Set Δt = 0.1 time units
        !>   \endcode
        subroutine set_Delta_t_homog(this,Delta_t)
            implicit none
            class(time_discr_homog_c) :: this           !< Uniform time discretization object (modified)
            real(kind=8), intent(in) :: Delta_t(:)      !< Time step values [T] (only first used)
            
            if (Delta_t(1) < 0d0) error stop "Time step must be positive"  !< Validate Δt > 0
            this%Delta_t = Delta_t(1)                   !< Assign first element to scalar member
        end subroutine
        
        !> \brief Set non-uniform time step values
        !> \param[in,out] this    Non-uniform time discretization object to modify
        !> \param[in]     Delta_t Array of time step values [T]
        !> \details 
        !>   For non-uniform stepping, each time step can differ.
        !>   Assigns entire array: Δt_i = Delta_t(i) for i = 1..N
        !>   
        !>   Validates that all Δt_i > 0 (strictly positive).
        !>   Program terminates if any step is non-positive.
        !>
        !>   Loop checks every element before assignment ensures data integrity.
        !>
        !>   Example:
        !>   \code{.f90}
        !>   real(kind=8) :: steps(3) = [0.1d0, 0.05d0, 0.2d0]
        !>   call time_discr%set_Delta_t(steps)
        !>   \endcode
        subroutine set_Delta_t_heterog(this,Delta_t)
            implicit none
            class(time_discr_heterog_c) :: this         !< Non-uniform time discretization object (modified)
            real(kind=8), intent(in) :: Delta_t(:)      !< Array of time step values [T]
            integer(kind=4) :: i                        !< Loop counter
            
            do i = 1, size(Delta_t)                     !< Loop over all time steps
                if (Delta_t(i) <= 0d0) error stop "Time steps must all be positive"  !< Validate each Δt_i > 0
            end do
            this%Delta_t = Delta_t                      !< Assign entire array to member
        end subroutine
        
        !> \brief Read uniform time discretization parameters from file
        !> \param[in,out] this     Uniform time discretization object to populate
        !> \param[in]     filename Path to input file
        !> \details 
        !>   Reads uniform time stepping parameters from formatted text file.
        !>   
        !>   Expected file format (3 lines):
        !>   ```
        !>   int_method    ! Line 1: Integration method (1-4)
        !>   Delta_t       ! Line 2: Time step size [T]
        !>   Num_time      ! Line 3: Number of time steps
        !>   ```
        !>
        !>   Validation performed:
        !>   - int_method ∈ [1,4]
        !>   - Delta_t > 0
        !>   - Num_time > 0
        !>
        !>   After reading, computes Final_time = Num_time * Delta_t
        !>
        !>   Uses unit 1 for file I/O (assumes unit is available).
        !>   File must exist and be readable.
        !>
        !>   Example file:
        !>   ```
        !>   3
        !>   0.1
        !>   1000
        !>   ```
        !>   Results in: Crank-Nicolson, Δt=0.1 time units, N=1000, t_f=100
        subroutine read_time_discr_homog(this,filename)
            implicit none
            class(time_discr_homog_c) :: this           !< Uniform time discretization object (modified)
            character(len=*), intent(in) :: filename    !< Input filename
            
            character(len=100) :: label                 !< Temporary string for reading lines
            integer :: i                                !< I/O status variable
            real(kind=8) :: theta_t_read, theta_r_read  !< Values read from file

            open(unit=1, file=filename, status='old', action='read')  !< Open file for reading
            do
                read(1,*) label                               !< Read line
                if (trim(label)=='TIME DISCRETISATION') then
                    read(1,*) theta_t_read                      !< Read time weighting factor for transport
                    call this%set_theta_t(theta_t_read)
                    
                    read(1,*) theta_r_read                      !< Read time weighting factor for reactions
                    call this%set_theta_r(theta_r_read)
                    
                    !> Derive int_method from theta_t
                    if (theta_t_read==0d0) then
                        call this%set_int_method(1)             !< Euler explicit
                    else if (theta_t_read==1d0) then
                        call this%set_int_method(2)             !< Euler fully implicit
                    else if (theta_t_read==5d-1) then
                        call this%set_int_method(3)             !< Crank-Nicolson
                    else
                        error stop "Theta must be 0, 0.5, or 1"
                    end if
                    
                    read(1,*) this%Delta_t                      !< Read time step from line 3
                    if (this%Delta_t < 0d0) then                !< Validate Δt > 0
                        error stop "Time step must be positive"
                    end if
                
                    read(1,*) this%Num_time                     !< Read number of steps from line 4
                    if (this%Num_time < 0) then                 !< Validate N > 0
                        error stop "Number of time steps must be positive"
                    end if
                else if (trim(label)=='end') then
                    exit                                        !< End of time discretization block
                else
                    continue
                end if
            end do
            close(1)                                    !< Close file
            
            call this%compute_Final_time()              !< Compute t_f = N * Δt
        end subroutine
        
        !> \brief Read non-uniform time discretization parameters from file
        !> \param[in,out] this     Non-uniform time discretization object to populate
        !> \param[in]     filename Path to input file
        !> \details 
        !>   Reads non-uniform time stepping via time breakpoints and subdivision counts.
        !>   
        !>   Expected file format:
        !>   ```
        !>   size_t_vec              ! Number of time breakpoints
        !>   t_1 t_2 ... t_n         ! Time values [T]
        !>   n_1 n_2 ... n_{n-1}     ! Subdivisions per interval
        !>   ```
        !>
        !>   Algorithm:
        !>   For each interval [t_i, t_{i+1}], create n_i uniform substeps:
        !>   \f[ \Delta t_{\text{interval}} = \frac{t_{i+1} - t_i}{n_i} \f]
        !>
        !>   Example file:
        !>   ```
        !>   3
        !>   0.0  1.0  10.0
        !>   10  5
        !>   ```
        !>   Creates:
        !>   - 10 steps of Δt=0.1 in [0,1]
        !>   - 5 steps of Δt=1.8 in [1,10]
        !>   Total: 15 time steps
        !>
        !>   After construction, computes:
        !>   - Final_time = sum(Delta_t)
        !>   - Num_time = size(Delta_t)
        !>
        !>   Uses unit 2 for file I/O (different from homog to avoid conflict).
        subroutine read_time_discr_heterog(this,filename)
            implicit none
            class(time_discr_heterog_c) :: this         !< Non-uniform time discretization object (modified)
            character(len=*), intent(in) :: filename    !< Input filename
            real(kind=8), allocatable :: t_vec(:)       !< Time breakpoints [T]
            integer(kind=4), allocatable :: n_vec(:)    !< Subdivisions per interval
            integer(kind=4) :: i, j, k, size_t_vec      !< Loop counters and size
            
            open(unit=2, file=filename, status='old', action='read')  !< Open file for reading
            
            read(2,*) size_t_vec                        !< Read number of time breakpoints
            allocate(t_vec(size_t_vec))                 !< Allocate breakpoint array
            read(2,*) t_vec                             !< Read time values [t_0, t_1, ..., t_n]
            
            allocate(n_vec(size_t_vec - 1))             !< Allocate subdivision counts (n-1 intervals)
            read(2,*) n_vec                             !< Read subdivisions per interval
            
            allocate(this%Delta_t(sum(n_vec)))          !< Allocate time step array (total steps)
            
            i = 1                                        !< Initialize time step array index
            k = 1                                        !< Initialize interval index
            do                                           !< Loop over intervals
                this%Delta_t(i) = (t_vec(k+1) - t_vec(k)) / n_vec(k)  !< Compute Δt for interval k
                if (n_vec(k) > 1) then                  !< If multiple substeps in interval
                    do j = 1, n_vec(k) - 1              !< Fill remaining substeps with same Δt
                        this%Delta_t(i + j) = this%Delta_t(i)  !< Copy Δt to subsequent elements
                    end do
                end if
                i = i + n_vec(k)                        !< Advance array index by substep count
                k = k + 1                               !< Move to next interval
                if (k == size_t_vec) exit               !< Exit when all intervals processed
            end do
            
            close(2)                                    !< Close file
            
            this%Final_time = sum(this%Delta_t)         !< Compute final time = Σ Δt_i
            this%Num_time = size(this%Delta_t)          !< Count total number of steps
        end subroutine
        
        !> \brief Get uniform time step value
        !> \param[in] this Uniform time discretization object
        !> \param[in] k    Optional time step index (ignored for uniform)
        !> \return Delta_t Uniform time step [T]
        !> \details 
        !>   Returns the constant time step value for uniform stepping.
        !>   
        !>   Parameter k is ignored (present for interface compatibility).
        !>   All time steps have same value in uniform case.
        !>
        !>   Example:
        !>   \code{.f90}
        !>   dt = time_discr%get_Delta_t()    ! Returns Δt
        !>   dt = time_discr%get_Delta_t(5)   ! Also returns same Δt
        !>   \endcode
        function get_Delta_t_homog(this,k) result(Delta_t)
            implicit none
            class(time_discr_homog_c) :: this           !< Uniform time discretization object
            integer(kind=4), intent(in), optional :: k  !< Time step index (ignored)
            real(kind=8) :: Delta_t                     !< Uniform time step [T]
            
            Delta_t = this%Delta_t                      !< Return scalar time step (k not used)
        end function
        
        !> \brief Get non-uniform time step value at specific index
        !> \param[in] this Non-uniform time discretization object
        !> \param[in] k    Optional time step index (default: 1)
        !> \return Delta_t Time step at index k [T]
        !> \details 
        !>   Returns the time step value at index k.
        !>   
        !>   If k provided: returns Delta_t(k)
        !>   If k omitted: returns Delta_t(1) (first time step)
        !>
        !>   Useful in time integration loops to get current step size.
        !>
        !>   Example:
        !>   \code{.f90}
        !>   dt_1 = time_discr%get_Delta_t()    ! First step
        !>   dt_5 = time_discr%get_Delta_t(5)   ! Fifth step
        !>   \endcode
        function get_Delta_t_heterog(this,k) result(Delta_t)
            implicit none
            class(time_discr_heterog_c) :: this         !< Non-uniform time discretization object
            integer(kind=4), intent(in), optional :: k  !< Time step index (optional)
            real(kind=8) :: Delta_t                     !< Time step at index k [T]
            
            if (present(k)) then                        !< Check if k provided
                Delta_t = this%Delta_t(k)               !< Return Δt at index k
            else                                        !< k not provided
                Delta_t = this%Delta_t(1)               !< Return first time step
            end if
        end function
        
        !> \brief Compute dimensionless uniform time step
        !> \param[in,out] this      Uniform time discretization object to update
        !> \param[in]     char_time Characteristic time scale [T]
        !> \details 
        !>   Computes dimensionless time step for non-dimensional analysis:
        !>   \f[ \Delta t_D = \frac{\Delta t}{t_c} \f]
        !>
        !>   Where:
        !>   - Δt: Physical time step [T]
        !>   - t_c: Characteristic time of problem [T]
        !>   - Δt_D: Dimensionless time step [-]
        !>
        !>   Result stored in this%Delta_t_D
        !>
        !>   Example:
        !>   \code{.f90}
        !>   t_c = L**2 / D  ! Diffusive characteristic time
        !>   call time_discr%compute_dimless_Delta_t(t_c)
        !>   \endcode
        subroutine compute_dimless_Delta_t_homog(this,char_time)
            implicit none
            class(time_discr_homog_c) :: this           !< Uniform time discretization object (modified)
            real(kind=8), intent(in) :: char_time       !< Characteristic time [T]
            
            this%Delta_t_D = this%Delta_t / char_time   !< Compute Δt_D = Δt / t_c
        end subroutine
        
        !> \brief Compute dimensionless non-uniform time steps
        !> \param[in,out] this      Non-uniform time discretization object to update
        !> \param[in]     char_time Characteristic time scale [T]
        !> \details 
        !>   Computes dimensionless time step array for non-dimensional analysis:
        !>   \f[ \Delta t_D(i) = \frac{\Delta t(i)}{t_c} \quad \forall i \f]
        !>
        !>   Fortran's array syntax: performs element-wise division
        !>   All elements divided by same characteristic time
        !>
        !>   Result stored in this%Delta_t_D(:)
        !>
        !>   Example:
        !>   \code{.f90}
        !>   t_c = L / v  ! Advective characteristic time
        !>   call time_discr%compute_dimless_Delta_t(t_c)
        !>   \endcode
        subroutine compute_dimless_Delta_t_heterog(this,char_time)
            implicit none
            class(time_discr_heterog_c) :: this         !< Non-uniform time discretization object (modified)
            real(kind=8), intent(in) :: char_time       !< Characteristic time [T]
            
            this%Delta_t_D = this%Delta_t / char_time   !< Element-wise division: Δt_D(:) = Δt(:) / t_c
        end subroutine
        
        !> \brief Get dimensionless uniform time step value
        !> \param[in] this Uniform time discretization object
        !> \param[in] k    Optional time step index (ignored for uniform)
        !> \return Delta_t_D Dimensionless uniform time step [-]
        !> \details 
        !>   Returns the constant dimensionless time step for uniform stepping.
        !>   
        !>   Parameter k is ignored (present for interface compatibility).
        !>   Must call compute_dimless_Delta_t first to populate Delta_t_D.
        !>
        !>   Returns Δt_D = Δt / t_c (computed earlier)
        !>
        !>   Example:
        !>   \code{.f90}
        !>   call time_discr%compute_dimless_Delta_t(t_c)
        !>   dt_d = time_discr%get_Delta_t_D()  ! Returns Δt_D
        !>   \endcode
        function get_Delta_t_D_homog(this,k) result(Delta_t_D)
            implicit none
            class(time_discr_homog_c) :: this           !< Uniform time discretization object
            integer(kind=4), intent(in), optional :: k  !< Time step index (ignored)
            real(kind=8) :: Delta_t_D                   !< Dimensionless uniform time step [-]
            
            Delta_t_D = this%Delta_t_D                  !< Return scalar dimensionless step (k not used)
        end function
        
        !> \brief Get dimensionless non-uniform time step at specific index
        !> \param[in] this Non-uniform time discretization object
        !> \param[in] k    Optional time step index (default: 1)
        !> \return Delta_t_D Dimensionless time step at index k [-]
        !> \details 
        !>   Returns the dimensionless time step value at index k.
        !>   
        !>   If k provided: returns Delta_t_D(k)
        !>   If k omitted: returns Delta_t_D(1) (first dimensionless step)
        !>
        !>   Must call compute_dimless_Delta_t first to populate Delta_t_D array.
        !>
        !>   Useful for checking CFL numbers or stability criteria.
        !>
        !>   Example:
        !>   \code{.f90}
        !>   call time_discr%compute_dimless_Delta_t(t_c)
        !>   dt_d_5 = time_discr%get_Delta_t_D(5)  ! Fifth dimless step
        !>   \endcode
        function get_Delta_t_D_heterog(this,k) result(Delta_t_D)
            implicit none
            class(time_discr_heterog_c) :: this         !< Non-uniform time discretization object
            integer(kind=4), intent(in), optional :: k  !< Time step index (optional)
            real(kind=8) :: Delta_t_D                   !< Dimensionless time step at k [-]
            
            if (present(k)) then                        !< Check if k provided
                Delta_t_D = this%Delta_t_D(k)           !< Return Δt_D at index k
            else                                        !< k not provided
                Delta_t_D = this%Delta_t_D(1)           !< Return first dimensionless step
            end if
        end function
        
        !> \brief Compute dimensionless final time
        !> \param[in,out] this Time discretization object to update
        !> \param[in]     t_c  Characteristic time scale [T]
        !> \details 
        !>   Computes dimensionless final simulation time:
        !>   \f[ t_f^D = \frac{t_f}{t_c} \f]
        !>
        !>   Where:
        !>   - t_f: Physical final time [T]
        !>   - t_c: Characteristic time scale [T]
        !>   - t_f^D: Dimensionless final time [-]
        !>
        !>   Used for non-dimensional analysis and scaling studies.
        !>   Result stored in this%Final_time_D
        !>
        !>   Prerequisite: Final_time must be set or computed first.
        !>
        !>   Example:
        !>   \code{.f90}
        !>   t_c = L**2 / D  ! Diffusive characteristic time
        !>   call time_discr%compute_Final_time_D(t_c)
        !>   print *, 'Dimensionless final time:', time_discr%Final_time_D
        !>   \endcode
        subroutine compute_Final_time_D(this,t_c)
            implicit none
            class(time_discr_c) :: this                 !< Time discretization object (modified)
            real(kind=8), intent(in) :: t_c             !< Characteristic time [T]
            
            this%Final_time_D = this%Final_time / t_c   !< Compute dimensionless final time
        end subroutine
        
        !> \brief Read WMA time discretisation from a "_WMA_discr.dat" file
        !>
        !> \details Reads the integration method for chemical reactions, the time step size,
        !>          and the number of time steps from the problem-specific WMA discretisation
        !>          file. The file format is:
        !>          - Title line (skipped)
        !>          - Separator line (skipped)
        !>          - 'TIME DISCRETISATION WMA' label
        !>          - Integration method (1: Euler Explicit, 2: Euler fully Implicit, 3: Crank-Nicolson)
        !>          - Time step size (days)
        !>          - Number of time steps
        !>          - Separator line (skipped)
        !>          - 'end'
        !>
        !> @param[in,out] this Homogeneous time discretisation object to configure
        !> @param[in] dir Problem directory path
        !> @param[in] root Root name (prefix) for input file
        !> @param[out] int_method_chem_reacts Integration method for chemical reactions
        subroutine read_time_discr_WMA(this,dir,root,int_method_chem_reacts)
            implicit none
            class(time_discr_homog_c) :: this
            character(len=*), intent(in) :: dir
            character(len=*), intent(in) :: root
            integer(kind=4), intent(out) :: int_method_chem_reacts
            
            integer(kind=4) :: unit
            character(len=256) :: label
            
            unit=58
            open(unit,file=dir//root//'_WMA_discr.dat',status='old',action='read')
            do
                read(unit,*) label
                if (trim(label).eq.'end') then
                    exit
                else if (index(label,'TIME').gt.0) then
                    read(unit,*) int_method_chem_reacts
                    if (int_method_chem_reacts == 1) then
                        call this%set_theta_r(0d0)
                    else if (int_method_chem_reacts == 2) then
                        call this%set_theta_r(1d0)
                    else if (int_method_chem_reacts == 3) then
                        call this%set_theta_r(0.5d0)
                    else
                        error stop "Integration method not implemented"
                    end if
                    this%int_method = int_method_chem_reacts
                    read(unit,*) this%Delta_t
                    if (this%Delta_t < 0d0) error stop "Time step must be positive"
                    read(unit,*) this%Num_time
                    if (this%Num_time < 0) error stop "Number of time steps must be positive"
                    call this%compute_Final_time()
                else
                    continue
                end if
            end do
            close(unit)
        end subroutine
        
end module time_discr_m