!> @file MRMT_m.f90
!> @brief Multi-Rate Mass Transfer (MRMT) module for modeling anomalous reactive transport
!> @details This module implements the Multi-Rate Mass Transfer (MRMT) model for simulating
!> solute reactive transport in heterogeneous porous media with non-equilibrium partitioning between
!> mobile and immobile domains. The MRMT approach represents physical or chemical heterogeneity
!> by dividing the medium into one mobile zone and multiple immobile zones, each characterized
!> by distinct exchange rates and volume fractions.
!>
!> @par MRMT Conceptual Model:
!> The porous medium is partitioned into:
!>   - One mobile zone: advective-dispersive transport
!>   - N immobile zones: diffusion-limited or sorption-controlled storage
!>   - First-order mass exchange between mobile and each immobile zone
!>
!> @par Governing Equations:
!> Mobile zone: φ_m ∂c_m/∂t = -∇·(q*c_m) + ∇·(φ_m*D*∇c_m) - Σ_j φ_j*α_j*(c_m - c_j) + φ_m*R_m
!> Immobile zone j: φ_j ∂c_j/∂t = α_j*φ_j*(c_m - c_j) + φ_j*R_j
!>
!> where:
!>   - φ_m = mobile porosity [-]
!>   - φ_j = φ_im * p_j = immobile zone j porosity [-]
!>   - p_j = probability (volume fraction) of zone j [-]
!>   - α_j = first-order exchange rate for zone j [1/T]
!>   - c_m = mobile concentration [M/L³]
!>   - c_j = immobile zone j concentration [M/L³]
!>   - R_m, R_j = reaction terms in mobile and immobile zones [M/(L³·T)]
!>
!> @see Haggerty & Gorelick (1995), Water Resources Research, 31(10), 2383-2400
!> @see Carrera et al. (1998), Water Resources Research, 34(12), 3303-3314
!> @author Jordi Petchamé-Guerrero
!> @date November 2025

module MRMT_m
    use PDE_model_m, only: PDE_model_c !< Import base PDE model class
    use arrays_m, only: tridiag_matrix_c !< Import tridiagonal matrix class for linear systems
    use mob_zone_m, only: mob_zone_c !< Import mobile zone class definition
    use imm_zone_m, only: imm_zone_c !< Import immobile zone class definition
    use PDE_transient_m, only: PDE_1D_transient_c !< Import transient PDE model class
    implicit none !< Enforce explicit variable declarations for type safety
    private !< Make all entities private by default
    
    !> @brief Multi-Rate Mass Transfer class for non-equilibrium transport modeling
    !> @details Extends the base PDE model class to include mobile/immobile zone partitioning
    !> with multiple exchange rates. This class represents a dual-domain or multi-continuum model
    !> where the mobile zone handles advective-dispersive transport and multiple immobile zones
    !> act as diffusion-limited or sorption-controlled storage domains.
    type, public, extends(PDE_model_c) :: MRMT_c
        type(mob_zone_c) :: mob_zone !< Mobile zone properties: porosity φ_m, concentrations, transport parameters [-]
        integer(kind=4) :: n_imm !< Number of immobile zones N (typically 1-10 zones for multi-rate models) [-]
        type(imm_zone_c), allocatable :: imm_zones(:) !< Array of immobile zones [1:n_imm], each with α_j, φ_j, P_j, c_j [-]
        real(kind=8) :: slope_mem_fct !< Slope of the memory function (dimensionless) [-]
        real(kind=8) :: res_time_max !< Maximum residence time among immobile zones [T]
        real(kind=8) :: imm_por !< Total immobile porosity φ_im (sum over all zones) [-]
    contains
        procedure :: set_n_imm !< Set number of immobile zones N
        procedure :: allocate_imm_zones !< Allocate memory for immobile zones array
        !procedure :: compute_A_mat_conc_mob !< Construct coefficient matrix A for mobile zone implicit system
        !procedure :: compute_b_vec_conc_mob !< Construct right-hand side vector b for mobile zone
        procedure :: compute_conc_imm_MRMT !< Update immobile zone concentrations at new time step
        !procedure :: solve_PDE_EI_Delta_t_homog_MRMT !< Solve MRMT system using Euler implicit with uniform time step
        procedure :: check_imm_zones !< Validate immobile zone parameters (probabilities sum to 1, etc.)
        !procedure :: run_PDE=>run_PDE_MRMT !< Execute complete MRMT simulation with file I/O
    end type

    type, public, extends(MRMT_c) :: MRMT_1D_trans_c
        !> 1D Transient MRMT model extending base MRMT class
        class(PDE_1D_transient_c), pointer :: PDE !< Pointer to 1D transient PDE object (for spatial discretization and transport operators)
    contains
        !> Additional 1D transient-specific methods can be added here
        procedure :: compute_A_mat_conc_mob !< Construct coefficient matrix A for mobile zone implicit system
        procedure :: compute_b_vec_conc_mob !< Construct right-hand side vector b for mobile zone
        procedure :: solve_PDE_EI_Delta_t_homog_MRMT_1D_trans !< Solve MRMT system using Euler implicit with uniform time step for 1D transient case
        procedure :: run_PDE=>run_PDE_MRMT_1D_trans !< Execute complete MRMT simulation with file I/O for 1D transient case
    end type
        !character(len=:), allocatable :: msg !< Error message describing the exception
    
    !> @brief Interface declarations for MRMT subroutines
    !> @details These interfaces define the calling signatures for MRMT-specific computation
    !> routines. The actual implementations are in separate source files for modularity.
    interface
        !> @brief Solve MRMT PDE system using Euler implicit with homogeneous (uniform) time step
        !> @details Advances the coupled mobile-immobile zone system through time using theta-weighted
        !> implicit time integration. Solves tridiagonal linear systems at each time step for mobile
        !> zone, then updates immobile zones.
        !> @param[in,out] this MRMT object with initial conditions and parameters
        !> @param[in] theta Time weighting factor θ ∈ [0,1]: 0=explicit, 0.5=Crank-Nicolson, 1=implicit [-]
        !> @param[in] Time_out Output times array [T]
        !> @param[out] output Concentration output array at specified times [M/L³]
        subroutine solve_PDE_EI_Delta_t_homog_MRMT_1D_trans(this,theta,Time_out,output)
            import MRMT_1D_trans_c !< Import MRMT class definition
            class(MRMT_1D_trans_c) :: this !< MRMT object to evolve in time [-]
            real(kind=8), intent(in) :: theta !< Time weighting factor θ [-]
            real(kind=8), intent(in) :: Time_out(:) !< Output times [T]
            real(kind=8), intent(out) :: output(:,:) !< Concentration output [M/L³]
        end subroutine
        
        !> @brief Construct coefficient matrix A for mobile zone implicit linear system
        !> @details Builds the tridiagonal matrix A accounting for transport (advection, dispersion)
        !> and mass exchange with all immobile zones. Matrix form: A*c_mob^(k+1) = b
        !> @param[in] this MRMT object with transport and exchange parameters
        !> @param[in] theta Time weighting factor θ ∈ [0,1] [-]
        !> @param[in] Delta_t Time step size Δt [T]
        !> @param[out] A_mat Tridiagonal coefficient matrix A [-]
        subroutine compute_A_mat_conc_mob(this,theta,Delta_t,A_mat)
            import MRMT_1D_trans_c, tridiag_matrix_c !< Import MRMT and matrix classes
            implicit none !< Enforce explicit declarations
            class(MRMT_1D_trans_c), intent(in) :: this !< MRMT object with parameters [-]
            real(kind=8), intent(in) :: theta !< Time weighting factor θ [-]
            real(kind=8), intent(in) :: Delta_t !< Time step Δt [T]
            class(tridiag_matrix_c), intent(out) :: A_mat !< Coefficient matrix A for linear system A*c_mob^(k+1)=b [-]
        end subroutine
        
        !> @brief Construct right-hand side vector b for mobile zone implicit system
        !> @details Builds vector b accounting for old time step concentrations (mobile and immobile),
        !> transport operators, and boundary conditions. System form: A*c_mob^(k+1) = b
        !> @param[in] this MRMT object with transport and exchange parameters
        !> @param[in] theta Time weighting factor θ ∈ [0,1] [-]
        !> @param[in] Delta_t Time step size Δt [T]
        !> @param[in] conc_mob_old Mobile zone concentrations at time k: c_mob^k [M/L³]
        !> @param[in] conc_imm_old Immobile zone concentrations at time k: c_imm^k [M/L³]
        !> @param[out] b Right-hand side vector for linear system A*c_mob^(k+1)=b [M/L³]
        subroutine compute_b_vec_conc_mob(this,theta,Delta_t,conc_mob_old,conc_imm_old,b_vec)
            import MRMT_1D_trans_c !< Import MRMT class definition
            implicit none !< Enforce explicit declarations
            class(MRMT_1D_trans_c), intent(in) :: this !< MRMT object with parameters [-]
            real(kind=8), intent(in) :: theta !< Time weighting factor θ [-]
            real(kind=8), intent(in) :: Delta_t !< Time step Δt [T]
            real(kind=8), intent(in) :: conc_mob_old(:) !< Mobile concentrations at time k: c_mob^k [M/L³]
            real(kind=8), intent(in) :: conc_imm_old(:) !< Immobile concentrations at time k: c_imm^k [M/L³]
            real(kind=8), intent(out) :: b_vec(:) !< Right-hand side vector b [M/L³]
        end subroutine
        
        !> @brief Update immobile zone concentrations at new time step
        !> @details Computes c_imm^(k+1) for all immobile zones using theta-weighted scheme and
        !> newly computed mobile zone concentration c_mob^(k+1). Each zone updates independently.
        !> @param[in] this MRMT object with exchange parameters
        !> @param[in] theta Time weighting factor θ ∈ [0,1] [-]
        !> @param[in] conc_imm_old Immobile concentrations at time k: c_imm^k [M/L³]
        !> @param[in] conc_mob_old Mobile concentration at time k: c_m^k [M/L³]
        !> @param[in] conc_mob_new Mobile concentration at time k+1: c_m^(k+1) [M/L³]
        !> @param[in] Delta_t Time step size Δt [T]
        !> @param[out] conc_imm_new Immobile concentrations at time k+1: c_imm^(k+1) [M/L³]
        subroutine compute_conc_imm_MRMT(this,theta,conc_imm_old,conc_mob_old,conc_mob_new,Delta_t,conc_imm_new)
            import MRMT_c !< Import MRMT class definition
            implicit none !< Enforce explicit declarations
            class(MRMT_c), intent(in) :: this !< MRMT object with parameters [-]
            real(kind=8), intent(in) :: theta !< Time weighting factor θ [-]
            real(kind=8), intent(in) :: conc_imm_old(:) !< Immobile concentrations at time k: c_imm^k [M/L³]
            real(kind=8), intent(in) :: conc_mob_old(:) !< Mobile concentration at time k: c_m^k [M/L³]
            real(kind=8), intent(in) :: conc_mob_new(:) !< Mobile concentration at time k+1: c_m^(k+1) [M/L³]
            real(kind=8), intent(in) :: Delta_t !< Time step Δt [T]
            real(kind=8), intent(out) :: conc_imm_new(:) !< Immobile concentrations at time k+1: c_imm^(k+1) [M/L³]
        end subroutine
    end interface !< End of interface declarations
    
contains !< Begin module procedure implementations
        
        !> @brief Set number of immobile zones
        !> @details Assigns the number of immobile zones N for the MRMT model. This must be called
        !> before allocate_imm_zones. Typical values range from 1 (dual-porosity) to ~10 (multi-rate).
        !> @param[in,out] this MRMT object to configure
        !> @param[in] n_imm Number of immobile zones N (must be > 0) [-]
        subroutine set_n_imm(this,n_imm)
            implicit none !< Enforce explicit declarations
            class(MRMT_c), intent(inout) :: this !< MRMT object to configure [-]
            integer(kind=4), intent(in) :: n_imm !< Number of immobile zones N (typically 1-10) [-]
            this%n_imm=n_imm !< Assign number of immobile zones to MRMT object [-]
        end subroutine
        
        !> @brief Allocate memory for immobile zones array
        !> @details Allocates the imm_zones array with size n_imm. Must be called after set_n_imm
        !> and before assigning individual zone properties (exchange rates, porosities, probabilities).
        !> @param[in,out] this MRMT object with n_imm already set
        subroutine allocate_imm_zones(this)
            implicit none !< Enforce explicit declarations
            class(MRMT_c), intent(inout) :: this !< MRMT object with n_imm defined [-]
            allocate(this%imm_zones(this%n_imm)) !< Allocate array for N immobile zones [1:n_imm]
        end subroutine
        
        !> @brief Validate immobile zone parameters
        !> @details Checks that immobile zone probabilities are physically meaningful:
        !>   1. Each P_j ∈ [0,1] (valid probability/volume fraction)
        !>   2. Σ P_j = 1 (probabilities sum to unity, accounting for all immobile volume)
        !> Stops execution with error message if validation fails.
        !> @param[in] this MRMT object with immobile zones configured
        !> @note Uses tolerance ε = 10^-12 for floating-point comparison of probability sum
        subroutine check_imm_zones(this)
            implicit none !< Enforce explicit declarations
            class(MRMT_c), intent(in) :: this !< MRMT object with immobile zones to validate [-]
            integer(kind=4) :: i !< Loop counter for iterating over immobile zones [-]
            real(kind=8) :: prob_tot !< Accumulator for total probability Σ P_j (should equal 1.0) [-]
            real(kind=8), parameter :: epsilon=1d-12 !< Tolerance for floating-point comparison (10^-12) [-]
            
            prob_tot = 0d0 !< Initialize probability accumulator [-]
            !> Loop over all immobile zones to validate and accumulate probabilities
            do i=1,this%n_imm !< Iterate through all N immobile zones
                if (this%imm_zones(i)%prob<0d0 .or. this%imm_zones(i)%prob>1d0) &
                    error stop "Probability must be in range [0,1]" !< Ensure 0 ≤ P_j ≤ 1
                prob_tot = prob_tot + this%imm_zones(i)%prob !< Accumulate probability sum [-]
            end do
            if (abs(prob_tot-1d0)>=epsilon) &
                error stop "Probabilities must sum to 1" !< Ensure |Σ P_j - 1| < ε
        end subroutine
        
        !> @brief Execute complete MRMT simulation with file I/O
        !> @details Main driver routine for MRMT simulations. Performs the following steps:
        !>   1. Initialize PDE problem (read input, set up spatial/temporal discretization)
        !>   2. Determine time integration method (Euler explicit θ=0, implicit θ=1, Crank-Nicolson θ=0.5)
        !>   3. Solve coupled MRMT system through time
        !>   4. Write mobile and immobile zone concentrations to output file
        !> @param[in,out] this MRMT object to simulate
        !> @param[in] root Root filename for input/output files (output: root_MRMT.out)
        subroutine run_PDE_MRMT_1D_trans(this,dir,root,mesh_type)
            implicit none !< Enforce explicit declarations
            class(MRMT_1D_trans_c), intent(inout) :: this !< MRMT object to simulate [-]
            character(len=*), intent(in) :: dir !< Full filename for input/output operations [-]
            character(len=*), intent(in) :: root !< Root filename for I/O (e.g., "problem" → "problem_MRMT.out")
            integer(kind=4), intent(in) :: mesh_type !< Integer code for mesh type (e.g., 1 for uniform, 2 for non-uniform)
            
            integer(kind=4) :: i !< Loop counter for output writing [-]
            real(kind=8) :: theta !< Time weighting factor θ determined by integration method [-]
            real(kind=8), allocatable :: MRMT_output(:,:) !< Output array for concentrations at target locations [M/L³]
            
            !> Initialize PDE problem: read input file, set up discretization, initial conditions
            call this%PDE%main_PDE(dir,root,mesh_type) !< Read input using root filename, initialize spatial and temporal discretization
            
            !> Type-check PDE to ensure it is a 1D transient problem (required for MRMT)
            !select type (PDE=>this%PDE) !< Dynamic type checking of PDE object
            !class is (PDE_1D_transient_c) !< MRMT requires 1D transient PDE formulation
                !> Determine time integration method from PDE configuration
                if (this%PDE%time_discr%int_method==1) then !< Method 1: Euler Explicit
                    theta=0d0 !< Explicit scheme: θ = 0 (forward Euler)
                else if (this%PDE%time_discr%int_method==3) then !< Method 3: Euler Implicit
                    theta=1d0 !< Fully implicit scheme: θ = 1 (backward Euler)
                else if (this%PDE%time_discr%int_method==4) then !< Method 4: Crank-Nicolson
                    theta=5d-1 !< Crank-Nicolson scheme: θ = 0.5 (trapezoidal rule)
                end if
                !> Allocate output array: (number of target locations) × (initial + final time)
                allocate(MRMT_output(this%PDE%spatial_discr%Num_targets,2)) !< Allocate for targets at t=0 and t=Final [M/L³]
                
                !> Solve MRMT system from t=0 to final time using selected time integration method
                call this%solve_PDE_EI_Delta_t_homog_MRMT_1D_trans(theta,[0d0,this%PDE%time_discr%Final_time],MRMT_output) !< Solve with time array [t_0, t_final]
                
                !> Write results to output file
                open(unit=46,file=root//'_MRMT.out',status='unknown') !< Open output file with MRMT suffix
                write(46,"(2x,'Mobile concentrations:',/)") !< Write header for mobile zone section
                do i=1,this%PDE%spatial_discr%Num_targets !< Loop over all target output locations
                    write(46,*) this%mob_zone%conc(i) !< Write mobile concentration at target i [M/L³]
                end do !< End mobile concentrations output
                write(46,"(/,2x,'Immobile concentrations:',/)") !< Write header for immobile zones section
                do i=1,this%n_imm !< Loop over all immobile zones
                    write(46,*) this%imm_zones(i)%conc !< Write immobile zone i concentration [M/L³]
                end do !< End immobile concentrations output
                close(46) !< Close output file
            !end select !< End PDE type selection
        end subroutine !< End of run_PDE_MRMT subroutine
end module !< End of MRMT_m module