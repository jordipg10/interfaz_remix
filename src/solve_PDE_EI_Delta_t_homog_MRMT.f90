!> @file solve_PDE_EI_Delta_t_homog_MRMT.f90
!> @brief Solve 1D transient PDE with homogeneous time step using Euler Implicit and MRMT
!> @details This subroutine solves the one-dimensional transient partial differential equation
!> for solute transport with Multi-Rate Mass Transfer (MRMT) using a theta-weighted implicit
!> Euler method with constant (homogeneous) time stepping.
!>
!> @par Algorithm Overview:
!> For each time step k:
!>   1. Construct coefficient matrix A (depends on theta, Δt, transport, and MRMT parameters)
!>   2. Construct RHS vector b (from old mobile and immobile concentrations)
!>   3. Solve linear system: A * c_mob^(k+1) = b using Thomas algorithm (tridiagonal solver)
!>   4. Update immobile zone concentrations: c_imm^(k+1) = f(c_imm^k, c_mob^k, c_mob^(k+1))
!>   5. Write results to binary output file
!>   6. Check if current time matches output time and store results
!>   7. Update old concentrations for next time step
!>
!> @par Temporal Discretization:
!> Uses theta-weighted scheme:
!>   - theta = 0: Fully explicit (forward Euler)
!>   - theta = 0.5: Crank-Nicolson (second-order accurate)
!>   - theta = 1: Fully implicit (backward Euler, unconditionally stable)
!>
!> @par Output:
!>   - Binary file "conc_binary_EI_MRMT.txt" contains all intermediate time steps
!>   - output array contains concentrations at requested Time_out values
!>   - Format: [mobile concentrations (n values), immobile concentrations (n_imm values)]
!>
!> @see compute_A_mat_conc_mob For coefficient matrix construction
!> @see compute_b_conc_mob For right-hand side vector construction
!> @see compute_conc_imm_MRMT For immobile zone concentration updates
!> @see Thomas For tridiagonal linear system solver
!> @author Generated documentation
!> @date November 2025

!> @brief Solve 1D transient PDE with MRMT using theta-weighted Euler Implicit method
!> @param[in,out] this MRMT object containing PDE, mobile zone, immobile zones, and transport matrix
!> @param[in] theta Time weighting factor θ ∈ [0,1]: 0=explicit, 0.5=Crank-Nicolson, 1=implicit [-]
!> @param[in] Time_out Array of output times at which to store results [T]
!> @param[out] output Concentration matrix: output(spatial_point, time_index) = [c_mob; c_imm] [M/L³]
subroutine solve_PDE_EI_Delta_t_homog_MRMT_1D_trans(this,theta,Time_out,output)
    use MRMT_m, only: MRMT_1D_trans_c !< Import Multi-Rate Mass Transfer class with mobile/immobile zone structure
    use metodos_sist_lin_m, only: Thomas, tridiag_matrix_c !< Import Thomas algorithm (tridiagonal solver) and matrix class
    use diffusion_transient_m, only: diffusion_1D_transient_c !< Import 1D transient diffusion class for PDE type checking
    use time_discr_m, only: time_discr_homog_c !< Import homogeneous time discretization class (constant Δt)
    implicit none !< Enforce explicit variable declarations for type safety
    
    !> Input/Output Variables
    class(MRMT_1D_trans_c) :: this !< MRMT object containing PDE (with spatial_discr, BCs, conc_init, time_discr), mob_zone, imm_zones, n_imm [-]
    real(kind=8), intent(in) :: theta !< Time weighting factor θ for temporal discretization: 0 ≤ θ ≤ 1 [-]
    real(kind=8), intent(in) :: Time_out(:) !< Array of output times at which to record concentrations [T]
    real(kind=8), intent(out) :: output(:,:) !< Output concentration matrix (n+n_imm × Num_output): output(:,i) at Time_out(i) [M/L³]

    !> Local Variables - Integers
    integer(kind=4) :: n !< Number of spatial discretization points (targets) in mobile zone [-]
    integer(kind=4) :: i !< Loop counter for iterating over spatial points or immobile zones [-]
    integer(kind=4) :: icol !< Column index for output array (tracks which Time_out we're at) [-]
    integer(kind=4) :: k !< Time step counter (1 to Num_time) [-]
    integer(kind=4) :: Num_output !< Number of output time points (size of Time_out array) [-]
    
    !> Local Variables - Reals
    real(kind=8) :: Time !< Current simulation time [T]
    real(kind=8), parameter :: epsilon=1d-9 !< Tolerance for time comparison (determines if Time ≈ Time_out) [-]
    real(kind=8), parameter :: tol_Thomas=1d-9 !< Convergence tolerance for Thomas algorithm (tridiagonal solver) [-]
    
    !> Local Variables - Allocatable Arrays
    real(kind=8), allocatable :: conc_mob_old(:) !< Mobile zone concentrations at time k: c_mob^k (n spatial points) [M/L³]
    real(kind=8), allocatable :: conc_mob_new(:) !< Mobile zone concentrations at time k+1: c_mob^(k+1) (n spatial points) [M/L³]
    real(kind=8), allocatable :: conc_imm_init(:) !< Initial immobile zone concentrations (n_imm zones, typically zero) [M/L³]
    real(kind=8), allocatable :: conc_imm_old(:) !< Immobile zone concentrations at time k: c_imm^k (n_imm zones) [M/L³]
    real(kind=8), allocatable :: conc_imm_new(:) !< Immobile zone concentrations at time k+1: c_imm^(k+1) (n_imm zones) [M/L³]
    real(kind=8), allocatable :: b(:) !< Right-hand side vector for linear system A*c_mob^(k+1)=b (n spatial points) [M/L³]
    
    !> Local Variables - Matrices
    type(tridiag_matrix_c) :: A_mat !< Coefficient matrix A for implicit linear system (tridiagonal: sub, diag, super) [-]
    type(tridiag_matrix_c) :: A_mat_ODE !< Transport operator matrix from PDE (used for eigenvalue analysis in commented code) [1/T]
    
    !> Get number of spatial discretization points from PDE spatial discretization
    n=this%PDE%spatial_discr%Num_targets !< Extract number of targets (spatial mesh points) from PDE object [-]
    
    !> Type-check PDE object to ensure it's a 1D transient diffusion problem
    select type (PDE=>this%PDE) !< Use type-select to access specific PDE type methods and attributes
    class is (diffusion_1D_transient_c) !< Only proceed if PDE is a 1D transient diffusion object
        !> Type-check time discretization to ensure homogeneous (constant) time stepping
        select type (time_discr=>PDE%time_discr) !< Access time discretization from PDE object
        type is (time_discr_homog_c) !< Only proceed if time discretization has constant Δt (homogeneous)
            
            !> Initialize mobile zone concentrations from PDE initial conditions
            conc_mob_old=PDE%conc_init !< Copy initial concentration distribution to mobile zone old concentrations [M/L³]
            
            !> Apply outlet boundary condition to last spatial point
            conc_mob_old(n)=PDE%BCs%conc_out !< Set downstream boundary concentration (Dirichlet BC at outlet) [M/L³]
            
            !> Configure MRMT model: set number of immobile zones equal to number of spatial points
            call this%set_n_imm(n) !< Each spatial target gets one associated immobile zone (1:1 correspondence) [-]
            
            !> Set mobile zone porosity (hard-coded value, should be parameterized in future)
            this%mob_zone%mob_por=4d-1 !< Assign mobile zone porosity φ_m = 0.4 (40% pore space) [-]
            
            !> Allocate all concentration and auxiliary arrays based on spatial/zone dimensions
            allocate(conc_mob_new(n),conc_imm_init(this%n_imm),conc_imm_old(this%n_imm),conc_imm_new(this%n_imm),b(n)) !< Allocate: n mobile points, n_imm immobile zones, n RHS vector [M/L³]
            
            !> Initialize immobile zone concentrations to zero (no initial mass in immobile zones)
            conc_imm_init=0d0 !< Set all initial immobile concentrations c_imm^0 = 0 [M/L³]
            
            !> Copy initial immobile concentrations to old concentrations for first time step
            conc_imm_old=conc_imm_init !< Initialize c_imm^0 = conc_imm_init for time stepping [M/L³]
            
            !> Allocate memory for immobile zones array in MRMT object
            call this%allocate_imm_zones() !< Allocate this%imm_zones(1:n_imm) array for storing immobile zone properties [-]
            
            !> Compute transport operator matrix for PDE (used for eigenvalue-based MRMT parameterization)
            call PDE%compute_A_mat_ODE_1D(A_mat_ODE) !< Construct transport matrix A_ODE (diffusion/advection operator) [1/T]
            
            !> Compute eigenvalues of transport operator (for potential MRMT rate distribution)
            call A_mat_ODE%compute_eigenvalues() !< Calculate eigenvalues λ_i of transport matrix (related to relaxation times) [1/T]
            
            !> Validate immobile zone parameters (check normalization, positivity, etc.)
            call this%check_imm_zones() !< Verify Σ P_i = 1, α_i > 0, φ_i > 0, etc. [-]
            
            !> Compute coefficient matrix A for implicit linear system (depends on theta, Δt, MRMT)
            call this%compute_A_mat_conc_mob(theta,time_discr%Delta_t,A_mat) !< Construct A matrix: A*c_mob^(k+1)=b [-]
            
            !> Get number of output time points from Time_out array size
            Num_output=size(Time_out) !< Extract length of Time_out array (number of times to record) [-]
            
            !> Open binary output file for writing all intermediate concentration profiles
            open(unit=0,file="conc_binary_EI_MRMT.txt",form="unformatted",access="sequential",status="unknown") !< Binary file for efficient storage of all time steps
            
            !> Initialize output column index to 1 (first requested output time)
            icol=1 !< Start at first column of output array (Time_out(1)) [-]
            
            !> Initialize simulation time to zero (start of simulation)
            Time=0 !< Set current time t = 0 [T]
            
            !> Check if initial time (t=0) matches first output time and record if so
            if (abs(Time-Time_out(icol))<epsilon) then !< Test if |t - Time_out(1)| < ε (essentially t ≈ Time_out(1))
                output(:,icol)=[conc_mob_old,conc_imm_old] !< Store initial concentrations: [c_mob^0; c_imm^0] [M/L³]
                icol=icol+1 !< Advance to next output time index [-]
            end if !< End initial time output check
            
            !> MAIN TIME-STEPPING LOOP: Advance solution from t=0 to t=T_final
            do k=1,time_discr%Num_time !< Loop over all time steps k = 1, 2, ..., Num_time
                
                !> Update current simulation time
                Time=k*time_discr%Delta_t !< Calculate current time: t^k = k * Δt [T]
                
                !> Write current time and concentrations to binary output file
                write(0) Time, [conc_mob_old,conc_imm_old] !< Binary write: (time, [c_mob^k; c_imm^k]) for all intermediate steps
                
                !> STEP 1: Construct right-hand side vector b for linear system
                !> b depends on old mobile and immobile concentrations, theta, and Δt
                call this%compute_b_vec_conc_mob(theta,time_discr%Delta_t,conc_mob_old,conc_imm_old,b) !< Compute b = X*c_mob^k - f [M/L³]
                
                !> STEP 2: Solve linear system A * c_mob^(k+1) = b using Thomas algorithm
                !> Thomas algorithm is optimal O(n) solver for tridiagonal systems
                call Thomas(A_mat,b,tol_Thomas,conc_mob_new) !< Solve for new mobile concentrations c_mob^(k+1) [M/L³]
                
                !> STEP 3: Update immobile zone concentrations using analytical MRMT formula
                !> c_imm^(k+1) = f(c_imm^k, c_mob^k, c_mob^(k+1), α, θ, Δt)
                call this%compute_conc_imm_MRMT(theta,conc_imm_old,conc_mob_old,conc_mob_new,time_discr%Delta_t,conc_imm_new) !< Compute c_imm^(k+1) from MRMT exchange [M/L³]
                
                !> Check if current time matches a requested output time
                if (abs(Time-Time_out(icol))<epsilon) then !< Test if |t^k - Time_out(icol)| < ε (time match within tolerance)
                    output(:,icol)=[conc_mob_new,conc_imm_new] !< Store concentrations at output time: [c_mob^k; c_imm^k] [M/L³]
                    icol=icol+1 !< Advance to next requested output time [-]
                    
                    !> Check if all requested output times have been recorded
                    if (icol>Num_output) then !< Test if we've exceeded the number of output times
                        write(*,*) "Reached Num_output" !< Inform user that all outputs have been recorded
                        exit !< Exit time loop early (all requested outputs obtained)
                    end if !< End early exit check
                end if !< End output time recording
                
                !> Update old concentrations for next time step (k → k+1)
                conc_mob_old=conc_mob_new !< Set c_mob^k := c_mob^(k+1) for next iteration [M/L³]
                conc_imm_old=conc_imm_new !< Set c_imm^k := c_imm^(k+1) for next iteration [M/L³]
                
                !> Compute total concentration production (for mass balance or analysis)
                call PDE%prod_total_conc(A_mat_ODE,Time) !< Calculate production rate or total mass at current time
            end do !< End main time-stepping loop
            
            !> Store final mobile zone concentration in MRMT object
            this%mob_zone%conc=conc_mob_new !< Copy final mobile zone concentration c_mob^final to MRMT mobile zone [M/L³]
            
            !> Store final immobile zone concentrations in MRMT object
            do i=1,this%n_imm !< Loop over all immobile zones to store final concentrations
                this%imm_zones(i)%conc=conc_imm_new(i) !< Copy final immobile zone i concentration c_imm_i^final to MRMT immobile zone i [M/L³]
            end do !< End loop over immobile zones
            
            !> Deallocate all temporary concentration arrays to free memory
            deallocate(conc_mob_old,conc_mob_new,conc_imm_old,conc_imm_new) !< Free memory for concentration arrays (no longer needed after simulation)
            
            !> Close binary output file
            close(0) !< Close "conc_binary_EI_MRMT.txt" after writing all time steps
            
        end select !< End time discretization type select (time_discr_homog_c)
    end select !< End PDE type select (diffusion_1D_transient_c)
end subroutine !< End of solve_PDE_EI_Delta_t_homog_MRMT subroutine 