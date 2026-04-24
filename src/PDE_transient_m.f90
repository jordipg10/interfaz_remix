!> \file PDE_transient_m.f90
!> \brief Transient (time-dependent) PDE module
!> \details
!>   Extends PDE_m to handle time-dependent partial differential equations.
!>   
!>   Governing equation:
!>   \f[
!>   \mathbf{F} \frac{dc}{dt} = \mathbf{T} c + \mathbf{R} c_{rech} + \mathbf{B} c_{bd}
!>   \f]
!>   
!>   Where:
!>   - \f$\mathbf{F}\f$: Storage matrix (diagonal) - represents accumulation
!>   - \f$c\f$: State variable (concentration, temperature, etc.)
!>   - \f$\mathbf{T}\f$: Transition (transport) matrix - advection & dispersion
!>   - \f$\mathbf{R}\f$: Recharge matrix - sink/source terms
!>   - \f$\mathbf{B}\f$: Boundary matrix - boundary conditions
!>   
!>   θ-method discretization:
!>   \f[
!>   \mathbf{A} c^{k+1} = \mathbf{X} c^k + \theta_t*\mathbf{Y}^{k+1} c_{rech}^{k+1} + (1-\theta_t)\mathbf{Y}^k c_{rech}^k + \theta_t*\mathbf{Z}^{k+1} c_{bd}^{k+1} + (1-\theta_t)*\mathbf{Z}^k c_{bd}^k
!>   \f]
!>   
!>   Where:
!>   - \f$\mathbf{A}\f$: Implicit linear system matrix (current time step)
!>   - \f$\mathbf{X}\f$: Explicit linear system matrix (previous time step)
!>   - \f$\mathbf{Y}\f$: Explicit linear system matrix for recharge
!>   - \f$\mathbf{Z}\f$: Explicit linear system matrix for boundary
!>   - θ: Implicitness parameter (0=explicit Euler, 0.5=Crank-Nicolson, 1=fully implicit Euler)
!>   
!>   Time discretization schemes:
!>   - Explicit Euler (θ=0): Conditionally stable, small time steps
!>   - Fully Implicit Euler (θ=1): Unconditionally stable, large time steps allowed
!>   - Crank-Nicolson (θ=0.5): Second-order accurate, unconditionally stable
!>   - RKF45: Adaptive Runge-Kutta-Fehlberg with error control
!>   
!>   Water Mixing Approach (WMA):
!>   - Interprets numerical solution as physical water mixing
!>   - mixing_ratios: Proportions of waters mixing at each node
!>   - Enables operator splitting for reactive transport
!>   - Separates transport (mixing) from reactions
!>
!> \author jordi Petchamé-Guerrero
!> \date October 2025

module PDE_transient_m
    use PDE_m, only: PDE_c, PDE_1D_c, PDE_2D_c                                     !< Base PDE module
    use time_discr_m, only: time_discr_c, time_discr_homog_c, time_discr_heterog_c  !< Time discretization
    use char_params_m, only: char_params_c                                   !< Characteristic parameters
    use arrays_m, only: diag_matrix_c, tridiag_matrix_c, pentadiag_matrix_c, real_array_c, int_array_c  !< Matrix types
    use BCs_m, only: BCs_1D_c, BCs_2D_c                             !< Boundary conditions
    use spatial_discr_2D_m, only: mesh_2D_Euler_homog_c  !< Needed for 2D band-size correction in allocate_trans_mat_2D_trans
    implicit none
    save
    private    
    !> \brief Transient 1D PDE subclass
    !> \details
    !>   Extends PDE_1D_c to handle time-dependent problems.
    !>   
    !>   Key additional features:
    !>   - Time discretization (homogeneous or heterogeneous time steps)
    !>   - Storage matrix F (accumulation term)
    !>   - θ-method matrices (A, X, Y, Z) for time integration
    !>   - Water Mixing Approach (WMA) for reactive transport
    !>   - Adaptive time stepping (RKF45)
    !>   
    !>   Matrix components:
    !>   - F_mat: Storage/accumulation (diagonal)
    !>   - X_mat: Explicit part from previous time (tridiagonal)
    !>   - Y_mat: Explicit recharge mixing (diagonal)
    !>   - Z_mat: Explicit boundary mixing (size 2 vector)
    !>   - A_mat: Implicit part for current time (tridiagonal)
    !>   - A_mat_inv: Inverse of A_mat for explicit methods
    !>   
    !>   Water Mixing Approach (WMA):
    !>   - mixing_ratios_conc: Concentrations of mixing waters
    !>   - mixing_ratios_R: Reaction amounts for mixing waters
    !>   - mix_conc_indices: Connectivity for water mixing
    !>   - Enables operator splitting: Transport → Reaction → Transport
    !>   - Preserves mass balance exactly
    !>   
    !>   Time integration methods:
    !>   - Explicit Euler (EE): Fast but conditionally stable
    !>   - Implicit Euler (EI): Slower but unconditionally stable
    !>   - RKF45: Adaptive with error control
    type, public, abstract, extends(PDE_c) :: PDE_transient_c
        class(time_discr_c), allocatable :: time_discr                            !< Time discretization (polymorphic)
        type(diag_matrix_c) :: F_mat                                          !< Storage matrix (diagonal)
        type(diag_matrix_c) :: F_mat_prev                                       !< Explicit domain mixing ratios (previous time step)
        type(diag_matrix_c) :: Y_mat                                          !< Explicit recharge mixing ratios
        type(diag_matrix_c) :: Y_mat_prev                                         !< Explicit recharge mixing ratios (previous time step)
        real(kind=8), allocatable :: Z_mat(:)                                 !< Explicit boundary mixing ratios
        real(kind=8), allocatable :: Z_mat_prev(:)                            !< Explicit boundary mixing ratios (previous time step)
        real(kind=8), allocatable :: A_mat_inv(:,:)                           !< Inverse of A_mat (for implicit methods)
        real(kind=8), allocatable :: f_vec(:)                                 !< Independent term: f = Y*c_rech + Z*c_bd
        type(real_array_c) :: mixing_ratios_conc                              !< WMA mixing ratios for concentrations (any dimension)
        type(real_array_c) :: mixing_ratios_R                                !< WMA mixing ratios for reaction amounts
        type(real_array_c) :: mixing_ratios_R_init                           !< WMA initial mixing ratios for reactions
        type(int_array_c) :: mix_conc_indices                            !< Indices of target waters mixing with each target
        type(int_array_c) :: mix_react_indices                        !< Indices of domain target waters only
        real(kind=8), allocatable :: mixing_ratios_mat_conc_mesh(:,:)          !< WMA mesh mixing ratios matrix (implicit Euler)
        real(kind=8), allocatable :: mixing_ratios_mat_conc_bd(:,:)           !< WMA boundary mixing ratios matrix (implicit Euler)
        real(kind=8), allocatable :: bd_mat_prev(:)          !< boundary matrix at previous time step (for time-dependent BCs)
        logical :: Lagr_flag = .false. !< Flag to indicate if Lagrangian tracking is enabled (for debugging and analysis)
    contains
    !> Set procedures
        procedure :: set_time_discr                                   !< Set time discretization
        procedure :: set_mix_react_indices                                 !< [COMMENTED] Set characteristic parameters
        procedure :: set_Lagr_flag
        procedure :: is_initialized                                   !< Check if object is fully initialized
    !> Allocate procedures
        procedure :: allocate_arrays_PDE=>allocate_arrays_PDE_trans  !< Allocate all transient arrays
        procedure :: allocate_F_mat                                   !< Allocate storage matrix
        procedure :: allocate_Y_mat                                   !< Allocate explicit recharge matrix
        procedure :: allocate_Z_mat                                   !< Allocate explicit boundary matrix
        !procedure :: allocate_A_mat                                   !< Allocate implicit domain matrix
        procedure :: allocate_A_mat_inv                               !< Allocate inverse of A_mat
        procedure :: allocate_f_vec                                   !< Allocate independent term vector
        procedure :: allocate_mixing_ratios                           !< Allocate WMA mixing ratios
        procedure :: allocate_mixing_ratios_mat_conc_mesh              !< Allocate WMA matrix
        procedure(compute_mixing_ratios_Delta_t_homog), public, deferred :: compute_mixing_ratios_Delta_t_homog !< Allocate WMA domain matrix
        procedure(allocate_A_mat), public, deferred :: allocate_A_mat !< Allocate WMA domain matrix
        procedure(allocate_X_mat), public, deferred :: allocate_X_mat !< Allocate WMA domain matrix
        procedure :: allocate_mixing_ratios_mat_conc_bd               !< Allocate WMA boundary matrix
        procedure :: allocate_mix_conc_indices                   !< Allocate WMA water indices
    !> Computation procedures
        procedure(compute_F_mat_PDE), public, deferred :: compute_F_mat_PDE   !< Compute storage matrix (deferred)
        !procedure :: compute_E_mat                                    !< Compute external matrix
        !procedure :: compute_X_mat                                    !< Compute explicit domain matrix
        procedure :: compute_Y_mat                                    !< Compute explicit recharge matrix
        procedure :: compute_Z_mat                                    !< Compute explicit boundary matrix
        !procedure :: compute_A_mat                                    !< Compute implicit domain matrix
        !procedure :: compute_lumped_A_mat                             !< Compute mass-lumped A matrix
        procedure :: compute_mix_ratios_R_opt4                       !< Compute mixing ratios option 4
        procedure :: compute_mix_ratio_R_opt4
        procedure :: compute_f_vec                                    !< Compute independent term
        !procedure :: compute_b_vec                           !< Compute RHS of linear system
        !procedure :: compute_A_mat_ODE                                !< Compute ODE matrix
        procedure :: compute_b_ODE                                    !< Compute ODE RHS
        !procedure :: solve_PDE_EE_Delta_t_homog                      !< [COMMENTED] Explicit Euler homog
        !procedure :: solve_PDE_EE_Delta_t_heterog                    !< [COMMENTED] Explicit Euler heterog
        !procedure :: solve_PDE_EI_Delta_t_homog                      !< [COMMENTED] Implicit Euler homog
        !procedure :: solve_PDE_RKF45                                 !< [COMMENTED] RKF45 adaptive
        procedure :: compute_k_RKF45                                  !< Compute RKF45 stages
    end type

    type, public, abstract, extends(PDE_transient_c) :: PDE_1D_transient_c
        !< Inherits all from PDE_transient_c and PDE_1D_c        
        type(BCs_1D_c) :: BCs                                                    !< Boundary conditions
        type(tridiag_matrix_c) :: trans_mat                                   !< Transition matrix T (tridiagonal)
        type(tridiag_matrix_c) :: trans_mat_prev                                       !< Transition matrix T at previous time step (tridiagonal)
        type(tridiag_matrix_c) :: X_mat                                       !< Explicit domain mixing ratios (previous time)
        type(tridiag_matrix_c) :: A_mat                                       !< Implicit domain mixing ratios (current time)
    contains
        procedure :: set_BCs_1D_trans                                    !< Set 1D boundary conditions
        procedure :: compute_X_mat_1D                                    !< Compute explicit domain matrix
        procedure :: compute_A_mat_1D   !< Compute source/sink term
        procedure :: compute_A_mat_ODE_1D
        procedure :: compute_E_mat_1D
        procedure :: allocate_X_mat=>allocate_X_mat_1D                                   !< Allocate explicit domain matrix
        procedure :: allocate_A_mat=>allocate_A_mat_1D                                   !< Allocate implicit domain matrix
        procedure :: allocate_trans_mat=>allocate_trans_mat_1D_trans                              !< Allocate transition matrix
        procedure :: allocate_bd_mat=>allocate_bd_mat_1D_trans                                  !< Allocate boundary matrix
        procedure :: compute_b_vec_1D                           !< Compute RHS of linear system
    end type

    type, public, abstract, extends(PDE_transient_c) :: PDE_2D_transient_c
        !< Inherits all from PDE_transient_c and PDE_2D_c        
        type(BCs_2D_c) :: BCs                                                    !< Boundary conditions
        type(pentadiag_matrix_c) :: trans_mat                                 !< Transition matrix T (pentadiagonal)
        type(pentadiag_matrix_c) :: trans_mat_prev                             !< Transition matrix T at previous time step (pentadiagonal)
        type(pentadiag_matrix_c) :: X_mat                                       !< Explicit domain mixing ratios (previous time)
        type(pentadiag_matrix_c) :: A_mat                                       !< Implicit domain mixing ratios (current time)
    contains
        procedure :: set_BCs_2D_trans
        !procedure :: solve_PDE=>solve_PDE_2D_stat                             !< Solve steady-state PDE
    !> Allocate procedures
        procedure :: allocate_trans_mat=>allocate_trans_mat_2D_trans                              !< Allocate transition matrix
        procedure :: allocate_bd_mat=>allocate_bd_mat_2D_trans                                  !< Allocate boundary matrix
        procedure :: allocate_Z_mat=>allocate_Z_mat_2D                                         !< Allocate explicit boundary matrix (2D)
        procedure :: compute_E_mat_2D
        procedure :: compute_X_mat_2D
        procedure :: compute_A_mat_2D
        procedure :: allocate_X_mat=>allocate_X_mat_2D                                   !< Allocate explicit domain matrix
        procedure :: allocate_A_mat=>allocate_A_mat_2D                                   !< Allocate implicit domain matrix
    end type
!*****************************************************************************************************************************
    !> \brief Abstract interface to compute storage matrix F
    !> \details
!>   Deferred procedure that must be implemented by concrete subclasses.
!>   
!>   The storage matrix F represents the accumulation term in:
!>   \f[ F \frac{dc}{dt} = Tc + Ec_{ext} + Bc_{bd} \f]
!>   
!>   For example:
!>   - Diffusion: F = diag(porosity * retardation)
!>   - Transport: F = diag(porosity)
!>   - MRMT: F accounts for mobile/immobile zones
!>   
!>   Implementation must populate this%F_mat
    abstract interface
        !> \brief Compute storage matrix F
        !> \param[in,out] this The PDE object (modified: F_mat is set)
        subroutine compute_F_mat_PDE(this)
            import PDE_transient_c
            implicit none
            class(PDE_transient_c) :: this
        end subroutine

        subroutine compute_mixing_ratios_Delta_t_homog(this)
           import PDE_transient_c
           implicit none
           class(PDE_transient_c) :: this
           !real(kind=8), intent(in) :: theta
           !type(diag_matrix_c), intent(out), optional :: A_mat_lumped
        end subroutine
        
        subroutine allocate_X_mat(this)
            import PDE_transient_c
            implicit none
            class(PDE_transient_c) :: this
        end subroutine
        
        subroutine allocate_A_mat(this)
            import PDE_transient_c
            implicit none
            class(PDE_transient_c) :: this
        end subroutine
        
    end interface
!*****************************************************************************************************************************
    !> \brief Concrete interfaces for transient PDE methods
    !> \details
    !>   These interfaces define specific computational procedures for:
    !>   - External source/sink matrix E
    !>   - θ-method matrices (A, X, Y, Z)
    !>   - Mass-lumped approximations
    !>   - ODE formulations
    !>   - RKF45 time stepping
    interface
        
        !> \brief Compute external source/sink matrix E
        !> \details
        !>   Computes E for external forcings (recharge, discharge, evaporation)
        !>   E represents spatially distributed sources/sinks
        !> \param[in,out] this The PDE object
        !> \param[out] E_mat External matrix (tridiagonal)
        !> \param[in] k Optional time step index (for time-dependent sources)
        subroutine compute_E_mat_1D(this,Delta_t,E_mat,E_mat_prev)
            import PDE_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(PDE_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: Delta_t
            type(tridiag_matrix_c), intent(out) :: E_mat
            type(tridiag_matrix_c), intent(out) :: E_mat_prev
            !integer(kind=4), intent(in), optional :: k
        end subroutine
        
        !> \brief Compute external source/sink matrix E
        !> \details
        !>   Computes E for external forcings (recharge, discharge, evaporation)
        !>   E represents spatially distributed sources/sinks
        !> \param[in,out] this The PDE object
        !> \param[out] E_mat External matrix (tridiagonal)
        !> \param[in] k Optional time step index (for time-dependent sources)
        subroutine compute_E_mat_2D(this,Delta_t,E_mat,E_mat_prev)
            import PDE_2D_transient_c
            import pentadiag_matrix_c
            implicit none
            class(PDE_2D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: Delta_t
            type(pentadiag_matrix_c), intent(out) :: E_mat
            type(pentadiag_matrix_c), intent(out) :: E_mat_prev
            !integer(kind=4), intent(in), optional :: k
        end subroutine
        
        
        
        !> \brief Compute mass-lumped approximation of A matrix
        !> \details
        !>   Mass lumping diagonalizes A for efficient explicit time stepping:
        !>   - Row-sum lumping: diagonal = sum of each row
        !>   - Preserves total mass
        !>   - Avoids matrix inversion
        !>   - Conditionally stable (Courant condition)
        !> \param[in] this The PDE object (const)
        !> \param[out] A_mat_lumped Diagonal approximation of A (must be allocated)
        ! subroutine compute_lumped_A_mat(this,A_mat_lumped)
        !     import PDE_1D_transient_c
        !     import diag_matrix_c
        !     implicit none
        !     class(PDE_1D_transient_c), intent(in) :: this
        !     type(diag_matrix_c), intent(out) :: A_mat_lumped !> must be allocated
        ! end subroutine
        
        
        
        !> \brief Compute independent term f in linear system
        !> \details
        !>   Computes f = Y*c_rech + Z*c_bd
        !>   where:
        !>   - Y: explicit recharge mixing matrix (diagonal)
        !>   - Z: explicit boundary mixing (size 2)
        !>   - Z: explicit boundary mixing (size 2)
        !>   - c_rech: recharge concentration
        !>   - c_bd: boundary concentrations (2 values)
        !>   
        !>   Result stored in this%f_vec
        !> \param[in,out] this The PDE object (modified: f_vec is set)
        !> \param[in] Delta_t Time step size
        subroutine compute_f_vec(this,Delta_t)
            import PDE_transient_c
            implicit none
            class(PDE_transient_c) :: this
            real(kind=8), intent(in) :: Delta_t
        end subroutine
        
        !> \brief Compute RHS vector b for linear system
        !> \details
        !>   Assembles the right-hand side for θ-method:
        !>   \f[ b = X c^k + f \f]
        !>   where:
        !>   - X = F + (1-θ)Δt(T+E): explicit mixing matrix
        !>   - c^k: concentration at previous time
        !>   - f: independent term (recharge + boundary)
        !> \param[in] this The PDE object (const)
        !> \param[in] theta Time integration parameter
        !> \param[in] conc_old Concentration vector at previous time
        !> \param[in,out] b RHS vector to be computed
        !> \param[in] k Optional time step index
        subroutine compute_b_vec_1D(this,theta,conc_old,b_vec,k)
            import PDE_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(PDE_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(out) :: b_vec(:)
            integer(kind=4), intent(in), optional :: k
        end subroutine
        
        !> \brief Compute ODE system matrix for RKF45
        !> \details
        !>   Transforms PDE into ODE form: dc/dt = A_ODE * c + b_ODE
        !>   Used by Runge-Kutta-Fehlberg 4-5 adaptive time stepping
        !>   \f[ A_{ODE} = F^{-1}(T + E) \f]
        !> \param[in] this The PDE object (const)
        !> \param[out] A_mat ODE system matrix (tridiagonal)
        subroutine compute_A_mat_ODE_1D(this,A_mat)
            import PDE_1D_transient_c
            import tridiag_matrix_c
            implicit none
            class(PDE_1D_transient_c), intent(in) :: this
            type(tridiag_matrix_c), intent(out) :: A_mat
        end subroutine
        
        
        !> \brief Compute ODE RHS vector for RKF45
        !> \details
        !>   Computes b_ODE = F^{-1} * f for ODE formulation
        !>   where f contains recharge and boundary contributions
        !> \param[in] this The PDE object (const)
        !> \return b ODE RHS vector
        function compute_b_ODE(this) result(b)
            import PDE_transient_c
            implicit none
            class(PDE_transient_c), intent(in) :: this
            real(kind=8), allocatable :: b(:)
        end function
        
        !> \brief Solve PDE with Explicit Euler (homogeneous Δt)
        !> \details
        !>   Explicit time integration: c^{k+1} = A_inv * X * c^k + A_inv * f
        !>   - Fast (no linear solve)
        !>   - Conditionally stable (requires small Δt)
        !>   - Homogeneous time steps
        !> \param[in,out] this The PDE object
        !> \param[in] Time_out Output times for solution
        !> \param[out] output Solution matrix (time x space)
        subroutine solve_PDE_EE_Delta_t_homog(this,Time_out,output)
            import PDE_transient_c
            class(PDE_transient_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        !> \brief Solve PDE with Explicit Euler (heterogeneous Δt)
        !> \details
        !>   Like solve_PDE_EE_Delta_t_homog but with adaptive time steps
        !>   Allows different Δt for each time interval
        !> \param[in,out] this The PDE object
        !> \param[in] Time_out Output times for solution
        !> \param[out] output Solution matrix (time x space)
        subroutine solve_PDE_EE_Delta_t_heterog(this,Time_out,output)
            import PDE_transient_c
            class(PDE_transient_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        !> \brief Solve PDE with θ-method (homogeneous Δt)
        !> \details
        !>   General θ-method: A*c^{k+1} = X*c^k + f
        !>   - θ=0: Explicit Euler (fast, conditionally stable)
        !>   - θ=1: Implicit Euler (slow, unconditionally stable)
        !>   - θ=0.5: Crank-Nicolson (2nd order accurate)
        !>   Requires linear solve at each time step
        !> \param[in,out] this The PDE object
        !> \param[in] theta Time integration parameter (0 to 1)
        !> \param[in] Time_out Output times for solution
        !> \param[out] output Solution matrix (time x space)
        subroutine solve_PDE_EI_Delta_t_homog(this,theta,Time_out,output)
            import PDE_transient_c
            class(PDE_transient_c) :: this
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        !> \brief Solve PDE with Runge-Kutta-Fehlberg 4-5 adaptive time stepping
        !> \details
        !>   RKF45 features:
        !>   - 4th and 5th order embedded RK pair
        !>   - Adaptive time step based on local error estimate
        !>   - Automatic error control via tolerance
        !>   - Efficient for smooth solutions
        !>   
        !>   Error estimate: ||c_RK5 - c_RK4||
        !>   Time step update: Δt_new = Δt_old * (tol / error)^(1/5)
        !> \param[in,out] this The PDE object
        !> \param[in] Delta_t_init Initial time step size
        !> \param[in] tolerance Error tolerance for adaptive stepping
        subroutine solve_PDE_RKF45(this,Delta_t_init,tolerance)
            import PDE_transient_c
            class(PDE_transient_c) :: this
            real(kind=8), intent(in) :: Delta_t_init
            real(kind=8), intent(in) :: tolerance
        end subroutine
        
        !> \brief Update time step for RKF45 based on error estimate
        !> \details
        !>   Adaptive time stepping formula:
        !>   \f[ \Delta t_{new} = \Delta t_{old} \left(\frac{\epsilon}{||c_{RK5} - c_{RK4}||}\right)^{1/5} \f]
        !>   where ε is the tolerance and the exponent 1/5 comes from 5th order accuracy
        !>   
        !>   Safety factors are typically applied to avoid excessive changes
        !> \param[in] Delta_t_old Current time step
        !> \param[in] tolerance Error tolerance
        !> \param[in] conc_RK4 4th order solution
        !> \param[in] conc_RK5 5th order solution
        !> \param[out] Delta_t_new Updated time step
        subroutine update_time_step_RKF45(Delta_t_old,tolerance,conc_RK4,conc_RK5,Delta_t_new)
            implicit none
            real(kind=8), intent(in) :: Delta_t_old
            real(kind=8), intent(in) :: tolerance
            real(kind=8), intent(in) :: conc_RK4(:)
            real(kind=8), intent(in) :: conc_RK5(:)
            real(kind=8), intent(out) :: Delta_t_new
        end subroutine
        
        !> \brief Compute RKF45 stage derivatives k1-k6
        !> \details
        !>   Computes the 6 Runge-Kutta stages for RKF45:
        !>   - k1 = f(t, c)
        !>   - k2 = f(t + 1/4 Δt, c + 1/4 k1 Δt)
        !>   - k3 = f(t + 3/8 Δt, c + 3/32 k1 Δt + 9/32 k2 Δt)
        !>   - k4 = f(t + 12/13 Δt, c + ...)
        !>   - k5 = f(t + Δt, c + ...)
        !>   - k6 = f(t + 1/2 Δt, c + ...)
        !>   
        !>   RK4 solution: c_RK4 = c + (25/216 k1 + 1408/2565 k3 + 2197/4104 k4 - 1/5 k5) Δt
        !>   RK5 solution: c_RK5 = c + (16/135 k1 + 6656/12825 k3 + 28561/56430 k4 - 9/50 k5 + 2/55 k6) Δt
        !> \param[in] this The PDE object (const)
        !> \param[in] Delta_t Time step size
        !> \param[in] conc_RK4 Current concentration (base for stages)
        !> \return k Matrix of stage derivatives (6 stages x n points)
        function compute_k_RKF45(this,Delta_t,conc_RK4) result(k)
            import PDE_transient_c
            implicit none
            class(PDE_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: conc_RK4(:)
            real(kind=8), allocatable :: k(:,:)
        end function
    end interface
!*****************************************************************************************************************************
    contains
    !******************************************************************************************************************************
        !> \brief Check if transient PDE object is fully initialized
        !> \details
        !>   Returns .true. only if both required pointer fields are associated
        !>   and the spatial discretization has been built (targets defined).
        !>   This guards against polymorphic assignment of an uninitialised object.
        !> \param[in] this The PDE object
        !> \return res .true. if spatial_discr and time_discr are valid
        logical function is_initialized(this)
            implicit none
            class(PDE_transient_c), intent(in) :: this
            is_initialized = .false.
            if (.not. allocated(this%spatial_discr)) return
            if (.not. this%spatial_discr%Num_targets_defined) return
            if (.not. allocated(this%time_discr)) return
            if (.not. allocated(this%F_mat%diag)) return
            is_initialized = .true.
        end function is_initialized

        !> \brief Set time discretization object
        !> \details
        !>   Associates a time discretization object with this PDE.
        !>   Time discretization can be:
        !>   - Homogeneous (constant Δt)
        !>   - Heterogeneous (variable Δt)
        !>   - Adaptive (RKF45)
        !> \param[in,out] this The PDE object
        !> \param[in] time_discr_obj Time discretization object (polymorphic, must remain in scope)
        subroutine set_time_discr(this,time_discr_obj)
            implicit none
            class(PDE_transient_c) :: this
            class(time_discr_c), intent(in) :: time_discr_obj
            if (allocated(this%time_discr)) deallocate(this%time_discr)
            allocate(this%time_discr, source=time_discr_obj)
        end subroutine
        
        subroutine set_BCs_1D_trans(this,BCs_1D)
           implicit none
           class(PDE_1D_transient_c) :: this
            type(BCs_1D_c), intent(in) :: BCs_1D
           this%BCs=BCs_1D
        end subroutine
        
        subroutine allocate_F_mat(this)
            implicit none
            class(PDE_transient_c) :: this
            !print *, "Allocating F_mat of size ", this%spatial_discr%Num_targets
            call this%F_mat%allocate_array(this%spatial_discr%Num_targets)
            call this%F_mat_prev%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
        subroutine allocate_X_mat_1D(this)
            implicit none
            class(PDE_1D_transient_c) :: this
            call this%X_mat%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
        subroutine allocate_X_mat_2D(this)
            implicit none
            class(PDE_2D_transient_c) :: this
            integer(kind=4) :: nx, ny, n_sub, n_sub2
            call this%X_mat%allocate_array(this%spatial_discr%Num_targets)
            select type (mesh=>this%spatial_discr)
            type is (mesh_2D_Euler_homog_c)
                nx     = mesh%Num_cells_x
                ny     = mesh%Num_cells_y
                n_sub  = ny * (nx - 1)
                n_sub2 = (ny - 1) * nx
                if (allocated(this%X_mat%sub))    deallocate(this%X_mat%sub)
                if (allocated(this%X_mat%super))  deallocate(this%X_mat%super)
                if (allocated(this%X_mat%sub2))   deallocate(this%X_mat%sub2)
                if (allocated(this%X_mat%super2)) deallocate(this%X_mat%super2)
                allocate(this%X_mat%sub(n_sub),    this%X_mat%super(n_sub))
                allocate(this%X_mat%sub2(n_sub2),  this%X_mat%super2(n_sub2))
                this%X_mat%sub    = 0.0d0
                this%X_mat%super  = 0.0d0
                this%X_mat%sub2   = 0.0d0
                this%X_mat%super2 = 0.0d0
            end select
        end subroutine
        
        subroutine allocate_Y_mat(this)
            implicit none
            class(PDE_transient_c) :: this
            call this%Y_mat%allocate_array(this%spatial_discr%Num_targets)
            call this%Y_mat_prev%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
        subroutine allocate_Z_mat(this)
            implicit none
            class(PDE_transient_c) :: this
            allocate(this%Z_mat(2))
            allocate(this%Z_mat_prev(2))
        end subroutine
        
        subroutine allocate_A_mat_1D(this)
            implicit none
            class(PDE_1D_transient_c) :: this
            call this%A_mat%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
        subroutine allocate_A_mat_2D(this)
            implicit none
            class(PDE_2D_transient_c) :: this
            integer(kind=4) :: nx, ny, n_sub, n_sub2
            call this%A_mat%allocate_array(this%spatial_discr%Num_targets)
            select type (mesh=>this%spatial_discr)
            type is (mesh_2D_Euler_homog_c)
                nx     = mesh%Num_cells_x
                ny     = mesh%Num_cells_y
                n_sub  = ny * (nx - 1)
                n_sub2 = (ny - 1) * nx
                if (allocated(this%A_mat%sub))    deallocate(this%A_mat%sub)
                if (allocated(this%A_mat%super))  deallocate(this%A_mat%super)
                if (allocated(this%A_mat%sub2))   deallocate(this%A_mat%sub2)
                if (allocated(this%A_mat%super2)) deallocate(this%A_mat%super2)
                allocate(this%A_mat%sub(n_sub),    this%A_mat%super(n_sub))
                allocate(this%A_mat%sub2(n_sub2),  this%A_mat%super2(n_sub2))
                this%A_mat%sub    = 0.0d0
                this%A_mat%super  = 0.0d0
                this%A_mat%sub2   = 0.0d0
                this%A_mat%super2 = 0.0d0
            end select
        end subroutine
        
        subroutine allocate_f_vec(this)
            implicit none
            class(PDE_transient_c) :: this
            allocate(this%f_vec(this%spatial_discr%Num_targets))
        end subroutine
        
        subroutine allocate_mixing_ratios(this)
            implicit none
            class(PDE_transient_c) :: this
            call this%mixing_ratios_conc%allocate_array(this%spatial_discr%Num_targets)
            call this%mixing_ratios_R%allocate_array(this%spatial_discr%Num_targets)
            call this%mixing_ratios_R_init%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
        !subroutine allocate_mixing_ratios_mat_conc_dom(this)
        !    implicit none
        !    class(PDE_1D_transient_c) :: this
        !    allocate(this%mixing_ratios_mat_conc_dom(this%spatial_discr%Num_targets-2*this%spatial_discr%targets_flag,this%spatial_discr%Num_targets-2*this%spatial_discr%targets_flag))
        !    !call this%mixing_ratios_R%allocate_array(this%spatial_discr%Num_targets)
        !    !call this%mixing_ratios_R_init%allocate_array(this%spatial_discr%Num_targets)
        !end subroutine
        
        subroutine allocate_mixing_ratios_mat_conc_mesh(this)
            implicit none
            class(PDE_transient_c) :: this
            allocate(this%mixing_ratios_mat_conc_mesh(this%spatial_discr%Num_targets,this%spatial_discr%Num_targets))
            this%mixing_ratios_mat_conc_mesh=0d0
            !call this%mixing_ratios_R%allocate_array(this%spatial_discr%Num_targets)
            !call this%mixing_ratios_R_init%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
        subroutine allocate_mixing_ratios_mat_conc_bd(this)
            implicit none
            class(PDE_transient_c) :: this
            allocate(this%mixing_ratios_mat_conc_bd(this%spatial_discr%Num_targets,2))
            this%mixing_ratios_mat_conc_bd=0d0
            !call this%mixing_ratios_R%allocate_array(this%spatial_discr%Num_targets)
            !call this%mixing_ratios_R_init%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
         subroutine allocate_A_mat_inv(this)
            implicit none
            class(PDE_transient_c) :: this
            allocate(this%A_mat_inv(this%spatial_discr%Num_targets,this%spatial_discr%Num_targets))
            !call this%mixing_ratios_conc%allocate_array(this%spatial_discr%Num_targets)
            !call this%mixing_ratios_R%allocate_array(this%spatial_discr%Num_targets)
            !call this%mixing_ratios_R_init%allocate_array(this%spatial_discr%Num_targets)
        end subroutine
        
        subroutine allocate_arrays_PDE_trans(this)
            implicit none
            class(PDE_transient_c) :: this
            
            ! Allocate parent class arrays (same as allocate_arrays_PDE_1D_stat)
            call this%allocate_trans_mat()
            call this%allocate_rech_mat()
            call this%allocate_bd_mat()
            !call this%allocate_bd_mat_prev()
            call this%allocate_source_term_PDE()
            
            ! Allocate transient-specific arrays
            call this%allocate_F_mat()
            !call this%allocate_X_mat()
            call this%allocate_Y_mat()
            call this%allocate_Z_mat()
            !call this%allocate_A_mat()
            call this%allocate_f_vec()

            select type (this)
            class is (PDE_transient_c)
                call this%allocate_X_mat()
                call this%allocate_A_mat()
            end select
        end subroutine

        
        subroutine allocate_mix_conc_indices(this)
            implicit none
            class(PDE_transient_c) :: this
            integer(kind=4) :: i
            call this%mix_conc_indices%allocate_array(this%mixing_ratios_conc%num_cols)
            !call this%mix_react_indices%allocate_array(this%mixing_ratios_conc%num_cols-&
            !    2*this%spatial_discr%targets_flag)
            do i=1,this%mix_conc_indices%num_cols
                call this%mix_conc_indices%cols(i)%allocate_vector(this%mixing_ratios_conc%cols(i)%dim+2) !> +2 for number of upstream and downstream waters, respectively
                !call this%mix_react_indices%cols(i)%allocate_vector(this%mix_react_indices%num_cols+2)
            end do
            ! do i=1,this%mix_react_indices%num_cols
            !     call this%mix_react_indices%cols(i)%allocate_vector(this%mix_react_indices%num_cols+2)
            ! end do
        end subroutine

        subroutine compute_X_mat_1D(this,theta,E_mat)
            !import PDE_1D_transient_c
            !import tridiag_matrix_c
            implicit none
            class(PDE_1D_transient_c) :: this
            real(kind=8), intent(in) :: theta
            class(tridiag_matrix_c), intent(in) :: E_mat
            
            integer(kind=4) :: n
            real(kind=8) :: B_norm_inf,B_norm_1
        
            !call this%X_mat%allocate_array(n)
            this%X_mat%sub=(1d0-theta)*E_mat%sub
            this%X_mat%diag=1d0+(1d0-theta)*E_mat%diag
            this%X_mat%super=(1d0-theta)*E_mat%super
            !print *, this%X_mat%sub
            !print *, this%X_mat%diag
            !print *, this%X_mat%super
        end subroutine

        subroutine compute_X_mat_2D(this,theta,E_mat)
            !import PDE_1D_transient_c
            !import pentadiag_matrix_c
            implicit none
            class(PDE_2D_transient_c) :: this
            real(kind=8), intent(in) :: theta
            class(pentadiag_matrix_c), intent(in) :: E_mat
            
            integer(kind=4) :: nx, ny, n_sub, n_sub2
            real(kind=8) :: B_norm_inf,B_norm_1
        
            !call this%X_mat%allocate_array(n)
            this%X_mat%sub=(1d0-theta)*E_mat%sub
            this%X_mat%sub2=(1d0-theta)*E_mat%sub2
            this%X_mat%diag=1d0+(1d0-theta)*E_mat%diag
            this%X_mat%super=(1d0-theta)*E_mat%super
            this%X_mat%super2=(1d0-theta)*E_mat%super2

            select type (mesh=>this%spatial_discr)
            type is (mesh_2D_Euler_homog_c)
                nx     = mesh%Num_cells_x
                ny     = mesh%Num_cells_y
                n_sub  = ny * (nx - 1)
                n_sub2 = (ny - 1) * nx
                if (allocated(this%X_mat%sub))    deallocate(this%X_mat%sub)
                if (allocated(this%X_mat%super))  deallocate(this%X_mat%super)
                if (allocated(this%X_mat%sub2))   deallocate(this%X_mat%sub2)
                if (allocated(this%X_mat%super2)) deallocate(this%X_mat%super2)
                allocate(this%X_mat%sub(n_sub),    this%X_mat%super(n_sub))
                allocate(this%X_mat%sub2(n_sub2),  this%X_mat%super2(n_sub2))
                this%X_mat%sub    = 0.0d0
                this%X_mat%super  = 0.0d0
                this%X_mat%sub2   = 0.0d0
                this%X_mat%super2 = 0.0d0
            end select
        end subroutine
        
        subroutine compute_Y_mat(this,k)
            implicit none
            class(PDE_transient_c) :: this
            !real(kind=8), intent(in) :: theta
            integer(kind=4), intent(in), optional :: k
            !class(tridiag_matrix_c), intent(in) :: E_mat
            
            integer(kind=4) :: n
            real(kind=8) :: B_norm_inf,B_norm_1
        
            select type (time_discr=>this%time_discr)
            type is (time_discr_homog_c)
                if (this%dimless) then
                    this%Y_mat%diag=this%rech_mat%diag*time_discr%Delta_t_D
                else
                    this%Y_mat%diag=this%rech_mat%diag*time_discr%Delta_t
                end if
                !E_mat%super=E_mat%super*time_discr%Delta_t
            type is (time_discr_heterog_c)
                if (this%dimless) then
                    this%Y_mat%diag=this%rech_mat%diag*time_discr%Delta_t_D(k)
                else
                    this%Y_mat%diag=this%rech_mat%diag*time_discr%Delta_t(k)
                end if
            end select
            !this%Y_mat_prev%diag=this%Y_mat%diag
            this%Y_mat%diag=this%Y_mat%diag/this%F_mat_prev%diag
            this%Y_mat_prev%diag=this%Y_mat%diag/this%F_mat%diag
            !call this%X_mat%allocate_array(n)
            !this%Y_mat%sub=(1d0-theta)*this%rech_mat%sub
            !this%Y_mat%diag=1d0+(1d0-theta)*this%rech_mat%diag
            !this%Y_mat%super=(1d0-theta)*E_mat%super
        end subroutine
       
        subroutine compute_Z_mat(this,k)
            !import PDE_1D_transient_c
            !import tridiag_matrix_c
            implicit none
            class(PDE_transient_c) :: this
            !real(kind=8), intent(in) :: theta
            integer(kind=4), intent(in), optional :: k
            !class(tridiag_matrix_c), intent(in) :: E_mat
            
            integer(kind=4) :: n
            real(kind=8) :: B_norm_inf,B_norm_1
        
            select type (time_discr=>this%time_discr)
            type is (time_discr_homog_c)
                if (this%dimless) then
                    this%Z_mat=this%bd_mat*time_discr%Delta_t_D
                    this%Z_mat_prev=this%bd_mat_prev*time_discr%Delta_t_D
                else
                    this%Z_mat=this%bd_mat*time_discr%Delta_t
                    this%Z_mat_prev=this%bd_mat_prev*time_discr%Delta_t
                end if
                !this%Z_mat=this%bd_mat*time_discr%Delta_t
            type is (time_discr_heterog_c)
                if (this%dimless) then
                    this%Z_mat=this%bd_mat*time_discr%Delta_t_D(k)
                    this%Z_mat_prev=this%bd_mat_prev*time_discr%Delta_t_D(k)
                else
                    this%Z_mat=this%bd_mat*time_discr%Delta_t(k)
                    this%Z_mat_prev=this%bd_mat_prev*time_discr%Delta_t(k)
                end if
            end select
            this%Z_mat(1)=this%Z_mat(1)/this%F_mat_prev%diag(1)
            this%Z_mat(2)=this%Z_mat(2)/this%F_mat_prev%diag(this%spatial_discr%num_targets)
            this%Z_mat_prev(1)=this%Z_mat_prev(1)/this%F_mat%diag(1)
            this%Z_mat_prev(2)=this%Z_mat_prev(2)/this%F_mat%diag(this%spatial_discr%num_targets)
            !this%Z_mat_prev=this%Z_mat
            !call this%X_mat%allocate_array(n)
            !this%Y_mat%sub=(1d0-theta)*this%rech_mat%sub
            !this%Y_mat%diag=1d0+(1d0-theta)*this%rech_mat%diag
            !this%Y_mat%super=(1d0-theta)*E_mat%super
        end subroutine

        function compute_mix_ratios_R_opt4(this) result(mix_ratios_R)
            implicit none
            class(PDE_transient_c) :: this !> object
            !real(kind=8), intent(in) :: theta !> theta parameter
            real(kind=8), allocatable :: mix_ratios_R(:) !> new mixing ratios
            
            integer(kind=4) :: i,j,num_up
            real(kind=8) :: sum
            
            allocate(mix_ratios_R(this%mixing_ratios_R%num_cols))
            !n_targets=this%spatial_discr%Num_targets
            !> Mixing ratios for reaction amounts in next time step at computation targets
            do i=1,this%mixing_ratios_R%num_cols
                sum=this%mixing_ratios_R%cols(i)%col_1(1) !> initial mixing ratio (first entry)
                !> loop over downstream waters
                do j=1,this%mix_react_indices%cols(i)%col_1(&
                    this%mix_react_indices%cols(i)%dim)
                    num_up=this%mix_react_indices%cols(i)%col_1(&
                        this%mix_react_indices%cols(i)%dim-1) !> number of upstream waters
                    sum=sum+this%mixing_ratios_R%cols(i)%col_1(1+num_up+j) !> we update sum with downstream waters
                end do
                mix_ratios_R(i)=sum !> final mixing ratio for target i
            end do
            !mix_ratios_R=mix_ratios_R*theta !> multiply by theta parameter
        end function

        function compute_mix_ratio_R_opt4(this,i) result(mix_ratio_R)
            implicit none
            class(PDE_transient_c) :: this !> object
            integer(kind=4), intent(in) :: i !> target index
            !real(kind=8), intent(in) :: theta !> theta parameter
            real(kind=8) :: mix_ratio_R !> new mixing ratio
            
            integer(kind=4) :: j,num_up
            real(kind=8) :: sum
            
            num_up=this%mix_react_indices%cols(i)%col_1(&
                this%mix_react_indices%cols(i)%dim-1) !> number of upstream waters
            !allocate(mix_ratio_R(this%mixing_ratios_R%num_cols))
            !n_targets=this%spatial_discr%Num_targets
            !> Mixing ratios for reaction amounts in next time step at computation targets
            !do i=1,this%mixing_ratios_R%num_cols
                sum=this%mixing_ratios_R%cols(i)%col_1(1) !> initial mixing ratio (first entry)
                !> loop over downstream waters
                do j=1,this%mix_react_indices%cols(i)%col_1(&
                    this%mix_react_indices%cols(i)%dim)
                    ! num_up=this%mix_react_indices%cols(i)%col_1(&
                    !     this%mix_react_indices%cols(i)%dim-1) !> number of upstream waters
                    sum=sum+this%mixing_ratios_R%cols(i)%col_1(1+num_up+j) !> we update sum with downstream waters
                end do
                mix_ratio_R=sum !> final mixing ratio for target i
                if (mix_ratio_R < 0d0 .or. mix_ratio_R > 1d0) then
                    error stop "Error in compute_mix_ratio_R_opt4: mixing ratio out of bounds"
                end if
            !end do
            !mix_ratio_R=mix_ratio_R*theta !> multiply with theta parameter
        end function

        subroutine allocate_trans_mat_1D_trans(this)
            implicit none
            class(PDE_1D_transient_c) :: this
            call this%trans_mat%allocate_array(this%spatial_discr%Num_targets)
            call this%trans_mat_prev%allocate_array(this%spatial_discr%Num_targets)
        end subroutine

        subroutine allocate_trans_mat_2D_trans(this)
            implicit none
            class(PDE_2D_transient_c) :: this
            integer(kind=4) :: nx, ny, n_sub, n_sub2
            call this%trans_mat%allocate_array(this%spatial_discr%Num_targets)
            call this%trans_mat_prev%allocate_array(this%spatial_discr%Num_targets)
            select type (mesh=>this%spatial_discr)
            type is (mesh_2D_Euler_homog_c)
                nx     = mesh%Num_cells_x
                ny     = mesh%Num_cells_y
                n_sub  = ny * (nx - 1)
                n_sub2 = (ny - 1) * nx
                !> --- trans_mat ---
                if (allocated(this%trans_mat%sub))    deallocate(this%trans_mat%sub)
                if (allocated(this%trans_mat%super))  deallocate(this%trans_mat%super)
                if (allocated(this%trans_mat%sub2))   deallocate(this%trans_mat%sub2)
                if (allocated(this%trans_mat%super2)) deallocate(this%trans_mat%super2)
                allocate(this%trans_mat%sub(n_sub),    this%trans_mat%super(n_sub))
                allocate(this%trans_mat%sub2(n_sub2),  this%trans_mat%super2(n_sub2))
                this%trans_mat%sub    = 0.0d0
                this%trans_mat%super  = 0.0d0
                this%trans_mat%sub2   = 0.0d0
                this%trans_mat%super2 = 0.0d0
                !> --- trans_mat_prev ---
                if (allocated(this%trans_mat_prev%sub))    deallocate(this%trans_mat_prev%sub)
                if (allocated(this%trans_mat_prev%super))  deallocate(this%trans_mat_prev%super)
                if (allocated(this%trans_mat_prev%sub2))   deallocate(this%trans_mat_prev%sub2)
                if (allocated(this%trans_mat_prev%super2)) deallocate(this%trans_mat_prev%super2)
                allocate(this%trans_mat_prev%sub(n_sub),    this%trans_mat_prev%super(n_sub))
                allocate(this%trans_mat_prev%sub2(n_sub2),  this%trans_mat_prev%super2(n_sub2))
                this%trans_mat_prev%sub    = 0.0d0
                this%trans_mat_prev%super  = 0.0d0
                this%trans_mat_prev%sub2   = 0.0d0
                this%trans_mat_prev%super2 = 0.0d0
            end select
        end subroutine

        subroutine allocate_bd_mat_1D_trans(this)
            implicit none
            class(PDE_1D_transient_c) :: this
            allocate(this%bd_mat(2))
            allocate(this%bd_mat_prev(2))
            this%bd_mat=0d0
            this%bd_mat_prev=this%bd_mat
        end subroutine

        subroutine allocate_bd_mat_2D_trans(this)
            implicit none
            class(PDE_2D_transient_c) :: this
            allocate(this%bd_mat(4))
            allocate(this%bd_mat_prev(4))
            this%bd_mat = 0d0
            this%bd_mat_prev = 0d0
        end subroutine

        subroutine allocate_Z_mat_2D(this)
            implicit none
            class(PDE_2D_transient_c) :: this
            allocate(this%Z_mat(4))
            allocate(this%Z_mat_prev(4))
        end subroutine

        subroutine set_BCs_2D_trans(this,BCs_2D)
            implicit none
            class(PDE_2D_transient_c) :: this
            type(BCs_2D_c), intent(in) :: BCs_2D
           this%BCs=BCs_2D
        end subroutine

        subroutine set_mix_react_indices(this)
            implicit none
            class(PDE_transient_c) :: this
            integer(kind=4) :: i, j, num_up_conc, num_up_domain, num_down_domain, num_up_ext
            
            call this%mix_react_indices%allocate_array(this%mixing_ratios_R%num_cols)
            do i=1,this%mix_react_indices%num_cols
                num_up_conc = this%mix_conc_indices%cols(i)%col_1( &
                    this%mix_conc_indices%cols(i)%dim-1)
                
                if (this%time_discr%int_method > 1) then
                    !> Implicit: all domain targets are coupled
                    num_up_domain = i - 1
                    num_down_domain = this%spatial_discr%Num_targets - i
                else
                    !> Explicit: only actual neighbors, all domain
                    num_up_domain = num_up_conc
                    num_down_domain = this%mix_conc_indices%cols(i)%col_1( &
                        this%mix_conc_indices%cols(i)%dim)
                end if
                num_up_ext = num_up_conc - num_up_domain
                
                call this%mix_react_indices%cols(i)%allocate_vector( &
                    1 + num_up_domain + num_down_domain + 2)
                
                !> Current water index
                this%mix_react_indices%cols(i)%col_1(1) = &
                    this%mix_conc_indices%cols(i)%col_1(1)
                !> Number of upstream and downstream domain waters
                this%mix_react_indices%cols(i)%col_1( &
                    this%mix_react_indices%cols(i)%dim-1) = num_up_domain
                this%mix_react_indices%cols(i)%col_1( &
                    this%mix_react_indices%cols(i)%dim) = num_down_domain
                
                !> Copy upstream domain water indices (skip external upstream entries)
                do j=1, num_up_domain
                    this%mix_react_indices%cols(i)%col_1(1+j) = &
                        this%mix_conc_indices%cols(i)%col_1(1 + num_up_ext + j)
                end do
                !> Copy downstream domain water indices (skip external downstream entries)
                do j=1, num_down_domain
                    this%mix_react_indices%cols(i)%col_1(1 + num_up_domain + j) = &
                        this%mix_conc_indices%cols(i)%col_1(1 + num_up_conc + j)
                end do
            end do
        end subroutine

        !> \brief Set Lagrangian/Eulerian method flag
        !> \param[in,out] this      Transport problem object
        !> \param[in]     Lagr_flag Method selector (.true.=Lagrangian, .false.=Eulerian)
        !> \details 
        !>   Controls which method is used for transport matrix computation.
        !>   
        !>   Lagrangian Method:
        !>   - Follows fluid particles
        !>   - Mass-conservative by construction
        !>   - Better for advection-dominated flows
        !>   
        !>   Eulerian Method:
        !>   - Fixed spatial grid
        !>   - Traditional finite difference/volume
        !>   - Better for diffusion-dominated flows
        subroutine set_Lagr_flag(this,Lagr_flag)
        implicit none
        class(PDE_transient_c) :: this                 !< Transport problem object (modified)
        logical, intent(in) :: Lagr_flag                        !< Method flag (.true.=Lagrangian)
        this%Lagr_flag=Lagr_flag                                !< Store flag
        end subroutine

        !> \brief Compute implicit domain matrix A for θ-method
        !> \details
        !>   Computes the implicit part A based on θ-method discretization:
        !>   \f[ A = F - \theta \Delta t (T + E) \f]
        !>   where:
        !>   - F: storage matrix (accumulation)
        !>   - T: transport/diffusion operator
        !>   - E: external source/sink matrix
        !>   - θ: time integration parameter (0=explicit, 1=implicit, 0.5=Crank-Nicolson)
        !>   
        !>   Result is stored in this%A_mat
        !> \param[in] this The PDE object (const)
        !> \param[in] theta Time integration parameter (0 to 1)
        !> \param[in] E_mat External matrix E (tridiagonal)
        subroutine compute_A_mat_1D(this,theta,E_mat)
            implicit none
            class(PDE_1D_transient_c), intent(inout) :: this
            real(kind=8), intent(in) :: theta
            class(tridiag_matrix_c), intent(in) :: E_mat

            this%A_mat%sub=-theta*E_mat%sub
            this%A_mat%diag=1d0-theta*E_mat%diag
            this%A_mat%super=-theta*E_mat%super
            !print *, this%A_mat%super
        end subroutine

        subroutine compute_A_mat_2D(this,theta,E_mat)
            implicit none
            class(PDE_2D_transient_c), intent(inout) :: this
            real(kind=8), intent(in) :: theta
            class(pentadiag_matrix_c), intent(in) :: E_mat

            this%A_mat%sub=-theta*E_mat%sub
            this%A_mat%sub2=-theta*E_mat%sub2
            this%A_mat%diag=1d0-theta*E_mat%diag
            this%A_mat%super=-theta*E_mat%super
            this%A_mat%super2=-theta*E_mat%super2
        end subroutine
end module