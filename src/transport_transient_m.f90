!> \file transport_transient_m.f90
!> \brief Time-dependent 1D advection-dispersion-reaction transport module
!> \details 
!>   Defines classes and procedures for solving transient 1D transport problems.
!>   Extends diffusion with advection and multiple boundary condition scenarios.
!>
!>   Governing Equation (1D Transient ADE):
!>   \f[
!>     \phi \frac{\partial c}{\partial t} = \frac{\partial}{\partial x}\left(D \frac{\partial c}{\partial x}\right) 
!>     - q\frac{\partial(c)}{\partial x} + r \left( c_r - c \right)
!>   \f]
!>   where:
!>   - φ: Porosity [-]
!>   - c: Concentration [C]
!>   - D: Dispersion coefficient [L²/T]
!>   - q: Darcy flux [L/T]
!>   - r: sink/source term [1/T]
!>
!>   Numerical Methods:
!>   1. **Eulerian Approach**: Fixed spatial grid, track concentration evolution
!>   2. **Lagrangian Approach**: Follow fluid particles, track mass transport
!>
!>   Time Integration Schemes:
!>   - Explicit Euler (EE): Simple, conditionally stable
!>   - Implicit Euler (EI): Unconditionally stable, requires linear solver
!>   - Crank-Nicolson: Second-order accurate
!>
!>   Boundary Conditions:
!>   - Dirichlet: Fixed concentration c = c₀
!>   - Neumann homogeneous: zero dispersive flux
!>   - Robin: prescribed mass flux -D∂c/∂x + qc = F
!>
!>   Applications:
!>   - Recharge/discharge scenarios
!>   - Evaporation problems
!>   - Reactive transport (coupled with chemistry)
!>   - Contaminant migration
!>
!>   Class Hierarchy:
!>   - transport_1D_transient_c extends diffusion_1D_transient_c
!>   - Adds advection physics and transport-specific methods
!>
!>   \see diffusion_transient_m, transport_properties_heterog_m, transport_stab_params_m
module transport_transient_m
    use diffusion_transient_m, only: diffusion_1D_transient_c, fund_sol_diff_eqn_1D, &
            diffusion_2D_transient_c
    use transport_properties_heterog_m, only: tpt_props_heterog_1D_c, &
            tpt_props_heterog_2D_c
    use transport_stab_params_m, only: stab_params_tpt_1D_c, &
            stab_params_tpt_2D_c
    use arrays_m, only: diag_matrix_c
    use time_discr_m, only: time_discr_homog_c
    implicit none
    save
    private
    !> \class transport_1D_transient_c
    !> \brief 1D time-dependent transport problem class
    !> \extends diffusion_1D_transient_c
    !> \details 
    !>   Encapsulates all data and methods for solving 1D transient ADE problems.
    !>   Extends diffusion with advection terms and transport-specific properties.
    type, public, extends(diffusion_1D_transient_c) :: transport_1D_transient_c
        !> \var tpt_props_heterog
        !> Transport properties object (flux, dispersion, porosity, source term)
        type(tpt_props_heterog_1D_c) :: tpt_props_heterog
        
        !> \var stab_params_tpt
        !> Stability parameters (Courant, Peclet, Fourier numbers)
        type(stab_params_tpt_1D_c) :: stab_params_tpt
        
        
    contains
        
        !> Set concentration reaction flags from source term
        procedure :: set_conc_r_flag=>set_conc_r_flag_tpt
        !> Compute flux matrix F for transport equation
        procedure :: compute_F_mat_PDE=>compute_F_mat_tpt
        !> Compute transient transport matrix (dispatches to Lagr/Euler)
        procedure :: compute_trans_mat_PDE=>compute_trans_mat_tpt_trans
        !> Compute transport matrix using Lagrangian method
        procedure :: compute_trans_mat_tpt_trans_Lagr
        !> Compute transport matrix using Eulerian method
        procedure :: compute_trans_mat_tpt_trans_Euler
        !> Set stability parameters object
        procedure :: set_stab_params_tpt
        !> Initialize transport problem from input files
        !procedure :: initialise_PDE=>initialise_transport_1D_transient
        !> Initialize for reactive transport (RT) coupling
        procedure :: initialise_PDE=>initialise_transport_1D_transient_RT
        !> Write transport solution to output files
        procedure :: write_PDE=>write_transport_1D_transient
        !> Set transport properties object
        procedure :: set_tpt_props_heterog_obj
        !> Compute mass balance error (PMF BC, evaporation)
        procedure :: mass_balance_error_ADE_trans_PMF_evap
        !> Compute mass balance error (Dirichlet BC, evaporation)
        procedure :: mass_balance_error_ADE_trans_Dirichlet_evap
        !> Compute mass balance error (PMF BC, recharge)
        procedure :: mass_balance_error_ADE_trans_PMF_recharge
        !> Compute mass balance error (Dirichlet BC, recharge)
        procedure :: mass_balance_error_ADE_trans_Dirichlet_recharge
        !> Compute mass balance error (PMF BC, discharge)
        procedure :: mass_balance_error_ADE_trans_PMF_discharge
        !> Compute mass balance error (Dirichlet BC, discharge)
        procedure :: mass_balance_error_ADE_trans_Dirichlet_discharge
        !> Check time step against stability criteria
        procedure :: check_Delta_t
        !> Read transport data for WMA (Water Mixing Analysis)
        procedure :: read_transport_data_WMA
        !> Write transport data for WMA
        procedure :: write_transport_data_WMA=>write_transport_data_WMA_1D
        !> Fundamental solution of 1D transport equation
        procedure :: fund_sol_tpt_eqn_1D
        !> Compute mixing ratios for homogeneous time step
        procedure :: compute_mixing_ratios_Delta_t_homog=>compute_mixing_ratios_Delta_t_homog_1D
        !> Update mixing ratios for homogeneous time step
        procedure :: update_mixing_ratios_Delta_t_homog
        !> Solve transport using explicit Euler with homogeneous Δt
        procedure :: solve_tpt_EE_Delta_t_homog
    end type
    
    type, public, extends(diffusion_2D_transient_c) :: transport_2D_transient_c
        type(tpt_props_heterog_2D_c) :: tpt_props_heterog
        type(stab_params_tpt_2D_c) :: stab_params_tpt
    contains
        !> Dispatch: selects Lagrangian or Eulerian based on Lagr_flag
        procedure :: compute_trans_mat_PDE=>compute_trans_mat_tpt_trans_2D_dispatch
        !> Compute 2D transport matrix using Lagrangian method (dispersion only)
        procedure :: compute_trans_mat_tpt_trans_Lagr_2D
        !> Compute 2D transport matrix using Eulerian method (advection + dispersion)
        procedure :: compute_trans_mat_tpt_trans_Euler_2D=>compute_trans_mat_tpt_trans_2D
        procedure :: compute_F_mat_PDE=>compute_F_mat_tpt_trans_2D
        procedure :: write_transport_data_WMA=>write_transport_data_WMA_2D
        procedure :: compute_mixing_ratios_Delta_t_homog=>compute_mixing_ratios_Delta_t_homog_2D
        procedure :: initialise_PDE=>initialise_transport_2D_transient_RT
        procedure :: set_tpt_props_heterog_obj_2D
        procedure :: set_stab_params_tpt_2D
    end type
    !> \brief Interface declarations for external subroutines
    !> \details 
    !>   Defines abstract interfaces for procedures implemented in separate files.
    !>   Allows polymorphic dispatch and modular code organization.
    interface
        !> \brief Compute flux matrix F for transport equation
        !> \param[in,out] this Transport problem object
        !> \details Computes F matrix including advection terms
        !subroutine compute_F_mat_tpt(this)
        !    import transport_1D_transient_c
        !    implicit none
        !    class(transport_1D_transient_c) :: this
        !end subroutine
        
        !> \brief Compute transport matrix using Lagrangian method
        !> \param[in,out] this Transport problem object
        !> \details Lagrangian approach: follow fluid particles
        subroutine compute_trans_mat_tpt_trans_Lagr(this)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c) :: this
        end subroutine
        
        !> \brief Compute transport matrix using Eulerian method
        !> \param[in,out] this Transport problem object
        !> \details Eulerian approach: fixed spatial grid
        subroutine compute_trans_mat_tpt_trans_Euler(this)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c) :: this
        end subroutine
        
        !> \brief Initialize 1D transient transport problem
        !> \param[in,out] this Transport problem object
        !> \param[in]     root Root name for input files
        !> \details Reads mesh, BCs, time discretization, transport properties
        subroutine initialise_transport_1D_transient(this,path,root,mesh_type)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c) :: this
            character(len=*), intent(in) :: path
            character(len=*), intent(in) :: root
            integer(kind=4), intent(in) :: mesh_type
        end subroutine

        subroutine initialise_transport_2D_transient_RT(this,path,root,mesh_type)
            import transport_2D_transient_c
            implicit none
            class(transport_2D_transient_c) :: this
            character(len=*), intent(in) :: path
            character(len=*), intent(in) :: root
            integer(kind=4), intent(in) :: mesh_type
        end subroutine
        
        !> \brief Initialize transport for reactive transport coupling
        !> \param[in,out] this      Transport problem object
        !> \param[in]     path      Directory path for input files
        !> \param[in]     root      Root name for input files
        !> \param[in]     mesh_type Mesh type (1=1D homog, 2=1D heterog, 3=radial)
        !> \details Sets up transport problem for coupling with geochemistry
        subroutine initialise_transport_1D_transient_RT(this,path,root,mesh_type)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c) :: this
            character(len=*), intent(in) :: path
            character(len=*), intent(in) :: root
            integer(kind=4), intent(in) :: mesh_type
        end subroutine
        
        !> \brief Write transport solution to output files
        !> \param[in] this     Transport problem object
        !> \param[in] root     Root name for output files
        !> \param[in] Time_out Output time values [T]
        !> \param[in] output   Concentration vs time output [C]
        !> \details Writes concentration profiles at specified times
        subroutine write_transport_1D_transient(this)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
        !     character(len=*), intent(in) :: root
        !     real(kind=8), intent(in) :: Time_out(:)
        !     real(kind=8), intent(in) :: output(:,:)
        end subroutine
        
        !> \brief Mass balance error for Dirichlet BC with recharge
        !> \param[in] this     Transport problem object
        !> \param[in] conc_old Concentration at previous time [C]
        !> \param[in] conc_new Concentration at current time [C]
        !> \param[in] Delta_t  Time step [T]
        !> \param[in] Delta_x  Cell size [L]
        !> \return mass_bal_err Mass balance error [-]
        function mass_balance_error_ADE_trans_Dirichlet_recharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        !> \brief Mass balance error for Dirichlet BC with discharge
        !> \param[in] this     Transport problem object
        !> \param[in] conc_old Concentration at previous time [C]
        !> \param[in] conc_new Concentration at current time [C]
        !> \param[in] Delta_t  Time step [T]
        !> \param[in] Delta_x  Cell size [L]
        !> \return mass_bal_err Mass balance error [-]
        function mass_balance_error_ADE_trans_Dirichlet_discharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        !> \brief Mass balance error for Dirichlet BC with evaporation
        !> \param[in] this     Transport problem object
        !> \param[in] conc_old Concentration at previous time [C]
        !> \param[in] conc_new Concentration at current time [C]
        !> \param[in] Delta_t  Time step [T]
        !> \param[in] Delta_x  Cell size [L]
        !> \return mass_bal_err Mass balance error [-]
        function mass_balance_error_ADE_trans_Dirichlet_evap(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        !> \brief Mass balance error for PMF (Prescribed Mass Flux) BC with recharge
        !> \param[in] this     Transport problem object
        !> \param[in] conc_old Concentration at previous time [C]
        !> \param[in] conc_new Concentration at current time [C]
        !> \param[in] Delta_t  Time step [T]
        !> \param[in] Delta_x  Cell size [L]
        !> \return mass_bal_err Mass balance error [-]
        function mass_balance_error_ADE_trans_PMF_recharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        !> \brief Mass balance error for PMF BC with discharge
        !> \param[in] this     Transport problem object
        !> \param[in] conc_old Concentration at previous time [C]
        !> \param[in] conc_new Concentration at current time [C]
        !> \param[in] Delta_t  Time step [T]
        !> \param[in] Delta_x  Cell size [L]
        !> \return mass_bal_err Mass balance error [-]
        function mass_balance_error_ADE_trans_PMF_discharge(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function
        
        !> \brief Mass balance error for PMF BC with evaporation
        !> \param[in] this     Transport problem object
        !> \param[in] conc_old Concentration at previous time [C]
        !> \param[in] conc_new Concentration at current time [C]
        !> \param[in] Delta_t  Time step [T]
        !> \param[in] Delta_x  Cell size [L]
        !> \return mass_bal_err Mass balance error [-]
        function mass_balance_error_ADE_trans_PMF_evap(this,conc_old,conc_new,Delta_t,Delta_x) result(mass_bal_err)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: conc_new(:)
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: Delta_x
            real(kind=8) :: mass_bal_err
        end function

        !> \brief Read transport data for Water Mixing Approach (WMA)
        !> \param[in,out] this Transport problem object
        !> \param[in]     root Root name for input files
        !> \details Reads transport-specific data for WMA applications
        subroutine read_transport_data_WMA(this,dir,root)
            import transport_1D_transient_c
            class(transport_1D_transient_c) :: this
            character(len=*), intent(in) :: dir
            character(len=*), intent(in) :: root
        end subroutine
        
        !> \brief Read spatial/temporal discretization for WMA
        !> \param[in,out] this Transport problem object
        !> \param[in]     unit File unit number
        !> \param[in]     root Root name for input files
        !> \details Reads mesh and time discretization parameters
        subroutine read_discretisation_WMA(this,unit,root)
            import transport_1D_transient_c
            class(transport_1D_transient_c) :: this
            integer(kind=4), intent(in) :: unit
            character(len=*), intent(in) :: root
        end subroutine
        
        !> \brief Write transport data for WMA to file
        !> \param[in] this Transport problem object
        !> \param[in] unit File unit number for output
        !> \details Writes transport state variables and parameters
        subroutine write_transport_data_WMA_1D(this,path,root)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c), intent(in) :: this
            character(len=*), intent(in) :: path
            character(len=*), intent(in) :: root
        end subroutine

        subroutine write_transport_data_WMA_2D(this,path,root)
            import transport_2D_transient_c
            implicit none
            class(transport_2D_transient_c), intent(in) :: this
            character(len=*), intent(in) :: path
            character(len=*), intent(in) :: root
        end subroutine
        
        !> \brief Compute mixing ratios for homogeneous time step
        !> \param[in,out] this         Transport problem object
        !> \param[out]    A_mat_lumped Optional lumped mass matrix
        !> \details 
        !>   Computes mixing ratios between cells for uniform time stepping.
        !>   Used in Lagrangian and Eulerian methods for mass transfer.
        subroutine compute_mixing_ratios_Delta_t_homog_1D(this)
            import transport_1D_transient_c
            implicit none
            class(transport_1D_transient_c) :: this
            !type(diag_matrix_c), intent(out), optional :: A_mat_lumped
        end subroutine

        subroutine compute_mixing_ratios_Delta_t_homog_2D(this)
            import transport_2D_transient_c
            implicit none
            class(transport_2D_transient_c) :: this
            !type(diag_matrix_c), intent(out), optional :: A_mat_lumped
        end subroutine
        
        !> \brief Solve 1D transient transport using explicit Euler
        !> \param[in,out] this     Transport problem object
        !> \param[in]     Time_out Output time values [T]
        !> \param[out]    output   Concentration vs time output [C]
        !> \details 
        !>   Solves 1D transient ADE with homogeneous time step using explicit Euler.
        !>   Suitable for Lagrangian or Eulerian methods.
        !>   
        !>   Algorithm:
        !>   - Explicit time integration: c^{n+1} = c^n + Δt·f(c^n)
        !>   - Conditional stability (Courant ≤ 1, Fourier ≤ 0.5)
        !>   - Intermediate results written to binary file
        subroutine solve_tpt_EE_Delta_t_homog(this,Time_out,output)
            import transport_1D_transient_c
            class(transport_1D_transient_c) :: this
            real(kind=8), intent(in) :: Time_out(:)
            real(kind=8), intent(out) :: output(:,:)
        end subroutine
        
        !> \brief Update mixing ratios for next time step
        !> \param[in,out] this Transport problem object
        !> \details Updates mixing ratios when properties or time step change
        subroutine update_mixing_ratios_Delta_t_homog(this,Delta_t)
        import transport_1D_transient_c
        class(transport_1D_transient_c) :: this
        real(kind=8), intent(in) :: Delta_t
        end subroutine

        subroutine compute_trans_mat_tpt_trans_2D(this)
            import transport_2D_transient_c
            implicit none
            class(transport_2D_transient_c) :: this
        end subroutine
        
        !> \brief Compute 2D transport matrix using Lagrangian method
        !> \param[in,out] this 2D Transport problem object
        !> \details Lagrangian approach: dispersion only, no advection (symmetric matrix)
        subroutine compute_trans_mat_tpt_trans_Lagr_2D(this)
            import transport_2D_transient_c
            implicit none
            class(transport_2D_transient_c) :: this
        end subroutine
    end interface
    
    contains
      
        !> \brief Set stability parameters object
        !> \param[in,out] this             Transport problem object
        !> \param[in]     stab_params_tpt  Stability parameters (Courant, Peclet, Fourier)
        !> \details 
        !>   Assigns pre-computed stability parameters to transport object.
        !>   Parameters must be computed externally using compute_stab_params_tpt.
        subroutine set_stab_params_tpt(this,stab_params_tpt)
            implicit none
            class(transport_1D_transient_c) :: this             !< Transport problem object (modified)
            class(stab_params_tpt_1D_c), intent(in) :: stab_params_tpt  !< Stability parameters to assign
            this%stab_params_tpt=stab_params_tpt                !< Copy stability parameters
        end subroutine

        !> \brief Set stability parameters object
        !> \param[in,out] this             Transport problem object
        !> \param[in]     stab_params_tpt  Stability parameters (Courant, Peclet, Fourier)
        !> \details 
        !>   Assigns pre-computed stability parameters to transport object.
        !>   Parameters must be computed externally using compute_stab_params_tpt.
        subroutine set_stab_params_tpt_2D(this,stab_params_tpt)
            implicit none
            class(transport_2D_transient_c) :: this             !< Transport problem object (modified)
            class(stab_params_tpt_2D_c), intent(in) :: stab_params_tpt  !< Stability parameters to assign
            this%stab_params_tpt=stab_params_tpt                !< Copy stability parameters
        end subroutine
        
        !> \brief Set transport properties object
        !> \param[in,out] this              Transport problem object
        !> \param[in]     tpt_props_heterog Transport properties (flux, dispersion, source)
        !> \details 
        !>   Assigns transport properties to problem object.
        !>   Properties must be initialized externally (read from file or computed).
        subroutine set_tpt_props_heterog_obj(this,tpt_props_heterog)
            implicit none
            class(transport_1D_transient_c) :: this             !< Transport problem object (modified)
            class(tpt_props_heterog_1D_c), intent(in) :: tpt_props_heterog  !< Transport properties to assign
            this%tpt_props_heterog=tpt_props_heterog            !< Copy transport properties
        end subroutine

        subroutine set_tpt_props_heterog_obj_2D(this,tpt_props_heterog)
            implicit none !> Enforce explicit variable declarations
            class(transport_2D_transient_c) :: this !> Transport equation object (modified in place)
            class(tpt_props_heterog_2D_c), intent(in)  :: tpt_props_heterog !> Input transport properties to assign
            this%tpt_props_heterog=tpt_props_heterog !> Copy transport properties to object's internal storage
        end subroutine
        
        !> \brief Set concentration reaction flags from source term
        !> \param[in,out] this Transport problem object
        !> \details 
        !>   Identifies cells with reactive source terms (r > 0).
        !>   Sets conc_r_flag(i) = 1 for reactive cells, 0 otherwise.
        !>   
        !>   Used for:
        !>   - Optimizing reactive transport (skip non-reactive cells)
        !>   - Coupling transport with geochemistry
        !>   - Tracking reactive zones
        subroutine set_conc_r_flag_tpt(this)
            implicit none
            class(transport_1D_transient_c) :: this             !< Transport problem object (modified)
            
            integer(kind=4) :: i                                !< Loop index
            
            allocate(this%conc_obj%conc_r_flag(this%spatial_discr%Num_targets))  !< Allocate flag array
            this%conc_obj%conc_r_flag=0                         !< Initialize all to zero (non-reactive)
            do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag  !< Loop over all cells
                if (this%tpt_props_heterog%source_term(i)>0) then  !< Check if source term is positive
                    this%conc_obj%conc_r_flag(i)=1              !< Mark as reactive
                end if
            end do
        end subroutine
        
        !> \brief Check time step against critical stability limit
        !> \param[in,out] this Transport problem object
        !> \details 
        !>   Verifies that Δt ≤ Δt_crit for numerical stability.
        !>   Only applies to homogeneous time discretization.
        !>   
        !>   Critical Time Step:
        !>   \f[
        !>     \Delta t_{crit} = \frac{\phi \Delta x^2}{2D}
        !>   \f]
        !>   
        !>   Current implementation: Warning only (no error stop)
        subroutine check_Delta_t(this)
            implicit none
            class(transport_1D_transient_c) :: this             !< Transport problem object
            select type (time=>this%time_discr)                 !< Polymorphic dispatch on time discretization
            type is (time_discr_homog_c)                        !< Case: homogeneous time stepping
                if (time%Delta_t>this%stab_params_tpt%Delta_t_crit) then  !< Check if Δt exceeds critical value
                    !print *, "Critical time step: ", this%stab_params_tpt%Delta_t_crit  !< Optional: print warning
                    !error stop "You must reduce time step to have stability"  !< Optional: halt execution
                end if
            end select
        end subroutine

        !> \brief Fundamental solution of 1D transport equation
        !> \param[in] this    Transport problem object
        !> \param[in] c0      Background concentration [C]
        !> \param[in] M       Injected mass [M]
        !> \param[in] Delta_x Cell size [L]
        !> \param[in] x       Spatial coordinate [L]
        !> \param[in] mu      Source location [L]
        !> \param[in] t       Time [T]
        !> \return conc       Concentration at (x,t) [C]
        !> \details 
        !>   Analytical fundamental solution for ADE with instantaneous point source.
        !>   
        !>   Transport Solution:
        !>   \f[
        !>     c(x,t) = c_0 + \exp\left(\frac{q}{2D}(x-\mu-\frac{qt}{2\phi})\right) \cdot G_{diff}(x,\mu,t)
        !>   \f]
        !>   where G_{diff} is the diffusion fundamental solution.
        !>   
        !>   Used for:
        !>   - Analytical verification of numerical schemes
        !>   - Green's function approach
        !>   - Benchmark problems
        function fund_sol_tpt_eqn_1D(this,c0,M,Delta_x,x,mu,t) result(conc)
            class(transport_1D_transient_c), intent(in) :: this  !< Transport problem object
            real(kind=8), intent(in) :: c0                      !< Background concentration [C]
            real(kind=8), intent(in) :: M                       !< Injected mass [M]
            real(kind=8), intent(in) :: Delta_x                 !< Cell size [L]
            real(kind=8), intent(in) :: x                       !< Spatial position [L]
            real(kind=8), intent(in) :: mu                      !< Source location [L]
            real(kind=8), intent(in) :: t                       !< Time [T]
            real(kind=8) :: conc                                !< Computed concentration [C]
            
            real(kind=8), parameter :: eps=1d-12                !< Tolerance for numerical precision

                conc=exp(5d-1*(this%tpt_props_heterog%flux_cent(1)/this%tpt_props_heterog%disp_cent(1))*&
                    (x-mu-t*this%tpt_props_heterog%flux_cent(1)/(2*this%tpt_props_heterog%porosity(1))))  !< Advection correction factor
                conc=conc*fund_sol_diff_eqn_1D(this,M,Delta_x,x,mu,t)+c0  !< Multiply by diffusion solution, add background
        end function
        
        
        
        !> \brief Compute transport matrix (dispatcher)
        !> \param[in,out] this Transport problem object
        !> \details 
        !>   Dispatches to appropriate method based on Lagr_flag.
        !>   
        !>   Dispatch Logic:
        !>   - If Lagr_flag = .true.: Call compute_trans_mat_tpt_trans_Lagr
        !>   - If Lagr_flag = .false.: Call compute_trans_mat_tpt_trans_Euler
        !>   
        !>   Provides polymorphic interface for method selection.
        subroutine compute_trans_mat_tpt_trans(this)
            implicit none
            class(transport_1D_transient_c) :: this             !< Transport problem object (modified)
            if (this%Lagr_flag) then                            !< Check method flag
                call this%compute_trans_mat_tpt_trans_Lagr()    !< Use Lagrangian method
            else                                                !< Otherwise
                call this%compute_trans_mat_tpt_trans_Euler()   !< Use Eulerian method
            end if
        end subroutine
        
        !> \brief 2D transport transition matrix dispatch (Lagrangian/Eulerian)
        !> \param[in,out] this 2D Transport problem object
        !> \details 
        !>   Dispatch Logic:
        !>   - If Lagr_flag = .true.: Call compute_trans_mat_tpt_trans_Lagr_2D (dispersion only, symmetric T)
        !>   - If Lagr_flag = .false.: Call compute_trans_mat_tpt_trans_Euler_2D (advection + dispersion)
        subroutine compute_trans_mat_tpt_trans_2D_dispatch(this)
            implicit none
            class(transport_2D_transient_c) :: this
            if (this%Lagr_flag) then
                call this%compute_trans_mat_tpt_trans_Lagr_2D()
            else
                call this%compute_trans_mat_tpt_trans_Euler_2D()
            end if
        end subroutine
        
        subroutine compute_F_mat_tpt_trans_2D(this)
            implicit none
            class(transport_2D_transient_c) :: this
            
            integer(kind=4) :: i,num_cells,num_nodes
            real(kind=8) :: r_i !> radius (for radial symmetry)
    
            num_cells=this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
            num_nodes=num_cells+1
            !print *, n
     
            if (this%spatial_discr%adapt_ref==1) then
                deallocate(this%F_mat%diag)
                call this%F_mat%allocate_array(this%spatial_discr%Num_targets)
            end if
    
            if (this%spatial_discr%targets_flag==0) then !> we take cell porosity
                this%F_mat%diag=this%tpt_props_heterog%porosity
            else !> we interpolate to nodes
                this%F_mat%diag(1)=this%tpt_props_heterog%porosity(1) !> left boundary
                do i=2,num_cells
                    this%F_mat%diag(i)=0.5d0*(this%tpt_props_heterog%porosity(i)+&
                        this%tpt_props_heterog%porosity(i-1)) !> internal nodes
                end do
                this%F_mat%diag(num_nodes)=this%tpt_props_heterog%porosity(num_cells) !> right boundary
            end if
            this%F_mat_prev%diag=this%F_mat%diag !> Store previous F matrix diagonal for time integration schemes
        end subroutine
        
        subroutine compute_F_mat_tpt(this) !> diagonal matrix
        !> F_{ii}=phi_i or phi_{i+1/2} depending on the spatial discretization
            implicit none
            class(transport_1D_transient_c) :: this
    
            integer(kind=4) :: i,num_cells,num_nodes
            real(kind=8) :: r_i !> radius (for radial symmetry)
    
            num_cells=this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
            num_nodes=num_cells+1
            !print *, n
     
            if (this%spatial_discr%adapt_ref==1) then
                deallocate(this%F_mat%diag)
                call this%F_mat%allocate_array(this%spatial_discr%Num_targets)
            end if
    
            if (this%spatial_discr%targets_flag==0) then !> we take cell porosity
                this%F_mat%diag=this%tpt_props_heterog%porosity
            else !> we interpolate to nodes
                this%F_mat%diag(1)=this%tpt_props_heterog%porosity(1) !> left boundary
                do i=2,num_cells
                    this%F_mat%diag(i)=0.5d0*(this%tpt_props_heterog%porosity(i)+&
                        this%tpt_props_heterog%porosity(i-1)) !> internal nodes
                end do
                this%F_mat%diag(num_nodes)=this%tpt_props_heterog%porosity(num_cells) !> right boundary
            end if
            this%F_mat_prev%diag=this%F_mat%diag !> Store previous F matrix diagonal for time integration schemes
        end subroutine
        
end module