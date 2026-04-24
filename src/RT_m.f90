!> \file RT_m.f90
!> \brief reactive transport module
!> \details
!>   This module provides the complete framework for coupled chemistry-transport simulations.
!>   
!>   Reactive transport couples:
!>   - Solute transport (advection-dispersion equation)
!>   - Chemical reactions (equilibrium and kinetic)
!>   - WMA for numerical efficiency
!>   
!>   Numerical methods available:
!>   - **Lagrangian methods**: Track particles moving with flow
!>     * Ideal for advection-dominated systems
!>     * No numerical dispersion
!>     * Natural mass conservation
!>   - **Eulerian methods**: Fixed spatial grid
!>     * Better for diffusion-dominated systems
!>     * Allows complex boundary conditions
!>     * Requires CFL stability condition
!>   
!>   Time integration options:
!>   - Stationary (steady-state) problems
!>   - Transient (time-dependent) problems
!>   - Explicit vs implicit chemistry integration
!>   
!>   WMA approach:
!>   1. Transport step: Conservative mixing of particles or grid cells
!>   2. Reaction step: Reactive mixing and chemical equilibrium calculations
!>   3. Repeat for each time step
!>   
!>   Key assumptions:
!>   - Isochrone mesh for Lagrangian method
!>   
!>   Applications:
!>   - Groundwater contamination
!>   - Column experiments
!>   - Weathering profiles
!>   - Biogeochemical gradients
!>   
!> \author Jordi Petchamé-Guerrero
!> \date October 2025
module RT_m                                                                  !< One-dimensional reactive transport module with Lagrangian and Eulerian methods
    use chemistry_m, only: chemistry_c
    use PDE_transient_m, only: PDE_transient_c
    use aqueous_chemistry_m, only: aqueous_chemistry_c                         !< Import complete chemistry management system
    use transport_transient_m, only: transport_1D_transient_c, transport_2D_transient_c                   !< Import transient transport classes
    use transport_m, only: transport_1D_c!, transport_2D_c                                      !< Import stationary transport class for steady-state problems
    use time_discr_m, only: time_discr_homog_c, time_discr_heterog_c                               !< Import time discretization classes
    use spatial_discr_m, only: spatial_discr_c                                  !< Import spatial discretization class
    use target_water_m, only: target_water_c !< Import target water class
    implicit none                                                              !< Require explicit declaration of all variables
    save                                                                       !< Preserve module variables between procedure calls
    private                                                                   !< Default visibility is private; expose only specified types/procedures
    public :: move_particles_stat_flux_1D, move_particles_stat_flux_EC_1D, move_particles_stat_flux_2D           !< Public types and procedures accessible outside the module
    !> \brief Base reactive transport class
    !> \details
    !>   Abstract superclass for all reactive transport implementations.
    !>   
    !>   This class provides the foundation for 1D reactive transport modeling by
    !>   combining chemistry and transport processes. It serves as the base for:
    !>   - RT_1D_stat_c: Steady-state problems
    !>   - RT_1D_transient_c: Time-dependent problems
    !>   
    !>   Key component:
    !>   - chemistry: Complete chemical system (species, reactions, zones)
    !>   
    !>   Solver methods:
    !>   - Lagrangian with steady flux (particles move at constant velocity)
    !>   - Lagrangian with transient flux (particles accelerate/decelerate)
    !>   - Eulerian with transient flux (fixed grid, time-varying flow)
    !>   
    !>   All methods assume ideal solution behavior and lumped chemistry.
    type, public, abstract :: RT_c                                                   !< Base reactive transport class - abstract superclass for all RT implementations
        type(chemistry_c) :: chemistry                                        !< Chemistry manager object containing all chemical processes and data
        type(target_water_c), allocatable :: target_waters(:)                                    !< Target water object for boundary conditions
        !type(target_water_c), allocatable :: target_waters_old(:)                                    !< Target water object for boundary conditions
        type(target_water_c), allocatable :: target_waters_init(:)                                    !< Target water object for boundary conditions
    contains                                                                   !< Type-bound procedures for reactive transport operations
        procedure :: set_target_waters                  !< Set positions of target waters
        procedure :: set_chemistry                                    !< Initialize and configure chemistry system
        !procedure :: solve_RT_ideal_lump_Lagr_stat_flux           !< Solve RT using Lagrangian method with steady flux (lumped chemistry)
        procedure :: write_RT                                     !< Write reactive transport results to output files
        procedure :: compute_eq_react_rates_stat                     !< Introduce new particles at injection points or boundaries
        procedure :: set_target_solids_time                        !< Read time discretization parameters from input files
    end type !< End base reactive transport class

    type, public, abstract, extends(RT_c) :: RT_transient_c                             !< Stationary (steady-state) reactive transport subclass
        real(kind=8) :: Delta_t_crit                                          !< Critical time step for numerical stability (CFL condition)
        integer(kind=4) :: int_method_chem_reacts                             !< Time integration method for chemical reactions: 1=explicit, 2=implicit                                                                   !< Additional procedures specific to stationary problems
        !procedure :: set_transport_stat                              !< Configure transport for steady-state conditions
    contains                                                                   !< Additional procedures specific to transient problems
        procedure :: set_int_method_chem_reacts                              !< Set chemical reaction integration method (explicit vs implicit)
        procedure :: set_transport_trans                              !< Set chemical reaction integration method (explicit vs implicit)
        procedure :: compute_Delta_t_crit_RT                         !< Compute critical time step for numerical stability
        procedure :: check_Delta_t_RT                        !< Move Lagrangian particles under steady flux conditions
        procedure(solve_RT_trans_ideal_iface), deferred :: solve_RT_ideal_lump_Lagr_stat_flux   !< Deferred: solve RT 1D Lagrangian transient flux
        procedure(solve_RT_trans_ideal_iface), deferred :: solve_RT_ideal_lump_Lagr_trans_flux   !< Deferred: solve RT 1D Lagrangian transient flux
        procedure(solve_RT_trans_ideal_iface), deferred :: solve_RT_ideal_lump_Euler_trans_flux  !< Deferred: solve RT 1D Eulerian transient flux
        procedure(solve_RT_trans_ideal_iface), deferred :: solve_RT_ideal_cons_Lagr_stat_flux    !< Deferred: solve RT 1D consistent Lagrangian stationary flux
        procedure(solve_RT_trans_ideal_iface), deferred :: solve_RT_ideal_cons_Euler_stat_flux   !< Deferred: solve RT 1D consistent Eulerian stationary flux
        procedure(solve_RT_trans_ideal_iface), deferred :: solve_RT_ideal_lump_Euler_stat_flux   !< Deferred: solve RT 1D lumped Eulerian stationary flux
    end type                                                                   !< End stationary reactive transport class
!***************************************************************************************************************************************************!
    !> \brief Stationary 1D reactive transport class
    !> \details
    !>   Extends RT_c for steady-state (time-independent) problems.
    !>   
    !>   Used when:
    !>   - Flow field is constant in time
    !>   - Seeking equilibrium concentration profiles
    !>   - Boundary conditions don't change
    !>   
    !>   Contains stationary transport object with:
    !>   - Steady velocity field
    !>   - Time-independent dispersion
    !>   - Fixed boundary conditions
    type, public,extends(RT_c) :: RT_1D_stat_c                             !< Stationary (steady-state) reactive transport subclass
        type(transport_1D_c) :: transport                                     !< Stationary transport object for steady-state flow and transport
    contains                                                                   !< Additional procedures specific to stationary problems
        procedure :: set_transport_stat                              !< Configure transport for steady-state conditions
    end type                                                                   !< End stationary reactive transport class
!***************************************************************************************************************************************************!
    !> \brief Transient 1D reactive transport class
    !> \details
    !>   Extends RT_c for time-dependent problems - most commonly used.
    !>   
    !>   Features:
    !>   - Time-varying flow and transport
    !>   - Adaptive time stepping
    !>   - CFL stability control
    !>   - Lagrangian particle tracking
    !>   
    !>   Critical time step (Delta_t_crit):
    !>   - Ensures numerical stability (CFL condition)
    !>   - Computed from flow velocity and grid spacing
    !>   - Must satisfy: Delta_t <= Delta_t_crit for explicit methods
    !>   
    !>   Integration methods for chemistry:
    !>   - 1: Explicit Euler (fast, conditionally stable)
    !>   - 2: Implicit Euler (slow, unconditionally stable)
    !>   - 3: Other methods (to be implemented)
    !>   
    !>   Lagrangian particle operations:
    !>   - move_particles_stat_flux: Update positions under steady flow
    !>   - introduce_particle: Add new particles at inflow boundary
    type, public,extends(RT_transient_c) :: RT_1D_transient_c                        !< Transient (time-dependent) 1D reactive transport subclass - most commonly used
        type(transport_1D_transient_c) :: transport                           !< Transient transport object for time-dependent flow and transport
        !real(kind=8) :: Delta_t_crit                                          !< Critical time step for numerical stability (CFL condition)
        !integer(kind=4) :: int_method_chem_reacts                             !< Time integration method for chemical reactions: 1=explicit, 2=implicit
    contains                                                                   !< Additional procedures specific to transient problems
        procedure :: move_particles_stat_flux_EC => move_particles_stat_flux_EC_1D  !< Set chemical reaction integration method (explicit vs implicit)
        !procedure :: set_transport_trans                             !< Configure transport for transient conditions
        procedure :: move_particles_stat_flux => move_particles_stat_flux_1D    !< Move Lagrangian particles under steady flux conditions
        procedure :: introduce_particle                             !< Introduce new particle at injection points or boundaries
        procedure :: solve_RT_ideal_lump_Lagr_stat_flux => solve_RT_ideal_lump_Lagr_stat_flux_1D           !< Override: solve RT 1D Lagrangian stationary flux
        procedure :: solve_RT_ideal_lump_Lagr_trans_flux=>solve_RT_ideal_lump_Lagr_trans_flux_1D          !< Override: solve RT 1D Lagrangian transient flux
        procedure :: solve_RT_ideal_lump_Euler_trans_flux=>solve_RT_ideal_lump_Euler_trans_flux_1D         !< Override: solve RT 1D Eulerian transient flux
        procedure :: solve_RT_ideal_cons_Lagr_stat_flux=>solve_RT_ideal_cons_Lagr_stat_flux_1D           !< Override: solve RT 1D consistent Lagrangian stationary flux
        procedure :: solve_RT_ideal_cons_Euler_stat_flux=>solve_RT_ideal_cons_Euler_stat_flux_1D          !< Override: solve RT 1D consistent Eulerian stationary flux
        procedure :: solve_RT_ideal_lump_Euler_stat_flux=>solve_RT_ideal_lump_Euler_stat_flux_1D          !< Override: solve RT 1D lumped Eulerian stationary flux
    end type                                                                   !< End transient 1D reactive transport class
!***************************************************************************************************************************************************!
    !> \brief Transient 2D reactive transport class
    !> \details
    !>   Extends RT_c for time-dependent problems - most commonly used.
    !>   
    !>   Features:
    !>   - Time-varying flow and transport
    !>   - Adaptive time stepping
    !>   - CFL stability control
    !>   - Lagrangian particle tracking
    !>   
    !>   Critical time step (Delta_t_crit):
    !>   - Ensures numerical stability (CFL condition)
    !>   - Computed from flow velocity and grid spacing
    !>   - Must satisfy: Delta_t <= Delta_t_crit for explicit methods
    !>   
    !>   Integration methods for chemistry:
    !>   - 1: Explicit Euler (fast, conditionally stable)
    !>   - 2: Implicit Euler (slow, unconditionally stable)
    !>   - 3: Other methods (to be implemented)
    !>   
    !>   Lagrangian particle operations:
    !>   - move_particles_stat_flux: Update positions under steady flow
    !>   - introduce_particle: Add new particles at inflow boundary
    type, public,extends(RT_transient_c) :: RT_2D_transient_c                        !< Transient (time-dependent) 2D reactive transport subclass - most commonly used
        type(transport_2D_transient_c) :: transport                           !< Transient transport object for time-dependent flow and transport
        !real(kind=8) :: Delta_t_crit                                          !< Critical time step for numerical stability (CFL condition)
        !integer(kind=4) :: int_method_chem_reacts                             !< Time integration method for chemical reactions: 1=explicit, 2=implicit
    contains                                                                   !< Additional procedures specific to transient problems
        !procedure :: move_particles_stat_flux_EC                              !< Set chemical reaction integration method (explicit vs implicit)
        !procedure :: set_transport_trans                             !< Configure transport for transient conditions
        procedure :: move_particles_stat_flux=>move_particles_stat_flux_2D                        !< Move Lagrangian particles under steady flux conditions
        procedure :: introduce_particles => introduce_particles_2D  !< Introduce new particle at injection points or boundaries
        procedure :: remap_mix_indices => remap_mix_indices_2D        !< Remap transport-local water indices to global chemistry water indices
        procedure :: solve_RT_ideal_lump_Lagr_stat_flux => solve_RT_ideal_lump_Lagr_stat_flux_2D         !< Get transport object for transient conditions
        procedure :: solve_RT_ideal_lump_Lagr_trans_flux => solve_RT_ideal_lump_Lagr_trans_flux_2D   !< 2D stub
        procedure :: solve_RT_ideal_lump_Euler_trans_flux => solve_RT_ideal_lump_Euler_trans_flux_2D !< 2D stub
        procedure :: solve_RT_ideal_cons_Lagr_stat_flux => solve_RT_ideal_cons_Lagr_stat_flux_2D     !< 2D stub
        procedure :: solve_RT_ideal_cons_Euler_stat_flux => solve_RT_ideal_cons_Euler_stat_flux_2D   !< 2D stub
        procedure :: solve_RT_ideal_lump_Euler_stat_flux => solve_RT_ideal_lump_Euler_stat_flux_2D   !< 2D stub
    end type                                                                   !< End transient reactive transport class
    !***************************************************************************************************************************************************!
    abstract interface
        subroutine solve_RT_trans_ideal_iface(this, dir, root)
            import RT_transient_c
            class(RT_transient_c) :: this
            character(len=*), intent(in) :: dir
            character(len=*), intent(in) :: root
        end subroutine
    end interface
    !***************************************************************************************************************************************************!
    !> \brief External procedure interfaces for reactive transport operations
    !> \details
    !>   These interfaces define external procedures for solving RT problems,
    !>   I/O operations, and time discretization management.
    interface
        !> \brief Solve RT using ideal lumped Eulerian method with transient flux
        !> \param[in,out] this RT_c object
        !> \param[in]     root Root name for output files
        !> \details
        !>   Fixed-grid Eulerian approach with time-varying flow.
        !>   Concentrations computed at fixed spatial points, ideal solution assumed.
        subroutine solve_RT_ideal_lump_Euler_trans_flux_1D(this,dir,root)
            import RT_1D_transient_c
            implicit none
            class(RT_1D_transient_c) :: this
            character(len=*), intent(in) :: dir !< Directory for output file [string]
            character(len=*), intent(in) :: root !> root name for output file
        end subroutine

        subroutine solve_RT_ideal_cons_Lagr_stat_flux_1D(this,dir,root)
            import RT_1D_transient_c
            implicit none
            class(RT_1D_transient_c) :: this
            character(len=*), intent(in) :: dir !< Directory for output file [string]
            character(len=*), intent(in) :: root !< Root name for output file (base filename without extension) [string]
        end subroutine
        
        !> \brief Solve RT using ideal lumped Lagrangian method with steady flux
        !> \param[in,out] this RT_c object
        !> \param[in]     root Root name for output files
        !> \details
        !>   Lagrangian particle tracking under steady-state flux conditions.
        !>   Particles move at constant velocity determined by flow field.
        !>   Ideal solution and lumped chemistry assumed.
        subroutine solve_RT_ideal_lump_Lagr_stat_flux_1D(this,dir,root)
            import RT_1D_transient_c
            implicit none
            class(RT_1D_transient_c) :: this
            character(len=*), intent(in) :: dir !> directory for output file
            character(len=*), intent(in) :: root !> root name for output file
        end subroutine

        subroutine solve_RT_ideal_lump_Lagr_stat_flux_2D(this,dir,root)
            import RT_2D_transient_c
            implicit none
            class(RT_2D_transient_c) :: this
            character(len=*), intent(in) :: dir !> directory for output file
            character(len=*), intent(in) :: root !> root name for output file
        end subroutine
        
        subroutine solve_RT_ideal_lump_Euler_stat_flux_1D(this,dir,root)
            import RT_1D_transient_c
            implicit none
            class(RT_1D_transient_c) :: this
            character(len=*), intent(in) :: dir !> directory for output file
            character(len=*), intent(in) :: root !> root name for output file
        end subroutine
        
        subroutine solve_RT_ideal_cons_Euler_stat_flux_1D(this,dir,root)
            import RT_1D_transient_c
            implicit none
            class(RT_1D_transient_c) :: this
            character(len=*), intent(in) :: dir !> directory for output file
            character(len=*), intent(in) :: root !> root name for output file
        end subroutine
        
        !> \brief Solve RT using ideal lumped Lagrangian method with transient flux
        !> \param[in,out] this RT_c object
        !> \param[in]     root Root name for output files
        !> \details
        !>   Lagrangian particle tracking under time-varying flux.
        !>   Particles accelerate/decelerate according to changing velocity field.
        !>   Ideal solution and lumped chemistry assumed.
        subroutine solve_RT_ideal_lump_Lagr_trans_flux_1D(this,dir,root)
            import RT_1D_transient_c
            implicit none
            class(RT_1D_transient_c) :: this
            character(len=*), intent(in) :: dir !> directory for output file
            character(len=*), intent(in) :: root !> root name for output file
        end subroutine
        
        !> \brief Write reactive transport results to output files
        !> \param[in] this    RT_c object
        !> \param[in] root    Root name for output files
        !> \param[in] path_py Optional path for Python visualization scripts
        !> \details
        !>   Writes concentration profiles, reaction rates, and other RT data.
        !>   Optionally generates Python visualization scripts.
        subroutine write_RT(this,dir,root,path_py)
            import RT_c
            implicit none
            class(RT_c), intent(in) :: this
            !integer(kind=4), intent(in) :: unit
            character(len=*), intent(in) :: dir !> directory for output files
            character(len=*), intent(in) :: root !> root name for output files
            character(len=*), intent(in), optional :: path_py !> path for python scripts
        end subroutine
        
        !> \brief Generate Python visualization scripts
        !> \param[in] this RT_c object
        !> \param[in] path Path for output files
        !> \details
        !>   Creates Python scripts for plotting concentration profiles,
        !>   reaction rates, and spatiotemporal evolution.
        !> \brief Compute critical time step for CFL stability
        !> \param[in,out] this RT_c object
        !> \details
        !>   Computes maximum stable time step based on CFL condition:
        !>   \f[
        !>     \Delta t_{\text{crit}} = \min_i \frac{\Delta x_i}{|v_i| + \sqrt{D_i/\Delta x_i}}
        !>   \f]
        !>   where \f$v_i\f$ is velocity, \f$D_i\f$ is dispersion, \f$\Delta x_i\f$ is cell width.
        !>   Updates this%Delta_t_crit internally.
        subroutine compute_Delta_t_crit_RT(this)
            import RT_transient_c
            implicit none
            class(RT_transient_c) :: this
        end subroutine
        
        !> \brief Read transport data
        !> \param[in,out] this RT_c object
        !> \param[in]     unit     File unit number
        !> \param[in]     file_tpt Transport data filename
        !> \details
        !>   Reads transport-related data: velocity profiles, dispersion coefficients,
        !>   boundary conditions, and mixing ratios.
       subroutine read_transport_data(this,unit,file_tpt,mixing_ratios)
            import RT_c
            class(RT_c) :: this
            integer(kind=4), intent(in) :: unit
            character(len=*), intent(in) :: file_tpt
       end subroutine
        
        !subroutine write_transport_data(this,unit)
        !    import RT_c
        !    class(RT_c) :: this
        !    integer(kind=4), intent(in) :: unit
        !end subroutine
       
       subroutine compute_eq_react_rates_stat(this)
            import RT_c
            implicit none
            class(RT_c) :: this
       end subroutine

    end interface
    
    contains
        !> \brief Set stationary transport object
        !> \param[in,out] this          RT_1D_stat_c object
        !> \param[in]     transport_obj Transport object for steady-state
        !> \details
        !>   Assigns the stationary transport object containing:
        !>   - Steady velocity field
        !>   - Time-independent dispersion coefficients
        !>   - Fixed boundary conditions
        subroutine set_transport_stat(this,transport_obj)
            implicit none
            class(RT_1D_stat_c) :: this
            class(transport_1D_c), intent(in) :: transport_obj
            this%transport=transport_obj
        end subroutine
        
        !> \brief Set transient transport object
        !> \param[in,out] this          RT_1D_transient_c object
        !> \param[in]     transport_obj Transport object for time-dependent flow
        !> \details
        !>   Assigns the transient transport object containing:
        !>   - Time-varying velocity field
        !>   - Dynamic dispersion coefficients
        !>   - Time-dependent boundary conditions
        subroutine set_transport_trans(this, transport_obj)
        implicit none
        class(RT_transient_c) :: this
        class(PDE_transient_c), intent(in) :: transport_obj

        select type (this)
        type is (RT_1D_transient_c)
            select type (transport_obj)
            class is (transport_1D_transient_c)
                if (.not. transport_obj%is_initialized()) &
                    error stop "set_transport_trans: transport_1D object is not initialized"
                this%transport = transport_obj
            class default
                error stop "Transport object type not compatible with RT_1D_transient_c"
            end select
        type is (RT_2D_transient_c)
            select type (transport_obj)
            class is (transport_2D_transient_c)
                if (.not. transport_obj%is_initialized()) &
                    error stop "set_transport_trans: transport_2D object is not initialized"
                this%transport = transport_obj
            class default
                error stop "Transport object type not compatible with RT_2D_transient_c"
            end select
        end select
    end subroutine
        
        !> \brief Set chemistry system
        !> \param[in,out] this     RT_c object
        !> \param[in]     chem_obj Chemistry object
        !> \details
        !>   Assigns the chemistry object containing:
        !>   - Chemical species definitions
        !>   - Reaction network (equilibrium + kinetic)
        !>   - Reactive zones
        !>   - Thermodynamic database
        subroutine set_chemistry(this,chem_obj)
            implicit none
            class(RT_c) :: this
            class(chemistry_c), intent(in) :: chem_obj
            this%chemistry=chem_obj
        end subroutine
      
        
        !> \brief Check if current time step satisfies stability criteria
        !> \param[in,out] this RT_c object
        !> \details
        !>   Validates that Delta_t <= Delta_t_crit to ensure CFL stability.
        !>   Aborts simulation with error message if condition violated.
        !>   
        !>   Checks:
        !>   - CFL number: \f$ \text{CFL} = \frac{v \Delta t}{\Delta x} \leq 1 \f$
        !>   - Diffusion number: \f$ D_n = \frac{D \Delta t}{\phi \Delta x^2} \leq 0.5 \f$
        subroutine check_Delta_t_RT(this)
            implicit none
            class(RT_transient_c) :: this    
            
            real(kind=8), parameter :: eps=1d-16
            
            !> Explicit method requires stability check
            select type (this)
            class is (RT_1D_transient_c)
                if (this%transport%time_discr%int_method.eq.1) then !> Check only for explicit methods
                    call this%compute_Delta_t_crit_RT()              !> Compute critical time step
                    if (abs(this%Delta_t_crit)<eps) then
                        continue !> Near-zero critical timestep, skip check
                    else
                        select type (time=>this%transport%time_discr)
                        type is (time_discr_homog_c)
                            if (time%Delta_t>this%Delta_t_crit) then !> CFL violation detected
                                print *, this%Delta_t_crit
                                error stop "Delta_t is larger than Delta_t_crit"
                            end if
                        end select
                    end if
                end if
            end select
        end subroutine
        
        !> \brief Set integration method for chemical reactions based on theta_r
        !> \param[in,out] this  RT_transient_c object
        !> \details
        !>   Reads theta_r from the transport time discretization and sets:
        !>   - theta_r == 0: Explicit Euler (int_method_chem_reacts = 1)
        !>   - theta_r /= 0: Implicit Euler (int_method_chem_reacts = 2)
        subroutine set_int_method_chem_reacts(this)
            implicit none
            class(RT_transient_c) :: this
            
            real(kind=8) :: theta_r
            
            select type (this)
            class is (RT_1D_transient_c)
                theta_r=this%transport%time_discr%theta_r
            class is (RT_2D_transient_c)
                theta_r=this%transport%time_discr%theta_r
            class default
                error stop "set_int_method_chem_reacts: unsupported RT type"
            end select
            
            if (theta_r==0d0) then
                this%int_method_chem_reacts=1 !> Explicit Euler
            else
                this%int_method_chem_reacts=2 !> Implicit Euler
            end if
        end subroutine

        !> \brief Move Lagrangian particles under steady flux in 1D reactive transport
        !> 
        !> \details
        !> This subroutine updates particle positions in a Lagrangian framework for transient
        !> reactive transport with steady-state flow. It performs a particle shifting operation
        !> where all particles advance one position downstream, creating space for new inflow
        !> particles and removing particles that exit at the outflow boundary.
        !> 
        !> **Physical Interpretation:**
        !> - Particles represent fluid parcels advecting through the domain
        !> - Each particle carries its own chemical composition (aqueous + solid phases)
        !> - Under steady flux, particles move at constant velocity between grid points
        !> - Position update follows: \f$ x_p(t+\Delta t) = x_p(t) + v \Delta t \f$
        !> 
        !> **Numerical Algorithm:**
        !> 1. **Save current state:** Store existing particle array and indices
        !> 2. **Shift boundary water:** Move inflow boundary water to first domain position
        !> 3. **Shift domain waters:** Move each domain water one position downstream (reverse loop)
        !> 4. **Update outflow:** Move last domain water to outflow boundary position
        !> 5. **Transfer solid chemistry:** Link each moved particle to its new target solid
        !> 6. **Remove exited particle:** Reallocate arrays with one fewer particle
        !> 7. **Reassign indices:** Update domain and boundary water index arrays
        !> 
        !> **Particle Movement Pattern:**
        !> ```
        !> Before:  [BD_in] [D1] [D2] [D3] ... [Dn] [BD_out]
        !>                    ↓    ↓    ↓         ↓      ↓
        !> After:            [D1] [D2] [D3] ... [Dn] [BD_out] (BD_in removed, new space created)
        !> ```
        !> 
        !> **Index Management:**
        !> - `old_bd_indices(1)`: Initial inflow boundary water index
        !> - `old_bd_indices(2)`: Initial outflow boundary water index  
        !> - `old_tar_wat_indices(i)`: Initial domain water indices (i=1 to num_target_waters)
        !> - After movement, inflow boundary index becomes first domain position
        !> - Last domain water becomes new outflow boundary water
        !> - Total particle count decreases by 1 (one exits at outflow)
        !> 
        !> **Solid Chemistry Coupling:**
        !> - Each particle is associated with a target solid via `tar%id`
        !> - When particle moves to new position, it inherits solid chemistry from that target
        !> - Accumulated reaction rates (`Rk_accum`) are transferred with the solid
        !> - Kinetic reaction indices are recalculated after position update
        !> 
        !> **Assumptions and Restrictions:**
        !> - Isochrone mesh: Particles initialized at regular grid intervals
        !> - Steady flux: Velocity field constant in time (transient chemistry only)
        !> - Outflow boundary: Particles can exit domain (removed from simulation)
        !> - Incompressible flow: Particle mass conserved during advection
        !> - One-to-one mapping: Each particle corresponds to exactly one target solid
        !> 
        !> **Error Handling:**
        !> - Terminates if particle moves to target_id = 0 (non-existing target)
        !> - Diagnostic prints show particle movement for debugging
        !> 
        !> \param[in,out] this Transient reactive transport object (RT_1D_transient_c)
        !>                     - `this%chemistry%waters`: Particle array (aqueous chemistry)
        !>                     - `this%chemistry%tar_wat_indices`: Domain water indices
        !>                     - `this%chemistry%bd_waters_indices`: Boundary water indices
        !>                     - `this%chemistry%num_waters`: Total particle count (decremented)
        !>                     - `this%chemistry%target_solids`: Associated solid phases
        !> \param[in]     k    Current time step index (used for diagnostics only)
        !> 
        !> \pre
        !> - Chemistry object must be initialized with target waters and solids
        !> - Domain and boundary indices must be allocated and valid
        !> - Each domain water must have valid solid_chemistry%tar%id > 0
        !> - Isochrone mesh structure must be maintained
        !> 
        !> \post
        !> - All particles shifted one position downstream
        !> - Total particle count reduced by 1
        !> - Inflow boundary water removed from particle array
        !> - Last domain water becomes new outflow boundary water
        !> - Domain indices updated to reflect new configuration
        !> - Solid chemistry properly linked to new particle positions
        !> 
        !> \warning
        !> - This subroutine assumes steady flux (transient chemistry with stationary flow)
        !> - Not suitable for time-varying velocity fields (use move_particles_trans_flux)
        !> - Particle removal at outflow may violate mass conservation if not balanced by inflow
        !> - Reverse loop (i = num_dom down to 1) is critical to avoid overwriting data
        !> 
        !> \note
        !> - Companion subroutine `introduce_particle` must be called to add new inflow particle
        !> - Diagnostic prints can be commented out for production runs
        !> - Target ID check could be disabled for performance if mesh is guaranteed valid
        !> 
        !> \see introduce_particle, set_solid_chemistry, set_indices_rk, copy_aq_chem
        subroutine move_particles_stat_flux_1D(this,k)
            !> Move particles to the new coordinates based on the current transport spatial discretisation
            !> This is used in the transient case to update particle positions after each time step.
            !> It also checks if the particles have left the domain and removes them from the chemistry object
            !> It assumes an isochrone mesh
            implicit none
            class(RT_1D_transient_c) :: this !> transient reactive transport object
            integer(kind=4), intent(in) :: k !> time step index
            
            integer(kind=4) :: target_id !> target id associated to the displaced particle
            integer(kind=4) :: num_exit !> number of particles that left the domain
            integer(kind=4) :: old_num_tar_wat !> old number of target waters
            real(kind=8), allocatable :: new_coords(:) !> new coordinates of the particles
            integer(kind=4) :: i !> loop index
            logical :: exit_flag !> flag to check if the particle has left the domain
            type(target_water_c), allocatable :: old_particles(:) !> old particles
            type(target_water_c), allocatable :: disp_particles(:) !> displaced particles
            integer(kind=4), allocatable :: old_tar_wat_indices(:) !> old particles indices
            integer(kind=4), allocatable :: old_bd_indices(:) !> old particles indices
            type(aqueous_chemistry_c) :: aq_chem_exit !> aqueous chemistry of the particle that left the domain
            
            num_exit=0 !> initialise number of particles that left the domain
            
            allocate(new_coords(size(this%transport%spatial_discr%targets(1)%coord))) !> allocate new coordinates array
            !allocate(old_particles(this%chemistry%num_target_waters)) !> allocate old particles array
            !allocate(disp_particles(this%chemistry%num_target_waters-1)) !> allocate old domain indices array
            allocate(old_tar_wat_indices(this%chemistry%num_target_waters)) !> allocate old domain indices array
            allocate(old_bd_indices(this%chemistry%num_bd_waters)) !> allocate old boundary indices array
            ! do i=1,this%chemistry%num_target_waters
            !     call old_particles(i)%copy_tar_wat(this%target_waters(i)) !> get old particles
            ! end do
            old_tar_wat_indices=this%chemistry%tar_wat_indices !> get old particle indices
            !old_num_tar_wat=this%chemistry%num_target_waters !> get old number of domain waters
            old_bd_indices=this%chemistry%bd_waters_indices !> get old boundary indices
            !call this%chemistry%allocate_tar_wat_indices(this%chemistry%num_target_waters-1) !> allocate new indices array
            ! this%chemistry%tar_wat_indices=&
            !     old_tar_wat_indices(1:this%chemistry%num_target_waters) !> set new domain indices
            ! this%chemistry%waters(old_bd_indices(1))%pos=&
            !     old_particles(old_tar_wat_indices(1))%pos !> we move the inflow boundary water to the position of the first domain water
            ! target_id=old_particles(old_tar_wat_indices(1))%solid_chemistry%get_tar_id() !> we get the id of the first domain target
            ! print *, "Particle ", this%chemistry%waters(old_bd_indices(1))%id, & 
            !         " moved to target ", target_id
            ! call this%chemistry%waters(old_bd_indices(1))%set_solid_chemistry(&
            !     this%chemistry%target_solids(target_id)) !> we set the solid chemistry of the displaced inflow boundary water
            ! call this%chemistry%waters(old_bd_indices(1))%set_indices_rk() !> we set the indices of the kinetic reactions to the displaced inflow boundary water
            ! this%chemistry%waters(old_bd_indices(1))%concentrations=&
            !     this%chemistry%wat_types(1)%concentrations !+ &
                !this%chemistry%waters(old_tar_wat_indices(1))%concentrations/5d0 !> we set the concentrations of the displaced inflow boundary water
            ! this%chemistry%waters(old_tar_wat_indices(this%chemistry%num_target_waters))%pos=&
            !     this%chemistry%waters(old_bd_indices(2))%pos !> we move the last domain water to the position of the outflow boundary water
            ! target_id=this%chemistry%waters(old_bd_indices(2))%solid_chemistry%tar%id !> we get the id of the last target
            ! print *, "Water ", this%chemistry%waters(old_tar_wat_indices(this%chemistry%num_target_waters))%id, & 
            !         " moved to target ", target_id
            ! print *, "Particle ", this%target_waters(this%chemistry%num_target_waters)%id, & 
            !     " left the domain "
            call aq_chem_exit%copy_aq_chem(this%target_waters(&
                this%chemistry%num_target_waters)%aq_chem) !> water that left the domain
            ! call this%chemistry%waters(old_tar_wat_indices(this%chemistry%num_target_waters))%set_solid_chemistry(&
            !     this%chemistry%target_solids(target_id)) !> we set the solid chemistry of the new outflow boundary water
            ! call this%chemistry%waters(old_tar_wat_indices(this%chemistry%num_target_waters))%set_indices_rk() !> we set the indices of the kinetic reactions to the new outflow boundary water
            
            do i=1,this%chemistry%num_target_waters-1 !> loop over all target waters (except first one) in reverse order
                call this%target_waters(&
                    this%chemistry%num_target_waters-i+1)%aq_chem%copy_aq_chem_Lagr(&
                    this%target_waters(this%chemistry%num_target_waters-i)%aq_chem) !> copy displaced particle
                call this%target_waters(&
                    this%chemistry%num_target_waters-i+1)%set_id(&
                    this%target_waters(this%chemistry%num_target_waters-i)%id) !> set particle id
                !> compute new position of the particle
                ! new_coords=this%target_waters(old_num_tar_wat-i+1)%pos_old
                ! target_id=old_particles(old_num_tar_wat-i+1)%aq_chem%solid_chemistry_old%tar%id
                ! if (target_id>=0 .and. target_id<=this%chemistry%num_target_solids) then
                !     !print *, "Rk accumulated in target ", target_id, ":", this%chemistry%target_solids(target_id)%Rk_accum
                    ! call this%target_waters(&
                    !     old_num_tar_wat-i+1)%aq_chem%set_solid_chemistry(&
                    !     this%target_waters(old_num_tar_wat-i+1)%solid_chemistry_old)
                !     !print *, "Rk accumulated in target ", target_id, ":", this%chemistry%waters(&
                !     !    old_tar_wat_indices(this%chemistry%num_target_waters-i))%solid_chemistry%Rk_accum
                ! else
                !     error stop "Particle moved to a non-existing target"
                !     !call this%transport%spatial_discr%check_exit(new_coords,exit_flag)
                !     !if (exit_flag) then
                !     !    if (this%chemistry%waters(this%chemistry%num_waters-i+1)%pos(1)<this%chemistry%waters(&
                !     !        this%chemistry%bd_waters_indices(2))%pos(1)) then
                !     !        print *, "New outflow particle:", this%chemistry%num_waters-i+1
                !     !    else
                !     !        print *, "Particle left the mesh at time step ", k
                !     !        print *, "Particle index: ", this%chemistry%num_waters-i+1
                !     !        num_exit=num_exit+1 !> increment number of particles that left the mesh
                !     !    end if
                !     !else
                !     !    error stop "Particle moved to a non-existing target"
                !     !end if
                ! end if
                ! call this%target_waters(&
                !     old_num_tar_wat-i+1)%set_pos(new_coords)
                !call this%chemistry%waters(old_tar_wat_indices(&
                !    old_num_tar_wat-i))%set_indices_rk()
                ! this%target_waters(old_num_tar_wat-i+1)%aq_chem%indices_rk=&
                !     this%target_waters(old_num_tar_wat-i+1)%aq_chem%indices_rk_old !> we set the indices of the kinetic reactions to the displaced particles
                ! this%target_waters(old_num_tar_wat-i+1)%aq_chem%ind_prim_species=&
                !     this%target_waters(old_num_tar_wat-i+1)%aq_chem%ind_prim_species_old !> we set the concentrations of the displaced particles
                ! this%target_waters(old_num_tar_wat-i+1)%aq_chem%ind_sec_species=&
                !     this%target_waters(old_num_tar_wat-i+1)%aq_chem%ind_sec_species_old !> we set the concentrations of the displaced particles
                ! print *, "Particle ", this%target_waters(&
                !     this%chemistry%num_target_waters-i+1)%id, & 
                !     " moved to target ", this%target_waters(this%chemistry%num_target_waters-i+1)%aq_chem%solid_chemistry%tar%id
            end do
            if (this%transport%BCs%labels(2) .eq. 2) then !> Neumann outflow BC
                call this%chemistry%waters(this%chemistry%bd_waters_indices(2))%copy_aq_chem(&
                    aq_chem_exit) !> Copy outflow particle
            !else if (this%transport%BCs%labels(2) .eq. 1) then !> Dirichlet outflow BC: boundary waters remain fixed
            end if 
            !call disp_particles(1)%copy_aq_chem(this%chemistry%waters(old_bd_indices(1))) !> get displaced inflow boundary water
            ! do i=1,this%chemistry%num_target_waters-1
            !     call disp_particles(i)%copy_aq_chem(this%chemistry%waters(&
            !         old_tar_wat_indices(i))) !> get displaced particles
            ! end do
            !print *, this%chemistry%waters_init(3)%id
            !call this%chemistry%allocate_waters(this%chemistry%num_waters-1) !> allocate new particles array
            
            !this%chemistry%tar_wat_indices(1)=old_bd_indices(1) !> set the new domain target water index
            !this%chemistry%bd_waters_indices(2)=old_tar_wat_indices(this%chemistry%num_target_waters) !> set new outflow boundary water index
            !this%chemistry%bd_waters_indices(1)=0 !> no inflow boundary water
            ! print *, "New domain indices: ", this%chemistry%tar_wat_indices
            ! print *, "New boundary indices: ", this%chemistry%bd_waters_indices
            ! do i=1,old_bd_indices(2)-1
            !     call this%chemistry%waters(i)%copy_aq_chem(old_particles(i)) !> set old particles to the chemistry object
            !     !print *, this%chemistry%waters(i)%id !> set old particles to the chemistry object
            ! end do
            ! do i=1,this%chemistry%num_target_waters
            !     call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%copy_aq_chem(&
            !         disp_particles(i)) !> set displaced particles to the domain target waters
            ! end do
            ! call this%chemistry%waters(this%chemistry%bd_waters_indices(2))%copy_aq_chem(&
            !     disp_particles(old_tar_wat_indices(this%chemistry%num_target_waters))) !> set boundary outflow water
            ! print *, "New domain indices: ", this%chemistry%tar_wat_indices
            ! !print *, "New boundary indices: ", this%chemistry%bd_waters_indices
            ! print *, "Total number of particles: ", this%chemistry%num_target_waters
            !if (num_exit>0) then
            !    call this%chemistry%allocate_waters(this%chemistry%num_waters-num_exit) !> allocate new particles array
            !    call this%chemistry%allocate_tar_wat_indices(this%chemistry%num_target_waters-num_exit) !> allocate new indices array
            !    this%chemistry%tar_wat_indices=old_tar_wat_indices(1:size(old_tar_wat_indices)-num_exit) !> set new indices
            !    !call this%chemistry%allocate_bd_wat_indices(this%chemistry%num_target_waters-num_exit) !> allocate new indices array
            !    this%chemistry%bd_waters_indices(2)=this%chemistry%num_waters !> set new boundary water index
            !    this%chemistry%waters=old_particles(1:this%chemistry%num_waters) !> set new particles
            !    !this%chemistry%waters(this%chemistry%bd_waters_indices)=old_particles(old_bd_indices) !> set old particles to the chemistry object
            !    !this%chemistry%waters(this%chemistry%tar_wat_indices)=old_particles(old_tar_wat_indices(1:size(old_tar_wat_indices)-num_exit)) !> set old particles to the chemistry object
            !    print *, "Total number of particles: ", this%chemistry%num_waters
            !end if
            deallocate(new_coords) !> deallocate new coordinates array
            !deallocate(old_particles) !> deallocate old particles array
            deallocate(old_tar_wat_indices) !> deallocate old domain indices array
            deallocate(old_bd_indices) !> deallocate old boundary indices array
        end subroutine
        !> 
        !> \details
        !> This subroutine updates particle positions in a Lagrangian framework for transient
        !> reactive transport with steady-state flow. It performs a particle shifting operation
        !> where all particles advance one position downstream, creating space for new inflow
        !> particles and removing particles that exit at the outflow boundary.
        !> 
        !> **Physical Interpretation:**
        !> - Particles represent fluid parcels advecting through the domain
        !> - Each particle carries its own chemical composition (aqueous + solid phases)
        !> - Under steady flux, particles move at constant velocity between grid points
        !> - Position update follows: \f$ x_p(t+\Delta t) = x_p(t) + v \Delta t \f$
        !> 
        !> **Numerical Algorithm:**
        !> 1. **Save current state:** Store existing particle array and indices
        !> 2. **Shift boundary water:** Move inflow boundary water to first domain position
        !> 3. **Shift domain waters:** Move each domain water one position downstream (reverse loop)
        !> 4. **Update outflow:** Move last domain water to outflow boundary position
        !> 5. **Transfer solid chemistry:** Link each moved particle to its new target solid
        !> 6. **Remove exited particle:** Reallocate arrays with one fewer particle
        !> 7. **Reassign indices:** Update domain and boundary water index arrays
        !> 
        !> **Particle Movement Pattern:**
        !> ```
        !> Before:  [BD_in] [D1] [D2] [D3] ... [Dn] [BD_out]
        !>                    ↓    ↓    ↓         ↓      ↓
        !> After:            [D1] [D2] [D3] ... [Dn] [BD_out] (BD_in removed, new space created)
        !> ```
        !> 
        !> **Index Management:**
        !> - `old_bd_indices(1)`: Initial inflow boundary water index
        !> - `old_bd_indices(2)`: Initial outflow boundary water index  
        !> - `old_tar_wat_indices(i)`: Initial domain water indices (i=1 to num_target_waters)
        !> - After movement, inflow boundary index becomes first domain position
        !> - Last domain water becomes new outflow boundary water
        !> - Total particle count decreases by 1 (one exits at outflow)
        !> 
        !> **Solid Chemistry Coupling:**
        !> - Each particle is associated with a target solid via `tar%id`
        !> - When particle moves to new position, it inherits solid chemistry from that target
        !> - Accumulated reaction rates (`Rk_accum`) are transferred with the solid
        !> - Kinetic reaction indices are recalculated after position update
        !> 
        !> **Assumptions and Restrictions:**
        !> - Isochrone mesh: Particles initialized at regular grid intervals
        !> - Steady flux: Velocity field constant in time (transient chemistry only)
        !> - Outflow boundary: Particles can exit domain (removed from simulation)
        !> - Incompressible flow: Particle mass conserved during advection
        !> - One-to-one mapping: Each particle corresponds to exactly one target solid
        !> 
        !> **Error Handling:**
        !> - Terminates if particle moves to target_id = 0 (non-existing target)
        !> - Diagnostic prints show particle movement for debugging
        !> 
        !> \param[in,out] this Transient reactive transport object (RT_1D_transient_c)
        !>                     - `this%chemistry%waters`: Particle array (aqueous chemistry)
        !>                     - `this%chemistry%tar_wat_indices`: Domain water indices
        !>                     - `this%chemistry%bd_waters_indices`: Boundary water indices
        !>                     - `this%chemistry%num_waters`: Total particle count (decremented)
        !>                     - `this%chemistry%target_solids`: Associated solid phases
        !> \param[in]     k    Current time step index (used for diagnostics only)
        !> 
        !> \pre
        !> - Chemistry object must be initialized with target waters and solids
        !> - Domain and boundary indices must be allocated and valid
        !> - Each domain water must have valid solid_chemistry%tar%id > 0
        !> - Isochrone mesh structure must be maintained
        !> 
        !> \post
        !> - All particles shifted one position downstream
        !> - Total particle count reduced by 1
        !> - Inflow boundary water removed from particle array
        !> - Last domain water becomes new outflow boundary water
        !> - Domain indices updated to reflect new configuration
        !> - Solid chemistry properly linked to new particle positions
        !> 
        !> \warning
        !> - This subroutine assumes steady flux (transient chemistry with stationary flow)
        !> - Not suitable for time-varying velocity fields (use move_particles_trans_flux)
        !> - Particle removal at outflow may violate mass conservation if not balanced by inflow
        !> - Reverse loop (i = num_dom down to 1) is critical to avoid overwriting data
        !> 
        !> \note
        !> - Companion subroutine `introduce_particle` must be called to add new inflow particle
        !> - Diagnostic prints can be commented out for production runs
        !> - Target ID check could be disabled for performance if mesh is guaranteed valid
        !> 
        !> \see introduce_particle, set_solid_chemistry, set_indices_rk, copy_aq_chem
        subroutine move_particles_stat_flux_2D(this,k)
            !> Displace target waters horizontally (left to right) in a 2D mesh.
            !> Each row is processed independently so that a particle in the rightmost
            !> column of row r never wraps into the leftmost column of row r+1.
            !>
            !> Target waters are in row-major order: index = (r-1)*Nx + c
            !>   r = row (1..Ny, bottom to top), c = column (1..Nx, left to right).
            !>
            !> For each row:
            !>   1. The particle at column Nx exits the domain.
            !>   2. Particles at columns c = Nx..2 receive the chemistry of column c-1.
            !>   3. Column 1 (leftmost) is left for introduce_particles to fill.
            !>
            !> Boundary waters ordering (from input data):
            !>   bd_waters_indices(1..num_left_bd_waters)                     : left boundary
            !>   bd_waters_indices(num_left_bd_waters+1..2*num_left_bd_waters): right boundary
            implicit none
            class(RT_2D_transient_c) :: this !> transient reactive transport object
            integer(kind=4), intent(in) :: k !> time step index

            integer(kind=4) :: Nx !> number of cells in x-direction (columns)
            integer(kind=4) :: Ny !> number of cells in y-direction (rows)
            integer(kind=4) :: r  !> row loop index (1 = bottom, Ny = top)
            integer(kind=4) :: c  !> column loop index
            integer(kind=4) :: idx_src !> linear index of source target water
            integer(kind=4) :: idx_dst !> linear index of destination target water
            type(aqueous_chemistry_c), allocatable :: aq_chem_exit(:) !> exit chemistry per row

            Nx = this%transport%spatial_discr%get_num_cells(1)
            Ny = this%transport%spatial_discr%get_num_cells(2)

            !> 1. Save chemistry of rightmost column (these particles exit the domain)
            allocate(aq_chem_exit(Ny))
            do r = 1, Ny
                idx_src = r * Nx !> rightmost cell of row r
                call aq_chem_exit(r)%copy_aq_chem(this%target_waters(idx_src)%aq_chem)
                ! print *, "Particle ", this%target_waters(idx_src)%id, &
                !     " left the domain (row ", r, ")"
            end do

            !> 2. Shift particles horizontally within each row (right to left to avoid overwriting)
            do r = 1, Ny
                do c = Nx, 2, -1
                    idx_dst = (r - 1) * Nx + c
                    idx_src = (r - 1) * Nx + c - 1
                    call this%target_waters(idx_dst)%aq_chem%copy_aq_chem_Lagr( &
                        this%target_waters(idx_src)%aq_chem)
                    call this%target_waters(idx_dst)%set_id( &
                        this%target_waters(idx_src)%id)
                    ! print *, "  Shift row", r, ": target_water(", idx_src, &
                    !     ") -> target_water(", idx_dst, ")", &
                    !     " solid_chem tar_id=", &
                    !     this%target_waters(idx_dst)%aq_chem%solid_chemistry%tar%id
                end do
                !> Column 1 of each row will be filled by introduce_particles
            end do

            !> 3. Update right-boundary outflow waters (Neumann BC)
            if (this%transport%BCs%labels(2) .eq. 2) then
                do r = 1, this%chemistry%num_right_bd_waters
                    ! print *, "  Neumann outflow: row", r, &
                    !     ", bd_water index =", &
                    !     this%chemistry%bd_waters_indices( &
                    !     this%chemistry%num_left_bd_waters + &
                    !     this%chemistry%num_right_bd_waters - r + 1)
                    call this%chemistry%waters( &
                        this%chemistry%bd_waters_indices( &
                        this%chemistry%num_left_bd_waters + &
                        this%chemistry%num_right_bd_waters - r + 1))%copy_aq_chem( &
                        aq_chem_exit(r))
                end do
            end if

            deallocate(aq_chem_exit)
        end subroutine

        subroutine move_particles_stat_flux_EC_1D(this,k)
            !> Move particles to the new coordinates based on the current transport spatial discretisation
            !> This is used in the transient case to update particle positions after each time step.
            !> It also checks if the particles have left the domain and removes them from the chemistry object
            !> It assumes an isochrone mesh
            implicit none
            class(RT_1D_transient_c) :: this !> transient reactive transport object
            integer(kind=4), intent(in) :: k !> time step index
            
            integer(kind=4) :: target_id !> target id associated to the displaced particle
            integer(kind=4) :: num_exit !> number of particles that left the domain
            integer(kind=4) :: old_num_tar_wat !> old number of target waters
            real(kind=8), allocatable :: new_coords(:) !> new coordinates of the particles
            integer(kind=4) :: i !> loop index
            logical :: exit_flag !> flag to check if the particle has left the domain
            type(target_water_c), allocatable :: old_particles(:) !> old particles
            type(target_water_c), allocatable :: disp_particles(:) !> displaced particles
            integer(kind=4), allocatable :: old_tar_wat_indices(:) !> old particles indices
            integer(kind=4), allocatable :: old_bd_indices(:) !> old particles indices
            type(aqueous_chemistry_c) :: aq_chem_exit !> aqueous chemistry of the particle that left the domain
            
            num_exit=0 !> initialise number of particles that left the domain
            
            allocate(new_coords(size(this%transport%spatial_discr%targets(1)%coord))) !> allocate new coordinates array
            !allocate(old_particles(this%chemistry%num_target_waters)) !> allocate old particles array
            !allocate(disp_particles(this%chemistry%num_target_waters-1)) !> allocate old domain indices array
            allocate(old_tar_wat_indices(this%chemistry%num_target_waters)) !> allocate old domain indices array
            allocate(old_bd_indices(this%chemistry%num_bd_waters)) !> allocate old boundary indices array
            ! do i=1,this%chemistry%num_target_waters
            !     call old_particles(i)%copy_tar_wat(this%target_waters(i)) !> get old particles
            ! end do
            old_tar_wat_indices=this%chemistry%tar_wat_indices !> get old particle indices
            !old_num_tar_wat=this%chemistry%num_target_waters !> get old number of domain waters
            old_bd_indices=this%chemistry%bd_waters_indices !> get old boundary indices
            !call this%chemistry%allocate_tar_wat_indices(this%chemistry%num_target_waters-1) !> allocate new indices array
            ! this%chemistry%tar_wat_indices=&
            !     old_tar_wat_indices(1:this%chemistry%num_target_waters) !> set new domain indices
            ! this%chemistry%waters(old_bd_indices(1))%pos=&
            !     old_particles(old_tar_wat_indices(1))%pos !> we move the inflow boundary water to the position of the first domain water
            ! target_id=old_particles(old_tar_wat_indices(1))%solid_chemistry%get_tar_id() !> we get the id of the first domain target
            ! print *, "Particle ", this%chemistry%waters(old_bd_indices(1))%id, & 
            !         " moved to target ", target_id
            ! call this%chemistry%waters(old_bd_indices(1))%set_solid_chemistry(&
            !     this%chemistry%target_solids(target_id)) !> we set the solid chemistry of the displaced inflow boundary water
            ! call this%chemistry%waters(old_bd_indices(1))%set_indices_rk() !> we set the indices of the kinetic reactions to the displaced inflow boundary water
            ! this%chemistry%waters(old_bd_indices(1))%concentrations=&
            !     this%chemistry%wat_types(1)%concentrations !+ &
                !this%chemistry%waters(old_tar_wat_indices(1))%concentrations/5d0 !> we set the concentrations of the displaced inflow boundary water
            ! this%chemistry%waters(old_tar_wat_indices(this%chemistry%num_target_waters))%pos=&
            !     this%chemistry%waters(old_bd_indices(2))%pos !> we move the last domain water to the position of the outflow boundary water
            ! target_id=this%chemistry%waters(old_bd_indices(2))%solid_chemistry%tar%id !> we get the id of the last target
            ! print *, "Water ", this%chemistry%waters(old_tar_wat_indices(this%chemistry%num_target_waters))%id, & 
            !         " moved to target ", target_id
            ! print *, "Particle ", this%target_waters(this%chemistry%num_target_waters)%id, & 
            !         " left the domain "
            call aq_chem_exit%copy_aq_chem(this%target_waters(&
                this%chemistry%num_target_waters)%aq_chem) !> water that left the domain
            ! call this%chemistry%waters(old_tar_wat_indices(this%chemistry%num_target_waters))%set_solid_chemistry(&
            !     this%chemistry%target_solids(target_id)) !> we set the solid chemistry of the new outflow boundary water
            ! call this%chemistry%waters(old_tar_wat_indices(this%chemistry%num_target_waters))%set_indices_rk() !> we set the indices of the kinetic reactions to the new outflow boundary water
            
            do i=1,this%chemistry%num_target_waters-1 !> loop over all target waters (except first one) in reverse order
                call this%target_waters(&
                    this%chemistry%num_target_waters-i+1)%aq_chem%copy_aq_chem_Lagr(&
                    this%target_waters(this%chemistry%num_target_waters-i)%aq_chem) !> copy displaced particle
                call this%target_waters(&
                    this%chemistry%num_target_waters-i+1)%set_id(&
                    this%target_waters(this%chemistry%num_target_waters-i)%id) !> set particle id
                !> compute new position of the particle
                ! new_coords=this%target_waters(old_num_tar_wat-i+1)%pos_old
                ! target_id=old_particles(old_num_tar_wat-i+1)%aq_chem%solid_chemistry_old%tar%id
                ! if (target_id>=0 .and. target_id<=this%chemistry%num_target_solids) then
                !     !print *, "Rk accumulated in target ", target_id, ":", this%chemistry%target_solids(target_id)%Rk_accum
                    ! call this%target_waters(&
                    !     old_num_tar_wat-i+1)%aq_chem%set_solid_chemistry(&
                    !     this%target_waters(old_num_tar_wat-i+1)%solid_chemistry_old)
                !     !print *, "Rk accumulated in target ", target_id, ":", this%chemistry%waters(&
                !     !    old_tar_wat_indices(this%chemistry%num_target_waters-i))%solid_chemistry%Rk_accum
                ! else
                !     error stop "Particle moved to a non-existing target"
                !     !call this%transport%spatial_discr%check_exit(new_coords,exit_flag)
                !     !if (exit_flag) then
                !     !    if (this%chemistry%waters(this%chemistry%num_waters-i+1)%pos(1)<this%chemistry%waters(&
                !     !        this%chemistry%bd_waters_indices(2))%pos(1)) then
                !     !        print *, "New outflow particle:", this%chemistry%num_waters-i+1
                !     !    else
                !     !        print *, "Particle left the mesh at time step ", k
                !     !        print *, "Particle index: ", this%chemistry%num_waters-i+1
                !     !        num_exit=num_exit+1 !> increment number of particles that left the mesh
                !     !    end if
                !     !else
                !     !    error stop "Particle moved to a non-existing target"
                !     !end if
                ! end if
                ! call this%target_waters(&
                !     old_num_tar_wat-i+1)%set_pos(new_coords)
                !call this%chemistry%waters(old_tar_wat_indices(&
                !    old_num_tar-wat-i))%set_indices_rk()
                ! this%target_waters(old_num_tar-wat-i+1)%aq_chem%indices_rk=&
                !     this%target_waters(old_num_tar-wat-i+1)%aq_chem%indices_rk_old !> we set the indices of the kinetic reactions to the displaced particles
                ! this%target_waters(old_num_tar-wat-i+1)%aq_chem%ind_prim_species=&
                !     this%target_waters(old_num_tar-wat-i+1)%aq_chem%ind_prim_species_old !> we set the concentrations of the displaced particles
                ! this%target_waters(old_num_tar-wat-i+1)%aq_chem%ind_sec_species=&
                !     this%target_waters(old_num_tar-wat-i+1)%aq_chem%ind_sec_species_old !> we set the concentrations of the displaced particles
                ! print *, "Particle ", this%target_waters(&
                !     this%chemistry%num_target_waters-i+1)%id, & 
                !     " moved to target ", this%target_waters(this%chemistry%num_target_waters-i+1)%aq_chem%solid_chemistry%tar%id
            end do
            !> Chapuza para segunda target water
            this%target_waters(2)%aq_chem%concentrations=&
                this%target_waters(2)%aq_chem%concentrations/2d0 + &
                this%chemistry%waters(this%chemistry%bd_waters_indices(1))%concentrations/2d0
            this%target_waters(2)%aq_chem%activities=&
                this%target_waters(2)%aq_chem%activities/2d0 + &
                this%chemistry%waters(this%chemistry%bd_waters_indices(1))%activities/2d0
            this%target_waters(2)%aq_chem%volume=2d0*this%target_waters(2)%aq_chem%volume !> we set the volume of the first domain water to the double to account for the inflow boundary water
            !> Chapuza para agua que sale del dominio
            
            if (this%transport%BCs%labels(2) > 1) then !> flux BC outflow
                call this%chemistry%waters(this%chemistry%bd_waters_indices(2))%copy_aq_chem(&
                    aq_chem_exit) !> Copy exit water
                this%chemistry%waters(this%chemistry%bd_waters_indices(2))%concentrations=&
                    this%chemistry%waters(this%chemistry%bd_waters_indices(2))%concentrations/2d0 + &
                    this%target_waters(this%chemistry%num_target_waters)%aq_chem%concentrations/2d0
                this%chemistry%waters(this%chemistry%bd_waters_indices(2))%activities=&
                    this%chemistry%waters(this%chemistry%bd_waters_indices(2))%activities/2d0 + &
                    this%target_waters(this%chemistry%num_target_waters)%aq_chem%activities/2d0
                this%chemistry%waters(this%chemistry%bd_waters_indices(2))%volume=&
                    2d0*this%chemistry%waters(this%chemistry%bd_waters_indices(2))%volume !> we set the volume of the first domain water to the double to account for the inflow boundary water
                !> Chapuza para ultima target water
                !this%target_waters(this%chemistry%num_target_waters)%concentrations=&
                !    this%target_waters(this%chemistry%num_target_waters)%concentrations/2d0
                !this%target_waters(this%chemistry%num_target_waters)%activities=&
                !    this%target_waters(this%chemistry%num_target_waters)%activities/2d0
                this%target_waters(this%chemistry%num_target_waters)%aq_chem%volume=&
                    5d-1*this%target_waters(this%chemistry%num_target_waters)%aq_chem%volume !> we set the volume of the first domain water to the double to account for the inflow boundary water
                !    this%transport%spatial_discr%final_point) !> Assign unique ID
                !call this%chemistry%waters(this%chemistry%bd_waters_indices(2))%set_solid_chemistry(&
                !    this%chemistry%target_solids(this%chemistry%num_target_solids)) !> Assign unique ID
                !call this%chemistry%waters(this%chemistry%bd_waters_indices(2))%set_indices_rk() !> we set the indices of the kinetic reactions to the new outflow boundary water
            else !> Dirichlet outflow BC
                ! call this%target_waters(this%chemistry%num_target_waters)%aq_chem%copy_aq_chem(&
                !     aq_chem_exit) !> Copy exit water
            end if
            !else if (this%transport%BCs%labels(2) .eq. 1) then !> Dirichlet outflow BC
            !    call this%chemistry%waters(this%chemistry%bd_waters_indices(2))%copy_aq_chem(&
                    !this%target_waters(this%chemistry%num_target_waters)%aq_chem) !> Copy outflow BC particle
            !     call this%chemistry%waters(old_tar_wat_indices(this%chemistry%num_target_waters))%copy_aq_chem(&
            !         old_particles(old_bd_indices(2))) !> Copy old outflow BC particle
            !end if 
            !> Chapuza para ultimo nodo
            !this%target_waters(this%chemistry%num_target_waters)%aq_chem%concentrations=&
            !    this%target_waters(this%chemistry%num_target_waters)%aq_chem%concentrations*2d0
            !this%target_waters(this%chemistry%num_target_waters)%aq_chem%activities=&
            !    this%target_waters(this%chemistry%num_target_waters)%aq_chem%activities*2d0
            !this%target_waters(this%chemistry%num_target_waters)%aq_chem%volume=&
            !    this%target_waters(this%chemistry%num_target_waters)%aq_chem%volume/2d0 !> we set the volume of the first domain water to the double to account for the inflow boundary water
            !> we set the concentrations of the inflow boundary water
            !call disp_particles(1)%copy_aq_chem(this%chemistry%waters(old_bd_indices(1))) !> get displaced inflow boundary water
            ! do i=1,this%chemistry%num_target_waters-1
            !     call disp_particles(i)%copy_aq_chem(this%chemistry%waters(&
            !         old_tar_wat_indices(i))) !> get displaced particles
            ! end do
            !print *, this%chemistry%waters_init(3)%id
            !call this%chemistry%allocate_waters(this%chemistry%num_waters-1) !> allocate new particles array
            
            !this%chemistry%tar_wat_indices(1)=old_bd_indices(1) !> set the new domain target water index
            !this%chemistry%bd_waters_indices(2)=old_tar_wat_indices(this%chemistry%num_target_waters) !> set new outflow boundary water index
            !this%chemistry%bd_waters_indices(1)=0 !> no inflow boundary water
            ! print *, "New domain indices: ", this%chemistry%tar_wat_indices
            ! print *, "New boundary indices: ", this%chemistry%bd_waters_indices
            ! do i=1,old_bd_indices(2)-1
            !     call this%chemistry%waters(i)%copy_aq_chem(old_particles(i)) !> set old particles to the chemistry object
            !     !print *, this%chemistry%waters(i)%id !> set old particles to the chemistry object
            ! end do
            ! do i=1,this%chemistry%num_target_waters
            !     call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%copy_aq_chem(&
            !         disp_particles(i)) !> set displaced particles to the domain target waters
            ! end do
            ! call this%chemistry%waters(this%chemistry%bd_waters_indices(2))%copy_aq_chem(&
            !     disp_particles(old_tar_wat_indices(this%chemistry%num_target_waters))) !> set boundary outflow water
            ! print *, "New domain indices: ", this%chemistry%tar_wat_indices
            ! !print *, "New boundary indices: ", this%chemistry%bd_waters_indices
            ! print *, "Total number of particles: ", this%chemistry%num_target_waters
            !if (num_exit>0) then
            !    call this%chemistry%allocate_waters(this%chemistry%num_waters-num_exit) !> allocate new particles array
            !    call this%chemistry%allocate_tar_wat_indices(this%chemistry%num_target_waters-num_exit) !> allocate new indices array
            !    this%chemistry%tar_wat_indices=old_tar_wat_indices(1:size(old_tar_wat_indices)-num_exit) !> set new indices
            !    !call this%chemistry%allocate_bd_wat_indices(this%chemistry%num_target_waters-num_exit) !> allocate new indices array
            !    this%chemistry%bd_waters_indices(2)=this%chemistry%num_waters !> set new boundary water index
            !    this%chemistry%waters=old_particles(1:this%chemistry%num_waters) !> set new particles
            !    !this%chemistry%waters(this%chemistry%bd_waters_indices)=old_particles(old_bd_indices) !> set old particles to the chemistry object
            !    !this%chemistry%waters(this%chemistry%tar_wat_indices)=old_particles(old_tar_wat_indices(1:size(old_tar_wat_indices)-num_exit)) !> set old particles to the chemistry object
            !    print *, "Total number of particles: ", this%chemistry%num_waters
            !end if
            deallocate(new_coords) !> deallocate new coordinates array
            !deallocate(old_particles) !> deallocate old particles array
            deallocate(old_tar_wat_indices) !> deallocate old domain indices array
            deallocate(old_bd_indices) !> deallocate old boundary indices array
        end subroutine

        !> \brief Introduce new particles at inflow boundary
        !> \param[in,out] this Transient RT object
        !> \param[in]     k    Time step index
        !> \details
        !>   Adds new Lagrangian particle at inflow boundary each time step:
        !>   
        !>   Algorithm:
        !>   1. Save current particle array
        !>   2. Reallocate with num_particles + 1
        !>   3. Copy old particles to new array
        !>   4. Initialize new particle at inflow position
        !>   5. Set chemistry and solid interactions for new particle
        !>   6. Update boundary and domain indices
        !>   
        !>   New particle properties:
        !>   - Position: Inflow boundary location
        !>   - Chemistry: From initial/boundary conditions
        !>   - Volume: Default cell volume
        !>   - ID: Unique identifier (num_particles + k)
        !>   
        !>   Used in Lagrangian methods to maintain particle count as
        !>   particles exit the domain at the outflow.
        subroutine introduce_particle(this,k)
            !> Introduce new particle at the inflow boundary
            !> This is used in the transient case to introduce new particles at each time step.
            implicit none
            class(RT_1D_transient_c) :: this !> transient reactive transport object
            integer(kind=4), intent(in) :: k !> time step index

            type(target_water_c), allocatable :: old_particles(:) !> old particles
            integer(kind=4), allocatable :: old_tar_wat_indices(:) !> old particles indices
            integer(kind=4), allocatable :: old_bd_indices(:) !> old particles indices
            integer(kind=4) :: i !> loop index

            ! allocate(old_particles(this%chemistry%num_target_waters)) !> allocate old particles array
            ! do i=1,this%chemistry%num_target_waters
            !     call old_particles(i)%copy_tar_wat(this%target_waters(&
            !         i)) !> get old particles
            ! end do
            !old_particles=this%chemistry%waters !> get old particles
            !old_tar_wat_indices=this%chemistry%tar_wat_indices !> get old domain indices
            !old_bd_indices=this%chemistry%bd_waters_indices !> get old boundary indices
            !call this%chemistry%allocate_waters(this%chemistry%num_waters+1) !> Expand particle array by 1
            !call this%chemistry%allocate_tar_wat_indices(this%chemistry%num_target_waters+1) !> allocate new indices array
            !this%chemistry%bd_waters_indices=this%chemistry%bd_waters_indices_init !> Restore initial boundary indices
            ! this%chemistry%tar_wat_indices(1)=this%chemistry%bd_waters_indices(1)+1 !> set new domain target water index
            ! this%chemistry%tar_wat_indices(2:)=old_tar_wat_indices+1 !> Restore initial domain indices
            !this%chemistry%tar_wat_indices(1:this%chemistry%num_target_waters-1)=old_tar_wat_indices !> set new indices
            !this%chemistry%tar_wat_indices(this%chemistry%num_target_waters)=this%chemistry%num_waters-1 !> set the last domain target water index
            ! do i=2,this%chemistry%num_target_waters
            !     call this%target_waters(i)%copy_tar_wat(&
            !         this%target_waters(i-1)) !> Copy particles
            ! end do
            ! if (sum(this%chemistry%bd_waters_indices)>this%chemistry%CV_params%zero) then
            call this%target_waters(1)%aq_chem%copy_aq_chem_Lagr(&
                this%chemistry%waters(this%chemistry%bd_waters_indices(1))) !> Copy inflow boundary particle
            if (this%transport%spatial_discr%targets_flag.eq.1) then !> edge-centered spatial discretisation
            !    this%target_waters(1)%aq_chem%concentrations=&
            !        this%target_waters(1)%aq_chem%concentrations*2d0
            !    this%target_waters(1)%aq_chem%activities=&
            !        this%target_waters(1)%aq_chem%activities*2d0
                this%target_waters(1)%aq_chem%volume=&
                    this%target_waters(1)%aq_chem%volume/2d0
            end if
            ! else
            !     call this%chemistry%waters(this%chemistry%tar_wat_indices(1))%copy_aq_chem(&
            !         this%chemistry%waters_init(this%chemistry%tar_wat_indices_init(1)))
            !     if (this%transport%BCs%labels(2).eq.1) then !> Dirichlet outflow BC
            !         call this%chemistry%waters(this%chemistry%tar_wat_indices(&
            !             this%chemistry%num_target_waters))%copy_aq_chem(&
            !             this%chemistry%waters_init(this%chemistry%tar_wat_indices_init(&
            !             this%chemistry%num_target_waters_init)))
            !         !this%chemistry%waters(this%chemistry%tar_wat_indices(1))%concentrations=&
            !         !    old_particles(old_tar_wat_indices(1))%concentrations*1.1d0 !> increase by 10%
            !     end if
            ! end if
            !> Copy old inflow BC particle
            call this%target_waters(1)%set_id(&
                this%chemistry%num_target_waters+k) !> Assign unique ID
            ! call this%target_waters(1)%set_pos(&
            !     this%transport%spatial_discr%targets(1+this%transport%spatial_discr%targets_flag)%coord) !> Set position of new particle
            ! call this%target_waters(1)%aq_chem%set_solid_chemistry(&
            !     this%chemistry%target_solids(1+this%transport%spatial_discr%targets_flag)) !> Set solid chemistry for new particle
            ! !call this%target_waters(this%chemistry%tar_wat_indices(1))%aq_chem%set_indices_rk() !> Set kinetic reaction indices for new particle
            ! this%target_waters(1)%aq_chem%indices_rk=&
            !     this%target_waters_init(1)%aq_chem%indices_rk !> Copy primary species indices
            ! this%target_waters(1)%aq_chem%ind_prim_species=&
            !     this%target_waters_init(1)%aq_chem%ind_prim_species !> Copy primary species indices
            ! this%target_waters(1)%aq_chem%ind_sec_species=&
            !     this%target_waters_init(1)%aq_chem%ind_sec_species !> Copy secondary species indices
            !> Set chemistry for outflow boundary particle (chapuza)
            ! if (this%transport%BCs%labels(2) .eq. 1) then !> Dirichlet outflow BC
            !     !print *, this%chemistry%waters_init(1)%solid_chemistry%tar%id
            !     !print *, this%chemistry%waters_init(2)%solid_chemistry%tar%id
            !     !print *, this%chemistry%waters_init(3)%solid_chemistry%tar%id
            !     call this%chemistry%waters(this%chemistry%bd_waters_indices(2))%copy_aq_chem(&
            !         this%chemistry%waters_init(this%chemistry%bd_waters_indices_init(2))) !> Copy initial outflow BC particle
            !     call this%chemistry%waters(this%chemistry%bd_waters_indices(2))%set_id(&
            !         old_particles(size(old_particles))%id) !> Assign unique ID
            ! else
            !     call this%chemistry%waters(this%chemistry%bd_waters_indices(2))%copy_aq_chem(&
            !         old_particles(old_bd_indices(2))) !> Copy old outflow BC particle
            ! end if 
            
            !call this%chemistry%waters(this%chemistry%bd_waters_indices_init(1))%set_volume() !> Set default volume for new particle
            !call this%chemistry%waters(this%chemistry%bd_waters_indices_init(1))%set_id(this%chemistry%num_waters+k) !> Assign unique ID
            ! print *, "New particle introduced at time step ", k
            ! print *, "Particle ID: ", this%target_waters(1)%id
            !print *, "Total number of particles: ", this%chemistry%num_target_waters
            !print *, "New domain indices: ", this%chemistry%tar_wat_indices
            !print *, "New boundary indices: ", this%chemistry%bd_waters_indices
            !deallocate(old_particles) !> deallocate old particles array
            !deallocate(old_tar_wat_indices) !> deallocate old domain indices array
            !deallocate(old_bd_indices) !> deallocate old boundary indices array
        end subroutine
        
        !> \brief Introduce new particles at inflow boundary
        !> \param[in,out] this Transient RT object
        !> \param[in]     k    Time step index
        !> \details
        !>   Adds new Lagrangian particle at inflow boundary each time step:
        !>   
        !>   Algorithm:
        !>   1. Save current particle array
        !>   2. Reallocate with num_particles + 1
        !>   3. Copy old particles to new array
        !>   4. Initialize new particle at inflow position
        !>   5. Set chemistry and solid interactions for new particle
        !>   6. Update boundary and domain indices
        !>   
        !>   New particle properties:
        !>   - Position: Inflow boundary location
        !>   - Chemistry: From initial/boundary conditions
        !>   - Volume: Default cell volume
        !>   - ID: Unique identifier (num_particles + k)
        !>   
        !>   Used in Lagrangian methods to maintain particle count as
        !>   particles exit the domain at the outflow.
        subroutine introduce_particles_2D(this,k)
            !> Introduce new target waters at the leftmost cells (column 1) of each row.
            !> Each new target water receives the chemistry of the corresponding left
            !> boundary water: bd_waters_indices(i) -> target_waters(1+(i-1)*Nx), i=1..Ny.
            implicit none
            class(RT_2D_transient_c) :: this !> transient reactive transport object
            integer(kind=4), intent(in) :: k !> time step index

            integer(kind=4) :: i      !> row loop index
            integer(kind=4) :: Nx     !> number of cells in x-direction (columns)
            integer(kind=4) :: Ny     !> number of cells in y-direction (rows)
            integer(kind=4) :: idx    !> linear index of leftmost cell in row i
            integer(kind=4) :: bd_idx !> reversed boundary water index

            Nx = this%transport%spatial_discr%get_num_cells(1)
            Ny = this%transport%spatial_discr%get_num_cells(2)

            do i = 1, this%chemistry%num_left_bd_waters
                idx = 1 + (i - 1) * Nx !> leftmost cell of row i
                bd_idx = this%chemistry%num_left_bd_waters - i + 1
                !> Copy chemistry from left boundary water into target_water at column 1
                call this%target_waters(idx)%aq_chem%copy_aq_chem_Lagr( &
                    this%chemistry%waters(this%chemistry%bd_waters_indices(bd_idx)))
                if (this%transport%spatial_discr%targets_flag .eq. 1) then !> edge-centered
                    this%target_waters(idx)%aq_chem%volume = &
                        this%target_waters(idx)%aq_chem%volume / 2d0
                end if
                !> Assign unique ID
                call this%target_waters(idx)%set_id( &
                    this%chemistry%num_target_waters + &
                    (k - 1) * this%transport%spatial_discr%get_num_cells(2) + i)
                ! print *, "  Introduce row", i, ": bd_water(", &
                !     this%chemistry%bd_waters_indices(bd_idx), &
                !     ") -> target_water(", idx, "), new id =", &
                !     this%target_waters(idx)%id, &
                !     " solid_chem tar_id=", &
                !     this%target_waters(idx)%aq_chem%solid_chemistry%tar%id
            end do
        end subroutine

        subroutine set_target_waters(this,mesh,Lagr_flag,BCs,dir,root)
            implicit none
            class(RT_c) :: this !> chemistry object
            class(spatial_discr_c), intent(in) :: mesh !> mesh object ( we assume it has targets defined )
            logical, intent(in) :: Lagr_flag !> flag indicating if Lagrangian method is used
            integer(kind=4), intent(in) :: BCs(:) !> flag indicating if Dirichlet BCs are used
            character(len=*), intent(in) :: dir !> directory for output file
            character(len=*), intent(in) :: root !> root name for output file
            integer(kind=4) :: i !> loop index for domain waters
            real(kind=8), allocatable :: pos_inf(:) !> spatial position of the inflow boundary target water
            real(kind=8), allocatable :: pos_out(:) !> spatial position of the outflow boundary target water
            integer(kind=4), allocatable :: aux_tar_wat_ind(:) !> spatial position of the outflow boundary target water
            type(aqueous_chemistry_c), allocatable :: aux_wat(:) !> auxiliary water array for Lagrangian method
            type(aqueous_chemistry_c), allocatable :: aux_wat_init(:) !> auxiliary water array for Lagrangian method
            !> Boundary waters
            !allocate(pos_inf(1)) !< Allocate 1D array for inflow position coordinate
            !allocate(pos_out(1)) !< Allocate 1D array for outflow position coordinate
            if (mesh%Num_targets_defined.eqv..false.) then !< Validate that this%transport%spatial_discr has targets defined before attempting to access them
                error stop "mesh does not have targets defined" !< Abort execution with error message if targets not defined
            !end if
            !select type (mesh)
            !    class is (mesh!else if (mesh%targets_flag==0 .and. mesh%dim==1) then !< Check if mesh uses cell-centered targets (flag=0)
            !    !pos_inf(1)=mesh%init_point(1)-0.5d0*mesh%get_cell_size(1) !< Position inflow water half a cell before first cell center (upstream boundary)
            !    !pos_out(1)=mesh%final_point(1)+0.5d0*mesh%get_cell_size(mesh%num_targets) !< Position outflow water half a cell after last cell center (downstream boundary)
            else if (mesh%targets_flag==0 .and. mesh%dim==2) then
                this%chemistry%num_left_bd_waters=mesh%get_num_cells(2)
                !this%chemistry%num_right_bd_waters=this%chemistry%num_left_bd_waters
                !pos_inf(1)=mesh%init_point(1)-0.5d0*mesh%get_cell_size(1) !< Position inflow water half a cell before first cell center (upstream boundary)
                !pos_out(1)=mesh%final_point(1)+0.5d0*mesh%get_cell_size(mesh%num_targets) !< Position outflow water half a cell after last cell center (downstream boundary)
            else if (mesh%targets_flag==0 .and. mesh%dim==1) then
                this%chemistry%num_left_bd_waters=1
                !this%chemistry%num_right_bd_waters=this%chemistry%num_left_bd_waters
                !pos_inf(1)=mesh%init_point(1)-0.5d0*mesh%get_cell_size(1) !< Position inflow water half a cell before first cell center (upstream boundary)
                !pos_out(1)=mesh%final_point(1)+0.5d0*mesh%get_cell_size(mesh%num_targets) !< Position outflow water half a cell after last cell center (downstream boundary)
            else if (mesh%targets_flag==1 .and. mesh%dim==1) then !< Check if mesh uses edge-centered targets (flag=1)
                !if (Lagr_flag) then
                !    this%chemistry%num_target_waters=this%chemistry%num_target_waters+&
                !        this%chemistry%num_bd_waters-1 !< Update number of target waters based on mesh targets (excluding boundary edges)
                !else
                this%chemistry%num_left_bd_waters=1
                !this%chemistry%num_right_bd_waters=this%chemistry%num_left_bd_waters
                if (BCs(1) .eq. 1 .and. BCs(2) .eq. 1) then !< Check if Dirichlet BC is used
                     this%chemistry%num_target_waters=this%chemistry%num_target_waters+&
                        this%chemistry%num_bd_waters !< Update number of target waters based on mesh targets (including boundary edges)
                !end if
                    aux_tar_wat_ind=this%chemistry%tar_wat_indices !< Store current number of target waters before adding boundary waters
                    call this%chemistry%allocate_tar_wat_indices() !< Reallocate target water indices array to accommodate new boundary waters
                    !do i=1,this%chemistry%num_bd_waters !< Loop over all target waters to shift existing indices
                    !if (Lagr_flag .eqv. .false.) then
                    !    this%chemistry%tar_wat_indices(1)=this%chemistry%bd_waters_indices(1) !< Assign new index for boundary water
                    !    !this%chemistry%tar_wat_indices(&
                    !    !    this%chemistry%num_target_waters)=this%chemistry%bd_waters_indices(this%chemistry%num_bd_waters) !< Assign new index for boundary water
                    !    do i=1,this%chemistry%num_target_waters-this%chemistry%num_bd_waters !< Loop over all target waters to update indices
                    !        this%chemistry%tar_wat_indices(1+i)=&
                    !            aux_tar_wat_ind(i) !< Shift existing interior target water indices by number of boundary waters
                    !    end do
                    !    this%chemistry%tar_wat_indices(&
                    !        this%chemistry%num_target_waters)=this%chemistry%bd_waters_indices(this%chemistry%num_bd_waters) !< Assign new index for boundary water
                    !else
                        allocate(aux_wat(this%chemistry%num_waters)) !< Allocate auxiliary water array to temporarily hold water data
                        allocate(aux_wat_init(this%chemistry%num_waters)) !< Allocate auxiliary water array to temporarily hold water data
                        do i=1,this%chemistry%num_waters !< Loop over all waters to find boundary water index
                            call aux_wat(i)%copy_aq_chem(this%chemistry%waters(i))
                            call aux_wat_init(i)%copy_aq_chem(this%chemistry%waters_init(i))
                        end do
                        call this%chemistry%allocate_waters(this%chemistry%num_waters+2) !< Expand water array by 1 to accommodate inflow water
                        call this%chemistry%allocate_waters_init(this%chemistry%num_waters) !< Expand water array by 1 to accommodate inflow water
                        call this%chemistry%waters(1)%copy_aq_chem(aux_wat(this%chemistry%bd_waters_indices(1))) !< Assign new index for inflow water
                        call this%chemistry%waters(1)%set_id(this%chemistry%num_waters-1) !< Assign new index for inflow water
                        call this%chemistry%waters_init(1)%copy_aq_chem(aux_wat_init(this%chemistry%bd_waters_indices(1))) !< Assign new index for inflow water
                        call this%chemistry%waters(this%chemistry%num_waters)%copy_aq_chem(&
                            aux_wat(this%chemistry%bd_waters_indices(2))) !< Assign new index for inflow water
                        call this%chemistry%waters(this%chemistry%num_waters)%set_id(&
                            this%chemistry%num_waters) !< Assign new index for inflow water
                        call this%chemistry%waters_init(this%chemistry%num_waters)%copy_aq_chem(&
                            aux_wat_init(this%chemistry%bd_waters_indices(2))) !< Assign new index for inflow water
                        !call this%chemistry%waters_init()%copy_aq_chem(aux_wat_init(this%chemistry%bd_waters_indices(1))) !< Assign new index for inflow water
                        !this%chemistry%waters(1)%concentrations=&
                        !    this%chemistry%waters(1)%concentrations/2d0 !< Adjust concentrations for new inflow water
                        !this%chemistry%waters(1)%activities=&
                        !    this%chemistry%waters(1)%activities/2d0 !< Adjust activities for new inflow water
                        do i=2,this%chemistry%num_waters-1 !< Loop over all waters to copy back data
                            call this%chemistry%waters(i)%copy_aq_chem(aux_wat(i-1))
                            call this%chemistry%waters_init(i)%copy_aq_chem(aux_wat_init(i-1))
                        end do
                        !do i=1,this%chemistry%num_target_waters !< Loop over all target waters to shift existing indices
                        !    this%chemistry%tar_wat_indices(i+1)=&
                        !        aux_tar_wat_ind(i)+1 !< Shift existing interior target water indices by number of boundary waters
                        !end do
                        do i=1,this%chemistry%num_target_waters-2 !< Loop over all target waters to update indices
                            this%chemistry%tar_wat_indices(i)=&
                                aux_tar_wat_ind(i) !< Shift existing interior target water indices by number of boundary waters
                        end do
                        this%chemistry%tar_wat_indices(&
                            this%chemistry%num_target_waters-1)=&
                            this%chemistry%bd_waters_indices(this%chemistry%num_bd_waters)
                        this%chemistry%bd_waters_indices(this%chemistry%num_bd_waters)=&
                            this%chemistry%num_waters !< Update outflow boundary water index
                        this%chemistry%tar_wat_indices(&
                            this%chemistry%num_target_waters)=&
                            this%chemistry%bd_waters_indices(this%chemistry%num_bd_waters)-1 !< Assign new index for boundary water
                        print *, "Updated target water indices: ", this%chemistry%tar_wat_indices
                        print *, "Updated boundary water indices: ", this%chemistry%bd_waters_indices
                        deallocate(aux_wat) !< Deallocate auxiliary water array after use
                        deallocate(aux_tar_wat_ind) !< Deallocate auxiliary target water indices array after use
                    !end if
                end if
                    !if (Lagr_flag) then
                    !    do i=1,this%chemistry%num_target_waters-this%chemistry%num_bd_waters+1 !< Loop over all target waters to update indices
                    !        this%chemistry%tar_wat_indices(i)=&
                    !            aux_tar_wat_ind(i) !< Shift existing interior target water indices by number of boundary waters
                    !    end do
                    !else
                    !    do i=1,this%chemistry%num_target_waters-this%chemistry%num_bd_waters !< Loop over all target waters to update indices
                    !        this%chemistry%tar_wat_indices(1+i)=&
                    !            aux_tar_wat_ind(i) !< Shift existing interior target water indices by number of boundary waters
                    !    end do
                    !end if
                    !if (i<=this%chemistry%num_bd_waters) then !< Check if current index corresponds to a boundary water
                    !    this%chemistry%tar_wat_indices(i)=i !< Assign new index for boundary water
                    !else !< Current index corresponds to an interior target water
                    !    this%chemistry%tar_wat_indices(i)=aux_tar_wat_ind(i-this%chemistry%num_bd_waters) + this%chemistry%num_bd_waters !< Shift existing interior target water indices by number of boundary waters
                    !end if !< End boundary vs interior water index branching
                !end do !< End loop over target waters to update indices
                this%chemistry%waters(this%chemistry%tar_wat_indices(&
                    1))%volume=&
                    this%chemistry%waters(this%chemistry%tar_wat_indices(&
                    2))%volume/2d0
                this%chemistry%waters(this%chemistry%tar_wat_indices(&
                    this%chemistry%num_target_waters))%volume=&
                    this%chemistry%waters(this%chemistry%tar_wat_indices(&
                    this%chemistry%num_target_waters-1))%volume/2d0
                !pos_inf(1)=mesh%init_point(1) !< Position inflow water at domain start edge
                !pos_out(1)=mesh%final_point(1) !< Position outflow water at domain end edge
            else !< mesh targets_flag has invalid value (not 0 or 1)
                error stop "mesh targets_flag not recognized" !< Abort execution with error message if flag value is invalid
            end if !< End mesh targets_flag branching
            !allocate(this%target_waters(this%chemistry%num_target_waters)) !< Allocate target waters array based on number of target waters in chemistry object
            !select type (mesh)
            !class is (mesh_2D_Euler_homog_c)
            !    this%chemistry%num_left_bd_waters=mesh%Num_cells_y
            this%chemistry%num_right_bd_waters=this%chemistry%num_left_bd_waters
            
            
            
            call this%chemistry%allocate_tar_wat_indices_init(this%chemistry%num_target_waters) !< Allocate initial state target waters array based on number of target waters in chemistry object
            this%chemistry%tar_wat_indices_init=this%chemistry%tar_wat_indices !< Store initial target water indices for reference
            allocate(this%target_waters(this%chemistry%num_target_waters)) !< Allocate target waters array based on number of target waters in chemistry object
            allocate(this%target_waters_init(this%chemistry%num_target_waters)) !< Allocate target waters array based on number of target waters in chemistry object
            !call this%chemistry%allocate_waters_init(this%chemistry%num_waters) !< Allocate initial state boundary waters array based on number of boundary waters in chemistry object
            
            !if (this%num_bd_waters==2) then !< Check if system has two boundary waters (standard inflow/outflow configuration)
                ! call this%waters(this%bd_waters_indices(1))%set_pos(pos_inf) !< Assign inflow position to first boundary target water (current state)
                ! call this%waters_init(this%bd_waters_indices(1))%set_pos(pos_inf) !< Assign inflow position to first boundary target water (initial state)
                ! call this%waters(this%bd_waters_indices(2))%set_pos(pos_out) !< Assign outflow position to second boundary target water (current state)
                ! call this%waters_init(this%bd_waters_indices(2))%set_pos(pos_out) !< Assign outflow position to second boundary target water (initial state)
                !if (mesh%targets_flag==0 .or. this%chemistry%num_waters==this%chemistry%num_target_waters) then !< Check if mesh uses cell-centered targets (flag=0)
                    !> Target waters
                    !print *, "[DEBUG set_target_waters] num_target_waters=", this%chemistry%num_target_waters
                    !print *, "[DEBUG set_target_waters] mesh%Num_targets=", mesh%Num_targets
                    !print *, "[DEBUG set_target_waters] allocated(mesh%targets)=", allocated(mesh%targets)
                    !if (allocated(mesh%targets)) print *, "[DEBUG set_target_waters] size(mesh%targets)=", size(mesh%targets)
                    do i=1,this%chemistry%num_target_waters !< Loop over all target waters (interior spatial points)
                        !print *, "[DEBUG set_target_waters] i=", i, " allocated(coord)=", allocated(mesh%targets(i)%coord)
                        if (allocated(mesh%targets(i)%coord)) then
                            !print *, "[DEBUG set_target_waters] coord=", mesh%targets(i)%coord
                        else
                            print *, "[DEBUG set_target_waters] ERROR: coord not allocated for target ", i
                            cycle
                        end if
                        !call this%target_waters(this%chemistry%tar_wat_indices(i))%set_id(i) !< Assign spatial coordinates from mesh target to domain target water i (current state), offset by +1 to account for edge/cell centering
                        call this%target_waters(i)%set_pos(mesh%targets(i)%coord) !< Assign spatial coordinates from mesh target to domain target water i (current state), offset by +1 to account for edge/cell centering
                        !call this%target_waters_init(i)%set_pos(mesh%targets(i)%coord) !< Assign spatial coordinates from mesh target to domain target water i (initial state), ensuring consistent positioning
                    end do !< End loop over target waters
                !else if (mesh%targets_flag==1) then !< Check if mesh uses edge-centered targets (flag=1)
                !    !> Target waters
                !    do i=1,this%chemistry%num_target_waters !< Loop over all target waters (interior spatial points)
                !        !call this%target_waters(this%chemistry%tar_wat_indices(i))%set_id(i) !< Assign spatial coordinates from mesh target to domain target water i (current state), offset by +1 to account for edge/cell centering
                !        call this%target_waters(i)%set_pos(mesh%targets(i+1)%coord) !< Assign spatial coordinates from mesh target to domain target water i (current state)
                !        !call this%target_waters_init(i)%set_pos(mesh%targets(i+1)%coord) !< Assign spatial coordinates from mesh target to domain target water i (initial state), ensuring consistent positioning
                !    end do !< End loop over target waters
                !else !< mesh targets_flag has invalid value (not 0 or 1)
                !    error stop "mesh targets_flag not recognized" !< Abort execution with error message if flag value is invalid
                !end if !< End mesh targets_flag branching within two boundary waters case
                ! !> Target waters
                ! do i=1,this%num_target_waters !< Loop over all target waters (interior spatial points)
                !     call this%waters(this%tar_wat_indices(i))%set_pos(mesh%targets(i+mesh%targets_flag)%coord) !< Assign spatial coordinates from mesh target to domain target water i (current state), offset by targets_flag to account for edge/cell centering
                !     call this%waters_init(this%tar_wat_indices(i))%set_pos(mesh%targets(i+mesh%targets_flag)%coord) !< Assign spatial coordinates from mesh target to domain target water i (initial state), ensuring consistent positioning
                ! end do !< End loop over target waters
            ! else if (this%num_bd_waters==0) then !< Check if system has zero boundary waters (closed system configuration)
            !     !call this%waters(this%bd_waters_indices(1))%set_pos(pos_inf) !< Assign inflow position to the single boundary target water (current state)
            !     !call this%waters_init(this%bd_waters_indices(1))%set_pos(pos_inf) !< Assign inflow position to the single boundary target water (initial state)
            !     !> Target waters
            !     do i=1,this%num_target_waters !< Loop over all target waters (interior spatial points)
            !         call this%waters(this%tar_wat_indices(i))%set_pos(mesh%targets(i)%coord) !< Assign spatial coordinates from mesh target to domain target water i (current state), offset by targets_flag
            !         call this%waters_init(this%tar_wat_indices(i))%set_pos(mesh%targets(i)%coord) !< Assign spatial coordinates from mesh target to domain target water i (initial state)
            !     end do !< End loop over target waters
            ! else !< Unsupported number of boundary waters configuration
            !     error stop "Number of boundary waters not supported in set_target_waters" !< Abort execution with error message for unsupported configuration
            ! end if !< End boundary waters configuration branching
                !if (all(this%chemistry%bd_waters_indices==0)) then !< Check if inflow boundary water index is valid (greater than zero)
                !    this%chemistry%bd_waters_indices(1)=1 !< Set inflow boundary water index to first water index
                !    !this%chemistry%bd_waters_indices_init(1)=1 !< Set inflow boundary water index to first water index
                !    this%chemistry%bd_waters_indices(2)=this%chemistry%num_waters !< Set outflow boundary water index to last water index
                !    !this%chemistry%bd_waters_indices_init(2)=this%chemistry%num_waters !< Set outflow boundary water index to last water index
                !end if
                    do i=1,this%chemistry%num_target_waters !< Loop over all target waters (interior spatial points)
                        call this%target_waters(i)%set_id(i) !< Assign spatial coordinates from mesh target to domain target water i (initial state), ensuring consistent positioning
                        call this%target_waters_init(i)%set_id(&
                            this%target_waters(i)%id) !< Assign spatial coordinates from mesh target to domain target water i (initial state), ensuring consistent positioning
                        ! call this%target_waters_old(i)%set_id(&
                        !     this%target_waters(i)%id) !< Assign spatial coordinates from mesh target to domain target water i (current state)
                        call this%target_waters(i)%set_aq_chem(this%chemistry%waters(&
                            this%chemistry%tar_wat_indices(i))) !< Assign spatial coordinates from mesh target to domain target water i (initial state), ensuring consistent positioning
                        call this%target_waters_init(i)%set_aq_chem(this%chemistry%waters_init(&
                            this%chemistry%tar_wat_indices_init(i))) !< Assign spatial coordinates from mesh target to domain target water i (initial state), ensuring consistent positioning
                        ! call this%target_waters_old(i)%set_aq_chem(this%chemistry%waters(&
                        !     this%chemistry%tar_wat_indices(i))) !< Assign spatial coordinates from mesh target to domain target water i (current state)
                        call this%target_waters_init(i)%set_pos(this%target_waters(i)%pos) !< Assign spatial coordinates from mesh target to domain target water i (initial state), ensuring consistent positioning
                        ! call this%target_waters_old(i)%set_pos(this%target_waters(i)%pos) !< Assign spatial coordinates from mesh target to domain target water i (initial state), ensuring consistent positioning
                    end do !< End loop over target waters
                    !deallocate(pos_inf) !< Deallocate inflow position array after use
                    !deallocate(pos_out) !< Deallocate outflow position array after use
                    !> Write initial aqueous concentrations to file (rows=targets, columns=species in natural order)
                    ! block
                    !     integer :: iunit, it, js
                    !     open(newunit=iunit, file=trim(dir)//trim(root)//'_conc_aq_init.dat', status='replace', action='write')
                    !     do it = 1, this%chemistry%num_target_waters
                    !         write(iunit, '(*(ES15.5))') &
                    !             (this%target_waters(it)%aq_chem%concentrations(js), &
                    !              js = 1, size(this%target_waters(it)%aq_chem%concentrations))
                    !     end do
                    !     close(iunit)
                    ! end block
                    
        end subroutine
        
        !> \brief Remap transport-local water indices to global chemistry water indices
        !> \details
        !>   During compute_mixing_ratios, mix_conc_indices and mix_react_indices store
        !>   transport-local indices following the convention:
        !>     [1=inflow_bd, 2..N+1=targets, N+2=outflow_bd]
        !>   But the actual indices into chemistry%waters(:) come from tar_wat_indices
        !>   and bd_waters_indices, which are read from the input file.
        !>   This subroutine remaps the transport-local indices to the actual global
        !>   water indices after chemistry has been initialized.
        subroutine remap_mix_indices_2D(this)
            use spatial_discr_2D_m, only: mesh_2D_Euler_homog_c
            implicit none
            class(RT_2D_transient_c) :: this
            integer(kind=4) :: i, j, v, k
            integer(kind=4) :: num_tar, num_inf, num_out
            integer(kind=4) :: Nx, Ny, row_i
            integer(kind=4) :: num_left_bd, num_right_bd
            
            num_tar = this%transport%spatial_discr%Num_targets
            num_inf = 1  !< num_inf_ext_mix_rat used in compute_mixing_ratios
            num_out = 1  !< num_out_ext_mix_rat used in compute_mixing_ratios
            num_left_bd = this%chemistry%num_left_bd_waters
            num_right_bd = this%chemistry%num_right_bd_waters
            
            !> Get grid dimensions
            select type (mesh => this%transport%spatial_discr)
            type is (mesh_2D_Euler_homog_c)
                Nx = mesh%Num_cells_x
                Ny = mesh%Num_cells_y
            class default
                error stop "remap_mix_indices: 2D homogeneous mesh expected"
            end select
            
            write(*,*) ''
            write(*,*) '================================================================================'
            write(*,*) 'DEBUG [remap_mix_indices]: ENTRY'
            write(*,*) '================================================================================'
            write(*,*) '  num_tar =', num_tar, '  Nx =', Nx, '  Ny =', Ny
            write(*,*) '  num_inf =', num_inf, '  num_out =', num_out
            write(*,*) '  num_bd_waters       =', this%chemistry%num_bd_waters
            write(*,*) '  num_left_bd_waters  =', num_left_bd
            write(*,*) '  num_right_bd_waters =', num_right_bd
            !write(*,*) '  bd_waters_indices   =', this%chemistry%bd_waters_indices
            !write(*,*) '  tar_wat_indices     =', this%chemistry%tar_wat_indices
            
            !> Remap mix_conc_indices
            write(*,*) ''
            write(*,*) '--- mix_conc_indices BEFORE remap ---'
            ! do i = 1, this%transport%mix_conc_indices%num_cols
            !     write(*,'(A,I4,A,*(I6))') '  col(', i, ') = ', &
            !         this%transport%mix_conc_indices%cols(i)%col_1
            ! end do
            
            do i = 1, this%transport%mix_conc_indices%num_cols
                row_i = (i - 1) / Nx + 1
                do j = 1, this%transport%mix_conc_indices%cols(i)%dim - 2
                    v = this%transport%mix_conc_indices%cols(i)%col_1(j)
                    k = v - num_inf
                    if (k >= 1 .and. k <= num_tar) then
                        !> Domain target water k -> tar_wat_indices(k)
                        this%transport%mix_conc_indices%cols(i)%col_1(j) = &
                            this%chemistry%tar_wat_indices(k)
                    else if (v == num_inf) then
                        !> Inflow (left) boundary for this row
                        !> bd_waters_indices ordered top-to-bottom: index 1=row Ny, index Ny=row 1
                        this%transport%mix_conc_indices%cols(i)%col_1(j) = &
                            this%chemistry%bd_waters_indices( &
                            num_left_bd - row_i + 1)
                    else if (v == num_tar + num_inf + num_out) then
                        !> Outflow (right) boundary for this row
                        !> Right boundaries follow left in same top-to-bottom order
                        this%transport%mix_conc_indices%cols(i)%col_1(j) = &
                            this%chemistry%bd_waters_indices( &
                            num_left_bd + num_right_bd - row_i + 1)
                    end if
                end do
            end do
            
            write(*,*) ''
            write(*,*) '--- mix_conc_indices AFTER remap ---'
            ! do i = 1, this%transport%mix_conc_indices%num_cols
            !     write(*,'(A,I4,A,*(I6))') '  col(', i, ') = ', &
            !         this%transport%mix_conc_indices%cols(i)%col_1
            ! end do
            
            !> Remap mix_react_indices (typically only domain water indices)
            write(*,*) ''
            write(*,*) '--- mix_react_indices BEFORE remap ---'
            ! do i = 1, this%transport%mix_react_indices%num_cols
            !     write(*,'(A,I4,A,*(I6))') '  col(', i, ') = ', &
            !         this%transport%mix_react_indices%cols(i)%col_1
            ! end do
            
            do i = 1, this%transport%mix_react_indices%num_cols
                row_i = (i - 1) / Nx + 1
                do j = 1, this%transport%mix_react_indices%cols(i)%dim - 2
                    v = this%transport%mix_react_indices%cols(i)%col_1(j)
                    k = v - num_inf
                    if (k >= 1 .and. k <= num_tar) then
                        this%transport%mix_react_indices%cols(i)%col_1(j) = &
                            this%chemistry%tar_wat_indices(k)
                    else if (v == num_inf) then
                        this%transport%mix_react_indices%cols(i)%col_1(j) = &
                            this%chemistry%bd_waters_indices( &
                            num_left_bd - row_i + 1)
                    else if (v == num_tar + num_inf + num_out) then
                        this%transport%mix_react_indices%cols(i)%col_1(j) = &
                            this%chemistry%bd_waters_indices( &
                            num_left_bd + num_right_bd - row_i + 1)
                    end if
                end do
            end do
            
            write(*,*) ''
            write(*,*) '--- mix_react_indices AFTER remap ---'
            ! do i = 1, this%transport%mix_react_indices%num_cols
            !     write(*,'(A,I4,A,*(I6))') '  col(', i, ') = ', &
            !         this%transport%mix_react_indices%cols(i)%col_1
            ! end do
            
            write(*,*) ''
            write(*,*) '================================================================================'
            write(*,*) 'DEBUG [remap_mix_indices]: EXIT'
            write(*,*) '================================================================================'
            write(*,*) ''
        end subroutine

        subroutine set_target_solids_time(this)
            implicit none
            class(RT_c) :: this !> chemistry object
            !real(kind=8), intent(in) :: Delta_t !> time step size for transient simulation
            !class(spatial_discr_c), intent(in) :: mesh !> mesh object ( we assume it has targets defined )
            integer(kind=4) :: i, j , Nx, Ny !> loop index for domain waters
            real(kind=8) :: Delta_t !> time step size for transient simulation
            !real(kind=8), allocatable :: pos_inf(:) !> spatial position of the inflow boundary target water
            !real(kind=8), allocatable :: pos_out(:) !> spatial position of the outflow boundary target water
            select type (this)
            class is (RT_2D_transient_c)
                if (this%transport%Lagr_flag) then
                    Delta_t=this%transport%time_discr%get_Delta_t() !< Update target water positions based on transport time discretization and time step size for Lagrangian method
                    Nx=this%transport%spatial_discr%get_num_cells(1)
                    Ny=this%transport%spatial_discr%get_num_cells(2)
                    do i=1,Ny !< Loop over all target waters to update positions
                        do j=1,Nx
                            call this%chemistry%target_solids((i-1)*Nx+j)%set_time(Delta_t*(2*j-1)/2d0) !< Update target solid arrival time assuming isochronal mesh
                            !print *, "Updated target solid ", (i-1)*Nx+j, " time to ", this%chemistry%target_solids((i-1)*Nx+j)%time
                        end do
                    end do
                    !call this%transport%spatial_discr%update_target_positions(this%target_waters, Delta_t) !< Update target water positions based on transport spatial discretization and time step size
                end if
            end select
                !call this%transport%spatial_discr%update_target_positions(this%target_waters, Delta_t) !< Update target water positions based on transport spatial discretization and time step size
            ! if (mesh%Num_targets_defined.eqv..false.) then !< Validate that this%transport%spatial_discr has targets defined before attempting to access them
            !     error stop "mesh does not have targets defined" !< Abort execution with error message if targets not defined
            ! end if
            ! allocate(pos_inf(1)) !< Allocate 1D array for inflow position coordinate
            ! allocate(pos_out(1)) !< Allocate 1D array for outflow position coordinate
            ! pos_inf(1)=mesh%init_point(1)-0.5d0*mesh%get_cell_size(1) !< Position inflow water half a cell before first cell center (upstream boundary)
            ! pos_out(1)=mesh%final_point(1)+0.5d0*mesh%get_cell_size(mesh%num_targets) !< Position outflow water half a cell after last cell center (downstream boundary)
            ! do i=1,this%chemistry%num_bd_waters !< Loop over all boundary waters to set positions
            !     call this%chemistry%waters(this%chemistry%bd_waters_indices(i))%set_pos(pos_inf) !< Assign inflow position to boundary target water i (current state)
            !     call this%chemistry%waters_init(this%chemistry%bd_waters_indices(i))%set_pos(pos_inf) !< Assign inflow position to boundary target water i (initial state)
            ! end do
            ! deallocate(pos_inf) !< Deallocate inflow position array after use
            ! deallocate(pos_out) !< Deallocate outflow position array after use
        end subroutine
        
        !function get_transport_obj_2D_trans(this) result(transport_obj)
        !    implicit none
        !    class(RT_2D_transient_c) :: this !> transient reactive transport object
        !    class(transport_2D_Euler_c), pointer :: get_transport_obj

        subroutine solve_RT_ideal_lump_Lagr_trans_flux_2D(this, dir, root)
            class(RT_2D_transient_c) :: this
            character(len=*), intent(in) :: dir, root
            error stop "solve_RT_ideal_lump_Lagr_trans_flux not applicable for 2D"
        end subroutine

        subroutine solve_RT_ideal_lump_Euler_trans_flux_2D(this, dir, root)
            class(RT_2D_transient_c) :: this
            character(len=*), intent(in) :: dir, root
            error stop "solve_RT_ideal_lump_Euler_trans_flux not applicable for 2D"
        end subroutine

        subroutine solve_RT_ideal_cons_Lagr_stat_flux_2D(this, dir, root)
            class(RT_2D_transient_c) :: this
            character(len=*), intent(in) :: dir, root
            error stop "solve_RT_ideal_cons_Lagr_stat_flux not applicable for 2D"
        end subroutine

        subroutine solve_RT_ideal_cons_Euler_stat_flux_2D(this, dir, root)
            class(RT_2D_transient_c) :: this
            character(len=*), intent(in) :: dir, root
            error stop "solve_RT_ideal_cons_Euler_stat_flux not applicable for 2D"
        end subroutine

        subroutine solve_RT_ideal_lump_Euler_stat_flux_2D(this, dir, root)
            class(RT_2D_transient_c) :: this
            character(len=*), intent(in) :: dir, root
            error stop "solve_RT_ideal_lump_Euler_stat_flux not applicable for 2D"
        end subroutine

end module RT_m