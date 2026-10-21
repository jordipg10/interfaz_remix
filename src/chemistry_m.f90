!> @file chemistry_m.f90
!> @brief Main chemistry module for reactive transport modeling.
!> @details This module defines the top-level chemistry_c class, which coordinates every
!>          chemical process in a reactive transport simulation. It manages water types,
!>          target waters (the per-cell aqueous states), solid phases, gas phases, mineral
!>          zones and reactive zones, and owns the chemical system (species, reactions and
!>          thermodynamic data).
!>
!>          Responsibilities include:
!>          - Reading chemical data (CHEPROO / PHREEQC / PFLOTRAN input formats).
!>          - Driving the Water Mixing Approach (WMA) reactive-mixing solvers, in both the
!>            lumped and consistent variants and with explicit or implicit chemistry.
!>          - Managing equilibrium speciation and kinetic reaction integration.
!>          - Writing chemical output (see write_chemistry).
!> @author Jordi Petchamé-Guerrero
!> @date 2025
module chemistry_m
    use aqueous_chemistry_m, only: aqueous_chemistry_c, compute_c_mix, &
        reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2, compute_c_mix_by_zones_lump, & !< Import aqueous chemistry class
        zone_water_c !< current-step within-zone sub-mixture store
    use solid_chemistry_m, only: solid_chemistry_c !< Import solid chemistry class
    use gas_chemistry_m, only: gas_chemistry_c !< Import gas chemistry class
    use chem_system_m, only: chem_system_c !< Import chemical system and convergence parameters
    use reactive_zone_m, only: reactive_zone_c, compare_react_zones !< Import reactive zone classes
    use mineral_zone_m, only: mineral_zone_c !< Import mineral and mineral zone classes
    use chem_out_options_m, only: chem_out_options_s !< Import output control structure
    use spatial_discr_m, only: spatial_discr_c !< Import 1D spatial discretization
    use time_discr_m, only: time_discr_c !< Import time discretization classes
    use vectors_m, only: inf_norm_vec_int, inf_norm_vec_real !< Import vector norm functions
    use CV_params_m, only: CV_params_s !< Import convergence parameters structure
    use arrays_m, only: int_array_c, real_array_c !< Import sparse matrix class
    implicit none !< Require explicit declaration of all variables
    private !< Make all entities private by default
    public :: chemistry_c !< Make chemistry class public
    public :: write_out_step_vars !< Shared per-output-step writer honouring chem_out_options%vars
    public :: check_surf_compl_act_sum !< DEBUG: per-step/per-cell surface-complex activity sum check
    public :: reset_wat_indices_zones !< Rebuild wat_indices_zones when the reactive-zone set/membership changed
    public :: interfaz_comps_arch_eq_kin, interfaz_comps_arch_eq, interfaz_esp_arch !< External interface subroutines for reactive mixing iterations
    public :: interfaz_comps_arch_eq_kin_T, interfaz_comps_arch_eq_T !< Transposed-input variants (rows: targets, columns: components)
    public :: interfaz_comps_arch_eq_kin_mix !< Conservative + reactive mixing variant (reads concentrations prior to mixing and the mixing ratios)
    !> @brief Main chemistry class - central coordinator for all chemical processes
    !> @details This class manages the complete chemical system including:
    !>          - Water types and target waters at spatial locations
    !>          - Solid phases (minerals, surfaces, biofilms)
    !>          - Gas phases and reactive zones
    !>          - Chemical reaction networks and thermodynamic data
    !>          - Numerical algorithms for reactive mixing
    !> @note This is the primary interface for all chemical calculations in RT simulations
    type :: chemistry_c
        !> @name Numerical Algorithm Control Parameters
        !> @{
        integer(kind=4) :: read_opt                                        !< Chemical data input format: 1=CHEPROO-based, 2=PHREEQC, 3=PFLOTRAN
        integer(kind=4) :: act_coeffs_model                                !< Activity coefficient model: Debye-Hückel, Davies, Pitzer, etc.
        logical :: lump_flag                                             !< Spatial lumping flag: .true.=apply lumping, .false.=full resolution
        logical :: mix_by_zones_flag                             !< Mixing flag: .true.=mix per reactive zone, .false.=mix directly over all waters regardless of reactive zone
        real(kind=8) :: est_prm !< [-] estimation parameter for downstream kinetic reactions (0 = no extrapolation)
        integer(kind=4) :: r_down_opt                                   !< Downstream reaction rates estimation method: 1,2,3,4 (different interpolation schemes) (see documentation)
        integer(kind=4) :: rk_avg_opt                                          !< Reaction rate averaging method: 1=concentration-based, 2=direct rate averaging
        integer(kind=4) :: Jac_opt                                             !< Jacobian computation model: 0=incremental coefficients, 1=analytical
        !> @}
        
        !> @name Core Chemical System
        !> @{
        type(chem_system_c) :: chem_syst                                       !< Complete chemical system: species, reactions, thermodynamic data, matrices
        type(CV_params_s) :: CV_params                                         !< Convergence parameters for iterative solvers
        type(chem_out_options_s) :: chem_out_options                           !< Chemistry output options and control variables
        !> @}
        
        !> @name Water Type Management
        !> @{
        integer(kind=4) :: num_wat_types=0                                     !< Total number of distinct water composition types (initial + boundary conditions)
        integer(kind=4) :: num_waters=0                            !< Number of aqueous components in chemical system
        integer(kind=4) :: num_target_waters=0                                 !< Current number of target water locations (evolves during simulation)
        integer(kind=4) :: num_waters_init=0                            !< Initial number of target waters (reference for comparison/restart)
        integer(kind=4) :: num_target_waters_init=0                        !< Initial number of target waters in computational domain (excludes boundaries)
        type(aqueous_chemistry_c), allocatable :: wat_types(:)                 !< Array of water type compositions (initial + boundary conditions)
        type(aqueous_chemistry_c), allocatable :: waters(:)             !< Current target water compositions at spatial locations
        type(aqueous_chemistry_c), allocatable :: waters_init(:)        !< Initial target water compositions (backup for restart/mass balance)
        type(zone_water_c), allocatable :: zone_waters(:,:)              !< Array of reactive zone waters
        integer(kind=4), allocatable :: tar_wat_indices(:)                 !< Indices of domain target waters (excludes boundaries)
        integer(kind=4), allocatable :: tar_wat_indices_init(:)            !< Initial domain target water indices (for restart)
        integer(kind=4), allocatable :: upstream_water_indices(:)             !< Index of closest upstream water for each target water (for transport coupling)
        integer(kind=4), allocatable :: downstream_water_indices(:)           !< Index of closest downstream water for each target water (for transport coupling)
        !> @}
        
        !> @name External Water Sources
        !> @{
        integer(kind=4) :: num_ext_waters=0                                    !< Number of external water sources entering domain from outside
        integer(kind=4) :: num_rech_waters=0                                   !< Number of recharge water sources (infiltration, injection wells)
        integer(kind=4) :: num_bd_waters=0                                     !< Number of boundary waters (inflow/outflow conditions)
        integer(kind=4) :: num_left_bd_waters=0                                     !< Number of boundary waters (inflow/outflow conditions)
        integer(kind=4) :: num_right_bd_waters=0                                     !< Number of boundary waters (inflow/outflow conditions)
        integer(kind=4), allocatable :: ext_waters_indices(:)                  !< Indices of external water sources in waters array
        integer(kind=4), allocatable :: rech_waters_indices(:)                 !< Indices of recharge water sources in waters array
        integer(kind=4), allocatable :: bd_waters_indices(:)                   !< Indices of boundary waters in waters array
        integer(kind=4), allocatable :: wat_type_indices(:)                    !< Water type index (into wat_types) per water; set in read_tar_wat_line
        !> @}
        
        !> @name Solid Phase Management
        !> @{
        integer(kind=4) :: num_target_solids=0                                 !< Number of target solids (≤ num_waters)
        integer(kind=4) :: num_target_solids_dummy=0                           !< Number of dummy target solids (for locations without solids)
        integer(kind=4) :: num_materials=0                                     !< Number of distinct material types (≤ num_target_solids)
        integer(kind=4) :: num_init_cat_exch_zones=0                           !< Number of initial cation exchange zones for sorption
        integer(kind=4) :: num_interfaces=0                                     !< Number of solid-solid interfaces
        type(solid_chemistry_c), allocatable :: materials(:)                   !< Array of distinct material types (solid compositions)
        type(solid_chemistry_c), allocatable :: init_cat_exch_zones(:)         !< Initial cation exchange zone definitions
        type(solid_chemistry_c), allocatable :: sol_types(:)                   !< Solid types associated with initial & boundary water types
        type(solid_chemistry_c), allocatable :: target_solids(:)               !< Current target solid compositions at spatial locations
        type(solid_chemistry_c), allocatable :: target_solids_init(:)          !< Initial target solid compositions (for restart/mass balance)
        type(solid_chemistry_c), allocatable :: target_solids_dummy(:)         !< Dummy target solids for locations without associated solids
        type(solid_chemistry_c), allocatable :: zone_solids(:,:)              !< Array of solid chemistry objects associated with zone_waters
        !> @}
        
        !> @name Gas Phase Management
        !> @{
        integer(kind=4) :: num_target_gases=0                                  !< Number of target gas phases
        integer(kind=4) :: num_gas_zones=0                                     !< Number of distinct gas zone types
        type(gas_chemistry_c), allocatable :: gas_zones(:)                     !< Array of distinct gas zone types
        type(gas_chemistry_c), allocatable :: gas_zones_wat_types(:)           !< Gas zones associated with water types
        type(gas_chemistry_c), allocatable :: target_gases(:)                  !< Current target gas compositions at spatial locations
        type(gas_chemistry_c), allocatable :: target_gases_init(:)             !< Initial target gas compositions (for restart)
        !> @}
        
        !> @name Reactive Zone Management
        !> @{
        integer(kind=4) :: num_reactive_zones=0                                !< Number of reactive zones (≤ num_target_solids)
        integer(kind=4) :: num_reactive_zones_dummy=0                          !< Number of dummy reactive zones (for locations without heterogeneous equilibrium reactions)
        type(reactive_zone_c), allocatable :: reactive_zones(:)                !< Array of reactive zone definitions (reaction networks)
        type(reactive_zone_c), allocatable :: react_zones_wat_types(:)         !< Reactive zones associated with water types
        type(reactive_zone_c), allocatable :: react_zones_dummy(:)             !< Dummy reactive zones for waters with no heterogeneous equilibrium reactions
        type(int_array_c) :: wat_indices_zones                             !< Per reactive zone (column j = zone id): indices of the target waters whose previous-step state (solid_chemistry_old) belongs to reactive_zones(j); the cell grouping in the Eulerian case, refreshed after each Lagrangian displacement (reset_wat_indices_zoness, which also classifies the Neumann outflow boundary waters holding the parcels that left the last cells); the within-zone groups of cross_zone_mix_core
        !> @name Mineral Zone Management
        !> @{
        integer(kind=4) :: num_mineral_zones=0                                 !< Number of mineral zones (≤ num_target_solids)
        type(mineral_zone_c), allocatable :: mineral_zones(:)                  !< Array of mineral zone definitions
        type(mineral_zone_c), allocatable :: min_zones_wat_types(:)            !< Mineral zones associated with water types (analogous to react_zones_wat_types)
        type(mineral_zone_c) :: min_zone_dummy                                 !< Dummy mineral zone for locations without minerals
        !> @}
    contains
        !> @name Configuration and Setup Methods
        !> @{
        procedure :: set_lump_flag                    !< Set spatial lumping flag for Water Mixing Algorithm
        procedure :: set_mix_by_zones_flag            !< Set mixing flag (per reactive zone vs. directly over all waters)
        procedure :: set_rk_avg_opt                   !< Set reaction rate averaging method
        procedure :: set_r_down_opt                  !< Set downstream reaction rate estimation method
        procedure :: set_est_prm                     !< Set estimation parameter for downstream kinetic reactions
        procedure :: set_chem_syst                    !< Set complete chemical system
        procedure :: set_num_wat_types                !< Set number of water types
        procedure :: set_num_materials                !< Set number of material types
        procedure :: set_num_waters                  !< Set number of target waters
        procedure :: set_num_target_waters               !< Set number of domain target waters
        procedure :: set_num_target_solids            !< Set number of target solids
        procedure :: set_num_ext_waters               !< Set number of external waters
        procedure :: set_ext_waters_indices           !< Set external water indices
        procedure :: set_num_rech_waters              !< Set number of recharge waters
        procedure :: set_num_bd_waters                !< Set number of boundary waters
        procedure :: set_waters_target_solids  !< Link target waters to target solids
        procedure :: set_target_solids_materials      !< Link target solids to materials
        procedure :: set_target_gases                 !< Set target gas phases
        procedure :: set_read_opt                     !< Set chemical data input format option
        procedure :: set_reactive_zones               !< Set reactive zone definitions
        procedure :: set_Jac_opt                      !< Set Jacobian computation method
        procedure :: set_target_solids_mesh           !< Link target solids to spatial mesh
        procedure :: set_upstream_water_indices       !< Set upstream water indices for transport coupling
        procedure :: set_downstream_water_indices     !< Set downstream water indices for transport coupling
        procedure :: set_num_interfaces               !< Set number of interfaces
        !> @}
        
        !> @name Query and Retrieval Methods
        !> @{
        procedure :: get_num_aq_comps_chem_syst       !< Get number of aqueous components in chemical system
        procedure :: get_num_aq_comps_tar_wat         !< Get number of aqueous components in a domain (target) water
        procedure :: get_num_wat_types                !< Get number of water types
        procedure :: get_conc_comp_wat_types          !< Get component concentrations of water types
        procedure :: get_num_aq_var_act_species       !< Get number of aqueous variable activity species
        procedure :: get_tar_sol_ind                  !< Get target solid index
        procedure :: get_tar_gas_ind                  !< Get target gas index
        !> @}
        
        !> @name Memory Management Methods
        !> @{
        procedure :: allocate_waters           !< Allocate target waters array
        procedure :: allocate_waters_init      !< Allocate initial target waters array
        procedure :: allocate_tar_wat_indices     !< Allocate domain target water indices
        procedure :: allocate_tar_wat_indices_init !< Allocate initial domain target water indices
        procedure :: allocate_ext_waters_indices      !< Allocate external water indices
        procedure :: allocate_rech_waters_indices     !< Allocate recharge water indices
        procedure :: allocate_bd_waters_indices       !< Allocate boundary water indices
        procedure :: allocate_target_solids           !< Allocate target solids array
        procedure :: allocate_target_solids_dummy     !< Allocate dummy target solids
        procedure :: allocate_target_gases            !< Allocate target gases array
        procedure :: allocate_reactive_zones          !< Allocate reactive zones array
        procedure :: allocate_react_zones_wat_types   !< Allocate reactive zones for water types
        procedure :: allocate_min_zones_wat_types     !< Allocate mineral zones for water types
        procedure :: allocate_mineral_zones           !< Allocate mineral zones array
        procedure :: allocate_materials               !< Allocate materials array
        procedure :: allocate_init_cat_exch_zones     !< Allocate initial cation exchange zones
        procedure :: allocate_wat_types               !< Allocate water types array
        procedure :: allocate_sol_types               !< Allocate solid types array
        procedure :: allocate_gas_zones               !< Allocate gas zones array
        procedure :: allocate_gas_zones_wat_types     !< Allocate gas zones for water types
        procedure :: allocate_reactive_zones_dummy          !< Allocate dummy reactive zones
        procedure :: allocate_zone_waters           !< Allocate zone waters array
        procedure :: allocate_zone_solids           !< Allocate zone solids array
        
        !> @}
        
        !> @name Input/Output Methods
        !> @{
        procedure :: read_waters_init          !< Read initial target water definitions
        procedure :: read_tar_wat_line                !< Read single target water line from input
        procedure :: read_tar_sol                     !< Read target solid definitions
        procedure :: read_tar_gas                     !< Read target gas definitions
        procedure :: read_init_min_zones_CHEPROO      !< Read initial mineral zones (CHEPROO format)
        procedure :: read_chemistry                   !< Read complete chemistry setup
        procedure :: read_chemistry_interface         !< Read chemistry setup for the reactive-mixing interface executable (optional CV_params, no comp/out options)
        procedure :: read_chemistry_CHEPROO           !< Read chemistry data in CHEPROO format
        procedure :: read_init_bd_wat_types_CHEPROO   !< Read initial boundary water types (CHEPROO)
        procedure :: read_init_cat_exch_zones_CHEPROO !< Read initial cation exchange zones (CHEPROO)
        procedure :: read_gas_zones_CHEPROO             !< Read initial gas zones (CHEPROO format)
        procedure :: read_comp_opts                   !< Read computation options from input
        procedure :: write_chemistry                  !< Write chemistry state to output
        procedure :: write_conc_comp_wat_types        !< Write component concentrations of water types
        procedure :: write_u_mix_init                   !< Write initial u_mix concentrations
        procedure :: write_u_init                       !< Write initial aqueous components and species to file
        !> @}
        
        !> @name Initialization Methods
        !> @{
        procedure :: initialise_chemistry             !< Initialize complete chemistry system
        !> @}
        
        !> @name Main Solver Methods
        !> @{
        procedure :: solve_reactive_mixing_lump              !< Main reactive mixing solver with spatial lumping
        procedure :: solve_reactive_mixing_ideal_cons        !< Ideal reactive mixing solver with consistency
        procedure :: solve_reactive_mixing_cons              !< General reactive mixing solver with consistency
        procedure :: solve_reactive_mixing_ideal_lump        !< Ideal reactive mixing solver with lumping
        procedure :: solve_reactive_mixing_ideal_lump_iter   !< Single time-step iteration of ideal lumped solver
        procedure :: solve_reactive_mixing         !< Reactive mixing with time-dependent boundary conditions
        procedure :: write_react_mix_mass_bal_dom   !< Write the domain-integrated reactive-mixing mass balance (per species) for one time step
        !> @}
        
        !> @name Linkage and Association Methods
        !> @{
        procedure :: link_waters_target_solids !< Link target waters to corresponding target solids
        procedure :: link_target_waters_target_gases  !< Link target waters to corresponding target gases
        procedure :: link_target_solids_reactive_zone !< Link target solids to reactive zones
        procedure :: link_target_gases_reactive_zone  !< Link target gases to reactive zones
        procedure :: set_wat_indices_zones !< Link waters to a reactive zone
        procedure :: link_target_waters_old_reaction_zones !< Group target waters by their previous-time-step reaction zones
        procedure :: link_waters_mineral_zone !< Link target waters to mineral zones
        procedure :: link_target_solids_mineral_zone !< Link target solids to mineral zones
        procedure :: associate_zone_waters_zone_solids
        procedure :: set_zone_waters !< Allocate and link the per-(target water, reactive zone) zone waters and zone solids
        !> @}
        
        !> @name Validation and Checking Methods
        !> @{
        procedure :: check_new_reactive_zones         !< Check for new reactive zones during simulation
        !> @}
        
        !> @name Utility and Helper Methods
        !> @{
        procedure :: compute_u_mix_init            !< Compute initial u_mix concentrations
        procedure :: WMA_iter_EE_eq_kin_ideal_lump  !< Water Mixing Approach iteration using EE method with equilibrium and kinetics (ideal, lumped)
        procedure :: compute_conc_bd_waters_Lagr   !< Compute boundary water concentrations using Lagrangian extrapolation
        procedure :: reorder_mix_react_indices      !< Reorder mixing reaction indices for upstream/downstream processing
        !> @}
    end type !< End of chemistry_c type definition
    
    interface
    !> @brief Solve reactive mixing using lumped approach
    !> @details Performs reactive mixing calculations using a lumped parameter approach for all domain and external waters.
    !> Uses provided mixing ratios and water indices, with time discretization and integration method for chemical reactions.
    !> @param this Chemistry object
    !> @param root Output file root
    !> @param mixing_ratios Mixing ratios matrix
    !> @param mix_conc_indices Indices of mixing waters
    !> @param time_discr Time discretization object
    !> @param int_method_chem_reacts Integration method for chemical reactions
    subroutine solve_reactive_mixing_lump(this,dir,root,mixing_ratios_conc,mixing_ratios_R,mix_conc_indices,&
        time_discr,theta_r)
            import chemistry_c                            !< Import chemistry class definition
            import real_array_c                           !< Import real array class for mixing ratios
            import int_array_c                            !< Import integer array class for indices
            import time_discr_c                           !< Import time discretization class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object instance
            character(len=*), intent(in) :: dir           !< Directory for output files
            character(len=*), intent(in) :: root          !< Root name for output files
            class(real_array_c), intent(in) :: mixing_ratios_conc !< Mixing ratios for concentrations
            class(real_array_c), intent(in) :: mixing_ratios_R !< Mixing ratios for reaction rates (mix_conc_indices layout)
            class(int_array_c), intent(in) :: mix_conc_indices !< Indices of mixing waters
            class(time_discr_c), intent(in) :: time_discr !< Time discretization object
            real(kind=8), intent(in) :: theta_r           !< Integration method for chemical reactions
        end subroutine
        
        !> @brief Solve ideal conservative reactive mixing
        !> @details Solves reactive mixing for ideal conservative systems, using concentration and reaction rate mixing ratios.
        !> Handles domain and external waters, with time discretization and integration method for chemical reactions.
        !> @param this Chemistry object
        !> @param root Output file root
        !> @param mixing_ratios_conc Mixing ratios for concentrations
        !> @param mixing_ratios_R Mixing ratios for reaction rates (mix_conc_indices layout, 0 for conc-only waters)
        !> @param mix_conc_indices Indices of mixing waters
        !> @param time_discr Time discretization object
        !> @param int_method_chem_reacts Integration method for chemical reactions
        subroutine solve_reactive_mixing_ideal_cons(this,dir,root,mixing_ratios_conc,&
            mixing_ratios_R,mix_conc_indices,&
            time_discr,theta_r)
            import chemistry_c                            !< Import chemistry class definition
            import real_array_c                           !< Import real array class for mixing ratios
            import int_array_c                            !< Import integer array class for indices
            import time_discr_c                           !< Import time discretization class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object instance
            character(len=*), intent(in) :: dir           !< Directory for output files
            character(len=*), intent(in) :: root          !< Root name for output files
            class(real_array_c), intent(in) :: mixing_ratios_conc !< Mixing ratios for concentrations
            class(real_array_c), intent(in) :: mixing_ratios_R !< Mixing ratios for reaction rates (mix_conc_indices layout)
            class(int_array_c), intent(in) :: mix_conc_indices !< Indices of mixing waters
            class(time_discr_c), intent(in) :: time_discr !< Time discretization object
            real(kind=8), intent(in) :: theta_r           !< Integration method for chemical reactions
        end subroutine

        !> @brief Solve conservative reactive mixing
        !> @details Solves reactive mixing for conservative systems, using concentration and reaction rate mixing ratios.
        !> Handles domain and external waters, with time discretization and integration method for chemical reactions.
        !> @param this Chemistry object
        !> @param root Output file root
        !> @param mixing_ratios_conc Mixing ratios for concentrations
        !> @param mixing_ratios_R Mixing ratios for reaction rates (mix_conc_indices layout, 0 for conc-only waters)
        !> @param mix_conc_indices Indices of mixing waters
        !> @param time_discr Time discretization object
        !> @param int_method_chem_reacts Integration method for chemical reactions
        subroutine solve_reactive_mixing_cons(this,dir,root,mixing_ratios_conc,mixing_ratios_R,mix_conc_indices,&
            time_discr,theta_r)
            import chemistry_c                            !< Import chemistry class definition
            import real_array_c                           !< Import real array class for mixing ratios
            import int_array_c                            !< Import integer array class for indices
            import time_discr_c                           !< Import time discretization class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object instance
            character(len=*), intent(in) :: dir           !< Directory for output files
            character(len=*), intent(in) :: root          !< Root name for output files
            class(real_array_c), intent(in) :: mixing_ratios_conc !< Mixing ratios for concentrations
            class(real_array_c), intent(in) :: mixing_ratios_R !< Mixing ratios for reaction rates (mix_conc_indices layout)
            class(int_array_c), intent(in) :: mix_conc_indices !< Indices of mixing waters
            class(time_discr_c), intent(in) :: time_discr !< Time discretization object
            real(kind=8), intent(in) :: theta_r           !< Integration method for chemical reactions
        end subroutine
            
        !> @brief Solve ideal lumped reactive mixing
        !> @details Performs reactive mixing calculations using an ideal lumped parameter approach for all domain and external waters.
        !> Uses provided mixing ratios and water indices, with time discretization and integration method for chemical reactions.
        !> @param this Chemistry object
        !> @param root Output file root
        !> @param mixing_ratios Mixing ratios matrix
        !> @param mix_conc_indices Indices of mixing waters
        !> @param time_discr Time discretization object
        !> @param int_method_chem_reacts Integration method for chemical reactions
        subroutine solve_reactive_mixing_ideal_lump(this,dir,root,mixing_ratios_conc,mixing_ratios_R,mix_conc_indices,&
                time_discr,theta_r)
            import chemistry_c                            !< Import chemistry class definition
            import real_array_c                           !< Import real array class for mixing ratios
            import int_array_c                            !< Import integer array class for indices
            import time_discr_c                           !< Import time discretization class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object instance
            character(len=*), intent(in) :: dir           !< Directory for output files
            character(len=*), intent(in) :: root          !< Root name for output files
            class(real_array_c), intent(in) :: mixing_ratios_conc !< Mixing ratios for concentrations
            class(real_array_c), intent(in) :: mixing_ratios_R !< Mixing ratios for reaction rates (mix_conc_indices layout)
            class(int_array_c), intent(in) :: mix_conc_indices !< Indices of mixing waters
            class(time_discr_c), intent(in) :: time_discr !< Time discretization object
            real(kind=8), intent(in) :: theta_r           !< Integration method for chemical reactions
        end subroutine
        
        
        !> @brief Solve reactive mixing with time-dependent boundary conditions
        !> @details Solves reactive mixing for systems with time-dependent boundary conditions, using provided mixing ratios, water indices, and transport parameters.
        !> Includes time and spatial discretization, integration method, and analytical solution for validation.
        !> @param this Chemistry object
        !> @param root Output file root
        !> @param unit File unit for output
        !> @param mixing_ratios Mixing ratios matrix
        !> @param mix_conc_indices Indices of mixing waters
        !> @param time_discr_tpt Time discretization object for transport
        !> @param int_method_chem_reacts Integration method for chemical reactions
        !> @param spatial_discr_tpt Spatial discretization object for transport
        !> @param D Diffusion coefficient
        !> @param q Flow rate
        !> @param phi Porosity
        !> @param anal_sol Analytical solution function
        subroutine solve_reactive_mixing_BCs_dep_t(this,root,unit,mixing_ratios,mix_conc_indices,time_discr_tpt,&
            theta_r,spatial_discr_tpt,D,q,phi,anal_sol)
            import chemistry_c                            !< Import chemistry class definition
            import spatial_discr_c                        !< Import spatial discretization class
            import time_discr_c                           !< Import time discretization class
            import real_array_c                           !< Import real array class for mixing ratios
            import int_array_c                            !< Import integer array class for indices
            implicit none                                 !< Require explicit variable declarations
        !> Arguments
            class(chemistry_c) :: this                    !< Chemistry object instance
            character(len=*), intent(in) :: root          !< Root name for output files
            integer(kind=4), intent(in) :: unit           !< File unit for output
            class(real_array_c), intent(in) :: mixing_ratios !< Mixing ratios matrix
            class(int_array_c), intent(in) :: mix_conc_indices !< Indices of mixing waters
            class(time_discr_c), intent(in) :: time_discr_tpt !< Time discretization for transport
            real(kind=8), intent(in) :: theta_r           !< Integration method for chemical reactions
            class(spatial_discr_c), intent(in) :: spatial_discr_tpt !< Spatial discretization for transport
            real(kind=8), intent(in) :: D                 !< Diffusion coefficient
            real(kind=8), intent(in) :: q                 !< Flow rate
            real(kind=8), intent(in) :: phi               !< Porosity
            real(kind=8), external :: anal_sol             !< Analytical solution function
        end subroutine
        
    !> @brief Read chemistry data from PHREEQC files
    !> @details Reads chemical system and database information from PHREEQC input and database files.
    !> Used for initializing chemistry objects from external geochemical databases.
    !> @param this Chemistry object
    !> @param path_inp Path to PHREEQC input file
    !> @param path_DB Path to PHREEQC database file
    !> @param filename Input filename
    subroutine read_chemistry_PHREEQC(this,path_inp,path_DB,filename)
            import chemistry_c                            !< Import chemistry class definition
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object to initialize
            character(len=*), intent(in) :: path_inp      !< Path to PHREEQC input file
            character(len=*), intent(in) :: path_DB       !< Path to PHREEQC database file
            character(len=*), intent(in) :: filename      !< Input filename
        end subroutine
        
    !> @brief Read chemistry data from CHEPROO files
    !> @details Reads chemical system and database information from CHEPROO input and database files.
    !> Used for initializing chemistry objects from external CHEPROO databases.
    !> @param this Chemistry object
    !> @param root Root name for files
    !> @param path_pb Path to problem files
    !> @param path_DB Path to CHEPROO database file
    !> @param unit_chem_syst_file File unit for chemical system file
    !> @param unit_loc_chem_file File unit for local chemistry file
    subroutine read_chemistry_CHEPROO(this,root,path_pb,path_DB,unit_chem_syst_file,unit_loc_chem_file)
            import chemistry_c                            !< Import chemistry class definition
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object to initialize
            character(len=*), intent(in) :: root          !< Root name for files
            character(len=*), intent(in) :: path_pb       !< Path to problem files
            character(len=*), intent(in) :: path_DB       !< Path to CHEPROO database file
            integer(kind=4), intent(in) :: unit_chem_syst_file !< File unit for chemical system
            integer(kind=4), intent(in) :: unit_loc_chem_file !< File unit for local chemistry
        end subroutine
        
    !> @brief Read complete chemistry setup from files
    !> @details Reads and initializes the complete chemistry configuration including chemical system,
    !>          water types, solid/gas phases, and reactive zones from input files.
    !> @param[in,out] this Chemistry object to initialize
    !> @param[in] root Root name for input/output files
    !> @param[in] path_pb Path to problem directory
    !> @param[in] path_DB Path to database files
        subroutine read_chemistry(this,root,path_pb,path_DB)
            import chemistry_c
            class(chemistry_c) :: this                    !< Chemistry object to initialize
            character(len=*), intent(in) :: root          !< Root name for input/output files
            character(len=*), intent(in) :: path_pb       !< Path to problem directory
            character(len=*), intent(in) :: path_DB       !< Path to database files
        end subroutine

    !> @brief Read chemistry setup for the reactive-mixing interface executable
    !> @details Same CHEPROO reading and reactive-zone linking as read_chemistry, but the
    !>          convergence parameters file is optional (type defaults otherwise) and the
    !>          computation / output options files are not read (the interface writes its own output).
    !> @param[in,out] this Chemistry object to initialize
    !> @param[in] root Root name for input/output files
    !> @param[in] path_pb Path to problem directory
    !> @param[in] path_DB Path to database files
    !> @param[in] num_tar Number of target (domain) waters expected from the mesh
        subroutine read_chemistry_interface(this,root,path_pb,path_DB,num_tar)
            import chemistry_c
            class(chemistry_c) :: this                    !< Chemistry object to initialize
            character(len=*), intent(in) :: root          !< Root name for input/output files
            character(len=*), intent(in) :: path_pb       !< Path to problem directory
            character(len=*), intent(in) :: path_DB       !< Path to database files
            integer(kind=4), intent(in) :: num_tar        !< Expected number of target (domain) waters from the mesh
        end subroutine
                
    !> @brief Read reactive zones in Lagrangian formulation from file
    !> @details Reads reactive zone definitions for Lagrangian reactive transport,
    !>          including zone assignments to spatial locations and chemical properties.
    !> @param[in,out] this Chemistry object to populate with reactive zones
    !> @param[in] unit File unit number for reading reactive zone data
        subroutine read_reactive_zones_Lagr(this,unit)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this                    !< Chemistry object to populate
            integer(kind=4), intent(in) :: unit           !< File unit number for input
        end subroutine
        
    !> @brief Link domain target waters to a specific reactive zone
    !> @details Finds all domain target water indices whose reactive zone matches reactive zone i.
    !>          Returns an allocatable array of water indices belonging to the zone.
    !> @param[in,out] this Chemistry object containing target waters and reactive zones
    !> @param[in] i Index of the reactive zone to query
    !> @param[out] tar_wat_indices Allocatable array of water indices in zone i
    !> @param[in] exit_wat_indices Optional: boundary waters holding parcels that left the domain (Neumann outflow), classified like the target waters
        subroutine set_wat_indices_zones(this)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this                    !< Chemistry object instance
            !integer(kind=4), intent(in) :: i              !< Reactive zone index to query
            !integer(kind=4), intent(out), allocatable :: tar_wat_indices_rz(:) !< Output: target water indices in this zone
            !integer(kind=4), intent(in), optional :: exit_wat_indices(:) !< Global indices of the exited-parcel boundary waters to classify too
        end subroutine
        
    !> @brief Group target waters by their previous-time-step reaction zone
    !> @details Groups all target waters into distinct reaction zones according to their
    !>          full previous-time-step reaction set (equilibrium + kinetic reactions).
    !> @param[in] this Chemistry object containing target waters and reactive zones
    !> @param[out] num_rz Number of distinct old reaction zones found
    !> @param[out] tar_wat_indices_rz Per old reaction zone (column): global target-water indices in that zone
        subroutine link_target_waters_old_reaction_zones(this,num_rz,tar_wat_indices_rz)
            import chemistry_c
            import int_array_c
            implicit none
            class(chemistry_c), intent(in) :: this        !< Chemistry object instance
            integer(kind=4), intent(out) :: num_rz        !< Output: number of old reaction zones
            class(int_array_c), intent(out) :: tar_wat_indices_rz !< Output: per old reaction zone, global target-water indices
        end subroutine
        
    !> @brief Link target waters to a specific mineral zone
    !> @details Finds all target water indices associated with mineral zone i.
    !>          Returns an allocatable array of water indices belonging to the mineral zone.
    !> @param[in,out] this Chemistry object containing waters and mineral zones
    !> @param[in] i Index of the mineral zone to query
    !> @param[out] wat_indices Allocatable array of water indices in mineral zone i
        subroutine link_waters_mineral_zone(this,i,wat_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this                    !< Chemistry object instance
            integer(kind=4), intent(in) :: i              !< Mineral zone index to query
            integer(kind=4), intent(out), allocatable :: wat_indices(:) !< Output: water indices in this mineral zone
        end subroutine
        
    !> @brief Link target solids to a specific mineral zone
    !> @details Finds all target solid indices associated with mineral zone i.
    !>          Returns an allocatable array of target solid indices belonging to the zone.
    !> @param[in] this Chemistry object containing target solids and mineral zones
    !> @param[in] i Index of the mineral zone to query
    !> @param[out] tar_sol_indices Allocatable array of target solid indices in mineral zone i
        subroutine link_target_solids_mineral_zone(this,i,tar_sol_indices)
            import chemistry_c
            implicit none
            class(chemistry_c), intent(in) :: this        !< Chemistry object instance
            integer(kind=4), intent(in) :: i              !< Mineral zone index to query
            integer(kind=4), intent(out), allocatable :: tar_sol_indices(:) !< Output: target solid indices in this mineral zone
        end subroutine
        
    !> @brief Link target waters to a specific target solid
    !> @details Finds all target water indices associated with target solid i.
    !>          Returns an allocatable array of target water indices linked to the solid.
    !> @param[in,out] this Chemistry object containing waters and solids
    !> @param[in] i Index of the target solid to query
    !> @param[out] tw_indices Allocatable array of target water indices linked to solid i
        subroutine link_target_waters_target_solid(this,i,tw_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this                    !< Chemistry object instance
            integer(kind=4), intent(in) :: i              !< Target solid index to query
            integer(kind=4), intent(out), allocatable :: tw_indices(:) !< Output: target water indices linked to this solid
        end subroutine
        
    !> @brief Link multiple target waters to multiple target solids
    !> @details Finds all water indices associated with a set of target solid indices.
    !>          Returns an allocatable array of water indices linked to the specified solids.
    !> @param[in,out] this Chemistry object containing waters and solids
    !> @param[in] tar_sol_indices Array of target solid indices to query
    !> @param[out] wat_indices Allocatable array of water indices linked to specified solids
        subroutine link_waters_target_solids(this,tar_sol_indices,wat_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this                    !< Chemistry object instance
            integer(kind=4), intent(in) :: tar_sol_indices(:) !< Target solid indices to query
            integer(kind=4), intent(out), allocatable :: wat_indices(:) !< Output: water indices linked to these solids
        end subroutine
        
    !> @brief Link target waters to target gas phases
    !> @details Finds all target water indices associated with a set of target gas indices.
    !>          Returns an allocatable array of target water indices linked to the specified gases.
    !> @param[in,out] this Chemistry object containing waters and gases
    !> @param[in] tar_gas_indices Array of target gas indices to query
    !> @param[in,out] tar_wat_indices Allocatable array of target water indices linked to gases
        subroutine link_target_waters_target_gases(this,tar_gas_indices,tar_wat_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this                    !< Chemistry object instance
            integer(kind=4), intent(in) :: tar_gas_indices(:) !< Target gas indices to query
            integer(kind=4), intent(inout), allocatable :: tar_wat_indices(:) !< Output: target water indices linked to these gases
        end subroutine
        
    !> @brief Link target solids to a specific reactive zone
    !> @details Finds all target solid indices associated with reactive zone i.
    !>          Returns an allocatable array of target solid indices belonging to the zone.
    !> @param[in,out] this Chemistry object containing solids and reactive zones
    !> @param[in] i Index of the reactive zone to query
    !> @param[out] tar_sol_indices Allocatable array of target solid indices in zone i
        subroutine link_target_solids_reactive_zone(this,i,tar_sol_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this                    !< Chemistry object instance
            integer(kind=4), intent(in) :: i              !< Reactive zone index to query
            integer(kind=4), intent(out), allocatable :: tar_sol_indices(:) !< Output: target solid indices in this zone
        end subroutine
        
    !> @brief Link target gases to a specific reactive zone
    !> @details Finds all target gas indices associated with reactive zone i.
    !>          Returns an allocatable array of target gas indices belonging to the zone.
    !> @param[in,out] this Chemistry object containing gases and reactive zones
    !> @param[in] i Index of the reactive zone to query
    !> @param[out] tar_gas_indices Allocatable array of target gas indices in zone i
        subroutine link_target_gases_reactive_zone(this,i,tar_gas_indices) 
            import chemistry_c
            implicit none
            class(chemistry_c) :: this                    !< Chemistry object instance
            integer(kind=4), intent(in) :: i              !< Reactive zone index to query
            integer(kind=4), intent(out), allocatable :: tar_gas_indices(:) !< Output: target gas indices in this zone
        end subroutine
        
    !> @brief Initialize domain waters from initial water types
    !> @details Sets up domain waters using an array of initial water type compositions.
    !>          Copies aqueous chemistry from water types to corresponding domain water locations.
    !> @param[in,out] this Chemistry object to initialize
    !> @param[in] initial_water_types Array of initial water type aqueous chemistry objects
        subroutine initialise_waters_dom(this,initial_water_types)
            import chemistry_c                            !< Import chemistry class definition
            import aqueous_chemistry_c                    !< Import aqueous chemistry class
            import chem_system_c                          !< Import chemical system class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object to initialize
            class(aqueous_chemistry_c), intent(in) :: initial_water_types(:) !< Array of initial water type compositions
        end subroutine
        
    !> @brief Read initial target waters from file
    !> @details Reads initial target water definitions from input files, including their
    !>          water type assignments, solid zone associations, and gas zone associations.
    !> @param[in,out] this Chemistry object to populate with initial waters
    !> @param[in] root Root name for input files
    !> @param[in] nsrz Number of solid reactive zones
    !> @param[in] ngrz Number of gas reactive zones
        subroutine read_waters_init(this,root,nsrz,ngrz)
            import chemistry_c                            !< Import chemistry class definition
            import aqueous_chemistry_c                    !< Import aqueous chemistry class
            import solid_chemistry_c                      !< Import solid chemistry class
            import gas_chemistry_c                        !< Import gas chemistry class
            import reactive_zone_c                        !< Import reactive zone class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object to populate
            character(len=*), intent(in) :: root          !< Root name for input files
            integer(kind=4), intent(in) :: nsrz           !< Number of solid reactive zones
            integer(kind=4), intent(in) :: ngrz           !< Number of gas reactive zones
        end subroutine
        
    !> @brief Read domain waters with iterative equilibrium speciation
    !> @details Reads domain water compositions and performs iterative speciation calculations
    !>          to achieve chemical equilibrium. Returns convergence status and iteration count.
    !> @param[in,out] this Chemistry object
    !> @param[in] unit File unit number for reading
    !> @param[in] init_water_types Array of initial water type compositions
    !> @param[in] bd_water_types Array of boundary water type compositions
    !> @param[in] init_sol_types Array of initial solid type chemistry objects
    !> @param[out] niter Number of iterations performed for convergence
    !> @param[out] CV_flag Convergence flag (.true. if converged, .false. otherwise)
        subroutine read_waters_dom_bis(this,unit,init_water_types,bd_water_types,init_sol_types,niter,CV_flag)
            import chemistry_c                            !< Import chemistry class definition
            import aqueous_chemistry_c                    !< Import aqueous chemistry class
            import solid_chemistry_c                      !< Import solid chemistry class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object
            integer(kind=4), intent(in) :: unit           !< File unit number for reading
            class(aqueous_chemistry_c), intent(in) :: init_water_types(:) !< Initial water type compositions
            class(aqueous_chemistry_c), intent(in) :: bd_water_types(:) !< Boundary water type compositions
            class(solid_chemistry_c), intent(in) :: init_sol_types(:) !< Initial solid type chemistry objects
            integer(kind=4), intent(out) :: niter         !< Number of iterations for convergence
            logical, intent(out) :: CV_flag               !< Convergence flag (.true.=converged)
        end subroutine
        
    !> @brief Initialize target solids array with a specific count
    !> @details Allocates and initializes target solid chemistry objects for a given number
    !>          of spatial locations in the reactive transport domain.
    !> @param[in,out] this Chemistry object to initialize
    !> @param[in] n Number of target solids to initialize
        subroutine initialise_target_solids(this,n)
            import chemistry_c                            !< Import chemistry class definition
            import solid_chemistry_c                      !< Import solid chemistry class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object to initialize
            integer(kind=4), intent(in) :: n              !< Number of target solids to create
        end subroutine
        
    !> @brief Initialize external waters from problem files
    !> @details Reads and initializes external (boundary/recharge) water definitions
    !>          from problem files using root path and database path.
    !> @param[in,out] this Chemistry object to populate with external waters
        subroutine initialise_ext_waters(this)
            import chemistry_c                            !< Import chemistry class definition
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object to populate
        end subroutine
        
    !> @brief Write chemistry state to output files
    !> @details Writes the current state of the chemistry system to output files,
    !>          including water compositions, solid phases, and gas phases.
    !>          Optionally writes for a specific time step.
    !> @param[in] this Chemistry object to output
    !> @param[in] dir Directory path for output files
    !> @param[in] root Root name for output files
    !> @param[in] time_step Optional time step number for time-series output
        subroutine write_chemistry(this,dir,root,time_step,Lagr_flag,num_time)
            import chemistry_c
            class(chemistry_c), intent(in) :: this        !< Chemistry object to output
            character(len=*), intent(in) :: dir           !< Directory for output files
            character(len=*), intent(in) :: root          !< Root name for output files
            integer(kind=4), intent(in), optional :: time_step !< Optional time step number
            logical, intent(in), optional :: Lagr_flag    !< Transport approach: .true.=Lagrangian, .false.=Eulerian (used in output filename)
            integer(kind=4), intent(in), optional :: num_time !< Total number of time steps in the run (appended to output filename)
        end subroutine
        
    !> @brief Read initial and boundary water types from CHEPROO input
    !> @details Reads initial and boundary water type definitions from CHEPROO format files,
    !>          optionally including gas phase chemistry for multiphase systems.
    !> @param[in,out] this Chemistry object to populate with water types
    !> @param[in] unit File unit number for reading
    !> @param[in] gas_species_chem Optional gas chemistry for multiphase initialization
        subroutine read_init_bd_wat_types_CHEPROO(this,unit)
            import chemistry_c                            !< Import chemistry class definition
            import aqueous_chemistry_c                    !< Import aqueous chemistry class
            import solid_chemistry_c                      !< Import solid chemistry class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object to populate
            integer(kind=4), intent(in) :: unit           !< File unit number for reading
        end subroutine

    !> @brief Read initial mineral zones from CHEPROO input
    !> @details Reads mineral zone definitions from CHEPROO format input files,
    !>          returning the number of mineral reactive zones found.
    !> @param[in,out] this Chemistry object to populate
    !> @param[in] unit File unit number for reading
    !> @param[out] nmrz Number of mineral reactive zones read from file
    !> @param[in] surf_chem Optional surface chemistry for adsorption zones
        subroutine read_init_min_zones_CHEPROO(this,unit,nmrz,surf_chem)
            import chemistry_c                            !< Import chemistry class definition
            import solid_chemistry_c                      !< Import solid chemistry class
            import reactive_zone_c                        !< Import reactive zone class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object to populate
            integer(kind=4), intent(in) :: unit           !< File unit number for reading
            integer(kind=4), intent(out) :: nmrz          !< Number of mineral reactive zones
            type(solid_chemistry_c), intent(in), optional :: surf_chem !< Optional surface chemistry
        end subroutine
        
    !> @brief Read initial mineral zones from CHEPROO input (alternative format)
    !> @details Alternative reader for initial mineral zones that returns allocatable
    !>          arrays of solid chemistry objects and optionally reactive zone objects.
    !> @param[in,out] this Chemistry object
    !> @param[in] unit File unit number for reading
    !> @param[out] init_min_zones Allocatable array of initial mineral zone solid chemistry objects
    !> @param[out] reactive_zones Optional allocatable array of reactive zone objects
        subroutine read_init_min_zones_CHEPROO_bis(this,unit,init_min_zones,reactive_zones)
            import chemistry_c                            !< Import chemistry class definition
            import solid_chemistry_c                      !< Import solid chemistry class
            import reactive_zone_c                        !< Import reactive zone class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object
            integer(kind=4), intent(in) :: unit           !< File unit number for reading
            type(solid_chemistry_c), intent(out), allocatable :: init_min_zones(:) !< Output: initial mineral zone objects
            type(reactive_zone_c), intent(out), allocatable, optional :: reactive_zones(:) !< Optional output: reactive zone objects
        end subroutine

    !> @brief Read initial cation exchange zones from CHEPROO input
    !> @details Reads cation exchange zone definitions from CHEPROO format input files,
    !>          returning the number of adsorption reactive zones found.
    !> @param[in,out] this Chemistry object to populate
    !> @param[in] unit File unit number for reading
    !> @param[out] ndrz Number of adsorption reactive zones read from file
        subroutine read_init_cat_exch_zones_CHEPROO(this,unit,ndrz)
            import chemistry_c                            !< Import chemistry class definition
            import solid_chemistry_c                      !< Import solid chemistry class
            import reactive_zone_c                        !< Import reactive zone class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object to populate
            integer(kind=4), intent(in) :: unit           !< File unit number for reading
            integer(kind=4), intent(out):: ndrz           !< Number of adsorption reactive zones
        end subroutine
        
    !> @brief Read gas boundary zones from CHEPROO input
    !> @details Reads gas boundary zone definitions from CHEPROO format input files.
    !>          Returns allocatable arrays of gas chemistry objects and optionally reactive zones.
    !> @param[in,out] this Chemistry object
    !> @param[in] unit File unit number for reading
    !> @param[out] gas_bd_zones Allocatable array of gas boundary zone chemistry objects
    !> @param[out] reactive_zones Optional allocatable array of reactive zone objects
        subroutine read_gas_bd_zones_CHEPROO(this,unit,gas_bd_zones,reactive_zones)
            import chemistry_c
            import gas_chemistry_c
            import reactive_zone_c
            implicit none
            class(chemistry_c) :: this                    !< Chemistry object
            integer(kind=4), intent(in) :: unit           !< File unit number for reading
            type(gas_chemistry_c), intent(out), allocatable :: gas_bd_zones(:) !< Output: gas boundary zone objects
            type(reactive_zone_c), intent(out), allocatable, optional :: reactive_zones(:) !< Optional output: reactive zone objects
        end subroutine
        
    !> @brief Read gas zones from CHEPROO input
    !> @details Reads gas zone definitions from CHEPROO format input files,
    !>          returning the number of gas reactive zones found.
    !> @param[in,out] this Chemistry object to populate
    !> @param[in] unit File unit number for reading
    !> @param[out] ngrz Number of gas reactive zones read from file
        subroutine read_gas_zones_CHEPROO(this,unit,ngrz)
            import chemistry_c                            !< Import chemistry class definition
            import gas_chemistry_c                        !< Import gas chemistry class
            import reactive_zone_c                        !< Import reactive zone class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object to populate
            integer(kind=4), intent(in) :: unit           !< File unit number for reading
            integer(kind=4), intent(out) :: ngrz          !< Number of gas reactive zones
        end subroutine
        
    !> @brief Interface for component concentrations (equilibrium + kinetic reactions) via file I/O
    !> @details Reads aqueous component concentrations from a file after conservative transport,
    !>          solves the coupled equilibrium+kinetic step implicitly per target water
    !>          (theta=1 backward Euler, lumped kinetic mixing ratio; the initial state has zero
    !>          kinetic rates, so no previous-step term exists) with an explicit-Euler
    !>          fallback for non-converging cells, and writes the resulting concentrations and
    !>          end-of-step kinetic reaction rates to the output file.
        subroutine interfaz_comps_arch_eq_kin(this,path,num_aq_comps,file_in,Delta_t,file_out)
            import chemistry_c
            class(chemistry_c) :: this                    !< Chemistry object
            character(len=*), intent(in) :: path          !< Path for input and output files
            integer(kind=4), intent(in) :: num_aq_comps   !< Number of aqueous components
            character(len=*), intent(in) :: file_in       !< Input filename with post-transport concentrations
            real(kind=8), intent(in) :: Delta_t           !< Time step size for reactive chemistry
            character(len=*), intent(in) :: file_out      !< Output filename for post-reaction concentrations
        end subroutine

    !> @brief Interface for component concentrations (equilibrium reactions only) via file I/O
    !> @details Reads aqueous component concentrations from a file after conservative transport,
    !>          performs reactive chemistry step (eq only, no kinetics), and writes updated concentrations to output file.
        subroutine interfaz_comps_arch_eq(this,path,num_aq_comps,file_in,Delta_t,file_out)
            import chemistry_c
            class(chemistry_c) :: this                    !< Chemistry object
            character(len=*), intent(in) :: path          !< Path for input and output files
            integer(kind=4), intent(in) :: num_aq_comps   !< Number of aqueous components
            character(len=*), intent(in) :: file_in       !< Input filename with post-transport concentrations
            real(kind=8), intent(in) :: Delta_t           !< Time step size for reactive chemistry
            character(len=*), intent(in) :: file_out      !< Output filename for post-reaction concentrations
        end subroutine

    !> @brief Interface for component concentrations (equilibrium + kinetic reactions) via file I/O - transposed input
    !> @details Same as interfaz_comps_arch_eq_kin (implicit theta=1 backward-Euler solve with
    !>          explicit-Euler fallback, writes end-of-step kinetic reaction rates) but file_in has
    !>          rows=targets and columns=components.
        subroutine interfaz_comps_arch_eq_kin_T(this,path,num_aq_comps,file_in,Delta_t,file_out)
            import chemistry_c
            class(chemistry_c) :: this                    !< Chemistry object
            character(len=*), intent(in) :: path          !< Path for input and output files
            integer(kind=4), intent(in) :: num_aq_comps   !< Number of aqueous components
            character(len=*), intent(in) :: file_in       !< Input filename with post-transport concentrations (rows: targets, columns: components)
            real(kind=8), intent(in) :: Delta_t           !< Time step size for reactive chemistry
            character(len=*), intent(in) :: file_out      !< Output filename for post-reaction concentrations
        end subroutine

    !> @brief Interface for component concentrations (equilibrium reactions only) via file I/O - transposed input
    !> @details Same as interfaz_comps_arch_eq but file_in has rows=targets and columns=components.
        subroutine interfaz_comps_arch_eq_T(this,path,num_aq_comps,file_in,Delta_t,file_out)
            import chemistry_c
            class(chemistry_c) :: this                    !< Chemistry object
            character(len=*), intent(in) :: path          !< Path for input and output files
            integer(kind=4), intent(in) :: num_aq_comps   !< Number of aqueous components
            character(len=*), intent(in) :: file_in       !< Input filename with post-transport concentrations (rows: targets, columns: components)
            real(kind=8), intent(in) :: Delta_t           !< Time step size for reactive chemistry
            character(len=*), intent(in) :: file_out      !< Output filename for post-reaction concentrations
        end subroutine

    !> @brief Interface for component concentrations (equilibrium + optional kinetic) with internal conservative mixing
    !> @details Reads aqueous component concentrations PRIOR to mixing and the mixing ratios, performs the
    !>          conservative mixing internally and then the reactive chemistry step, and writes the updated
    !>          concentrations to the output file. file_in has rows=components and columns=waters; file_mix
    !>          stores the mixing ratios in row-major order (rows=target waters, columns=waters) and its
    !>          column ordering must match file_in.
        subroutine interfaz_comps_arch_eq_kin_mix(this,path,num_aq_comps,file_in,file_mix,Delta_t,file_out)
            import chemistry_c
            class(chemistry_c) :: this                    !< Chemistry object
            character(len=*), intent(in) :: path          !< Path for input and output files
            integer(kind=4), intent(in) :: num_aq_comps   !< Number of aqueous components
            character(len=*), intent(in) :: file_in       !< Input filename with concentrations PRIOR to mixing (rows: components, columns: waters)
            character(len=*), intent(in) :: file_mix      !< Input filename with mixing ratios in row-major order (rows: target waters, columns: waters)
            real(kind=8), intent(in) :: Delta_t           !< Time step size for reactive chemistry
            character(len=*), intent(in) :: file_out      !< Output filename for post-reaction concentrations
        end subroutine

    !> @brief Interface for aqueous species concentrations via file I/O (no equilibrium reactions)
    !> @details Reads aqueous species concentrations from a file after conservative transport,
    !>          performs reactive chemistry step (kinetic only), and writes updated concentrations to output file.
        subroutine interfaz_esp_arch(this,path,num_aq_comps,file_in,Delta_t,file_out)
            import chemistry_c
            class(chemistry_c) :: this                    !< Chemistry object
            character(len=*), intent(in) :: path          !< Path for input and output files
            integer(kind=4), intent(in) :: num_aq_comps   !< Number of aqueous species
            character(len=*), intent(in) :: file_in       !< Input filename with post-transport concentrations
            real(kind=8), intent(in) :: Delta_t           !< Time step
            character(len=*), intent(in) :: file_out      !< Output filename for post-reaction concentrations
        end subroutine
        
    !> @brief Interface for component concentrations via variables
    !> @details Takes aqueous component concentrations after conservative transport (u_mix),
    !>          performs reactive chemistry step, and returns updated concentrations (u_new).
    !> @param[in,out] this Chemistry object
    !> @param[in] u_mix Post-transport concentrations [mol/L]
    !> @param[in] Delta_t Time step size for reactive chemistry
    !> @param[out] u_new Post-reaction concentrations [mol/L]
            subroutine interfaz_comps_vars(this,u_mix,Delta_t,u_new)
            import chemistry_c
            class(chemistry_c) :: this                    !< Chemistry object
            real(kind=8), intent(in) :: u_mix(:)          !< Post-transport concentrations [mol/L]
            real(kind=8), intent(in) :: Delta_t           !< Time step size for reactive chemistry
            real(kind=8), intent(out) :: u_new(:)         !< Post-reaction concentrations [mol/L]
            end subroutine

    !> @brief Read target solid definitions from file
    !> @details Reads target solid chemistry definitions from input files,
    !>          associating solids with spatial locations and reactive zones.
    !> @param[in,out] this Chemistry object to populate with target solids
    !> @param[in] root Root name for input files
    !> @param[in] nsrz Number of solid reactive zones
    !> @param[in] ngrz Number of gas reactive zones
            subroutine read_tar_sol(this,root,nsrz,ngrz)
            import chemistry_c
            class(chemistry_c) :: this                    !< Chemistry object to populate
            character(len=*), intent(in) :: root          !< Root name for input files
            integer(kind=4), intent(in) :: nsrz           !< Number of solid reactive zones
            integer(kind=4), intent(in) :: ngrz           !< Number of gas reactive zones
            end subroutine
            
    !> @brief Read target gas definitions from file
    !> @details Reads target gas phase chemistry definitions from input files,
    !>          associating gas phases with spatial locations and reactive zones.
    !> @param[in,out] this Chemistry object to populate with target gases
    !> @param[in] root Root name for input files
    !> @param[in] ngrz Number of gas reactive zones
            subroutine read_tar_gas(this,root,ngrz)
            import chemistry_c
            class(chemistry_c) :: this                    !< Chemistry object to populate
            character(len=*), intent(in) :: root          !< Root name for input files
            integer(kind=4), intent(in) :: ngrz           !< Number of gas reactive zones
            end subroutine

    !> @brief Initialize complete chemistry system
    !> @details Reads chemical system, initial/boundary conditions, and reactive zone
    !>          definitions, then allocates and populates all chemistry data structures.
    !> @param[in,out] this Chemistry object to initialize
    !> @param[in] root Root name for input/output files
    !> @param[in] path_pb Path to problem directory
    !> @param[in] path_DB Path to database files
        subroutine initialise_chemistry(this,root,path_pb,path_DB)
            import chemistry_c                            !< Import chemistry class definition
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object to initialize
            character(len=*), intent(in) :: root          !< Root name for input/output files
            character(len=*), intent(in) :: path_pb       !< Path to problem directory
            character(len=*), intent(in) :: path_DB       !< Path to database files
        end subroutine

    !> @brief Single iteration of the ideal lumped reactive mixing solver
    !> @details Performs one iteration of the Water Mixing Algorithm using the ideal
    !>          activity model with spatial lumping. Processes a batch of target waters
    !>          and delegates to the provided procedure pointer for the actual
    !>          equilibrium-kinetic chemical reaction step.
    !> @param[in,out] this Chemistry object coordinating the simulation
    !> @param[in] k Current time step index
    !> @param[in] Delta_t Time step size
    !> @param[in] theta_r Integration method parameter for chemical reactions
    !> @param[in] num_tar_wat Number of target waters in this batch
    !> @param[in] ind_tar_wat Indices of target waters in this batch
    !> @param[in] n_p_cache Number of primary species per target water
    !> @param[in] n_v_cache Number of variable-activity species per target water
    !> @param[in] n_eq_cache Number of equilibrium reactions per target water
    !> @param[in] has_minerals_flag Flag indicating presence of minerals per target water
    !> @param[in] has_gases_flag Flag indicating presence of gases per target water
    !> @param[in] mix_conc_indices Mixing water indices for concentrations
    !> @param[in] mixing_ratios_conc Mixing ratios for concentrations
    !> @param[in] lumped_lambdas Lumped eigenvalues for spatial lumping
    !> @param[in] all_conc_old Old concentrations for all waters
    !> @param[in] all_ind_aq_sp Aqueous species indices for all waters
    !> @param[in,out] conc_comp Component concentrations (updated in place)
    !> @param[in] p_solver Procedure pointer to equilibrium-kinetic solver
        subroutine solve_reactive_mixing_ideal_lump_iter(this, k, Delta_t, theta_r, time, &
            num_tar_wat, ind_tar_wat, &
            n_p_cache, n_v_cache, n_eq_cache, &
            has_minerals_flag, has_gases_flag, &
            mix_conc_indices, mixing_ratios_conc, &
            mixing_ratios_R, &
            lumped_lambdas, &
            conc_comp, iunit_mb, iunit_celldbg, &
            p_solver, p_mixer)
            import chemistry_c                            !< Import chemistry class definition
            import real_array_c                           !< Import real array container type
            import int_array_c                            !< Import integer array container type
            import reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2 !< Import solver procedure interface
            import compute_c_mix_by_zones_lump
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object coordinating the simulation
            integer(kind=4), intent(in) :: k              !< Current time step index
            real(kind=8), intent(in) :: Delta_t           !< Time step size
            real(kind=8), intent(in) :: theta_r           !< Integration method parameter for chemical reactions
            real(kind=8), intent(in) :: time              !< Current simulation time (mass-balance header)
            integer(kind=4), intent(in) :: num_tar_wat    !< Number of target waters in this batch
            integer(kind=4), intent(in) :: ind_tar_wat(:) !< Indices of target waters in this batch
            integer(kind=4), intent(in) :: n_p_cache(:)   !< Number of primary species per target water
            integer(kind=4), intent(in) :: n_v_cache(:)   !< Number of variable-activity species per target water
            integer(kind=4), intent(in) :: n_eq_cache(:)  !< Number of equilibrium reactions per target water
            logical, intent(in) :: has_minerals_flag(:)   !< Flag indicating presence of minerals per target water
            logical, intent(in) :: has_gases_flag(:)      !< Flag indicating presence of gases per target water
            class(int_array_c), intent(in) :: mix_conc_indices !< Mixing water indices for concentrations
            class(real_array_c), intent(in) :: mixing_ratios_conc !< Mixing ratios for concentrations
            class(real_array_c), intent(in) :: mixing_ratios_R !< Reaction mixing ratios (mix_conc_indices layout)
            real(kind=8), intent(in) :: lumped_lambdas(:) !< Lumped eigenvalues for spatial lumping
            real(kind=8), intent(inout) :: conc_comp(:)   !< Component concentrations (updated in place)
            integer(kind=4), intent(in) :: iunit_mb       !< Mass-balance output file unit (already open)
            integer(kind=4), intent(in) :: iunit_celldbg  !< Per-cell [CELL-DBG] trace file unit (already open)
            procedure(reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2) :: p_solver !< Procedure pointer to equilibrium-kinetic solver
            procedure(compute_c_mix_by_zones_lump), pointer, intent(in) :: p_mixer
        end subroutine

    end interface
    
    contains
    
        !> @brief Set spatial lumping flag for Water Mixing Algorithm
        !> @details Controls whether spatial lumping is applied in the Water Mixing Algorithm
        !>          to reduce computational cost by grouping similar spatial locations
        !> @param[in] this Chemistry object
        !> @param[in] lump_flag Spatial lumping flag (.true.=apply lumping, .false.=full resolution)
        subroutine set_lump_flag(this,lump_flag)
            class(chemistry_c) :: this                    !< Chemistry object instance
            logical, intent(in) :: lump_flag              !< Spatial lumping control flag
            this%lump_flag=lump_flag                      !< Assign input flag to object attribute
        end subroutine
        
        !> @brief Set mixing flag for the Water Mixing Algorithm
        !> @details Controls whether waters are mixed per reactive zone (grouping the mixing
        !>          waters by their associated reactive zone, as in compute_c_mix_by_zones) or
        !>          directly over all mixing waters regardless of their reactive zone.
        !> @param[in] this Chemistry object
        !> @param[in] mix_by_zones_flag Mixing control flag (.true.=per reactive zone, .false.=direct)
        subroutine set_mix_by_zones_flag(this,mix_by_zones_flag)
            class(chemistry_c) :: this                    !< Chemistry object instance
            logical, intent(in) :: mix_by_zones_flag      !< Mixing control flag
            this%mix_by_zones_flag=mix_by_zones_flag      !< Assign input flag to object attribute
        end subroutine
        
        !> @brief Set reaction rate averaging method
        !> @details Controls how reaction rates are averaged in the Water Mixing Algorithm:
        !>          - 1: Concentration-based averaging
        !>          - 2: Direct rate averaging
        !> @param[in] this Chemistry object
        !> @param[in] rk_avg_opt Rate averaging option (1 or 2)
        subroutine set_rk_avg_opt(this,rk_avg_opt)
            class(chemistry_c) :: this                    !< Chemistry object instance
            integer(kind=4), intent(in) :: rk_avg_opt     !< Rate averaging option [1-2]
            if (rk_avg_opt<1 .or. rk_avg_opt>2) then                      !< Validate rate averaging option is in valid range [1-2]
                error stop "Chemistry attribute 'rk_avg_opt' must be 1 or 2"  !< Terminate program with error message for invalid input
            else                                                               !< Input is valid
                this%rk_avg_opt=rk_avg_opt                                     !< Assign validated rate averaging option to object attribute
            end if                                                             !< End validation block
        end subroutine

        !> @brief Set downstream reaction rate estimation method
        !> @details Controls the interpolation scheme for estimating reaction rates
        !>          at downstream locations in the Water Mixing Algorithm:
        !>          - 1, 2, 3, 4: Different interpolation schemes
        !> @param[in] this Chemistry object  
        !> @param[in] r_down_opt Downstream rate estimation option (1-4)
        subroutine set_r_down_opt(this,r_down_opt)
            class(chemistry_c) :: this                    !< Chemistry object instance
            integer(kind=4), intent(in) :: r_down_opt    !< Downstream rate option [1-4]
            if (r_down_opt<1 .or. r_down_opt>4) then                        !< Validate downstream rate option is in valid range [1-4]
                error stop "Chemistry attribute 'r_down_opt' must be 1, 2, 3 or 4" !< Terminate program with error message for invalid input
            else                                                                   !< Input is valid
                this%r_down_opt=r_down_opt                                       !< Assign validated downstream rate option to object attribute
            end if                                                                 !< End validation block
        end subroutine
        
        !> @brief Set chemical data input format option
        !> @details Specifies the format for reading chemical input data:
        !>          - 1: CHEPROO format (fully implemented)
        !>          - 2: PHREEQC format (not fully implemented)
        !>          - 3: PFLOTRAN format (not implemented)
        !> @param[in] this Chemistry object
        !> @param[in] option Input format option (currently only 1 is supported)
        subroutine set_read_opt(this,option)
            class(chemistry_c) :: this                    !< Chemistry object instance
            integer(kind=4), intent(in) :: option         !< Input format option [1-3]
            if (option<1 .or. option>3) then                                    !< Check if option is outside implemented range [1-3]
                error stop "Chemistry input option not implemented yet"            !< Terminate for unimplemented option
            else if (option>1) then                                              !< Check if option is PHREEQC (2) or PFLOTRAN (3)
                error stop "Chemistry input option not fully implemented yet"      !< Terminate for partially implemented options
            else                                                                 !< Option is 1 (CHEPROO format - fully supported)
                this%read_opt=option                                             !< Assign validated input format option to object attribute
            end if                                                               !< End validation block
        end subroutine
        
        !> @brief Set Jacobian computation method
        !> @details Controls how Jacobian matrices are computed for Newton-Raphson methods:
        !>          - 0: Incremental coefficients (finite differences)
        !>          - 1: Analytical Jacobian (exact derivatives)
        !> @param[in] this Chemistry object
        !> @param[in] Jac_opt Jacobian computation option (0 or 1)
        subroutine set_Jac_opt(this,Jac_opt)
            implicit none
            class(chemistry_c) :: this                    !< Chemistry object instance
            integer(kind=4), intent(in) :: Jac_opt        !< Jacobian computation method [0-1]
            if (Jac_opt<0 .or. Jac_opt>1) then                        !< Validate Jacobian option is in valid range [0-1]
                error stop "Chemistry attribute 'Jac_opt' must be 0 or 1" !< Terminate program with error message for invalid input
            else                                                           !< Input is valid
                this%Jac_opt=Jac_opt                                       !< Assign validated Jacobian computation option to object attribute
            end if                                                         !< End validation block
        end subroutine
        
        !> @brief Set complete chemical system
        !> @details Assigns a complete chemical system object containing all species,
        !>          reactions, thermodynamic data, and stoichiometric matrices
        !> @param[in] this Chemistry object
        !> @param[in] chem_syst_obj Complete chemical system object
        subroutine set_chem_syst(this,chem_syst_obj)
            implicit none
            class(chemistry_c) :: this                        !< Chemistry object instance
            type(chem_system_c), intent(in) :: chem_syst_obj  !< Chemical system object to assign
            this%chem_syst=chem_syst_obj                       !< Copy complete chemical system to object attribute
        end subroutine
        
        !> @brief Set number of water types
        !> @details Sets the total number of distinct water composition types,
        !>          including initial conditions and boundary water types
        !> @param[in] this Chemistry object
        !> @param[in] num_wat_types Number of water types (must be > 0)
        subroutine set_num_wat_types(this,num_wat_types)
            implicit none
            class(chemistry_c) :: this                        !< Chemistry object instance
            integer(kind=4), intent(in) :: num_wat_types      !< Number of water types [>0]
            if (num_wat_types<1) then                                               !< Validate number of water types is positive
                error stop "Chemistry attribute 'num_wat_types' must be greater than 0" !< Terminate program if invalid (must have at least 1 water type)
            else                                                                        !< Input is valid
                this%num_wat_types=num_wat_types                                        !< Assign validated number of water types to object attribute
            end if                                                                      !< End validation block
        end subroutine
        
        !> @brief Set number of target waters
        !> @details Sets the number of target waters within the computational domain,
        !>          excluding boundary waters. If not provided, calculates automatically
        !>          as total target waters minus external waters
        !> @param[in] this Chemistry object
        !> @param[in] num_target_waters Number of domain target waters (optional)
        subroutine set_num_target_waters(this,num_target_waters)
            implicit none
            class(chemistry_c) :: this                           !< Chemistry object instance
            integer(kind=4), intent(in), optional :: num_target_waters !< Number of domain target waters
            if (present(num_target_waters)) then                                                     !< Check if optional parameter was provided
                this%num_target_waters=num_target_waters                                             !< Use explicitly provided value for domain target waters count
            else                                                                                       !< Optional parameter not provided
                this%num_target_waters=this%num_waters-this%num_ext_waters                 !< Calculate domain waters by subtracting external waters from total
            end if                                                                                     !< End conditional assignment
        end subroutine
        
        !> @brief Set number of target solids
        !> @details Sets the number of target solid phases at spatial locations
        !> @param[in] this Chemistry object
        !> @param[in] num_target_solids Number of target solids
        subroutine set_num_target_solids(this,num_target_solids)
            implicit none
            class(chemistry_c) :: this                           !< Chemistry object instance
            integer(kind=4), intent(in) :: num_target_solids     !< Number of target solids to set
            this%num_target_solids=num_target_solids             !< Assign number of target solids to object attribute
        end subroutine
        
        !> @brief Set number of external (boundary) waters
        !> @details Sets the total number of external water types at boundaries,
        !>          used for boundary condition specification in reactive transport
        !> @param[in] this Chemistry object instance
        !> @param[in] num_ext_waters Number of external water types
        subroutine set_num_ext_waters(this,num_ext_waters)
            implicit none
            class(chemistry_c) :: this               !< Chemistry object instance
            integer(kind=4), intent(in) :: num_ext_waters !< Number of external water types
            this%num_ext_waters=num_ext_waters   !< Assign external waters count
        end subroutine
        
        !> @brief Set number of boundary water types
        !> @details Sets the number of distinct boundary water types for boundary conditions.
        !>          These are typically Dirichlet or flux boundary waters
        !> @param[in] this Chemistry object instance
        !> @param[in] num_bd_waters Number of boundary water types
        subroutine set_num_bd_waters(this,num_bd_waters)
            implicit none
            class(chemistry_c) :: this                !< Chemistry object instance
            integer(kind=4), intent(in) :: num_bd_waters !< Number of boundary waters
            this%num_bd_waters=num_bd_waters      !< Assign boundary waters count
        end subroutine
        
        !> @brief Set number of recharge water types
        !> @details Sets the number of distinct recharge (infiltration) water types,
        !>          used for modeling groundwater recharge with varying chemistry
        !> @param[in] this Chemistry object instance
        !> @param[in] num_rech_waters Number of recharge water types
        subroutine set_num_rech_waters(this,num_rech_waters)
            implicit none
            class(chemistry_c) :: this                    !< Chemistry object instance
            integer(kind=4), intent(in) :: num_rech_waters !< Number of recharge water types
            this%num_rech_waters=num_rech_waters          !< Assign recharge waters count
        end subroutine

        !> @brief Set array of target solids
        !> @details Assigns the array of target solid chemistry objects and updates count.
        !>          Each target solid represents mineralogical composition at a spatial location
        !> @param[in] this Chemistry object instance
        !> @param[in] target_solids Array of solid chemistry objects
        subroutine set_target_solids(this,target_solids)
            implicit none
            class(chemistry_c) :: this                       !< Chemistry object instance
            type(solid_chemistry_c), intent(in) :: target_solids(:) !< Array of solid chemistry objects
            this%target_solids=target_solids                 !< Perform deep copy of target solids array to object attribute
            this%num_target_solids=size(target_solids)       !< Update count by querying size of input array
        end subroutine
        
        !> @brief Set array of target gases
        !> @details Assigns the array of target gas chemistry objects and updates count.
        !>          Each target gas represents gas phase composition at a spatial location
        !> @param[in] this Chemistry object instance
        !> @param[in] target_gases Array of gas chemistry objects
        subroutine set_target_gases(this,target_gases)
            implicit none
            class(chemistry_c) :: this                     !< Chemistry object instance
            type(gas_chemistry_c), intent(in) :: target_gases(:) !< Array of gas chemistry objects
            this%target_gases=target_gases                 !< Perform deep copy of target gases array to object attribute
            this%num_target_gases=size(target_gases)       !< Update count by querying size of input array
        end subroutine
        
       !> @brief Allocate target waters array
       !> @details Allocates memory for the waters array. If already allocated,
       !>          deallocates first to prevent memory leaks. Optionally updates count
       !> @param[in] this Chemistry object instance
       !> @param[in] num_tar_wat Optional number of target waters (updates attribute if present)
       subroutine allocate_waters(this,num_tar_wat)
            implicit none
            class(chemistry_c) :: this                       !< Chemistry object instance
            integer(kind=4), intent(in), optional :: num_tar_wat !< Optional target waters count
            if (present(num_tar_wat)) then                   !< Check if optional count parameter was provided by caller
                this%num_waters=num_tar_wat           !< Update object attribute with provided count value
            end if                                           !< Otherwise use existing count value
            if (allocated(this%waters)) then          !< Check allocation status to prevent double allocation
                deallocate(this%waters)               !< Free existing memory to prevent memory leak
            end if                                           !< Array now safely deallocated or was never allocated
            allocate(this%waters(this%num_waters)) !< Allocate array with size from object attribute
        end subroutine
        
        !> @brief Allocate domain target water indices arrays
        !> @details Allocates memory for indices of target waters within the computational domain
        !>          (excluding boundary waters). Creates both current and initial index arrays.
        !>          Deallocates first if already allocated to prevent memory leaks
        !> @param[in] this Chemistry object instance
        !> @param[in] num_tar_wat Optional number of domain target waters
        subroutine allocate_tar_wat_indices(this,num_tar_wat)
            implicit none
            class(chemistry_c) :: this                           !< Chemistry object instance
            integer(kind=4), intent(in), optional :: num_tar_wat !< Optional domain target waters count
            if (present(num_tar_wat)) then                   !< Check if optional count parameter was provided by caller
                this%num_target_waters=num_tar_wat       !< Update object attribute with provided domain waters count
            end if                                               !< Otherwise use existing count value
            if (allocated(this%tar_wat_indices)) then        !< Check allocation status of current indices array
                deallocate(this%tar_wat_indices)             !< Free existing memory for current indices
            end if                                               !< Current indices array now safely deallocated
            !> Allocate both current and initial domain indices arrays with same size
            allocate(this%tar_wat_indices(this%num_target_waters))
        end subroutine

        !> @brief Allocate initial domain target water indices array
        !> @details Allocates memory for the initial indices of target waters within the
        !>          computational domain (excluding boundaries). Used for restart and mass balance.
        !>          Deallocates first if already allocated to prevent memory leaks.
        !> @param[in,out] this Chemistry object instance
        !> @param[in] num_target_waters_init Optional number of initial domain target waters
        subroutine allocate_tar_wat_indices_init(this,num_target_waters_init)
            implicit none
            class(chemistry_c) :: this                           !< Chemistry object instance
            integer(kind=4), intent(in), optional :: num_target_waters_init !< Optional domain target waters count
            if (present(num_target_waters_init)) then                   !< Check if optional count parameter was provided by caller
                this%num_target_waters_init=num_target_waters_init       !< Update object attribute with provided domain waters count
            end if                                               !< Otherwise use existing count value
            if (allocated(this%tar_wat_indices_init)) then        !< Check allocation status of current indices array
                deallocate(this%tar_wat_indices_init)             !< Free existing memory for current indices
            end if                                               !< Current indices array now safely deallocated
            !> Allocate both current and initial domain indices arrays with same size
            allocate(this%tar_wat_indices_init(this%num_target_waters_init))
        end subroutine

        !> @brief Allocate external waters indices array
        !> @details Allocates memory for indices of external (boundary) water types.
        !>          External waters are used to impose boundary conditions in reactive transport.
        !>          Deallocates first if already allocated to prevent memory leaks
        !> @param[in] this Chemistry object instance
        !> @param[in] num_ext_wat Optional number of external waters
        subroutine allocate_ext_waters_indices(this,num_ext_wat)
            implicit none
            class(chemistry_c) :: this                       !< Chemistry object instance
            integer(kind=4), intent(in), optional :: num_ext_wat !< Optional external waters count
            if (present(num_ext_wat)) then                   !< Check if optional count parameter was provided by caller
                this%num_ext_waters=num_ext_wat              !< Update object attribute with provided external waters count
            end if                                           !< Otherwise use existing count value
            if (allocated(this%ext_waters_indices)) then     !< Check allocation status to prevent double allocation
                deallocate(this%ext_waters_indices)          !< Free existing memory to prevent memory leak
            end if                                           !< Array now safely deallocated or was never allocated
            allocate(this%ext_waters_indices(this%num_ext_waters)) !< Allocate array with size from object attribute
        end subroutine
        
        !> @brief Allocate target solids arrays
        !> @details Allocates memory for target_solids and target_solids_init arrays.
        !>          Target solids store mineralogical composition at each spatial location.
        !>          Deallocates first if already allocated to prevent memory leaks
        !> @param[in] this Chemistry object instance
        !> @param[in] n Optional number of target solids
        subroutine allocate_target_solids(this,n)
            implicit none
            class(chemistry_c) :: this                !< Chemistry object instance
            integer(kind=4), intent(in), optional :: n !< Optional target solids count
            if (present(n)) then                      !< Check if optional count parameter was provided by caller
                this%num_target_solids=n              !< Update object attribute with provided solids count
            end if                                    !< Otherwise use existing count value
            if (allocated(this%target_solids)) then   !< Check allocation status of current solids array
                deallocate(this%target_solids,this%target_solids_init) !< Free memory for both current and initial arrays simultaneously
            end if                                    !< Both arrays now safely deallocated or were never allocated
            !> Allocate both current and initial target solids arrays with same size
            allocate(this%target_solids(this%num_target_solids),this%target_solids_init(this%num_target_solids)) !< Single statement allocates both arrays for efficiency
        end subroutine

        subroutine allocate_target_solids_dummy(this)
            !> Allocate dummy array for target solids
            !> Used when no target solids are defined, but some routines need to access this attribute
            implicit none
            class(chemistry_c) :: this                                !< Chemistry object instance
            if (allocated(this%target_solids_dummy)) then             !< Check if dummy array already allocated
                deallocate(this%target_solids_dummy)                  !< Free existing memory to prevent leak
            end if                                                    !< Array now safely deallocated
            this%num_target_solids_dummy=this%num_gas_zones+1 !< Set dummy solids count to num_gas_zones+1 to ensure non-zero dimension even with no gas zones
            allocate(this%target_solids_dummy(this%num_target_solids_dummy)) !< Allocate with size=num_gas_zones+1 to ensure non-zero dimension even with no gas zones
        end subroutine

        subroutine allocate_reactive_zones_dummy(this)
            !> Allocate dummy array for target solids
            !> Used when no target solids are defined, but some routines need to access this attribute
            implicit none
            class(chemistry_c) :: this                                !< Chemistry object instance
            if (allocated(this%react_zones_dummy)) then             !< Check if dummy array already allocated
                deallocate(this%react_zones_dummy)                  !< Free existing memory to prevent leak
            end if                                                    !< Array now safely deallocated
            allocate(this%react_zones_dummy(this%num_target_solids_dummy)) !< Allocate with size=num_gas_zones+1 to ensure non-zero dimension even with no gas zones
        end subroutine
        
        
        !> @brief Allocate target gases arrays
        !> @details Allocates memory for target_gases and target_gases_init arrays.
        !>          Target gases store gas phase composition at each spatial location.
        !>          Does not deallocate first (assumes first allocation)
        !> @param[in] this Chemistry object instance
        !> @param[in] n Optional number of target gases
        subroutine allocate_target_gases(this,n)
            implicit none
            class(chemistry_c) :: this                !< Chemistry object instance
            integer(kind=4), intent(in), optional :: n !< Optional target gases count
            if (present(n)) then                      !< Check if optional count parameter was provided by caller
                this%num_target_gases=n               !< Update object attribute with provided gases count
            end if                                    !< Otherwise use existing count value
            !> Allocate both current and initial target gases arrays with same size
            allocate(this%target_gases(this%num_target_gases),this%target_gases_init(this%num_target_gases)) !< Single statement allocates both arrays for efficiency, no deallocation check (assumes first allocation)
        end subroutine
        
        !> @brief Allocate recharge waters indices array
        !> @details Allocates memory for indices of recharge (infiltration) water types.
        !>          Recharge waters represent time-varying infiltration chemistry
        !> @param[in] this Chemistry object instance
        !> @param[in] num_rech_waters Optional number of recharge water types
        subroutine allocate_rech_waters_indices(this,num_rech_waters)
            implicit none
            class(chemistry_c) :: this                           !< Chemistry object instance
            integer(kind=4), intent(in), optional :: num_rech_waters !< Optional recharge waters count
            if (present(num_rech_waters)) then                   !< Check if optional count parameter was provided by caller
                this%num_rech_waters=num_rech_waters             !< Update object attribute with provided recharge waters count
            end if                                               !< Otherwise use existing count value
            allocate(this%rech_waters_indices(this%num_rech_waters)) !< Allocate integer array to store recharge water indices with size from object attribute
        end subroutine
        
        !> @brief Allocate boundary waters indices arrays
        !> @details Allocates memory for indices of boundary water types used in boundary conditions.
        !>          Creates both current and initial index arrays for temporal BC tracking
        !> @param[in] this Chemistry object instance
        !> @param[in] num_bd_waters Optional number of boundary water types
        subroutine allocate_bd_waters_indices(this,num_bd_waters)
            implicit none                                        !< Require explicit variable declarations
            class(chemistry_c) :: this                           !< Chemistry object instance
            integer(kind=4), intent(in) :: num_bd_waters         !< Boundary waters count
            if (mod(num_bd_waters,2) .ne. 0 .or. num_bd_waters <= 0) then !< Validate even positive integer
                error stop "Chemistry attribute 'num_bd_waters' must be an even positive integer" !< Terminate for invalid input
            end if                                               !< End validation block
                this%num_bd_waters=num_bd_waters                 !< Update object attribute with provided boundary waters count
            allocate(this%bd_waters_indices(this%num_bd_waters))  !< Allocate boundary waters indices array
            this%bd_waters_indices=0 !< Initialize current indices array to zero
        end subroutine
        
        !> @brief Set array of target waters
        !> @details Assigns the array of target aqueous chemistry objects and updates count.
        !>          Target waters represent aqueous chemistry at each spatial location
        !> @param[in] this Chemistry object instance
        !> @param[in] waters Array of aqueous chemistry objects
        subroutine set_waters(this,waters)
            implicit none
            class(chemistry_c) :: this                            !< Chemistry object instance
            class(aqueous_chemistry_c), intent(in) :: waters(:) !< Array of aqueous chemistry objects
            this%waters=waters                      !< Perform deep copy of entire target waters array to object attribute
            this%num_waters=size(waters)            !< Automatically update count by querying size of input array
        end subroutine
        
       !> @brief Allocate reactive zones array
       !> @details Allocates memory for the reactive_zones array. Reactive zones group
       !>          spatial locations with identical chemical system and reaction network.
       !>          Deallocates first if already allocated to prevent memory leaks
       !> @param[in] this Chemistry object instance
       !> @param[in] num_reactive_zones Optional number of reactive zones
       subroutine allocate_reactive_zones(this,num_reactive_zones)
            implicit none
            class(chemistry_c) :: this                               !< Chemistry object instance
            integer(kind=4), intent(in), optional :: num_reactive_zones !< Optional reactive zones count
            if (present(num_reactive_zones)) then                    !< Check if optional count parameter was provided by caller
                this%num_reactive_zones=num_reactive_zones           !< Update object attribute with provided reactive zones count
            end if                                                   !< Otherwise use existing count value
            if (allocated(this%reactive_zones)) then                 !< Check if already allocated
                deallocate(this%reactive_zones)                      !< Deallocate to prevent memory leak
            end if
            allocate(this%reactive_zones(this%num_reactive_zones))   !< Allocate with current count
       end subroutine
       
       !> @brief Allocate mineral zones array
       !> @details Allocates memory for the mineral_zones array. Mineral zones group
       !>          spatial locations with identical mineralogical composition and reactions.
       !>          Deallocates first if already allocated to prevent memory leaks
       !> @param[in] this Chemistry object instance
       !> @param[in] num_mineral_zones Optional number of mineral zones
       subroutine allocate_mineral_zones(this,num_mineral_zones)
            implicit none
            class(chemistry_c) :: this                              !< Chemistry object instance
            integer(kind=4), intent(in), optional :: num_mineral_zones !< Optional mineral zones count
            if (present(num_mineral_zones)) then                    !< Check if count parameter provided
                this%num_mineral_zones=num_mineral_zones            !< Update mineral zones count
            end if
            if (allocated(this%mineral_zones)) then                 !< Check if already allocated
                deallocate(this%mineral_zones)                      !< Deallocate to prevent memory leak
            end if
            allocate(this%mineral_zones(this%num_mineral_zones))    !< Allocate with current count
       end subroutine
       
       !> @brief Set reactive zones array with automatic detection
       !> @details Sets reactive zones either from provided array or automatically detects
       !>          unique reactive zones by comparing target solids/gases chemistry.
       !>          Uses pairwise comparison algorithm to identify unique reactive zones
       !> @param[in] this Chemistry object instance
       !> @param[in] reactive_zones Optional array of reactive zone objects (if provided, directly assign)
       subroutine set_reactive_zones(this,reactive_zones)
            implicit none
            class(chemistry_c) :: this                              !< Chemistry object instance
            type(reactive_zone_c), intent(in), optional :: reactive_zones(:) !< Optional reactive zones array
            integer(kind=4) :: i,j,l                                !< Loop indices: i,j for comparison, l for unique zones
            logical :: flag                                         !< Comparison flag (true if zones match)
            integer(kind=4), allocatable :: rz_indices(:)           !< Marker array for unique reactive zones
            if (present(reactive_zones)) then                       !< Direct assignment if provided
                this%reactive_zones=reactive_zones                  !< Assign reactive zones array
                this%num_reactive_zones=size(reactive_zones)        !< Update count from array size
            else if (this%num_target_solids>0) then                 !< Auto-detect from target solids
                i=1                                                 !< Initialize first comparison index
                j=2                                                 !< Initialize second comparison index
                this%num_reactive_zones=1                           !< Start with at least one zone
                allocate(rz_indices(this%num_target_solids))        !< Allocate marker array
                rz_indices=0                                        !< Initialize all markers to zero
                do                                                  !< Pairwise comparison loop
                        !> Compare reactive zones of target solids i and j
                        call compare_react_zones(this%target_solids(i)%reactive_zone,this%target_solids(j)%reactive_zone,flag)
                        if (flag.eqv..true.) then                   !< Zones match - continue comparison
                            if (i<this%num_target_solids-1) then    !< More elements to compare
                                i=i+1                               !< Advance to next pair
                                j=i+1
                            else                                    !< Reached end of array
                                exit                                !< Exit comparison loop
                            end if
                        else if (j<this%num_target_solids) then     !< Zones differ - try next j
                            j=j+1                                   !< Increment second index
                        else if (i<this%num_target_solids-1) then   !< No match found for i
                            this%num_reactive_zones=this%num_reactive_zones+1 !< Increment unique zone count
                            rz_indices(i)=1                         !< Mark this solid as unique zone
                            i=i+1                                   !< Advance to next comparison
                            j=i+1
                        else                                        !< Final element reached
                            exit                                    !< Exit comparison loop
                        end if
                end do
                rz_indices(j)=1                                     !< Mark last unique solid
                call this%allocate_reactive_zones()                 !< Allocate with determined count
                l=1                                                 !< Initialize unique zone index
                do i=1,this%num_target_solids                       !< Assign unique reactive zones
                    if (rz_indices(i)==1) then                      !< This solid marks a unique zone
                        call this%reactive_zones(l)%copy_react_zone(this%target_solids(i)%reactive_zone) !< Copy zone
                        l=l+1                                       !< Increment zone counter
                    end if
                end do
            else if (this%num_target_gases>0) then                  !< Auto-detect from target gases
                i=1                                                 !< Initialize first comparison index
                j=2                                                 !< Initialize second comparison index
                this%num_reactive_zones=1                           !< Start with at least one zone
                allocate(rz_indices(this%num_target_gases))         !< Allocate marker array for gases
                rz_indices=0                                        !< Initialize all markers to zero
                do                                                  !< Pairwise comparison loop for gases
                        !> Compare reactive zones of target gases i and j
                        call compare_react_zones(this%target_gases(i)%reactive_zone,this%target_gases(j)%reactive_zone,flag)
                        if (flag.eqv..true.) then                   !< Zones match - continue comparison
                            if (i<this%num_target_gases-1) then     !< More elements to compare
                                i=i+1                               !< Advance to next pair
                                j=i+1
                            else                                    !< Reached end of array
                                exit                                !< Exit comparison loop
                            end if
                        else if (j<this%num_target_gases) then      !< Zones differ - try next j
                            j=j+1                                   !< Increment second index
                        else if (i<this%num_target_gases-1) then    !< No match found for i
                            this%num_reactive_zones=this%num_reactive_zones+1 !< Increment unique zone count
                            rz_indices(i)=1                         !< Mark this gas as unique zone
                            i=i+1                                   !< Advance to next comparison
                            j=i+1
                        else                                        !< Final element reached
                            exit                                    !< Exit comparison loop
                        end if
                end do
                rz_indices(j)=1                                     !< Mark last unique gas
                call this%allocate_reactive_zones()                 !< Allocate with determined count
                l=1                                                 !< Initialize unique zone index
                do i=1,this%num_target_gases                        !< Assign unique reactive zones from gases
                    if (rz_indices(i)==1) then                      !< This gas marks a unique zone
                        this%reactive_zones(l)=this%target_gases(i)%reactive_zone !< Copy gas reactive zone
                        l=l+1                                       !< Increment zone counter
                    end if
                end do
            else                                                    !< No solids or gases - empty system
                call this%allocate_reactive_zones(0)                !< Allocate zero reactive zones
            end if
            !> The zone id is its position in this%reactive_zones (used to index the array)
            do i=1,this%num_reactive_zones
                this%reactive_zones(i)%id=i
            end do
       end subroutine
       


        !> @brief Read and initialize a target water from file line data
        !> @details Similar to loop_read_tar_wat_init but used when reading from files with
        !>          pre-defined target solid and gas indices. Performs the following operations:
        !>          1. Validates target water index bounds
        !>          2. Assigns water type to target water and initial target water
        !>          3. Links solid and gas chemistry using provided indices
        !>          4. Detects new reactive zones and computes speciation algorithm
        !>          5. Handles species swapping for numerical stability
        !>          6. Computes stoichiometric matrices and reaction rates
        !>          Used in file-based initialization where chemistry objects are pre-built
        !> @param[in,out] this Chemistry object instance
        !> @param[in] flag_ext Flag indicating if target water is external
        !> @param[in] iszn Solid zone index
        !> @param[in] igzn Gas reactive zone index
        !> @param[in] tar_wat_ind Target water index being initialized
        !> @param[in] wtype Water type index
        !> @param[in] tar_sol_ind Target solid index (0 if no solid, uses dummy)
        !> @param[in] tar_gas_ind Target gas index (0 if no gas)
        !> @param[in,out] aux_iszn Auxiliary solid zone index from previous iteration
        !> @param[in,out] aux_igzn Auxiliary gas zone index from previous iteration
        subroutine read_tar_wat_line(this,flag_ext,iszn,igzn,tar_wat_ind,wtype,tar_sol_ind,tar_gas_ind,aux_iszn,aux_igzn)
            class(chemistry_c) :: this                    !< Chemistry object instance
            logical, intent(in) :: flag_ext               !< Flag to indicate if water is external
            integer(kind=4), intent(in) :: iszn           !< Index of solid zone
            integer(kind=4), intent(in) :: igzn           !< Index of gas zone
            integer(kind=4), intent(in) :: tar_wat_ind    !< Target water index
            integer(kind=4), intent(in) :: wtype          !< Water type index
            integer(kind=4), intent(in) :: tar_sol_ind    !< Target solid index
            integer(kind=4), intent(in) :: tar_gas_ind    !< Target gas index
            integer(kind=4), intent(inout) :: aux_iszn    !< Auxiliary solid zone index from previous iteration
            integer(kind=4), intent(inout) :: aux_igzn    !< Auxiliary gas zone index from previous iteration

            integer(kind=4) :: n_p_aq                     !< Number of aqueous primary species
            integer(kind=4) :: n2v_aq                     !< Number of secondary aqueous variable activity species
            integer(kind=4) :: j                          !< Loop index for swap detection
            logical :: flag_Se                            !< Flag to swap species indices for numerical stability
            integer(kind=4), allocatable :: swap(:),aux_swap(:) !< Indices of species to swap
            real(kind=8), allocatable :: rk(:)            !< Reaction rates array (used for initialization)
            type(solid_chemistry_c), target :: aux_solid_chem !< Auxiliary solid chemistry object (unused)
            
            
            allocate(swap(2),aux_swap(2))                 !< TEMPORARY FIX: allocate swap arrays for species index swapping

            !> Validate target water index bounds
            !> Validate target water index bounds
            if (tar_wat_ind<1 .or. tar_wat_ind>this%num_waters) then
                error stop "Target water index out of bounds"
            end if

            !> Record the water type per water so set_wat_indices_zones can group
            !> mineral-free external/boundary waters by shared water type.
            if (.not. allocated(this%wat_type_indices)) allocate(this%wat_type_indices(this%num_waters))
            this%wat_type_indices(tar_wat_ind)=wtype

            !> STEP 1: Assign aqueous chemistry from water type to current and initial target waters
            !> Copy the composition (species concentrations, activities, etc.) from the specified water type
            call this%waters(tar_wat_ind)%copy_aq_chem(this%wat_types(wtype)) !< set current target water composition
            call this%waters_init(tar_wat_ind)%copy_aq_chem(this%wat_types(wtype)) !< set initial target water composition
            !> Set unique identifiers for target waters
            call this%waters(tar_wat_ind)%set_id(tar_wat_ind) !< assign unique ID to current target water
            call this%waters_init(tar_wat_ind)%set_id(tar_wat_ind) !< assign unique ID to initial target water
            
            !> STEP 2: Link solid and gas chemistry based on provided indices
            !> This establishes the multiphase chemistry coupling for reactive transport
            if (tar_sol_ind>0) then !< Case A: We have a target solid (mineral phase present)
                if (tar_gas_ind>0) then !< Case A1: Both solid and gas phases present
                    !> Link gas chemistry to both current and initial target waters
                    call this%waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_gas_ind)) !< link current water to gas chemistry
                    call this%waters_init(tar_wat_ind)%set_gas_chemistry(this%target_gases_init(tar_gas_ind)) !< link initial water to gas chemistry
                !else !< Case A2: Solid phase only (no gas)
                end if
                !> Link solid chemistry to both current and initial target waters
                call this%waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(tar_sol_ind)) !< link current water to solid chemistry
                call this%waters(tar_wat_ind)%set_solid_chemistry_old(this%target_solids(tar_sol_ind)) !< link current water to solid chemistry
                call this%waters_init(tar_wat_ind)%set_solid_chemistry(this%target_solids_init(tar_sol_ind)) !< link initial water to solid chemistry
            else !< Case B: No solid zone (tar_sol_ind <= 0)
                if (tar_gas_ind>0) then !< Case B1: Gas phase only (no solid minerals)
                    !> Link gas chemistry to both current and initial target waters
                    call this%waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_gas_ind))
                    call this%waters_init(tar_wat_ind)%set_gas_chemistry(this%target_gases_init(tar_gas_ind))
                    !> TEMPORARY FIX: Use dummy solid chemistry object with reactive zone from gas
                    !> This maintains consistency in code structure even without actual mineral phase
                    call this%target_solids_dummy(1+igzn)%set_reactive_zone(this%target_gases(tar_gas_ind)%reactive_zone) !< link dummy solid to gas reactive zone
                    call this%waters(tar_wat_ind)%set_solid_chemistry(this%target_solids_dummy(1+igzn)) !< set dummy solid for current water
                    call this%waters(tar_wat_ind)%set_solid_chemistry_old(this%target_solids_dummy(1+igzn)) !< set dummy solid for current water
                    call this%waters_init(tar_wat_ind)%set_solid_chemistry(this%target_solids_dummy(1+igzn)) !< set default dummy for initial water
                else
                    !> Case B2: No target solid or gas phase
                    if (flag_ext .and. associated(this%wat_types(wtype)%solid_chemistry)) then
                        !> External (boundary/recharge) water: keep the solid chemistry of its WATER TYPE
                        !> (sol_types entry, whose reactive zone is the water type's react_zones_wat_types
                        !> zone with its minerals/components), instead of the mineral-free dummy. This is
                        !> the zone the mixer expects for inflow waters, and the component basis in which
                        !> the boundary water is reported.
                        call this%waters(tar_wat_ind)%set_solid_chemistry(this%wat_types(wtype)%solid_chemistry)
                        call this%waters(tar_wat_ind)%set_solid_chemistry_old(this%wat_types(wtype)%solid_chemistry)
                        call this%waters_init(tar_wat_ind)%set_solid_chemistry(this%wat_types(wtype)%solid_chemistry)
                    else
                        !> Domain water without solids (pure aqueous system): dummy solid chemistry
                        call this%waters(tar_wat_ind)%set_solid_chemistry(this%target_solids_dummy(1)) !< set dummy solid for current water
                        call this%waters(tar_wat_ind)%set_solid_chemistry_old(this%target_solids_dummy(1)) !< set dummy solid for current water
                        call this%waters_init(tar_wat_ind)%set_solid_chemistry(this%target_solids_dummy(1)) !< set default dummy for initial water
                    end if
                end if
            end if
            !> STEP 2b: The domain cell's cation-exchange equilibrium is (re)solved AFTER the
            !> per-water logK is built (STEP 4b below), once the now-known CEC is available.
            !> The water type's exchange was speciated before the surface-adsorption zone was
            !> read (CEC=0) and so cannot seed a physical free site here; solving it below with
            !> the merged reactive zone's CEC gives a free site and adsorbed cations that satisfy
            !> the site balance (conc_free + sum z*conc_ads = CEC).
            n_p_aq=this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species !< number of primary aqueous species
            n2v_aq=this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species !< number of secondary aqueous variable activity species
            !> STEP 3: Detect reactive zone changes and compute speciation algorithm arrays
            !> This step checks if the current target water has a different reactive zone than
            !> the previous one. Reactive zones group target waters with identical chemistry
            !> to optimize computations. If a new zone is detected, speciation arrays must be recomputed.
            if (n_p_aq>0) then !< Proceed only if there are aqueous species to process
                !> Check if reactive zone has changed from previous target water
                !> Reactive zone change detection based on solid zone AND gas zone indices
                if (aux_iszn==0 .or. aux_iszn/=iszn .or. aux_igzn/=igzn) then !< NEW REACTIVE ZONE DETECTED
                    !> When reactive zone changes, must recalculate species indices and speciation arrays
                    
                    !> Set species indices for primary, secondary, and variable activity species
                    call this%waters(tar_wat_ind)%set_ind_species()          !< set indices for current target water
                    call this%waters_init(tar_wat_ind)%set_ind_species()     !< set indices for initial target water
                    
                    !> NOTE: Do NOT call set_ind_var_act_species() here.
                    !> Multiple reactive zones may share the same object through pointer aliasing.
                    !> compute_speciation_alg_arrays may swap stoich_mat columns in-place for the
                    !> first water, and subsequent waters see the already-swapped matrix (flag_Se=F).
                    !> Preserving ind_var_act_species allows detecting that swap below.
                    !> Compute speciation algorithm arrays. The aqueous chemistry owns both the
                    !> solid (reactive zone) and gas chemistry, so it builds the modified equilibrium
                    !> constants K' = K/(prod a_cst) (constant-activity gas partial pressures fold in;
                    !> minerals have a=1). compute_speciation_alg_arrays builds ONLY the structural,
                    !> K-independent arrays on the (possibly shared) reactive zone and returns K after
                    !> any reaction swap; the transformed constants logK_tilde/logK_star are stored
                    !> PER-WATER (waters and waters_init are distinct aqueous chemistry objects that
                    !> may have different partial pressures, so each builds its own).
                    !> Structural, K-independent arrays on the (possibly shared) reactive zone.
                    !> Rebuilt only when the reactive zone changes (they depend only on its
                    !> stoichiometry, not on the individual water's temperature or partial pressures).
                    call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
                        flag_Se,swap)
                    
                    !> STEP 4: Handle species index swapping for numerical stability
                    !> The compute_speciation_alg_arrays routine may flag that two species should
                    !> swap positions to improve numerical conditioning of the linear system
                    !> If flag_Se=F, check if the shared reactive_zone already has a
                    !> swap applied (from a previous water that triggered the column swap).
                    if (.not. flag_Se) then
                        do j = 1, this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
                            if (this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%ind_var_act_species(j) /= j) then
                                flag_Se = .true.
                                swap(1) = j
                                swap(2) = this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%ind_var_act_species(j)
                                exit
                            end if
                        end do
                    end if
                    if (flag_Se.eqv..true.) then !< Swap operation required
                        !> Save original indices before swapping
                        aux_swap(1)=this%waters(tar_wat_ind)%ind_var_act_species(swap(1))
                        aux_swap(2)=this%waters(tar_wat_ind)%ind_var_act_species(swap(2))
                        !> Perform swap in current target water
                        this%waters(tar_wat_ind)%ind_var_act_species(swap(1))=aux_swap(2)
                        this%waters(tar_wat_ind)%ind_var_act_species(swap(2))=aux_swap(1)
                        this%waters(tar_wat_ind)%ind_prim_species=&
                            this%waters(tar_wat_ind)%ind_var_act_species(1:n_p_aq) !< update primary species indices accordingly
                        this%waters(tar_wat_ind)%ind_sec_species(1:n2v_aq)=&
                            this%waters(tar_wat_ind)%ind_var_act_species(n_p_aq+1:n_p_aq+n2v_aq) !< update secondary species indices accordingly
                        !> Perform same swap in initial target water to maintain consistency
                        this%waters_init(tar_wat_ind)%ind_var_act_species(swap(1))=aux_swap(2)
                        this%waters_init(tar_wat_ind)%ind_var_act_species(swap(2))=aux_swap(1)
                        this%waters_init(tar_wat_ind)%ind_prim_species=&
                            this%waters_init(tar_wat_ind)%ind_var_act_species(1:n_p_aq) !< update primary species indices accordingly
                        this%waters_init(tar_wat_ind)%ind_sec_species(1:n2v_aq)=&
                            this%waters_init(tar_wat_ind)%ind_var_act_species(n_p_aq+1:n_p_aq+n2v_aq) !< update secondary species indices accordingly
                    end if
                else !< SAME REACTIVE ZONE as previous target water
                    !> Optimization: reuse species indices from previous water (no need to recompute)
                    this%waters(tar_wat_ind)%ind_var_act_species=this%waters(tar_wat_ind-1)%ind_var_act_species
                    this%waters_init(tar_wat_ind)%ind_var_act_species=&
                        this%waters_init(tar_wat_ind-1)%ind_var_act_species
                    this%waters(tar_wat_ind)%ind_prim_species=&
                        this%waters(tar_wat_ind-1)%ind_prim_species
                    this%waters_init(tar_wat_ind)%ind_prim_species=&
                        this%waters_init(tar_wat_ind-1)%ind_prim_species
                    this%waters(tar_wat_ind)%ind_sec_species=&
                        this%waters(tar_wat_ind-1)%ind_sec_species
                    this%waters_init(tar_wat_ind)%ind_sec_species=&
                        this%waters_init(tar_wat_ind-1)%ind_sec_species
                end if
                !> Per-water temperature-dependent equilibrium constants (logK) and transformed
                !> constants (logK_tilde/logK_star), built for EVERY target water regardless of
                !> whether the reactive zone changed: logK depends on the water's own temperature
                !> and logK_tilde/logK_star on its own constant-activity gas partial pressures
                !> (folded in by get_eq_csts_aq_chem). The structural arrays they rely on live on
                !> the shared reactive zone and were (re)built above only when it changed.
                block
                    call this%waters(tar_wat_ind)%compute_logK_PHREEQC_aq_chem()
                    call this%waters_init(tar_wat_ind)%compute_logK_PHREEQC_aq_chem()
                    !> Transformed constants: each reads its own logK and folds in its own
                    !> constant-activity gas partial pressures (minerals have a=1). Select
                    !> logK_star (component elimination) or logK_tilde (all secondary species)
                    !> to match this reactive zone's active speciation algebra.
                    if (this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%flag_comp) then
                        call this%waters(tar_wat_ind)%compute_logK_star() !> this water's transformed constants (logK_star)
                        call this%waters_init(tar_wat_ind)%compute_logK_star() !> initial water's transformed constants (logK_star)
                    else
                        call this%waters(tar_wat_ind)%compute_logK_tilde_aq_chem() !> this water's transformed constants (logK_tilde)
                        call this%waters_init(tar_wat_ind)%compute_logK_tilde_aq_chem() !> initial water's transformed constants (logK_tilde)
                    end if
                end block
                !> STEP 4b: (Re)solve this domain cell's cation-exchange equilibrium with the
                !> now-known CEC. The water type's exchange was speciated before the surface
                !> adsorption zone was read (CEC=0), collapsing the free site and adsorbed
                !> cations to ~0; here each domain (non-external) target water re-solves the
                !> ideal Gaines-Thomas exchange from its aqueous cation concentrations and the
                !> merged reactive zone's CEC, giving a physical free site and adsorbed cations.
                !> The current and initial snapshots own separate target solids, so solve each.
                !> External (recharge/boundary) waters share these solids and are not re-solved.
                if (.not. flag_ext) then
                    call solve_cell_exchange(this%waters(tar_wat_ind))
                    call solve_cell_exchange(this%waters_init(tar_wat_ind))
                end if
        !> STEP 5: Compute stoichiometric transformation matrix U_SkT_prod for kinetic reactions
        !> TEMPORARY FIX: This section handles mineral zones with fewer kinetic minerals than
        !> the full chemical system, computing transformation matrices for subset of minerals
                if (this%waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin<&
                    this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin) then
                    !> Partial mineral zone: compute U_SkT_prod for subset of kinetic minerals present in this zone
                    call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod(&
                        this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                        this%waters(tar_wat_ind)%solid_chemistry%mineral_zone%ind_min_chem_syst(1:&
                        this%waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin))
                else
                    !> Full mineral zone: compute U_SkT_prod for all kinetic minerals in chemical system
                    call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod()
                end if
            end if
            !> STEP 6: Allocate and initialize reaction rates, amounts and their contributions to mass balances. This is done for both the current and initial target waters, ensuring that the reaction rates are properly set up for subsequent calculations. The reaction rates are computed based on the current state of the aqueous chemistry, and the initial rates are stored for reference. The pore volume for each water is also set, which is important for scaling the reaction rates to the volume of the water phase
            !> Set up reaction rate arrays and compute initial rates for non-external waters
            call this%waters(tar_wat_ind)%allocate_rxn_rates_and_amounts()      !< allocate reaction rates array for current water
            call this%waters_init(tar_wat_ind)%allocate_rxn_rates_and_amounts() !< allocate reaction rates array for initial water
            call this%waters(tar_wat_ind)%allocate_contrib_rxn_amnt()          !< allocate contribution reaction amounts for current water
            call this%waters_init(tar_wat_ind)%allocate_contrib_rxn_amnt()     !< allocate contribution reaction amounts for initial water
            call this%waters(tar_wat_ind)%allocate_conc_mix()          !< allocate contribution reaction rates for current water
            call this%waters_init(tar_wat_ind)%allocate_conc_mix()     !< allocate contribution reaction rates for initial water
            call this%waters(tar_wat_ind)%set_indices_rk()               !< set indices mapping for reaction rates (current)
            call this%waters_init(tar_wat_ind)%set_indices_rk()          !< set indices mapping for reaction rates (initial)
            call this%waters(tar_wat_ind)%set_volume()                   !< set pore volume for current target water
            call this%waters_init(tar_wat_ind)%set_volume()              !< set pore volume for initial target water
            
            !> Compute initial reaction rates. Only for the DOMAIN target waters: an external (boundary or
            !> recharge) water is a fixed inflow composition that does not react, so its kinetic rates
            !> stay at the zero set by allocate_rxn_rates_and_amounts (the mixers read rk_old of every
            !> stencil water; a nonzero rate on an external water would inject a spurious previous-step
            !> kinetic term, e.g. a huge calcite dissolution rate on an acid inflow water).
            ! if (.not. flag_ext) then
            !     allocate(rk(this%waters(tar_wat_ind)%indices_rk%num_cols)) !< allocate temporary reaction rates vector
            !     call this%waters(tar_wat_ind)%compute_rxn_rate_kin(rk)           !< compute current reaction rates
            !     call this%waters(tar_wat_ind)%set_rk_old(rk)      !< compute initial reaction rates
            !     call this%waters_init(tar_wat_ind)%compute_rxn_rate_kin(rk)      !< compute initial reaction rates
            !     call this%waters_init(tar_wat_ind)%set_rk_old(rk)          !< set old rates for initial water
            !     deallocate(rk) !< clean up temporary array
            ! end if
            call this%waters(tar_wat_ind)%compute_saturation_mins() !< initial saturation indices of all minerals (kinetic and equilibrium) for the t=0 output
            call this%waters_init(tar_wat_ind)%compute_saturation_mins() !< initial saturation indices of all minerals (kinetic and equilibrium) for the t=0 output
            !> Update auxiliary zone indices for next iteration (optimization for detecting zone changes)
            if (tar_sol_ind>0) then
                aux_iszn=this%target_solids(tar_sol_ind)%id !< store current solid zone ID for comparison with next water
            end if
            if (tar_gas_ind>0) then
                aux_igzn=this%target_gases(tar_gas_ind)%id !< store current gas zone ID for comparison with next water
            end if
            
            !> Clean up temporary arrays
            deallocate(swap,aux_swap) !< deallocate species swap index arrays
        end subroutine


        !> @brief (Re)solve one domain cell's ideal cation-exchange (Gaines-Thomas) equilibrium
        !> @details Equilibrates the exchanger of a single target water using the now-known CEC
        !>          of its (merged) reactive zone. The water-type exchange speciation runs before
        !>          the surface-adsorption zone is read (CEC=0), so the free site and adsorbed
        !>          cations collapse to ~0; this routine re-solves them from the aqueous cation
        !>          concentrations so they satisfy both the exchange mass-action law and the site
        !>          balance conc_free + sum_i z_i*conc_ads_i = CEC. Mirrors the ideal surface
        !>          speciation used for water types with only direct constraints
        !>          (read_wat_type_CHEPROO Branch 2). No-op for waters without an exchanger.
        !> @param[in,out] aq Target water whose exchanger is (re)equilibrated
        subroutine solve_cell_exchange(aq)
            class(aqueous_chemistry_c), intent(inout) :: aq !< target water to (re)equilibrate
            real(kind=8), allocatable :: act_ads_ig(:)      !< initial guess for adsorbed cation activities (site fractions)
            integer(kind=4) :: nec                          !< number of exchange cations
            integer(kind=4) :: fer,ler                      !< first/last exchange reaction index in this water's logK
            integer(kind=4) :: niter                        !< Newton iterations performed
            integer(kind=4) :: oc                           !< offset of the surface-complex block in the solid arrays
            logical :: cv_flag                              !< convergence flag
            if (.not. associated(aq%solid_chemistry)) return
            if (aq%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl<=0) return
            if (aq%solid_chemistry%CEC<=0d0) return
            nec=aq%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
            !> Surface complexes occupy the trailing num_surf_compl solid slots (after minerals).
            oc=aq%solid_chemistry%num_solids-aq%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl
            !> Free site (x-, valence 1): Gaines-Thomas log10 gamma = log10(z)-log10(CEC) = -log10(CEC),
            !> so its concentration a/10^(log gamma) = (1-sum a)*CEC shares the molar (equivalent) basis
            !> of the adsorbed cations, making conc_free + sum_i z_i*conc_ads_i = CEC hold exactly. The
            !> adsorbed-cation coefficients (trailing slots) were already set from the CEC when the
            !> surface zone was read.
            aq%solid_chemistry%log_act_coeffs(oc+1)=-log10(aq%solid_chemistry%CEC)
            !> The exchange reactions are the last equilibrium reactions before the gases.
            ler=aq%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions- &
                aq%solid_chemistry%reactive_zone%gas_phase%num_gases_eq
            fer=ler-nec+1
            !> Ideal model: aqueous activities equal concentrations. exch_cat_indices are positions
            !> in the aq_phase species list, so map them to state-variable positions with
            !> indices_aq_species before indexing this water's concentrations (get_c1's reduced
            !> primary ordering differs and is shorter once minerals eliminate aqueous species).
            call aq%compute_activities_aq()
            block
                integer(kind=4) :: dbg_naq
                dbg_naq = size(aq%solid_chemistry%reactive_zone%chem_syst%aq_phase%aq_species)
                write(0,'(A,7(1X,I0),1X,L1)') ' [SCE-DBG] nsc,nsol,nec,oc,naqsp,sizeIAS,sizeconc,csassoc=', &
                    aq%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl, &
                    aq%solid_chemistry%num_solids, nec, oc, dbg_naq, size(aq%indices_aq_species), &
                    size(aq%concentrations), associated(aq%solid_chemistry%reactive_zone%chem_syst)
                write(0,'(A,*(1X,I0))') ' [SCE-DBG] exch_cat_indices=', &
                    aq%solid_chemistry%reactive_zone%cat_exch_zone%exch_cat_indices
                write(0,'(A,2(1X,I0),1X,I0)') ' [SCE-DBG] fer,ler,sizelogK=', fer, ler, size(aq%logK)
            end block
            if (nec==2) then
                !> Two exchange cations: closed-form binary solution.
                call aq%solid_chemistry%compute_conc_surf_ideal_bin( &
                    aq%concentrations(aq%indices_aq_species( &
                    aq%solid_chemistry%reactive_zone%cat_exch_zone%exch_cat_indices)), &
                    aq%logK(fer:ler))
            else
                !> More than two exchange cations: Newton-Raphson in log10 site fractions. Uniform
                !> initial guess (each < 1, summing to < 1) keeps the first residual finite.
                allocate(act_ads_ig(nec))
                act_ads_ig=1d0/(nec+1d0)
                call aq%solid_chemistry%compute_conc_surf_ideal_Newton( &
                    aq%concentrations(aq%indices_aq_species( &
                    aq%solid_chemistry%reactive_zone%cat_exch_zone%exch_cat_indices)), &
                    act_ads_ig,niter,cv_flag,aq%logK(fer:ler))
                deallocate(act_ads_ig)
            end if
        end subroutine


        !> @brief Read computation options from file
        !> @details Reads computation configuration from a file named {root}_comp_opts.dat,
        !>          including Jacobian option, lumping flag, downstream reaction rate method,
        !>          and reaction rate averaging method.
        !> @param[in,out] this Chemistry object to configure
        !> @param[in] dir Directory containing the computation options file
        !> @param[in] root Root name for the computation options file
        !> @param[in] unit File unit number for reading
        subroutine read_comp_opts(this,dir,root,unit)
            implicit none
            class(chemistry_c) :: this                    !< Chemistry object to configure
            character(len=*), intent(in) :: dir           !< Directory of the file name
            character(len=*), intent(in) :: root          !< Root of the file name
            integer(kind=4), intent(in) :: unit           !< File unit number

            integer(kind=4) :: opt                        !< Computation option read from file
            logical :: flag                               !< Lumping flag read from file
            character(len=256) :: label                   !< Label string read from file
            real(kind=8) :: mu                            !< Estimation parameter (unused)
            !> Open file with computation options
            open(unit,file=dir//root//'_comp_opts.dat',status='old',action='read') !< Open computation options file for reading
            do                                           !< Loop over file lines until end marker
                read(unit,*) label                       !< Read label from current file line
                if (label=='end') then                   !< Check for end-of-file marker
                    exit                                 !< Exit reading loop
                else if (label=='COMPUTATION OPTIONS') then !< Check for computation options section header
                    read(unit,*) opt                      !< Read Jacobian computation option from file
                    call this%set_Jac_opt(opt)            !< Set Jacobian computation method
                    read(unit,*) flag                     !< Read lumping flag from file
                    call this%set_lump_flag(flag)         !< Set spatial lumping flag
                    read(unit,*) opt                      !< Read downstream reaction rate option from file
                    call this%set_r_down_opt(opt)         !< Set downstream reaction rate estimation method
                    read(unit,*) opt                      !< Read reaction rate averaging option from file
                    call this%set_rk_avg_opt(opt)         !< Set reaction rate averaging method
                    read(unit,*) flag                     !< Read mixing flag from file
                    call this%set_mix_by_zones_flag(flag) !< Set mixing flag (per reactive zone vs. direct)
                else                                     !< Unknown label encountered
                    continue                             !< Skip to next line
                end if                                   !< End label parsing
            end do                                       !< End file reading loop
            close(unit)                                  !< Close computation options file
        end subroutine
        
        !> @brief Set the number of waters
        !> @details Sets the total number of target waters in the chemistry system.
        !> Target waters are specific water compositions used in reactive transport simulations.
        !> This value must be consistent with the allocated waters array.
        !> @param this Chemistry object
        !> @param num_wat Number of waters (must be positive)
        subroutine set_num_waters(this,num_wat)
        implicit none                                     !< Require explicit variable declarations
        class(chemistry_c) :: this                        !< Chemistry object instance
        integer(kind=4), intent(in) :: num_wat            !< Number of waters to set
        if (num_wat<1) then                               !< Validate number of waters is positive
            error stop "Number of waters must be greater than 0" !< Terminate with error for invalid input
        end if                                            !< End validation block
        this%num_waters=num_wat                           !< Assign validated number of waters to object attribute
        end subroutine
        
        !> @brief Allocate memory for water types array
        !> @details Safely allocates memory for the water types array, ensuring proper deallocation
        !> of any existing array. Validates that the number of water types is positive before
        !> allocation. Sets the number of water types in the chemistry object.
        !> @param this Chemistry object
        !> @param num_wat_types Number of water types to allocate (must be > 0)
        subroutine allocate_wat_types(this,num_wat_types)
        implicit none                                     !< Require explicit variable declarations
        class(chemistry_c) :: this                        !< Chemistry object to modify
        integer(kind=4), intent(in) :: num_wat_types      !< Number of water types to allocate
        if (allocated(this%wat_types)) deallocate(this%wat_types) !< Deallocate existing water types array if previously allocated
        if (num_wat_types<1) then                         !< Validate number of water types is positive
            error stop "Number of water types must be greater than 0" !< Terminate with error for invalid input
        end if                                            !< End validation block
        call this%set_num_wat_types(num_wat_types)         !< Set the number of water types attribute
        allocate(this%wat_types(this%num_wat_types))       !< Allocate water types array with validated size
        end subroutine
        
        !> @brief Associate target waters with target solids
        !> @details Links each target water with its corresponding solid chemistry object
        !> based on the provided indices. If an index is 0 or negative, assigns the default
        !> solid chemistry object (first target solid). This establishes the solid-aqueous
        !> chemistry coupling for reactive transport calculations.
        !> @param this Chemistry object
        !> @param ind_tar_sol Array of target solid indices (dimension: number of target waters)
        subroutine set_waters_target_solids(this,ind_tar_sol)
        implicit none                                     !< Require explicit variable declarations
        class(chemistry_c) :: this                        !< Chemistry object instance
        integer(kind=4), intent(in) :: ind_tar_sol(:)     !< Indices of target solids (dim=num_target_waters)
        
        integer(kind=4) :: i                              !< Loop index for iterating over waters
        
        do i=1,this%num_waters                            !< Loop over all target waters to link solids
            if (ind_tar_sol(i)>0) then                    !< Check if this water has a valid target solid index
                call this%waters(i)%set_solid_chemistry(this%target_solids(ind_tar_sol(i))) !< Link water to specified target solid chemistry
            else                                          !< No valid solid index (0 or negative)
                call this%waters(i)%set_solid_chemistry(this%target_solids(1)) !< Link to default (first) target solid as fallback
            end if                                        !< End conditional solid assignment
        end do                                            !< End loop over target waters
        end subroutine
        
        !> @brief Associate target solids with material indices
        !> @details Links target solid chemistry objects with material definitions based on
        !> provided material indices. This appears to contain a logic error as it sets
        !> solid chemistry on waters instead of target_solids. Should be reviewed.
        !> @param this Chemistry object  
        !> @param ind_materials Array of material indices (dimension: number of target solids)
        !> @warning This subroutine may contain a bug - check target assignment logic
        subroutine set_target_solids_materials(this,ind_materials)
        implicit none                                     !< Require explicit variable declarations
        class(chemistry_c) :: this                        !< Chemistry object instance
        integer(kind=4), intent(in) :: ind_materials(:)   !< Indices of materials (dim=num_target_solids)

        integer(kind=4) :: i                              !< Loop index for iterating over target solids

        do i=1,this%num_target_solids
            if (ind_materials(i)>0) then                                                        !< Check if this target water has an associated material index
                call this%waters(i)%set_solid_chemistry(this%target_solids(ind_materials(i))) !< Link target water to corresponding target solid chemistry using material index
            else                                                                                !< No material index specified for this target water
                call this%waters(i)%set_solid_chemistry(this%target_solids(1))           !< Link to first (default) target solid chemistry object as fallback
            end if                                                                              !< End conditional solid chemistry assignment
        end do                                                                                  !< End loop over all target waters
        end subroutine
        
        !> @brief Get the number of aqueous components in the chemical system
        !> @details Retrieves the number of primary aqueous species (components) from the
        !> chemical system's speciation algorithm. This represents the minimum set of
        !> independent aqueous species needed to describe the system composition.
        !> @param this Chemistry object
        !> @return num_aq_comps Number of aqueous components in the chemical system
        function get_num_aq_comps_chem_syst(this) result(num_aq_comps)
        implicit none                                     !< Require explicit variable declarations
        class(chemistry_c) :: this                        !< Chemistry object instance
        integer(kind=4) :: num_aq_comps                   !< Number of aqueous components in chemical system
        num_aq_comps=this%chem_syst%speciation_alg%num_aq_prim_species !< Retrieve count from chemical system's speciation algorithm object
        end function
        
        !> @brief Get the number of aqueous components in a domain (target) water
        !> @details Retrieves the number of primary aqueous species (components) from the
        !> speciation algebra of the reactive zone of a domain water, which may differ from the
        !> chemical system's count when constant-activity minerals eliminate components.
        !> @param[in] this Chemistry object
        !> @param[in] ind_tw Optional index into `tar_wat_indices` selecting the target water (default: first)
        !> @return num_aq_comps Number of aqueous components in the selected target water
        function get_num_aq_comps_tar_wat(this,ind_tw) result(num_aq_comps)
        implicit none                                     !< Require explicit variable declarations
        class(chemistry_c), intent(in) :: this            !< Chemistry object instance
        integer(kind=4), intent(in), optional :: ind_tw   !< Optional index into tar_wat_indices
        integer(kind=4) :: num_aq_comps                   !< Number of aqueous components in the target water
        integer(kind=4) :: idx                            !< Resolved position in tar_wat_indices
        if (present(ind_tw)) then                         !< Validate caller-provided index
            if (ind_tw<1 .or. ind_tw>size(this%tar_wat_indices)) then
                error stop "get_num_aq_comps_tar_wat: ind_tw out of range"
            end if
            idx=ind_tw                                    !< Use caller-provided index
        else
            if (.not. allocated(this%tar_wat_indices) .or. size(this%tar_wat_indices)<1) then
                error stop "get_num_aq_comps_tar_wat: no target waters available"
            end if
            idx=1                                         !< Default to the first target water
        end if
        num_aq_comps=this%waters(this%tar_wat_indices(idx))%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species !< Retrieve count from the local reactive zone's speciation algorithm
        end function
        
        !> @brief Get number of aqueous variable activity species in a target water
        !> @details Retrieves the number of aqueous species with variable activity coefficients
        !> for a specific target water or the first target water if no index is provided.
        !> Assumes all target waters have the same number of variable activity species.
        !> @param this Chemistry object
        !> @param ind_tw Optional index of target water (must be within valid range)
        !> @return n_v_aq Number of aqueous variable activity species
        function get_num_aq_var_act_species(this,ind_tw) result(n_v_aq)
        implicit none                                     !< Require explicit variable declarations
        class(chemistry_c) :: this                        !< Chemistry object instance
        integer(kind=4), intent(in), optional :: ind_tw   !< Optional index of target water
        integer(kind=4) :: n_v_aq                         !< Number of aqueous variable activity species
        if (present(ind_tw)) then                                                 !< Check if optional target water index was provided
            if (ind_tw>0 .and. ind_tw<=this%num_waters) then               !< Validate index is within valid range [1, num_waters]
                n_v_aq=this%reactive_zones(ind_tw)%speciation_alg%num_aq_var_act_species !< Retrieve count from specified reactive zone's speciation algorithm
            else                                                                  !< Index is out of valid bounds
                error stop "Index of target water out of bounds"                  !< Terminate with error message for invalid index
            end if                                                                !< End index validation
        else                                                                      !< No index provided, use default (first target water)
            n_v_aq=this%waters(1)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species !< Retrieve from first target water, assuming uniform species count across all waters
        end if                                                                    !< End conditional retrieval
        end function
        
        !> @brief Get number of water types
        !> @details Returns the total number of water types defined in the chemistry system.
        !> Water types represent distinct aqueous compositions used as building blocks
        !> for reactive transport simulations.
        !> @param this Chemistry object
        !> @return num_wat_types Number of water types in the system
        function get_num_wat_types(this) result(num_wat_types)
        implicit none                                     !< Require explicit variable declarations
        class(chemistry_c) :: this                        !< Chemistry object instance
        integer(kind=4) :: num_wat_types                  !< Number of water types in the system
        num_wat_types=this%num_wat_types                  !< Return water types count directly from object attribute
        end function
        
        !> @brief Get concentrations of aqueous components for all water types
        !> @details Retrieves the aqueous component concentrations for all defined water types
        !> in the chemistry system. Returns a 2D array where each column represents a water type
        !> and each row represents an aqueous component concentration. Assumes all water types
        !> have the same number of aqueous components.
        !> @param this Chemistry object
        !> @return conc_comp_wat_types 2D array of component concentrations (components × water types)
        function get_conc_comp_wat_types(this) result(conc_comp_wat_types)
        implicit none                                     !< Require explicit variable declarations
        class(chemistry_c) :: this                        !< Chemistry object instance
        real(kind=8), allocatable :: conc_comp_wat_types(:,:) !< Concentrations of aqueous components of water types

        integer(kind=4) :: num_aq_comps                   !< Number of aqueous components in the chemical system
        integer(kind=4) :: i                              !< Loop index

        num_aq_comps=this%get_num_aq_comps_chem_syst() !< Call getter function to retrieve number of aqueous components
        allocate(conc_comp_wat_types(num_aq_comps,this%num_wat_types)) !< Allocate 2D array with dimensions [components × water types]
        do i=1,this%num_wat_types                                       !< Loop over all defined water types
            conc_comp_wat_types(:,i)=this%wat_types(i)%get_u_aq()       !< Extract aqueous component concentrations for water type i and store in column i
        end do                                                          !< End water types loop
        end function

        !> @brief Write aqueous species and components (with concentrations) for all water types
        !> @details For each water type this writes: (1) the aqueous variable activity species with
        !> their concentrations, and (2) the aqueous components expressed as linear combinations of
        !> those species (from the water type's component matrix) together with their component
        !> concentrations. Each water type is written in its own basis (component definitions and
        !> counts may differ between types).
        !> @param this Chemistry object
        !> @param path File path for output
        !> @param filename Output filename
        subroutine write_conc_comp_wat_types(this,path,filename)
        implicit none                                     !< Require explicit variable declarations
        class(chemistry_c) :: this                        !< Chemistry object instance
        character(len=*), intent(in) :: path              !< Path to the output file
        character(len=*), intent(in) :: filename          !< Output filename

        real(kind=8), allocatable :: comp_mat_ref(:,:)    !< aqueous component matrix of the RESIDENT water type (rows=components, cols=species)
        real(kind=8), allocatable :: u_ref(:)             !< a water type's component concentrations in the resident basis
        real(kind=8), allocatable :: c_i(:)               !< a water type's species concentrations reordered to the resident species order
        real(kind=8) :: coef                              !< component-matrix coefficient
        real(kind=8), parameter :: tol_coef=1d-12         !< tolerance to treat a coefficient as zero
        integer(kind=4) :: i,ci,j,jj,nvg,iref,ncomp_ref,nsp_ref !< loop and size indices
        character(len=512) :: expr                        !< right-hand side of a component definition
        character(len=512), allocatable :: comp_expr(:)   !< component definitions (built once from the resident matrix)
        character(len=64) :: term,sp_name,coef_str        !< auxiliary strings
        character(len=64), allocatable :: ref_sp_names(:) !< resident water's variable activity species names (component-matrix column order)
        logical :: first                                  !< flag for first term in a component definition

        open(unit=10, file=path//filename, status='unknown', action='write', form='formatted') !< Open output file for writing water type compositions
        write(10,"(2x,'Aqueous species and components (with concentrations) of the water types')")

        !> Reference basis = the RESIDENT/initial water type (the type of the first domain target water).
        !> ALL water types' component concentrations are expressed in this water's component matrix,
        !> so external/boundary components are comparable to the domain (not each water's own basis).
        iref=1
        if (allocated(this%wat_type_indices) .and. this%num_target_waters>0) iref=this%wat_type_indices(1)
        comp_mat_ref=this%wat_types(iref)%solid_chemistry%reactive_zone%speciation_alg%comp_mat_aq !< resident component matrix
        ncomp_ref=size(comp_mat_ref,1) !< number of components in the resident basis
        nsp_ref=size(comp_mat_ref,2)   !< number of species (columns) in the resident basis
        allocate(ref_sp_names(nsp_ref))
        do j=1,nsp_ref
            ref_sp_names(j)=this%wat_types(iref)%aq_phase%aq_species(this%wat_types(iref)%indices_aq_phase(&
                this%wat_types(iref)%ind_var_act_species(j)))%name !< resident species j (component-matrix column j)
        end do
        write(10,"(2x,'Component basis: water type ',A,' (',I0,' components)')") &
            trim(this%wat_types(iref)%name), ncomp_ref

        !> Build the component definitions (u_i = linear combination of species) once, from the resident matrix
        allocate(comp_expr(ncomp_ref))
        do ci=1,ncomp_ref
            expr=''
            first=.true.
            do j=1,nsp_ref
                coef=comp_mat_ref(ci,j)
                if (abs(coef)<tol_coef) cycle !< skip zero coefficients
                sp_name=ref_sp_names(j)
                if (first) then
                    if (coef>0d0) then
                        if (abs(coef-1d0)<tol_coef) then
                            term=trim(sp_name)
                        else
                            write(coef_str,'(F0.4)') coef
                            term=trim(adjustl(coef_str))//' '//trim(sp_name)
                        end if
                    else
                        if (abs(coef+1d0)<tol_coef) then
                            term='- '//trim(sp_name)
                        else
                            write(coef_str,'(F0.4)') abs(coef)
                            term='- '//trim(adjustl(coef_str))//' '//trim(sp_name)
                        end if
                    end if
                    first=.false.
                else
                    if (coef>0d0) then
                        if (abs(coef-1d0)<tol_coef) then
                            term=' + '//trim(sp_name)
                        else
                            write(coef_str,'(F0.4)') coef
                            term=' + '//trim(adjustl(coef_str))//' '//trim(sp_name)
                        end if
                    else
                        if (abs(coef+1d0)<tol_coef) then
                            term=' - '//trim(sp_name)
                        else
                            write(coef_str,'(F0.4)') abs(coef)
                            term=' - '//trim(adjustl(coef_str))//' '//trim(sp_name)
                        end if
                    end if
                end if
                expr=trim(expr)//trim(term) !< append term (leading sign space preserved)
            end do
            comp_expr(ci)=expr
        end do

        do i=1,this%num_wat_types !< Loop over all water types
            nvg=min(this%wat_types(i)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species, &
                    size(this%wat_types(i)%ind_var_act_species)) !< number of aqueous variable activity species (guarded)

            write(10,"(/,2x,'=== Water type: ',A,' ===')") trim(this%wat_types(i)%name) !< water type name header

            !> Aqueous variable activity species with their concentrations
            write(10,"(4x,'Aqueous variable activity species (name : concentration):')")
            do j=1,nvg
                sp_name=this%wat_types(i)%aq_phase%aq_species(this%wat_types(i)%indices_aq_phase(&
                    this%wat_types(i)%ind_var_act_species(j)))%name !< name of variable activity species j
                write(10,"(6x,I0,': ',A,' : ',ES15.5)") j, trim(sp_name), &
                    this%wat_types(i)%concentrations(this%wat_types(i)%ind_var_act_species(j))
            end do

            !> Aqueous components in the RESIDENT water's basis: u_i = comp_mat_resident * c
            !> (c reordered by species name to the resident species order; absent species contribute 0)
            write(10,"(4x,'Aqueous components (u_i = linear combination of species : concentration):')")
            allocate(c_i(nsp_ref))
            do j=1,nsp_ref
                c_i(j)=0d0
                do jj=1,nvg
                    if (trim(this%wat_types(i)%aq_phase%aq_species(this%wat_types(i)%indices_aq_phase(&
                        this%wat_types(i)%ind_var_act_species(jj)))%name)==trim(ref_sp_names(j))) then
                        c_i(j)=this%wat_types(i)%concentrations(this%wat_types(i)%ind_var_act_species(jj))
                        exit
                    end if
                end do
            end do
            u_ref=matmul(comp_mat_ref,c_i) !< water type i's composition in the resident component basis
            do ci=1,ncomp_ref
                write(10,"(6x,'u',I0,' = ',A,' : ',ES15.5)") ci, trim(comp_expr(ci)), u_ref(ci)
            end do
            deallocate(c_i,u_ref)

            !> Gas partial pressures as computed by the code during speciation
            !> (gas_chemistry%activities holds P_i filled by compute_partial_pressures)
            if (associated(this%wat_types(i)%gas_chemistry)) then
                write(10,"(4x,'Gas partial pressures (name : P [bar], computed by the code):')")
                do j=1,this%wat_types(i)%gas_chemistry%reactive_zone%gas_phase%num_species
                    write(10,"(6x,I0,': ',A,' : ',ES15.5)") j, &
                        trim(this%wat_types(i)%gas_chemistry%reactive_zone%gas_phase%gases(j)%name), &
                        this%wat_types(i)%gas_chemistry%activities(j)
                end do
            end if
        end do !< End loop over water types
        deallocate(comp_mat_ref,ref_sp_names,comp_expr) !< free reference-basis arrays
        close(10) !< Close output file after writing all water type data
        end subroutine
        
        !> @brief Associate target solids with spatial mesh targets
        !> @details Links each target solid chemistry object with its corresponding spatial target
        !> from the mesh discretization. This establishes the spatial-chemical coupling needed
        !> for reactive transport simulations by connecting chemical properties to mesh elements.
        !> @param this Chemistry object
        !> @param mesh Spatial discretization object with defined targets
        subroutine set_target_solids_mesh(this,mesh)
        class(chemistry_c) :: this                        !< Chemistry object instance
        class(spatial_discr_c), intent(in) :: mesh        !< Mesh object (assumes targets are defined)
        integer(kind=4) :: i                              !< Loop index
        do i=1,this%num_target_solids !< Loop over all target solid chemistry objects to establish spatial coupling
            call this%target_solids(i)%set_target(mesh%targets(i)) !< Associate target solid i with corresponding spatial mesh target i, linking chemical properties to mesh element
            call this%target_solids_init(i)%set_target(mesh%targets(i)) !< Associate target solid i with corresponding spatial mesh target i, linking chemical properties to mesh element
        end do !< End loop over target solids
        end subroutine
        
        !> @brief Allocate memory for gas zones array
        !> @details Safely allocates memory for the gas zones array, deallocating any existing
        !> array first. Optionally sets the number of gas zones if provided as parameter.
        !> Gas zones represent distinct gas phase compositions in the system.
        !> @param this Chemistry object
        !> @param num_gas_zones Optional number of gas zones to allocate
        subroutine allocate_gas_zones(this,num_gas_zones)
        implicit none
        class(chemistry_c) :: this !< Chemistry object to modify
        integer(kind=4), intent(in), optional :: num_gas_zones !< Optional parameter to set number of gas zones
        if (present(num_gas_zones)) then !< Check if optional parameter was provided by caller
            this%num_gas_zones=num_gas_zones !< Update number of gas zones attribute from optional parameter
        end if !< End optional parameter handling
        if (allocated(this%gas_zones)) deallocate(this%gas_zones) !< Deallocate existing gas zones array if previously allocated to prevent memory leak
        allocate(this%gas_zones(this%num_gas_zones))              !< Allocate gas zones array with current count
        end subroutine

        !> @brief Allocate memory for initial target waters array
        !> @details Safely allocates memory for the initial target waters array, deallocating
        !> any existing array first. Optionally sets the number of initial target waters if
        !> provided as parameter. Initial target waters define starting compositions for
        !> reactive transport simulations.
        !> @param this Chemistry object
        !> @param num_waters_init Optional number of initial target waters to allocate
        subroutine allocate_waters_init(this,num_waters_init)
        implicit none
        class(chemistry_c) :: this !< Chemistry object to modify
        integer(kind=4), intent(in), optional :: num_waters_init !< Optional parameter to set number of initial target waters
        if (present(num_waters_init)) then !< Check if optional parameter was provided by caller
            this%num_waters_init=num_waters_init !< Update number of initial target waters attribute from optional parameter
        end if !< End optional parameter handling
        if (allocated(this%waters_init)) deallocate(this%waters_init) !< Deallocate existing initial target waters array if previously allocated to prevent memory leak
        allocate(this%waters_init(this%num_waters_init)) !< Allocate initial target waters array with current number of initial target waters
        end subroutine
        
        !> @brief Compute u_mix for initial target waters from mixing
        !> @details Calculates the u_mix (transformed component concentrations) for initial
        !> target waters based on mixing of multiple waters with specified mixing ratios.
        !> The algorithm first computes c_mix from the mixing waters, then transforms
        !> to u_mix coordinates. This is essential for initializing reactive transport
        !> simulations with mixed water compositions.
        !> @param this Chemistry object
        !> @param mix_conc_indices Matrix containing indices of target waters that mix
        !> @param mixing_ratios Matrix of mixing ratios for concentrations
        !> @return u_mix_init 2D array of u_mix values (components × domain target waters)
        function compute_u_mix_init(this,mix_conc_indices,mixing_ratios) result(u_mix_init)
        class(chemistry_c) :: this
        class(int_array_c), intent(in) :: mix_conc_indices !< matrix that contains indices of target waters that mix with each target water
        class(real_array_c), intent(in) :: mixing_ratios !< mixing ratios matrix for concentrations
        real(kind=8), allocatable :: u_mix_init(:,:) !< u_mix of initial target waters
        
        integer(kind=4) :: i !< loop index
        integer(kind=4) :: j_mix !< loop index for mixing waters
        integer(kind=4) :: num_mix_loc !< number of mixing waters
        real(kind=8), allocatable :: c_mix(:) !< c_mix of mixing waters
        real(kind=8), allocatable :: conc_old_mix(:,:) !< conc_old from each mixing water
        integer(kind=4), allocatable :: ind_aq_sp_mix(:,:) !< indices_aq_species from each mixing water
        ! Body of compute_u_mix_init
        allocate(u_mix_init(this%get_num_aq_comps_chem_syst(),this%num_target_waters)) !< Allocate result array: rows = aqueous components, columns = domain target waters
        do i=1,this%num_target_waters !< Loop over all domain target waters to compute transformed concentrations from mixing
            allocate(c_mix(this%waters_init(this%tar_wat_indices(i)&
                )%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species)) !< Allocate c_mix with size = number of variable activity aqueous species in this water's reactive zone
            num_mix_loc=mix_conc_indices%cols(i)%dim-3
            allocate(conc_old_mix(size(this%waters_init(mix_conc_indices%cols(i)%col_1(1))%conc_old), num_mix_loc))
            allocate(ind_aq_sp_mix(size(this%waters_init(mix_conc_indices%cols(i)%col_1(1))%indices_aq_species), num_mix_loc))
            do j_mix=1,num_mix_loc
                conc_old_mix(:,j_mix)=this%waters_init(mix_conc_indices%cols(i)%col_1(j_mix+1))%conc_old
                ind_aq_sp_mix(:,j_mix)=this%waters_init(mix_conc_indices%cols(i)%col_1(j_mix+1))%indices_aq_species
            end do
            call compute_c_mix(this%waters_init(mix_conc_indices%cols(i)%col_1(1)),&
                conc_old_mix,ind_aq_sp_mix,&
                mixing_ratios%cols(i)%col_1,c_mix) !< Compute c_mix (transformed concentrations) by mixing waters using provided indices and ratios for domain water i
            deallocate(conc_old_mix,ind_aq_sp_mix)
            u_mix_init(:,i)=this%waters_init(this%tar_wat_indices(i))%compute_u_mix(c_mix) !< Transform c_mix to u_mix coordinate system for domain water i using chemical system's transformation
            deallocate(c_mix) !< Free temporary c_mix array before next iteration
        end do !< End loop over domain target waters
        end function compute_u_mix_init
        
        !> @brief Write u_mix initial values to file
        !> @details Outputs the u_mix (transformed component concentrations) for initial target
        !> waters to a formatted file. Each row represents an aqueous component and each column
        !> represents a domain target water. Values are written in scientific notation for
        !> numerical precision and analysis purposes.
        !> @param this Chemistry object
        !> @param path Output file path
        !> @param filename Output filename
        !> @param u_mix_init 2D array of u_mix values to write
        subroutine write_u_mix_init(this,path,filename,u_mix_init)
        class(chemistry_c) :: this !< chemistry object
        character(len=*), intent(in) :: path !< path to write the file
        character(len=*), intent(in) :: filename !< name of the file
        real(kind=8), intent(in) :: u_mix_init(:,:) !< u_mix of initial target waters
        
        integer(kind=4) :: i !< loop index
        
        open(unit=20,file=path//filename,status='replace',action='write') !< Open file for writing, replacing if it exists
        do i=1,size(u_mix_init,1) !< Loop over all aqueous components (rows of u_mix_init matrix)
            write(20,"(*(ES15.5))") u_mix_init(i,:) !< Write row i (one component across all domain target waters) in scientific notation
        end do !< End loop over components
        close(20) !< Close output file after writing complete u_mix matrix
        end subroutine
        
        !> @brief Write initial aqueous components and species description to file
        !> @details Generates a human-readable summary file (e.g. u_init.dat) describing the
        !> aqueous component definitions in terms of variable activity species, the list of
        !> aqueous variable activity species, and the initial aqueous component and species
        !> concentrations both in the domain (target waters) and in the external (boundary/recharge)
        !> waters. The format reproduces the layout expected by the REMIX interface.
        !> @param this Chemistry object
        !> @param path Output file path
        !> @param filename Output filename
        subroutine write_u_init(this,path,filename)
        class(chemistry_c) :: this !< chemistry object
        character(len=*), intent(in) :: path !< path to write the file
        character(len=*), intent(in) :: filename !< name of the file

        integer(kind=4) :: i,ci,j,jj,wi,iref,ncomp,nsp_ref,n_tar,n_ext !< loop and size indices
        real(kind=8), parameter :: tol_coef=1d-12 !< tolerance to treat a component-matrix coefficient as zero
        real(kind=8), allocatable :: comp_mat_ref(:,:) !< aqueous component matrix of the resident/initial water type
        real(kind=8), allocatable :: c_i(:) !< a water's species concentrations reordered to the resident species order
        real(kind=8), allocatable :: ug_dom(:,:),cg_dom(:,:),ug_ext(:,:),cg_ext(:,:) !< concentration tables (resident basis)
        character(len=64), allocatable :: ref_sp_names(:) !< resident water type variable activity species names (comp-matrix column order)
        character(len=512) :: expr !< right-hand side of a component definition
        character(len=64) :: term,coef_str !< auxiliary strings
        real(kind=8) :: coef !< component-matrix coefficient
        logical :: first !< flag for first term in a component definition

        !> Reference/resident basis = the INITIAL (resident) water TYPE of the first domain target
        !> water. All component definitions, the species list and every concentration column (domain
        !> and external) are expressed in this single basis (species reordered to the resident order,
        !> components u = comp_mat_ref * c); waters are written in order (no per-reactive-zone grouping).
        iref=1
        if (allocated(this%wat_type_indices) .and. this%num_target_waters>0) iref=this%wat_type_indices(1)
        comp_mat_ref=this%wat_types(iref)%solid_chemistry%reactive_zone%speciation_alg%comp_mat_aq
        ncomp=size(comp_mat_ref,1)
        nsp_ref=size(comp_mat_ref,2)
        allocate(ref_sp_names(nsp_ref))
        do j=1,nsp_ref
            ref_sp_names(j)=this%wat_types(iref)%aq_phase%aq_species(this%wat_types(iref)%indices_aq_phase(&
                this%wat_types(iref)%ind_var_act_species(j)))%name
        end do
        n_tar=size(this%tar_wat_indices)
        n_ext=size(this%ext_waters_indices)

        open(unit=20,file=path//filename,status='replace',action='write') !< open output file, replacing if it exists

        !> Number of aqueous components (resident water type basis)
        write(20,"(2x,'Number of aqueous components: ',I0)") ncomp

        !> Variable activity species (resident basis order): primaries (rows 1..ncomp) then secondaries
        write(20,*) !< blank line between sections
        write(20,"(2x,'Aqueous variable activity species (rows 1 to ',I0,': primary; rows ',I0,' to ',I0,': secondary):')") &
            ncomp, ncomp+1, nsp_ref
        do j=1,nsp_ref
            write(20,"(4x,I0,': ',A)") j, trim(ref_sp_names(j))
        end do

        !> Component definitions (u_i = linear combination of variable activity species), resident basis
        write(20,*) !< blank line between sections
        write(20,"(2x,'Aqueous components:')")
        do ci=1,ncomp
            expr=''
            first=.true.
            do j=1,nsp_ref
                coef=comp_mat_ref(ci,j)
                if (abs(coef)<tol_coef) cycle
                if (first) then
                    if (coef>0d0) then
                        if (abs(coef-1d0)<tol_coef) then
                            term=trim(ref_sp_names(j))
                        else
                            write(coef_str,'(F0.4)') coef
                            term=trim(adjustl(coef_str))//' '//trim(ref_sp_names(j))
                        end if
                    else
                        if (abs(coef+1d0)<tol_coef) then
                            term='- '//trim(ref_sp_names(j))
                        else
                            write(coef_str,'(F0.4)') abs(coef)
                            term='- '//trim(adjustl(coef_str))//' '//trim(ref_sp_names(j))
                        end if
                    end if
                    first=.false.
                else
                    if (coef>0d0) then
                        if (abs(coef-1d0)<tol_coef) then
                            term=' + '//trim(ref_sp_names(j))
                        else
                            write(coef_str,'(F0.4)') coef
                            term=' + '//trim(adjustl(coef_str))//' '//trim(ref_sp_names(j))
                        end if
                    else
                        if (abs(coef+1d0)<tol_coef) then
                            term=' - '//trim(ref_sp_names(j))
                        else
                            write(coef_str,'(F0.4)') abs(coef)
                            term=' - '//trim(adjustl(coef_str))//' '//trim(ref_sp_names(j))
                        end if
                    end if
                end if
                expr=trim(expr)//trim(term)
            end do
            write(20,"(4x,'u',I0,' = ',A)") ci, trim(expr)
        end do

        !> Domain (target waters): concentrations expressed in the resident basis (species reordered
        !> to the resident order; components u = comp_mat_ref * c). Columns = target waters, in order.
        allocate(ug_dom(ncomp,n_tar),cg_dom(nsp_ref,n_tar),c_i(nsp_ref))
        do i=1,n_tar
            wi=this%tar_wat_indices(i)
            c_i=0d0
            do j=1,nsp_ref
                do jj=1,min(this%waters(wi)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species, &
                            size(this%waters(wi)%ind_var_act_species))
                    if (trim(this%waters(wi)%aq_phase%aq_species(this%waters(wi)%indices_aq_phase(&
                        this%waters(wi)%ind_var_act_species(jj)))%name)==trim(ref_sp_names(j))) then
                        c_i(j)=this%waters(wi)%concentrations(this%waters(wi)%ind_var_act_species(jj))
                        exit
                    end if
                end do
            end do
            cg_dom(:,i)=c_i
            ug_dom(:,i)=matmul(comp_mat_ref,c_i)
        end do
        write(20,*) !< blank line between sections
        write(20,"(2x,'Aqueous variable activity species concentrations in the domain (rows = species, cols = target waters):')")
        do j=1,nsp_ref
            write(20,"(2x,*(ES15.5))") cg_dom(j,:)
        end do
        write(20,*) !< blank line between sections
        write(20,"(2x,'Aqueous component concentrations in the domain (rows = components, cols = target waters):')")
        do ci=1,ncomp
            write(20,"(2x,*(ES15.5))") ug_dom(ci,:)
        end do

        !> External (boundary/recharge) waters: same resident-basis reordering. Columns = external waters, in order.
        allocate(ug_ext(ncomp,n_ext),cg_ext(nsp_ref,n_ext))
        do i=1,n_ext
            wi=this%ext_waters_indices(i)
            c_i=0d0
            do j=1,nsp_ref
                do jj=1,min(this%waters(wi)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species, &
                            size(this%waters(wi)%ind_var_act_species))
                    if (trim(this%waters(wi)%aq_phase%aq_species(this%waters(wi)%indices_aq_phase(&
                        this%waters(wi)%ind_var_act_species(jj)))%name)==trim(ref_sp_names(j))) then
                        c_i(j)=this%waters(wi)%concentrations(this%waters(wi)%ind_var_act_species(jj))
                        exit
                    end if
                end do
            end do
            cg_ext(:,i)=c_i
            ug_ext(:,i)=matmul(comp_mat_ref,c_i)
        end do
        write(20,*) !< blank line between sections
        write(20,"(2x,'Aqueous variable activity species concentrations of external waters (rows = species, cols = external waters):')")
        write(20,"(2x,*(A24))") (trim(this%waters(this%ext_waters_indices(i))%name), i=1,n_ext)
        do j=1,nsp_ref
            write(20,"(2x,*(ES24.10))") cg_ext(j,:)
        end do
        write(20,*) !< blank line between sections
        write(20,"(2x,'Aqueous component concentrations of external waters (rows = components, cols = external waters):')")
        write(20,"(2x,*(A24))") (trim(this%waters(this%ext_waters_indices(i))%name), i=1,n_ext)
        do ci=1,ncomp
            write(20,"(2x,*(ES24.10))") ug_ext(ci,:)
        end do

        close(20) !< close output file

        end subroutine
        
        !> @brief Set the number of materials in the system
        !> @details Sets the total number of material objects in the chemistry system.
        !> Materials represent distinct solid compositions with associated chemical properties.
        !> Validates that the number is non-negative before assignment.
        !> @param this Chemistry object
        !> @param num_materials Number of materials (must be ≥ 0)
        subroutine set_num_materials(this,num_materials)
        class(chemistry_c) :: this                        !< Chemistry object instance
        integer(kind=4), intent(in) :: num_materials      !< Number of materials
        if (num_materials<0) then !< Validate that number of materials is non-negative (physical constraint)
            error stop "Number of materials must be non-negative" !< Abort execution if validation fails with descriptive error message
        end if !< End validation block
        this%num_materials=num_materials !< Assign validated number of materials to chemistry object attribute
        end subroutine
        
        !> @brief Allocate memory for materials array
        !> @details Allocates memory for the materials array based on the previously set
        !> number of materials. Each material object will contain solid chemistry properties
        !> and associations with target solids in the reactive transport system.
        !> @param this Chemistry object
        subroutine allocate_materials(this)
        class(chemistry_c) :: this                        !< Chemistry object instance
        allocate(this%materials(this%num_materials))       !< Allocate materials array with size equal to the number of materials set previously
        end subroutine
        
        !> @brief Get indices of target solids associated with a material
        !> @details Finds all target solid indices that are associated with a specific material ID.
        !> This function scans through all target solids and returns an array containing the
        !> indices of those that match the provided material index. Essential for linking
        !> material properties to spatial locations in reactive transport simulations.
        !> @param this Chemistry object
        !> @param imat Index of material to search for
        !> @return tar_sol_ind Array of target solid indices associated with the material
        function get_tar_sol_ind(this,imat) result(tar_sol_ind) !< function to get indices of target solids associated to a material
        class(chemistry_c) :: this !< chemistry object
        integer(kind=4), intent(in) :: imat !< index of material
        integer(kind=4), allocatable :: tar_sol_ind(:) !< indices of target solids associated to material
        integer(kind=4) :: i !< loop index for scanning target solids
        integer(kind=4) :: j !< loop index for filling result array
        integer(kind=4) :: dim !< dimension of tar_sol_ind (count of matching target solids)
        dim=0 !< Initialize dimension counter to zero before first scan
        do i=1,this%num_target_solids !< First pass: count how many target solids have material ID matching imat
            if (this%target_solids(i)%id==imat) then !< Check if current target solid belongs to requested material
                dim=dim+1 !< Increment counter for each matching target solid
            end if !< End material ID match check
        end do !< End first pass counting loop
        allocate(tar_sol_ind(dim)) !< Allocate result array with exact size needed (number of matching target solids)
        !i=1
        j=1 !< Initialize index for filling result array
        do i=1,this%num_target_solids !< Second pass: populate result array with indices of matching target solids
            if (this%target_solids(i)%id==imat) then !< Check if current target solid belongs to requested material
                tar_sol_ind(j)=i !< Store index i of matching target solid in result array at position j
                if (j<dim) then !< Check if more matching target solids remain to be found
                    j=j+1 !< Advance result array index to next position
                else !< All matching target solids have been found
                    exit !< Exit loop early once result array is completely filled
                end if !< End check for remaining matches
            end if !< End material ID match check
        end do !< End second pass population loop
        end function
        
        !> @brief Get indices of target gases associated with a gas zone
        !> @details Finds all target gas indices that are associated with a specific gas zone ID.
        !> This function scans through all target gases and returns an array containing the
        !> indices of those that match the provided gas zone index. Used for linking gas
        !> zone properties to spatial locations in multiphase reactive transport.
        !> @param this Chemistry object
        !> @param igzn Index of gas zone to search for
        !> @return tar_gas_ind Array of target gas indices associated with the gas zone
        function get_tar_gas_ind(this,igzn) result(tar_gas_ind) !< function to get indices of target solids associated to a gas zone
        class(chemistry_c) :: this !< chemistry object
        integer(kind=4), intent(in) :: igzn !< index of gas zone to search for
        integer(kind=4), allocatable :: tar_gas_ind(:) !< indices of target gases associated to gas zone (result array)
        integer(kind=4) :: i,j !< loop indices (i for scanning, j for filling result)
        integer(kind=4) :: dim !< dimension of tar_gas_ind (count of matching target gases)
        dim=0 !< Initialize dimension counter to zero before first scan
        do i=1,this%num_target_gases !< First pass: count how many target gases have gas zone ID matching igzn
            if (this%target_gases(i)%id==igzn) then !< Check if current target gas belongs to requested gas zone
                dim=dim+1 !< Increment counter for each matching target gas
            end if !< End gas zone ID match check
        end do !< End first pass counting loop
        allocate(tar_gas_ind(dim)) !< Allocate result array with exact size needed (number of matching target gases)
        !i=1
        j=1 !< Initialize index for filling result array
        do i=1,this%num_target_gases !< Second pass: populate result array with indices of matching target gases
            if (this%target_gases(i)%id==igzn) then !< Check if current target gas belongs to requested gas zone
                tar_gas_ind(j)=i !< Store index i of matching target gas in result array at position j
                if (j<dim) then !< Check if more matching target gases remain to be found
                    j=j+1 !< Advance result array index to next position
                else !< All matching target gases have been found
                    exit !< Exit loop early once result array is completely filled
                end if !< End check for remaining matches
            end if !< End gas zone ID match check
        end do !< End second pass population loop
        end function

        !> @brief Allocate memory for solid types array
        !> @details Safely allocates memory for the solid types array, assuming a bijective
        !> relationship between water types and solid types. Deallocates any existing array
        !> before allocation. Each solid type corresponds to a specific water type composition.
        !> @param this Chemistry object
        subroutine allocate_sol_types(this)
        implicit none
        class(chemistry_c) :: this !< Chemistry object to modify
        if (allocated(this%sol_types)) deallocate(this%sol_types) !< Deallocate existing solid types array if previously allocated to prevent memory leak
        allocate(this%sol_types(this%num_wat_types)) !< Allocate solid types array with size equal to number of water types (assuming 1-to-1 correspondence)
        end subroutine allocate_sol_types

        !> @brief Allocate memory for reactive zones associated with water types
        !> @details Safely allocates memory for the reactive zones array associated with water types,
        !> assuming a bijective relationship between water types and reactive zones. Each water type
        !> corresponds to a specific reactive zone with defined chemical reactions and speciation.
        !> @param this Chemistry object
        subroutine allocate_react_zones_wat_types(this)
        implicit none
        class(chemistry_c) :: this !< Chemistry object to modify
        if (allocated(this%react_zones_wat_types)) deallocate(this%react_zones_wat_types) !< Deallocate existing reactive zones array if previously allocated to prevent memory leak
        allocate(this%react_zones_wat_types(this%num_wat_types)) !< Allocate reactive zones array with size equal to number of water types (assuming 1-to-1 correspondence)
        end subroutine allocate_react_zones_wat_types

        !> @brief Allocate memory for mineral zones associated with water types
        !> @details Safely allocates the mineral zones array associated with water types,
        !> analogous to allocate_react_zones_wat_types. Deallocates any existing array first.
        !> @param this Chemistry object
        subroutine allocate_min_zones_wat_types(this)
        implicit none
        class(chemistry_c) :: this !< Chemistry object to modify
        if (allocated(this%min_zones_wat_types)) deallocate(this%min_zones_wat_types) !< Deallocate existing mineral zones array if previously allocated to prevent memory leak
        allocate(this%min_zones_wat_types(this%num_wat_types)) !< Allocate mineral zones array with size equal to number of water types
        end subroutine allocate_min_zones_wat_types

        !> @brief Allocate memory for gas zones associated with water types
        !> @details Safely allocates memory for the gas zones array associated with water types,
        !> assuming a bijective relationship between water types and gas zones. Each water type
        !> corresponds to a specific gas zone with defined gas phase composition and equilibrium.
        !> @param this Chemistry object
        subroutine allocate_gas_zones_wat_types(this)
        implicit none
        class(chemistry_c) :: this !< Chemistry object to modify
        if (allocated(this%gas_zones_wat_types)) deallocate(this%gas_zones_wat_types) !< Deallocate existing gas zones array if previously allocated to prevent memory leak
        allocate(this%gas_zones_wat_types(this%num_wat_types)) !< Allocate gas zones array with size equal to number of water types (assuming 1-to-1 correspondence)
        end subroutine allocate_gas_zones_wat_types



        !> @brief Set indices of external waters (boundary and recharge)
        !> @details Populates the ext_waters_indices array with indices of boundary waters
        !> followed by recharge waters. This is used to identify all external water sources
        !> in the simulation domain for boundary condition and recharge handling.
        !> @param this Chemistry object
        subroutine set_ext_waters_indices(this)
            implicit none
            class(chemistry_c) :: this !< Chemistry object to modify
            integer(kind=4) :: i !< Loop index for populating external waters indices
            !> First: boundary waters
            do i=1,this%num_bd_waters !< Loop over boundary waters to copy their indices to external waters array first
                this%ext_waters_indices(i)=this%bd_waters_indices(i) !< Copy boundary water index i to first section of external waters indices array
            end do !< End boundary waters loop
            !> Second: recharge waters
            do i=1,this%num_rech_waters !< Loop over recharge waters to append their indices after boundary waters
                this%ext_waters_indices(this%num_bd_waters+i)=this%rech_waters_indices(i) !< Copy recharge water index i to second section of external waters indices array, offset by number of boundary waters
            end do !< End recharge waters loop
        end subroutine


        !> @brief Allocate memory for initial cation exchange zones
        !> @details Safely allocates memory for the initial cation exchange zones array, deallocating
        !> any existing array first. Optionally sets the number of initial cation exchange zones if
        !> provided as parameter. Used for initializing cation exchange properties in the system.
        !> @param this Chemistry object
        !> @param num_init_cat_exch_zones Optional number of initial cation exchange zones to allocate
        subroutine allocate_init_cat_exch_zones(this,num_init_cat_exch_zones)
            implicit none
            class(chemistry_c) :: this !< Chemistry object to modify
            integer(kind=4), intent(in), optional :: num_init_cat_exch_zones !< Optional parameter to set number of initial cation exchange zones
            if (present(num_init_cat_exch_zones)) then !< Check if optional parameter was provided by caller
                this%num_init_cat_exch_zones=num_init_cat_exch_zones !< Update number of initial cation exchange zones attribute from optional parameter
            end if !< End optional parameter handling
            if (allocated(this%init_cat_exch_zones)) deallocate(this%init_cat_exch_zones) !< Deallocate existing initial cation exchange zones array if previously allocated to prevent memory leak
            allocate(this%init_cat_exch_zones(this%num_init_cat_exch_zones)) !< Allocate initial cation exchange zones array with current number of zones
        end subroutine

        !> @brief Water mixing algorithm with explicit Euler, equilibrium/kinetic reactions, ideal activity, lumped formulation
        !> @details Performs mixing iteration using the Weighted Mixing Algorithm (WMA)
        !> with the following characteristics:
        !> - **Temporal discretization**: Explicit Euler (EE) for transport step
        !> - **Chemical reactions**: Both equilibrium (instantaneous) and kinetic (rate-limited)
        !> - **Activity model**: Ideal solution (activity coefficients = 1)
        !> - **Formulation**: Lumped (reduced computational complexity)
        !>
        !> Algorithm steps for each domain target water:
        !> 1. Retrieve old primary species concentrations (c1_old) from previous time step
        !> 2. Compute transformed concentrations (c_mix) by mixing waters using provided indices and lambda ratios
        !> 3. Solve coupled transport-reaction system using explicit Euler with equilibrium-kinetic split
        !> 4. Return updated aqueous component (u_new) and species (c_new) concentrations
        !>
        !> This is a high-level coordinator that loops over all domain target waters and delegates
        !> the actual water mixing computation to individual target water objects.
        !>
        !> @param[in,out] this Chemistry object coordinating the mixing algorithm
        !> @param[in] mix_wat_indices Matrix of integer arrays containing indices of waters that mix for each target water
        !> @param[in] lambdas Matrix of real arrays containing mixing ratios (weights) for weighted mixing algorithm
        !> @param[in] Delta_t Time step size [time units] for explicit Euler temporal discretization
        !> @param[out] u_new New aqueous component concentrations [mol/L] - 2D array (components × domain waters)
        !> @param[out] c_new New species concentrations [mol/L] - 2D array (species × domain waters)
        subroutine WMA_iter_EE_eq_kin_ideal_lump(this,mix_wat_indices,lambdas,Delta_t,u_new,c_new)
            class(chemistry_c) :: this !< Chemistry object coordinating water mixing algorithm
            type(int_array_c), intent(in) :: mix_wat_indices !< Matrix of indices indicating which waters mix for each target water
            type(real_array_c), intent(in) :: lambdas !< Matrix of mixing ratios (lambdas) for weighted mixing algorithm (WMA)
            real(kind=8), intent(in) :: Delta_t !< Time step size for explicit Euler temporal discretization
            real(kind=8), intent(out) :: u_new(:,:) !< Output: new aqueous component concentrations (components × domain waters)
            real(kind=8), intent(out) :: c_new(:,:) !< Output: new species concentrations (species × domain waters)

            integer(kind=4) :: i !< Loop index for iterating over domain target waters
            integer(kind=4) :: j_mix !< Loop index for mixing waters
            integer(kind=4) :: num_mix_loc !< Number of mixing waters
            real(kind=8) :: lambda_r !< Placeholder for mixing ratio
            real(kind=8), allocatable :: c_mix(:) !< Transformed concentrations from mixing waters (temporary)
            real(kind=8), allocatable :: c1_old(:) !< Old primary species concentrations for current target water
            real(kind=8), allocatable :: conc_old_mix(:,:) !< conc_old from each mixing water
            integer(kind=4), allocatable :: ind_aq_sp_mix(:,:) !< indices_aq_species from each mixing water

            do i=1,this%num_target_waters !< Loop over all target waters to apply water mixing with reactions
                c1_old=this%waters(this%tar_wat_indices(i))%get_c1() !< Retrieve old primary species concentrations for domain water i from previous time step
                allocate(c_mix(this%waters(this%tar_wat_indices(i))%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)) !< Allocate c_mix with size = number of variable activity species in this water's reactive zone
                num_mix_loc=mix_wat_indices%cols(i)%dim-3
                allocate(conc_old_mix(size(this%waters(mix_wat_indices%cols(i)%col_1(1))%conc_old), num_mix_loc))
                allocate(ind_aq_sp_mix(size(this%waters(mix_wat_indices%cols(i)%col_1(1))%indices_aq_species), num_mix_loc))
                do j_mix=1,num_mix_loc
                    conc_old_mix(:,j_mix)=this%waters(mix_wat_indices%cols(i)%col_1(j_mix+1))%conc_old
                    ind_aq_sp_mix(:,j_mix)=this%waters(mix_wat_indices%cols(i)%col_1(j_mix+1))%indices_aq_species
                end do
                call compute_c_mix(this%waters(mix_wat_indices%cols(i)%col_1(1)),&
                    conc_old_mix,ind_aq_sp_mix,&
                    lambdas%cols(i)%col_1,c_mix) !< Compute mixed transformed concentrations by applying weighted mixing to old target waters using provided indices and lambda ratios
                deallocate(conc_old_mix,ind_aq_sp_mix)
                lambda_r=lambdas%cols(i)%col_1(1) !< Extract mixing ratio (lambda_r) for current domain water i from lambdas matrix
                call this%waters(this%tar_wat_indices(i))%reactive_mixing_iter_EE_eq_kin_ideal( &
                    c1_old,c_mix,lambda_r,Delta_t,0d0,c_new(:,i),u_new(:,i)) !< Solve water mixing with equilibrium and kinetic reactions using explicit Euler, ideal activity model, and lumped formulation for domain water i
                deallocate(c_mix,c1_old) !< Free temporary arrays before next iteration
            end do !< End loop over domain target waters
        end subroutine WMA_iter_EE_eq_kin_ideal_lump

        !> @brief Check for new reactive zones due to concentration changes
        !> @details Monitors whether concentration changes (e.g., mineral dissolution to zero)
        !>          have created new reactive zones. Scans target solids for zero-concentration
        !>          non-flowing species and identifies unique patterns to create new zones.
        !>          Uses pairwise comparison to detect distinct zero-concentration patterns.
        !> @param[in,out] this Chemistry object with reactive zones to check/update
        !> @param[in] i Index of the reactive zone to check for splitting
        !> @param[in] tolerance Tolerance for detecting zero concentrations of non-flowing species
        subroutine check_new_reactive_zones(this,i,tolerance)
            class(chemistry_c) :: this                   !< Chemistry object with reactive zones
            integer(kind=4), intent(in)  :: i             !< Reactive zone index to check
            real(kind=4), intent(in) :: tolerance         !< Tolerance for zero concentration detection
    
            integer(kind=4) :: j,k,l,p,r,s,m,ii,kk,ll,old_num_react_zones,num_new_react_zones,num_sol_conc_zero,num_nf_sp_new,num_tar_sol,& !< Local loop indices and counters
            num_rep_rows                                  !< Number of repeated (duplicate) rows in zero-flag matrix
            integer(kind=4), allocatable :: tar_sol_indices(:),old_conc_zero_flag(:,:),new_conc_zero_flag(:,:),ind_sol_conc_zero(:),& !< Allocatable working arrays
            ind_conc_zero(:,:),ind_rep_rows(:)            !< Indices of zero-concentration species and repeated rows
            type(reactive_zone_c), allocatable :: new_react_zones(:),old_react_zones(:) !< Temporary reactive zone arrays
    
            old_react_zones=this%reactive_zones           !< Save current reactive zones before modification
            old_num_react_zones=this%num_reactive_zones   !< Save current count of reactive zones
    
            call this%link_target_solids_reactive_zone(i,tar_sol_indices) !< Get target solid indices associated with reactive zone i
            num_tar_sol=size(tar_sol_indices)             !< Count target solids in this reactive zone
            num_sol_conc_zero=0                           !< Initialize counter for solids with zero-concentration species
    
            !> Count target solids with at least one zero-concentration non-flowing species
            do j=1,num_tar_sol                           !< Loop over target solids in this reactive zone
                do k=1,this%reactive_zones(i)%num_non_flow_species !< Loop over non-flowing species in this zone
                    if (abs(this%target_solids(tar_sol_indices(j))%concentrations(k))<tolerance) then !< Check if species concentration is below tolerance
                        num_sol_conc_zero=num_sol_conc_zero+1 !< Increment counter for solids with zero concentrations
                        exit                             !< Exit inner loop (only need to find one zero species per solid)
                    end if                               !< End concentration check
                end do                                   !< End loop over non-flowing species
            end do                                       !< End loop over target solids
            allocate(ind_sol_conc_zero(num_sol_conc_zero)) !< Allocate array for indices of solids with zero concentrations
            allocate(old_conc_zero_flag(num_sol_conc_zero,this%reactive_zones(i)%num_non_flow_species)) !< Allocate binary flag matrix [solids x species]
            l=0                                          !< Initialize index for filling ind_sol_conc_zero array
            do j=1,num_tar_sol                           !< Second pass: populate array of zero-concentration solid indices
                do k=1,this%reactive_zones(i)%num_non_flow_species !< Loop over non-flowing species
                    if (abs(this%target_solids(tar_sol_indices(j))%concentrations(k))<tolerance) then !< Check for zero concentration
                        l=l+1                            !< Advance index in result array
                        ind_sol_conc_zero(l)=tar_sol_indices(j) !< Store target solid index with zero concentration
                        if (l==num_sol_conc_zero) exit   !< All zero-concentration solids found, exit early
                    end if                               !< End concentration check
                end do                                   !< End loop over species
            end do                                       !< End loop over target solids
            do l=1,num_sol_conc_zero                     !< Build binary flag matrix: 1=zero concentration, 0=nonzero
                do k=1,this%reactive_zones(i)%num_non_flow_species !< Loop over non-flowing species for each zero-concentration solid
                    if (abs(this%target_solids(ind_sol_conc_zero(l))%concentrations(k))<tolerance) then !< Check if species k is zero in solid l
                        old_conc_zero_flag(l,k)=1        !< Mark species as zero (depleted)
                    else                                 !< Species concentration is above tolerance
                        old_conc_zero_flag(l,k)=0        !< Mark species as nonzero (present)
                    end if                               !< End concentration check
                end do                                   !< End loop over species
            end do                                       !< End loop over zero-concentration solids
            num_rep_rows=0                               !< Initialize counter for duplicate zero-pattern rows
            do l=1,num_sol_conc_zero-1                   !< Pairwise comparison: find duplicate zero-concentration patterns
                do m=l+1,num_sol_conc_zero               !< Compare row l with all subsequent rows
                    if (inf_norm_vec_int(old_conc_zero_flag(l,:)-old_conc_zero_flag(m,:))<tolerance) then !< Check if patterns are identical (infinity norm < tolerance)
                        num_rep_rows=num_rep_rows+1      !< Increment duplicate row counter
                        exit                             !< Found a match, no need to compare further
                    end if                               !< End pattern comparison
                end do                                   !< End inner comparison loop
            end do                                       !< End outer comparison loop
            allocate(ind_rep_rows(num_rep_rows))          !< Allocate array for indices of duplicate rows
            p=0                                          !< Initialize index for filling ind_rep_rows
            do l=1,num_sol_conc_zero-1                   !< Second pass: identify specific duplicate row indices
                do m=l+1,num_sol_conc_zero               !< Compare row l with subsequent rows
                    if (inf_norm_vec_int(old_conc_zero_flag(l,:)-old_conc_zero_flag(m,:))<tolerance) then !< Check for duplicate pattern
                        p=p+1                            !< Advance index in result array
                        ind_rep_rows(p)=l                !< Store index of duplicate row
                        exit                             !< Found a match, move to next row
                    end if                               !< End pattern comparison
                end do                                   !< End inner comparison loop
            end do                                       !< End outer comparison loop
            num_new_react_zones=num_sol_conc_zero-num_rep_rows !< Number of unique new reactive zones = total - duplicates
            allocate(new_react_zones(num_new_react_zones)) !< Allocate array for new unique reactive zones
            allocate(new_conc_zero_flag(num_new_react_zones,this%reactive_zones(i)%num_non_flow_species)) !< Allocate flag matrix for new unique zones
            p=1                                          !< Initialize pointer into ind_rep_rows array
            r=0                                          !< Initialize counter for unique (non-duplicate) rows
            do l=1,num_sol_conc_zero                     !< Filter out duplicate rows to get unique zero-concentration patterns
                if (l/=ind_rep_rows(p)) then              !< Current row is not a duplicate
                    r=r+1                                 !< Increment unique zone counter
                    new_conc_zero_flag(r,:)=old_conc_zero_flag(l,:) !< Copy unique zero-pattern to new flag matrix
                    if (p<num_rep_rows) then              !< Check if more duplicate indices remain
                        p=p+1                             !< Advance pointer to next duplicate index
                    else                                  !< All duplicates have been skipped
                        exit                             !< Exit filtering loop
                    end if                               !< End duplicate index check
                end if                                   !< End duplicate row check
            end do                                       !< End filtering loop
            do s=1,num_sol_conc_zero-ind_rep_rows(num_rep_rows) !< Copy remaining rows after the last duplicate index
                r=r+1                                    !< Increment unique zone counter
                new_conc_zero_flag(r,:)=old_conc_zero_flag(ind_rep_rows(num_rep_rows)+s,:) !< Copy tail rows to new flag matrix
            end do                                       !< End tail copy loop
            num_nf_sp_new=0                              !< Initialize non-flowing species count for new zones
            do ii=1,num_new_react_zones                  !< Loop over each new reactive zone to initialize its properties
                num_nf_sp_new=this%reactive_zones(i)%num_non_flow_species-sum(new_conc_zero_flag(ii,:)) !< Count remaining (nonzero) non-flowing species for this new zone
                call new_react_zones(ii)%allocate_ind_non_flow_species(num_nf_sp_new) !< Allocate non-flowing species index array for new zone
                call new_react_zones(ii)%set_chem_syst_react_zone(this%chem_syst) !< Link new zone to the chemical system
                ll=0                                     !< Initialize index for populating non-flowing species indices
                do kk=1,this%reactive_zones(i)%num_non_flow_species !< Loop over all non-flowing species in parent zone
                    if (new_conc_zero_flag(ii,kk)==0) then !< Species is still present (nonzero concentration)
                        ll=ll+1                          !< Increment non-flowing species index
                        new_react_zones(ii)%ind_non_flow_species(ll)=this%reactive_zones(i)%ind_non_flow_species(kk) !< Copy species index from parent zone
                    end if                               !< End nonzero species check
                end do                                   !< End loop over parent zone species
            end do                                       !< End loop over new reactive zones
        
            deallocate(this%reactive_zones)               !< Free old reactive zones array before reallocation
            call this%allocate_reactive_zones(old_num_react_zones+num_new_react_zones) !< Reallocate with space for both old and new zones
        
            do l=1,old_num_react_zones                   !< Restore original reactive zones into expanded array
                this%reactive_zones(l)=old_react_zones(l) !< Copy old zone l back to position l
            end do                                       !< End restoration loop
            do l=1,num_new_react_zones                   !< Append new reactive zones after the old ones
                this%reactive_zones(old_num_react_zones+l)=new_react_zones(l) !< Place new zone l at position old_count+l
            end do                                       !< End append loop
            do l=1,this%num_reactive_zones               !< Keep each zone id equal to its index in this%reactive_zones
                this%reactive_zones(l)%id=l
            end do
        
            deallocate(tar_sol_indices,old_conc_zero_flag,new_conc_zero_flag,ind_sol_conc_zero,new_react_zones,old_react_zones) !< Free all temporary working arrays
        end subroutine
        
        !> @brief Main reactive mixing solver dispatcher
        !> @details Selects and calls the appropriate reactive mixing solver based on the
        !>          activity coefficient model and spatial lumping flag configuration:
        !>          - act_coeffs_model=0 + lump=.true.:  solve_reactive_mixing_ideal_lump
        !>          - act_coeffs_model=0 + lump=.false.: solve_reactive_mixing_ideal_cons
        !>          - act_coeffs_model>0 + lump=.true.:  solve_reactive_mixing_lump (not yet implemented)
        !>          - act_coeffs_model>0 + lump=.false.: solve_reactive_mixing_cons (not yet implemented)
        !> @param[in,out] this Chemistry object coordinating the simulation
        !> @param[in] dir_pb Problem directory path
        !> @param[in] root Root name for output files
        !> @param[in] mixing_ratios_conc Mixing ratios for concentrations
        !> @param[in] mixing_ratios_R Mixing ratios for reaction rates (same layout as mixing_ratios_conc; 0 for conc-only waters)
        !> @param[in] mix_conc_indices Matrix of indices of target waters that mix
        !> @param[in] time_discr Time discretization object
        !> @param[in] theta_r Integration method parameter for chemical reactions
        subroutine solve_reactive_mixing(this,dir_pb,root,mixing_ratios_conc,mixing_ratios_R, &
            mix_conc_indices,time_discr,theta_r)
            class(chemistry_c) :: this                    !< Chemistry object coordinating the simulation
            character(len=*), intent(in) :: dir_pb        !< Problem directory path
            character(len=*), intent(in) :: root          !< Root name for output files
            class(real_array_c), intent(in) :: mixing_ratios_conc !< Mixing ratios for concentrations
            class(real_array_c), intent(in) :: mixing_ratios_R !< Mixing ratios for reaction rates (mix_conc_indices layout, 0 for conc-only waters)
            class(int_array_c), intent(in) :: mix_conc_indices !< Matrix of indices of target waters that mix
            class(time_discr_c), intent(in) :: time_discr !< Time discretization object
            real(kind=8), intent(in) :: theta_r           !< Integration method parameter for chemical reactions
            
            procedure(solve_reactive_mixing_ideal_cons), pointer :: p_solver=>null() !< Procedure pointer to selected reactive mixing solver
            
            if (this%act_coeffs_model==0 .and. this%lump_flag) then !< Check for ideal activity with spatial lumping
                p_solver=>solve_reactive_mixing_ideal_lump !< Select ideal lumped solver
            else if (this%act_coeffs_model==0 .and. (.not. this%lump_flag)) then !< Check for ideal activity without lumping
                p_solver=>solve_reactive_mixing_ideal_cons !< Select ideal consistent solver
            else if (this%act_coeffs_model>0 .and. this%lump_flag) then !< Check for non-ideal activity (extended Debye-Hückel, Pitzer, etc.) with lumping
                p_solver=>solve_reactive_mixing_lump    
            else if (this%act_coeffs_model>0 .and. (.not. this%lump_flag)) then !< Check for non-ideal activity without lumping (most accurate, slowest)
                p_solver=>solve_reactive_mixing_cons !< Select non-ideal consistent solver (uniform 8-argument signature)
            else                                                               !< Invalid combination of act_coeffs_model and lump_flag
                error stop "Error: act_coeffs_model and lump_flag not compatible" !< Terminate with error: check chemistry configuration
            end if
            call p_solver(this,dir_pb,root,mixing_ratios_conc,& !< Call the selected reactive mixing solver procedure
                mixing_ratios_R,mix_conc_indices,& !< Pass reaction mixing ratios and water indices
                time_discr,theta_r) !< Pass time discretization and integration method
        end subroutine

        !> @brief Writes one time step of the domain-integrated reactive-mixing mass balance (per variable-activity species).
        !> @details idx_ref is the reference target water (homogeneous domain: it supplies the species names/flags).
        !> For each species: LHS (storage), RHS = transport + reaction, residual = LHS - RHS (all already summed over cells).
        !> A residual above mb_tol*max(1,|LHS|,|RHS|) is a hard failure (error stop) after printing the offending species
        !> and its LHS/transport/reaction split: a violated balance means the reactive-mixing step is wrong, not a warning.
        subroutine write_react_mix_mass_bal_dom(this,idx_ref,iunit,time,n_v,lhs,transp,react)
            class(chemistry_c), intent(in) :: this !< Chemistry object instance
            integer(kind=4), intent(in) :: idx_ref !< reference target-water index (species names/flags)
            integer(kind=4), intent(in) :: iunit !< mass-balance output file unit
            real(kind=8), intent(in) :: time !< current simulation time
            integer(kind=4), intent(in) :: n_v !< number of variable activity species
            real(kind=8), intent(in) :: lhs(:),transp(:),react(:) !< domain-summed per-species terms
            integer(kind=4) :: j,nmax
            real(kind=8) :: rhs,residual,mb_tol,mb_scale
            character(len=:), allocatable :: sp_name
            write(iunit,"(2x,'t = ',ES15.5)") time
            write(iunit,"(4x,A15,3A20)") 'species','LHS','RHS','residual'
            !> Species names/flags come from the aqueous phase via indices_aq_phase o ind_var_act_species;
            !> chem_syst%species is reordered by the speciation algebra, so indexing it directly with
            !> ind_var_act_species mislabels and drops species (e.g. ca+2 for gypsum).
            mb_tol=sqrt(this%CV_params%abs_tol) !< mass-balance residuals must be tighter than sqrt(abs_tol)
            associate(w => this%waters(idx_ref))
            nmax=min(n_v,size(w%ind_var_act_species))
            do j=1,nmax
                if (w%aq_phase%aq_species(w%indices_aq_phase(w%ind_var_act_species(j)))%cst_act_flag .eqv. .false.) then
                    rhs=transp(j)+react(j)
                    residual=lhs(j)-rhs
                    sp_name=trim(w%aq_phase%aq_species(w%indices_aq_phase(w%ind_var_act_species(j)))%name)
                    write(iunit,"(4x,A15,3ES20.10)") sp_name,lhs(j),rhs,residual
                    !> Scale-relative violation test (absolute floor at 1) so the guard is meaningful in
                    !> any solver regardless of the balance magnitude: |residual| > mb_tol*max(1,|LHS|,|RHS|).
                    mb_scale=max(abs(lhs(j)),abs(rhs))
                    if (abs(residual)>mb_tol*max(1d0,mb_scale)) then
                        write(*,"(A,ES15.5,A,A,A,ES15.5,A,ES15.5,A,ES15.5,A)") "ERROR: reactive-mixing mass balance violated at t=",time,&
                            ", species '",sp_name,"': |residual|=",abs(residual)," > tol=",mb_tol*max(1d0,mb_scale),&
                            " (|LHS|=",abs(lhs(j)),")"
                        write(*,"(A,I0,A,ES15.5,A,ES15.5,A,ES15.5)") "   [MB-SPLIT] j=",j," lhs=",lhs(j), &
                            " transp=",transp(j)," react=",react(j)
                        error stop "Reactive-mixing mass balance violated (residual exceeds tolerance)"
                    end if
                end if
            end do
            end associate
            !> Companion component-space balance (u = U*c): apply the component matrix U to the
            !> species terms. Equilibrium reactions cancel under U (U*Se_eq^T = 0), so its reaction
            !> term is kinetic-only. Written to a sibling "<name>_comp.out", opened lazily here and
            !> kept open across time steps (re-found each call via INQUIRE, no extra arguments needed).
            block
                integer(kind=4) :: ncomp,iunit_comp,jc
                logical :: comp_open
                character(len=512) :: fname
                character(len=:), allocatable :: comp_fname
                character(len=15) :: clabel
                real(kind=8), allocatable :: lhs_c(:),transp_c(:),react_c(:)
                real(kind=8) :: rhs_c,resid_c
                associate(sa => this%waters(idx_ref)%solid_chemistry%reactive_zone%speciation_alg)
                ncomp=size(sa%comp_mat,1)
                lhs_c=matmul(sa%comp_mat,lhs(1:n_v))
                transp_c=matmul(sa%comp_mat,transp(1:n_v))
                react_c=matmul(sa%comp_mat,react(1:n_v))
                end associate
                inquire(unit=iunit,name=fname)
                comp_fname=trim(fname)
                comp_fname=comp_fname(1:len_trim(comp_fname)-4)//'_comp.out'
                inquire(file=comp_fname,opened=comp_open,number=iunit_comp)
                if (.not. comp_open) then
                    open(newunit=iunit_comp,file=comp_fname,status='replace',action='write')
                    write(iunit_comp,*) 'Domain-integrated reactive-mixing mass balance (per aqueous component u = U*c)'
                    write(iunit_comp,*) 'Per time step & component: LHS = U*sum V*(c^{k+1}-c^k)/Dt, RHS = U*(transport+reaction), residual = LHS - RHS'
                end if
                write(iunit_comp,"(2x,'t = ',ES15.5)") time
                write(iunit_comp,"(4x,A15,3A20)") 'component','LHS','RHS','residual'
                do jc=1,ncomp
                    rhs_c=transp_c(jc)+react_c(jc)
                    resid_c=lhs_c(jc)-rhs_c
                    write(clabel,'(A,I0)') 'u',jc
                    write(iunit_comp,"(4x,A15,3ES20.10)") clabel,lhs_c(jc),rhs_c,resid_c
                end do
            end block
        end subroutine

        !> @brief Set upstream water indices using flow direction and grid topology
        !> @details For each target water, identifies the closest upstream water as the
        !>          mixing contributor whose displacement is most anti-parallel to the
        !>          local flow velocity (most negative dot product).
        !>          Waters are assumed ordered left-to-right, bottom-to-top on a grid
        !>          with nx columns. For 1D, set nx = num_waters and vy = 0.
        subroutine set_upstream_water_indices(this, mix_conc_indices, nx, vx, vy)
            class(chemistry_c) :: this                   !< Chemistry object instance
            class(int_array_c), intent(in) :: mix_conc_indices !< Mixing water indices for each water
            integer(kind=4), intent(in) :: nx             !< Number of waters per row in grid (for 1D: num_waters)
            real(kind=8), intent(in) :: vx(:)             !< x-velocity at each target water (dim: num_target_waters)
            real(kind=8), intent(in) :: vy(:)             !< y-velocity at each target water (dim: num_target_waters)
            integer(kind=4) :: i, j, k, w, num_mix, ix, iy, jx, jy, max_self !< Local loop indices and grid coordinates
            real(kind=8) :: dx, dy, dot_prod, best_dot    !< Displacement components, dot product, and best match tracker
            integer(kind=4), allocatable :: self_to_tpos(:) !< Inverse map: self-index -> target position
            if (allocated(this%upstream_water_indices)) deallocate(this%upstream_water_indices) !< Safe deallocation before reallocation
            allocate(this%upstream_water_indices(this%num_target_waters)) !< Allocate upstream indices array
            !> Build inverse map from mix_conc_indices self-index to target position.
            !> mix_conc_indices%cols(k)%col_1(1) is the self-index of target k,
            !> which may be a transport-local index (1D) or a global water index (2D after remap).
            max_self = 0                                  !< Initialize maximum self-index to zero
            do i = 1, this%num_target_waters             !< First pass: find maximum self-index across all target waters
                w = mix_conc_indices%cols(i)%col_1(1)     !< Get self-index for target water i
                if (w > max_self) max_self = w            !< Update maximum self-index if current is larger
            end do                                       !< End first pass
            allocate(self_to_tpos(max_self))              !< Allocate inverse map array with size = max self-index
            self_to_tpos = 0                              !< Initialize all mappings to zero (no target position)
            do i = 1, this%num_target_waters             !< Second pass: build inverse map from self-index to position
                self_to_tpos(mix_conc_indices%cols(i)%col_1(1)) = i !< Map self-index of target i to position i
            end do                                       !< End inverse map construction
            do i=1,this%num_target_waters                 !< Main loop: find upstream water for each target water
                num_mix=mix_conc_indices%cols(i)%dim-3    !< Number of mixing contributors for target water i
                ix=mod(i-1,nx)+1                          !< Compute column position of target water i on grid
                iy=(i-1)/nx+1                             !< Compute row position of target water i on grid
                best_dot=-huge(1d0)                       !< Initialize best dot product to -infinity (closest upstream = least negative)
                this%upstream_water_indices(i)=0          !< Default: 0 = no upstream target (boundary water)
                do j=1,num_mix                            !< Loop over mixing contributors of target water i
                    w=mix_conc_indices%cols(i)%col_1(j+1) !< Get self-index of mixing contributor j
                    if (w < 1 .or. w > max_self) cycle    !< Skip invalid self-indices
                    k=self_to_tpos(w)                     !< Map self-index to target position
                    if (k == 0) then                      !< Self-index not in target list (boundary water)
                        !> w is a boundary water: use grid position 0 (upstream of target 1)
                        jx=0                              !< Boundary water grid column = 0
                        jy=iy                             !< Boundary water grid row = same as target
                    else                                  !< Self-index is a target water
                        jx=mod(k-1,nx)+1                  !< Compute column of mixing contributor on grid
                        jy=(k-1)/nx+1                     !< Compute row of mixing contributor on grid
                    end if                                !< End boundary vs target check
                    dx=dble(jx-ix)                        !< Compute x-displacement from target to contributor
                    dy=dble(jy-iy)                        !< Compute y-displacement from target to contributor
                    dot_prod=dx*vx(i)+dy*vy(i)            !< Compute dot product of displacement with flow velocity
                    if (dot_prod<0d0 .and. dot_prod>best_dot) then !< Check if contributor is upstream and closer than current best
                        best_dot=dot_prod                 !< Update best dot product
                        this%upstream_water_indices(i)=k   !< Update upstream index (0 for boundary, 1..N for target)
                    end if                                !< End upstream candidate check
                end do                                    !< End loop over mixing contributors
            end do                                        !< End main loop over target waters
            deallocate(self_to_tpos)                       !< Free inverse map array
        end subroutine set_upstream_water_indices

        !> @brief Set downstream water indices using flow direction and grid topology
        !> @details For each target water, identifies the closest downstream water as the
        !>          water (that receives from it) whose displacement is most parallel
        !>          to the local flow velocity (most positive dot product).
        !>          Waters are assumed ordered left-to-right, bottom-to-top on a grid
        !>          with nx columns. For 1D, set nx = num_target_waters and vy = 0.
        subroutine set_downstream_water_indices(this, mix_conc_indices, nx, vx, vy)
            class(chemistry_c) :: this                   !< Chemistry object instance
            class(int_array_c), intent(in) :: mix_conc_indices !< Mixing water indices for each water
            integer(kind=4), intent(in) :: nx             !< Number of waters per row in grid
            real(kind=8), intent(in) :: vx(:)             !< x-velocity at each target water (dim: num_target_waters)
            real(kind=8), intent(in) :: vy(:)             !< y-velocity at each target water (dim: num_target_waters)
            integer(kind=4) :: i, j, k, num_mix, ix, iy, jx, jy, self_i !< Local loop indices and grid coordinates
            real(kind=8) :: dx, dy, dot_prod, best_dot    !< Displacement components, dot product, and best match tracker
            if (allocated(this%downstream_water_indices)) deallocate(this%downstream_water_indices) !< Safe deallocation before reallocation
            allocate(this%downstream_water_indices(this%num_target_waters)) !< Allocate downstream indices array
            do i=1,this%num_target_waters                 !< Main loop: find downstream water for each target water
                self_i=mix_conc_indices%cols(i)%col_1(1)  !< Get self-index (transport-local or global) for target water i
                ix=mod(i-1,nx)+1                          !< Compute column position of target water i on grid
                iy=(i-1)/nx+1                             !< Compute row position of target water i on grid
                best_dot=huge(1d0)                        !< Initialize best dot product to +infinity (closest downstream = least positive)
                this%downstream_water_indices(i)=this%tar_wat_indices(i) !< Default: self (global water index, no downstream neighbor)
                !> Scan all other target waters to find which has i as a mixing contributor
                do j=1,this%num_target_waters             !< Loop over all target waters to find downstream neighbors
                    if (j==i) cycle                       !< Skip self-comparison
                    num_mix=mix_conc_indices%cols(j)%dim-3 !< Number of mixing contributors for target water j
                    do k=1,num_mix                        !< Loop over mixing contributors of water j
                        if (mix_conc_indices%cols(j)%col_1(k+1)==self_i) then !< Check if water i contributes to water j
                            jx=mod(j-1,nx)+1              !< Compute column of downstream candidate on grid
                            jy=(j-1)/nx+1                 !< Compute row of downstream candidate on grid
                            dx=dble(jx-ix)                !< Compute x-displacement from target to candidate
                            dy=dble(jy-iy)                !< Compute y-displacement from target to candidate
                            dot_prod=dx*vx(i)+dy*vy(i)    !< Compute dot product of displacement with flow velocity
                            if (dot_prod>0d0 .and. dot_prod<best_dot) then !< Check if candidate is downstream and closer than current best
                                best_dot=dot_prod         !< Update best dot product
                                this%downstream_water_indices(i)=this%tar_wat_indices(j) !< Update downstream index
                            end if                        !< End downstream candidate check
                            exit                          !< self_i appears at most once in j's mix list
                        end if                            !< End contributor match check
                    end do                                !< End loop over mixing contributors
                end do                                    !< End loop over all target waters
            end do                                        !< End main loop over target waters
            !> Assign outflow boundary water as downstream for target waters with no downstream neighbor
            do i=1,this%num_target_waters                 !< Loop to assign boundary water for targets without downstream neighbor
                if (this%downstream_water_indices(i)==this%tar_wat_indices(i)) then !< Check if still pointing to self (no downstream found)
                    this%downstream_water_indices(i)=this%bd_waters_indices(this%num_bd_waters) !< Assign outflow boundary water as downstream
                end if                                    !< End self-check
            end do                                        !< End boundary assignment loop
        end subroutine set_downstream_water_indices

        !> @brief Compute boundary water concentrations using Lagrangian extrapolation
        !> @details Computes boundary water concentrations by linear extrapolation from
        !>          water types and nearest target waters. Uses the formula:
        !>          c_bd = 2 * c_type - c_nearest_target
        !>          This ensures continuity of concentration gradients at boundaries.
        !> @param[in,out] this Chemistry object containing boundary and target waters
        subroutine compute_conc_bd_waters_Lagr(this)
            class(chemistry_c) :: this                    !< Chemistry object instance

            this%waters(this%bd_waters_indices(1)&
                    )%concentrations=2d0*this%wat_types(1)%concentrations - & !< Left boundary: extrapolate from first water type and first target water
                    this%waters(this%tar_wat_indices(1))%concentrations !< c_bd_left = 2*c_type_1 - c_target_1
                this%waters(this%bd_waters_indices(2)&
                    )%concentrations=2d0*this%wat_types(this%num_wat_types& !< Right boundary: extrapolate from last water type and last target water
                    )%concentrations - &
                    this%waters(this%tar_wat_indices(&
                    this%num_target_waters))%concentrations !< c_bd_right = 2*c_type_N - c_target_N
        end subroutine compute_conc_bd_waters_Lagr

        !> @brief Reorder mixing reaction indices for upstream/downstream processing
        !> @details Reorders the mixing reaction indices array to prioritize upstream and
        !>          downstream target waters. This ordering is used by the reactive mixing
        !>          solver to improve convergence of the iterative algorithm.
        !> @param[in,out] this Chemistry object
        !> @param[in] ind_tw Index of target water in the waters attribute
        !> @param[in] old_mix_react_indices Original mixing reaction indices array
        !> @param[in] num_tar_wat_up Number of upstream target waters
        !> @param[in] num_tar_wat_down Number of downstream target waters
        !> @param[out] new_mix_react_indices Reordered mixing reaction indices array
        subroutine reorder_mix_react_indices(this,ind_tw,&
            old_mix_react_indices,num_tar_wat_up,num_tar_wat_down,&
            new_mix_react_indices)
            class(chemistry_c) :: this                    !< Chemistry object instance
            integer(kind=4), intent(in) :: ind_tw         !< Index of target water in waters attribute
            integer(kind=4), intent(in) :: old_mix_react_indices(:) !< Original mixing reaction indices array
            integer(kind=4), intent(in) :: num_tar_wat_up !< Number of upstream target waters
            integer(kind=4), intent(in) :: num_tar_wat_down !< Number of downstream target waters
            integer(kind=4), intent(out) :: new_mix_react_indices(:) !< Output: reordered mixing reaction indices
            integer(kind=4) :: i,num_tar_wat              !< Loop index and total number of target waters
            num_tar_wat=size(old_mix_react_indices)       !< Get total number of mixing target waters from input array size
        end subroutine reorder_mix_react_indices

        !> @brief Set estimation parameter for downstream kinetic reactions
        !> @details Sets the estimation parameter used for extrapolating reaction rates
        !>          at downstream locations. Must be non-negative.
        !>          - est_prm = 0: No extrapolation (use local rates only)
        !>          - est_prm > 0: Weighted extrapolation from upstream rates
        !> @param[in,out] this Chemistry object to configure
        !> @param[in] est_prm Estimation parameter [-] (must be >= 0)
        subroutine set_est_prm(this, est_prm)
            class(chemistry_c) :: this                    !< Chemistry object instance
            real(kind=8), intent(in) :: est_prm           !< Estimation parameter [-] (must be >= 0)
            if (est_prm<0d0) then                         !< Validate estimation parameter is non-negative
                error stop "Error: est_prm must be non-negative" !< Terminate if negative
            end if                                        !< End validation
            this%est_prm=est_prm                           !< Assign validated estimation parameter
        end subroutine set_est_prm

        !> @brief Writes one output-step block honouring chem_out_options%vars.
        !> @details For each selected variable type ('conc','react_rate','min-to-pore') writes a
        !>          labelled table: a header row ('Water index' + column names taken from the
        !>          reference output water) followed by one row per requested output water.
        !>          react_rate columns = [aqueous eq | mineral/exch eq | aqueous kin | mineral kin]
        !>          (end-of-step instantaneous rates). min-to-pore columns = mineral-to-pore volume ratios.
        subroutine write_out_step_vars(chem, unit, time)
            class(chemistry_c), intent(in) :: chem          !< Chemistry object (waters + output options)
            integer(kind=4), intent(in) :: unit             !< Output file unit
            real(kind=8), intent(in) :: time                !< Simulation time to label the block
            integer(kind=4) :: iv, ii, j, jr, m, g, idx, ridx
            integer(kind=4) :: nre, nm, nrct, nvf, nk
            integer(kind=4) :: n_v, n_v_w, n_aq_prim, n_prim, n_prim_w, n_lin, n_redox, n_ref_out
            integer(kind=4) :: n_vas, nterm, maxlen, nsc, off
            integer(kind=4) :: im !< mix-variable loop: 1 = conc_mix (conservative mix), 2 = conc_hat (speciation input)
            logical :: has_mix_var !< requested flag of the mix variable im
            logical :: mix_var_alloc !< the water's mix variable im is allocated
            integer(kind=4) :: mix_var_size !< its size
            logical :: has_conc, has_rate, has_vf, has_conc_mix, has_conc_hat, has_rxn_amnt, has_pH, has_sat
            logical, save :: header_written = .false.  !< file title/description written once (first call of the run)
            real(kind=8) :: val, coeff
            real(kind=8), allocatable :: rowvals(:), conc_nc(:), conc_comp(:)
            real(kind=8), allocatable :: rct_fact(:) !< abs_cum_react_amnt: factor of each column (mineral molecular weight [kg/mol] for the mineral columns, 1 for the rest, written in mol)
            character(len=256), allocatable :: rct_names(:), vf_names(:), ref_names(:), wat_names(:), vas_names(:)
            character(len=1024) :: comb_line
            character(len=1024), allocatable :: comb_lines(:)
            character(len=64) :: term_str

            !> Decide which variable types to write (defaults to all if none selected)
            has_conc=.false.; has_rate=.false.; has_vf=.false.; has_conc_mix=.false.; has_conc_hat=.false.; has_rxn_amnt=.false.
            has_pH=.false.; has_sat=.false.
            if (allocated(chem%chem_out_options%vars) .and. chem%chem_out_options%num_vars>0) then
                do iv=1,chem%chem_out_options%num_vars
                    select case (trim(chem%chem_out_options%vars(iv)))
                    case ('conc');       has_conc=.true.
                    case ('react_rate'); has_rate=.true.
                    case ('min-to-pore');  has_vf=.true.
                    case ('conc_mix'); has_conc_mix=.true.
                    case ('conc_hat'); has_conc_hat=.true.
                    case ('react_amnt'); has_rxn_amnt=.true.
                    case ('pH'); has_pH=.true.
                    case ('saturation'); has_sat=.true.
                    end select
                end do
            else
                !> no variable list: same default as read_chem_out_options (conc_hat, not conc_mix)
                has_conc=.true.; has_rate=.true.; has_vf=.true.; has_conc_mix=.false.; has_conc_hat=.true.; has_rxn_amnt=.true.
                has_pH=.true.; has_sat=.true.
            end if
            !> pH is only defined when h+ is in the aqueous phase (compute_pH is a no-op otherwise)
            if (chem%chem_syst%aq_phase%ind_prot==0) has_pH=.false.
            !> the saturation block is written only for systems with kinetic minerals
            if (chem%chem_syst%num_minerals_kin==0) has_sat=.false.

            ridx = chem%chem_out_options%ind_waters(1)   !< reference water: defines column headers

            !> File title and description, written once before the first time-step block
            if (.not. header_written) then
                write(unit,"(2x,'REMIX - reactive transport chemistry output')")
                write(unit,"(2x,'For each output time step, one labelled block per requested variable:')")
                if (chem%chem_syst%cat_exch_zone%num_surf_compl>0) then
                    write(unit,"(2x,'  conc       : aqueous variable-activity species concentrations, then all component concentrations U*c, then surface-complex concentrations (all output waters)')")
                else
                    write(unit,"(2x,'  conc       : aqueous variable-activity species concentrations, then all component concentrations U*c (all output waters)')")
                end if
                write(unit,"(2x,'  conc_mix   : aqueous variable-activity species concentrations of the conservative mix (previous-step concentrations mixed; for the by-zones mixers including the previous-step kinetic term), then all component concentrations U*c_mix (target waters only, t>0)')")
                write(unit,"(2x,'  conc_hat   : aqueous variable-activity species concentrations of the mixed waters entering the Newton speciation (mix + previous-step reaction terms of all mixing waters + current-step reaction terms of the mixing waters except the water itself), then all component concentrations U*c_hat (target waters only, t>0)')")
                write(unit,"(2x,'  react_rate : equilibrium/kinetic reaction rates (target waters only)')")
                write(unit,"(2x,'  min-to-pore : mineral-to-pore volume ratios V_min/V_pore (target waters only)')")
                write(unit,"(2x,'  abs_cum_react_amnt : absolute cumulative reaction amounts (target waters only); mineral columns in kg (amount in mol times the molecular weight of the mineral, 0 if the database gives none), the other columns in mol')")
                if (chem%chem_syst%aq_phase%ind_prot/=0) then
                    write(unit,"(2x,'  pH        : pH of the solution (all output waters)')")
                end if
                if (chem%chem_syst%num_minerals_kin>0) then
                    write(unit,"(2x,'  saturation : saturation index Omega = IAP/K of each kinetic mineral of the chemical system at the end of the step (all output waters; 1 = equilibrium if never evaluated)')")
                end if
                write(unit,"(2x,'  Each block: header row (water index + species/component names) followed by one row per requested output water')")
                header_written = .true.
            end if
            write(unit,"(/,2x,'t = ',ES26.16,/)") time

            !> --- Concentrations (variable-activity aqueous species) ---
            !> Column IDENTITY is the species NAME, not its position: each water's own
            !> ind_var_act_species order can differ from the reference water's (per-reactive-zone
            !> swap for numerical stability, and different secondary species across zones, e.g.
            !> calcite zone eliminates co3-2 while anhydrite zone eliminates so4-2). Printing by
            !> raw position under the reference header silently SWAPS/mislabels those columns.
            if (has_conc) then
                write(unit,"(2x,'conc aq var act sp')")
                write(unit,"(A18)", advance='no') 'Water index'
                associate(rw => chem%waters(ridx))
                n_v = rw%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
                allocate(ref_names(n_v))
                n_ref_out = 0
                do j=1,n_v
                    if (.not. rw%aq_phase%aq_species(rw%indices_aq_phase(rw%ind_var_act_species(j)))%cst_act_flag) then
                        n_ref_out = n_ref_out+1
                        ref_names(n_ref_out) = rw%aq_phase%aq_species(rw%indices_aq_phase(rw%ind_var_act_species(j)))%name
                        write(unit,"(A26)", advance='no') trim(ref_names(n_ref_out))
                    end if
                end do
                end associate
                write(unit,*)
                do ii=1,chem%chem_out_options%num_waters
                    idx = chem%chem_out_options%ind_waters(ii)
                    associate(w => chem%waters(idx))
                    conc_nc = w%get_conc_var_act()
                    n_v_w = w%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
                    write(unit,"(I18)", advance='no') idx
                    do jr=1,n_ref_out
                        val = 0d0
                        do j=1,n_v_w
                            if (.not. w%aq_phase%aq_species(w%indices_aq_phase(w%ind_var_act_species(j)))%cst_act_flag) then
                                if (trim(w%aq_phase%aq_species(w%indices_aq_phase(w%ind_var_act_species(j)))%name)==&
                                    trim(ref_names(jr))) then
                                    !> conc_nc places aqueous secondaries AFTER the solid/gas primaries
                                    !> (free site) at num_prim_species+i, not contiguously after the aqueous
                                    !> primaries; map j accordingly so exchange systems don't report the free
                                    !> site in place of an eliminated aqueous secondary (e.g. calcite's ca+2).
                                    if (j<=w%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species) then
                                        m = j
                                    else
                                        m = w%solid_chemistry%reactive_zone%speciation_alg%num_prim_species + &
                                            (j - w%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)
                                    end if
                                    if (m<=size(conc_nc)) val = conc_nc(m)
                                    exit
                                end if
                            end if
                        end do
                        write(unit,"(ES26.16)", advance='no') val
                    end do
                    write(unit,*)
                    end associate
                    deallocate(conc_nc)
                end do
                deallocate(ref_names)
            end if

            if (has_pH) then
                write(unit,"(2x,'pH')")
                write(unit,"(A18)", advance='no') 'Water index'
                associate(rw => chem%waters(ridx))
                n_v = rw%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
                allocate(ref_names(n_v))
                n_ref_out = 0
                do j=1,n_v
                    if (.not. rw%aq_phase%aq_species(rw%indices_aq_phase(rw%ind_var_act_species(j)))%cst_act_flag) then
                        n_ref_out = n_ref_out+1
                        ref_names(n_ref_out) = rw%aq_phase%aq_species(rw%indices_aq_phase(rw%ind_var_act_species(j)))%name
                        !write(unit,"(A20)", advance='no') trim(ref_names(n_ref_out))
                    end if
                end do
                end associate
                write(unit,*)
                do ii=1,chem%chem_out_options%num_waters
                    idx = chem%chem_out_options%ind_waters(ii)
                    associate(w => chem%waters(idx))
                    conc_nc = w%get_conc_var_act()
                    n_v_w = w%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
                    write(unit, "(I18, ES26.16)") idx, w%pH
                    ! do jr=1,n_ref_out
                    !     val = 0d0
                    !     ! do j=1,n_v_w
                    !     !     if (.not. w%aq_phase%aq_species(w%indices_aq_phase(w%ind_var_act_species(j)))%cst_act_flag) then
                    !     !         if (trim(w%aq_phase%aq_species(w%indices_aq_phase(w%ind_var_act_species(j)))%name)==&
                    !     !             trim(ref_names(jr))) then
                    !     !             !> conc_nc places aqueous secondaries AFTER the solid/gas primaries
                    !     !             !> (free site) at num_prim_species+i, not contiguously after the aqueous
                    !     !             !> primaries; map j accordingly so exchange systems don't report the free
                    !     !             !> site in place of an eliminated aqueous secondary (e.g. calcite's ca+2).
                    !     !             if (j<=w%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species) then
                    !     !                 m = j
                    !     !             else
                    !     !                 m = w%solid_chemistry%reactive_zone%speciation_alg%num_prim_species + &
                    !     !                     (j - w%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)
                    !     !             end if
                    !     !             if (m<=size(conc_nc)) val = conc_nc(m)
                    !     !             exit
                    !     !         end if
                    !     !     end if
                    !     ! end do
                    !     write(unit,"(ES20.10)", advance='no') val
                    ! end do
                    !write(unit,*)
                    end associate
                    deallocate(conc_nc)
                end do
                deallocate(ref_names)
            end if

            !> --- Saturation indices Omega = IAP/K(T) of the KINETIC minerals ---
            !> Written only when the chemical system has kinetic minerals (has_sat); the columns are the
            !> kinetic minerals of the GLOBAL chem_syst list (its first num_minerals_kin minerals), resolved
            !> by NAME as in react_rate. Each water reports the kinetic entries of its saturation array
            !> (compute_saturation_mins at the end of the step, from the rate law, also for kinetic
            !> minerals absent from the water's mineral zone); the equilibrium-mineral entries of the array
            !> are kept but not written. Waters never evaluated (e.g. boundary waters) print 1
            !> (equilibrium, the attribute's initial value).
            if (has_sat) then
                nk = chem%chem_syst%num_minerals_kin
                write(unit,"(2x,'saturation kin min')")
                write(unit,"(A18)", advance='no') 'Water index'
                do j=1,nk
                    write(unit,"(A26)", advance='no') trim(strip_trailing2(chem%chem_syst%minerals(j)%name))
                end do
                write(unit,*)
                allocate(rowvals(nk))
                do ii=1,chem%chem_out_options%num_waters
                    idx = chem%chem_out_options%ind_waters(ii)
                    rowvals=1d0 !> never evaluated: equilibrium (as the attribute's initial value)
                    associate(w => chem%waters(idx))
                    if (allocated(w%saturation_mins) .and. associated(w%solid_chemistry)) then
                        associate(wcs => w%solid_chemistry%reactive_zone%chem_syst)
                        do m=1,min(size(w%saturation_mins),wcs%num_minerals_kin) !> kinetic minerals come first in chem_syst%minerals
                            g = find_mineral_col(wcs%minerals(m)%name)
                            if (g>=1 .and. g<=nk) rowvals(g) = w%saturation_mins(m)
                        end do
                        end associate
                    end if
                    end associate
                    write(unit,"(I18)", advance='no') idx
                    do j=1,nk
                        write(unit,"(ES26.16)", advance='no') rowvals(j)
                    end do
                    write(unit,*)
                end do
                deallocate(rowvals)
            end if

            !> --- Aqueous variable-activity species concentrations of the two mix stages, one block each ---
            !> im=1: conc_mix, the conservative mix (compute_c_mix_cons_global / cross_zone_mix_core).
            !> im=2: conc_hat, the mixed waters entering the Newton speciation (the speciation input c_hat
            !>       stored by the reactive mixing solvers via set_conc_hat).
            !> Both attributes are in ind_var_act_species aqueous order so attr(j) matches var-act species j
            !> directly. Column identity is the species NAME (as in the conc block). Waters with no stored
            !> value print zeros. External (boundary/recharge) waters are skipped (they are not mixed), and
            !> nothing is written at the initial time step (t=0).
            do im=1,2
                if (im==1) then
                    has_mix_var = has_conc_mix
                else
                    has_mix_var = has_conc_hat
                end if
                if (.not. (has_mix_var .and. time>0d0)) cycle
                if (im==1) then
                    write(unit,"(2x,'conc mix aq var act sp')")
                else
                    write(unit,"(2x,'conc hat aq var act sp')")
                end if
                write(unit,"(A18)", advance='no') 'Water index'
                associate(rw => chem%waters(ridx))
                n_v = rw%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
                allocate(ref_names(n_v))
                n_ref_out = 0
                do j=1,n_v
                    if (.not. rw%aq_phase%aq_species(rw%indices_aq_phase(rw%ind_var_act_species(j)))%cst_act_flag) then
                        n_ref_out = n_ref_out+1
                        ref_names(n_ref_out) = rw%aq_phase%aq_species(rw%indices_aq_phase(rw%ind_var_act_species(j)))%name
                        write(unit,"(A26)", advance='no') trim(ref_names(n_ref_out))
                    end if
                end do
                end associate
                write(unit,*)
                do ii=1,chem%chem_out_options%num_waters
                    idx = chem%chem_out_options%ind_waters(ii)
                    if (allocated(chem%ext_waters_indices)) then
                        if (any(chem%ext_waters_indices==idx)) cycle  !< external waters are not mixed
                    end if
                    associate(w => chem%waters(idx))
                    n_v_w = w%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
                    if (im==1) then
                        mix_var_alloc = allocated(w%conc_mix)
                    else
                        mix_var_alloc = allocated(w%conc_hat)
                    end if
                    mix_var_size = 0
                    if (mix_var_alloc) then
                        if (im==1) then
                            mix_var_size = size(w%conc_mix)
                        else
                            mix_var_size = size(w%conc_hat)
                        end if
                    end if
                    write(unit,"(I18)", advance='no') idx
                    do jr=1,n_ref_out
                        val = 0d0
                        if (mix_var_alloc) then
                            do j=1,n_v_w
                                if (.not. w%aq_phase%aq_species(w%indices_aq_phase(w%ind_var_act_species(j)))%cst_act_flag) then
                                    if (trim(w%aq_phase%aq_species(w%indices_aq_phase(w%ind_var_act_species(j)))%name)==&
                                        trim(ref_names(jr))) then
                                        if (j<=mix_var_size) then
                                            if (im==1) then
                                                val = w%conc_mix(j)
                                            else
                                                val = w%conc_hat(j)
                                            end if
                                        end if
                                        exit
                                    end if
                                end if
                            end do
                        end if
                        write(unit,"(ES26.16)", advance='no') val
                    end do
                    write(unit,*)
                    end associate
                end do
                deallocate(ref_names)
            end do

            !> --- Component (total) concentrations U*c, printed PER WATER ---
            !> The component basis is ZONE-dependent: the number of components, their names, and
            !> their definition (rows of comp_mat) all change from one reactive zone to another
            !> (e.g. the calcite zone has u = {cl-, ca+2-co3-2(p), so4-2} while the anhydrite zone
            !> has u = {cl-, ca+2-so4-2, co3-2(p)}). A single shared header would therefore mislabel
            !> or drop columns. Instead, for each output water we print its OWN components: for each
            !> component we show its NAME (this water's primary species), its definition as a linear
            !> combination of that water's variable-activity species (u_i = sum_j comp_mat(i,j)*c_j),
            !> and its value u_i = (comp_mat*c_nc)_i. Only nonzero coefficients are shown.
            !> Boundary/recharge waters carry the reactive zone of their water type (see read_tar_wat_line),
            !> so they are printed in that zone's component basis like any other water.
            if (has_conc) then
                write(unit,"(2x,'conc comp')")
                do ii=1,chem%chem_out_options%num_waters
                    idx = chem%chem_out_options%ind_waters(ii)
                    associate(w => chem%waters(idx), &
                              sa => chem%waters(idx)%solid_chemistry%reactive_zone%speciation_alg)
                    conc_comp = w%compute_conc_comp(w%get_conc_var_act())
                    n_prim_w = sa%num_prim_species
                    n_vas    = sa%num_var_act_species
                    allocate(ref_names(n_prim_w))
                    call prim_species_names(w, n_prim_w, ref_names)
                    allocate(vas_names(n_vas))
                    call var_act_species_names(w, n_vas, vas_names)
                    write(unit,"(4x,'water ',I0)") idx
                    allocate(comb_lines(n_prim_w))
                    maxlen = 1
                    do jr=1,n_prim_w
                        comb_line = ''
                        nterm = 0
                        do j=1,n_vas
                            coeff = sa%comp_mat(jr,j)
                            if (abs(coeff) < 1d-12) cycle
                            if (abs(coeff-1d0) < 1d-12) then
                                write(term_str,'(A)') trim(vas_names(j))
                                if (nterm==0) then
                                    comb_line = trim(adjustl(term_str))
                                else
                                    comb_line = trim(comb_line)//' + '//trim(adjustl(term_str))
                                end if
                            else if (abs(coeff+1d0) < 1d-12) then
                                if (nterm==0) then
                                    comb_line = '- '//trim(vas_names(j))
                                else
                                    comb_line = trim(comb_line)//' - '//trim(vas_names(j))
                                end if
                            else
                                !> one decimal, and take the sign from the separator (never print
                                !> a doubled '+ +'): print |coeff| and choose ' + ' / ' - '
                                write(term_str,'(F0.1,1x,A)') abs(coeff), trim(vas_names(j))
                                if (nterm==0) then
                                    if (coeff < 0d0) then
                                        comb_line = '- '//trim(adjustl(term_str))
                                    else
                                        comb_line = trim(adjustl(term_str))
                                    end if
                                else if (coeff < 0d0) then
                                    comb_line = trim(comb_line)//' - '//trim(adjustl(term_str))
                                else
                                    comb_line = trim(comb_line)//' + '//trim(adjustl(term_str))
                                end if
                            end if
                            nterm = nterm + 1
                        end do
                        comb_lines(jr) = comb_line
                        maxlen = max(maxlen, len_trim(comb_line))
                    end do
                    !> left-justify every combination to the longest one in this water block so the
                    !> value column stays aligned while no expression is ever truncated
                    do jr=1,n_prim_w
                        val = 0d0
                        if (jr<=size(conc_comp)) val = conc_comp(jr)
                        write(unit,"(6x,A24,' = ',A,' = ',ES26.16)") &
                            trim(ref_names(jr)), comb_lines(jr)(1:maxlen), val
                    end do
                    deallocate(comb_lines)
                    end associate
                    deallocate(conc_comp)
                    deallocate(ref_names)
                    if (allocated(vas_names)) deallocate(vas_names)
                end do
            end if

            !> --- Component concentrations of the two mix stages, u = U*c, per water, one block each ---
            !> im=1: 'conc comp mix' = U*c_mix (conservative mix); im=2: 'conc comp hat' = U*c_hat (speciation
            !> input). Same per-water component basis as the conc comp block (all components, incl. the
            !> exchange free site). The mix vector is rebuilt in the full conc_nc layout as the mixer/solver
            !> forms it: aqueous slots from the attribute ([aq prim | aq sec]), non-aqueous variable-activity
            !> slots from the water's previous-step state (get_conc_old_var_act_new_rz). Skipped for external waters and at t=0.
            do im=1,2
                if (im==1) then
                    has_mix_var = has_conc_mix
                else
                    has_mix_var = has_conc_hat
                end if
                if (.not. (has_mix_var .and. time>0d0)) cycle
                if (im==1) then
                    write(unit,"(2x,'conc comp mix')")
                else
                    write(unit,"(2x,'conc comp hat')")
                end if
                do ii=1,chem%chem_out_options%num_waters
                    idx = chem%chem_out_options%ind_waters(ii)
                    if (allocated(chem%ext_waters_indices)) then
                        if (any(chem%ext_waters_indices==idx)) cycle
                    end if
                    associate(w => chem%waters(idx), &
                              sa => chem%waters(idx)%solid_chemistry%reactive_zone%speciation_alg)
                    if (im==1) then
                        mix_var_alloc = allocated(w%conc_mix)
                    else
                        mix_var_alloc = allocated(w%conc_hat)
                    end if
                    if (.not. mix_var_alloc) cycle
                    if (im==1) then
                        mix_var_size = size(w%conc_mix)
                    else
                        mix_var_size = size(w%conc_hat)
                    end if
                    n_aq_prim = sa%num_aq_prim_species
                    n_prim    = sa%num_prim_species
                    n_prim_w  = sa%num_prim_species
                    n_vas     = sa%num_var_act_species
                    if (mix_var_size/=sa%num_aq_var_act_species) cycle
                    conc_nc = w%get_conc_old_var_act_new_rz()
                    if (im==1) then
                        conc_nc(1:n_aq_prim) = w%conc_mix(1:n_aq_prim)
                        conc_nc(n_prim+1:n_prim+sa%num_aq_sec_var_act_species) = &
                            w%conc_mix(n_aq_prim+1:sa%num_aq_var_act_species)
                    else
                        conc_nc(1:n_aq_prim) = w%conc_hat(1:n_aq_prim)
                        conc_nc(n_prim+1:n_prim+sa%num_aq_sec_var_act_species) = &
                            w%conc_hat(n_aq_prim+1:sa%num_aq_var_act_species)
                    end if
                    conc_comp = matmul(sa%comp_mat,conc_nc)
                    allocate(ref_names(n_prim_w))
                    call prim_species_names(w, n_prim_w, ref_names)
                    allocate(vas_names(n_vas))
                    call var_act_species_names(w, n_vas, vas_names)
                    write(unit,"(4x,'water ',I0)") idx
                    allocate(comb_lines(n_prim_w))
                    maxlen = 1
                    do jr=1,n_prim_w
                        comb_lines(jr) = comb_expr(sa%comp_mat(jr,:), vas_names)
                        maxlen = max(maxlen, len_trim(comb_lines(jr)))
                    end do
                    do jr=1,n_prim_w
                        write(unit,"(6x,A24,' = ',A,' = ',ES26.16)") &
                            trim(ref_names(jr)), comb_lines(jr)(1:maxlen), conc_comp(jr)
                    end do
                    deallocate(comb_lines, vas_names, ref_names, conc_comp, conc_nc)
                    end associate
                end do
            end do

            !> --- Surface-complex (exchange) concentrations ---
            !> Adsorbed exchange species are the LAST num_surf_compl entries of each water's
            !> solid_chemistry%concentrations (after the minerals), surf_compl(1) being the free
            !> site x-; the exchanger zone is shared, so one chem_syst header applies to all rows.
            if (has_conc .and. chem%chem_syst%cat_exch_zone%num_surf_compl>0 .and. allocated(chem%target_solids)) then
                nsc = chem%chem_syst%cat_exch_zone%num_surf_compl
                write(unit,"(2x,'conc surf compl')")
                write(unit,"(A18)", advance='no') 'Target solid'
                do j=1,nsc
                    write(unit,"(A26)", advance='no') trim(chem%chem_syst%cat_exch_zone%surf_compl(j)%name)
                end do
                write(unit,*)
                do ii=1,chem%num_target_solids
                    if (.not. associated(chem%target_solids(ii)%reactive_zone)) cycle
                    if (chem%target_solids(ii)%reactive_zone%cat_exch_zone%num_surf_compl<=0) cycle
                    associate(sc => chem%target_solids(ii))
                    off = sc%num_solids - sc%reactive_zone%cat_exch_zone%num_surf_compl
                    write(unit,"(I18)", advance='no') ii
                    do j=1,nsc
                        val = 0d0
                        if (off+j>=1 .and. off+j<=size(sc%concentrations)) val = sc%concentrations(off+j)
                        write(unit,"(ES26.16)", advance='no') val
                    end do
                    write(unit,*)
                    end associate
                end do
            end if

            !> --- Reaction rates (end-of-step instantaneous) ---
            !> Mineral columns are indexed by the GLOBAL chem_syst mineral list (chem%chem_syst%minerals),
            !> NOT by the reference water's local (per-zone) mineral order: different reactive/mineral
            !> zones can have DIFFERENT minerals (e.g. calcite zone vs anhydrite zone), so a fixed local
            !> offset silently mislabels/mixes rates across zones. Each water maps its own local eq/kin
            !> mineral rate into the correct global column via ind_mins_chem_syst/ind_min_chem_syst.
            if (has_rate .and. associated(chem%waters(ridx)%solid_chemistry)) then
                associate(rw => chem%waters(ridx), rz => chem%waters(ridx)%solid_chemistry%reactive_zone, &
                          cs => chem%waters(ridx)%solid_chemistry%reactive_zone%chem_syst)
                n_v = size(rw%ind_var_act_species)
                n_aq_prim = rz%speciation_alg%num_aq_prim_species
                nre=0; if (allocated(rw%rxn_rate_eq)) nre = size(rw%rxn_rate_eq)
                nm = chem%chem_syst%num_minerals   !< GLOBAL mineral count, eq AND kin (find_mineral_col indexes this list)
                n_lin = cs%num_lin_kin_reacts
                n_redox = cs%num_redox_kin_reacts
                nrct = nre+nm+n_lin+n_redox
                if (nrct>0) then
                    allocate(rct_names(nrct))
                    do j=1,nre
                        if (n_aq_prim+j<=n_v) then
                            rct_names(j) = rw%aq_phase%aq_species(rw%indices_aq_phase(rw%ind_var_act_species(n_aq_prim+j)))%name
                        else
                            rct_names(j) = cs%eq_reacts(j)%name
                        end if
                    end do
                    do j=1,nm
                        rct_names(nre+j) = strip_trailing2(chem%chem_syst%minerals(j)%name)
                    end do
                    do j=1,n_lin
                        rct_names(nre+nm+j) = cs%lin_kin_reacts(j)%name
                    end do
                    do j=1,n_redox
                        rct_names(nre+nm+n_lin+j) = cs%redox_kin_reacts(j)%name
                    end do
                end if
                end associate
                if (nrct>0) then
                    write(unit,"(2x,'react rate')")
                    write(unit,"(A18)", advance='no') 'Target solid'
                    do j=1,nrct
                        write(unit,"(A26)", advance='no') trim(rct_names(j))
                    end do
                    write(unit,*)
                    allocate(rowvals(nrct))
                    do ii=1,chem%chem_out_options%num_waters
                        idx = chem%chem_out_options%ind_waters(ii)
                        if (idx > chem%num_target_waters) cycle  !< react_rate only for domain target waters
                        rowvals=0d0
                        associate(w => chem%waters(idx))
                        if (nre>0 .and. allocated(w%rxn_rate_eq)) then
                            m=min(nre,size(w%rxn_rate_eq)); rowvals(1:m)=w%rxn_rate_eq(1:m)
                        end if
                        if (associated(w%solid_chemistry)) then
                            associate(wz => w%solid_chemistry%reactive_zone)
                            if (allocated(w%solid_chemistry%rxn_rate_eq)) then
                                do j=1,wz%num_minerals
                                    if (j<=size(w%solid_chemistry%rxn_rate_eq)) then
                                        g = find_mineral_col(wz%chem_syst%minerals(wz%ind_mins_chem_syst(j))%name)
                                        if (g>=1) rowvals(nre+g)=rowvals(nre+g)+w%solid_chemistry%rxn_rate_eq(j)
                                    end if
                                end do
                            end if
                            end associate
                            if (associated(w%solid_chemistry%mineral_zone) .and. allocated(w%solid_chemistry%rxn_rate_kin)) then
                                associate(wmz => w%solid_chemistry%mineral_zone)
                                do j=1,wmz%num_minerals_kin
                                    if (j<=size(w%solid_chemistry%rxn_rate_kin)) then
                                        g = find_mineral_col(wmz%chem_syst%minerals(wmz%ind_min_chem_syst(j))%name)
                                        if (g>=1) rowvals(nre+g)=rowvals(nre+g)+w%solid_chemistry%rxn_rate_kin(j)
                                    end if
                                end do
                                end associate
                            end if
                        end if
                        if ((n_lin+n_redox)>0 .and. allocated(w%rxn_rate_kin)) then
                            m=min(n_lin+n_redox,size(w%rxn_rate_kin)); rowvals(nre+nm+1:nre+nm+m)=w%rxn_rate_kin(1:m)
                        end if
                        end associate
                        write(unit,"(I18)", advance='no') idx
                        do j=1,nrct
                            write(unit,"(ES26.16)", advance='no') rowvals(j)
                        end do
                        write(unit,*)
                    end do
                    deallocate(rowvals)
                    deallocate(rct_names)
                end if
            end if

            !> --- Absolute cumulative reaction amounts (integrated over time) ---
            !> Mineral columns are indexed by the GLOBAL chem_syst mineral list (chem%chem_syst%minerals), 
            !! NOT by the reference water's local (per-zone) mineral order: different reactive/mineral
            !! zones can have DIFFERENT minerals (e.g. calcite zone vs anhydrite zone), so a fixed local
            !! offset silently mislabels/mixes amounts across zones. Each water maps its own local eq/kin mineral amount 
            !! into the correct global column via ind_mins_chem_syst/ind_min_chem_syst.
            !! Units: the mineral columns are written in kg (reaction amount [mol] times the mineral's molecular weight
            !! [kg/mol], an optional trailing value of its database line; 0 kg when the database gives none), the aqueous
            !! equilibrium and the linear/redox kinetic columns in mol (master25 carries no molecular weight for the
            !! aqueous complexes); each column header carries its unit.
            if (has_rxn_amnt .and. associated(chem%waters(ridx)%solid_chemistry)) then
                associate(rw => chem%waters(ridx), rz => chem%waters(ridx)%solid_chemistry%reactive_zone, &
                          cs => chem%waters(ridx)%solid_chemistry%reactive_zone%chem_syst)
                n_v = size(rw%ind_var_act_species)
                n_aq_prim = rz%speciation_alg%num_aq_prim_species
                nre=0; if (allocated(rw%rxn_rate_eq)) nre = size(rw%rxn_rate_eq)
                nm = chem%chem_syst%num_minerals   !< GLOBAL mineral count, eq AND kin (find_mineral_col indexes this list)
                n_lin = cs%num_lin_kin_reacts
                n_redox = cs%num_redox_kin_reacts
                nrct = nre+nm+n_lin+n_redox
                if (nrct>0) then
                    allocate(rct_names(nrct),rct_fact(nrct))
                    rct_fact=1d0 !< mol
                    do j=1,nre
                        if (n_aq_prim+j<=n_v) then
                            rct_names(j) = trim(rw%aq_phase%aq_species(rw%indices_aq_phase(rw%ind_var_act_species(n_aq_prim+j)))%name)//&
                                ' [mol]'
                        else
                            rct_names(j) = trim(cs%eq_reacts(j)%name)//' [mol]'
                        end if
                    end do
                    do j=1,nm
                        rct_names(nre+j) = trim(strip_trailing2(chem%chem_syst%minerals(j)%name))//' [kg]'
                        rct_fact(nre+j) = chem%chem_syst%minerals(j)%mineral%molecular_weight !< kg/mol (0 if not in the database)
                    end do
                    do j=1,n_lin
                        rct_names(nre+nm+j) = trim(cs%lin_kin_reacts(j)%name)//' [mol]'
                    end do
                    do j=1,n_redox
                        rct_names(nre+nm+n_lin+j) = trim(cs%redox_kin_reacts(j)%name)//' [mol]'
                    end do
                end if
                end associate
                if (nrct>0) then
                    write(unit,"(2x,'abs cum react amount')")
                    write(unit,"(A18)", advance='no') 'Target solid'
                    do j=1,nrct
                        write(unit,"(A26)", advance='no') trim(rct_names(j))
                    end do
                    write(unit,*)
                    allocate(rowvals(nrct))
                    do ii=1,chem%chem_out_options%num_waters
                        idx = chem%chem_out_options%ind_waters(ii)
                        if (idx > chem%num_target_waters) cycle  !< react_rate only for domain target waters
                        rowvals=0d0
                        associate(w => chem%waters(idx))
                        if (nre>0 .and. allocated(w%rxn_rate_eq)) then
                            m=min(nre,size(w%rxn_rate_eq)); rowvals(1:m)=abs(w%cum_rxn_amnt_eq(1:m))
                        end if
                        if (associated(w%solid_chemistry)) then
                            associate(wz => w%solid_chemistry%reactive_zone)
                            if (allocated(w%solid_chemistry%rxn_rate_eq)) then
                                do j=1,wz%num_minerals
                                    if (j<=size(w%solid_chemistry%rxn_rate_eq)) then
                                        g = find_mineral_col(wz%chem_syst%minerals(wz%ind_mins_chem_syst(j))%name)
                                        if (g>=1) rowvals(nre+g)=rowvals(nre+g)+abs(w%solid_chemistry%cum_rxn_amnt_eq(j))
                                    end if
                                end do
                            end if
                            end associate
                            if (associated(w%solid_chemistry%mineral_zone) .and. allocated(w%solid_chemistry%rxn_rate_kin)) then
                                associate(wmz => w%solid_chemistry%mineral_zone)
                                do j=1,wmz%num_minerals_kin
                                    if (j<=size(w%solid_chemistry%rxn_rate_kin)) then
                                        g = find_mineral_col(wmz%chem_syst%minerals(wmz%ind_min_chem_syst(j))%name)
                                        if (g>=1) rowvals(nre+g)=rowvals(nre+g)+abs(w%solid_chemistry%cum_rxn_amnt_kin(j))
                                    end if
                                end do
                                end associate
                            end if
                        end if
                        if ((n_lin+n_redox)>0 .and. allocated(w%rxn_rate_kin)) then
                            m=min(n_lin+n_redox,size(w%rxn_rate_kin)); rowvals(nre+nm+1:nre+nm+m)=abs(w%cum_rxn_amnt_kin(1:m))
                        end if
                        end associate
                        write(unit,"(I18)", advance='no') idx
                        do j=1,nrct
                            write(unit,"(ES26.16)", advance='no') rowvals(j)*rct_fact(j) !< mineral columns: mol -> kg
                        end do
                        write(unit,*)
                    end do
                    deallocate(rowvals)
                    deallocate(rct_names,rct_fact)
                end if
            end if

            !> --- Mineral-to-pore volume ratios V_min/V_pore (>= 0, may exceed 1) ---
            !> Columns are the GLOBAL chem_syst mineral list, same rationale as react_rate above.
            if (has_vf .and. chem%chem_syst%num_minerals>0) then
                nvf = chem%chem_syst%num_minerals
                allocate(vf_names(nvf))
                do j=1,nvf
                    vf_names(j) = strip_trailing2(chem%chem_syst%minerals(j)%name)
                end do
                write(unit,"(2x,'min-to-pore vol ratio')")
                write(unit,"(A18)", advance='no') 'Target solid'
                do j=1,nvf
                    write(unit,"(A26)", advance='no') trim(vf_names(j))
                end do
                write(unit,*)
                allocate(rowvals(nvf))
                do ii=1,chem%chem_out_options%num_waters
                    idx = chem%chem_out_options%ind_waters(ii)
                    if (idx > chem%num_target_waters) cycle  !< min-to-pore only for domain target waters
                    rowvals=0d0
                    associate(w => chem%waters(idx))
                    if (associated(w%solid_chemistry)) then
                        if (allocated(w%solid_chemistry%min_pore_vol_ratios) .and. associated(w%solid_chemistry%mineral_zone)) then
                            associate(wmz => w%solid_chemistry%mineral_zone)
                            nk = wmz%num_minerals_kin
                            do j=1,nk
                                if (j<=size(w%solid_chemistry%min_pore_vol_ratios)) then
                                    g = find_mineral_col(wmz%chem_syst%minerals(wmz%ind_min_chem_syst(j))%name)
                                    if (g>=1) rowvals(g) = w%solid_chemistry%min_pore_vol_ratios(j)
                                end if
                            end do
                            end associate
                            associate(wz => w%solid_chemistry%reactive_zone)
                            do j=1,wz%num_minerals
                                if (nk+j<=size(w%solid_chemistry%min_pore_vol_ratios)) then
                                    g = find_mineral_col(wz%chem_syst%minerals(wz%ind_mins_chem_syst(j))%name)
                                    if (g>=1) rowvals(g) = w%solid_chemistry%min_pore_vol_ratios(nk+j)
                                end if
                            end do
                            end associate
                        end if
                    end if
                    end associate
                    write(unit,"(I18)", advance='no') idx
                    do j=1,nvf
                        write(unit,"(ES26.16)", advance='no') rowvals(j)
                    end do
                    write(unit,*)
                end do
                deallocate(rowvals, vf_names)
            end if

            !> Within-zone sub-mixtures of this step go to their own sibling file (see write_zone_waters)
            if (allocated(chem%zone_waters) .and. time>0d0) call write_zone_waters()
        contains
            !> Resolves a mineral name to its column in the GLOBAL chem%chem_syst mineral list.
            !> NEEDED because ind_mins_chem_syst/ind_min_chem_syst index each water's OWN
            !> reactive_zone/mineral_zone chem_syst pointer, which need not share the same
            !> mineral numbering as chem%chem_syst (name is the only reliable common key).
            function find_mineral_col(mname) result(g)
                character(len=*), intent(in) :: mname
                integer(kind=4) :: g, k
                g = 0
                do k=1,chem%chem_syst%num_minerals
                    if (trim(chem%chem_syst%minerals(k)%name)==trim(mname)) then
                        g = k
                        exit
                    end if
                end do
            end function find_mineral_col

            !> Fills names(1:np) with the NAMES of the np PRIMARY species (components) of water w,
            !> matching the ordering of conc_nc / comp_mat: aqueous primary species first
            !> (1..num_aq_prim_species), then the solid/gas primary species. The only non-aqueous
            !> primary species in these models is the cation-exchange free site x- (surf_compl(1)),
            !> placed last (position num_prim_species) as in get_conc_var_act; any other slot falls back
            !> to a generic 'comp<j>' label.
            subroutine prim_species_names(w, np, names)
                class(aqueous_chemistry_c), intent(in) :: w
                integer(kind=4), intent(in) :: np
                character(len=*), intent(out) :: names(:)
                integer(kind=4) :: jj, naqp
                naqp = w%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                do jj=1,np
                    if (jj<=naqp) then
                        names(jj) = w%aq_phase%aq_species(w%indices_aq_phase(w%ind_var_act_species(jj)))%name
                    else if (w%solid_chemistry%reactive_zone%speciation_alg%flag_cat_exch .and. jj==np) then
                        names(jj) = w%solid_chemistry%reactive_zone%cat_exch_zone%surf_compl(1)%name
                    else
                        write(names(jj),'(A,I0)') 'comp', jj
                    end if
                end do
            end subroutine prim_species_names

            !> Fills names(1:nv) with the NAMES of the nv VARIABLE-ACTIVITY species of water w,
            !> matching the column ordering of comp_mat / conc_nc (see get_conc_var_act):
            !>   1..num_aq_prim_species        : aqueous primary species
            !>   num_prim_species              : cation-exchange free site x- (if flag_cat_exch)
            !>   num_prim_species+1 .. +Nsec   : aqueous secondary variable-activity species
            !>   then                          : surface complexes (adsorbed exchange cations)
            !> Any slot that cannot be identified falls back to a generic 'sp<j>' label.
            subroutine var_act_species_names(w, nv, names)
                class(aqueous_chemistry_c), intent(in) :: w
                integer(kind=4), intent(in) :: nv
                character(len=*), intent(out) :: names(:)
                integer(kind=4) :: jj, naqp, naqsec, nprim, nexch
                associate(sa => w%solid_chemistry%reactive_zone%speciation_alg)
                naqp   = sa%num_aq_prim_species
                naqsec = sa%num_aq_sec_var_act_species
                nprim  = sa%num_prim_species
                do jj=1,nv
                    write(names(jj),'(A,I0)') 'sp', jj   !< default fallback label
                end do
                !> aqueous primary species (positions 1..naqp)
                do jj=1,min(naqp,nv)
                    names(jj) = w%aq_phase%aq_species(w%indices_aq_phase(w%ind_var_act_species(jj)))%name
                end do
                !> aqueous secondary variable-activity species (positions nprim+1 .. nprim+naqsec)
                do jj=1,naqsec
                    if (nprim+jj<=nv) then
                        names(nprim+jj) = &
                            w%aq_phase%aq_species(w%indices_aq_phase(w%ind_var_act_species(naqp+jj)))%name
                    end if
                end do
                !> cation-exchange species (free site + adsorbed cations), if present
                if (sa%flag_cat_exch) then
                    if (nprim>=1 .and. nprim<=nv) then
                        names(nprim) = w%solid_chemistry%reactive_zone%cat_exch_zone%surf_compl(1)%name
                    end if
                    nexch = w%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
                    do jj=1,nexch
                        if (nprim+naqsec+jj<=nv) then
                            names(nprim+naqsec+jj) = &
                                w%solid_chemistry%reactive_zone%cat_exch_zone%surf_compl(1+jj)%name
                        end if
                    end do
                end if
                end associate
            end subroutine var_act_species_names

            !> Removes a trailing '2' from a mineral display name (e.g. 'calcite2' -> 'calcite').
            function strip_trailing2(name) result(disp)
                character(len=*), intent(in) :: name
                character(len=256) :: disp
                integer(kind=4) :: L
                disp = name
                L = len_trim(disp)
                if (L>=1) then
                    if (disp(L:L)=='2') disp(L:L) = ' '
                end if
            end function strip_trailing2

            !> Linear combination string "a + 2.0 b - c" from one component-matrix row (zero coefficients omitted).
            function comb_expr(coeffs, names) result(expr)
                real(kind=8), intent(in) :: coeffs(:)
                character(len=*), intent(in) :: names(:)
                character(len=1024) :: expr
                character(len=64) :: term
                integer(kind=4) :: jj, nt
                expr = ''
                nt = 0
                do jj=1,size(coeffs)
                    if (abs(coeffs(jj)) < 1d-12) cycle
                    if (abs(abs(coeffs(jj))-1d0) < 1d-12) then
                        term = trim(names(jj))
                    else
                        write(term,'(F0.1,1x,A)') abs(coeffs(jj)), trim(names(jj))
                    end if
                    if (nt==0) then
                        if (coeffs(jj) < 0d0) then
                            expr = '- '//trim(adjustl(term))
                        else
                            expr = trim(adjustl(term))
                        end if
                    else if (coeffs(jj) < 0d0) then
                        expr = trim(expr)//' - '//trim(adjustl(term))
                    else
                        expr = trim(expr)//' + '//trim(adjustl(term))
                    end if
                    nt = nt + 1
                end do
            end function comb_expr

            !> Writes the current step's within-zone sub-mixtures, chem%zone_waters(i,z) (one zone water per
            !> domain target water i and reactive zone z, filled by cross_zone_mix_core, aqueous_chemistry_m),
            !> to the sibling file "<root><tags>_zone_waters.out", found/opened through INQUIRE on the main
            !> output unit like the component mass-balance file. A zone water belongs to the current step iff
            !> lambda_c>0: cross_zone_mix_core resets lambda_c on every zone water of the target before mixing
            !> and sets it only on the zones it builds. The zone's speciation algebra and the sub-mixture's
            !> mineral/exchange equilibrium rates are reached through the zone water's solid chemistry
            !> (chem%zone_solids(i,z), linked once by associate_zone_waters_zone_solids). Same variable blocks
            !> as the main file where they exist for a temporary water (conc, pH, conc_hat, react_rate:
            !> equilibrium only); min-to-pore and cumulative amounts do not. Only records of requested output
            !> waters are written, in output-water order then by zone id.
            subroutine write_zone_waters()
                integer(kind=4) :: r, i, z, jj, jr2, zunit, nzw, nzz, nsel
                integer(kind=4) :: n_ref_z, n_prim_z, n_aq_prim_z, n_aq_sec_z, n_vas_z, maxlen_z
                logical :: zopen
                integer(kind=4), allocatable :: sel_i(:), sel_z(:) !< (row, zone) of the selected zone waters
                character(len=512) :: fname
                character(len=:), allocatable :: zfname
                character(len=256), allocatable :: znames(:), zprim_names(:), zvas_names(:)
                character(len=1024), allocatable :: zcomb(:)
                real(kind=8), allocatable :: zconc_nc(:), zconc_comp(:)
                real(kind=8) :: zval
                logical :: hat_ok !< conc_hat of the zone water is laid out like its ind_var_act_species ([aq prim | aq sec])

                !> current-step zone waters (lambda_c>0) of the requested output waters, matched by the global
                !> target-water index stored on each zone water (all built zone waters of a row share it)
                nzw = size(chem%zone_waters,1)
                nzz = size(chem%zone_waters,2)
                allocate(sel_i(nzw*nzz), sel_z(nzw*nzz))
                nsel = 0
                do jj=1,chem%chem_out_options%num_waters
                    do i=1,nzw
                        if (.not. any(chem%zone_waters(i,:)%lambda_c>0d0 .and. &
                                      chem%zone_waters(i,:)%tar_wat_idx==chem%chem_out_options%ind_waters(jj))) cycle
                        do z=1,nzz
                            if (chem%zone_waters(i,z)%lambda_c<=0d0) cycle
                            if (.not. associated(chem%zone_waters(i,z)%aq_chem%solid_chemistry)) cycle
                            nsel = nsel+1
                            sel_i(nsel) = i
                            sel_z(nsel) = z
                        end do
                    end do
                end do
                if (nsel==0) then
                    deallocate(sel_i, sel_z)
                    return
                end if

                inquire(unit=unit,name=fname)
                zfname = trim(fname)
                if (len(zfname)>4) then
                    if (zfname(len(zfname)-3:)=='.out') zfname = zfname(1:len(zfname)-4)
                end if
                if (len(zfname)>7) then
                    if (zfname(len(zfname)-6:)=='_output') zfname = zfname(1:len(zfname)-7)
                end if
                zfname = zfname//'_zone_waters.out'
                inquire(file=zfname,opened=zopen,number=zunit)
                if (.not. zopen) then
                    open(newunit=zunit,file=zfname,status='replace',action='write')
                    write(zunit,"(2x,'REMIX - within-zone sub-mixtures of the cross-zone mixing (compute_c_mix_by_zones)')")
                    write(zunit,"(2x,'One temporary water per (target water, reactive zone of its mixing stencil) and output time step;')")
                    write(zunit,"(2x,'only target waters whose mixing stencil spans several reactive zones produce them.')")
                    write(zunit,"(2x,'Row labels: Water index = target water being mixed, Reactive zone = zone id of the sub-mixture,')")
                    write(zunit,"(2x,'            lambda_c = total concentration mixing ratio of that zone (its cross-zone weight)')")
                    write(zunit,"(2x,'  conc_eq       : speciated sub-mixture c_eq_z (aqueous variable-activity species),')")
                    write(zunit,"(2x,'                  then its component concentrations U_z*c_eq_z')")
                    write(zunit,"(2x,'  pH            : pH of the speciated sub-mixture')")
                    write(zunit,"(2x,'  conc_hat      : within-zone mix c_hat_z entering the speciation (conservative mix c_mix_z')")
                    write(zunit,"(2x,'                  plus the limited previous-step kinetic term), then its component concentrations U_z*c_hat_z')")
                    write(zunit,"(2x,'  react_rate    : equilibrium reaction rates of the sub-mixture (aqueous, mineral, exchange);')")
                    write(zunit,"(2x,'                  kinetic rates, min-to-pore and cumulative amounts are not defined for a temporary water')")
                end if
                write(zunit,"(/,2x,'t = ',ES26.16,/)") time

                !> reference-water species names define the species columns (as in the main file)
                associate(rw => chem%waters(ridx))
                allocate(znames(size(rw%ind_var_act_species)))
                n_ref_z = 0
                do jj=1,size(rw%ind_var_act_species)
                    if (.not. rw%aq_phase%aq_species(rw%indices_aq_phase(rw%ind_var_act_species(jj)))%cst_act_flag) then
                        n_ref_z = n_ref_z+1
                        znames(n_ref_z) = rw%aq_phase%aq_species(rw%indices_aq_phase(rw%ind_var_act_species(jj)))%name
                    end if
                end do
                end associate

                if (has_conc) then
                    write(zunit,"(2x,'conc eq aq var act sp')")
                    write(zunit,"(A18,A14,A26)", advance='no') 'Water index','Reactive zone','lambda_c'
                    do jr2=1,n_ref_z
                        write(zunit,"(A26)", advance='no') trim(znames(jr2))
                    end do
                    write(zunit,*)
                    do r=1,nsel
                        associate(zr => chem%zone_waters(sel_i(r),sel_z(r)), w => chem%zone_waters(sel_i(r),sel_z(r))%aq_chem)
                        write(zunit,"(I18,I14,ES26.16)", advance='no') zr%tar_wat_idx, zr%zone_id, zr%lambda_c
                        do jr2=1,n_ref_z
                            zval = 0d0
                            do jj=1,size(w%ind_var_act_species)
                                if (trim(w%aq_phase%aq_species(w%indices_aq_phase(w%ind_var_act_species(jj)))%name)== &
                                    trim(znames(jr2))) then
                                    zval = w%concentrations(w%ind_var_act_species(jj))
                                    exit
                                end if
                            end do
                            write(zunit,"(ES26.16)", advance='no') zval
                        end do
                        write(zunit,*)
                        end associate
                    end do
                end if

                if (has_pH) then
                    write(zunit,"(2x,'pH')")
                    write(zunit,"(A18,A14,A26,A26)") 'Water index','Reactive zone','lambda_c','pH'
                    do r=1,nsel
                        associate(zr => chem%zone_waters(sel_i(r),sel_z(r)))
                        write(zunit,"(I18,I14,ES26.16,ES26.16)") zr%tar_wat_idx, zr%zone_id, zr%lambda_c, zr%aq_chem%pH
                        end associate
                    end do
                end if

                if (has_conc_hat) then
                    write(zunit,"(2x,'conc hat aq var act sp')")
                    write(zunit,"(A18,A14,A26)", advance='no') 'Water index','Reactive zone','lambda_c'
                    do jr2=1,n_ref_z
                        write(zunit,"(A26)", advance='no') trim(znames(jr2))
                    end do
                    write(zunit,*)
                    do r=1,nsel
                        associate(zr => chem%zone_waters(sel_i(r),sel_z(r)), w => chem%zone_waters(sel_i(r),sel_z(r))%aq_chem)
                        write(zunit,"(I18,I14,ES26.16)", advance='no') zr%tar_wat_idx, zr%zone_id, zr%lambda_c
                        !> conc_hat and ind_var_act_species share the [aq prim | aq sec var act] layout of the zone's
                        !> aqueous variable-activity species (both sized num_aq_var_act_species), so position jj of
                        !> one names the species of position jj of the other. Anything else is a bug upstream
                        !> (mix_within_zone_core writes conc_hat in exactly this layout): refuse to print a
                        !> misaligned row instead of silently writing zeros.
                        hat_ok=.false.
                        if (allocated(w%conc_hat)) hat_ok=(size(w%conc_hat)==size(w%ind_var_act_species))
                        if (allocated(w%conc_hat) .and. .not. hat_ok) then
                            print '(A,I0,A,I0,A,I0,A,I0)', "[ZONE-WATERS] ERROR: target water ", zr%tar_wat_idx, &
                                " zone ", zr%zone_id, ": size(conc_hat)=", size(w%conc_hat), &
                                " /= size(ind_var_act_species)=", size(w%ind_var_act_species)
                            error stop "[ZONE-WATERS] zone water conc_hat not laid out like its aqueous variable-activity species"
                        end if
                        do jr2=1,n_ref_z
                            zval = 0d0
                            if (hat_ok) then
                                do jj=1,size(w%ind_var_act_species)
                                    if (trim(w%aq_phase%aq_species(w%indices_aq_phase(w%ind_var_act_species(jj)))%name)== &
                                        trim(znames(jr2))) then
                                        zval = w%conc_hat(jj)
                                        exit
                                    end if
                                end do
                            end if
                            write(zunit,"(ES26.16)", advance='no') zval
                        end do
                        write(zunit,*)
                        end associate
                    end do
                end if

                !> Components in each sub-mixture's OWN zone basis (as the per-water conc comp block). The
                !> conc_nc layout is rebuilt from the aqueous species only: a sub-mixture carries no
                !> exchanger state, so any non-aqueous variable-activity slot stays 0.
                if (has_conc) then
                    write(zunit,"(2x,'conc comp')")
                    do r=1,nsel
                        associate(zr => chem%zone_waters(sel_i(r),sel_z(r)), w => chem%zone_waters(sel_i(r),sel_z(r))%aq_chem, &
                                  sa => chem%zone_waters(sel_i(r),sel_z(r))%aq_chem%solid_chemistry%reactive_zone%speciation_alg)
                        n_aq_prim_z = sa%num_aq_prim_species
                        n_aq_sec_z  = sa%num_aq_sec_var_act_species
                        n_prim_z    = sa%num_prim_species
                        n_vas_z     = sa%num_var_act_species
                        allocate(zconc_nc(n_vas_z)); zconc_nc = 0d0
                        do jj=1,n_aq_prim_z
                            zconc_nc(jj) = w%concentrations(w%ind_var_act_species(jj))
                        end do
                        do jj=1,n_aq_sec_z
                            zconc_nc(n_prim_z+jj) = w%concentrations(w%ind_var_act_species(n_aq_prim_z+jj))
                        end do
                        zconc_comp = matmul(sa%comp_mat,zconc_nc)
                        allocate(zprim_names(n_prim_z), zvas_names(n_vas_z), zcomb(n_prim_z))
                        call prim_species_names(w, n_prim_z, zprim_names)
                        call var_act_species_names(w, n_vas_z, zvas_names)
                        write(zunit,"(4x,'water ',I0,'  zone ',I0,'  lambda_c ',ES14.6)") zr%tar_wat_idx, zr%zone_id, zr%lambda_c
                        maxlen_z = 1
                        do jr2=1,n_prim_z
                            zcomb(jr2) = comb_expr(sa%comp_mat(jr2,:), zvas_names)
                            maxlen_z = max(maxlen_z, len_trim(zcomb(jr2)))
                        end do
                        do jr2=1,n_prim_z
                            write(zunit,"(6x,A24,' = ',A,' = ',ES26.16)") &
                                trim(zprim_names(jr2)), zcomb(jr2)(1:maxlen_z), zconc_comp(jr2)
                        end do
                        deallocate(zconc_nc, zconc_comp, zprim_names, zvas_names, zcomb)
                        end associate
                    end do
                end if

                if (has_conc_hat) then
                    write(zunit,"(2x,'conc comp hat')")
                    do r=1,nsel
                        associate(zr => chem%zone_waters(sel_i(r),sel_z(r)), w => chem%zone_waters(sel_i(r),sel_z(r))%aq_chem, &
                                  sa => chem%zone_waters(sel_i(r),sel_z(r))%aq_chem%solid_chemistry%reactive_zone%speciation_alg)
                        if (.not. allocated(w%conc_hat)) cycle
                        if (size(w%conc_hat)/=sa%num_aq_var_act_species) cycle
                        n_aq_prim_z = sa%num_aq_prim_species
                        n_aq_sec_z  = sa%num_aq_sec_var_act_species
                        n_prim_z    = sa%num_prim_species
                        n_vas_z     = sa%num_var_act_species
                        allocate(zconc_nc(n_vas_z)); zconc_nc = 0d0
                        zconc_nc(1:n_aq_prim_z) = w%conc_hat(1:n_aq_prim_z)
                        zconc_nc(n_prim_z+1:n_prim_z+n_aq_sec_z) = w%conc_hat(n_aq_prim_z+1:n_aq_prim_z+n_aq_sec_z)
                        zconc_comp = matmul(sa%comp_mat,zconc_nc)
                        allocate(zprim_names(n_prim_z), zvas_names(n_vas_z), zcomb(n_prim_z))
                        call prim_species_names(w, n_prim_z, zprim_names)
                        call var_act_species_names(w, n_vas_z, zvas_names)
                        write(zunit,"(4x,'water ',I0,'  zone ',I0,'  lambda_c ',ES14.6)") zr%tar_wat_idx, zr%zone_id, zr%lambda_c
                        maxlen_z = 1
                        do jr2=1,n_prim_z
                            zcomb(jr2) = comb_expr(sa%comp_mat(jr2,:), zvas_names)
                            maxlen_z = max(maxlen_z, len_trim(zcomb(jr2)))
                        end do
                        do jr2=1,n_prim_z
                            write(zunit,"(6x,A24,' = ',A,' = ',ES26.16)") &
                                trim(zprim_names(jr2)), zcomb(jr2)(1:maxlen_z), zconc_comp(jr2)
                        end do
                        deallocate(zconc_nc, zconc_comp, zprim_names, zvas_names, zcomb)
                        end associate
                    end do
                end if

                !> Equilibrium reaction rates per sub-mixture, labelled in its OWN zone: aqueous reactions by
                !> the secondary species they eliminate there, then the zone's minerals and exchange cations.
                if (has_rate) then
                    write(zunit,"(2x,'react rate')")
                    do r=1,nsel
                        associate(zr => chem%zone_waters(sel_i(r),sel_z(r)), w => chem%zone_waters(sel_i(r),sel_z(r))%aq_chem, &
                                  zs => chem%zone_waters(sel_i(r),sel_z(r))%aq_chem%solid_chemistry, &
                                  rz => chem%zone_waters(sel_i(r),sel_z(r))%aq_chem%solid_chemistry%reactive_zone)
                        write(zunit,"(4x,'water ',I0,'  zone ',I0,'  lambda_c ',ES14.6)") zr%tar_wat_idx, zr%zone_id, zr%lambda_c
                        n_aq_prim_z = rz%speciation_alg%num_aq_prim_species
                        if (allocated(w%rxn_rate_eq)) then
                            do jj=1,size(w%rxn_rate_eq)
                                if (n_aq_prim_z+jj<=size(w%ind_var_act_species)) then
                                    write(zunit,"(6x,A24,' = ',ES26.16)") &
                                        trim(w%aq_phase%aq_species(w%indices_aq_phase(w%ind_var_act_species(n_aq_prim_z+jj)))%name), &
                                        w%rxn_rate_eq(jj)
                                else if (allocated(rz%chem_syst%eq_reacts)) then
                                    if (jj<=size(rz%chem_syst%eq_reacts)) write(zunit,"(6x,A24,' = ',ES26.16)") &
                                        trim(rz%chem_syst%eq_reacts(jj)%name), w%rxn_rate_eq(jj)
                                end if
                            end do
                        end if
                        if (allocated(zs%rxn_rate_eq)) then
                            do jj=1,rz%num_minerals
                                if (jj<=size(zs%rxn_rate_eq)) write(zunit,"(6x,A24,' = ',ES26.16)") &
                                    trim(strip_trailing2(rz%chem_syst%minerals(rz%ind_mins_chem_syst(jj))%name)), zs%rxn_rate_eq(jj)
                            end do
                            do jj=1,rz%cat_exch_zone%num_exch_cats
                                if (rz%num_minerals+jj<=size(zs%rxn_rate_eq)) write(zunit,"(6x,A24,' = ',ES26.16)") &
                                    trim(rz%cat_exch_zone%surf_compl(1+jj)%name), zs%rxn_rate_eq(rz%num_minerals+jj)
                            end do
                        end if
                        end associate
                    end do
                end if
                flush(zunit)
                deallocate(znames, sel_i, sel_z)
            end subroutine write_zone_waters
        end subroutine write_out_step_vars

        !> @brief Check the surface-complex site-balance invariants at every target water.
        !> @details Under the Gaines-Thomas convention each surface-complex activity is an
        !>          equivalent fraction, so at every time step and every exchanging target
        !>          water TWO invariants must hold: Sum(activities) = 1 and the adsorbed
        !>          equivalents Sum(z_i*c_i) = CEC. Both are recomputed from the CURRENT solid
        !>          concentrations (the stored activities array is not resynced by the reactive
        !>          step) and enforced (error stop if either deviates beyond a relative
        !>          tolerance). A per-step [SURF-DBG] line reports the worst deviations.
        subroutine check_surf_compl_act_sum(chem, k, time)
            class(chemistry_c), intent(in) :: chem   !< Chemistry object (target waters + solids)
            integer(kind=4), intent(in) :: k         !< Time step index
            real(kind=8), intent(in) :: time         !< Simulation time
            integer(kind=4) :: i, nsc, nex, ns, off, ic, z, ncheck, worst_act, worst_cec
            real(kind=8) :: equiv_sum, act_sum, cec, dev_act, dev_cec, max_dev_act, max_dev_cec
            real(kind=8), parameter :: tol = 1d-6    !< relative tolerance for the site-balance invariants
            ncheck = 0; max_dev_act = 0d0; max_dev_cec = 0d0; worst_act = 0; worst_cec = 0
            do i=1,chem%num_target_waters
                if (.not. associated(chem%waters(i)%solid_chemistry)) cycle
                if (.not. associated(chem%waters(i)%solid_chemistry%reactive_zone)) cycle
                if (chem%waters(i)%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl <= 0) cycle
                associate(sc => chem%waters(i)%solid_chemistry)
                nsc = sc%reactive_zone%cat_exch_zone%num_surf_compl
                nex = sc%reactive_zone%cat_exch_zone%num_exch_cats
                ns  = sc%num_solids
                off = ns - nsc
                cec = sc%CEC
                !> adsorbed equivalents Sum(z_i c_i): free site x- (z=1) plus adsorbed cations
                !> (z=valence), recomputed from the CURRENT solid concentrations.
                equiv_sum = sc%concentrations(off+1)                        !< free site x- (valence 1)
                do ic=1,nex
                    z = sc%reactive_zone%chem_syst%aq_phase%aq_species( &
                        sc%reactive_zone%cat_exch_zone%exch_cat_indices(ic))%valence
                    equiv_sum = equiv_sum + z*sc%concentrations(off+1+ic)
                end do
                !> Gaines-Thomas activities are the equivalent fractions a_i = z_i c_i/CEC,
                !> so Sum(a_i) = equiv_sum/CEC. Invariants: Sum(a_i)=1 and equiv_sum=CEC.
                act_sum = equiv_sum
                if (cec>0d0) act_sum = equiv_sum/cec
                dev_act = abs(act_sum - 1d0)                                !< |Sum(activities) - 1|
                dev_cec = abs(equiv_sum - cec)/max(abs(cec),1d-30)         !< |Sum(equivalents) - CEC| (relative)
                if (dev_act>max_dev_act) then; max_dev_act = dev_act; worst_act = i; end if
                if (dev_cec>max_dev_cec) then; max_dev_cec = dev_cec; worst_cec = i; end if
                ncheck = ncheck + 1
                if (dev_act>tol .or. dev_cec>tol) then
                    write(*,'(A,I6,A,I6,A,ES13.5,A,ES13.5,A,ES13.5)') &
                        ' [SURF-DBG] VIOLATION step ', k, ' target water ', i, &
                        ' sum(act)=', act_sum, ' sum(equiv)=', equiv_sum, ' CEC=', cec
                    error stop "Surface-complex site balance violated: sum(activities)/=1 or sum(equivalents)/=CEC"
                end if
                end associate
            end do
            if (ncheck>0) then
                write(*,'(A,I6,A,ES13.5,A,I6,A,ES12.4,A,ES12.4)') &
                    ' [SURF-DBG] step ', k, ' t=', time, ' exch cells=', ncheck, &
                    ' max|sum(act)-1|=', max_dev_act, ' max|sum(equiv)-CEC|/CEC=', max_dev_cec
            end if
        end subroutine check_surf_compl_act_sum

        !> @brief Relink target waters to reactive zones if the zone set or membership changed
        !> @details Reactive-zone mixing step vii: rebuilds chem%wat_indices_zones from the
        !>          reactive-zone ids of the waters' previous-step state (solid_chemistry_old, via
        !>          set_wat_indices_zones) only when it differs from the stored grouping.
        !>          The Lagrangian solvers call it after each displacement so the grouping follows
        !>          the parcels (zone of the cell each parcel comes from). With a Neumann outflow
        !>          they pass exit_wat_indices, the outflow boundary waters holding the parcels that
        !>          left the last cells: their previous-step concentrations were equilibrated in
        !>          those cells' zones, so they are classified like the target waters (and must then
        !>          NOT be listed among the mixer's external waters).
        subroutine reset_wat_indices_zones(chem)
            class(chemistry_c), intent(inout) :: chem
            !integer(kind=4), intent(in), optional :: exit_wat_indices(:) !< global indices of the exited-parcel (Neumann outflow) boundary waters
            integer(kind=4) :: i
            integer(kind=4), allocatable :: twiz(:)
            logical :: changed
            changed=(chem%wat_indices_zones%num_cols/=chem%num_reactive_zones)
            ! if (.not. changed) then
            !     ! do i=1,chem%num_reactive_zones
            !     !     call chem%set_wat_indices_zones(i,twiz)
            !     !     if (size(twiz)/=chem%wat_indices_zones%cols(i)%dim) then
            !     !         changed=.true.
            !     !     else if (any(twiz/=chem%wat_indices_zones%cols(i)%col_1)) then
            !     !         changed=.true.
            !     !     end if
            !     !     if (changed) exit
            !     ! end do
            ! end if
            ! if (.not. changed) return
            ! call chem%wat_indices_zones%allocate_array(chem%num_reactive_zones)
            ! do i=1,chem%num_reactive_zones
            !     call chem%set_wat_indices_zones(i,twiz)
            !     call chem%wat_indices_zones%cols(i)%allocate_vector(size(twiz))
            !     chem%wat_indices_zones%cols(i)%col_1=twiz
            ! end do
        end subroutine reset_wat_indices_zones

        subroutine set_num_interfaces(this,num_int)
            class(chemistry_c), intent(inout) :: this
            integer(kind=4), intent(in) :: num_int
            if (num_int<0) then
                error stop "Number of interfaces cannot be negative"
            end if
            this%num_interfaces=num_int
        end subroutine set_num_interfaces

        subroutine allocate_zone_waters(this)
            class(chemistry_c), intent(inout) :: this
            if (allocated(this%zone_waters)) then
                deallocate(this%zone_waters)
            end if
            allocate(this%zone_waters(this%num_target_solids,this%num_reactive_zones))
        end subroutine allocate_zone_waters

        subroutine allocate_zone_solids(this)
            class(chemistry_c), intent(inout) :: this
            if (allocated(this%zone_solids)) then
                deallocate(this%zone_solids)
            end if
            allocate(this%zone_solids(this%num_target_solids,this%num_reactive_zones))
        end subroutine allocate_zone_solids

        subroutine associate_zone_waters_zone_solids(this)
            class(chemistry_c), intent(inout) :: this
            integer(kind=4) :: i !> target water index
            integer(kind=4) :: j !> reactive zone index
            if (.not. allocated(this%zone_waters)) then
                error stop "Zone waters are not allocated"
            else if (.not. allocated(this%zone_solids)) then
                error stop "Zone solids are not allocated"
            else
                do i=1,this%num_target_solids
                    !if (this%target_solids(i)%is_interface)
                    do j=1,this%num_reactive_zones
                        call this%zone_solids(i,j)%set_reactive_zone(this%reactive_zones(j))
                        call this%zone_solids(i,j)%set_mineral_zone(this%min_zone_dummy)
                        !call this%zone_solids(i,j)%initialise_solid_chemistry()
                        call this%zone_waters(i,j)%aq_chem%set_solid_chemistry(this%zone_solids(i,j))
                        call this%zone_waters(i,j)%aq_chem%set_solid_chemistry_old()
                    end do
                end do
            end if
        end subroutine associate_zone_waters_zone_solids

        !> @brief Allocates the zone waters and their scratch zone solids (one per domain target water and
        !>        reactive zone) and links them (zone_waters(i,j)%aq_chem%solid_chemistry => zone_solids(i,j)).
        !> @details Must be called once num_target_waters and the reactive_zones array are FINAL: after
        !>          read_chemistry (which may append reactive zones) and after set_target_waters (RT_m), which
        !>          may add the boundary waters to the target waters. The solvers index zone_waters(idx,:) by
        !>          the target water index, so an array sized earlier (num_target_waters still 0) is written out
        !>          of bounds. Idempotent: previous arrays are released.
        subroutine set_zone_waters(this)
            class(chemistry_c), intent(inout) :: this
            if (this%num_target_solids<=0) error stop "set_zone_waters: number of target solids not set yet"
            if (.not. allocated(this%reactive_zones)) error stop "set_zone_waters: reactive zones not allocated"
            if (size(this%reactive_zones)<this%num_reactive_zones) error stop "set_zone_waters: fewer reactive zones than num_reactive_zones"
            call this%allocate_zone_waters()
            call this%allocate_zone_solids()
            call this%associate_zone_waters_zone_solids()
        end subroutine set_zone_waters

end module chemistry_m !< End of chemistry_m module
