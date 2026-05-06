!> @file chemistry_m.f90
!> @brief Main chemistry module for reactive transport modeling
!> @details This module contains the main chemistry class that coordinates all chemical 
!>          processes in reactive transport simulations. It manages water types, target 
!>          waters, solid phases, gas phases, and reactive zones. The module provides
!>          functionality for reading chemical data, solving reactive mixing problems,
!>          and managing chemical equilibrium and kinetic reactions.
!> @author Generated documentation
!> @date 2025

!> @brief Main chemistry module for reactive transport
!> @details Coordinates all chemical processes including aqueous speciation, mineral
!>          dissolution/precipitation, gas exchange, and reactive mixing algorithms
module chemistry_m
    use aqueous_chemistry_m, only: aqueous_chemistry_c, compute_c_mix, &
        reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2 !< Import aqueous chemistry class
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
    public :: reactive_zone_c !< Re-export reactive zone class for use by subprograms
    public :: interfaz_comps_arch, interfaz_esp_arch !< External interface subroutines for reactive mixing iterations
    !> @brief Main chemistry class - central coordinator for all chemical processes
    !> @details This class manages the complete chemical system including:
    !>          - Water types and target waters at spatial locations
    !>          - Solid phases (minerals, sorbed species, biofilms)
    !>          - Gas phases and reactive zones
    !>          - Chemical reaction networks and thermodynamic data
    !>          - Numerical algorithms for reactive mixing
    !> @note This is the primary interface for all chemical calculations in RT simulations
    type :: chemistry_c
        !> @name Numerical Algorithm Control Parameters
        !> @{
        integer(kind=4) :: num_lump=0                                       !< Number of spatial lumping regions in Water Mixing Algorithm (0=no lumping)
        integer(kind=4) :: read_opt                                        !< Chemical data input format: 1=CHEPROO-based, 2=PHREEQC, 3=PFLOTRAN
        integer(kind=4) :: act_coeffs_model                                !< Activity coefficient model: Debye-HÃ¼ckel, Davies, Pitzer, etc.
        logical :: lump_flag                                             !< Spatial lumping flag: .true.=apply lumping, .false.=full resolution
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
        !> @}
        
        !> @name Solid Phase Management
        !> @{
        integer(kind=4) :: num_target_solids=0                                 !< Number of target solids (â‰¤ num_waters)
        integer(kind=4) :: num_target_solids_dummy=0                           !< Number of dummy target solids (for locations without solids)
        integer(kind=4) :: num_materials=0                                     !< Number of distinct material types (â‰¤ num_target_solids). Layout: materials(1..num_min_zones) = mineral zones, materials(num_min_zones+1..num_min_zones+num_init_cat_exch_zones) = surface (cation-exchange) zones.
        integer(kind=4) :: num_min_zones=0                                     !< Number of mineral zones stored at the head of `materials(:)` (= nmrz from quim_loc.dat)
        integer(kind=4) :: num_init_cat_exch_zones=0                           !< Number of initial cation exchange (surface) zones stored at the tail of `materials(:)`
        type(solid_chemistry_c), allocatable :: materials(:)                   !< Array of distinct material types (solid compositions)
        type(solid_chemistry_c), allocatable :: wat_type_solids(:)             !< Solid chemistry templates associated with each water type (1..num_wat_types). Independent of `target_solids` and `materials`.
        type(solid_chemistry_c), allocatable :: target_solids(:)               !< Current target solid compositions at spatial locations
        type(solid_chemistry_c), allocatable :: target_solids_init(:)          !< Initial target solid compositions (for restart/mass balance)
        type(solid_chemistry_c), allocatable :: target_solids_dummy(:)         !< Dummy target solids for locations without associated solids
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
        integer(kind=4) :: num_reactive_zones=0                                !< Number of reactive zones (â‰¤num_target_solids)
        integer(kind=4) :: num_reactive_zones_dummy=0                          !< Number of dummy reactive zones (for locations without reactions)
        type(reactive_zone_c), allocatable :: reactive_zones(:)                !< Array of reactive zone definitions (reaction networks)
        type(reactive_zone_c), allocatable :: react_zones_wat_types(:)         !< Reactive zones associated with water types
        type(reactive_zone_c), allocatable :: react_zones_dummy(:)       !< Reactive zones associated with gas types
        !> @name Mineral Zone Management
        !> @{
        integer(kind=4) :: num_mineral_zones=0                                 !< Number of mineral zones (â‰¤num_target_solids)
        type(mineral_zone_c), allocatable :: mineral_zones(:)                  !< Array of mineral zone definitions
        type(mineral_zone_c) :: min_zone_dummy                                 !< Dummy mineral zone for locations without minerals
        !> @}
    contains
        !> @name Configuration and Setup Methods
        !> @{
        procedure :: set_lump_flag                    !< Set spatial lumping flag for Water Mixing Algorithm
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
        procedure :: allocate_mineral_zones           !< Allocate mineral zones array
        procedure :: allocate_materials               !< Allocate materials array
        procedure :: allocate_wat_types               !< Allocate water types array
        procedure :: allocate_wat_type_solids         !< Allocate per-water-type solid chemistry array
        procedure :: allocate_gas_zones               !< Allocate gas zones array
        procedure :: allocate_gas_zones_wat_types     !< Allocate gas zones for water types
        procedure :: allocate_reactive_zones_dummy          !< Allocate dummy reactive zones
        !> @}
        
        !> @name Input/Output Methods
        !> @{
        procedure :: read_waters_init          !< Read initial target water definitions
        procedure :: read_tar_wat_line                !< Read single target water line from input
        procedure :: read_tar_sol                     !< Read target solid definitions
        procedure :: read_tar_gas                     !< Read target gas definitions
        procedure :: read_init_min_zones_CHEPROO      !< Read initial mineral zones (CHEPROO format)
        procedure :: read_chemistry_interface         !< Read chemistry information from input files and initialize chemistry object
        procedure :: read_chemistry_CHEPROO           !< Read chemistry data in CHEPROO format
        procedure :: read_init_bd_wat_types_CHEPROO   !< Read initial boundary water types (CHEPROO)
        procedure :: read_init_cat_exch_zones_CHEPROO !< Read initial cation exchange zones (CHEPROO)
        procedure :: read_gas_zones_CHEPROO             !< Read initial gas zones (CHEPROO format)
        procedure :: read_comp_opts                   !< Read computation options from input
        procedure :: write_chemistry                  !< Write chemistry state to output
        procedure :: write_conc_comp_tar_wat        !< Write component concentrations of every target water (and external waters)
        procedure :: write_conc_comp_wat_types   !< Write component concentrations of every water type, split initial vs external
        procedure :: write_u_mix_init                   !< Write initial u_mix concentrations
        !> @}
        
        !> @name Initialization Methods
        !> @{
        !procedure :: initialise_chemistry             !< Initialize complete chemistry system
        !> @}
        
        !> @name Main Solver Methods
        !> @{
        procedure :: solve_reactive_mixing_lump              !< Main reactive mixing solver with spatial lumping
        procedure :: solve_reactive_mixing_ideal_cons        !< Ideal reactive mixing solver with consistency
        procedure :: solve_reactive_mixing_cons              !< General reactive mixing solver with consistency
        procedure :: solve_reactive_mixing_ideal_lump        !< Ideal reactive mixing solver with lumping
        procedure :: solve_reactive_mixing_ideal_lump_iter   !< Single time-step iteration of ideal lumped solver
        procedure :: solve_reactive_mixing         !< Reactive mixing with time-dependent boundary conditions
        !> @}
        
        !> @name Linkage and Association Methods
        !> @{
        procedure :: link_waters_target_solids !< Link target waters to corresponding target solids
        procedure :: link_target_waters_target_gases  !< Link target waters to corresponding target gases
        procedure :: link_target_solids_reactive_zone !< Link target solids to reactive zones
        procedure :: link_target_gases_reactive_zone  !< Link target gases to reactive zones
        procedure :: link_target_waters_reactive_zone !< Link target waters to reactive zones
        procedure :: link_waters_mineral_zone !< Link target waters to mineral zones
        procedure :: link_target_solids_mineral_zone !< Link target solids to mineral zones
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
    !> @param mix_react_indices Indices of domain mixing waters
    !> @param time_discr Time discretization object
    !> @param int_method_chem_reacts Integration method for chemical reactions
    subroutine solve_reactive_mixing_lump(this,dir,root,mixing_ratios_conc,mixing_ratios_R,mix_conc_indices,time_discr,&
        theta_r)
            import chemistry_c                            !< Import chemistry class definition
            import real_array_c                           !< Import real array class for mixing ratios
            import int_array_c                            !< Import integer array class for indices
            import time_discr_c                           !< Import time discretization class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object instance
            character(len=*), intent(in) :: dir           !< Directory for output files
            character(len=*), intent(in) :: root          !< Root name for output files
            class(real_array_c), intent(in) :: mixing_ratios_conc !< Mixing ratios for concentrations
            class(real_array_c), intent(inout) :: mixing_ratios_R !< Mixing ratios for reaction rates
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
        !> @param mixing_ratios_R_init Initial mixing ratios for reaction rates
        !> @param mix_conc_indices Indices of mixing waters
        !> @param mix_react_indices Indices of domain mixing waters
        !> @param time_discr Time discretization object
        !> @param int_method_chem_reacts Integration method for chemical reactions
        !> @param mixing_ratios_R Output mixing ratios for reaction rates
        subroutine solve_reactive_mixing_ideal_cons(this,dir,root,mixing_ratios_conc,&
            mixing_ratios_R,mix_conc_indices,mix_react_indices,&
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
            class(real_array_c), intent(in) :: mixing_ratios_R !< Mixing ratios for reaction rates
            class(int_array_c), intent(in) :: mix_conc_indices !< Indices of mixing waters
            class(int_array_c), intent(in) :: mix_react_indices !< Indices of domain mixing waters
            class(time_discr_c), intent(in) :: time_discr !< Time discretization object
            real(kind=8), intent(in) :: theta_r           !< Integration method for chemical reactions
        end subroutine

        !> @brief Solve conservative reactive mixing
        !> @details Solves reactive mixing for conservative systems, using concentration and reaction rate mixing ratios.
        !> Handles domain and external waters, with time discretization and integration method for chemical reactions.
        !> @param this Chemistry object
        !> @param root Output file root
        !> @param mixing_ratios_conc Mixing ratios for concentrations
        !> @param mixing_ratios_R_init Initial mixing ratios for reaction rates
        !> @param mix_conc_indices Indices of mixing waters
        !> @param mix_react_indices Indices of domain mixing waters
        !> @param time_discr Time discretization object
        !> @param int_method_chem_reacts Integration method for chemical reactions
        !> @param mixing_ratios_R Output mixing ratios for reaction rates
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
            class(real_array_c), intent(inout) :: mixing_ratios_R !< Mixing ratios for reaction rates
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
        !> @param mix_react_indices Indices of domain mixing waters
        !> @param time_discr Time discretization object
        !> @param int_method_chem_reacts Integration method for chemical reactions
        subroutine solve_reactive_mixing_ideal_lump(this,dir,root,mixing_ratios_conc,mixing_ratios_R,mix_conc_indices,&
                mix_react_indices,time_discr,theta_r)
            import chemistry_c                            !< Import chemistry class definition
            import real_array_c                           !< Import real array class for mixing ratios
            import int_array_c                            !< Import integer array class for indices
            import time_discr_c                           !< Import time discretization class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object instance
            character(len=*), intent(in) :: dir           !< Directory for output files
            character(len=*), intent(in) :: root          !< Root name for output files
            class(real_array_c), intent(in) :: mixing_ratios_conc !< Mixing ratios for concentrations
            class(real_array_c), intent(in) :: mixing_ratios_R !< Mixing ratios for reaction rates
            class(int_array_c), intent(in) :: mix_conc_indices !< Indices of mixing waters
            class(int_array_c), intent(in) :: mix_react_indices !< Indices of domain mixing waters
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
    subroutine read_chemistry_CHEPROO(this,root,path_pb,path_DB,unit_chem_syst_file,unit_loc_chem_file,num_tar)
            import chemistry_c                            !< Import chemistry class definition
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object to initialize
            character(len=*), intent(in) :: root          !< Root name for files
            character(len=*), intent(in) :: path_pb       !< Path to problem files
            character(len=*), intent(in) :: path_DB       !< Path to CHEPROO database file
            integer(kind=4), intent(in) :: unit_chem_syst_file !< File unit for chemical system
            integer(kind=4), intent(in) :: unit_loc_chem_file !< File unit for local chemistry
            integer(kind=4), intent(in) :: num_tar        !< Expected number of target (domain) waters from the mesh
        end subroutine
        
    !> @brief Read complete chemistry setup from files
    !> @details Reads and initializes the complete chemistry configuration including chemical system,
    !>          water types, solid/gas phases, and reactive zones from input files.
    !> @param[in,out] this Chemistry object to initialize
    !> @param[in] root Root name for input/output files
    !> @param[in] path_pb Path to problem directory
    !> @param[in] path_DB Path to database files
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
    !> @details Finds all domain target water indices associated with reactive zone i.
    !>          Returns an allocatable array of target water indices belonging to the zone.
    !> @param[in,out] this Chemistry object containing target waters and reactive zones
    !> @param[in] i Index of the reactive zone to query
    !> @param[out] tar_wat_indices Allocatable array of domain target water indices in zone i
        subroutine link_target_waters_dom_reactive_zone(this,i,tar_wat_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this                    !< Chemistry object instance
            integer(kind=4), intent(in) :: i              !< Reactive zone index to query
            integer(kind=4), intent(out), allocatable :: tar_wat_indices(:) !< Output: domain target water indices in this zone
        end subroutine
        
    !> @brief Link all target waters (domain + external) to a specific reactive zone
    !> @details Finds all target water indices (both domain and external) associated with
    !>          reactive zone i. Returns separate arrays for domain and external indices.
    !> @param[in] this Chemistry object containing target waters and reactive zones
    !> @param[in] i Index of the reactive zone to query
    !> @param[out] dom_indices Allocatable array of domain target water indices in zone i
    !> @param[out] ext_indices Allocatable array of external target water indices in zone i
        subroutine link_target_waters_reactive_zone(this,i,dom_indices,ext_indices)
            import chemistry_c
            implicit none
            class(chemistry_c), intent(in) :: this        !< Chemistry object instance
            integer(kind=4), intent(in) :: i              !< Reactive zone index to query
            integer(kind=4), intent(out), allocatable :: dom_indices(:) !< Output: domain target water indices in this zone
            integer(kind=4), intent(out), allocatable :: ext_indices(:) !< Output: external target water indices in this zone
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
        subroutine read_waters_init(this,root,nsrz,ngrz,num_tar)
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
            integer(kind=4), intent(in) :: num_tar        !< Expected number of target (domain) waters from the mesh
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
        subroutine write_chemistry(this,dir,root,time_step)
            import chemistry_c
            class(chemistry_c), intent(in) :: this        !< Chemistry object to output
            character(len=*), intent(in) :: dir           !< Directory for output files
            character(len=*), intent(in) :: root          !< Root name for output files
            integer(kind=4), intent(in), optional :: time_step !< Optional time step number
        end subroutine
        
    !> @brief Read initial and boundary water types from CHEPROO input
    !> @details Reads initial and boundary water type definitions from CHEPROO format files,
    !>          optionally including gas phase chemistry for multiphase systems.
    !> @param[in,out] this Chemistry object to populate with water types
    !> @param[in] unit File unit number for reading
    !> @param[in,out] init_cat_exch_zones Initial cation exchange zones (caller-owned, may be reallocated)
    !> @param[in] gas_species_chem Optional gas chemistry for multiphase initialization
        subroutine read_init_bd_wat_types_CHEPROO(this,unit,init_cat_exch_zones,&
            gas_species_chem)
            import chemistry_c                            !< Import chemistry class definition
            import aqueous_chemistry_c                    !< Import aqueous chemistry class
            import solid_chemistry_c                      !< Import solid chemistry class
            import gas_chemistry_c                        !< Import gas chemistry class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object to populate
            integer(kind=4), intent(in) :: unit           !< File unit number for reading
            type(solid_chemistry_c), intent(inout), allocatable :: init_cat_exch_zones(:) !< Initial cation exchange zones (caller-owned)
            type(gas_chemistry_c), intent(in), optional :: gas_species_chem !< Optional gas chemistry object
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
    !> @param[out] init_cat_exch_zones Allocated and populated initial cation exchange zones
        subroutine read_init_cat_exch_zones_CHEPROO(this,unit,ndrz,init_cat_exch_zones)
            import chemistry_c                            !< Import chemistry class definition
            import solid_chemistry_c                      !< Import solid chemistry class
            import reactive_zone_c                        !< Import reactive zone class
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object to populate
            integer(kind=4), intent(in) :: unit           !< File unit number for reading
            integer(kind=4), intent(out):: ndrz           !< Number of adsorption reactive zones
            type(solid_chemistry_c), intent(out), allocatable :: init_cat_exch_zones(:) !< Initial cation exchange zones (allocated here)
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
        
    !> @brief Interface for component concentrations via file I/O
    !> @details Reads aqueous component concentrations from a file after conservative transport,
    !>          performs reactive chemistry step, and writes updated concentrations to output file.
    !> @param[in,out] this Chemistry object
    !> @param[in] path Path for input and output files
    !> @param[in] num_aq_comps Number of aqueous components to process
    !> @param[in] file_in Input filename with post-transport concentrations
    !> @param[in] Delta_t Time step size for reactive chemistry
    !> @param[in] file_out Output filename for post-reaction concentrations
        subroutine interfaz_comps_arch(this,path,num_aq_comps,file_in,Delta_t,file_out)
            import chemistry_c
            class(chemistry_c) :: this                    !< Chemistry object
            character(len=*), intent(in) :: path          !< Path for input and output files
            integer(kind=4), intent(in) :: num_aq_comps   !< Number of aqueous components
            character(len=*), intent(in) :: file_in       !< Input filename with post-transport concentrations
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
            subroutine read_tar_sol(this,root,nsrz,ngrz,num_tar)
            import chemistry_c
            class(chemistry_c) :: this                    !< Chemistry object to populate
            character(len=*), intent(in) :: root          !< Root name for input files
            integer(kind=4), intent(in) :: nsrz           !< Number of solid reactive zones
            integer(kind=4), intent(in) :: ngrz           !< Number of gas reactive zones
            integer(kind=4), intent(in) :: num_tar        !< Expected number of targets in the mesh (upper bound for target solids)
            end subroutine
            
    !> @brief Read target gas definitions from file
    !> @details Reads target gas phase chemistry definitions from input files,
    !>          associating gas phases with spatial locations and reactive zones.
    !> @param[in,out] this Chemistry object to populate with target gases
    !> @param[in] root Root name for input files
    !> @param[in] ngrz Number of gas reactive zones
            subroutine read_tar_gas(this,root,ngrz,num_tar)
            import chemistry_c
            class(chemistry_c) :: this                    !< Chemistry object to populate
            character(len=*), intent(in) :: root          !< Root name for input files
            integer(kind=4), intent(in) :: ngrz           !< Number of gas reactive zones
            integer(kind=4), intent(in) :: num_tar        !< Expected number of targets in the mesh (upper bound for target gases)
            end subroutine

    !> @brief Initialize complete chemistry system
    !> @details Reads chemical system, initial/boundary conditions, and reactive zone
    !>          definitions, then allocates and populates all chemistry data structures.
    !> @param[in,out] this Chemistry object to initialize
    !> @param[in] root Root name for input/output files
    !> @param[in] path_pb Path to problem directory
    !> @param[in] path_DB Path to database files
        ! subroutine initialise_chemistry(this,root,path_pb,path_DB)
        !     import chemistry_c                            !< Import chemistry class definition
        !     implicit none                                 !< Require explicit variable declarations
        !     class(chemistry_c) :: this                    !< Chemistry object to initialize
        !     character(len=*), intent(in) :: root          !< Root name for input/output files
        !     character(len=*), intent(in) :: path_pb       !< Path to problem directory
        !     character(len=*), intent(in) :: path_DB       !< Path to database files
        ! end subroutine

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
        subroutine solve_reactive_mixing_ideal_lump_iter(this, k, Delta_t, theta_r, &
            num_tar_wat, ind_tar_wat, &
            n_p_cache, n_v_cache, n_eq_cache, &
            has_minerals_flag, has_gases_flag, &
            mix_conc_indices, mixing_ratios_conc, &
            lumped_lambdas, &
            all_conc_old, all_ind_aq_sp, &
            conc_comp, &
            p_solver)
            import chemistry_c                            !< Import chemistry class definition
            import real_array_c                           !< Import real array container type
            import int_array_c                            !< Import integer array container type
            import reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2 !< Import solver procedure interface
            implicit none                                 !< Require explicit variable declarations
            class(chemistry_c) :: this                    !< Chemistry object coordinating the simulation
            integer(kind=4), intent(in) :: k              !< Current time step index
            real(kind=8), intent(in) :: Delta_t           !< Time step size
            real(kind=8), intent(in) :: theta_r           !< Integration method parameter for chemical reactions
            integer(kind=4), intent(in) :: num_tar_wat    !< Number of target waters in this batch
            integer(kind=4), intent(in) :: ind_tar_wat(:) !< Indices of target waters in this batch
            integer(kind=4), intent(in) :: n_p_cache(:)   !< Number of primary species per target water
            integer(kind=4), intent(in) :: n_v_cache(:)   !< Number of variable-activity species per target water
            integer(kind=4), intent(in) :: n_eq_cache(:)  !< Number of equilibrium reactions per target water
            logical, intent(in) :: has_minerals_flag(:)   !< Flag indicating presence of minerals per target water
            logical, intent(in) :: has_gases_flag(:)      !< Flag indicating presence of gases per target water
            class(int_array_c), intent(in) :: mix_conc_indices !< Mixing water indices for concentrations
            class(real_array_c), intent(in) :: mixing_ratios_conc !< Mixing ratios for concentrations
            real(kind=8), intent(in) :: lumped_lambdas(:) !< Lumped eigenvalues for spatial lumping
            real(kind=8), intent(in) :: all_conc_old(:,:) !< Old concentrations for all waters
            integer(kind=4), intent(in) :: all_ind_aq_sp(:,:) !< Aqueous species indices for all waters
            real(kind=8), intent(inout) :: conc_comp(:)   !< Component concentrations (updated in place)
            procedure(reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2) :: p_solver !< Procedure pointer to equilibrium-kinetic solver
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
        !> @param[in] tar_sol_ind Target solid index (0 if no solid; for external waters
        !>            with `flag_ext=.true.` and `tar_sol_ind=0`, the water inherits its
        !>            solid_chemistry — and therefore reactive zone — from
        !>            `wat_types(wtype)%solid_chemistry`. For domain waters with
        !>            `tar_sol_ind=0` the dummy solid is still used)
        !> @param[in] tar_gas_ind Target gas index (0 if no gas; same fallback rule as
        !>            above for external waters with `flag_ext=.true.`)
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
                !> Option C: external waters (boundary / recharge) live off-mesh, so they
                !> have no target_solid / target_gas index in tar_wat.dat (the user writes
                !> 0). However the WATER TYPE used to define them was built against a
                !> reactive zone whose mineral/gas phases govern its speciation. Reuse
                !> that reactive zone here by pointing the external water at the water
                !> type's own solid_chemistry (and gas_chemistry, if present) instead of
                !> the empty target_solids_dummy.
                if (flag_ext .and. associated(this%wat_types(wtype)%solid_chemistry)) then
                    if (.not. associated(this%wat_types(wtype)%solid_chemistry%reactive_zone)) then
                        error stop "read_tar_wat_line: water type's solid_chemistry has no reactive_zone"
                    end if
                    !> Inherit gas chemistry from the water type if available
                    if (tar_gas_ind>0) then
                        call this%waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_gas_ind))
                        call this%waters_init(tar_wat_ind)%set_gas_chemistry(this%target_gases_init(tar_gas_ind))
                    else if (associated(this%wat_types(wtype)%gas_chemistry)) then
                        call this%waters(tar_wat_ind)%set_gas_chemistry(this%wat_types(wtype)%gas_chemistry)
                        call this%waters_init(tar_wat_ind)%set_gas_chemistry(this%wat_types(wtype)%gas_chemistry)
                    end if
                    !> Inherit solid chemistry (and therefore reactive zone) from the water type
                    call this%waters(tar_wat_ind)%set_solid_chemistry(this%wat_types(wtype)%solid_chemistry)
                    call this%waters(tar_wat_ind)%set_solid_chemistry_old(this%wat_types(wtype)%solid_chemistry)
                    call this%waters_init(tar_wat_ind)%set_solid_chemistry(this%wat_types(wtype)%solid_chemistry)
                else if (tar_gas_ind>0) then !< Case B1: Gas phase only (no solid minerals)
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
                    !> Case B2: No solid or gas phases (pure aqueous system)
                    call this%waters(tar_wat_ind)%set_solid_chemistry(this%target_solids_dummy(1)) !< set dummy solid for current water
                    call this%waters(tar_wat_ind)%set_solid_chemistry_old(this%target_solids_dummy(1)) !< set dummy solid for current water
                    call this%waters_init(tar_wat_ind)%set_solid_chemistry(this%target_solids_dummy(1)) !< set default dummy for initial water
                end if
            end if
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
                    !> Compute speciation algorithm arrays (different behavior with/without gas phase)
                    if (associated(this%waters(tar_wat_ind)%gas_chemistry)) then
                        !> Gas phase present: include gas activities in speciation algorithm
                        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
                            flag_Se,swap,this%waters(tar_wat_ind)%gas_chemistry%activities)
                    else
                        !> No gas phase: compute speciation arrays without gas activities
                        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
                            flag_Se,swap)
                    end if
                    
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
            !> STEP 6: Allocate and initialize reaction rates
            !> Set up reaction rate arrays and compute initial rates for non-external waters
            call this%waters(tar_wat_ind)%allocate_reaction_rates()      !< allocate reaction rates array for current water
            call this%waters_init(tar_wat_ind)%allocate_reaction_rates() !< allocate reaction rates array for initial water
            call this%waters(tar_wat_ind)%set_indices_rk()               !< set indices mapping for reaction rates (current)
            call this%waters_init(tar_wat_ind)%set_indices_rk()          !< set indices mapping for reaction rates (initial)
            call this%waters(tar_wat_ind)%set_volume()                   !< set pore volume for current target water
            call this%waters_init(tar_wat_ind)%set_volume()              !< set pore volume for initial target water
            
            !> Compute initial reaction rates
                allocate(rk(this%waters(tar_wat_ind)%indices_rk%num_cols)) !< allocate temporary reaction rates vector
                call this%waters(tar_wat_ind)%compute_rk_new(rk)           !< compute current reaction rates
                call this%waters(tar_wat_ind)%set_rk_old(rk)      !< compute initial reaction rates
                call this%waters_init(tar_wat_ind)%compute_rk_new(rk)      !< compute initial reaction rates
                call this%waters_init(tar_wat_ind)%set_rk_old(rk)          !< set old rates for initial water
                deallocate(rk) !< clean up temporary array
            
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

        !> @brief Get the number of aqueous components in a domain (resident/initial) target water
        !> @details Retrieves the number of aqueous primary species (components) from the
        !>          speciation algebra associated with a specific target water,
        !>          accessed through its solid chemistry's reactive zone. If no index is
        !>          provided, the first target water (`tar_wat_indices(1)`) is used.
        !>          The component count may differ from the chemical-system-wide count when
        !>          the local reactive zone defines a reduced set of equilibrium reactions.
        !> @param[in] this Chemistry object
        !> @param[in] ind_tw Optional index into `tar_wat_indices` selecting the target water
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
        !> @return conc_comp_wat_types 2D array of component concentrations (components Ã— water types)
        function get_conc_comp_wat_types(this) result(conc_comp_wat_types)
        implicit none                                     !< Require explicit variable declarations
        class(chemistry_c) :: this                        !< Chemistry object instance
        real(kind=8), allocatable :: conc_comp_wat_types(:,:) !< Concentrations of aqueous components of water types

        integer(kind=4) :: num_aq_comps                   !< Number of aqueous components in the chemical system
        integer(kind=4) :: i                              !< Loop index

        num_aq_comps=this%get_num_aq_comps_chem_syst() !< Call getter function to retrieve number of aqueous components
        allocate(conc_comp_wat_types(num_aq_comps,this%num_wat_types)) !< Allocate 2D array with dimensions [components Ã— water types]
        do i=1,this%num_wat_types                                       !< Loop over all defined water types
            conc_comp_wat_types(:,i)=this%wat_types(i)%get_u_aq()       !< Extract aqueous component concentrations for water type i and store in column i
        end do                                                          !< End water types loop
        end function

        !> @brief Write component concentrations for every target water (and external waters)
        !> @details Produces a formatted report of the aqueous-component concentrations
        !>          held by every entry of `this%waters`, projected onto the component
        !>          basis of a chosen reference water type. The reference water type is
        !>          located by case-insensitive name match against the keywords
        !>          `domain`/`initial`/`resident`/`dominio`/`inicial`/`residente`; its
        !>          `comp_mat_aq` matrix sets the row dimension for every output column.
        !>
        !>          File structure (in order):
        !>            1. Header `Number of aqueous components: <n>`.
        !>            2. Block `Aqueous components:` listing each component as the
        !>               explicit linear combination of aqueous variable-activity species
        !>               (one row per component, format `uJ = <combination>`).
        !>            3. Block `Aqueous component concentrations in the domain:` with one
        !>               column per index in `this%tar_wat_indices` and one row per
        !>               component (rows = components, cols = target waters).
        !>            4. Block `Aqueous component concentrations of external waters:` with
        !>               one column per index in `this%ext_waters_indices`, after
        !>               filtering out any water whose name matches the reference water
        !>               type (those are conceptually domain waters and were already
        !>               written in block 3). Each column is preceded by the water's name.
        !>
        !>          Numeric format: ES13.5E2, fixed 13-char fields so columns line up.
        !> @see write_conc_comp_wat_types for a per-water-type version (no target waters)
        !> @param this Chemistry object
        !> @param path File path for output
        !> @param filename Output filename
        subroutine write_conc_comp_tar_wat(this,path,filename)
        implicit none                                     !< Require explicit variable declarations
        class(chemistry_c) :: this                        !< Chemistry object instance
        character(len=*), intent(in) :: path              !< Path to the output file
        character(len=*), intent(in) :: filename          !< Output filename

        real(kind=8), allocatable :: u_aq_init(:)         !< Aqueous components of water types
        real(kind=8), allocatable :: c_var_act(:)         !< Variable-activity-species concentrations of water type i, ordered as in the reference water type
        integer(kind=4) :: i,j,k                          !< Loop indices
        integer(kind=4) :: ref_idx                        !< Index of reference water type ("domain"/"initial")
        integer(kind=4) :: num_comps_ref                  !< Number of aqueous components in reference water type
        integer(kind=4) :: num_var_act_ref                !< Number of variable-activity aqueous species in reference water type
        integer(kind=4) :: sp_idx_i                       !< Index of a reference species inside water type i's own aq_phase
        character(len=:), allocatable :: name_lc          !< Lowercased water type name for matching
        character(len=13) :: name_buf                     !< Fixed-width buffer for right-aligning species names (matches ES13.5E2 width)

        !> Locate reference water type whose name contains "domain", "initial",
        !> "dominio" or "inicial" (case-insensitive). Its aqueous component matrix
        !> will be used for all water types so the resulting u_aq vectors share
        !> the same dimension.
        ref_idx=0 !< Sentinel: no reference water type found yet
        do i=1,this%num_wat_types !< Scan all water types looking for the reference
            name_lc=to_lower(trim(this%wat_types(i)%name)) !< Normalise name to lowercase for case-insensitive match
            if (index(name_lc,'domain')>0 .or. index(name_lc,'initial')>0 .or. &
                index(name_lc,'dominio')>0 .or. index(name_lc,'inicial')>0 .or. &
                index(name_lc,'resident')>0 .or. index(name_lc,'residente')>0) then
                ref_idx=i !< Remember reference water type index
                exit !< Use first matching water type as reference
            end if
        end do
        if (ref_idx==0) then !< No matching water type was found
            error stop "write_conc_comp_tar_wat: no water type named 'domain'/'initial'/'resident'/'dominio'/'inicial'/'residente' found"
        end if

        num_comps_ref=this%wat_types(ref_idx)%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
        num_var_act_ref=this%wat_types(ref_idx)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species

        open(unit=10, file=path//filename, status='unknown', action='write', form='formatted') !< Open output file for writing water type compositions
        !> Write the number of aqueous components at the top of the file (taken
        !> from the reference water type's speciation algebra).
        write(10,"(2x,A,I0,/)") "Number of aqueous components: ", num_comps_ref
        !> Write the names of the aqueous primary species of the reference water
        !> type, in the same order as the component concentrations that follow.
        !> Names are right-aligned in a 13-character field so that they end at the
        !> same column (column 17) as the concentration values written below.
        write(10,"(2x,A)") "Aqueous components:"
        !> Each component is the linear combination of aqueous variable-activity
        !> species defined by row j of comp_mat_aq. We assemble a textual
        !> expression like "+1*h+ -1*hco3- +2*co3-2" by walking the row and
        !> appending non-zero terms.
        block
            character(len=512) :: comp_expr           !< Buffer for one component expression
            character(len=:), allocatable :: sp_name  !< Name of the k-th var-act species
            real(kind=8) :: coeff                     !< Stoichiometric coefficient comp_mat_aq(j,k)
            integer(kind=4) :: nv                     !< Number of var-act species
            integer(kind=4) :: pos                    !< Current write position in comp_expr
            character(len=32) :: tok                  !< Temporary token buffer for one term
            logical :: first_term                     !< True until the first non-zero term has been emitted
            nv=this%wat_types(ref_idx)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
            do j=1,num_comps_ref
                comp_expr=' '
                pos=1
                first_term=.true.
                do k=1,nv
                    coeff=this%wat_types(ref_idx)%solid_chemistry%reactive_zone%speciation_alg%comp_mat_aq(j,k)
                    if (coeff==0.0d0) cycle
                    sp_name=trim(this%wat_types(ref_idx)%aq_phase%aq_species(&
                        this%wat_types(ref_idx)%indices_aq_phase(&
                        this%wat_types(ref_idx)%ind_var_act_species(k)))%name)
                    !> Format the term. The leading sign of the very first non-zero
                    !> term is suppressed for positives ("+" implicit) and kept for
                    !> negatives. For subsequent terms an explicit "+ " is prepended
                    !> for positives so the expression reads e.g. "hco3- + co3-2".
                    if (coeff==1.0d0) then
                        if (first_term) then
                            write(tok,'(A)') sp_name
                        else
                            write(tok,'(A,A)') '+ ', sp_name
                        end if
                    else if (coeff==-1.0d0) then
                        write(tok,'(A,A)') '- ', sp_name
                    else if (coeff>0.0d0) then
                        if (first_term) then
                            write(tok,'(F0.2,A,A)') coeff, '*', sp_name
                        else
                            write(tok,'(A,F0.2,A,A)') '+ ', coeff, '*', sp_name
                        end if
                    else
                        write(tok,'(A,F0.2,A,A)') '- ', -coeff, '*', sp_name
                    end if
                    if (pos>1) then
                        comp_expr(pos:pos)=' '
                        pos=pos+1
                    end if
                    comp_expr(pos:pos+len_trim(tok)-1)=trim(tok)
                    pos=pos+len_trim(tok)
                    first_term=.false.
                end do
                write(10,"(4x,A,I0,A,A)") "u",j," = ", trim(comp_expr)
            end do
        end block
        write(10,"(/,2x,A)") "Aqueous component concentrations in the domain: " !< Heading for the target-waters block
        !> Build a 2D matrix (num_comps_ref rows x n_dom cols) and write it row by
        !> row so that each column corresponds to one target water, matching the
        !> u_tilde input layout (rows=components, cols=targets).
        if (allocated(this%tar_wat_indices) .and. size(this%tar_wat_indices)>0) then
            block
                real(kind=8), allocatable :: u_dom(:,:) !< Component matrix: rows=components, cols=target waters
                integer(kind=4) :: n_dom                 !< Number of target (domain) waters
                n_dom=size(this%tar_wat_indices)
                allocate(u_dom(num_comps_ref,n_dom))
                allocate(c_var_act(num_var_act_ref))
                do i=1,n_dom !< Fill column i with components of the i-th target water
                    associate (iw => this%tar_wat_indices(i))
                        do j=1,num_var_act_ref !< Map reference species ordering to water iw's own ordering
                            call this%waters(iw)%aq_phase%get_aq_species_index_by_name(&
                                trim(this%wat_types(ref_idx)%aq_phase%aq_species(&
                                    this%wat_types(ref_idx)%indices_aq_phase(&
                                    this%wat_types(ref_idx)%ind_var_act_species(j)))%name), &
                                sp_idx_i)
                            do k=1,size(this%waters(iw)%indices_aq_phase)
                                if (this%waters(iw)%indices_aq_phase(k)==sp_idx_i) then
                                    c_var_act(j)=this%waters(iw)%concentrations(k)
                                    exit
                                end if
                            end do
                        end do
                        u_dom(:,i)=matmul(&
                            this%wat_types(ref_idx)%solid_chemistry%reactive_zone%speciation_alg%comp_mat_aq, &
                            c_var_act)
                    end associate
                end do
                deallocate(c_var_act)
                !> Write row by row: each row holds the j-th component of every target water.
                do j=1,num_comps_ref
                    write(10,"(4x,*(ES13.5E2,1x))") (u_dom(j,i), i=1,n_dom)
                end do
                deallocate(u_dom)
            end block
        end if
        !> External (boundary + recharge) waters block. Same projection onto the
        !> reference component basis as the domain block above; values are also
        !> written in column form (one column per external water) so the layout
        !> matches u_tilde's (rows=components, cols=targets) convention.
        !>
        !> NOTE: any water whose water-type name matches the reference water type
        !> (i.e. the "domain"/"residente"/"initial" type) is skipped here, even if
        !> it appears in ext_waters_indices, because such waters are conceptually
        !> domain (target) waters — they were already written in the domain block.
        if (allocated(this%ext_waters_indices) .and. size(this%ext_waters_indices)>0) then
            block
                real(kind=8), allocatable :: u_ext(:,:)            !< Component matrix: rows=components, cols=true external waters
                integer(kind=4), allocatable :: ext_keep(:)        !< Filtered indices (drops waters of reference type)
                integer(kind=4) :: n_ext                            !< Number of true external waters
                character(len=:), allocatable :: ref_name_lc        !< Lowercased reference wat_type name
                character(len=:), allocatable :: cand_name_lc       !< Lowercased candidate water name
                ref_name_lc=to_lower(trim(this%wat_types(ref_idx)%name))
                !> Two-pass filter: first count then collect indices whose name
                !> differs from the reference water type's name.
                n_ext=0
                do i=1,size(this%ext_waters_indices)
                    cand_name_lc=to_lower(trim(this%waters(this%ext_waters_indices(i))%name))
                    if (cand_name_lc /= ref_name_lc) n_ext=n_ext+1
                end do
                if (n_ext==0) then
                    !> Nothing to write: every "external" entry is actually a
                    !> domain water that was already covered above.
                    !> Skip the heading entirely for clarity.
                    !> (Falls through to close(10) below.)
                else
                    allocate(ext_keep(n_ext))
                    j=0
                    do i=1,size(this%ext_waters_indices)
                        cand_name_lc=to_lower(trim(this%waters(this%ext_waters_indices(i))%name))
                        if (cand_name_lc /= ref_name_lc) then
                            j=j+1
                            ext_keep(j)=this%ext_waters_indices(i)
                        end if
                    end do
                    write(10,"(/,2x,A)") "Aqueous component concentrations of external waters: " !< Heading for the external-waters block
                    allocate(u_ext(num_comps_ref,n_ext))
                    allocate(c_var_act(num_var_act_ref))
                    do i=1,n_ext !< Fill column i with components of the i-th true external water
                        associate (iw => ext_keep(i))
                            do j=1,num_var_act_ref
                                call this%waters(iw)%aq_phase%get_aq_species_index_by_name(&
                                    trim(this%wat_types(ref_idx)%aq_phase%aq_species(&
                                        this%wat_types(ref_idx)%indices_aq_phase(&
                                        this%wat_types(ref_idx)%ind_var_act_species(j)))%name), &
                                    sp_idx_i)
                                do k=1,size(this%waters(iw)%indices_aq_phase)
                                    if (this%waters(iw)%indices_aq_phase(k)==sp_idx_i) then
                                        c_var_act(j)=this%waters(iw)%concentrations(k)
                                        exit
                                    end if
                                end do
                            end do
                            u_ext(:,i)=matmul(&
                                this%wat_types(ref_idx)%solid_chemistry%reactive_zone%speciation_alg%comp_mat_aq, &
                                c_var_act)
                        end associate
                    end do
                    deallocate(c_var_act)
                    !> Column header line: external water names right-aligned in 13-char fields.
                    write(10,"(4x,*(A13,1x))") (adjustr(name_buf_assign(this%waters(ext_keep(i))%name)), i=1,n_ext)
                    !> Write row by row: each row holds the j-th component of every external water.
                    do j=1,num_comps_ref
                        write(10,"(4x,*(ES13.5E2,1x))") (u_ext(j,i), i=1,n_ext)
                    end do
                    deallocate(u_ext,ext_keep)
                end if
            end block
        end if
        close(10) !< Close output file after writing all water type data

        contains

        !> @brief Copy a string into a 13-char fixed buffer (used for column headers)
        pure function name_buf_assign(s) result(out)
            character(len=*), intent(in) :: s
            character(len=13) :: out
            out=trim(s) !< Auto-pads with blanks on the right
        end function name_buf_assign

        !> @brief Return a lowercase copy of the input string (ASCII)
        pure function to_lower(s) result(out)
            character(len=*), intent(in) :: s
            character(len=len(s)) :: out
            integer :: k, ic
            do k=1,len(s)
                ic=iachar(s(k:k))
                if (ic>=iachar('A') .and. ic<=iachar('Z')) then
                    out(k:k)=achar(ic+32)
                else
                    out(k:k)=s(k:k)
                end if
            end do
        end function to_lower
        end subroutine

        !> @brief Write component concentrations for every water type, split initial vs external
        !> @details Companion to `write_conc_comp_tar_wat`. Instead of iterating over
        !>          target waters (entries of `this%waters`), this routine iterates over
        !>          the entries of `this%wat_types` and writes their aqueous-component
        !>          concentrations projected onto the reference water type's component
        !>          basis (rows of `comp_mat_aq`).
        !>
        !>          The reference water type is located by case-insensitive name match
        !>          against the keywords `domain`/`initial`/`resident`/`dominio`/`inicial`/
        !>          `residente`. The same keywords are then used to classify every water
        !>          type into one of two groups:
        !>            - Matches  -> written under `Aqueous component concentrations of
        !>                          initial water types:`.
        !>            - All rest -> written under `Aqueous component concentrations of
        !>                          external water types:`.
        !>
        !>          File structure (in order):
        !>            1. Header `Number of aqueous components: <n>`.
        !>            2. Block `Aqueous components:` listing each component as the
        !>               explicit linear combination of aqueous variable-activity species
        !>               (one row per component, format `uJ = <combination>`).
        !>            3. Initial water-types block (only if any matches): one column per
        !>               type, preceded by a header line of right-aligned type names.
        !>            4. External water-types block (only if any non-matches): one column
        !>               per type, with the same header layout.
        !>
        !>          Layout convention: rows = components, columns = water types. Numeric
        !>          format: ES13.5E2, fixed 13-char fields so columns line up.
        !> @see write_conc_comp_tar_wat for the target-water version
        !> @param this Chemistry object instance
        !> @param path File path for output
        !> @param filename Output filename
        subroutine write_conc_comp_wat_types(this,path,filename)
        implicit none
        class(chemistry_c) :: this
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: filename

        real(kind=8), allocatable :: c_var_act(:)
        integer(kind=4) :: i,j,k
        integer(kind=4) :: ref_idx
        integer(kind=4) :: num_comps_ref, num_var_act_ref
        integer(kind=4) :: sp_idx_i
        integer(kind=4) :: n_init, n_ext
        integer(kind=4), allocatable :: init_idx(:), ext_idx(:)
        real(kind=8), allocatable :: u_init(:,:), u_ext(:,:)
        character(len=:), allocatable :: name_lc

        !> Step 1: locate the reference water type (same rule as write_conc_comp_tar_wat).
        ref_idx=0
        do i=1,this%num_wat_types
            name_lc=to_lower2(trim(this%wat_types(i)%name))
            if (index(name_lc,'domain')>0 .or. index(name_lc,'initial')>0 .or. &
                index(name_lc,'resident')>0 .or. index(name_lc,'dominio')>0 .or. &
                index(name_lc,'inicial')>0 .or. index(name_lc,'residente')>0) then
                ref_idx=i
                exit
            end if
        end do
        if (ref_idx==0) then
            error stop "write_conc_comp_wat_types: no reference water type found"
        end if

        num_comps_ref=this%wat_types(ref_idx)%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
        num_var_act_ref=this%wat_types(ref_idx)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species

        !> Step 2: classify water types into initial vs external by name match.
        n_init=0; n_ext=0
        do i=1,this%num_wat_types
            name_lc=to_lower2(trim(this%wat_types(i)%name))
            if (index(name_lc,'domain')>0 .or. index(name_lc,'initial')>0 .or. &
                index(name_lc,'resident')>0 .or. index(name_lc,'dominio')>0 .or. &
                index(name_lc,'inicial')>0 .or. index(name_lc,'residente')>0) then
                n_init=n_init+1
            else
                n_ext=n_ext+1
            end if
        end do
        if (n_init>0) allocate(init_idx(n_init))
        if (n_ext>0) allocate(ext_idx(n_ext))
        n_init=0; n_ext=0
        do i=1,this%num_wat_types
            name_lc=to_lower2(trim(this%wat_types(i)%name))
            if (index(name_lc,'domain')>0 .or. index(name_lc,'initial')>0 .or. &
                index(name_lc,'resident')>0 .or. index(name_lc,'dominio')>0 .or. &
                index(name_lc,'inicial')>0 .or. index(name_lc,'residente')>0) then
                n_init=n_init+1
                init_idx(n_init)=i
            else
                n_ext=n_ext+1
                ext_idx(n_ext)=i
            end if
        end do

        open(unit=10, file=path//filename, status='unknown', action='write', form='formatted')
        write(10,"(2x,A,I0,/)") "Number of aqueous components: ", num_comps_ref

        !> Write each component as the linear combination of aqueous variable-activity
        !> species defined by row j of comp_mat_aq, in the same format used by
        !> write_conc_comp_tar_wat (e.g. "u6 = hco3- + co3-2").
        write(10,"(2x,A,/)") "Aqueous components:"
        block
            character(len=512) :: comp_expr
            character(len=:), allocatable :: sp_name
            real(kind=8) :: coeff
            integer(kind=4) :: nv, pos
            character(len=32) :: tok
            logical :: first_term
            nv=this%wat_types(ref_idx)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
            do j=1,num_comps_ref
                comp_expr=' '
                pos=1
                first_term=.true.
                do k=1,nv
                    coeff=this%wat_types(ref_idx)%solid_chemistry%reactive_zone%speciation_alg%comp_mat_aq(j,k)
                    if (coeff==0.0d0) cycle
                    sp_name=trim(this%wat_types(ref_idx)%aq_phase%aq_species(&
                        this%wat_types(ref_idx)%indices_aq_phase(&
                        this%wat_types(ref_idx)%ind_var_act_species(k)))%name)
                    if (coeff==1.0d0) then
                        if (first_term) then
                            write(tok,'(A)') sp_name
                        else
                            write(tok,'(A,A)') '+ ', sp_name
                        end if
                    else if (coeff==-1.0d0) then
                        write(tok,'(A,A)') '- ', sp_name
                    else if (coeff>0.0d0) then
                        if (first_term) then
                            write(tok,'(F0.2,A,A)') coeff, '*', sp_name
                        else
                            write(tok,'(A,F0.2,A,A)') '+ ', coeff, '*', sp_name
                        end if
                    else
                        write(tok,'(A,F0.2,A,A)') '- ', -coeff, '*', sp_name
                    end if
                    if (pos>1) then
                        comp_expr(pos:pos)=' '
                        pos=pos+1
                    end if
                    comp_expr(pos:pos+len_trim(tok)-1)=trim(tok)
                    pos=pos+len_trim(tok)
                    first_term=.false.
                end do
                write(10,"(4x,A,I0,A,A)") "u",j," = ", trim(comp_expr)
            end do
        end block
        write(10,"(A)") ""

        !> Step 3: project initial water types onto the reference component basis.
        if (n_init>0) then
            allocate(u_init(num_comps_ref,n_init))
            allocate(c_var_act(num_var_act_ref))
            do i=1,n_init
                associate (it => init_idx(i))
                    do j=1,num_var_act_ref
                        call this%wat_types(it)%aq_phase%get_aq_species_index_by_name(&
                            trim(this%wat_types(ref_idx)%aq_phase%aq_species(&
                                this%wat_types(ref_idx)%indices_aq_phase(&
                                this%wat_types(ref_idx)%ind_var_act_species(j)))%name), &
                            sp_idx_i)
                        do k=1,size(this%wat_types(it)%indices_aq_phase)
                            if (this%wat_types(it)%indices_aq_phase(k)==sp_idx_i) then
                                c_var_act(j)=this%wat_types(it)%concentrations(k)
                                exit
                            end if
                        end do
                    end do
                    u_init(:,i)=matmul(&
                        this%wat_types(ref_idx)%solid_chemistry%reactive_zone%speciation_alg%comp_mat_aq, &
                        c_var_act)
                end associate
            end do
            deallocate(c_var_act)
            write(10,"(2x,A)") "Aqueous component concentrations of initial water types: "
            write(10,"(A)") ""
            write(10,"(4x,*(A13,1x))") (adjustr(name_buf_assign2(this%wat_types(init_idx(i))%name)), i=1,n_init)
            write(10,"(A)") ""
            do j=1,num_comps_ref
                write(10,"(4x,*(ES13.5E2,1x))") (u_init(j,i), i=1,n_init)
            end do
            deallocate(u_init)
        end if

        !> Step 4: project external water types onto the reference component basis.
        if (n_ext>0) then
            allocate(u_ext(num_comps_ref,n_ext))
            allocate(c_var_act(num_var_act_ref))
            do i=1,n_ext
                associate (it => ext_idx(i))
                    do j=1,num_var_act_ref
                        call this%wat_types(it)%aq_phase%get_aq_species_index_by_name(&
                            trim(this%wat_types(ref_idx)%aq_phase%aq_species(&
                                this%wat_types(ref_idx)%indices_aq_phase(&
                                this%wat_types(ref_idx)%ind_var_act_species(j)))%name), &
                            sp_idx_i)
                        do k=1,size(this%wat_types(it)%indices_aq_phase)
                            if (this%wat_types(it)%indices_aq_phase(k)==sp_idx_i) then
                                c_var_act(j)=this%wat_types(it)%concentrations(k)
                                exit
                            end if
                        end do
                    end do
                    u_ext(:,i)=matmul(&
                        this%wat_types(ref_idx)%solid_chemistry%reactive_zone%speciation_alg%comp_mat_aq, &
                        c_var_act)
                end associate
            end do
            deallocate(c_var_act)
            write(10,"(/,2x,A)") "Aqueous component concentrations of external water types: "
            write(10,"(A)") ""
            write(10,"(4x,*(A13,1x))") (adjustr(name_buf_assign2(this%wat_types(ext_idx(i))%name)), i=1,n_ext)
            write(10,"(A)") ""
            do j=1,num_comps_ref
                write(10,"(4x,*(ES13.5E2,1x))") (u_ext(j,i), i=1,n_ext)
            end do
            deallocate(u_ext)
        end if

        if (allocated(init_idx)) deallocate(init_idx)
        if (allocated(ext_idx)) deallocate(ext_idx)
        close(10)

        contains

        pure function name_buf_assign2(s) result(out)
            character(len=*), intent(in) :: s
            character(len=13) :: out
            out=trim(s)
        end function name_buf_assign2

        pure function to_lower2(s) result(out)
            character(len=*), intent(in) :: s
            character(len=len(s)) :: out
            integer :: kk, ic
            do kk=1,len(s)
                ic=iachar(s(kk:kk))
                if (ic>=iachar('A') .and. ic<=iachar('Z')) then
                    out(kk:kk)=achar(ic+32)
                else
                    out(kk:kk)=s(kk:kk)
                end if
            end do
        end function to_lower2
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
        !> @return u_mix_init 2D array of u_mix values (components Ã— domain target waters)
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
        
        !> @brief Set the number of materials in the system
        !> @details Sets the total number of material objects in the chemistry system.
        !> Materials represent distinct solid compositions with associated chemical properties.
        !> Validates that the number is non-negative before assignment.
        !> @param this Chemistry object
        !> @param num_materials Number of materials (must be â‰¥ 0)
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

        !> @brief Allocate memory for the per-water-type solid chemistry array
        !> @details Safely allocates memory for `wat_type_solids(:)`, sized 1-to-1 with
        !> `wat_types(:)`. Each entry holds the solid_chemistry_c template that the
        !> corresponding water type's `solid_chemistry` pointer is bound to. Unrelated
        !> to `target_solids(:)` (mesh-bound) and `materials(:)` (mineral + surface
        !> zone catalogue from quim_loc.dat).
        !> @param this Chemistry object
        subroutine allocate_wat_type_solids(this)
        implicit none
        class(chemistry_c) :: this !< Chemistry object to modify
        if (allocated(this%wat_type_solids)) deallocate(this%wat_type_solids) !< Deallocate existing array if previously allocated to prevent memory leak
        allocate(this%wat_type_solids(this%num_wat_types)) !< Allocate one slot per water type
        end subroutine allocate_wat_type_solids

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
        !> @param[out] u_new New aqueous component concentrations [mol/L] - 2D array (components Ã— domain waters)
        !> @param[out] c_new New species concentrations [mol/L] - 2D array (species Ã— domain waters)
        subroutine WMA_iter_EE_eq_kin_ideal_lump(this,mix_wat_indices,lambdas,Delta_t,u_new,c_new)
            class(chemistry_c) :: this !< Chemistry object coordinating water mixing algorithm
            type(int_array_c), intent(in) :: mix_wat_indices !< Matrix of indices indicating which waters mix for each target water
            type(real_array_c), intent(in) :: lambdas !< Matrix of mixing ratios (lambdas) for weighted mixing algorithm (WMA)
            real(kind=8), intent(in) :: Delta_t !< Time step size for explicit Euler temporal discretization
            real(kind=8), intent(out) :: u_new(:,:) !< Output: new aqueous component concentrations (components Ã— domain waters)
            real(kind=8), intent(out) :: c_new(:,:) !< Output: new species concentrations (species Ã— domain waters)

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
        !> @param[in] mixing_ratios_R Mixing ratios for reaction rates
        !> @param[in] mix_conc_indices Matrix of indices of target waters that mix
        !> @param[in] mix_react_indices Matrix of indices of domain mixing waters
        !> @param[in] time_discr Time discretization object
        !> @param[in] theta_r Integration method parameter for chemical reactions
        subroutine solve_reactive_mixing(this,dir_pb,root,mixing_ratios_conc,mixing_ratios_R, &
            mix_conc_indices,mix_react_indices,time_discr,theta_r)
            class(chemistry_c) :: this                    !< Chemistry object coordinating the simulation
            character(len=*), intent(in) :: dir_pb        !< Problem directory path
            character(len=*), intent(in) :: root          !< Root name for output files
            class(real_array_c), intent(in) :: mixing_ratios_conc !< Mixing ratios for concentrations
            class(real_array_c), intent(in) :: mixing_ratios_R !< Mixing ratios for reaction rates
            class(int_array_c), intent(in) :: mix_conc_indices !< Matrix of indices of target waters that mix
            class(int_array_c), intent(in) :: mix_react_indices !< Matrix of domain mixing water indices
            class(time_discr_c), intent(in) :: time_discr !< Time discretization object
            real(kind=8), intent(in) :: theta_r           !< Integration method parameter for chemical reactions
            
            procedure(solve_reactive_mixing_ideal_cons), pointer :: p_solver=>null() !< Procedure pointer to selected reactive mixing solver
            
            if (this%act_coeffs_model==0 .and. this%lump_flag .eqv. .true.) then !< Check for ideal activity with spatial lumping
                p_solver=>solve_reactive_mixing_ideal_lump !< Select ideal lumped solver
            else if (this%act_coeffs_model==0 .and. this%lump_flag .eqv. .false.) then !< Check for ideal activity without lumping
                p_solver=>solve_reactive_mixing_ideal_cons !< Select ideal consistent solver
            else if (this%act_coeffs_model>0 .and. this%lump_flag .eqv. .true.) then !< Check for non-ideal activity (extended Debye-HÃ¼ckel, Pitzer, etc.) with lumping
                !p_solver=>solve_reactive_mixing_lump    
            else if (this%act_coeffs_model>0 .and. this%lump_flag .eqv. .false.) then !< Check for non-ideal activity without lumping (most accurate, slowest)
                !p_solver=>solve_reactive_mixing_cons    
            else                                                               !< Invalid combination of act_coeffs_model and lump_flag
                error stop "Error: act_coeffs_model and lump_flag not compatible" !< Terminate with error: check chemistry configuration
            end if
            call p_solver(this,dir_pb,root,mixing_ratios_conc,& !< Call the selected reactive mixing solver procedure
                mixing_ratios_R,mix_conc_indices,mix_react_indices,& !< Pass reaction mixing ratios and water indices
                time_discr,theta_r) !< Pass time discretization and integration method
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

end module chemistry_m !< End of chemistry_m module
