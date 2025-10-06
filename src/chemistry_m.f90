!> Chemistry class (main class): contains all chemical information and solves reactive mixing
module chemistry_m
    use aqueous_chemistry_m, only: aqueous_chemistry_c, solid_chemistry_c, gas_chemistry_c, int_array_c, real_array_c, &
        inf_norm_vec_real, LU_lin_syst, compute_c_tilde
    use reactive_zone_m, only: reactive_zone_c, compare_react_zones, chem_system_c, CV_params_s, gas_phase_c, gas_c, &
        aq_phase_c
    use mineral_zone_m, only: mineral_zone_c, mineral_c
    use chem_out_options_m, only: chem_out_options_t
    use spatial_discr_1D_m, only: spatial_discr_c, mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c
    use time_discr_m, only: time_discr_c, time_discr_homog_c, time_discr_heterog_c
    implicit none
    save
    type, public :: chemistry_c
        integer(kind=4) :: num_lump=0 !> number of lumpings in consistent WMA option
        integer(kind=4) :: read_opt !> option for reading chemical data
        integer(kind=4) :: act_coeffs_model !> model to compute activity coefficients
        logical :: lump_flag !> TRUE if lumping applied on reaction term, FALSE otherwise
        integer(kind=4) :: cons_opt !> consistent WMA option, either 1 (explicit) or 2 (upstream to downstream) (see documentation for more details)
        integer(kind=4) :: rk_down_opt !> estimation of downstream waters rk, either 1, 2 or 3 (see documentation for more details)
        integer(kind=4) :: rk_avg_opt !> average rk, either 1 (concentrations) or 2 (reaction rates) (see documentation for more details)
        type(chem_system_c) :: chem_syst !> chemical system object
        integer(kind=4) :: num_wat_types=0 !> number of water types
        integer(kind=4) :: num_target_waters=0 !> number of target waters
        integer(kind=4) :: num_target_waters_init=0 !> number of target waters initially
        type(aqueous_chemistry_c), allocatable :: wat_types(:) !> initial & boundary water types
        type(aqueous_chemistry_c), allocatable :: target_waters(:) !> target waters initial
        type(aqueous_chemistry_c), allocatable :: target_waters_init(:) !> target waters initial
        integer(kind=4) :: num_target_waters_dom=0 !> number of initial target waters
        integer(kind=4) :: num_ext_waters=0 !> number of external waters
        integer(kind=4) :: num_rech_waters=0 !> number of recharge waters
        integer(kind=4) :: num_bd_waters=0 !> number of boundary waters
        integer(kind=4), allocatable :: ext_waters_indices(:) !> external waters indices
        integer(kind=4), allocatable :: rech_waters_indices(:) !> recharge waters indices
        integer(kind=4), allocatable :: bd_waters_indices(:) !> boundary waters indices
        integer(kind=4), allocatable :: dom_tar_wat_indices(:) !> domain target waters indices
        integer(kind=4) :: num_target_solids=0 !> number of target solids (<= num_target_waters)
        integer(kind=4) :: num_target_solids_dom=0 !> number of target solids init (<= num_target_waters_dom)
        integer(kind=4) :: num_materials=0 !> number of materials (<= num_target_solids)
        type(solid_chemistry_c), allocatable :: materials(:) !> materials (solid types)
        type(solid_chemistry_c), allocatable :: sol_types(:) !> solid types associated to initial & boundary water types (chapuza)
        type(solid_chemistry_c), allocatable :: target_solids(:) !> target solids
        type(solid_chemistry_c), allocatable :: target_solids_init(:) !> initial target solids
        integer(kind=4) :: num_target_gases=0 !> number of target gases
        integer(kind=4) :: num_gas_zones=0 !> number of gas zones
        type(gas_chemistry_c), allocatable :: gas_zones(:) !> gas zones
        type(gas_chemistry_c), allocatable :: gas_zones_wat_types(:) !> gas zones associated to water types (chapuza)
        type(gas_chemistry_c), allocatable :: target_gases(:) !> target gases
        type(gas_chemistry_c), allocatable :: target_gases_init(:) !> target gases init
        integer(kind=4) :: num_reactive_zones=0 !> number of reactive zones (<=num_target_solids)
        type(reactive_zone_c), allocatable :: reactive_zones(:) !> reactive zones
        type(reactive_zone_c), allocatable :: react_zones_wat_types(:) !> reactive zones associated to water types (chapuza)
        integer(kind=4) :: num_mineral_zones=0 !> number of mineral zones (<=num_target_solids)
        type(mineral_zone_c), allocatable :: mineral_zones(:) !> reactive zones
        integer(kind=4) :: Jac_opt !> model to compute Jacobians (0: incremental coefficnets, 1: analytical)
        type(CV_params_s) :: CV_params !> parameters for convergence
        type(chem_out_options_t) :: chem_out_options !> chemistry output options variable
    contains
    !> Set
        procedure, public :: set_lump_flag
        procedure, public :: set_cons_opt
        procedure, public :: set_rk_avg_opt
        procedure, public :: set_rk_down_opt
        procedure, public :: set_chem_syst
        procedure, public :: set_num_wat_types
        procedure, public :: set_num_materials
        procedure, public :: set_num_tar_wat
        procedure, public :: set_num_tar_wat_dom
        procedure, public :: set_num_target_solids
        procedure, public :: set_num_ext_waters
        procedure, public :: set_num_rech_waters
        procedure, public :: set_num_bd_waters
        procedure, public :: set_target_waters_wat_types
        procedure, public :: set_target_waters_target_solids
        procedure, public :: set_target_solids_materials
        procedure, public :: set_target_gases
        procedure, public :: set_read_opt
        procedure, public :: set_reactive_zones
        procedure, public :: set_Jac_opt
        procedure, public :: set_target_solids_mesh !> chapuza
    !> Get
        procedure, public :: get_num_aq_comps_chem_syst
        procedure, public :: get_num_aq_comps_dom
        procedure, public :: get_num_wat_types
        procedure, public :: get_conc_comp_wat_types
        procedure, public :: get_num_aq_var_act_species_dom
        procedure, public :: get_tar_sol_ind
        procedure, public :: get_tar_gas_ind
    !> Allocate
        procedure, public :: allocate_target_waters
        procedure, public :: allocate_target_waters_init
        procedure, public :: allocate_dom_tar_wat_indices
        procedure, public :: allocate_rech_waters_indices
        procedure, public :: allocate_bd_waters_indices
        procedure, public :: allocate_target_solids
        procedure, public :: allocate_target_gases
        procedure, public :: allocate_reactive_zones
        procedure, public :: allocate_react_zones_wat_types
        procedure, public :: allocate_mineral_zones
        procedure, public :: allocate_materials
        procedure, public :: allocate_wat_types
        procedure, public :: allocate_sol_types
        procedure, public :: allocate_gas_zones
        procedure, public :: allocate_gas_zones_wat_types
    !> Read
        procedure, public :: read_target_waters_init
        procedure, public :: read_tar_wat_line
        procedure, public :: read_tar_sol
        procedure, public :: read_tar_gas
        procedure, public :: read_init_min_zones_CHEPROO
        procedure, public :: read_chemistry
        procedure, public :: read_chemistry_CHEPROO
        procedure, public :: read_init_bd_wat_types_CHEPROO
        procedure, public :: read_init_cat_exch_zones_CHEPROO
        procedure, public :: read_init_gas_zones_CHEPROO
        procedure, public :: read_chem_opts
    !> Initialisation
        procedure, public :: initialise_chemistry
    !> Write
        procedure, public :: write_chemistry
        procedure, public :: write_conc_comp_wat_types
        procedure, public :: write_u_tilde_init
    !> Solve
        procedure, public :: solve_reactive_mixing_lump !> main solver
        procedure, public :: solve_reactive_mixing_ideal_cons !> main solver
        procedure, public :: solve_reactive_mixing_cons !> main solver
        procedure, public :: solve_reactive_mixing_ideal_lump !> main solver
        !procedure, public :: solve_reactive_mixing_bis !> main solver
        procedure, public :: solve_reactive_mixing_BCs_dep_t !> main solver
    !> Link
        procedure, public :: link_target_waters_target_solids
        procedure, public :: link_target_waters_target_gases
        procedure, public :: link_target_solids_reactive_zone
        procedure, public :: link_target_gases_reactive_zone
        procedure, public :: link_target_waters_reactive_zone
        procedure, public :: link_target_waters_mineral_zone
    !> Check
        procedure, public :: check_new_reactive_zones
    !> Others
        procedure, public :: loop_read_tar_wat_init !> nombre muy malo
    !> Mixing
        procedure, public :: interfaz_comps_arch
        procedure, public :: interfaz_comps_vars
        procedure, public :: interfaz_esp_arch
        procedure, public :: compute_u_tilde_init
    end type
    
    
    interface

    subroutine interfaz_esp_arch(this,path,num_aq_comps,file_in,Delta_t,file_out)
        !> This subroutine solves reactive mixing iteration when there are no equilibrium reactions
        import chemistry_c
        implicit none
        class(chemistry_c) :: this !> chemistry object
        character(len=*), intent(in) :: path !> path for input and output files
        integer(kind=4), intent(in) :: num_aq_comps !> number of aqueous species (= number of aqueous components)
        character(len=*), intent(in) :: file_in !> name of file containing aqueous species concentrations after solving conservative transport
        real(kind=8), intent(in) :: Delta_t !> time step
        character(len=*), intent(in) :: file_out !> name of file containing aqueous species concentrations after solving reactive mixing
    end subroutine
        !subroutine solve_reactive_mixing_bis(this,root,mixing_ratios,mixing_waters_indices,F_mat,time_discr,&
        !    int_method_chem_reacts)
        !    import chemistry_c
        !    import real_array_c
        !    import int_array_c
        !    !import diag_matrix_c
        !    import time_discr_c
        !    implicit none
        !    class(chemistry_c) :: this
        !    character(len=*), intent(in) :: root
        !    !integer(kind=4), intent(in) :: unit
        !    class(real_array_c), intent(in) :: mixing_ratios
        !    class(int_array_c), intent(in) :: mixing_waters_indices
        !    !class(diag_matrix_c), intent(in) :: F_mat !> storage matrix
        !    real(kind=8), intent(in) :: F_mat(:) !> storage matrix (diagonal)
        !    class(time_discr_c), intent(in) :: time_discr !> time discretisation object
        !    integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method for chemical reactions
        !end subroutine
        
        subroutine solve_reactive_mixing_lump(this,root,mixing_ratios,mixing_waters_indices,mixing_waters_indices_dom,time_discr,&
                int_method_chem_reacts)
            import chemistry_c
            import real_array_c
            import int_array_c
            import time_discr_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root
            !integer(kind=4), intent(in) :: unit
            class(real_array_c), intent(in) :: mixing_ratios
            class(int_array_c), intent(in) :: mixing_waters_indices
            class(int_array_c), intent(in) :: mixing_waters_indices_dom
            !class(diag_matrix_c), intent(in) :: F_mat !> storage matrix
            !real(kind=8), intent(in) :: F_mat(:) !> storage matrix (diagonal)
            class(time_discr_c), intent(in) :: time_discr !> time discretisation object
            integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method for chemical reactions
        end subroutine
        
        subroutine solve_reactive_mixing_ideal_cons(this,root,mixing_ratios_conc,mixing_ratios_Rk_init,mixing_waters_indices,&
                mixing_waters_indices_dom,&
            time_discr,int_method_chem_reacts,mixing_ratios_Rk)
            import chemistry_c
            import real_array_c
            import int_array_c
            !import diag_matrix_c
            import time_discr_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root
            !integer(kind=4), intent(in) :: unit
            class(real_array_c), intent(in) :: mixing_ratios_conc
            class(real_array_c), intent(in) :: mixing_ratios_Rk_init
            class(int_array_c), intent(in) :: mixing_waters_indices
            class(int_array_c), intent(in) :: mixing_waters_indices_dom !> matrix that contains indices of domain target waters that mix with each target water
            !class(diag_matrix_c), intent(in) :: F_mat !> storage matrix
            !real(kind=8), intent(in) :: F_mat(:) !> storage matrix (diagonal)
            class(time_discr_c), intent(in) :: time_discr !> time discretisation object
            integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method for chemical reactions
            class(real_array_c), intent(inout) :: mixing_ratios_Rk
        end subroutine

        subroutine solve_reactive_mixing_cons(this,root,mixing_ratios_conc,mixing_ratios_Rk_init,mixing_waters_indices,&
            mixing_waters_indices_dom,time_discr,&
            int_method_chem_reacts,mixing_ratios_Rk)
            import chemistry_c
            import real_array_c
            import int_array_c
            !import diag_matrix_c
            import time_discr_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root
            !integer(kind=4), intent(in) :: unit
            class(real_array_c), intent(in) :: mixing_ratios_conc
            class(real_array_c), intent(in) :: mixing_ratios_Rk_init
            class(int_array_c), intent(in) :: mixing_waters_indices
            class(int_array_c), intent(in) :: mixing_waters_indices_dom
            !class(diag_matrix_c), intent(in) :: F_mat !> storage matrix
            !real(kind=8), intent(in) :: F_mat(:) !> storage matrix (diagonal)
            class(time_discr_c), intent(in) :: time_discr !> time discretisation object
            integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method for chemical reactions
            class(real_array_c), intent(inout) :: mixing_ratios_Rk
        end subroutine
            
        subroutine solve_reactive_mixing_ideal_lump(this,root,mixing_ratios,mixing_waters_indices,mixing_waters_indices_dom,&
                time_discr,&
            int_method_chem_reacts)
            import chemistry_c
            import real_array_c
            import int_array_c
            import time_discr_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root
            !integer(kind=4), intent(in) :: unit
            class(real_array_c), intent(in) :: mixing_ratios
            class(int_array_c), intent(in) :: mixing_waters_indices
            class(int_array_c), intent(in) :: mixing_waters_indices_dom
            !class(diag_matrix_c), intent(in) :: F_mat !> storage matrix
            !real(kind=8), intent(in) :: F_mat(:) !> storage matrix (diagonal)
            class(time_discr_c), intent(in) :: time_discr !> time discretisation object
            integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method for chemical reactions
        end subroutine
        
        
        subroutine solve_reactive_mixing_BCs_dep_t(this,root,unit,mixing_ratios,mixing_waters_indices,time_discr_tpt,&
            int_method_chem_reacts,spatial_discr_tpt,D,q,phi,anal_sol)
            import chemistry_c
            import spatial_discr_c
            import time_discr_c
            import real_array_c
            import int_array_c
            implicit none
        !> Arguments
            class(chemistry_c) :: this !> chemistry object
            character(len=*), intent(in) :: root
            integer(kind=4), intent(in) :: unit
            class(real_array_c), intent(in) :: mixing_ratios !> mixing ratios matrix
            class(int_array_c), intent(in) :: mixing_waters_indices !> matrix that contains indices of target waters that mix with each target water
            class(time_discr_c), intent(in) :: time_discr_tpt !> time discretisation object (used to solve transport)
            integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method for chemical reactions
            class(spatial_discr_c), intent(in) :: spatial_discr_tpt !> spatial discretisation object (used to solve transport)
            real(kind=8), intent(in) :: D
            real(kind=8), intent(in) :: q
            real(kind=8), intent(in) :: phi
            real(kind=8), external :: anal_sol
        end subroutine
        
        ! subroutine read_chemistry_PHREEQC(this,path_inp,path_DB,filename)
        !     import chemistry_c
        !     implicit none
        !     class(chemistry_c) :: this
        !     character(len=*), intent(in) :: path_inp
        !     character(len=*), intent(in) :: path_DB
        !     character(len=*), intent(in) :: filename
        ! end subroutine
        
        subroutine read_chemistry_CHEPROO(this,root,path_pb,path_DB,unit_chem_syst_file,unit_loc_chem_file)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root
            character(len=*), intent(in) :: path_pb
            character(len=*), intent(in) :: path_DB
            integer(kind=4), intent(in) :: unit_chem_syst_file
            !character(len=*), intent(in) :: chem_syst_file
            integer(kind=4), intent(in) :: unit_loc_chem_file
            !character(len=*), intent(in) :: loc_chem_file
            !integer(kind=4), intent(in) :: unit_tar_wat_file
            !character(len=*), intent(in) :: target_waters_dom_file
            !integer(kind=4), intent(in) :: unit_output_file
            !character(len=*), intent(in) :: output_file
        end subroutine
        
        subroutine read_chemistry(this,root,path_pb,path_DB)!,unit_chem_opts_file,unit_chem_syst_file,unit_loc_chem_file,&
            !unit_target_waters_dom_file,unit_output_file)
            import chemistry_c
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root !> root for input/output files
            character(len=*), intent(in) :: path_pb !> path to problem
            character(len=*), intent(in) :: path_DB !> path to database files
            !integer(kind=4), intent(in) :: unit_chem_opts_file
            !integer(kind=4), intent(in) :: unit_chem_syst_file
            !!character(len=*), intent(in) :: chem_syst_file
            !integer(kind=4), intent(in) :: unit_loc_chem_file
            !!character(len=*), intent(in) :: loc_chem_file
            !integer(kind=4), intent(in) :: unit_target_waters_dom_file
            !!character(len=*), intent(in) :: target_waters_dom_file
            !integer(kind=4), intent(in) :: unit_output_file
            !!character(len=*), intent(in) :: output_file
        end subroutine
                
        subroutine read_reactive_zones_Lagr(this,unit)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit
        end subroutine
        
        subroutine link_target_waters_dom_reactive_zone(this,i,tar_wat_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: i
            integer(kind=4), intent(out), allocatable :: tar_wat_indices(:) 
        end subroutine
        
        subroutine link_target_waters_reactive_zone(this,i,dom_indices,ext_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: i
            integer(kind=4), intent(out), allocatable :: dom_indices(:) 
            integer(kind=4), intent(out), allocatable :: ext_indices(:) 
        end subroutine
        
        subroutine link_target_waters_mineral_zone(this,i,dom_indices,ext_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: i
            integer(kind=4), intent(out), allocatable :: dom_indices(:) 
            integer(kind=4), intent(out), allocatable :: ext_indices(:) 
        end subroutine
        
        subroutine link_target_waters_target_solid(this,i,tw_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: i
            integer(kind=4), intent(out), allocatable :: tw_indices(:) 
        end subroutine
        
        subroutine link_target_waters_target_solids(this,tar_sol_indices,tar_wat_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: tar_sol_indices(:) !> target solid indices
            integer(kind=4), intent(inout), allocatable :: tar_wat_indices(:) 
        end subroutine
        
        subroutine link_target_waters_target_gases(this,tar_gas_indices,tar_wat_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: tar_gas_indices(:) !> target gas indices
            integer(kind=4), intent(inout), allocatable :: tar_wat_indices(:) 
        end subroutine
        
        subroutine link_target_solids_reactive_zone(this,i,tar_sol_indices)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: i
            integer(kind=4), intent(out), allocatable :: tar_sol_indices(:)
        end subroutine
        
        subroutine link_target_gases_reactive_zone(this,i,tar_gas_indices) 
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: i
            integer(kind=4), intent(out), allocatable :: tar_gas_indices(:)
        end subroutine
        
       

        
        subroutine initialise_target_waters_dom(this,initial_water_types)
            import chemistry_c
            import aqueous_chemistry_c
            import chem_system_c
            implicit none
            class(chemistry_c) :: this
            class(aqueous_chemistry_c), intent(in) :: initial_water_types(:)
        end subroutine
        
        subroutine read_target_waters_init(this,root,nsrz,ngrz)
            import chemistry_c
            import aqueous_chemistry_c
            import solid_chemistry_c
            import gas_chemistry_c
            import reactive_zone_c
            !import aq_phase_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root
            !integer(kind=4), intent(in) :: unit !> file
            !class(aqueous_chemistry_c), intent(in) :: init_water_types(:)
            !class(aqueous_chemistry_c), intent(in) :: bd_water_types(:)
            !type(aqueous_chemistry_c), intent(in) :: water_types(:)
            !type(solid_chemistry_c), intent(in) :: init_sol_types(:)
            !type(gas_chemistry_c), intent(in) :: init_gas_types(:)
            integer(kind=4), intent(in) :: nsrz !> number of solid reactive zones
            integer(kind=4), intent(in) :: ngrz !> number of gas reactive zones
            !type(reactive_zone_c), intent(in) :: react_zones(:) !> all possible reactive zones
            !integer(kind=4), intent(out) :: niter !> number of iterations
            !logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
            !class(aq_phase_c), intent(out), optional :: aq_phase_new
        end subroutine
        
        subroutine read_target_waters_dom_bis(this,unit,init_water_types,bd_water_types,init_sol_types,niter,CV_flag)
            import chemistry_c
            import aqueous_chemistry_c
            import solid_chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit !> file
            class(aqueous_chemistry_c), intent(in) :: init_water_types(:)
            class(aqueous_chemistry_c), intent(in) :: bd_water_types(:)
            class(solid_chemistry_c), intent(in) :: init_sol_types(:)
            !real(kind=8), intent(in) :: tolerance
            !real(kind=8), intent(in) :: rel_tolerance
            !real(kind=8), intent(in) :: control_factor
            !integer(kind=4), intent(in) :: niter_max
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine
       
        
        subroutine initialise_target_solids(this,n)
            import chemistry_c
            import solid_chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: n
        end subroutine
        
        subroutine initialise_ext_waters(this)
            import chemistry_c
            import chem_system_c
            implicit none
            class(chemistry_c) :: this
        end subroutine
        
       
        
        subroutine initialise_chemistry(this,path_DB,root,unit_chem_syst_file,unit_loc_chem_file,unit_target_waters_dom_file,&
            unit_output_file)
            import chemistry_c
            import real_array_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root
            character(len=*), intent(in) :: path_DB
            integer(kind=4), intent(in) :: unit_chem_syst_file
            !character(len=*), intent(in) :: chem_syst_file
            integer(kind=4), intent(in) :: unit_loc_chem_file
            !character(len=*), intent(in) :: loc_chem_file
            integer(kind=4), intent(in) :: unit_target_waters_dom_file
            !character(len=*), intent(in) :: target_waters_dom_file
            integer(kind=4), intent(in) :: unit_output_file
            !character(len=*), intent(in) :: output_file
        end subroutine
        
        subroutine write_chemistry(this,unit,Lagr_flag)
            import chemistry_c
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit
            logical, intent(in) :: Lagr_flag
        end subroutine
        
        subroutine check_new_reactive_zones(this,i,tolerance)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in)  :: i !> reactive zone index
            real(kind=8), intent(in) :: tolerance !> for concentration of non-flowing species
        end subroutine
        
       
        
        subroutine read_init_bd_wat_types_CHEPROO(this,unit,init_cat_exch_zones,&
            gas_chem)
            import chemistry_c
            import aqueous_chemistry_c
            import solid_chemistry_c
            import gas_chemistry_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit !> file
            type(solid_chemistry_c), intent(inout) :: init_cat_exch_zones(:)
            type(gas_chemistry_c), intent(in), optional :: gas_chem !> chapuza
        end subroutine
        
        subroutine read_init_min_zones_CHEPROO(this,unit,num_rz,surf_chem)
            import chemistry_c
            import solid_chemistry_c
            import reactive_zone_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit !> file
            !type(solid_chemistry_c), intent(out), allocatable :: init_min_zones(:)
            integer(kind=4), intent(out) :: num_rz !> number of mineral reactive zones
            type(solid_chemistry_c), intent(in), optional :: surf_chem
        end subroutine
        
        subroutine read_init_min_zones_CHEPROO_bis(this,unit,init_min_zones,reactive_zones)
            import chemistry_c
            import solid_chemistry_c
            import reactive_zone_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit !> file
            type(solid_chemistry_c), intent(out), allocatable :: init_min_zones(:)
            type(reactive_zone_c), intent(out), allocatable, optional :: reactive_zones(:)
        end subroutine
        
        subroutine read_init_cat_exch_zones_CHEPROO(this,unit,init_cat_exch_zones)
            import chemistry_c
            import solid_chemistry_c
            import reactive_zone_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit !> file
            type(solid_chemistry_c), intent(out), allocatable :: init_cat_exch_zones(:)
            !type(reactive_zone_c), intent(inout), allocatable, optional :: reactive_zones(:)
        end subroutine
        
        subroutine read_gas_bd_zones_CHEPROO(this,unit,gas_bd_zones,reactive_zones)
            import chemistry_c
            import gas_chemistry_c
            import reactive_zone_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit !> file
            type(gas_chemistry_c), intent(out), allocatable :: gas_bd_zones(:)
            type(reactive_zone_c), intent(out), allocatable, optional :: reactive_zones(:)
        end subroutine
        
        subroutine read_init_gas_zones_CHEPROO(this,unit,ngrz)
            import chemistry_c
            import gas_chemistry_c
            import reactive_zone_c
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: unit !> file
            !type(gas_chemistry_c), intent(out), allocatable :: gas_zones(:)
            !type(reactive_zone_c), intent(out), allocatable, optional :: reactive_zones(:)
            integer(kind=4), intent(out) :: ngrz !> number of gas reactive zones
        end subroutine

        subroutine interfaz_comps_arch(this,path,num_aq_comps,file_in,Delta_t,file_out)
            import chemistry_c
            class(chemistry_c) :: this !> chemistry object
            character(len=*), intent(in) :: path !> path for input and output files
            integer(kind=4), intent(in) :: num_aq_comps !> number of aqueous components
            character(len=*), intent(in) :: file_in !> name of file containing component concentrations after solving conservative transport
            !integer(kind=4), intent(in) :: unit_in !> file unit
            real(kind=8), intent(in) :: Delta_t !> time step
            character(len=*), intent(in) :: file_out !> name of file containing component concentrations after solving conservative transport
            !integer(kind=4), intent(in) :: unit_out !> file unit
        end subroutine
        
            subroutine interfaz_comps_vars(this,u_tilde,Delta_t,u_new)
            import chemistry_c
            class(chemistry_c) :: this
            real(kind=8), intent(in) :: u_tilde(:) !> concentrations after solving conservative transport
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(out) :: u_new(:) !> concentrations after solving reactive mixing
            end subroutine

            subroutine read_tar_sol(this,root,nsrz,ngrz)
            import chemistry_c
            !import solid_chemistry_c
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root !> root name of the file
            !type(solid_chemistry_c), intent(in) :: init_sol_types(:) !> initial solid zones
            integer(kind=4), intent(in) :: nsrz !> number of solid reactive zones
            integer(kind=4), intent(in) :: ngrz !> number of gas reactive zones
            end subroutine
            
            subroutine read_tar_gas(this,root,ngrz)
            import chemistry_c
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root !> root name of the file
            integer(kind=4), intent(in) :: ngrz !> number of gas reactive zones
            end subroutine

    end interface
    
    contains
    
        subroutine set_lump_flag(this,lump_flag)
        class(chemistry_c) :: this
        logical, intent(in) :: lump_flag
        this%lump_flag=lump_flag
        end subroutine
        
        subroutine set_cons_opt(this,cons_opt)
        class(chemistry_c) :: this
        integer(kind=4), intent(in) :: cons_opt
        if (cons_opt<1 .or. cons_opt>2) then
            error stop "Chemistry attribute 'cons_opt' must be 1 or 2"
        else
            this%cons_opt=cons_opt
        end if
        end subroutine
            
        subroutine set_rk_avg_opt(this,rk_avg_opt)
        class(chemistry_c) :: this
        integer(kind=4), intent(in) :: rk_avg_opt
        if (rk_avg_opt<1 .or. rk_avg_opt>2) then
            error stop "Chemistry attribute 'rk_avg_opt' must be 1 or 2"
        else
            this%rk_avg_opt=rk_avg_opt
        end if
        end subroutine

        subroutine set_rk_down_opt(this,rk_down_opt)
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: rk_down_opt
            if (rk_down_opt<1 .or. rk_down_opt>4) then
                error stop "Chemistry attribute 'rk_down_opt' must be 1, 2, 3 or 4"
            else
                this%rk_down_opt=rk_down_opt
            end if
            end subroutine
        
        subroutine set_read_opt(this,option)
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: option
            if (option<1 .or. option>3) then
                error stop "Chemistry input option not implemented yet"
            else if (option>1) then
                error stop "Chemistry input option not fully implemented yet"
            else
                this%read_opt=option
            end if
        end subroutine
        
        subroutine set_Jac_opt(this,Jac_opt)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: Jac_opt
            if (Jac_opt<0 .or. Jac_opt>1) then
                error stop "Chemistry attribute 'Jac_opt' must be 0 or 1"
            else
                this%Jac_opt=Jac_opt
            end if
        end subroutine
        
        subroutine set_chem_syst(this,chem_syst_obj)
            implicit none
            class(chemistry_c) :: this
            type(chem_system_c), intent(in) :: chem_syst_obj
            this%chem_syst=chem_syst_obj
        end subroutine
        
        subroutine set_num_wat_types(this,num_wat_types)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: num_wat_types
            if (num_wat_types<1) then
                error stop "Chemistry attribute 'num_wat_types' must be greater than 0"
            else
                this%num_wat_types=num_wat_types
            end if
        end subroutine
        
        subroutine set_num_tar_wat_dom(this,num_tar_wat_dom)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: num_tar_wat_dom
            if (present(num_tar_wat_dom)) then
                this%num_target_waters_dom=num_tar_wat_dom
            else
                this%num_target_waters_dom=this%num_target_waters-this%num_ext_waters
            end if
        end subroutine
        
        subroutine set_num_target_solids(this,num_target_solids)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: num_target_solids
            this%num_target_solids=num_target_solids
        end subroutine
        
        subroutine set_num_ext_waters(this,num_ext_waters)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: num_ext_waters
            !if (present(num_ext_waters)) then
                this%num_ext_waters=num_ext_waters
            !else
            !    this%num_ext_waters=this%num_target_waters_dom
            !end if
        end subroutine
        
        subroutine set_num_bd_waters(this,num_bd_waters)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: num_bd_waters
            !if (present(num_bd_waters)) then
                this%num_bd_waters=num_bd_waters
            !else
            !    this%num_bd_waters=this%num_target_waters_dom
            !end if
        end subroutine
        
        subroutine set_num_rech_waters(this,num_rech_waters)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in) :: num_rech_waters
            this%num_rech_waters=num_rech_waters
        end subroutine
        
        subroutine set_target_solids(this,target_solids)
            implicit none
            class(chemistry_c) :: this
            type(solid_chemistry_c), intent(in) :: target_solids(:)
            this%target_solids=target_solids
            this%num_target_solids=size(target_solids)
        end subroutine
        
        subroutine set_target_gases(this,target_gases)
            implicit none
            class(chemistry_c) :: this
            type(gas_chemistry_c), intent(in) :: target_gases(:)
            this%target_gases=target_gases
            this%num_target_gases=size(target_gases)
        end subroutine
        
       subroutine allocate_target_waters(this,num_tar_wat)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: num_tar_wat
            if (present(num_tar_wat)) then
                this%num_target_waters=num_tar_wat
            end if
            if (allocated(this%target_waters)) then
                deallocate(this%target_waters)
            end if
            allocate(this%target_waters(this%num_target_waters))
        end subroutine
        
        subroutine allocate_dom_tar_wat_indices(this,num_tar_wat_dom)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: num_tar_wat_dom
            if (present(num_tar_wat_dom)) then
                this%num_target_waters_dom=num_tar_wat_dom
            end if
            if (allocated(this%dom_tar_wat_indices)) then
                deallocate(this%dom_tar_wat_indices)
            end if
            allocate(this%dom_tar_wat_indices(this%num_target_waters_dom))
        end subroutine
        
        subroutine allocate_target_solids(this,n)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: n
            if (present(n)) then
                this%num_target_solids=n
            end if
            if (allocated(this%target_solids)) then
                deallocate(this%target_solids,this%target_solids_init)
            end if
            allocate(this%target_solids(this%num_target_solids),this%target_solids_init(this%num_target_solids))
        end subroutine
        
        
        subroutine allocate_target_gases(this,n)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: n
            if (present(n)) then
                this%num_target_gases=n
            end if
            allocate(this%target_gases(this%num_target_gases),this%target_gases_init(this%num_target_gases))
        end subroutine
        
        subroutine allocate_rech_waters_indices(this,num_rech_waters)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: num_rech_waters
            if (present(num_rech_waters)) then
                this%num_rech_waters=num_rech_waters
            end if
            allocate(this%rech_waters_indices(this%num_rech_waters))
        end subroutine
        
        subroutine allocate_bd_waters_indices(this,num_bd_waters)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: num_bd_waters
            if (present(num_bd_waters)) then
                this%num_bd_waters=num_bd_waters
            end if
            allocate(this%bd_waters_indices(this%num_bd_waters))
        end subroutine
        
       
        
        subroutine set_target_waters(this,target_waters)
            implicit none
            class(chemistry_c) :: this
            type(aqueous_chemistry_c), intent(in) :: target_waters(:)
            this%target_waters=target_waters
            this%num_target_waters=size(target_waters)
        end subroutine
        
        !subroutine set_target_waters_dom(this,target_waters_dom)
        !    implicit none
        !    class(chemistry_c) :: this
        !    class(aqueous_chemistry_c), intent(in) :: target_waters_dom(:)
        !    this%target_waters()=target_waters_dom
        !    this%num_target_waters_dom=size(target_waters_dom)
        !end subroutine
        
       subroutine allocate_reactive_zones(this,num_reactive_zones)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: num_reactive_zones
            if (present(num_reactive_zones)) then
                this%num_reactive_zones=num_reactive_zones
            end if
            if (allocated(this%reactive_zones)) then
                deallocate(this%reactive_zones)
            end if
            allocate(this%reactive_zones(this%num_reactive_zones))
       end subroutine
       
       subroutine allocate_mineral_zones(this,num_mineral_zones)
            implicit none
            class(chemistry_c) :: this
            integer(kind=4), intent(in), optional :: num_mineral_zones
            if (present(num_mineral_zones)) then
                this%num_mineral_zones=num_mineral_zones
            end if
            if (allocated(this%mineral_zones)) then
                deallocate(this%mineral_zones)
            end if
            allocate(this%mineral_zones(this%num_mineral_zones))
       end subroutine
       
       subroutine set_reactive_zones(this,reactive_zones)
            implicit none
            class(chemistry_c) :: this
            type(reactive_zone_c), intent(in), optional :: reactive_zones(:)
            integer(kind=4) :: i,j,l
            logical :: flag
            integer(kind=4), allocatable :: rz_indices(:)
            if (present(reactive_zones)) then
                this%reactive_zones=reactive_zones
                this%num_reactive_zones=size(reactive_zones)
            else if (this%num_target_solids>0) then
                i=1
                j=2
                this%num_reactive_zones=1
                allocate(rz_indices(this%num_target_solids))
                rz_indices=0
                do
                        call compare_react_zones(this%target_solids(i)%reactive_zone,this%target_solids(j)%reactive_zone,flag)
                        if (flag.eqv..true.) then
                            if (i<this%num_target_solids-1) then
                                i=i+1
                                j=i+1
                            else
                                exit
                            end if
                        else if (j<this%num_target_solids) then
                            j=j+1
                        else if (i<this%num_target_solids-1) then
                            this%num_reactive_zones=this%num_reactive_zones+1
                            rz_indices(i)=1
                            i=i+1
                            j=i+1
                        else
                            exit
                        end if
                end do
                rz_indices(j)=1
                call this%allocate_reactive_zones()
                l=1
                do i=1,this%num_target_solids
                    if (rz_indices(i)==1) then
                        call this%reactive_zones(l)%assign_react_zone(this%target_solids(i)%reactive_zone)
                        l=l+1
                    end if
                end do
            else if (this%num_target_gases>0) then
                i=1
                j=2
                this%num_reactive_zones=1
                allocate(rz_indices(this%num_target_gases))
                rz_indices=0
                do
                        call compare_react_zones(this%target_gases(i)%reactive_zone,this%target_gases(j)%reactive_zone,flag)
                        if (flag.eqv..true.) then
                            if (i<this%num_target_gases-1) then
                                i=i+1
                                j=i+1
                            else
                                exit
                            end if
                        else if (j<this%num_target_gases) then
                            j=j+1
                        else if (i<this%num_target_gases-1) then
                            this%num_reactive_zones=this%num_reactive_zones+1
                            rz_indices(i)=1
                            i=i+1
                            j=i+1
                        else
                            exit
                        end if
                end do
                rz_indices(j)=1
                call this%allocate_reactive_zones()
                l=1
                do i=1,this%num_target_gases
                    if (rz_indices(i)==1) then
                        this%reactive_zones(l)=this%target_gases(i)%reactive_zone
                        l=l+1
                    end if
                end do
            else
                call this%allocate_reactive_zones(0)
            end if
       end subroutine
       

            subroutine loop_read_tar_wat_init(this,flag_ext,water_types,init_sol_types,init_gas_types,nsrz,ngrz,tar_wat_ind,wtype,&
            istype,igzn,aux_istype,aux_igzn,solid_chem_def)
            class(chemistry_c) :: this !> chemistry object
            logical, intent(in) :: flag_ext !> flag to indicate if target water is external
            type(aqueous_chemistry_c), intent(in) :: water_types(:) !> initial water types
            type(solid_chemistry_c), intent(in) :: init_sol_types(:) !> initial solid types
            type(gas_chemistry_c), intent(in) :: init_gas_types(:) !> initial gas types
            !type(reactive_zone_c), intent(in) :: react_zone !> default reactive zone object
            !type(mineral_zone_c), intent(in) :: min_zone !> default mineral zone object
            integer(kind=4), intent(in) :: nsrz !> number of solid reactive zones
            integer(kind=4), intent(in) :: ngrz !> number of gas reactive zones
            integer(kind=4), intent(in) :: tar_wat_ind !> target water index
            integer(kind=4), intent(in) :: wtype !> water type index
            integer(kind=4), intent(in) :: istype !> solid type index
            integer(kind=4), intent(in) :: igzn !> gas zone index
            integer(kind=4), intent(inout) :: aux_istype !> auxiliary solid type index
            integer(kind=4), intent(inout) :: aux_igzn !> auxiliary gas zone index
            type(solid_chemistry_c), intent(in) :: solid_chem_def !> default solid chemistry object

            logical :: flag_Se !> flag to swap indices of species
            integer(kind=4), allocatable :: swap(:),aux_swap(:) !> indices of species to swap
            real(kind=8), allocatable :: rk(:) !> NO NECESARIA
            type(solid_chemistry_c), target :: aux_solid_chem !> auxiliary solid chemistry object
            allocate(swap(2),aux_swap(2)) !> chapuza    
            
            this%target_waters(tar_wat_ind)=water_types(wtype) !> we set the target water
            if (istype>0) then !> we have a solid zone
                this%target_solids(tar_wat_ind)=init_sol_types(istype) !> we set the target solid
                if (igzn>0) then !> we have a gas zone
                    call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngrz+istype)) !> we set the reactive zone of the target solid
                    this%target_gases(tar_wat_ind)=init_gas_types(igzn) !> we set the target gas
                    call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngrz+istype)) !> we set the reactive zone of the target gas
                    call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind)) !> we set the gas chemistry pointer
                else
                    call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngrz+istype)) !> we set the reactive zone of the target solid
                end if
                !call this%target_waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(tar_wat_ind)) !> we set the solid chemistry pointer
                call this%target_waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(1)) !> chapuza para probar Lagrangiano
            else
                if (igzn>0) then !> we assume we only have one gas reactive zone
                    this%target_gases(tar_wat_ind)=init_gas_types(igzn)
                    call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngrz)) !> chapuza
                    call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
                    call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngrz)) !> chapuza
                    call this%target_solids(tar_wat_ind)%set_mineral_zone(solid_chem_def%mineral_zone) !> default mineral zone
                else
                    this%target_solids(tar_wat_ind)=solid_chem_def !> we set the default solid chemistry object
                end if
                call this%target_waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(tar_wat_ind)) !> we set the solid chemistry pointer
            end if
            !> we check if there is a new reactive zone
            if (aux_istype==0 .or. aux_istype/=istype .or. aux_igzn/=igzn) then !> we assume target waters are grouped by their reactive zones
                call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_speciation_alg_dimensions(.true.) !< esto creo que no es necesario
                call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_ind_eq_reacts() !> chapuza
                call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_stoich_mat_react_zone() !> chapuza
                call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_ind_gases_stoich_mat() !> chapuza
                call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_ind_mins_stoich_mat() !> chapuza
                call this%target_waters(tar_wat_ind)%set_ind_species()
                if (associated(this%target_waters(tar_wat_ind)%gas_chemistry)) then
                    call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
                        flag_Se,swap,this%target_waters(tar_wat_ind)%gas_chemistry%activities)
                else
                    call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
                        flag_Se,swap)
                end if
                if (flag_Se.eqv..true.) then !> we swap indices of species
                    aux_swap(1)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))
                    aux_swap(2)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))
                    this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))=aux_swap(2)
                    this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))=aux_swap(1)
                end if
            else if (aux_istype>0 .or. aux_igzn>0) then !> indices remain the same because reactive zone is the same
                this%target_waters(tar_wat_ind)%ind_var_act_species=this%target_waters(tar_wat_ind-1)%ind_var_act_species
                !this%target_waters(tar_wat_ind)%ind_sec_species=this%target_waters(tar_wat_ind-1)%ind_sec_species
            end if
        !> Chapuza
            !if (associated(this%target_waters(tar_wat_ind)%solid_chemistry%mineral_zone)) then
                if (this%target_waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin<&
                    this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin) then
                    call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod(&
                        this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                        this%target_waters(tar_wat_ind)%solid_chemistry%mineral_zone%ind_min_chem_syst(1:&
                        this%target_waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin))
                else
                    call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod()
                end if
            !end if
            call this%target_waters(tar_wat_ind)%allocate_reaction_rates()
            call this%target_waters(tar_wat_ind)%set_indices_rk()
            if (flag_ext .eqv. .false.) then !> external waters do not have kinetic reaction rates
                allocate(rk(this%target_waters(tar_wat_ind)%indices_rk%num_cols))
                call this%target_waters(tar_wat_ind)%compute_rk(rk)
                !!call this%target_waters(tar_wat_ind)%update_rk_old()
            end if
            deallocate(swap,aux_swap)
            end subroutine

            subroutine read_tar_wat_line(this,flag_ext,iszn,igzn,tar_wat_ind,wtype,tar_sol_ind,tar_gas_ind,aux_iszn,aux_igzn)
            class(chemistry_c) :: this !> chemistry object
            logical, intent(in) :: flag_ext !> flag to indicate if target water is external
            integer(kind=4), intent(in) :: iszn !> index of solid zone
            integer(kind=4), intent(in) :: igzn !> index of gas reactive zone
            integer(kind=4), intent(in) :: tar_wat_ind !> target water index
            integer(kind=4), intent(in) :: wtype !> water type index
            integer(kind=4), intent(in) :: tar_sol_ind !> target solid index
            integer(kind=4), intent(in) :: tar_gas_ind !> target gas index
            integer(kind=4), intent(inout) :: aux_iszn !> auxiliary solid zone index
            integer(kind=4), intent(inout) :: aux_igzn !> auxiliary gas zone index
            !type(solid_chemistry_c), intent(in) :: solid_chem_def !> default solid chemistry object

            logical :: flag_Se !> flag to swap indices of species
            integer(kind=4), allocatable :: swap(:),aux_swap(:) !> indices of species to swap
            real(kind=8), allocatable :: rk(:) !> NO NECESARIA
            type(solid_chemistry_c), target :: aux_solid_chem !> auxiliary solid chemistry object
            allocate(swap(2),aux_swap(2)) !> chapuza

            if (tar_wat_ind<1 .or. tar_wat_ind>this%num_target_waters) then
                error stop "Target water index out of bounds"
            end if
            
            this%target_waters(tar_wat_ind)=this%wat_types(wtype) !> we set the target water
            call this%target_waters(tar_wat_ind)%set_id(tar_wat_ind) !> we set the target water id
            if (tar_sol_ind>0) then !> we have a target solid
                !this%target_solids(tar_wat_ind)=init_sol_types(istype) !> we set the target solid
                if (tar_gas_ind>0) then !> we have a target gas
                    !call this%target_solids(tar_sol_ind)%set_reactive_zone(this%reactive_zones(ngrz+istype)) !> we set the reactive zone of the target solid
                    !this%target_gases(tar_wat_ind)=this%gas_zones(igzn) !> we set the target gas
                    !call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngrz+istype)) !> we set the reactive zone of the target gas
                    call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_gas_ind)) !> we set the gas chemistry pointer
                else
                    !nullify(this%target_waters(tar_wat_ind)%gas_chemistry) !> (chapuza) we nullify the gas chemistry pointer
                end if
                !call this%target_waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(tar_wat_ind)) !> we set the solid chemistry pointer
                call this%target_waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(tar_sol_ind)) !> we set the solid chemistry pointer
            else !> we do not have a solid zone
                if (tar_gas_ind>0) then !> we assume we only have gases
                    !this%target_gases(tar_wat_ind)=this%gas_zones(igzn)
                    !call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngrz)) !> chapuza
                    call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_gas_ind))
                    !call this%target_solids(1)%set_reactive_zone(this%reactive_zones(ngrz)) !> chapuza
                    !call this%target_solids(1)%set_mineral_zone(solid_chem_def%mineral_zone) !> default mineral zone
                else
                    !this%target_solids(1)=solid_chem_def !> we set the default solid chemistry object
                end if
                if (flag_ext .and. tar_wat_ind<2) then !> chapuza
                    call this%target_waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(1)) !> we set the default solid chemistry pointer (chapuza)
                else if (flag_ext .or. tar_wat_ind>this%num_target_solids) then
                    call this%target_waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(this%num_target_solids)) !> we set the solid chemistry pointer
                else if (tar_wat_ind<3 .and. this%num_bd_waters>0) then !> we assume target waters are grouped by their solid chemistry
                    call this%target_waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(1)) !> we set the solid chemistry pointer
                else if (this%num_bd_waters>0) then
                    call this%target_waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(tar_wat_ind-1)) !> we set the default solid chemistry pointer (chapuza)
                end if
            end if
            !> we check if there is a new reactive zone
            !if (tar_sol_ind>0) then !> we have a solid reactive zone
                ! if (aux_iszn>0 .and. aux_iszn==this%target_solids(tar_sol_ind)%id) then
                !     continue !> we have the same reactive zone as previous target solid
                ! else !> we have a new solid reactive zone
                !     call this%target_waters(tar_wat_ind)%set_ind_species() !> we set the indices of species
                !     if (associated(this%target_waters(tar_wat_ind)%gas_chemistry)) then
                !         call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
                !             flag_Se,swap,this%target_waters(tar_wat_ind)%gas_chemistry%activities)
                !     else
                !         call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
                !             flag_Se,swap)
                !     end if
                !     if (flag_Se.eqv..true.) then !> we swap indices of species
                !         aux_swap(1)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))
                !         aux_swap(2)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))
                !         this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))=aux_swap(2)
                !         this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))=aux_swap(1)
                !     end if
                ! end if
                if (aux_iszn==0 .or. aux_iszn/=iszn .or. aux_igzn/=igzn) then !> we have a new reactive zone
                    call this%target_waters(tar_wat_ind)%set_ind_species()
                    if (associated(this%target_waters(tar_wat_ind)%gas_chemistry)) then
                        call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
                            flag_Se,swap,this%target_waters(tar_wat_ind)%gas_chemistry%activities)
                    else
                        call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
                            flag_Se,swap)
                    end if
                    if (flag_Se.eqv..true.) then !> we swap indices of species
                        aux_swap(1)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))
                        aux_swap(2)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))
                        this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))=aux_swap(2)
                        this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))=aux_swap(1)
                    end if
            ! else if (tar_gas_ind>0) then !> we only have a gas reactive zone
            !     if (aux_igzn>0 .and. aux_igzn==this%target_gases(tar_gas_ind)%id) then
            !         continue !> we have the same reactive zone as previous target gas
            !     else !> we have a new gas reactive zone
            !         call this%target_waters(tar_wat_ind)%set_ind_species() !> we set the indices of species
            !         if (associated(this%target_waters(tar_wat_ind)%gas_chemistry)) then
            !             call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
            !                 flag_Se,swap,this%target_waters(tar_wat_ind)%gas_chemistry%activities)
            !         else
            !             call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
            !                 flag_Se,swap)
            !         end if
            !         if (flag_Se.eqv..true.) then !> we swap indices of species
            !             aux_swap(1)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))
            !             aux_swap(2)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))
            !             this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))=aux_swap(2)
            !             this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))=aux_swap(1)
            !         end if
            !     end if
                else !> indices remain the same because reactive zone is the same
                    this%target_waters(tar_wat_ind)%ind_var_act_species=this%target_waters(tar_wat_ind-1)%ind_var_act_species
                    !this%target_waters(tar_wat_ind)%ind_sec_species=this%target_waters(tar_wat_ind-1)%ind_sec_species
                end if
            !end if
        !> Chapuza
            !if (associated(this%target_waters(tar_wat_ind)%solid_chemistry%mineral_zone)) then
                if (this%target_waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin<&
                    this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin) then
                    call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod(&
                        this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                        this%target_waters(tar_wat_ind)%solid_chemistry%mineral_zone%ind_min_chem_syst(1:&
                        this%target_waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin))
                else
                    call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod()
                end if
            !end if
            call this%target_waters(tar_wat_ind)%allocate_reaction_rates() !> we allocate reaction rates
            call this%target_waters(tar_wat_ind)%set_indices_rk() !> we set the indices for reaction rates
            call this%target_waters(tar_wat_ind)%set_volume() !> we set the volume of the target water
            if (flag_ext .eqv. .false.) then !> external waters do not have kinetic reaction rates
                allocate(rk(this%target_waters(tar_wat_ind)%indices_rk%num_cols))
                call this%target_waters(tar_wat_ind)%compute_rk(rk)
                !!call this%target_waters(tar_wat_ind)%update_rk_old()
                !!call this%target_waters(tar_wat_ind)%solid_chemistry%update_rk_old()
                deallocate(rk)
            end if
            if (tar_sol_ind>0) then
                aux_iszn=this%target_solids(tar_sol_ind)%id
            end if
            if (tar_gas_ind>0) then
                aux_igzn=this%target_gases(tar_gas_ind)%id
            end if
            deallocate(swap,aux_swap)
        end subroutine


        subroutine read_chem_opts(this,root,unit)
        !> read chemical options from file
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root !> root of the file name
            integer(kind=4), intent(in) :: unit !> file

            integer(kind=4) :: opt !> chemical option
            logical :: flag !> lumping flag
            character(len=256) :: label !> label to read
        !> We open file with chemical options
            open(unit,file=root//'_chem_opts.dat',status='old',action='read')
        do 
            read(unit,*) label !> read label
            if (label=='end') then
                exit !> end of file
            else if (label=='CHEMICAL OPTIONS') then
                read(unit,*) opt !> read input option
                call this%set_read_opt(opt) !> set read option
                read(unit,*) opt !> read label
                call this%set_Jac_opt(opt) !> set Jacobian option
                read(unit,*) flag !> read lumping flag
                call this%set_lump_flag(flag) !> set lumping flag
                read(unit,*) opt !> read consistent WMA option
                call this%set_cons_opt(opt) !> set consistent WMA option
                read(unit,*) opt !> read rk downstream waters option
                call this%set_rk_down_opt(opt) !> set rk downstream waters option
                read(unit,*) opt !> read rk average option
                call this%set_rk_avg_opt(opt) !> set rk average option
            else 
                continue !> continue reading file
            end if
        end do
        close(unit)
        end subroutine
        
        subroutine set_num_tar_wat(this,num_tar_wat)
        implicit none
        class(chemistry_c) :: this
        integer(kind=4), intent(in) :: num_tar_wat
        this%num_target_waters=num_tar_wat
        end subroutine
        
        subroutine set_target_waters_wat_types(this,ind_wat_types)
        implicit none
        class(chemistry_c) :: this
        integer(kind=4), intent(in) :: ind_wat_types(:) !> indices of water types (dim=n� target waters)
        integer(kind=4) :: i !> loop index
        do i=1,this%num_target_waters
            this%target_waters(i)=this%wat_types(ind_wat_types(i))
        end do
        end subroutine
        
        subroutine allocate_wat_types(this,num_wat_types)
        implicit none
        class(chemistry_c) :: this
        integer(kind=4), intent(in) :: num_wat_types !> number of water types
        if (allocated(this%wat_types)) deallocate(this%wat_types)
        if (num_wat_types<1) then
            error stop "Number of water types must be greater than 0"
        end if
        call this%set_num_wat_types(num_wat_types) !> we set the number of water types
        allocate(this%wat_types(this%num_wat_types))
        !this%num_target_waters=n
        end subroutine
        
        subroutine set_target_waters_target_solids(this,ind_tar_sol)
        implicit none
        class(chemistry_c) :: this
        integer(kind=4), intent(in) :: ind_tar_sol(:) !> indices of target solids (dim=n� target waters)
        
        integer(kind=4) :: i !> loop index
        
        do i=1,this%num_target_waters
            if (ind_tar_sol(i)>0) then
                call this%target_waters(i)%set_solid_chemistry(this%target_solids(ind_tar_sol(i)))
            else
                call this%target_waters(i)%set_solid_chemistry(this%target_solids(1)) !> we set the default solid chemistry object
            end if
        end do
        end subroutine
        
        subroutine set_target_solids_materials(this,ind_materials)
        implicit none
        class(chemistry_c) :: this
        integer(kind=4), intent(in) :: ind_materials(:) !> indices of materials (dim=n� target solids)

        integer(kind=4) :: i !> loop index

        do i=1,this%num_target_solids
            if (ind_materials(i)>0) then
                call this%target_waters(i)%set_solid_chemistry(this%target_solids(ind_materials(i)))
            else
                call this%target_waters(i)%set_solid_chemistry(this%target_solids(1)) !> we set the default solid chemistry object
            end if
        end do
        end subroutine
        
        function get_num_aq_comps_chem_syst(this) result(num_aq_comps)
        !> Gets the number of aqueous components in the chemical system
        implicit none
        class(chemistry_c) :: this !> chemistry object
        integer(kind=4) :: num_aq_comps !> number of aqueous components in the chemical system
        num_aq_comps=this%chem_syst%speciation_alg%num_aq_prim_species !> 
        end function

        function get_num_aq_comps_dom(this) result(num_aq_comps)
        !> Gets the number of aqueous components in the domain
        !! We assume all target waters in the domain have the same number of aqueous components
        implicit none
        class(chemistry_c) :: this !> chemistry object
        integer(kind=4) :: num_aq_comps !> number of aqueous components in the domain
        num_aq_comps=this%target_waters(this%dom_tar_wat_indices(1)&
            )%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species !> we assume all target waters in the domain have the same number of aqueous components
        end function

        function get_num_aq_var_act_species_dom(this) result(n_nc_aq)
        !> Gets number of aqueous variable activity species in the domain
        !! We assume all target waters in the domain have the same number of aqueous variable activity species
        implicit none
        class(chemistry_c) :: this !> chemistry object
        integer(kind=4) :: n_nc_aq !> number of aqueous variable activity species in the domain
        n_nc_aq=this%target_waters(this%dom_tar_wat_indices(1))%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species !> we assume all target waters have the same number of variable activity species
        end function
        
        function get_num_wat_types(this) result(num_wat_types)
        implicit none
        class(chemistry_c) :: this
        integer(kind=4) :: num_wat_types !> number of water types
        num_wat_types=this%num_wat_types
        end function
        
        function get_conc_comp_wat_types(this) result(conc_comp_wat_types)
        !> Gets the concentrations of aqueous components of water types
        !! We assume all water types have the same number of aqueous components
        implicit none
        class(chemistry_c) :: this
        real(kind=8), allocatable :: conc_comp_wat_types(:,:) !> concentrations of aqueous components of water types

        integer(kind=4) :: num_aq_comps !> number of aqueous components in the chemical system
        integer(kind=4) :: i !> loop index

        num_aq_comps=this%get_num_aq_comps_chem_syst() !> we get the number of aqueous components in the chemical system
        allocate(conc_comp_wat_types(num_aq_comps,this%num_wat_types))
        do i=1,this%num_wat_types
            conc_comp_wat_types(:,i)=this%wat_types(i)%get_u_aq() !> we get the aqueous components of each water type
        end do
        end function

        subroutine write_conc_comp_wat_types(this,path,filename)
        !> Writes the concentrations of aqueous components of water types to a file
        implicit none
        !> Arguments
        class(chemistry_c) :: this !> chemistry object
        character(len=*), intent(in) :: path !> path to the file
        character(len=*), intent(in) :: filename !> file name
        !> Local variables
        real(kind=8), allocatable :: u_aq_init(:) !> aqueous components of water types
        integer(kind=4) :: i,j !> loop indices
        !> Body of the subroutine
        open(unit=10, file=path//filename, status='unknown', action='write', form='formatted')
        write(10,*) 'Concentrations of aqueous components of water types:'
        do i=1,this%num_wat_types
            u_aq_init=this%wat_types(i)%get_u_aq() !> we get the aqueous components of each water type
            write(10,"(/,2x,A15,/)") trim(this%wat_types(i)%name)//':' !> we write the name of the water type
            do j=1,this%wat_types(i)%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                write(10,"(10x,ES15.5)") u_aq_init(j) !> we write the aqueous components of each water type
            end do
            deallocate(u_aq_init)
        end do
        close(10)
        end subroutine
        
        subroutine set_target_solids_mesh(this,mesh)
        class(chemistry_c) :: this !> chemistry object
        class(spatial_discr_c), intent(in) :: mesh !> mesh object ( we assume it has targets defined )
        integer(kind=4) :: i !> loop index
        do i=1,this%num_target_solids
            call this%target_solids(i)%set_target(mesh%targets(i)) !> we set the target associated to each target solid
            call this%target_solids(i)%set_id(mesh%targets(i)%id) !> we set the same ID as the target
        end do
        end subroutine
        
        subroutine allocate_gas_zones(this,num_gas_zones)
        implicit none
        class(chemistry_c) :: this
        integer(kind=4), intent(in), optional :: num_gas_zones
        if (present(num_gas_zones)) then
            this%num_gas_zones=num_gas_zones
        end if
        if (allocated(this%gas_zones)) deallocate(this%gas_zones)
        allocate(this%gas_zones(this%num_gas_zones))!,this%target_gases_init(this%num_target_gases))
        end subroutine

        subroutine allocate_target_waters_init(this,num_target_waters_init)
        implicit none
        class(chemistry_c) :: this
        integer(kind=4), intent(in), optional :: num_target_waters_init
        if (present(num_target_waters_init)) then
            this%num_target_waters_init=num_target_waters_init
        end if
        if (allocated(this%target_waters_init)) deallocate(this%target_waters_init)
        allocate(this%target_waters_init(this%num_target_waters_init))
        end subroutine
        
        function compute_u_tilde_init(this,mixing_waters_indices,mixing_ratios) result(u_tilde_init)
        ! Variables
        class(chemistry_c) :: this
        !type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> waters that mix with each other
        !real(kind=8), intent(in) :: mixing_ratios_conc(:)  !> same order as "mix_waters"
        class(int_array_c), intent(in) :: mixing_waters_indices !> matrix that contains indices of target waters that mix with each target water
        class(real_array_c), intent(in) :: mixing_ratios !> mixing ratios matrix for concentrations
        real(kind=8), allocatable :: u_tilde_init(:,:) !> u_tilde of initial target waters
        
        integer(kind=4) :: i !> loop index
        real(kind=8), allocatable :: c_tilde(:) !> c_tilde of mixing waters
        ! Body of compute_u_tilde_init
        allocate(u_tilde_init(this%get_num_aq_comps_chem_syst(),this%num_target_waters_dom))
        do i=1,this%num_target_waters_dom !> we compute c_tilde for each domain target water
            allocate(c_tilde(this%target_waters_init(this%dom_tar_wat_indices(i)&
                )%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species)) !> we allocate c_tilde
            call compute_c_tilde(this%target_waters_init(mixing_waters_indices%cols(i)%col_1(&
                1:mixing_waters_indices%cols(i)%dim-2)),&
                mixing_ratios%cols(i)%col_1,c_tilde) !> we compute c_tilde for each domain target water
            u_tilde_init(:,i)=this%target_waters_init(this%dom_tar_wat_indices(i))%compute_u_tilde(c_tilde) !> we compute u_tilde from c_tilde
            deallocate(c_tilde)
        end do
        end function compute_u_tilde_init
        
        subroutine write_u_tilde_init(this,path,filename,u_tilde_init)
        class(chemistry_c) :: this !> chemistry object
        character(len=*), intent(in) :: path !> path to write the file
        character(len=*), intent(in) :: filename !> name of the file
        real(kind=8), intent(in) :: u_tilde_init(:,:) !> u_tilde of initial target waters
        
        integer(kind=4) :: i !> loop index
        
        open(unit=20,file=path//filename,status='replace',action='write')
        !call this%target_waters_init(1)%compute_u_tilde_init
        do i=1,size(u_tilde_init,1)
            write(20,"(*(ES15.5))") u_tilde_init(i,:)
        end do
        close(20)
        end subroutine
        
        subroutine set_num_materials(this,num_materials)
        class(chemistry_c) :: this !> chemistry object
        integer(kind=4), intent(in) :: num_materials !> number of materials
        if (num_materials<0) then
            error stop "Number of materials must be non-negative"
        end if
        this%num_materials=num_materials
        end subroutine
        
        subroutine allocate_materials(this)
        class(chemistry_c) :: this !> chemistry object
        allocate(this%materials(this%num_materials)) !> we allocate materials
        end subroutine
        
        
        function get_tar_sol_ind(this,imat) result(tar_sol_ind) !> function to get indices of target solids associated to a material
        class(chemistry_c) :: this !> chemistry object
        integer(kind=4), intent(in) :: imat !> index of material
        integer(kind=4), allocatable :: tar_sol_ind(:) !> indices of target solids associated to material
        integer(kind=4) :: i !> loop index
        integer(kind=4) :: j !> loop index
        integer(kind=4) :: dim !> dimension of tar_sol_ind
        dim=0 !> we initialize dimension
        do i=1,this%num_target_solids
            if (this%target_solids(i)%id==imat) then
                dim=dim+1 !> we increase dimension
            end if
        end do
        allocate(tar_sol_ind(dim))
        !i=1
        j=1 !> index for tar_sol_ind
        do i=1,this%num_target_solids
            if (this%target_solids(i)%id==imat) then
                tar_sol_ind(j)=i !> we set the index
                if (j<dim) then
                    j=j+1 !> we increase j
                else
                    exit !> we exit the loop
                end if
            end if
        end do
        end function
        
        function get_tar_gas_ind(this,igzn) result(tar_gas_ind) !> function to get indices of target solids associated to a gas zone
        class(chemistry_c) :: this !> chemistry object
        integer(kind=4), intent(in) :: igzn !> index of gas zone
        integer(kind=4), allocatable :: tar_gas_ind(:) !> indices of target gases associated to gas zone
        integer(kind=4) :: i,j !> loop indices
        integer(kind=4) :: dim !> dimension of tar_gas_ind
        dim=0 !> we initialize dimension
        do i=1,this%num_target_gases
            if (this%target_gases(i)%id==igzn) then
                dim=dim+1 !> we increase dimension
            end if
        end do
        allocate(tar_gas_ind(dim))
        !i=1
        j=1 !> index for tar_gas_ind
        do i=1,this%num_target_gases
            if (this%target_gases(i)%id==igzn) then
                tar_gas_ind(j)=i !> we set the index
                if (j<dim) then
                    j=j+1 !> we increase j
                else
                    exit !> we exit the loop
                end if
            end if
        end do
        end function

        subroutine allocate_sol_types(this)
        implicit none
        class(chemistry_c) :: this
        if (allocated(this%sol_types)) deallocate(this%sol_types)
        allocate(this%sol_types(this%num_wat_types)) !> we assume bijection between water types and solid types
        end subroutine allocate_sol_types

        subroutine allocate_react_zones_wat_types(this)
        implicit none
        class(chemistry_c) :: this
        if (allocated(this%react_zones_wat_types)) deallocate(this%react_zones_wat_types)
        allocate(this%react_zones_wat_types(this%num_wat_types)) !> we assume bijection between water types and reactive zones
        end subroutine allocate_react_zones_wat_types

        subroutine allocate_gas_zones_wat_types(this)
        implicit none
        class(chemistry_c) :: this
        if (allocated(this%gas_zones_wat_types)) deallocate(this%gas_zones_wat_types)
        allocate(this%gas_zones_wat_types(this%num_wat_types)) !> we assume bijection between water types and gas zones
        end subroutine allocate_gas_zones_wat_types

end module chemistry_m