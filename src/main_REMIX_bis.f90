!> \file main_REMIX_bis.f90
!> \brief Main program entry point for REMIX reactive transport simulator (file-driven configuration)
!> \details This is an alternative driver program for the REMIX (Reactive Transport Code) simulator.
!> It is identical to main_REMIX.f90 except for the configuration section: instead of hard-coding the
!> database path and selecting a problem via an integer ID and if-else block, it reads the
!> configuration directly from an input file.
!>
!> It handles the complete workflow for reactive transport modeling:
!> 1. Problem configuration (read dir_DB, dir_pb, root, mesh_type, is_2D from an input file)
!> 2. Transport initialization (spatial discretization, boundary conditions, mixing ratios)
!> 3. Chemistry initialization (thermodynamic database, reactions, target waters)
!> 4. Reactive transport solution (coupling transport and chemistry)
!> 5. Output generation (results files, Python visualization scripts)
!>
!> The code supports various solution approaches:
!> - Lagrangian vs. Eulerian transport
!> - Ideal vs. non-ideal activity coefficients
!> - Lumped vs. consistent mixing for reactions
!> - Explicit vs. implicit time integration
!> - Stationary vs. transient flow fields
!>
!> Usage: REMIX.exe <directory of main_REMIX_input.dat>
!>   The directory containing the configuration file 'main_REMIX_input.dat' is passed as the
!>   first command-line argument. The configuration file format is:
!>   Line 1: dir_DB    - path to the thermodynamic database directory (quoted string)
!>   Line 2: dir_pb    - path to the problem directory (quoted string)
!>   Line 3: root      - base filename prefix for input/output files (quoted string)
!>   Line 4: mesh_type - spatial mesh type; meaning depends on is_2D. 1D: 1=uniform Cartesian, 2=non-uniform, 3=radial. 2D: 1=uniform Cartesian (only uniform implemented)
!>   Line 5: is_2D     - .true. for 2D problems, .false. for 1D problems
!>   Line 6: option_lambdas - mixing ratios option (0=compute from transport matrices, 1=read from file)
!>   Line 7: option_tpt    - transport method (0=Eulerian, 1=Lagrangian)
!>
!> Output files generated:
!> - Concentration profiles: root_conc_*.dat
!> - Mineral volumes: root_minerals_*.dat
!> - Python scripts: root_plot.py for visualization
!>
!> \author Jordi Petchamé-Guerrero
!> \date November 2025

program main                                                          !< Main program entry point for REMIX reactive transport simulator
    use RT_m, only: RT_transient_c, RT_1D_transient_c, RT_2D_transient_c  !< Import reactive transport classes: base polymorphic class and 1D/2D concrete types
    use PDE_transient_m, only: PDE_transient_c                                  !< Import base PDE transient class used for polymorphic transport object declaration
    use transport_transient_m, only: transport_1D_transient_c, transport_2D_transient_c !< Import 1D and 2D transient transport concrete types for allocation
    use chemistry_m, only: chemistry_c                           !< Import chemistry class containing chemical system, waters, and solver methods
    use spatial_discr_2D_m, only: mesh_2D_Euler_homog_c          !< Import 2D homogeneous Eulerian mesh type for downstream water index computation
    use time_discr_m, only: time_discr_homog_c                   !< Import homogeneous time discretisation for WMA reading
    implicit none                                                                !< Require explicit declaration of all variables (Fortran safety feature to prevent undeclared variable bugs)
!> @name Polymorphic Objects
!> @brief Allocated at runtime based on problem dimensionality (1D or 2D)
!> @{
    class(RT_transient_c), pointer :: my_RT => null()                           !< Reactive transport object (pointer avoids gfortran SIGSEGV during auto-deallocation at end program) [-]
    class(PDE_transient_c), pointer :: my_tpt => null()                         !< Transport-only PDE object (pointer avoids gfortran SIGSEGV during auto-deallocation at end program) [-]
    type(time_discr_homog_c) :: time_discr_WMA                                  !< Homogeneous time discretisation object for WMA reading [-]
!> @}
!> @name Control Parameters
!> @brief User-configurable parameters that control problem selection and solver behaviour
!> @{
    integer(kind=4) :: option_tpt                                                !< Transport method flag: 0=Eulerian (fixed grid), 1=Lagrangian (particle tracking) [-]
    integer(kind=4) :: option_lambdas                                            !< Mixing ratios option: 0=compute lambdas from transport matrices, 1=read precomputed lambdas from file [-]
    integer(kind=4) :: mesh_type                                                !< Spatial mesh type (meaning depends on is_2D). 1D: 1=uniform Cartesian, 2=non-uniform, 3=radial. 2D: 1=uniform Cartesian (only uniform implemented) [-]
    integer(kind=4) :: int_method_chem                                          !< Time integration method for chemistry: 1=Euler explicit (EE), 2=Euler implicit (EI) [-]
    logical :: is_2D                                                            !< Flag for 2D problems: .true. if problem is 2D (Paulina cases), .false. for 1D problems [-]
    logical :: stat_flag                                                        !< Stationary flow flag: .true. if velocity field is constant in time, .false. for transient flow [-]
    logical :: apply_time_scaling                                               !< User flag: .true. to apply explicit time scaling to kinetic mineral reaction rates, .false. to disable [-]
    integer(kind=4) :: i_min                                                    !< Loop counter over kinetic mineral reactions for time-scaling configuration [-]
!> @}

!> @name File Path Variables
!> @brief Define directory locations and filename roots for input/output files
!> @{
    character(len=256) :: root                                                  !< Base filename prefix for input/output files (e.g., "calcite" produces calcite_chem_system.dat, etc.) [string]
    character(len=256) :: dir_DB                                               !< File system path to thermodynamic database directory containing equilibrium constants and species data [string]
    character(len=:), allocatable :: dir_DB_trim                               !< Trimmed (no trailing whitespace) database path for exact string operations in file I/O [string]
    character(len=256) :: dir_pb                                               !< File system path to specific problem directory containing all input files for the selected test case [string]
    character(len=:), allocatable :: dir_pb_trim                               !< Trimmed problem directory path for memory-efficient string operations in file I/O [string]
    character(len=:), allocatable :: root_trim                                  !< Trimmed root filename for consistent file naming across all input/output operations [string]
!> @}
!> @name Input File Variables
!> @brief Variables used to read the configuration from an input file
!> @{
    integer(kind=4), parameter :: unit_input = 10                               !< File unit number used to read the configuration input file [-]
    integer(kind=4), parameter :: num_header_lines = 4                          !< Number of heading/comment lines to skip at the top of the input file [-]
    integer(kind=4) :: i_header                                                  !< Loop counter used to skip heading lines [-]
    integer(kind=4) :: arg_len                                                   !< Length of the command-line argument holding the input directory [-]
    integer(kind=4) :: arg_status                                                !< Status flag returned when retrieving the command-line argument [-]
    character(len=256) :: header_line                                           !< Temporary buffer used to read and discard heading lines [string]
    character(len=256) :: input_file                                            !< Name of the configuration input file [string]
    character(len=256) :: input_dir                                             !< Directory containing the configuration input file [string]
    character(len=512) :: input_path                                            !< Full path to the configuration input file [string]
!> @}
!****************************************************************************************************************************
!> @section config_section Configuration Section
!> @brief Read the configuration (dir_DB, dir_pb, root, mesh_type, is_2D) from an input file
!> @{
    input_file = 'main_REMIX_input.dat'                                        !< Default name of the configuration input file [string]
    call get_command_argument(1, input_dir, arg_len, arg_status)              !< Read the input-file directory from the first command-line argument [-]
    if (arg_status /= 0 .or. arg_len == 0) then                                !< Check that a directory argument was actually supplied
        write(*,*) 'Usage: REMIX.exe <directory of main_REMIX_input.dat>'      !< Inform the user about the expected command-line usage [-]
        error stop 'Missing command-line argument: directory of the configuration input file' !< Terminate: no input directory provided [-]
    end if                                                                     !< End of argument-presence check
    input_dir = trim(adjustl(input_dir))                                       !< Remove surrounding whitespace from the supplied directory [string]
    if (len_trim(input_dir)>0 .and. input_dir(len_trim(input_dir):len_trim(input_dir))/='/' .and. &
        input_dir(len_trim(input_dir):len_trim(input_dir))/='\') then         !< Append a path separator if the directory does not end with one
        input_dir = trim(input_dir)//'\'                                       !< Add trailing backslash separator [string]
    end if                                                                     !< End of separator check
    input_path = trim(input_dir)//trim(input_file)                            !< Build the full path to the configuration input file [string]
    open(unit=unit_input, file=trim(input_path), status='old', action='read') !< Open the configuration input file for reading [-]
    do i_header = 1, num_header_lines                                          !< Skip the heading/comment lines at the top of the input file
        read(unit_input,'(A)') header_line                                     !< Read and discard one heading line (whole record) [-]
    end do                                                                     !< End of heading-skipping loop
    read(unit_input,*) dir_DB                                                  !< Read absolute path to thermodynamic database directory [string]
    read(unit_input,*) dir_pb                                                  !< Read absolute path to specific problem directory [string]
    read(unit_input,*) root                                                    !< Read root filename prefix for all input/output files [string]
    read(unit_input,*) mesh_type                                               !< Read spatial mesh type (1D: 1=uniform Cartesian, 2=non-uniform, 3=radial; 2D: 1=uniform Cartesian, only uniform implemented) [-]
    read(unit_input,*) is_2D                                                   !< Read 2D problem flag (.true.=2D, .false.=1D) [logical]
    read(unit_input,*) option_lambdas                                          !< Read mixing ratios option (0=compute from transport matrices, 1=read from file) [-]
    read(unit_input,*) option_tpt                                              !< Read transport method flag (0=Eulerian, 1=Lagrangian) [-]
    close(unit_input)                                                          !< Close the configuration input file [-]
    dir_DB_trim = trim(dir_DB)                                                !< Remove trailing whitespace from database path using Fortran intrinsic trim() for exact string matching in file open() calls [string]
!> @}
    dir_pb_trim=trim(dir_pb)                                                  !< Remove trailing whitespace from problem directory path using Fortran intrinsic trim() for string operations [string]
    root_trim=trim(root)                                                        !< Remove trailing whitespace from root filename using Fortran intrinsic trim() for consistent file naming [string]
!> @section alloc_section Polymorphic Object Allocation
!> @brief Allocate RT and transport objects to correct concrete type based on dimensionality
    if (is_2D) then                                                            !< Check if current problem is a 2D case
        allocate(RT_2D_transient_c :: my_RT)                                   !< Allocate reactive transport object as 2D transient type (RT_2D_transient_c) for 2D domain [-]
        allocate(transport_2D_transient_c :: my_tpt)                           !< Allocate transport PDE object as 2D transient type (transport_2D_transient_c) for 2D mesh [-]
    else                                                                       !< Otherwise current problem is a 1D case
        allocate(RT_1D_transient_c :: my_RT)                                   !< Allocate reactive transport object as 1D transient type (RT_1D_transient_c) for 1D domain [-]
        allocate(transport_1D_transient_c :: my_tpt)                           !< Allocate transport PDE object as 1D transient type (transport_1D_transient_c) for 1D mesh [-]
    end if                                                                     !< End of dimensionality-based allocation block
!> @section transport_init Transport Initialization Section
!> @brief Initialize transport parameters, spatial discretisation, BCs, and compute mixing ratios
    if (option_tpt.eq.0) then                                                  !< Check if Eulerian transport method was selected (option_tpt=0)
        call my_tpt%set_Lagr_flag(.false.)                                     !< Set Lagrangian flag to .false. on the transport object (use Eulerian fixed-grid approach) [-]
    else if (option_tpt.eq.1) then                                             !< Check if Lagrangian transport method was selected (option_tpt=1)
        call my_tpt%set_Lagr_flag(.true.)                                      !< Set Lagrangian flag to .true. on the transport object (use particle tracking approach) [-]
    else                                                                       !< Invalid transport option: neither 0 (Eulerian) nor 1 (Lagrangian)
        error stop "This option is not implemented yet"                        !< Terminate program with error if transport option is unrecognised [-]
    end if                                                                     !< End of transport method selection block
    if (option_lambdas.eq.0) then                                              !< Check if mixing ratios should be computed from scratch (option_lambdas=0)
        call my_tpt%initialise_PDE(dir_pb_trim,root_trim,mesh_type)            !< Initialise PDE: read spatial discretisation, BCs, time discretisation, and transport properties from input files [-]
        call my_tpt%allocate_arrays_PDE()                                      !< Allocate internal arrays for transport matrices (tridiagonal system, mixing ratios, boundary vectors) [-]
        select type (tpt => my_tpt)                                            !< Polymorphic dispatch: select concrete type of transport object for type-specific operations
        type is (transport_1D_transient_c)                                     !< Match 1D transient transport type (allocate_conc_init is 1D-specific)
            call tpt%allocate_conc_init()                                      !< Allocate initial concentration arrays for 1D transport (conc_init, conc_old, conc_old_old) [-]
        end select                                                             !< End of type-specific transport allocation block
        call my_tpt%compute_mixing_ratios_Delta_t_homog()                      !< Compute mixing ratios matrices (lambdas) from transport discretisation for homogeneous time step [-]
        if (my_tpt%Lagr_flag) then                                             !< Lagrangian transport: use '_Lagr' tag in the mixing ratios filenames [-]
            call my_tpt%write_mixing_ratios_conc_rows(dir_pb_trim//root_trim//'_Lagr_mixing_ratios_conc.out') !< Dense concentration mixing ratios by rows: col 1=inflow, cols 2..N+1=domain, col N+2=outflow [-]
            call my_tpt%write_mixing_ratios_R_rows(dir_pb_trim//root_trim//'_Lagr_mixing_ratios_R.out')       !< Dense reaction mixing ratios by rows (row i of A^{-1} over the N domain targets) [-]
        else                                                                   !< Eulerian transport: use '_Euler' tag in the mixing ratios filenames [-]
            call my_tpt%write_mixing_ratios_conc_rows(dir_pb_trim//root_trim//'_Euler_mixing_ratios_conc.out') !< Dense concentration mixing ratios by rows: col 1=inflow, cols 2..N+1=domain, col N+2=outflow [-]
            call my_tpt%write_mixing_ratios_R_rows(dir_pb_trim//root_trim//'_Euler_mixing_ratios_R.out')       !< Dense reaction mixing ratios by rows (row i of A^{-1} over the N domain targets) [-]
        end if                                                                 !< End of transport-type filename tag selection [-]
        call my_RT%set_transport_trans(my_tpt)                                 !< Copy transport object into reactive transport object (sets my_RT%transport from my_tpt) [-]
        !call my_RT%set_int_method_chem_reacts()                               !< Set chemistry integration method to 2 (Euler Implicit) — only implemented option [-]
    else if (option_lambdas.eq.1) then                                         !< Check if mixing ratios should be read from precomputed file (option_lambdas=1)
        call my_RT%set_transport_trans(my_tpt)                                 !< Copy transport object into reactive transport object before reading transport data [-]
        call time_discr_WMA%read_time_discr_WMA(dir_pb_trim,root_trim,int_method_chem) !< Read WMA time discretisation from file [-]
        select type (rt => my_RT)                                              !< Polymorphic dispatch: select concrete RT type for 1D-specific read operation
        type is (RT_1D_transient_c)                                            !< Match 1D transient reactive transport type
            call rt%transport%set_time_discr(time_discr_WMA)                   !< Assign WMA time discretisation to transport object [-]
            call rt%transport%read_transport_data_WMA(dir_pb_trim,root_trim)   !< Read precomputed WMA (Water Mixing Approach) transport data: mixing ratios, indices from file [-]
        end select                                                             !< End of type-specific transport data reading block
    else                                                                       !< Invalid mixing ratios option: neither 0 (compute) nor 1 (read from file)
        error stop "This option is not implemented yet"                        !< Terminate program with error if mixing ratios option is unrecognised [-]
    end if                                                                     !< End of mixing ratios computation/reading block
!> @section stat_flag_section Extract Stationary Flow Flag
!> @brief Determine if the flow field is stationary or transient via polymorphic dispatch
    select type (rt => my_RT)                                                  !< Polymorphic dispatch: select concrete RT type to access transport properties
    type is (RT_1D_transient_c)                                                !< Match 1D transient reactive transport type
        stat_flag = rt%transport%tpt_props_heterog%stat_flag                   !< Extract stationary flow flag from 1D transport properties (.true.=constant velocity) [logical]
    type is (RT_2D_transient_c)                                                !< Match 2D transient reactive transport type
        stat_flag = rt%transport%tpt_props_heterog%stat_flag                   !< Extract stationary flow flag from 2D transport properties (.true.=constant velocity) [logical]
    end select                                                                 !< End of stationary flag extraction block
!> @section chemistry_init Chemistry Initialization Section
!> @brief Read chemical system definition, set reaction parameters, and initialise target waters
    call my_RT%chemistry%read_chemistry(root_trim,dir_pb_trim,dir_DB_trim)     !< Read complete chemistry setup: chemical system (species, reactions, stoichiometry), thermodynamic database, initial/boundary waters, output options [-]
    apply_time_scaling=.false.                                                  !< Set time-scaling flag — USER MUST MODIFY THIS VALUE TO ENABLE/DISABLE EXPLICIT TIME SCALING OF MINERAL REACTION RATES [-]
    do i_min=1,my_RT%chemistry%chem_syst%num_minerals_kin                       !< Loop over all kinetic mineral reactions in the chemical system
        call my_RT%chemistry%chem_syst%min_kin_reacts(i_min)%set_time_scaling(apply_time_scaling) !< Configure explicit time-scaling flag for each kinetic mineral reaction [-]
    end do                                                                      !< End of time-scaling configuration loop
    !> chemistry integration method now lives in transport%time_discr%int_method_chem (set from theta_r during time-discretisation read) [-]
    select type (rt => my_RT)                                                  !< Polymorphic dispatch: select concrete RT type for dimension-specific chemistry initialisation
    type is (RT_1D_transient_c)                                                !< Match 1D transient reactive transport type
        call rt%chemistry%set_target_solids_mesh(rt%transport%spatial_discr)   !< Assign mineral assemblages (solid chemistry) to spatial mesh cells based on reactive zone mapping [-]
        call rt%set_target_waters(rt%transport%spatial_discr,rt%transport%Lagr_flag,rt%transport%BCs%labels,dir_pb_trim,root_trim) !< Initialise target waters at each mesh cell: set initial concentrations, link to reactive zones, configure BCs [-]
        if (rt%transport%Lagr_flag) then                                       !< Check if Lagrangian transport is active (particle tracking requires upstream/downstream indices)
            block                                                              !< Begin local scope block for temporary array allocation (Fortran 2008 feature)
                real(kind=8), allocatable :: vy_zero(:)                        !< Temporary array for zero y-velocity component (1D has no transverse velocity) [m/s]
                allocate(vy_zero(rt%chemistry%num_target_waters))              !< Allocate zero y-velocity array sized to number of target waters [-]
                vy_zero = 0d0                                                  !< Set all y-velocity components to zero (1D problem has no transverse flow) [m/s]
                call rt%chemistry%set_downstream_water_indices(rt%transport%mix_conc_indices, & !< Compute downstream water indices for each target water based on flux direction
                    rt%chemistry%num_target_waters, rt%transport%tpt_props_heterog%flux_cent, vy_zero) !< Pass mixing indices, number of targets, x-flux at cell centres, and zero y-flux [-]
                call rt%chemistry%set_upstream_water_indices(rt%transport%mix_conc_indices, & !< Compute upstream water indices for each target water based on flux direction (for initial guess in Newton)
                    rt%chemistry%num_target_waters, rt%transport%tpt_props_heterog%flux_cent, vy_zero) !< Pass mixing indices, number of targets, x-flux at cell centres, and zero y-flux [-]
            end block                                                          !< End local scope block (vy_zero is automatically deallocated here)
        end if                                                                 !< End of Lagrangian-specific upstream/downstream index computation for 1D
    type is (RT_2D_transient_c)                                                !< Match 2D transient reactive transport type
        call rt%chemistry%set_target_solids_mesh(rt%transport%spatial_discr)   !< Assign mineral assemblages (solid chemistry) to 2D spatial mesh cells based on reactive zone mapping [-]
        call rt%set_target_waters(rt%transport%spatial_discr,rt%transport%Lagr_flag,rt%transport%BCs%labels,dir_pb_trim,root_trim) !< Initialise target waters at each 2D mesh cell: set initial concentrations, link to reactive zones, configure BCs [-]
        call rt%remap_mix_indices()                                            !< Remap transport-local water indices to global chemistry water indices for 2D [-]
        if (rt%transport%Lagr_flag) then                                       !< Check if Lagrangian transport is active in 2D (particle tracking requires upstream/downstream indices)
            select type (mesh => rt%transport%spatial_discr)                   !< Polymorphic dispatch on spatial discretisation to access 2D mesh-specific data (number of cells in x)
            type is (mesh_2D_Euler_homog_c)                                    !< Match 2D homogeneous Eulerian mesh type (structured rectangular grid)
                call rt%chemistry%set_downstream_water_indices(rt%transport%mix_conc_indices, & !< Compute downstream water indices for 2D based on both x and y flux directions
                    mesh%Num_cells_x, rt%transport%tpt_props_heterog%flux_cent, rt%transport%tpt_props_heterog%flux_cent_y) !< Pass mixing indices, number of x-cells, x-flux at centres, y-flux at centres [-]
                call rt%chemistry%set_upstream_water_indices(rt%transport%mix_conc_indices, & !< Compute upstream water indices for 2D based on both x and y flux directions (for Newton initial guess)
                    mesh%Num_cells_x, rt%transport%tpt_props_heterog%flux_cent, rt%transport%tpt_props_heterog%flux_cent_y) !< Pass mixing indices, number of x-cells, x-flux at centres, y-flux at centres [-]
            class default                                                      !< No match: unsupported 2D mesh type encountered
                error stop "Unsupported 2D mesh type for downstream water indices" !< Terminate with error if mesh type is not mesh_2D_Euler_homog_c [-]
            end select                                                         !< End of 2D mesh type selection block
            call rt%set_target_solids_time()                                   !< Set each target solid's Lagrangian travel time from the injection (left) boundary; required by kinetic rate time-scaling (time**(-m)) to avoid division by zero at time=0 [-]
        end if                                                                 !< End of Lagrangian-specific upstream/downstream index computation for 2D
    end select                                                                 !< End of dimension-specific chemistry initialisation block
    call my_RT%chemistry%write_u_init(dir_pb_trim,root_trim//'_u_init.out')    !< Write initial aqueous components and species (grouped by reactive zone) to root_u_init.out [-]
    !> Finalize output time steps before dispatch: last element = Num_time (expand empty request) [-]
    select type (rt => my_RT)
    type is (RT_1D_transient_c)
        call rt%chemistry%chem_out_options%finalize_time_steps(rt%transport%time_discr%Num_time)
    type is (RT_2D_transient_c)
        call rt%chemistry%chem_out_options%finalize_time_steps(rt%transport%time_discr%Num_time)
    end select
!> @section reactive_transport_solver Reactive Transport Solution Section
!> @brief Solver dispatch: selects the appropriate solve_RT_* method based on:
!> - Stationary vs. transient flow (stat_flag)
!> - Lagrangian vs. Eulerian transport (Lagr_flag)
!> - Ideal vs. non-ideal activity coefficients (act_coeffs_model)
!> - Lumped vs. consistent mixing ratios (lump_flag)
    if (stat_flag) then                                                        !< Check if flow field is stationary (constant velocity in time)
        if (my_tpt%Lagr_flag) then                                             !< Check if Lagrangian particle tracking is active within stationary flow
            if (my_RT%chemistry%act_coeffs_model==0 .and. (my_RT%chemistry%lump_flag .eqv. .true.)) then !< Check for ideal activity coefficients (gamma=1) AND lumped mixing ratios
                call my_RT%solve_RT_ideal_lump_Lagr_stat_flux(dir_pb_trim,root_trim) !< Solve RT with: ideal speciation, lumped WMA, Lagrangian tracking, stationary flux [-]
            else if (my_RT%chemistry%act_coeffs_model==0 .and. (my_RT%chemistry%lump_flag .eqv. .false.)) then !< Check for ideal activity coefficients (gamma=1) AND consistent (non-lumped) mixing ratios
                call my_RT%solve_RT_ideal_cons_Lagr_stat_flux(dir_pb_trim,root_trim) !< Solve RT with: ideal speciation, consistent WMA, Lagrangian tracking, stationary flux [-]
            else if (my_RT%chemistry%act_coeffs_model>0 .and. (my_RT%chemistry%lump_flag .eqv. .true.)) then !< Check for non-ideal activity coefficients AND lumped mixing ratios
                call my_RT%solve_RT_lump_Lagr_stat_flux(dir_pb_trim,root_trim) !< Solve RT with: non-ideal speciation, lumped WMA, Lagrangian tracking, stationary flux [-]
            end if                                                             !< End of activity model / lumping flag selection for Lagrangian stationary flow
        else                                                                   !< Eulerian transport with stationary flow
            select type (rt => my_RT)                                          !< Polymorphic dispatch: select concrete RT type to access transport members for Eulerian solver
            type is (RT_1D_transient_c)                                        !< Match 1D transient reactive transport type
                call rt%chemistry%solve_reactive_mixing(dir_pb_trim,root_trim,rt%transport%mixing_ratios_conc,& !< Solve 1D Eulerian reactive mixing: pass concentration mixing ratios matrix
                   rt%transport%mixing_ratios_R,rt%transport%mix_conc_indices,& !< Pass reaction mixing ratios (same layout as the concentration ones) and mixing water indices
                   rt%transport%time_discr,&                                   !< Pass time discretisation (Num_time, Delta_t sequence)
                   rt%transport%time_discr%theta_r)                            !< Pass theta_r time weighting factor for chemical reactions [-]
            type is (RT_2D_transient_c)                                        !< Match 2D transient reactive transport type
                call rt%chemistry%solve_reactive_mixing(dir_pb_trim,root_trim,rt%transport%mixing_ratios_conc,& !< Solve 2D Eulerian reactive mixing: pass concentration mixing ratios matrix
                   rt%transport%mixing_ratios_R,rt%transport%mix_conc_indices,& !< Pass reaction mixing ratios (same layout as the concentration ones) and mixing water indices
                   rt%transport%time_discr,&                                   !< Pass time discretisation (Num_time, Delta_t sequence)
                   rt%transport%time_discr%theta_r)                            !< Pass theta_r time weighting factor for chemical reactions [-]
            end select                                                         !< End of polymorphic dispatch for Eulerian stationary solver
        end if                                                                 !< End of Lagrangian vs. Eulerian selection for stationary flow
    else                                                                       !< Transient flow field (time-varying velocity)
        if (my_RT%chemistry%act_coeffs_model==0 .and. my_RT%chemistry%lump_flag .eqv. .true.) then !< Check for ideal activity coefficients AND lumped mixing ratios in transient flow
            if (my_tpt%Lagr_flag) then                                         !< Check if Lagrangian transport is active in transient flow
                call my_RT%solve_RT_ideal_lump_Lagr_trans_flux(dir_pb_trim,root_trim) !< Solve RT with: ideal speciation, lumped WMA, Lagrangian tracking, transient flux [-]
            else                                                               !< Eulerian transport with transient flow
                call my_RT%solve_RT_ideal_lump_Euler_trans_flux(dir_pb_trim,root_trim) !< Solve RT with: ideal speciation, lumped WMA, Eulerian grid, transient flux [-]
            end if                                                             !< End of Lagrangian vs. Eulerian selection for transient flow
        else                                                                   !< Non-ideal activity coefficients or consistent mixing with transient flow (not yet implemented)
            error stop "RT solver not implemented yet"                         !< Terminate with error: transient flow with non-ideal or consistent mixing is not available [-]
        end if                                                                 !< End of activity model / lumping flag selection for transient flow
    end if                                                                     !< End of stationary vs. transient flow solver dispatch
!> @section output_results Output Results Section
!> @brief Write simulation results to output files and generate visualisation scripts
    call my_RT%write_RT(dir_pb_trim,root_trim)                                 !< Write reactive transport results: concentration profiles, mineral volumes, and Python plotting scripts to problem directory [-]
    write(*,*) 'End of the simulation'                                         !< Print termination message to standard output confirming successful completion of the simulation [-]
end program                                                                    !< End of main program (my_RT and my_tpt are pointers so no automatic deallocation occurs)
