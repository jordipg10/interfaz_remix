program main_interfaz
    use RT_1D_m
    implicit none
    !> Objects
    type(RT_1D_transient_c) :: my_RT_trans !> 1D transient reactive transport class
    type(transport_1D_transient_c) :: my_tpt_trans !> 1D transient transport class
    type(chemistry_c) :: my_chem !> chemistry class
    !> Variables
    integer(kind=4) :: opc_u_tilde !> option to compute or read u_tilde
    integer(kind=4) :: tpt_opt !> transport option (0: compute lambdas, 1: read lambdas)
    integer(kind=4) :: int_method_chem !> integration method for chemical reactions (1: Euler explicit, 2: Euler fully implicit, 3: Crank-Nicolson)
    integer(kind=4) :: num_comps !> number of components
    integer(kind=4) :: flag !> loop flag
    integer(kind=4) :: mesh_type !> type of mesh (1: 1D homogeneous, 2: 1D heterogeneous, 3: radial)
    integer(kind=4) :: num_can_vec !> number of canonical vectors
    integer(kind=4) :: num_non_can_vec !> number of canonical vectors
    integer(kind=4) :: n_nc_aq !> number of variable activity aqueous species
    integer(kind=4), allocatable :: ind_can_vec(:) !> indices of canonical vectors in mixing ratios
    integer(kind=4), allocatable :: ind_non_can_vec(:) !> indices of non-canonical vectors in mixing ratios
    real(kind=8), allocatable :: u_tilde_init(:,:) !> u_tilde at initial time step
    real(kind=8) :: Delta_t !> time step
    character(len=100) :: dir_DB,dir_pb, root_files,file_u_tilde,file_u_new
    character(len=:), allocatable :: dir_DB_trimmed, dir_pb_trimmed, root_files_trimmed, file_u_tilde_trimmed, file_u_new_trimmed
    !logical :: flag !> loop flag
    !> Process
    write(*,*) "Esta es la interfaz para resolver iteraciones de mezcla reactiva en 1D con el WMA usando Euler explícito."
    write(*,*) "Directorio de las bases de datos: "
    read(*,*) dir_DB
    dir_DB_trimmed = trim(dir_DB)
    write(*,*) "Directorio del problema: "
    read(*,*) dir_pb
    dir_pb_trimmed = trim(dir_pb)
    write(*,*) "Root de los archivos de entrada y salida: "
    read(*,*) root_files
    root_files_trimmed = trim(root_files)
    !> we read chemistry
    call my_chem%read_chemistry(dir_pb_trimmed//root_files_trimmed,dir_DB_trimmed)
    write(*,*) "Nombre del archivo donde quieres que escriba las concentraciones después de la mezcla reactiva?"
    read(*,*) file_u_new
    file_u_new_trimmed = trim(file_u_new)
    write(*,*) "Quieres que calcule las concentraciones después de resolver el transporte conservativo para el primer paso de tiempo (0) o que las lea de un archivo (1)?"
    read(*,*) opc_u_tilde
    if (opc_u_tilde.eq.0) then
        write(*,*) "De acuerdo. Vamos a calcular las concentraciones después de resolver el transporte conservativo para el primer paso de tiempo."
        write(*,*) "Nombre del archivo donde quieres que escriba estas concentraciones? (este archivo estará en el mismo directorio)."
        read(*,*) file_u_tilde
        write(*,*) "Quieres que calcule yo los lambdas (0) o que los lea de un archivo (1)?"
        read(*,*) tpt_opt
        if (tpt_opt.eq.0) then !> compute lambdas
        !> we read transport data, BCs and discretisations
            !write(*,*) "De acuerdo. Vamos a inicializar el transporte."
            write(*,*) "Tipo de malla ? (1: 1D homogénea, 2: 1D heterogénea, 3: radial):"
            read(*,*) mesh_type !> must be written by the user
            !> we initialise transport object
            call my_tpt_trans%initialise_transport_1D_transient_RT(dir_pb_trimmed//root_files_trimmed,mesh_type)
            !> we allocate transport arrays
            call my_tpt_trans%allocate_arrays_PDE_1D()
            !> we compute transport arrays, including mixing ratios, and we impose BCs
            call my_tpt_trans%compute_mixing_ratios_Delta_t_homog() !> missing the case of heterogenous time steps
            !> we set transport attribute in reactive transport object
            call my_RT_trans%set_transport_trans(my_tpt_trans)
            !> we choose and set integration method for chemical reactions
            !! 1: Euler explicit, 2: Euler fully implicit, 3: Crank-Nicolson
            !write(*,*) "Elegir método de integración temporal para las reacciones químicas (1: Euler explicito, & 
            !    2: Euler totalmente implicito, 3: Crank-Nicolson):"
            !read(*,*) int_method_chem !> must be written by the user
            call my_RT_trans%set_int_method_chem_reacts(1)
        else if (tpt_opt.eq.1) then !> read lambdas
        !> we set transport attribute in reactive transport object
            call my_RT_trans%set_transport_trans(my_tpt_trans) !> esto es un create en realidad
        !> we read temporal discretisation
            call my_RT_trans%read_time_discretisation(dir_pb_trimmed//root_files_trimmed)
        !> we read transport data for WMA
            call my_RT_trans%transport%read_transport_data_WMA(dir_pb_trimmed//root_files_trimmed)
        else
            error stop "Esta opción no está implementada todavía"
        end if
        !call my_tpt_trans%mixing_ratios_conc%get_can_vec(my_chem%CV_params%abs_tol,num_can_vec,ind_can_vec,num_non_can_vec,ind_non_can_vec)
        !n_nc_aq=my_chem%get_num_aq_var_act_species() !> number of variable activity aqueous species
        !allocate(c_tilde(n_nc_aq,my_chem%num_target_waters_dom)) !> we allocate c_tilde
        !do i=1,my_chem%num_target_waters_dom !> we compute c_tilde for each domain target water
        !    call compute_c_tilde(my_chem%target_waters_init(mixing_waters_indices%cols(i))%col_1(&
        !        1:my_tpt_trans%mixing_waters_indices%cols(i)%dim-2),&
        !        my_tpt_trans%mixing_ratios_conc%cols(i)%col_1,c_tilde(:,i)) !> we compute c_tilde for each domain target water
        !end do
        num_comps=my_chem%get_num_aq_comps() !> we get number of aqueous components
        Delta_t=my_tpt_trans%time_discr%get_Delta_t() !> we get time step
        u_tilde_init=my_chem%compute_u_tilde_init(my_tpt_trans%mixing_waters_indices,my_tpt_trans%mixing_ratios_conc) !> we compute u_tilde for first time step
        call my_chem%write_u_tilde_init(dir_pb_trimmed,trim(file_u_tilde),u_tilde_init) !> we write u_tilde_init to a file
    else if (opc_u_tilde.eq.1) then
        write(*,*) "De acuerdo. El archivo debe estar en el directorio del problema. Recuerda que el numero de filas es el numero de componentes y el numero de columnas es el numero de targets."
        write(*,*) "Nombre del archivo que contiene estas concentraciones?"
        read(*,*) file_u_tilde
        write(*,*) "Numero de componentes?"
        read(*,*) num_comps
        write(*,*) "Paso de tiempo?"
        read(*,*) Delta_t
    end if
    !flag = 1 !> to enter the loop
    file_u_tilde_trimmed = trim(file_u_tilde)
    write(*,*) "Procedemos al bucle de mezcla reactiva. Tendrás que actualizar el archivo con las concentraciones después del transporte conservativo en cada iteración."
    do
        call my_chem%interfaz_comps_arch(dir_pb_trimmed,num_comps,file_u_tilde_trimmed,Delta_t,file_u_new_trimmed)
        write(*,*) "Quieres hacer otra iteración de mezcla reactiva? (1: sí, 0: no)"
        read(*,*) flag
        if (flag.eq.0) then
            exit
        else if (flag.ne.1) then
            error stop "Opción no válida. Tiene que ser 1 o 0."
        !else
            !write(*,*) "Recuerda que el archivo con las concentraciones después de resolver el transporte conservativo debe estar en el directorio del problema."
        end if
    end do
    !> we set chemistry attribute in reactive transport object
    call my_RT_trans%set_chemistry(my_chem) !> we set chemistry attribute in reactive transport object
    write(*,*) "El programa ha terminado."
end program main_interfaz