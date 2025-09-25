program main_interfaz
    use RT_1D_m
    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: ieee_exceptions
    implicit none
    !> Objects
    type(RT_1D_transient_c) :: my_RT_trans !> 1D transient reactive transport class
    type(transport_1D_transient_c) :: my_tpt_trans !> 1D transient transport class
    type(chemistry_c) :: my_chem !> chemistry class
    !> Variables
    integer(kind=4) :: opc_u_tilde !> option to compute or read u_tilde
    integer(kind=4) :: tpt_opt !> transport option (0: compute lambdas, 1: read lambdas)
    integer(kind=4) :: int_method_chem !> integration method for chemical reactions (1: Euler explicit, 2: Euler fully implicit, 3: Crank-Nicolson)
    integer(kind=4) :: num_aq_comps !> number of aqueous components
    integer(kind=4) :: flag !> loop flag
    integer(kind=4) :: flag_Delta_t !> whether time step is constant (1) or variable (0)
    integer(kind=4) :: mesh_type !> type of mesh (1: 1D homogeneous, 2: 1D heterogeneous, 3: radial)
    integer(kind=4) :: num_can_vec !> number of canonical vectors
    integer(kind=4) :: num_non_can_vec !> number of canonical vectors
    integer(kind=4) :: n_nc_aq !> number of variable activity aqueous species
    integer(kind=4), allocatable :: ind_can_vec(:) !> indices of canonical vectors in mixing ratios
    integer(kind=4), allocatable :: ind_non_can_vec(:) !> indices of non-canonical vectors in mixing ratios
    real(kind=8), allocatable :: u_tilde_init(:,:) !> u_tilde at initial time step
    real(kind=8) :: Delta_t !> time step
    character(len=100) :: dir_DB,dir_pb, root_files,file_u_tilde,file_u_new !> directories and file names
    character(len=:), allocatable :: dir_DB_trimmed, dir_pb_trimmed, root_files_trimmed, file_u_tilde_trimmed, file_u_new_trimmed !> trimmed strings
    integer :: ios !> I/O status for safe reads
    logical :: has_ieee !> whether IEEE support is available
    !> Pre-Process
    has_ieee = ieee_support_standard()
    if (has_ieee) then
        call ieee_set_underflow_mode(gradual=.false.)
        call ieee_set_flag(ieee_underflow, .false.)
        call ieee_set_flag(ieee_inexact, .false.)
        call ieee_set_flag(ieee_overflow, .false.)
        call ieee_set_flag(ieee_divide_by_zero, .false.)
        call ieee_set_flag(ieee_invalid, .false.)
    end if
    !> Process
    write(*,*) "Esta es la interfaz para resolver iteraciones de mezcla reactiva en 1D con el WMA usando Euler explicito."
    write(*,*) "Directorio de las bases de datos: "
    read(*,*, iostat=ios) dir_DB
    if (ios /= 0) then
        write(*,*) 'Error/EOF leyendo dir_DB. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
    end if
    dir_DB_trimmed = trim(dir_DB)
    write(*,*) "Directorio del problema: "
    read(*,*, iostat=ios) dir_pb
    if (ios /= 0) then
        write(*,*) 'Error/EOF leyendo dir_pb. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
    end if
    dir_pb_trimmed = trim(dir_pb)
    write(*,*) "Root de los archivos de entrada y salida: "
    read(*,*, iostat=ios) root_files
    if (ios /= 0) then
        write(*,*) 'Error/EOF leyendo root_files. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
    end if
    root_files_trimmed = trim(root_files)
    !> we read chemistry
    call my_chem%read_chemistry(root_files_trimmed,dir_pb_trimmed,dir_DB_trimmed)
    write(*,*) "Nombre del archivo donde quieres que escriba las concentraciones despues de la mezcla reactiva?"
    read(*,*, iostat=ios) file_u_new
    if (ios /= 0) then
        write(*,*) 'Error/EOF leyendo file_u_new. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
    end if
    file_u_new_trimmed = trim(file_u_new)
    write(*,*) "Quieres que calcule las concentraciones despues de resolver el transporte conservativo para el primer paso de & 
        tiempo (0) o que las lea de un archivo (1)?"
    read(*,*, iostat=ios) opc_u_tilde
    if (ios /= 0) then
        write(*,*) 'Error/EOF leyendo opc_u_tilde. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
    end if
    if (opc_u_tilde.eq.0) then !> compute u_tilde
        write(*,*) "De acuerdo. Vamos a calcular las concentraciones despues de resolver el transporte conservativo para el primer &
            paso de tiempo."
        write(*,*) "Nombre del archivo donde quieres que escriba estas concentraciones? (este archivo estara en el mismo &
            directorio)."
        read(*,*, iostat=ios) file_u_tilde
        if (ios /= 0) then
            write(*,*) 'Error/EOF leyendo file_u_tilde. Ejecuta en una terminal interactiva o redirige desde fort.5'; call & 
                safe_stop(1)
        end if
        write(*,*) "Quieres que calcule yo los lambdas (0) o que los lea de un archivo (1)?"
        read(*,*, iostat=ios) tpt_opt
        if (ios /= 0) then
            write(*,*) 'Error/EOF leyendo tpt_opt. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
        end if
        if (tpt_opt.eq.0) then !> compute lambdas
            !> we read transport data, BCs and discretisations
            write(*,*) "Tipo de malla ? (1: 1D homogenea, 2: 1D heterogenea, 3: radial):"
            read(*,*, iostat=ios) mesh_type !> must be written by the user
            if (ios /= 0) then
                write(*,*) 'Error/EOF leyendo mesh_type. Ejecuta en una terminal interactiva o redirige desde fort.5'; call & 
                    safe_stop(1)
            end if
            !> we initialise transport object
            call my_tpt_trans%initialise_transport_1D_transient_RT(dir_pb_trimmed//root_files_trimmed,mesh_type)
            !> we allocate transport arrays
            call my_tpt_trans%allocate_arrays_PDE_1D()
            !> we compute transport arrays, including mixing ratios, and we impose BCs
            call my_tpt_trans%compute_mixing_ratios_Delta_t_homog() !> missing the case of heterogenous time steps
            !> we set transport attribute in reactive transport object
            call my_RT_trans%set_transport_trans(my_tpt_trans)
            !> we set Euler explicit as integration method for chemical reactions
            call my_RT_trans%set_int_method_chem_reacts(1)
        else if (tpt_opt.eq.1) then !> read lambdas
            !> we set transport attribute in reactive transport object
            call my_RT_trans%set_transport_trans(my_tpt_trans) !> esto es un create en realidad
            !> we read temporal discretisation
            call my_RT_trans%read_time_discretisation(dir_pb_trimmed//root_files_trimmed)
            !> we read transport data for WMA
            call my_RT_trans%transport%read_transport_data_WMA(dir_pb_trimmed//root_files_trimmed)
        else
            error stop "Esta opcion no esta implementada todavia"
        end if
        Delta_t=my_tpt_trans%time_discr%get_Delta_t() !> we get initial time step
        select type (time_discr=>my_tpt_trans%time_discr)
        type is (time_discr_homog_c)
            flag_Delta_t=1 !> we set flag to constant time step
        type is (time_discr_heterog_c)
            flag_Delta_t=0 !> we set flag to variable time step
        end select
        u_tilde_init=my_chem%compute_u_tilde_init(my_tpt_trans%mixing_waters_indices,my_tpt_trans%mixing_ratios_conc) !> we compute u_tilde for first time step
        call my_chem%write_u_tilde_init(dir_pb_trimmed,trim(file_u_tilde),u_tilde_init) !> we write u_tilde_init to a file
    else if (opc_u_tilde.eq.1) then !> read u_tilde
        write(*,*) "De acuerdo. El archivo debe estar en el directorio del problema. Recuerda que el numero de filas es el numero &
            de componentes y el numero de columnas es el numero de targets."
        write(*,*) "Nombre del archivo que contiene estas concentraciones?"
        read(*,*, iostat=ios) file_u_tilde !> file with u_tilde
        if (ios /= 0) then
            write(*,*) 'Error/EOF leyendo file_u_tilde. Ejecuta en una terminal interactiva o redirige desde fort.5'; call & 
                safe_stop(1)
        end if
        write(*,*) "Paso de tiempo inicial?"
        read(*,*, iostat=ios) Delta_t !> initial time step
        if (ios /= 0) then
            write(*,*) 'Error/EOF leyendo Delta_t. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
        else if (Delta_t <= 0d0) then
            error stop "El paso de tiempo debe ser mayor que cero"
        end if
        write(*,*) "Paso de tiempo constante? (1: si, 0: no)"
        read(*,*, iostat=ios) flag_Delta_t !> whether time step is constant (1) or variable (0)
        if (ios /= 0) then
            write(*,*) 'Error/EOF leyendo flag_Delta_t. Ejecuta en una terminal interactiva o redirige desde fort.5'; call & 
                safe_stop(1)
        end if
    end if
    file_u_tilde_trimmed = trim(file_u_tilde) !> we trim file name
    write(*,*) "Procedemos al bucle de mezcla reactiva. Tendras que actualizar el archivo con las concentraciones de & 
        componentes obtenidas despues de resolver el transporte conservativo en cada iteracion."
    num_aq_comps=my_chem%get_num_aq_comps() !> we get number of aqueous components
    if (flag_Delta_t.eq.1) then !> constant time step
        do
            call my_chem%interfaz_comps_arch(dir_pb_trimmed,num_aq_comps,file_u_tilde_trimmed,Delta_t,file_u_new_trimmed)
            write(*,*) "Quieres hacer otra iteracion de mezcla reactiva? (1: si, 0: no)"
            read(*,*, iostat=ios) flag
            if (ios /= 0) then
                write(*,*) 'Error/EOF leyendo flag. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
            end if
            if (flag.eq.0) then
                exit
            else if (flag.ne.1) then
                error stop "Opcion no valida. Tiene que ser 1 o 0."
            end if
        end do
    else if (flag_Delta_t.eq.0) then !> variable time step
        do
            call my_chem%interfaz_comps_arch(dir_pb_trimmed,num_aq_comps,file_u_tilde_trimmed,Delta_t,file_u_new_trimmed)
            write(*,*) "Quieres hacer otra iteracion de mezcla reactiva? (1: si, 0: no)"
            read(*,*, iostat=ios) flag
            if (ios /= 0) then
                write(*,*) 'Error/EOF leyendo flag. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
            end if
            if (flag.eq.0) then
                exit
            else if (flag.ne.1) then
                error stop "Opcion no valida. Tiene que ser 1 o 0."
            else
                write(*,*) "Nuevo paso de tiempo?"
                read(*,*, iostat=ios) Delta_t
                if (ios /= 0) then
                    write(*,*) 'Error/EOF leyendo Delta_t. Ejecuta en una terminal interactiva o redirige desde fort.5'; call & 
                        safe_stop(1)
                else if (Delta_t <= 0d0) then
                    error stop "El paso de tiempo debe ser mayor que cero"
                end if
            end if
        end do
    end if
    !> Post-Process
    if (has_ieee) call clear_ieee_flags()
    write(*,*) "El programa ha terminado."
contains
    subroutine clear_ieee_flags()
        if (ieee_support_standard()) then
            call ieee_set_flag(ieee_underflow, .false.)
            call ieee_set_flag(ieee_inexact, .false.)
            call ieee_set_flag(ieee_overflow, .false.)
            call ieee_set_flag(ieee_divide_by_zero, .false.)
            call ieee_set_flag(ieee_invalid, .false.)
        end if
    end subroutine clear_ieee_flags

    subroutine safe_stop(code)
        integer, intent(in) :: code
        call clear_ieee_flags()
        stop code
    end subroutine safe_stop
end program main_interfaz