program main_interfaz
    use chemistry_m
    use, intrinsic :: ieee_arithmetic
    use, intrinsic :: ieee_exceptions
    implicit none
    !> Objects
    type(chemistry_c) :: my_chem !> chemistry class
    !> Variables
    integer(kind=4) :: int_method_chem !> integration method for chemical reactions (1: Euler explicit, 2: Euler fully implicit, 3: Crank-Nicolson)
    integer(kind=4) :: num_st_var !> number of state variables (components or species)
    integer(kind=4) :: flag !> loop flag
    integer(kind=4) :: flag_Delta_t !> whether time step is constant (1) or variable (0)
    integer(kind=4), allocatable :: ind_can_vec(:) !> indices of canonical vectors in mixing ratios
    integer(kind=4), allocatable :: ind_non_can_vec(:) !> indices of non-canonical vectors in mixing ratios
    real(kind=8), allocatable :: u_tilde_init(:,:) !> u_tilde at initial time step
    real(kind=8) :: Delta_t !> time step
    character(len=100) :: dir_DB,dir_pb, root_files,file_u_tilde,file_u_new,file_u_wat_types !> directories and file names
    character(len=:), allocatable :: dir_DB_trimmed, dir_pb_trimmed, root_files_trimmed, file_u_tilde_trimmed, file_u_new_trimmed,&
        file_u_wat_types_trimmed !> trimmed strings
    integer :: ios !> I/O status for safe reads
    logical :: has_ieee !> whether IEEE support is available
    procedure(interfaz_comps_arch), pointer :: p_interfaz=>null() !> procedure pointer to the interface subroutine
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
    !> write file name
    write(*,*) "Nombre del archivo donde quieres que escriba las concentraciones de los tipos de agua iniciales y externas?"
    read(*,*, iostat=ios) file_u_wat_types !> file with initial and external water types
    if (ios /= 0) then
        write(*,*) 'Error/EOF leyendo file_u_wat_types. Ejecuta en una terminal interactiva o redirige desde fort.5'; call & 
            safe_stop(1)
    end if
    file_u_wat_types_trimmed = trim(file_u_wat_types)
    !> we write concentrations of initial and external water types
    call my_chem%write_conc_comp_wat_types(dir_pb_trimmed,file_u_wat_types_trimmed)
    !> write file name 
    write(*,*) "Nombre del archivo donde quieres que escriba las concentraciones despues de la mezcla reactiva?"
    read(*,*, iostat=ios) file_u_new
    if (ios /= 0) then
        write(*,*) 'Error/EOF leyendo file_u_new. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
    end if
    file_u_new_trimmed = trim(file_u_new)
    write(*,*) "Nombre del archivo que contiene las concentraciones despues de resolver una iteracion de transporte conservativo? & 
        IMPORTANTE: El archivo debe estar en el directorio del problema, el numero de filas tiene que ser el numero &
            de componentes y el numero de columnas el numero de targets."
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
    file_u_tilde_trimmed = trim(file_u_tilde) !> we trim file name
    write(*,*) "Procedemos al bucle de mezcla reactiva. Tendras que actualizar el archivo con las concentraciones & 
        obtenidas despues de resolver el transporte conservativo en cada iteracion."
    !> We choose interface based on whether there are equilibrium reactions or not, using procedure pointers
    if (my_chem%target_waters(my_chem%dom_tar_wat_indices(1))%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions==0) then
        p_interfaz=>interfaz_esp_arch
        num_st_var=my_chem%chem_syst%speciation_alg%num_aq_var_act_species !> we get number of aqueous variable activity species
    else
        p_interfaz=>interfaz_comps_arch
        num_st_var=my_chem%get_num_aq_comps_dom() !> we get number of aqueous components
    end if
    if (flag_Delta_t.eq.1) then !> constant time step
        do
            call p_interfaz(my_chem,dir_pb_trimmed,num_st_var,file_u_tilde_trimmed,Delta_t,file_u_new_trimmed)
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
            call p_interfaz(my_chem,dir_pb_trimmed,num_st_var,file_u_tilde_trimmed,Delta_t,file_u_new_trimmed)
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