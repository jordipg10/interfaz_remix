!> @file main_interfaz.f90
!> @brief Interactive driver for 1D reactive mixing iterations using the WMA
!>        (Water Mixing Approach) with explicit Euler time integration.
!> @details Reads chemistry/database paths and problem inputs from the user,
!>          loads the chemical system, writes initial/external water type
!>          component concentrations, and runs a loop of reactive mixing
!>          iterations using either a constant or variable time step.
program main_interfaz
    !> Import chemistry container type and the two interface procedures
    !> (one for systems with equilibrium reactions, one without).
    use chemistry_m, only: chemistry_c, interfaz_comps_arch, interfaz_esp_arch
    !> Standard intrinsic module providing IEEE arithmetic helpers.
    use, intrinsic :: ieee_arithmetic
    !> Standard intrinsic module providing IEEE floating-point exception flags.
    use, intrinsic :: ieee_exceptions
    !> Disable implicit typing: every variable must be declared explicitly.
    implicit none
    !> Objects
    type(chemistry_c) :: my_chem !>< Chemistry container holding chemical system, waters and reactive zones.
    !> Variables
    integer(kind=4) :: int_method_chem !>< Integration method for chemical reactions (1: Euler explicit, 2: Euler fully implicit, 3: Crank-Nicolson).
    integer(kind=4) :: num_aq_comps !>< Number of aqueous components in the chemical system.
    integer(kind=4) :: flag !>< Loop continuation flag (1: continue, 0: exit).
    integer(kind=4) :: flag_Delta_t !>< Whether the time step is constant (1) or variable (0).
    integer(kind=4), allocatable :: ind_can_vec(:) !>< Indices of canonical vectors in mixing ratios.
    integer(kind=4), allocatable :: ind_non_can_vec(:) !>< Indices of non-canonical vectors in mixing ratios.
    real(kind=8), allocatable :: u_tilde_init(:,:) !>< Component concentrations u_tilde at the initial time step.
    real(kind=8) :: Delta_t !>< Current time step value.
    character(len=100) :: dir_DB,dir_pb, root_files,file_u_tilde,file_u_new,file_u_wat_types !>< Raw fixed-length directory and file name inputs.
    character(len=:), allocatable :: dir_DB_trimmed, dir_pb_trimmed, root_files_trimmed, file_u_tilde_trimmed, file_u_new_trimmed,&
        file_u_wat_types_trimmed !>< Trimmed deferred-length versions of the input strings.
    integer :: ios !>< I/O status code returned by safe read statements.
    logical :: has_ieee !>< True if the runtime supports the IEEE intrinsic standard.
    procedure(interfaz_comps_arch), pointer :: p_interfaz=>null() !>< Procedure pointer dispatching to the appropriate reactive-mixing interface.
    !> Pre-Process
    !> Query whether IEEE intrinsic support is available in this build.
    has_ieee = ieee_support_standard()
    !> Configure IEEE behaviour and clear all relevant exception flags.
    if (has_ieee) then
        !> Disable gradual underflow (use abrupt underflow to zero).
        call ieee_set_underflow_mode(gradual=.false.)
        !> Clear underflow exception flag.
        call ieee_set_flag(ieee_underflow, .false.)
        !> Clear inexact-result exception flag.
        call ieee_set_flag(ieee_inexact, .false.)
        !> Clear overflow exception flag.
        call ieee_set_flag(ieee_overflow, .false.)
        !> Clear divide-by-zero exception flag.
        call ieee_set_flag(ieee_divide_by_zero, .false.)
        !> Clear invalid-operation exception flag.
        call ieee_set_flag(ieee_invalid, .false.)
    end if
    !> Process
    !> Banner describing what the interactive program does.
    write(*,*) "Esta es la interfaz para resolver la parte reactiva de una iteración de mezcla reactiva en 1D usando &
        Euler explícito."
    !> Prompt the user for the database directory.
    write(*,*) "Directorio de las bases de datos: "
    !> Read the database directory from standard input.
    read(*,*, iostat=ios) dir_DB
    !> If the read failed (e.g. EOF when running non-interactively), abort gracefully.
    if (ios /= 0) then
        write(*,*) 'Error/EOF leyendo dir_DB. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
    end if
    !> Trim trailing blanks from the database directory string.
    dir_DB_trimmed = trim(dir_DB)
    !> Prompt the user for the problem directory.
    write(*,*) "Directorio del problema: "
    !> Read the problem directory from standard input.
    read(*,*, iostat=ios) dir_pb
    !> Abort gracefully on read error or EOF.
    if (ios /= 0) then
        write(*,*) 'Error/EOF leyendo dir_pb. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
    end if
    !> Trim trailing blanks from the problem directory string.
    dir_pb_trimmed = trim(dir_pb)
    !> Prompt the user for the common root of input/output filenames.
    write(*,*) "Root de los archivos de entrada y salida: "
    !> Read the file root from standard input.
    read(*,*, iostat=ios) root_files
    !> Abort gracefully on read error or EOF.
    if (ios /= 0) then
        write(*,*) 'Error/EOF leyendo root_files. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
    end if
    !> Trim trailing blanks from the file root string.
    root_files_trimmed = trim(root_files)
    !> we read chemistry
    !> Load the chemical system, water types and reactive zones from disk.
    call my_chem%read_chemistry_interface(root_files_trimmed,dir_pb_trimmed,dir_DB_trimmed)
    !> write file name
    !> Prompt for the output file name for initial/external water-type concentrations.
    write(*,*) "Nombre del archivo donde quieres que escriba las concentraciones de las componentes acuosas &
        de los tipos de agua iniciales y externas?"
    !> Read the output filename for water-type concentrations.
    read(*,*, iostat=ios) file_u_wat_types !> file with initial and external water types
    !> Abort gracefully on read error or EOF.
    if (ios /= 0) then
        write(*,*) 'Error/EOF leyendo file_u_wat_types. Ejecuta en una terminal interactiva o redirige desde fort.5'; call & 
            safe_stop(1)
    end if
    !> Trim trailing blanks from the water-types output filename.
    file_u_wat_types_trimmed = trim(file_u_wat_types)
    !> we write concentrations of initial and external water types
    !> Write the initial and external water-type component concentrations to the chosen file.
    call my_chem%write_conc_comp_wat_types(dir_pb_trimmed,file_u_wat_types_trimmed)
    !> write file name 
    !> Prompt for the output file name for post-reactive-mixing concentrations.
    write(*,*) "Nombre del archivo donde quieres que escriba las concentraciones despues de la mezcla reactiva?"
    !> Read the output filename for post-mixing concentrations.
    read(*,*, iostat=ios) file_u_new
    !> Abort gracefully on read error or EOF.
    if (ios /= 0) then
        write(*,*) 'Error/EOF leyendo file_u_new. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
    end if
    !> Trim trailing blanks from the post-mixing output filename.
    file_u_new_trimmed = trim(file_u_new)
    !> Prompt for the input file containing u_tilde (post conservative-transport concentrations).
    write(*,*) "Nombre del archivo que contiene las concentraciones de componentes acuosas despues de resolver una iteracion de &
        transporte conservativo? & 
        IMPORTANTE: El archivo debe estar en el directorio del problema, el numero de filas tiene que ser el numero &
        de componentes y el numero de columnas el numero de targets."
    !> Read the input filename providing u_tilde.
    read(*,*, iostat=ios) file_u_tilde !> file with u_tilde
    !> Abort gracefully on read error or EOF.
    if (ios /= 0) then
        write(*,*) 'Error/EOF leyendo file_u_tilde. Ejecuta en una terminal interactiva o redirige desde fort.5'; call & 
            safe_stop(1)
    end if
    !> Prompt for the initial time step value.
    write(*,*) "Paso de tiempo inicial?"
    !> Read the initial time step.
    read(*,*, iostat=ios) Delta_t !> initial time step
    !> Validate the read: report I/O failure or non-positive time step.
    if (ios /= 0) then
        write(*,*) 'Error/EOF leyendo Delta_t. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
    else if (Delta_t <= 0d0) then
        !> Time step must be strictly positive.
        error stop "El paso de tiempo debe ser mayor que cero"
    end if
    !> Prompt the user to choose constant or variable time step mode.
    write(*,*) "Paso de tiempo constante? (1: si, 0: no)"
    !> Read the constant/variable time step flag.
    read(*,*, iostat=ios) flag_Delta_t !> whether time step is constant (1) or variable (0)
    !> Abort gracefully on read error or EOF.
    if (ios /= 0) then
        write(*,*) 'Error/EOF leyendo flag_Delta_t. Ejecuta en una terminal interactiva o redirige desde fort.5'; call & 
            safe_stop(1)
    end if
    !> Trim trailing blanks from the u_tilde input filename.
    file_u_tilde_trimmed = trim(file_u_tilde) !> we trim file name
    !> We get number of aqueous components in the first target water,
    !> which is needed to choose the interface and to read the input file with u_tilde.
    !> Query the first target water for its number of aqueous components.
    num_aq_comps=my_chem%get_num_aq_comps_tar_wat()
    !> Echo the number of aqueous components for diagnostics.
    print *, "Numero de componentes acuosas: ", num_aq_comps
    !> We choose interface based on whether there are equilibrium reactions or not, using procedure pointers
    !> If the first target water has no equilibrium reactions, use the species-based interface.
    if (my_chem%waters(my_chem%tar_wat_indices(1))%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions==0) then
        !> Bind the procedure pointer to the no-equilibrium variant.
        p_interfaz=>interfaz_esp_arch !> no equilibrium reactions
    else
        !> Bind the procedure pointer to the components-based variant (with equilibrium).
        p_interfaz=>interfaz_comps_arch !> with equilibrium reactions
    end if
    !> Inform the user that the reactive mixing loop is about to start.
    write(*,*) "Procedemos al bucle de mezcla reactiva. Tendrás que actualizar el archivo con las concentraciones & 
        obtenidas despues de resolver el transporte conservativo en cada iteración."
    !> Branch on time-step mode: constant vs variable.
    if (flag_Delta_t.eq.1) then !> constant time step
        !> Constant-time-step iteration loop.
        do
            !> Perform one reactive mixing iteration using the selected interface.
            call p_interfaz(my_chem,dir_pb_trimmed,num_aq_comps,file_u_tilde_trimmed,Delta_t,file_u_new_trimmed)
            !> Ask the user whether to perform another iteration.
            write(*,*) "Quieres hacer otra iteración de mezcla reactiva? (1: si, 0: no)"
            !> Read the loop continuation flag.
            read(*,*, iostat=ios) flag !> loop flag
            !> Abort gracefully on read error or EOF.
            if (ios /= 0) then
                write(*,*) 'Error/EOF leyendo flag. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
            end if
            !> Exit on 0; reject any value other than 0 or 1.
            if (flag.eq.0) then
                exit
            else if (flag.ne.1) then
                error stop "Opción no válida. Tiene que ser 1 o 0."
            end if
        end do
    else if (flag_Delta_t.eq.0) then !> variable time step
        !> Variable-time-step iteration loop.
        do
            !> Perform one reactive mixing iteration using the selected interface.
            call p_interfaz(my_chem,dir_pb_trimmed,num_aq_comps,file_u_tilde_trimmed,Delta_t,file_u_new_trimmed)
            !> Ask the user whether to perform another iteration.
            write(*,*) "Quieres hacer otra iteración de mezcla reactiva? (1: si, 0: no)"
            !> Read the loop continuation flag.
            read(*,*, iostat=ios) flag !> loop flag
            !> Abort gracefully on read error or EOF.
            if (ios /= 0) then
                write(*,*) 'Error/EOF leyendo flag. Ejecuta en una terminal interactiva o redirige desde fort.5'; call safe_stop(1)
            end if
            !> Exit on 0, reject invalid values, otherwise prompt for a new time step.
            if (flag.eq.0) then
                exit
            else if (flag.ne.1) then
                error stop "Opción no válida. Tiene que ser 1 o 0."
            else
                !> Prompt the user for the next time step.
                write(*,*) "Nuevo paso de tiempo?"
                !> Read the new time step.
                read(*,*, iostat=ios) Delta_t !> new time step
                !> Validate the read: report I/O failure or non-positive time step.
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
    !> Clear any pending IEEE flags before normal termination.
    if (has_ieee) call clear_ieee_flags()
    !> Notify the user that execution finished successfully.
    write(*,*) "El programa ha terminado."
contains
    !> @brief Clear all IEEE floating-point exception flags if IEEE support is available.
    !> @details Used both before normal termination and inside @ref safe_stop to avoid
    !>          spurious end-of-program warnings about pending floating-point exceptions.
    subroutine clear_ieee_flags()
        !> Only act if the IEEE intrinsic standard is supported by the runtime.
        if (ieee_support_standard()) then
            !> Clear underflow flag.
            call ieee_set_flag(ieee_underflow, .false.)
            !> Clear inexact flag.
            call ieee_set_flag(ieee_inexact, .false.)
            !> Clear overflow flag.
            call ieee_set_flag(ieee_overflow, .false.)
            !> Clear divide-by-zero flag.
            call ieee_set_flag(ieee_divide_by_zero, .false.)
            !> Clear invalid-operation flag.
            call ieee_set_flag(ieee_invalid, .false.)
        end if
    end subroutine clear_ieee_flags

    !> @brief Stop the program with a given exit code, first clearing IEEE flags.
    !> @param[in] code Integer exit/stop code passed to the Fortran @c stop statement.
    subroutine safe_stop(code)
        !> Exit code propagated to the operating system.
        integer, intent(in) :: code
        !> Clear pending IEEE flags so the runtime does not emit warnings.
        call clear_ieee_flags()
        !> Terminate execution with the requested code.
        stop code
    end subroutine safe_stop
end program main_interfaz