program main_interfaz
    use RT_1D_m, only: RT_1D_transient_c, transport_1D_transient_c, chemistry_c
    implicit none
!> Objects
    type(RT_1D_transient_c) :: my_RT_trans !> 1D transient reactive transport class
    type(transport_1D_transient_c) :: my_tpt_trans !> 1D transient transport class
    type(chemistry_c) :: my_chem !> chemistry class
!> Variables
    integer(kind=4) :: num_comps !> number of components
    integer(kind=4) :: flag !> flag for reactive mixing loop
    integer(kind=4) :: option_tpt !> flag for lambdas
    integer(kind=4) :: int_method_chem !> integration method for chemical reactions
    real(kind=8) :: Delta_t !> time step
    character(len=256) :: file_in !> input file
    character(len=256) :: file_out !> output file
    character(len=256) :: path_inp !> path for input files
    character(len=:), allocatable :: path_inp_trim !> path for input files trimmed
    character(len=256) :: path_DB !> path for databases
    character(len=:), allocatable :: path_DB_trim !> path for databases trimmed
    character(len=256) :: path_pb !> path of problem to solve
    character(len=:), allocatable :: path_pb_trim !> path of problem to solve trimmed
    character(len=256) :: root_pb !> root of problem to solve
    character(len=:), allocatable :: root_pb_trim !> root of problem to solve trimmed
!****************************************************************************************************************************
    write(*,*) 'Esta es la interfaz para resolver la mecla reactiva con el WMA usando Euler explícito'
!> Name of path containing chemical and transport information
    write(*,*) "Directorio bases de datos:"
    read(*,*) path_DB !> must be written by the user
    path_DB_trim = trim(path_DB)
!> Choose problem
    write(*,*) "Directorio del problema que se quiere resolver:"
    read(*,*) path_pb !> must be written by the user
    path_pb_trim=trim(path_pb)
    write(*,*) "Root del problema que se quiere resolver (p.ej. yeso_eq)"
    read(*,*) root_pb !> must be written by the user
    root_pb_trim=trim(root_pb)
!> Initialise transport
    write(*,*) "Elegir opcion transporte (0: calcular lambdas, 1: leer lambdas):"
    read(*,*) option_tpt !> must be written by the user
    if (option_tpt.eq.0) then !> compute lambdas
    !> we read transport data, BCs and discretisations
        !> in the explicit case, we also compute stability parameters
            call my_tpt_trans%initialise_transport_1D_transient_RT(path_pb_trim//root_pb_trim)
        !> we allocate transport arrays
            call my_tpt_trans%allocate_arrays_PDE_1D()
            call my_tpt_trans%allocate_conc()
        !> we compute transport arrays, including mixing ratios, and we impose BCs
            call my_tpt_trans%compute_mixing_ratios_Delta_t_homog() !> missing the case of heterogenous time steps
        !> we set transport attribute in reactive transport object
            call my_RT_trans%set_transport_trans(my_tpt_trans)
        !> we choose and set integration method for chemical reactions
            !! 1: Euler explicit, 2: Euler fully implicit, 3: Crank-Nicolson
            !write(*,*) "Elegir metodo integracion temporal para reacciones quimicas (1: Euler explicito, 2: Euler totalmente implicito, 3: Crank-Nicolson):"
            !read(*,*) int_method_chem !> must be written by the user
            int_method_chem = 1 !> we set default method for chemical reactions
            call my_RT_trans%set_int_method_chem_reacts(int_method_chem)
    else if (option_tpt.eq.1) then !> read lambdas
    !> we set transport attribute in reactive transport object
        call my_RT_trans%set_transport_trans(my_tpt_trans) !> esto es un create en realidad
    !> we read temporal discretisation
        call my_RT_trans%read_time_discretisation(path_pb_trim//root_pb_trim)
    !> we read transport data for WMA
        call my_RT_trans%transport%read_transport_data_WMA(path_pb_trim//root_pb_trim)
    else
        error stop "Esta opcion no esta implementada todavia"
    end if
!> we read chemistry
    call my_chem%read_chemistry(path_pb_trim//root_pb_trim,path_DB_trim)
    write(*,*) 'Introduzca los siguientes datos:'
    write(*,*) 'Número de componentes acuosas:'
    read(*,*) num_comps !> read number of components
    write(*,*) 'Paso de tiempo (asumimos que es uniforme):'
    read(*,*) Delta_t !> read time step
    write(*,*) 'Nombre del fichero de lectura (concentraciones despues de resolver transporte conservativo). IMPORTANTE: las filas deben representar especies acuosas y las columnas targets.'
    read(*,*) file_in !> read input file name
    write(*,*) 'Nombre del fichero de escritura (concentraciones despues de resolver la mezcla reactiva):'
    read(*,*) file_out !> read output file name
    flag=1
    do while (flag==1)
        call my_chem%interfaz_comps_arch(path_pb_trim,num_comps,file_in,Delta_t,file_out)
        write(*,*) 'Introduzca 1 para iterar, cualquier otra tecla para salir del bucle:'
        read(*,*) flag !> read flag to continue or exit loop
        !file_in=file_out !> we set the input file for the next iteration to be the output file of the previous one
        !print *, file_in, file_out
    end do
end program main_interfaz