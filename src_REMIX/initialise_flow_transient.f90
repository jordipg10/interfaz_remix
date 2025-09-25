!> Transient flow equation
subroutine initialise_flow_transient(this,root)
    use flow_transient_m, only: flow_transient_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c, spatial_discr_c
    use time_discr_m, only: time_discr_homog_c, time_discr_heterog_c, time_discr_c
    use flow_props_heterog_m, only: flow_props_heterog_c, props_c, flow_props_heterog_conf_c
    use BCs_m, only: BCs_t
    use stab_params_flow_m, only: stab_params_flow_c, stab_params_c
    use char_params_flow_m, only: char_params_flow_c, char_params_c
    use matrices_m, only: tridiag_matrix_c
    use BCs_subroutines_m, only: Dirichlet_BCs_PDE, Neumann_homog_BCs, Robin_Neumann_homog_BCs
    implicit none

    !> Variables
    class(flow_transient_c) :: this
    character(len=*), intent(in) :: root
    
    
    class(props_c), pointer :: my_props=>null()
    type(flow_props_heterog_conf_c), target :: my_props_flow
    class(spatial_discr_c), pointer :: my_mesh=>null()
    type(spatial_discr_rad_c), target :: my_radial_mesh
    class(time_discr_c), pointer :: my_time_discr=>null()
    type(time_discr_homog_c), target :: my_homog_time_discr
    type(time_discr_heterog_c), target :: my_heterog_time_discr
    type(BCs_t) :: my_BCs
    type(stab_params_flow_c) :: my_stab_params
    class(char_params_c), pointer :: my_char_params=>null()
    type(char_params_flow_c), target :: my_char_params_flow
    
    real(kind=8) :: Final_time, measure, Delta_r_0, length, L2_norm_vi, a, h_c
    real(kind=8) :: porosity, flux, velocity, dispersion, q0,rho,Dirichlet_BC
    real(kind=8) :: retardo
    real(kind=8), allocatable :: source_term(:),porosity_vec(:), flux_vec(:), velocity_vec(:), dispersion_vec(:)
    real(kind=8) :: theta,sigma
    real(kind=8), allocatable :: A_sub(:),A_diag(:),A_super(:),A_mat_Dirichlet(:,:),T_sub(:),T_diag(:),T_super(:),lambda(:),&
    lambda_star(:),eigenvectors(:,:),eigenvectors_star(:,:),z0(:),A_lambda(:,:,:),A_lambda_v(:)
    real(kind=8) :: beta,courant,Gamma_sub,Gamma_diag,Gamma_super,radius
    integer(kind=4) :: scheme,int_method,Num_time,opcion_BCs,parameters_flag,i,j,Num_cells,eqn_flag,dim,niter_max,niter,&
    Dirichlet_BC_location,nK,targets_flag,info
    integer(kind=4) :: BCs(2)    
    real(kind=8), allocatable :: K(:,:),Delta_t(:)
    real(kind=8), allocatable :: h0(:),c_e(:),Delta_x_vec(:),g(:),b(:),y0(:),y(:), Time_out(:),Delta_r(:),Delta_r_D(:)
    real(kind=8), allocatable :: get_my_parameters(:), Dirichlet_BCs(:)
    real(kind=8), allocatable :: d(:),e(:)
    logical :: dimless
    real(kind=8), parameter :: tol=1d-5, pi=4d0*atan(1d0)
    character(len=200) :: filename
    
    type(tridiag_matrix_c) :: E_mat
!****************************************************************************************************************************************************
    call this%set_sol_method(1) !> numerical solution method
!> Dimensionless form flag
    dimless=.true.
    call this%set_dimless_flag(dimless)
!> Boundary conditions
    call my_BCs%read_BCs(root//"_BCs.dat")
    if (my_BCs%BCs_label(1).eq.1 .and. my_BCs%BCs_label(2).eq.1) then
        call my_BCs%read_Dirichlet_BCs_head(root//"_Dirichlet_BCs.dat")
    end if
    call this%set_BCs(my_BCs)
 !> Radial mesh
    !my_mesh=>my_radial_mesh
    allocate(spatial_discr_rad_c :: my_mesh)
    call my_mesh%read_mesh(root//"_discr_esp_rad.dat")
    !if (this%dimless) then
    !    call my_mesh%compute_dimless_mesh()
    !end if
    call this%set_spatial_discr(my_mesh)
    Num_cells=this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
!> time discretisation
    !my_time_discr=>my_homog_time_discr
    allocate(time_discr_homog_c :: my_time_discr)
    !allocate(Delta_t(1)) !> chapuza
    !Delta_t=1d-1*minval(my_radial_mesh%Delta_r)**2
    !Final_time=5d-1
    !int_method=1                                                !> 1: Euler explicit
    !                                                            !> 2: Euler semi-implicit
    !                                                            !> 3: Euler fully implicit
    !                                                            !> 4: Crank-Nicolson
    !                                                            !> 5: RKF45
    !call my_time_discr%set_Delta_t(Delta_t)
    !call my_time_discr%set_int_method(int_method)
    !call my_time_discr%set_Final_time(Final_time)
    !call my_time_discr%compute_Num_time()
    call my_time_discr%read_time_discr(root//"_discr_temp.dat")
    call this%set_time_discr(my_time_discr)
!****************************************************************************************************************************************************
!> flow properties
    call my_props_flow%read_props(root,this%spatial_discr)
    call my_props_flow%set_source_term_flag(this%BCs)
    call this%set_flow_props_heterog(my_props_flow)
!****************************************************************************************************************************************************
!> Stability parameters
    call my_stab_params%compute_stab_params(this%flow_props_heterog,my_mesh,this%time_discr%get_Delta_t())
    call this%set_stab_params_flow(my_stab_params)
!****************************************************************************************************************************************************
!> Characteristic parameters
    !allocate(char_params_flow_c :: my_char_params)
    my_char_params=>my_char_params_flow
    !!> chapuza
    !if (this%BCs%head_inf>this%BCs%head_out) then
    !    h_c=this%BCs%head_inf
    !else
    !    h_c=this%BCs%head_out
    !end if
    call my_char_params_flow%set_char_head(this%BCs)
    call my_char_params_flow%compute_char_params(this%flow_props_heterog,my_mesh)
    call this%set_char_params_flow(my_char_params_flow)
    if (this%dimless) then
        call this%spatial_discr%compute_dimless_mesh(my_char_params_flow%char_length)
        call this%time_discr%compute_dimless_Delta_t(my_char_params_flow%char_time)
        call this%time_discr%compute_Final_time_D(my_char_params_flow%char_time)
        call this%stab_params_flow%compute_Delta_t_D_crit(my_char_params_flow%char_time)
        call this%BCs%compute_dimless_BCs(this%char_params_flow%char_head)
    end if
!****************************************************************************************************************************************************
!> Initial concentration
    call this%allocate_head()
    allocate(h0(Num_cells))
    h0=0d0 !> Initialize head to zero
    call this%set_head_init(h0)
!****************************************************************************************************************************************************
!> We compute arrays
    call this%allocate_arrays_PDE_1D()
    call this%compute_trans_mat_PDE()
    call this%compute_rech_mat_PDE()
    call this%compute_source_term_PDE()
    call this%compute_F_mat_PDE()
!> We impose BCs
    if (this%BCs%BCs_label(1).eq.1 .and. this%BCs%BCs_label(2).eq.1) then
        call Dirichlet_BCs_PDE(this)
    else if (this%BCs%BCs_label(1).eq.2 .and. this%BCs%BCs_label(2).eq.2) then
        call Neumann_homog_BCs(this)
    else if (this%BCs%BCs_label(1).eq.3 .and. this%BCs%BCs_label(2).eq.2) then
        call Robin_Neumann_homog_BCs(this)
    else
        error stop "Boundary conditions not implemented yet"
    end if
!> We compute arrays for linear system
    if (this%time_discr%int_method.eq.1) then !> Euler explicit
        theta=0d0
    else if (this%time_discr%int_method.eq.2) then !> Euler fully implicit
        theta=1d0
    else if (this%time_discr%int_method.eq.3) then !> Crank-Nicolson
        theta=0.5
    end if
    call this%compute_E_mat(E_mat)
    call this%compute_X_mat(theta,E_mat)
    call this%compute_A_mat(theta,E_mat)
    call this%compute_Y_mat()
    call this%compute_Z_mat()
    call this%compute_f_vec()
!****************************************************************************************************************************************************
!> Post-process
    nullify(my_mesh,my_time_discr)
end subroutine
