!> Transient diffusion equation
subroutine initialise_diffusion_transient(this, root)
    use diffusion_transient_m, only: diffusion_1D_transient_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c
    use spatial_discr_m, only: spatial_discr_c
    use time_discr_m, only: time_discr_homog_c, time_discr_heterog_c, time_discr_c
    use diff_props_heterog_m, only: diff_props_heterog_1D_c
    use properties_m, only: props_c
    use BCs_m, only: BCs_1D_c
    use diff_stab_params_m, only: stab_params_diff_c
    use stability_parameters_m, only: stab_params_c
    use char_params_diff_m, only: char_params_diff_c
    !use eigenvalues_eigenvectors_m
    !use prob_dens_fcts_m
    implicit none

    !> Variables
    class(diffusion_1D_transient_c) :: this
    character(len=*), intent(in) :: root !> root name for output files

    class(props_c), pointer :: my_props=>null()
    type(diff_props_heterog_1D_c), target :: my_props_diff
    class(spatial_discr_c), pointer :: my_mesh=>null()
    type(spatial_discr_rad_c), target :: my_radial_mesh
    class(time_discr_c), pointer :: my_time_discr=>null()
    type(time_discr_homog_c), target :: my_homog_time_discr
    type(time_discr_heterog_c), target :: my_heterog_time_discr
    type(BCs_1D_c) :: my_BCs
    class(stab_params_c), pointer :: my_stab_params=>null()
    type(stab_params_diff_c), target :: my_stab_params_diff
    type(char_params_diff_c), target :: my_char_params
    
    real(kind=8) :: Final_time, measure, Delta_r_0, length, L2_norm_vi, a
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
    real(kind=8), allocatable :: c0(:),c_e(:),Delta_x_vec(:),g(:),b(:),y0(:),y(:), Time_out(:),Delta_r(:),Delta_r_D(:)
    real(kind=8), allocatable :: get_my_parameters(:), Dirichlet_BCs(:)
    real(kind=8), allocatable :: d(:),e(:)
    logical :: dimless
    real(kind=8), parameter :: tol=1d-5, pi=4d0*atan(1d0)
    character(len=200) :: filename
!****************************************************************************************************************************************************
!> Dimensionless form flag
    dimless=.false.
    this%dimless=dimless
!> Boundary conditions
    call my_BCs%read_BCs("BCs.dat")
    if (my_BCs%labels(1).eq.1 .and. my_BCs%labels(2).eq.1) then
        call my_BCs%read_Dirichlet_BCs_conc("Dirichlet_BCs.dat")
    end if
    call this%set_BCs_1D_trans(my_BCs)
 !> Uniform radial mesh
    !my_mesh=>my_radial_mesh
    allocate(spatial_discr_rad_c :: my_mesh)
    call my_mesh%read_mesh("Delta_r_homog.dat")
    call this%set_spatial_discr(my_mesh)
    Num_cells=this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
!> Uniform time discretisation
    !my_time_discr=>my_homog_time_discr
    allocate(time_discr_homog_c :: my_time_discr)
    allocate(Delta_t(1)) !> chapuza
    Delta_t=1d-1*minval(my_radial_mesh%Delta_r)**2
    Final_time=5d-1
    int_method=1                                                !> 1: Euler explicit
                                                                !> 2: Euler semi-implicit
                                                                !> 3: Euler fully implicit
                                                                !> 4: Crank-Nicolson
                                                                !> 5: RKF45
    call my_time_discr%set_Delta_t(Delta_t)
    call my_time_discr%set_int_method(int_method)
    call my_time_discr%set_Final_time(Final_time)
    call my_time_discr%compute_Num_time()
    call this%set_time_discr(my_time_discr)
!****************************************************************************************************************************************************
!> Diffusion properties
    call my_props_diff%read_props("diff_props.dat",this%spatial_discr)
    call my_props_diff%set_source_term_flag(this%BCs)
    call this%set_diff_props_heterog(my_props_diff)
!****************************************************************************************************************************************************
!> Stability parameters
    call my_stab_params_diff%compute_stab_params_diff_1D(this%diff_props_heterog,my_mesh,Delta_t(1))
    call this%set_stab_params_diff(my_stab_params_diff)
!****************************************************************************************************************************************************
!> External concentration
    allocate(c_e(Num_cells))
    c_e=0d0
    call this%conc_obj%set_conc_ext(c_e)
!> Initial concentration
    allocate(c0(Num_cells))
    c0=0d0
    call this%set_conc_init(c0)
!****************************************************************************************************************************************************
!> Matrices y vectors
    call this%allocate_trans_mat()
    call this%compute_trans_mat_PDE()
    call this%compute_source_term_PDE()
!****************************************************************************************************************************************************
!> Post-process
    nullify(my_mesh,my_time_discr)
end subroutine
