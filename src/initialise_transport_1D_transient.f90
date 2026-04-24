subroutine initialise_transport_1D_transient(this,dir,root)
    use spatial_discr_m, only: spatial_discr_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c
    use time_discr_m, only: time_discr_homog_c, time_discr_heterog_c, time_discr_c
    use vectors_m, only: inf_norm_vec_real
    use transport_transient_m, only: transport_1D_transient_c
    use transport_stab_params_m, only: stab_params_tpt_1D_c
    use transport_properties_heterog_m, only: tpt_props_heterog_1D_c
    use char_params_tpt_m, only: char_params_tpt_c
    use BCs_m, only: BCs_1D_c
    use BCs_subroutines_m, only: Dirichlet_BCs_PDE, Neumann_homog_BCs, Robin_Neumann_homog_BCs
    use arrays_m, only: tridiag_matrix_c
    implicit none

    class(transport_1D_transient_c) :: this
    character(len=*), intent(in) :: dir !> root name for input files
    character(len=*), intent(in) :: root
    !character(len=*), intent(in) :: path
    !character(len=*), intent(in) :: file_BCs
    !character(len=*), intent(in) :: file_spatial_discr
    !character(len=*), intent(in) :: file_time_discr
    !character(len=*), intent(in) :: file_tpt_props
    
    type(tpt_props_heterog_1D_c) :: my_props_tpt
    class(spatial_discr_c), pointer :: my_mesh=>null()
    type(mesh_1D_Euler_homog_c), target :: my_homog_mesh
    type(mesh_1D_Euler_heterog_c), target :: my_heterog_mesh
    class(time_discr_c), pointer :: my_time_discr=>null()
    type(time_discr_homog_c), target :: my_homog_time_discr
    type(time_discr_heterog_c), target :: my_heterog_time_discr
    type(BCs_1D_c) :: my_BCs
    type(stab_params_tpt_1D_c) :: my_stab_params_tpt
    type(char_params_tpt_c) :: my_char_params_tpt
    
    type(tridiag_matrix_c) :: E_mat,E_mat_prev
    
    real(kind=8) :: q0,Delta_x,theta,measure,Final_time,x_1,x_2,x_3,x
    real(kind=8), allocatable :: c0(:),c_e(:),source_term_vec(:),porosity_vec(:),dispersion_vec(:),flux_vec(:),flux_coeffs(:)
    real(kind=8), allocatable :: d(:),e(:),Delta_t(:)
    integer(kind=4) :: parameters_flag,i,Num_cells,Num_time,info,n,adapt_ref_flag,scheme,int_method,r_flag,flux_ord,half_num_tar
    real(kind=8), parameter :: pi=4d0*atan(1d0), eps=1d-12, epsilon_x=1d-2, epsilon_t=1d-6
    character(len=200) :: filename
    logical :: evap,dimless
!****************************************************************************************************************************************************
!> Dimensionless form flag
    dimless=.false. !> esto habria que leerlo
    this%dimless=dimless
!> Mesh (chapuza)
    !my_mesh=>my_homog_mesh
    !allocate(mesh_1D_Euler_homog_c :: my_mesh) !> chapuza
    allocate(spatial_discr_rad_c :: my_mesh) !> chapuza
    !call my_mesh%read_mesh(root//'_discr_esp.dat')
    call my_mesh%read_mesh(dir//root//"_discr_esp_rad.dat")
    call this%set_spatial_discr(my_mesh)
!> Boundary conditions
    call my_BCs%read_BCs(dir//root//"_BCs.dat")
    if (my_BCs%labels(1).eq.1 .and. my_BCs%labels(2).eq.1 .and. this%spatial_discr%targets_flag.eq.0) then
        call my_BCs%read_Dirichlet_BCs_conc(dir//root//"_Dirichlet_BCs.dat")
        call my_BCs%read_caudal_inf(dir//root//"_flux_inflow.dat")
    else if (my_BCs%labels(1).eq.3) then
        call my_BCs%read_Robin_BC_inflow(dir//root//"_Robin_BC_inflow.dat")
    end if
    call this%set_BCs_1D_trans(my_BCs)
    call this%BCs%compute_flux_inf(this%spatial_discr)
!> Uniform time discretisation
    !my_time_discr=>my_homog_time_discr
    allocate(time_discr_homog_c :: my_time_discr)
    call my_time_discr%read_time_discr(dir//root//"_discr_temp.dat")
    call this%set_time_discr(my_time_discr)
!****************************************************************************************************************************************************
!> Transport properties
    call my_props_tpt%read_props(dir//root//"_tpt_props.dat",this%spatial_discr)
    if (my_props_tpt%source_term_order.eq.0) then !> constant source term
        if (my_props_tpt%cst_flux_flag .eqv. .true.) then !> flux is constant
            call this%BCs%set_cst_flux_boundary(my_props_tpt%flux_int(1))
        else !> flux is non-constant
            select type (mesh=>this%spatial_discr)
            class is (spatial_discr_rad_c)
                call my_props_tpt%compute_flux_rad(mesh,this%BCs%caudal_inf)
            class default
                call my_props_tpt%compute_flux_lin(this%BCs%flux_inf,mesh)
            end select
            !print *, my_props_tpt%flux
        end if
    else if (my_props_tpt%source_term_order>0) then !> flux is polynomic
    !> chapuza
        open(unit=2,file=dir//root//"_flux_coeffs.dat",status='old',action='read')
        read(2,*) flux_ord
        allocate(flux_coeffs(flux_ord+1))
        read(2,*) flux_coeffs
        close(2)
        call my_props_tpt%set_source_term_order(flux_ord-1)
        call my_props_tpt%compute_flux_nonlin(flux_coeffs,this%spatial_discr)
        call my_props_tpt%compute_source_term(this%spatial_discr,flux_coeffs)
    end if
    call my_props_tpt%set_source_term_flag(this%BCs)
    call my_props_tpt%compute_dispersion_1D(this%spatial_discr%scheme)
    call this%set_tpt_props_heterog_obj(my_props_tpt)
!****************************************************************************************************************************************************
!> Stability parameters
    call my_stab_params_tpt%compute_stab_params_tpt_1D(this%tpt_props_heterog,my_mesh,my_time_discr%get_Delta_t())
    call this%set_stab_params_tpt(my_stab_params_tpt)
!****************************************************************************************************************************************************
!> Critical time step test
    !select type (time=>this%time_discr)
    !type is (time_discr_homog_c)
    !>    call time%set_Delta_t_homog(this%stab_params_tpt%Delta_t_crit-epsilon_t)
    !>    call time%compute_Num_time()
    !>    call this%stab_params_tpt%compute_stab_params(this%tpt_props_heterog,my_homog_mesh%Delta_x,time%Delta_t)
    !end select
!****************************************************************************************************************************************************
!> External concentration
    allocate(c_e(this%spatial_discr%Num_targets))
    c_e=0d0
    call this%conc_obj%set_conc_ext(c_e)
    call this%set_conc_r_flag()
!> Initial concentration
    allocate(c0(this%spatial_discr%Num_targets))
    half_num_tar=nint(this%spatial_discr%Num_targets/2d0)
    if (this%BCs%labels(1)<3 .and. this%BCs%labels(2)<3) then
        c0(1:half_num_tar)=1d0
        c0(half_num_tar+1:half_num_tar)=0d0
        !c0=0d0
    else
        c0=0d0
    end if
    if (this%spatial_discr%targets_flag.eq.1) then
        if (this%BCs%labels(1).eq.1) then !> Dirichlet inflow
            c0(1)=this%BCs%conc_inf
        end if
        if (this%BCs%labels(2).eq.1) then !> Dirichlet outflow
            c0(Num_cells)=this%BCs%conc_out
        end if
    end if
    call this%set_conc_init(c0)
!****************************************************************************************************************************************************
!> We compute transport arrays
    call this%allocate_arrays_PDE()
    call this%compute_trans_mat_PDE()
    call this%compute_rech_mat_PDE()
    call this%compute_source_term_PDE()
    call this%compute_F_mat_PDE()
    !> We compute arrays for linear system
    if (this%time_discr%int_method.eq.1) then !> Euler explicit
        theta=0d0
    else if (this%time_discr%int_method.eq.2) then !> Euler fully implicit
        theta=1d0
    else if (this%time_discr%int_method.eq.3) then !> Crank-Nicolson
        theta=0.5
    end if
    call this%compute_E_mat_1D(this%time_discr%get_Delta_t(),E_mat,E_mat_prev)
    call this%compute_X_mat_1D(theta,E_mat_prev)
    call this%compute_A_mat_1D(theta,E_mat)
    call this%compute_Y_mat()
    call this%compute_Z_mat()
    call this%compute_f_vec(this%time_discr%get_Delta_t())
!> We impose BCs
    if (this%BCs%labels(1).eq.1 .and. this%BCs%labels(2).eq.1) then
        call Dirichlet_BCs_PDE(this)
    else if (this%BCs%labels(1).eq.2 .and. this%BCs%labels(2).eq.2) then
        call Neumann_homog_BCs(this)
    else if (this%BCs%labels(1).eq.3 .and. this%BCs%labels(2).eq.2) then
        call Robin_Neumann_homog_BCs(this)
    else
        error stop "Boundary conditions not implemented yet"
    end if
!****************************************************************************************************************************************************
    nullify(my_mesh,my_time_discr)
end subroutine