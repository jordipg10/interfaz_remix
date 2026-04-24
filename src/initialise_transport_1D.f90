subroutine initialise_transport_1D(this,dir,root)
    use BCs_subroutines_m, only: Dirichlet_BCs_PDE, Neumann_homog_BCs, Robin_Neumann_homog_BCs
    use transport_m, only: transport_1D_c
    use transport_properties_heterog_m, only: tpt_props_heterog_1D_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c
    use spatial_discr_m, only: spatial_discr_c
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c
    use BCs_m, only: BCs_1D_c
    implicit none

    class(transport_1D_c) :: this
    character(len=*), intent(in) :: dir !> root name for input files
    character(len=*), intent(in) :: root !> root name for input files
    
    type(tpt_props_heterog_1D_c) :: my_props_tpt
    class(spatial_discr_c), pointer :: my_mesh=>null()
    type(mesh_1D_Euler_homog_c), target :: my_homog_mesh
    type(mesh_1D_Euler_heterog_c), target :: my_heterog_mesh
    type(BCs_1D_c) :: my_BCs
    
    real(kind=8) :: q0,Delta_t,Delta_x,theta,measure
    real(kind=8), allocatable :: c0(:),c_e(:),source_term_vec(:),porosity_vec(:),dispersion_vec(:),flux_vec(:),Dirichlet_BCs(:)
    real(kind=8), allocatable :: d(:),e(:),flux_coeffs(:)
    integer(kind=4) :: parameters_flag,i,Num_cells,Num_time,info,n,adapt_ref_flag,BCs,opcion_BCs,scheme,flux_ord
    real(kind=8), parameter :: pi=4d0*atan(1d0), epsilon_x=1d-2, epsilon_t=1d-4
    character(len=200) :: filename
    logical :: evap,dimless
!****************************************************************************************************************************************************
!> Dimensionless form flag
    dimless=.false.
    this%dimless=dimless
!> Boundary conditions
    call my_BCs%read_BCs(dir//root//"_BCs.dat")
    if (my_BCs%labels(1).eq.1 .and. my_BCs%labels(2).eq.1) then
        call my_BCs%read_Dirichlet_BCs_conc(dir//root//"_Dirichlet_BCs.dat")
        call my_BCs%read_caudal_inf(dir//root//"_flux_inflow.dat")
    else if (my_BCs%labels(1).eq.3) then
        call my_BCs%read_Robin_BC_inflow(dir//root//"_Robin_BC_inflow.dat")
    end if
    call this%set_BCs_1D(my_BCs)
 !> Mesh
    my_mesh=>my_homog_mesh
    call my_mesh%read_mesh(dir//root//"_Delta_x_homog.dat")
    call this%set_spatial_discr(my_mesh)
    Num_cells=this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
!> Transport properties
    call my_props_tpt%read_props(dir//root//"_tpt_props.dat",this%spatial_discr)
    if (my_props_tpt%source_term_order.eq.0) then
        call my_props_tpt%compute_flux_lin(this%BCs%flux_inf,this%spatial_discr)
    else
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
    call this%set_tpt_props_heterog_obj(my_props_tpt)
!****************************************************************************************************************************************************
!> External concentration
    allocate(c_e(Num_cells))
    c_e=0d0
    call this%set_conc_ext(c_e)
    call this%set_conc_r_flag()
!****************************************************************************************************************************************************
!> Matrices y vectors
    call this%allocate_trans_mat()
    call this%compute_trans_mat_PDE()
    call this%compute_source_term_PDE()
!> Condiciones contorno
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
    nullify(my_mesh)
end subroutine