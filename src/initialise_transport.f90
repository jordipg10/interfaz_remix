!subroutine initialise_transport(this)
!>    use BCs_subroutines_m
!>    use transport_stab_params_homog_m
!>    use stab_params_viena_m
!>    use properties_viena_m
!>    use transport_m
!>    !use eigenvalues_eigenvectors_m
!>    !use prob_dens_fcts_m
!>    implicit none
!
!>    class(transport_1D_C) :: this
!>    
!>    class(props_c), pointer :: my_props=>null()
!>    type(tpt_props_homog_c), target :: my_props_homog
!>    type(tpt_props_heterog_c), target :: my_props_heterog
!>    type(props_viena_c), target :: my_props_viena
!>    class(spatial_discr_c), pointer :: my_mesh=>null()
!>    type(mesh_1D_Euler_homog_c), target :: my_homog_mesh
!>    type(mesh_1D_Euler_heterog_c), target :: my_heterog_mesh
!>    class(time_discr_c), pointer :: my_time_discr=>null()
!>    type(time_discr_homog_c), target :: my_homog_time_discr
!>    type(time_discr_heterog_c), target :: my_heterog_time_discr
!>    type(BCs_t) :: my_BCs
!>    class(stab_params_c), pointer :: my_stab_params=>null()
!>    type(stab_params_tpt_homog_c), target :: my_stab_params_homog
!>    type(stab_params_viena_c), target :: my_stab_params_viena
!>    type(char_params_tpt_c), target :: my_char_params
!>    
!>    real(kind=8) :: q0,Delta_t,Delta_x,theta
!>    real(kind=8), allocatable :: c0(:),c_e(:),source_term_vec(:),porosity_vec(:),dispersion_vec(:)
!>    real(kind=8), allocatable :: d(:),e(:)
!>    integer(kind=4) :: parameters_flag,i,Num_cells,Num_time,info,n
!>    real(kind=8), parameter :: pi=4d0*atan(1d0), epsilon_x=1d-2, epsilon_t=1d-4
!>    character(len=200) :: filename
!>    
!>    !external :: dsterf
!!****************************************************************************************************************************************************
!!> Boundary conditions
!>    !BCs=2
!>    !evap=.false. !> Evaporation
!>    !call my_BCs%set_BCs_label(BCs)
!>    !call my_BCs%set_evap(evap)
!>    call my_BCs%read_BCs("BCs.dat")
!>    if (my_BCs%BCs_label.eqv.1) then
!>        call my_BCs%read_Dirichlet_BCs("Dirichlet_BCs_tpt.dat")
!>    end if
!>    call this%set_BCs(my_BCs)
!> !> Mesh
!>    my_mesh=>my_homog_mesh
!>    !Num_targets=10
!>    !Delta_x=measure/Num_targets
!>    call my_homog_mesh%read_mesh("Delta_x_homog.dat")
!>    !call my_mesh%set_Num_targets(Num_targets)
!>    !call my_homog_mesh%compute_measure()
!>    !call my_mesh%set_Delta_x_homog(Delta_x)
!>    call this%set_spatial_discr(my_homog_mesh)
!>    Num_cells=this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
!>    !> Time discretisation
!>    !my_time_discr=>my_homog_time_discr
!>    !Delta_t=1d-4
!>    !Final_time=1d-1
!>    !call my_homog_time_discr%read_time_discr("Delta_t_homog.dat")
!>    !call my_homog_time_discr%set_Delta_t_homog(Delta_t)
!>    !call my_homog_time_discr%set_Final_time(Final_time)
!>    !call my_homog_time_discr%compute_Num_time()
!>    !call this%set_time_discr(my_homog_time_discr)
!>    !write(*,*) "Transient transport with homogeneous(1) or heterogeneous(2) parameters?"
!>    !read(*,*) parameters_flag
!>    parameters_flag=3
!>    if (parameters_flag.eqv.1) then
!>    !> Homogeneous properties
!>        !!> Mesh
!>        !my_mesh=>my_homog_mesh
!>        !!Num_targets=10
!>        !!Delta_x=measure/Num_targets
!>        !call my_homog_mesh%read_mesh("Delta_x_homog.txt")
!>        !!call my_mesh%set_Num_targets(Num_targets)
!>        !call my_homog_mesh%compute_measure()
!>        !!call my_mesh%set_Delta_x_homog(Delta_x)
!>        !call this%set_spatial_discr(my_homog_mesh)
!>        !Num_targets=this%spatial_discr%Num_targets
!>        !!> Time discretisation
!>        !my_time_discr=>my_homog_time_discr
!>        !!Delta_t=1d-4
!>        !!Final_time=1d-1
!>        !call my_homog_time_discr%read_time_discr("Delta_t_homog.txt")
!>        !!call my_homog_time_discr%set_Delta_t_homog(Delta_t)
!>        !!call my_homog_time_discr%set_Final_time(Final_time)
!>        !call my_homog_time_discr%compute_Num_time()
!>        !call this%set_time_discr(my_homog_time_discr)
!>        !> Initial properties
!>        my_props=>my_props_homog
!>        !allocate(source_term(this%spatial_discr%Num_targets))
!>        !source_term=0d0
!>        call my_props_homog%read_source_term("source_term_homog.dat",this%spatial_discr)
!>        call my_props_homog%set_source_term_flag(this%BCs)
!>        call my_props_homog%read_props("tpt_props_homog.dat")
!>        call my_props_homog%compute_retardo("retardo.dat")
!>        !print *, my_props_homog%retardo        
!>        call this%set_props(my_props_homog)
!>        !print *, this%props%source_term
!>        !> Stability parameters
!>        !my_stab_params=>my_stab_params_homog
!>        !call my_stab_params_homog%compute_stab_params(this%props,this%spatial_discr,this%time_discr)
!>        !call this%set_stab_params(my_stab_params_homog)
!>        !!> Characteristic parameters
!>        !call my_char_params%compute_char_params(this%props)
!>        !call this%set_char_params(my_char_params)
!>        !print *, this%props%get_props()
!>    else if (parameters_flag.eqv.2) then
!>    !> Heterogeneous properties
!>        !!> Mesh
!>        !my_mesh=>my_heterog_mesh
!>        !call my_heterog_mesh%read_mesh("Delta_x_heterog.txt")
!>        !call my_heterog_mesh%compute_measure()
!>        !call this%set_spatial_discr(my_heterog_mesh)
!>        !Num_targets=this%spatial_discr%Num_targets
!>        !!> Time discretisation
!>        !my_time_discr=>my_heterog_time_discr
!>        !call my_heterog_time_discr%read_time_discr("Delta_t_heterog.txt")
!>        !call my_heterog_time_discr%compute_Num_time()
!>        !call this%set_time_discr(my_heterog_time_discr)
!>        !> Initial properties
!>        my_props=>my_props_heterog
!>        allocate(source_term_vec(Num_cells))
!>        source_term_vec=1d1
!>        call my_props_heterog%set_source_term(source_term_vec)
!>        call my_props_heterog%set_source_term_flag(this%BCs)
!>        !print *, my_props
!>        allocate(porosity_vec(Num_cells),dispersion_vec(Num_cells))
!>        porosity_vec=1d0
!>        q0=2d1
!>        dispersion_vec=4d-1
!>        
!>        if (this%spatial_discr%scheme.eqv.1) then
!>            call my_props_heterog%interpolate_flux('flux_interfaces.dat',this%spatial_discr)
!>        else if (this%spatial_discr%scheme.eqv.2) then
!>            call my_props_heterog%compute_flux_lin(q0,this%spatial_discr)
!>        end if
!>        call my_props_heterog%set_tpt_props_heterog(porosity_vec,dispersion_vec) 
!>        call my_props_heterog%compute_retardo("retardo.dat")
!>        print *, my_props_heterog%flux
!>        call this%set_props(my_props_heterog)
!>        !> Stability parameters
!>        !my_stab_params=>my_stab_params_heterog
!>        !call my_stab_params_heterog%compute_stab_params(this%props,this%spatial_discr,this%time_discr)
!>        !call this%set_stab_params(my_stab_params_heterog)
!>        !print *, my_stab_params_heterog%alpha(1:5,1)*my_heterog_props%flux(1:5), my_stab_params_heterog%beta(1:5,1)
!>        !print *, this%props%get_props(3)
!>    else if (parameters_flag.eqv.3) then !> Viena
!>        my_props=>my_props_viena
!>        this%spatial_discr%Num_targets=1
!>        this%spatial_discr%targets_flag=0
!>        Num_cells=this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
!>        call my_props_viena%read_props('viena_parametros.dat')
!>        !print *, "Time weighting factor for mixing ratio?"
!>        !read(*,*) my_props_viena%theta
!>        !call my_props_viena%compute_lambda(this%time_discr%get_Delta_t(1),theta)
!>        call my_props_viena%compute_source_term_viena()
!>        call my_props_viena%set_source_term_flag(this%BCs)
!>        call this%set_props(my_props_viena)
!>        !my_stab_params=>my_stab_params_viena
!>        !call my_stab_params_viena%compute_stab_params(this%props,this%spatial_discr,this%time_discr)
!>        !call this%set_stab_params(my_stab_params_viena)
!>    else
!>        error stop "Parameters not defined yet"
!>    end if
!>    
!!****************************************************************************************************************************************************
!!> Critical mesh size
!>    !Delta_x=this%stab_params%Delta_x_crit
!>    !Delta_x=this%stab_params%Delta_x_crit+epsilon_x
!>    !Delta_x=this%stab_params%Delta_x_crit-epsilon_x
!>    !call my_homog_mesh%set_Delta_x_homog(Delta_x)
!>    !call my_homog_mesh%compute_measure()
!>    !call this%set_spatial_discr(my_homog_mesh)
!>    !call my_stab_params_homog%compute_stab_params(this%props,this%spatial_discr,this%time_discr)
!>    !call this%set_stab_params(my_stab_params_homog)
!!> Critical time step
!>    !Delta_t=this%stab_params%Delta_t_crit
!>    !Delta_t=this%stab_params%Delta_t_crit+epsilon_t
!>    !Delta_t=this%stab_params%Delta_t_crit-epsilon_t
!>    !Num_time=5
!>    !call my_homog_time_discr%set_Delta_t_homog(Delta_t)
!>    !call my_homog_time_discr%set_Num_time(Num_time)
!>    !call my_homog_time_discr%compute_Final_time()
!>    !call this%set_time_discr(my_homog_time_discr)
!!> Initial concentration
!>    !allocate(c0(Num_cells))
!>    !i=nint(Num_cells/2)
!>    !c0(1:i)=1
!>    !c0(i+1:Num_cells)=0
!>    !c0(1)=1d0
!>    !c0(2:Num_targets)=0d0
!>    !if (my_radial_mesh%Dirichlet_BC_location.eqv.0) then
!>    !>    c0(1)=g(1)
!>    !else
!>    !>    c0(Num_targets)=g(Num_targets)
!>    !end if
!>    !c0=1d0
!>    !!c0(i+1)=5d-1
!>    !c0(i+1:Num_targets-2)=1d0 !> c0: pulse injection
!>    !c0(Num_targets-1:Num_targets)=c0(1:i)
!>    !c0=1d0
!>    !sigma=my_homog_mesh%Delta_x*Num_targets/5d0 !> standard deviation
!>    !c0=Gaussian_pdf(sigma,my_homog_mesh)
!>    !call this%set_conc_init(c0)
!!> External concentration
!>    !allocate(c_e(Num_cells))
!>    !c_e=1d0
!>    !call this%set_conc_ext(c_e)
!>    !call this%set_conc_star_flag()
!!> Discretisation
!>    !write(*,*) "Scheme? (1:CFD, 2:EFD, 3:Upwind)"
!>    !read(*,*) scheme
!>    !if (parameters_flag.eqv.1) then
!>    !>    scheme=1
!>    !else
!>    !>    scheme=2
!>    !end if
!>    !call this%spatial_discr%set_scheme(scheme)
!>    !write(*,*) "Integration method? (1:Euler explicit, 2:Euler implicit, 3:Crank-Nicolson)"
!>    !read(*,*) int_method
!>    !int_method=1
!>    !call this%time_discr%set_int_method(int_method)
!!****************************************************************************************************************************************************
!!> Matrices
!>    call this%compute_trans_mat_PDE()
!>    print *, this%trans_mat%diag
!>    call this%compute_source_term_PDE()
!>    if (this%BCs%BCs_label.eqv.1) then
!>        call Dirichlet_BCs_PDE(this)
!>    else if (this%BCs%BCs_label.eqv.2) then
!>        call Neumann_homog_BCs(this)
!>    else
!>        error stop "Boundary conditions not implemented yet"
!>    end if
!>    !call this%compute_F_mat_PDE()
!>    !call this%mixing_ratios%allocate_array(Num_cells)
!!****************************************************************************************************************************************************
!!> Prueba Lapack
!>    !d=dble([1,2,3,4,5])
!>    !e=dble([2,4,1,3])
!>    !n=size(d)
!>    !call dsterf(n,d,e,info)
!>    !print *, d,e
!!****************************************************************************************************************************************************
!>    nullify(my_mesh,my_props)
!end subroutine