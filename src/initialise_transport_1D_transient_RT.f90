!> \file initialise_transport_1D_transient_RT.f90
!> \brief Initialize 1D transient transport object for reactive transport (RT) simulation
!>
!> \details
!> This subroutine initializes a 1D transient transport object specifically configured for 
!> reactive transport simulations. It performs the following operations:
!> 
!> **Initialization Steps:**
!> 1. Set dimensionless form flag
!> 2. Read and configure spatial discretization (mesh)
!>    - Option 1: 1D homogeneous Eulerian mesh
!>    - Option 2: 1D heterogeneous Eulerian mesh
!>    - Option 3: Radial mesh (for cylindrical/spherical coordinates)
!> 3. Read and set boundary conditions (BCs)
!>    - Dirichlet, Neumann, Robin, or prescribed mass flux (PMF)
!>    - Read inflow flux data
!> 4. Read and set time discretization parameters
!>    - Time step size Δt
!>    - Final time
!>    - Time integration scheme (Euler explicit/implicit, RKF45)
!> 5. Read and configure transport properties
!>    - Porosity φ
!>    - Dispersion coefficient D
!>    - Darcy flux q
!>    - Source/sink terms
!> 6. Compute flux field (constant, linear, or polynomial)
!> 7. Compute stability parameters (Courant number, Peclet number)
!>
!> **Mesh Types:**
!> - mesh_type = 1: 1D homogeneous (uniform properties)
!> - mesh_type = 2: 1D heterogeneous (spatially variable properties)
!> - mesh_type = 3: Radial (cylindrical/spherical symmetry)
!>
!> **Flux Computation:**
!> - source_term_order = 0: Constant flux
!> - source_term_order > 0: Polynomial flux (coefficients from file)
!>
!> \param[in,out] this Transport 1D transient object to be initialized
!> \param[in] path Directory path for input and output files
!> \param[in] root Root name for input and output files
!> \param[in] mesh_type Mesh type selector: 1 = 1D homogeneous, 2 = 1D heterogeneous, 3 = radial
!>
!> \use BCs_m Provides BCs_1D_c type for boundary conditions
!> \use time_discr_m Provides time discretization classes
!> \use transport_stab_params_m Provides stability parameters for transport
!> \use spatial_discr_1D_m Provides 1D spatial discretization classes
!> \use spatial_discr_rad_m Provides radial spatial discretization
!> \use transport_transient_m Provides transport_1D_transient_c class
!> \use transport_properties_heterog_m Provides heterogeneous transport properties
!> \use char_params_m Provides characteristic parameters
!> \use char_params_tpt_m Provides transport characteristic parameters
!> \use vectors_m Provides vector utilities (infinity norm)

subroutine initialise_transport_1D_transient_RT(this,path,root,mesh_type)
    use spatial_discr_m, only: spatial_discr_c !> Import spatial discretization base class
    use BCs_m, only: BCs_1D_c !> Import boundary conditions type
    use time_discr_m, only: time_discr_homog_c, time_discr_heterog_c, time_discr_c !> Import time discretization classes (homogeneous, heterogeneous, base class)
    use transport_stab_params_m, only: stab_params_tpt_1D_c !> Import stability parameters type for transport
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c !> Import 1D mesh classes (homogeneous, heterogeneous, base class)
    use spatial_discr_2D_m, only: mesh_2D_Euler_homog_c !> Import 2D homogeneous Eulerian mesh class (for potential future use)
    use spatial_discr_rad_m, only: spatial_discr_rad_c !> Import radial spatial discretization class
    use transport_transient_m, only: transport_1D_transient_c !> Import 1D transient transport class
    use transport_properties_heterog_m, only: tpt_props_heterog_1D_c !> Import heterogeneous transport properties class
    use char_params_m, only: char_params_c !> Import characteristic parameters base class
    use char_params_tpt_m, only: char_params_tpt_c !> Import transport characteristic parameters class
    use vectors_m, only: inf_norm_vec_real !> Import infinity norm function for real vectors
    implicit none !> Enforce explicit variable declaration (no implicit typing)

    class(transport_1D_transient_c) :: this !> Transport 1D transient object to be initialized (contains mesh, time discretization, BCs, properties)
    character(len=*), intent(in) :: path !> Directory path for input and output files (e.g., "./input/")
    character(len=*), intent(in) :: root !> Root name for input and output files (e.g., "problem1")
    integer(kind=4), intent(in) :: mesh_type !> Mesh type selector: 1 = 1D homogeneous, 2 = 1D heterogeneous, 3 = radial [-]
    
    type(tpt_props_heterog_1D_c) :: my_props_tpt !> Heterogeneous transport properties object (porosity, dispersion, flux)
    class(spatial_discr_c), pointer :: my_mesh=>null() !> Pointer to spatial discretization base class (points to specific mesh type)
    type(mesh_1D_Euler_homog_c), target :: my_homog_mesh !> 1D homogeneous Eulerian mesh object (uniform properties)
    type(mesh_2D_Euler_homog_c), target :: my_2D_homog_mesh !> 2D homogeneous Eulerian mesh object (for potential future use)
    type(mesh_1D_Euler_heterog_c), target :: my_heterog_mesh !> 1D heterogeneous Eulerian mesh object (spatially variable properties)
    type(spatial_discr_rad_c), target :: my_rad_mesh !> Radial spatial discretization object (cylindrical/spherical symmetry)
    class(time_discr_c), pointer :: my_time_discr=>null() !> Pointer to time discretization base class (points to specific time discretization type)
    type(time_discr_homog_c), target :: my_homog_time_discr !> Homogeneous time discretization object (uniform time steps)
    type(time_discr_heterog_c), target :: my_heterog_time_discr !> Heterogeneous time discretization object (variable time steps)
    type(BCs_1D_c) :: my_BCs !> Boundary conditions object (contains BC types, values, and flux data)
    type(stab_params_tpt_1D_c) :: my_stab_params_tpt !> Transport stability parameters object (Courant number, Peclet number, critical time step)
    type(char_params_tpt_c) :: my_char_params_tpt !> Transport characteristic parameters object (dimensionless numbers)
    
    real(kind=8) :: q0,Delta_t,Delta_x,theta,measure,Final_time,x_1,x_2,x_3,x,phi !> Auxiliary real variables: q0=flux [L/T], Δt=time step [T], Δx=space step [L], θ=time weighting [-], measure [L], Final_time [T], x positions [L], φ=porosity [-]
    real(kind=8), allocatable :: c0(:),c_e(:),source_term_vec(:),porosity_vec(:),dispersion_vec(:),flux_vec(:),flux_coeffs(:) !> Allocatable arrays: c0=initial concentration [C], c_e=equilibrium concentration [C], source_term_vec [C/T], porosity_vec [-], dispersion_vec [L²/T], flux_vec [L/T], flux_coeffs (polynomial coefficients)
    real(kind=8), allocatable :: d(:),e(:) !> Auxiliary arrays for matrix construction: d=diagonal, e=off-diagonal
    integer(kind=4) :: parameters_flag,i,Num_cells,Num_time,info,n,adapt_ref_flag,scheme,int_method,r_flag,flux_ord,half_num_tar !> Integer variables: flags, indices, counters, array sizes [-]
    real(kind=8), parameter :: pi=4d0*atan(1d0), eps=1d-12, epsilon_x=1d-2, epsilon_t=1d-6 !> Mathematical and tolerance constants: π, machine epsilon, spatial tolerance, temporal tolerance [-]
    character(len=200) :: filename,mesh_file !> Character strings for file names
    logical :: evap,dimless !> Logical flags: evap=evaporation flag, dimless=dimensionless form flag
!****************************************************************************************************************************************************
!> Dimensionless form flag
    dimless=.false. !> Set dimensionless form flag to FALSE (use dimensional equations, should be read from input file)
    this%dimless=dimless !> Store dimensionless flag in transport object
!> Mesh (chapuza)
    mesh_file=trim(path//root//'_discr_esp.dat') !> Construct mesh filename by concatenating path, root, and '_discr_esp.dat' (spatial discretization file)
    if (mesh_type.eq.1) then !> Check if mesh type is 1 (1D homogeneous)
        allocate(mesh_1D_Euler_homog_c :: my_mesh) !> Allocate pointer to 1D homogeneous Eulerian mesh class
    else if (mesh_type.eq.2) then !> Check if mesh type is 2 (1D heterogeneous)
        allocate(mesh_1D_Euler_heterog_c :: my_mesh) !> Allocate pointer to 1D heterogeneous Eulerian mesh class
    else if (mesh_type.eq.3) then !> Check if mesh type is 3 (radial)
        allocate(spatial_discr_rad_c :: my_mesh) !> Allocate pointer to radial spatial discretization class
        mesh_file=trim(path//root//'_discr_esp_rad.dat') !> Override mesh filename for radial case (use '_discr_esp_rad.dat')
        !write(*,*) "Porosity?" !> Print prompt for porosity input (interactive)
        !read(*,*) phi !> Read porosity value from user input [-]
        phi=1d-4 !> Set porosity to 1d-4 (hardcoded workaround for push-pull test) [-]
    else if (mesh_type.eq.4) then !> Check if mesh type is 4 (2D homogeneous)
        allocate(mesh_2D_Euler_homog_c :: my_mesh) !> Allocate pointer to 2D homogeneous spatial discretization class
    else !> If mesh type is not 1, 2, or 3
        error stop "Error: mesh_type must be 1 (1D homog), 2 (1D heterog), 3 (radial) or 4 (2D homog)" !> Terminate program with error message
    end if !> End mesh type selection
    call my_mesh%read_mesh(mesh_file,phi) !> Read mesh data from file (node positions, cell sizes, connectivity, porosity φ)
    call this%set_spatial_discr(my_mesh) !> Store mesh object in transport object's spatial discretization component
!> Boundary conditions
    call my_BCs%read_BCs(path//root//'_BCs.dat') !> Read boundary condition types and values from file (Dirichlet, Neumann, Robin, PMF)
    !if (my_BCs%labels(1).eq.1 .and. my_BCs%labels(2).eq.1 .and. this%spatial_discr%targets_flag.eq.0) then !> (commented out) Check if both BCs are Dirichlet and no internal targets
        call my_BCs%read_caudal_inf(path//root//'_caudal_inf.dat') !> Read inflow flux (caudal) data from file q(t) [L³/T or L/T depending on dimension]
    ! else if (my_BCs%labels(1).eq.3) then !> (commented out) Check if left BC is Robin type
    !     call my_BCs%read_Robin_BC_inflow(path//root//'_Robin_BC_inflow.dat') !> (commented out) Read Robin BC inflow data from file
    !end if !> (commented out) End BC type check
    call this%set_BCs_1D_trans(my_BCs) !> Store boundary conditions object in transport object
    call this%BCs%compute_flux_inf(this%spatial_discr) !> Compute inflow flux value(s) from boundary condition data [L/T]
!> Time discretization setup
    !my_time_discr=>my_homog_time_discr !> Alternative: use pointer assignment instead of allocation (commented out)
    allocate(time_discr_homog_c :: my_time_discr) !> Allocate pointer to homogeneous time discretization class (uniform time steps) [-]
    call my_time_discr%read_time_discr(path//root//'_discr_temp.dat') !> Read time discretization from file: initial time t₀, final time t_f, time step Δt [T]
    call this%set_time_discr(my_time_discr) !> Assign time discretization object to transport%time_discr member [-]
!****************************************************************************************************************************************************
!> Transport properties configuration
    call my_props_tpt%read_props(path//root,this%spatial_discr) !> Read transport properties from file: porosity φ, dispersion coefficient D, boundary fluxes [L], [L²/T], [L/T]
    if (my_props_tpt%source_term_order.eq.0) then !> Check if source term is constant (order 0 polynomial) [-]
        if (my_props_tpt%cst_flux_flag .eqv. .true.) then !> Check if Darcy flux is spatially constant [L/T]
            call this%BCs%set_cst_flux_boundary(my_props_tpt%flux_int(1)) !> Set constant flux boundary condition with integrated flux value [L³/T]
        else !> Flux is spatially variable (non-constant)
            select type (mesh=>this%spatial_discr) !> Polymorphic type selection to determine flux computation method
            class is (spatial_discr_rad_c) !> If mesh is radial (cylindrical/spherical symmetry)
                call my_props_tpt%compute_flux_rad(mesh,this%BCs%caudal_inf) !> Compute radial flux field q(r) from flow rate and radial geometry [L/T]
            class default !> If mesh is 1D Cartesian (homogeneous or heterogeneous)
                call my_props_tpt%compute_flux_lin(this%BCs%flux_inf,mesh) !> Compute linear flux field q(x) = q_inf (uniform Darcy flux) [L/T]
            end select
            !print *, my_props_tpt%flux !> Debugging output (commented out): display computed flux field
        end if
    else if (my_props_tpt%source_term_order>0) then !> Source term is polynomial (order > 0): f(x) = ∑ᵢ aᵢxⁱ
    !> Temporary workaround solution (should be refactored)
        open(unit=2,file=path//root//"_flux_coeffs.dat",status='old',action='read') !> Open file containing polynomial flux coefficients [-]
        read(2,*) flux_ord !> Read polynomial order (degree + 1, number of coefficients) [-]
        allocate(flux_coeffs(flux_ord+1)) !> Allocate array for polynomial coefficients [varies by term]
        read(2,*) flux_coeffs !> Read polynomial coefficients: [a₀, a₁, ..., aₙ] for flux polynomial
        close(2) !> Close flux coefficients file
        call my_props_tpt%set_source_term_order(flux_ord-1) !> Set source term polynomial degree (order = degree) [-]
        call my_props_tpt%compute_flux_nonlin(flux_coeffs,this%spatial_discr) !> Compute nonlinear (polynomial) flux field from coefficients q(x) = ∑ᵢ aᵢxⁱ [L/T]
        call my_props_tpt%compute_source_term(this%spatial_discr,flux_coeffs) !> Compute source term field from flux: f = -∇·q = -dq/dx [C/T or 1/T]
    end if
    call my_props_tpt%set_source_term_flag(this%BCs) !> Set flag indicating presence of source/sink term based on boundary conditions [-]
    call my_props_tpt%compute_dispersion_1D(this%spatial_discr%scheme) !> Compute dispersion at cell interfaces from cell-centered values using numerical scheme [L²/T]
    call this%set_tpt_props_heterog_obj(my_props_tpt) !> Assign transport properties object to transport%tpt_props_heterog member [-]
!****************************************************************************************************************************************************
!> Stability parameters computation
    call my_stab_params_tpt%compute_stab_params_tpt_1D(this%tpt_props_heterog,my_mesh,my_time_discr%get_Delta_t()) !> Compute Courant number Co = qΔt/Δx and Peclet number Pe = qΔx/D for stability analysis [-]
    call this%set_stab_params_tpt(my_stab_params_tpt) !> Assign stability parameters object to transport%stab_params_tpt member [-]
!**************************************************************************************************************************************************
    nullify(my_mesh,my_time_discr) !> Nullify pointers to avoid dangling references (memory managed by transport object)
end subroutine