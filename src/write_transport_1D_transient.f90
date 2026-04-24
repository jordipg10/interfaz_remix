!> \file write_transport_1D_transient.f90
!> \brief Writes comprehensive output for 1D transient transport simulation results
!> \details This subroutine writes all transport equation data and results to a formatted output file,
!> including:
!> - Mesh and discretization parameters (mesh size, volume, number of targets)
!> - Boundary conditions (Dirichlet, Neumann, Robin)
!> - Spatial discretization scheme (Traditional FD, Petchamé-Carrera, Upwind)
!> - Time discretization (time step, final time, number of steps)
!> - Integration method (Euler explicit/implicit, Crank-Nicolson, RKF45)
!> - Transport properties (porosity, flux, dispersion)
!> - Stability parameters (critical time step, Courant number, Peclet number)
!> - System matrices (A, B) and source vector (f)
!> - Concentration results at specified output times
!> - MRMT (Multi-Rate Mass Transfer) zone data if applicable
!> 
!> The governing equation is: F*dc/dt = T*c + g, discretized as A*c^(k+1) = B*c^k + f
!> 
!> Output format uses scientific notation (ES15.5) for real values and includes stability
!> warnings when Courant or Peclet conditions are violated.

!> \brief Write 1D transient transport simulation data and results to output file
!> \details Comprehensive output subroutine that writes all transport equation information including
!> discretization parameters, transport properties, stability analysis, system matrices, and
!> concentration solutions at specified output times. Handles both Eulerian and MRMT formulations.
!> 
!> \param[in] this transport_1D_transient_c object containing all transport simulation data
!> \param[in] root Character string root name for output files
!> \param[in] Time_out Real array of output time points [T]
!> \param[in] output Real 2D array of concentration results [C] (rows=cells, columns=time points)
!> 
!> \use transport_transient_m Transport transient module providing transport_1D_transient_c class
!> \use spatial_discr_1D_m Spatial discretization module for 1D mesh definitions
!> \use time_discr_m Time discretization module for time stepping information
!> Writes data and results of 1D transient transport equation
subroutine write_transport_1D_transient(this,root,Time_out,output)
    use transport_transient_m, only: transport_1D_transient_c !> Import 1D transient transport class
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c !> Import 1D Eulerian homogeneous mesh class
    use time_discr_m, only: time_discr_homog_c !> Import homogeneous time discretization class
    implicit none
    !> Variables
    class(transport_1D_transient_c), intent(in) :: this !> 1D transient transport object containing all simulation data
    character(len=*), intent(in) :: root !> root name for output files
    real(kind=8), intent(in) :: Time_out(:) !> array of output time points [T]
    real(kind=8), intent(in) :: output(:,:) !> concentration output array [C] (rows=spatial cells, columns=time points)

    integer(kind=4) :: Num_output !> number of output time points [-]
    real(kind=8), allocatable :: stab_params(:) !> stability parameters array (unused)
    integer(kind=4) :: i,j,k,n,n_flux !> loop indices and counters [-]
    character(len=256) :: file_out !> output file name string
    real(kind=8) :: Delta_t=0d0 !> time step size [T]
    
    n=this%spatial_discr%Num_targets !> get number of spatial targets (cells or interfaces) [-]
    Num_output=size(Time_out) !> get number of output time points [-]
    !if (this%spatial_discr%scheme.eq.1) then
    !    if (this%dimless.eqv..true.) then
    !        write(file_out,root//'transport_1D_trans_CFDS_adim.out')
    !    else
    !        write(file_out,"('transport_1D_trans_CFDS.out')")
    !    end if
    !else if (this%spatial_discr%scheme.eq.2) then
    !    if (this%dimless.eqv..true.) then
    !        write(file_out,"('transport_1D_trans_IFDS_adim.out')")
    !    else
    !        write(file_out,root//'transport_1D_trans_IFDS.out')
    !    end if
    !end if
    open(unit=1,file=root//'.out',status='unknown') !> open output file with root name
    write(1,"(2x,'Equation:',5x,'F*dc/dt = T*c + g',/)") !> write governing equation header
    write(1,"(2x,'Mesh volume:',F15.5/)") this%spatial_discr%measure !> write total mesh volume [L³]
    write(1,"(2x,'Number of targets:',I5/)") this%spatial_discr%Num_targets!-this%spatial_discr%targets_flag !> write number of spatial targets [-]
    if (this%spatial_discr%targets_flag.eq.0) then !> check if targets are cell-centered
        write(1,"(2x,'Targets type: Cells',/)") !> write cell-centered targets label
    else !> targets are interface-based
        write(1,"(2x,'Targets type: Interfaces',/)") !> write interface targets label
    end if
    select type (mesh=>this%spatial_discr) !> polymorphic type selection for mesh
    type is (mesh_1D_Euler_homog_c) !> if mesh is 1D Eulerian homogeneous
        write(1,"(2x,'Mesh size:',F15.5/)") mesh%Delta_x !> write uniform cell size Δx [L]
    end select
    if (this%BCs%labels(1).eq.1 .and. this%BCs%labels(2).eq.1) then !> check if both boundaries are Dirichlet
        write(1,"(2x,'Boundary conditions:',10x,'Dirichlet',10x,2F15.5,/)") this%BCs%conc_inf, this%BCs%conc_out !> write Dirichlet BC concentrations [C]
    else if (this%BCs%labels(1).eq.2 .and. this%BCs%labels(2).eq.2) then !> check if both boundaries are Neumann
        write(1,"(2x,'Boundary conditions:',10x,'Neumann homogeneous',/)") !> write Neumann homogeneous BC label
    else if (this%BCs%labels(1).eq.3 .and. this%BCs%labels(2).eq.2) then !> check if Robin inflow and Neumann outflow
        write(1,"(2x,'Boundary conditions:',10x,'Robin inflow, Neumann homogeneous outflow',/)") !> write mixed BC label
    end if
    if (this%spatial_discr%scheme.eq.1) then !> check if traditional finite difference scheme
        write(1,"(2x,'Scheme:',10x,'Traditional FD',/)") !> write traditional FD scheme label
    else if (this%spatial_discr%scheme.eq.2) then !> check if Petchamé-Carrera scheme
        write(1,"(2x,'Scheme:',10x,'Petcham�-Carrera (2024)',/)") !> write Petchamé-Carrera scheme label
    else if (this%spatial_discr%scheme.eq.3) then !> check if upwind scheme
        write(1,"(2x,'Scheme:',10x,'Upwind',/)") !> write upwind scheme label
    end if
    select type (time=>this%time_discr) !> polymorphic type selection for time discretization
    type is (time_discr_homog_c) !> if time discretization is homogeneous
        Delta_t=this%time_discr%get_Delta_t() !> get uniform time step size Δt [T]
        write(1,"(2x,'Time step:',ES15.5/)") Delta_t !> write time step Δt [T]
    end select
    write(1,"(2x,'Final time:',ES15.5/)") this%time_discr%Final_time !> write simulation final time t_f [T]
    write(1,"(2x,'Number of time steps:',I10/)") this%time_discr%Num_time !> write number of time steps N_t [-]
    if (this%time_discr%int_method.eq.1) then !> check if Euler explicit method
        write(1,"(2x,'Integration method:',10x,'Euler explicit',/)") !> write Euler explicit label
    else if (this%time_discr%int_method.eq.2) then !> check if Euler fully implicit method
        write(1,"(2x,'Integration method:',10x,'Euler fully implicit',/)") !> write Euler fully implicit label
    else if (this%time_discr%int_method.eq.3) then !> check if Crank-Nicolson method
        write(1,"(2x,'Integration method:',10x,'Crank-Nicolson',/)") !> write Crank-Nicolson label
    else if (this%time_discr%int_method.eq.4) then !> check if Runge-Kutta-Fehlberg method
        write(1,"(2x,'Integration method:',10x,'RKF45',/)") !> write RKF45 adaptive method label
    end if
    write(1,"(2x,'Linear system:',5x,'A*c^(k+1) = B*c^k + f',/)") !> write discretized linear system equation
    if (this%dimless.eqv..false.) then !> check if dimensional form (not dimensionless)
        if (this%tpt_props_heterog%homog_flag.eqv..true.) then !> check if transport properties are homogeneous
            write(1,"(2x,'Transport properties:',10x,'Homogeneous',/)") !> write homogeneous properties label
        else !> transport properties are heterogeneous
            write(1,"(2x,'Transport properties:',10x,'Heterogeneous',/)") !> write heterogeneous properties label
        end if
        !write(1,"(2x,'Properties:'/)")
        write(1,"(10x,'Porosity:',/)") !> write porosity section header
        do i=1,n-this%spatial_discr%targets_flag !> loop over spatial cells
            write(1,"(20x,ES15.5)") this%tpt_props_heterog%porosity(i) !> write porosity φ [-]
        end do
        write(1,"(/,10x,'Flux at interfaces:',/)") !> write flux section header
        do i=1,size(this%tpt_props_heterog%flux_int) !> loop over interfaces
            write(1,"(20x,ES15.5)") this%tpt_props_heterog%flux_int(i) !> write Darcy flux q [L/T]
        end do
        write(1,"(10x,'Dispersion at interfaces:',/)") !> write dispersion section header
        do i=1,size(this%tpt_props_heterog%disp_int) !> loop over interfaces
            write(1,"(20x,ES15.5)") this%tpt_props_heterog%disp_int(i) !> write dispersion coefficient D [L²/T]
        end do
    end if
    if (this%tpt_props_heterog%homog_flag.eqv..true.) then !> check if properties are homogeneous for stability analysis
        write(1,"(/,2x,'Stability parameters:',/)") !> write stability parameters section header
        write(1,"(10x,'Critical time step:',/)") !> write critical time step subsection header
        write(1,"(20x,ES15.5)") this%stab_params_tpt%Delta_t_crit !> write critical time step Δt_crit [T]
        if (Delta_t>0d0 .and. Delta_t>this%stab_params_tpt%Delta_t_crit .and. this%time_discr%int_method<3) then !> check if time step violates stability for explicit methods
            write(1,"(/,20x,'You must reduce time step to have stability')") !> write stability warning message
        end if
        write(1,"(/,10x,'Courant:',/)") !> write Courant number header
        write(1,"(20x,ES15.5)") this%stab_params_tpt%Courant !> write Courant number Co = q*Δt/(φ*Δx) [-]
        if (this%stab_params_tpt%Courant>1d0) then !> check if Courant condition is violated
            write(1,"(/,20x,'Courant condition violated',/)") !> write Courant violation warning
        end if
        write(1,"(/,10x,'Cell Peclet:',/)") !> write Peclet number header
        write(1,"(20x,ES15.5)") this%stab_params_tpt%cell_Peclet !> write Peclet number Pe = q*Δx/D [-]
        if (this%stab_params_tpt%cell_Peclet>2d0) then !> check if Peclet condition is violated
            write(1,"(/,20x,'Peclet condition violated',/)") !> write Peclet violation warning
        end if
    end if
    !write(1,"(2x,'F:'/)") 
    !do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
    !>    write(1,"(2x,F15.5)") this%F_mat%diag(i)
    !end do
    !write(1,"(/,2x,'Transition matrix T (with BCs):'/)") 
    !write(1,"(17x,2F15.5)") this%trans_mat%diag(1), this%trans_mat%super(1)    
    !do i=2,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag-1
    !>    write(1,"(2x,3F15.5)") this%trans_mat%sub(i-1), this%trans_mat%diag(i), this%trans_mat%super(i)
    !end do
    !write(1,"(2x,2F15.5/)") this%trans_mat%sub(this%spatial_discr%Num_targets-1), this%trans_mat%diag(this%spatial_discr%Num_targets)
    write(1,"(/,2x,'Matrix A:',/)") !> write matrix A section header
    write(1,"(17x,2F15.5)") this%A_mat%diag(1), this%A_mat%super(1) !> write first row of tridiagonal matrix A (diagonal and superdiagonal) [-]
    do i=2,this%spatial_discr%Num_targets-1 !> loop over interior rows of matrix A
        write(1,"(2x,3F15.5)") this%A_mat%sub(i-1), this%A_mat%diag(i), this%A_mat%super(i) !> write row i of matrix A (subdiagonal, diagonal, superdiagonal) [-]
    end do
    write(1,"(2x,2F15.5/)") this%A_mat%sub(this%spatial_discr%Num_targets-1), this%A_mat%diag(this%spatial_discr%Num_targets) !> write last row of matrix A (subdiagonal and diagonal) [-]
    write(1,"(/,2x,'Matrix B:',/)") !> write matrix B section header
    write(1,"(17x,2F15.5)") this%X_mat%diag(1), this%X_mat%super(1) !> write first row of tridiagonal matrix B (diagonal and superdiagonal) [-]
    do i=2,this%spatial_discr%Num_targets-1 !> loop over interior rows of matrix B
        write(1,"(2x,3F15.5)") this%X_mat%sub(i-1), this%X_mat%diag(i), this%X_mat%super(i) !> write row i of matrix B (subdiagonal, diagonal, superdiagonal) [-]
    end do
    write(1,"(2x,2F15.5/)") this%X_mat%sub(this%spatial_discr%Num_targets-1), this%X_mat%diag(this%spatial_discr%Num_targets) !> write last row of matrix B (subdiagonal and diagonal) [-]
    write(1,"(/,2x,'Vector f:',/)") !> write source/boundary vector f section header
    do i=1,this%spatial_discr%Num_targets !> loop over all spatial targets
        write(1,"(2x,F15.5)") this%f_vec(i) !> write source/boundary vector component f_i [C/T]
    end do
   
    if (size(output,1).eq.this%spatial_discr%Num_targets-this%spatial_discr%targets_flag) then !> check if output is for single mobile zone (not MRMT)
        if (this%time_discr%int_method<4) then !> check if using fixed time step methods (not adaptive RKF45)
            write(1,"(/,2x,'Cell',*(ES20.5)/)") (Time_out(k), k=1,Num_output) !> write header with output times [T]
            do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag !> loop over spatial cells
                write(1,"(2x,I4,*(F20.5))") i,(output(i,k), k=1,Num_output) !> write cell index and concentration at each output time [C]
            end do
        else !> using adaptive RKF45 method
            write(1,"(/,2x,'Cell',2ES20.5/)") Time_out(1), this%time_discr%Final_time !> write header with initial and final times [T]
            do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag !> loop over spatial cells
                write(1,"(2x,I5,2ES20.10)") i, this%conc_init(i), this%conc_obj%conc(i) !> write cell index, initial concentration, and final concentration [C]
            end do
        end if
    else !> output is for MRMT (Multi-Rate Mass Transfer) with mobile and immobile zones
        write(1,"(/,2x,'Mobile zone:',/)") !> write mobile zone section header
        write(1,"(10x,'Cell',3ES20.5/)") (Time_out(k), k=1,Num_output) !> write mobile zone header with output times [T]
        do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag !> loop over spatial cells
            write(1,"(10x,I4,3F20.5)") i,(output(i,k), k=1,Num_output) !> write cell index and mobile zone concentration at each output time [C]
        end do
        write(1,"(/,2x,'Immobile zones:',/)") !> write immobile zones section header
        do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag !> loop over spatial cells for immobile zones
           !> write(1,"(10x,I4,3F20.5)") i,(output(this%spatial_discr%Num_targets-this%spatial_discr%targets_flag+i,k), k=1,Num_output) !> would write immobile zone concentrations [C]
        end do
    end if
    rewind(1) !> rewind file to beginning
    close(1) !> close output file
end subroutine write_transport_1D_transient