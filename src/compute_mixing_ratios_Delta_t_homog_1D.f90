!> \file compute_mixing_ratios_Delta_t_homog_1D.f90
!> \brief Computes mixing ratios matrix for reactive transport with uniform time stepping
!> \details This subroutine computes the mixing ratios that describe how target waters mix
!> with boundary waters, neighboring waters, and recharge waters during one time step.
!> For implicit methods, it computes and inverts the system matrix A to obtain mixing ratios.
!> For explicit methods, it directly uses the transport matrix X.
!>
!> Mathematical formulation for implicit methods:
!> \f[
!> \text{Solve: } A \cdot c^{n+1} = X \cdot c^n + Z \cdot c_{bd} + Y \cdot c_{rech} + f
!> \f]
!> where mixing ratios are extracted from \f$ A^{-1} \f$ and the right-hand side matrices.
!>
!> \param[in,out] this Transport object containing spatial/temporal discretization and BCs
!> \param[out] A_mat_lumped Optional lumped system matrix for mass balance computations [-]

subroutine compute_mixing_ratios_Delta_t_homog_1D(this)
    use BCs_subroutines_m !> Module with boundary condition subroutines (Dirichlet, Neumann, Robin, mixed BCs)
    use transport_transient_m, only: transport_1D_transient_c !> Import transport class and matrix types
    use arrays_m, only: diag_matrix_c,tridiag_matrix_c, copy_real_array !> Import matrix classes for diagonal and tridiagonal matrices
    use time_discr_m, only: time_discr_homog_c !> Import homogeneous time discretization class
    use metodos_sist_lin_m, only: compute_inverse_tridiag_matrix !> Import subroutine for inverting tridiagonal matrices
    implicit none !> Enforce explicit variable declarations
    
    class(transport_1D_transient_c) :: this !> Transport object with spatial/temporal discretization, BCs, and solution arrays [-]
    !type(diag_matrix_c), intent(out), optional :: A_mat_lumped !> Optional lumped system matrix for mass-lumped formulations [-]
    
    integer(kind=4) :: i,j !> Loop indices for targets and neighbors [-]
    integer(kind=4) :: num_mix_rat !> Total number of mixing ratios in each target (domain + boundary + recharge) [-]
    integer(kind=4) :: num_inf_ext_mix_rat !> Number of inflow external mixing ratios outside the mesh (0 or 1) [-]
    integer(kind=4) :: num_out_ext_mix_rat !> Number of outflow boundary mixing ratios outside the mesh (0 or 1) [-]
    integer(kind=4) :: num_rech_mix_rat !> Number of recharge mixing ratios in each target (currently set to 0) [-]
    real(kind=8) :: theta,Delta_t !> Time integration parameter θ (0=explicit, 1=fully implicit, 0.5=Crank-Nicolson) [-] and time step size [T]
    real(kind=8), parameter :: tol_inv=1d-12 !> Tolerance for matrix inversion convergence [-]
    real(kind=8), allocatable :: A_mat_inv(:,:) !> Inverse of the system matrix A (not currently used, for potential alternative implementations) [-]
    real(kind=8) :: mix_sum !> Sum of mixing ratios for mass balance check [-]
    
    type(tridiag_matrix_c) :: E_mat,E_mat_prev,X_mat_T !> Tridiagonal matrices: E = transport evolution matrix, X_mat_T = transpose of X matrix [-]
    
!> Step 1: Allocate mixing ratio arrays for concentrations and reaction amounts
    call this%mixing_ratios_conc%allocate_array(this%spatial_discr%Num_targets) !> Allocate array for concentration mixing ratios (one column per target) [-]
    call this%mixing_ratios_R%allocate_array(this%spatial_discr%Num_targets) !> Allocate array for kinetic reaction amount mixing ratios (one column per target) [-]
    !call this%mixing_ratios_R_init%allocate_array(this%spatial_discr%Num_targets) !> Allocate array for initial kinetic reaction amount mixing ratios (stored for reference) [-]
    
!> Step 2: Determine number of boundary mixing ratios based on boundary condition types
    num_mix_rat=0 !> Initialize total number of mixing ratios in each target to zero [-]
    !if (this%spatial_discr%targets_flag==1) then
    !     if (this%BCs%labels(1).eq.1) then
    !         num_inf_ext_mix_rat=0
    !     else
    !         num_inf_ext_mix_rat=1
    !     end if
    !     if (this%BCs%labels(2).eq.1) then
    !         num_out_ext_mix_rat=0
    !     else
    !         num_out_ext_mix_rat=1
    !     end if
    !    !num_inf_ext_mix_rat=0 !> One inflow boundary mixing ratio for cell targets (BCs handled differently) [-]
    !    !num_out_ext_mix_rat=0 !> One outflow boundary mixing ratio for cell targets (BCs handled differently) [-]
    !else
        !if (this%Lagr_flag) then
        !    num_inf_ext_mix_rat=1 !> No inflow external mixing ratio for edge targets (BCs handled differently) [-]
        !else
            num_inf_ext_mix_rat=1 !> No inflow external mixing ratio for edge targets (BCs handled differently) [-]
            num_out_ext_mix_rat=1 !> No inflow external mixing ratio for edge targets (BCs handled differently) [-]
        !end if
        !num_inf_ext_mix_rat=0 !> No inflow external mixing ratio for edge targets (BCs handled differently) [-]
        !num_out_ext_mix_rat=0 !> No outflow external mixing ratio for edge targets (BCs handled differently) [-]
        !if (this%BCs%labels(1).eq.3) then
        !    num_inf_ext_mix_rat=1
        !else
        !    num_inf_ext_mix_rat=0
        !end if
        !if (this%BCs%labels(2).eq.3) then
        !    num_out_ext_mix_rat=1
        !else
        !    num_out_ext_mix_rat=0
        !end if
    !end if
        !if (this%BCs%labels(1).eq.1 .and. this%targets_flag.eq.0) then !> Check if inflow BC is Dirichlet (label=1) and targets are cell-centered (flag=0)
        !    num_inf_ext_mix_rat=1 !> Dirichlet inflow BC with cell-centered targets requires 1 boundary mixing ratio [-]
        !else if (this%BCs%labels(1).eq.1 .and. this%targets_flag.eq.1) then !> Commented out: inflow Dirichlet BC with edge-centered targets (would require different treatment)
        !    num_inf_ext_mix_rat=0 !> Commented out: edge-centered Dirichlet BC would need zero boundary mixing ratios
        !else if (this%BCs%labels(1)>1 .and. this%BCs%labels(1)<4) then !> Check if inflow BC is mass flux (label=2) or dispersive flux (label=3)
        !    num_inf_ext_mix_rat=1 !> Neumann or Robin inflow BC requires 1 boundary mixing ratio [-]
        !else !> Any other inflow BC type (label≥4 or unrecognized)
        !    num_inf_ext_mix_rat=0 !> No inflow boundary mixing ratio needed [-]
        !end if
        !if (this%BCs%labels(2).eq.1 .and. this%targets_flag.eq.0) then !> Check if outflow BC is Dirichlet (label=1) and targets are cell-centered (flag=0)
        !    num_out_ext_mix_rat=1 !> Dirichlet outflow BC with cell-centered targets requires 1 boundary mixing ratio [-]
        !else if (this%BCs%labels(2).eq.1 .and. this%targets_flag.eq.1) then !> Commented out: outflow Dirichlet BC with edge-centered targets (would require different treatment)
        !    num_out_ext_mix_rat=0 !> Commented out: edge-centered Dirichlet BC would need zero boundary mixing ratios
        !else if (this%BCs%labels(2)>1 .and. this%BCs%labels(2)<4) then !> Check if outflow BC is mass flux (label=2) or dispersive flux (label=3)
        !    num_out_ext_mix_rat=1 !> Neumann or Robin outflow BC requires 1 boundary mixing ratio [-]
        !else !> Any other outflow BC type
        !    num_out_ext_mix_rat=0 !> No outflow boundary mixing ratio needed [-]
        !end if
    !else !> Alternative: targets are edges (commented out)
    !    num_inf_ext_mix_rat=0 !> Commented out: no inflow boundary mixing ratio for edge-centered targets (BCs handled differently)
    !    num_out_ext_mix_rat=0 !> Commented out: no outflow boundary mixing ratio for edge-centered targets (BCs handled differently)
    !end if
    ! select type (this) !> Commented out: type-check to determine recharge term presence based on transport properties
    ! class is (transport_1D_transient_c) !> Commented out: verify object is 1D transient transport class
    !     if (this%tpt_props_heterog%cst_flux_flag .eqv. .false. .and. this%spat) then !> Commented out: check if flux is non-constant and spatial flag is true
    !         num_rech_mix_rat=1 !> Commented out: non-constant flux requires recharge term mixing ratio
    !     else !> Commented out: constant flux or other condition
    !         num_rech_mix_rat=0 !> Commented out: no recharge mixing ratio needed
    !     end if
    ! end select
    num_rech_mix_rat=0 !> Temporary workaround: set number of recharge mixing ratios to zero (no recharge) [-]
    num_mix_rat=num_mix_rat+num_inf_ext_mix_rat+num_out_ext_mix_rat+num_rech_mix_rat !> Compute total number of boundary and recharge mixing ratios [-]
    
!> Step 3: Extract time step size and verify time discretization is homogeneous
    select type (time_discr=>this%time_discr) !> Type-select on time discretization object to access specific attributes
    class is (time_discr_homog_c) !> Check if time discretization is homogeneous (uniform Δt)
        Delta_t=time_discr%Delta_t !> Extract time step size from homogeneous time discretization [T]
    class default !> If time discretization is not homogeneous
        error stop "This subroutine is only applied if time discretisation is homogeneous" !> Abort execution with error message
    end select
    !Delta_t=this%time_discr%get_Delta_t() !> Alternative: call getter method to extract Δt (commented out)
!> Step 4: Compute PDE matrices (transport, recharge, source, and mass matrices)
    call this%compute_trans_mat_PDE() !> Compute transport matrix T for advection-dispersion equation [-]
    call this%compute_rech_mat_PDE() !> Compute recharge matrix Y for source/sink terms [-]
    call this%compute_source_term_PDE() !> Compute source term vector f for external sources [-]
    call this%compute_F_mat_PDE() !> Compute mass matrix F for transient accumulation term [-]
    !> We impose the boundary conditions
    if (this%BCs%labels(1).eq.1 .and. this%BCs%labels(2).eq.1) then !> Check if both inflow and outflow BCs are Dirichlet (labels both equal 1)
        call Dirichlet_BCs_PDE(this) !> Apply Dirichlet boundary conditions to PDE matrices (modify T, F, f to enforce prescribed concentrations)
        !this%mixing_ratios_conc%cols(1)%dim=1+num_rech_mix_rat+num_inf_ext_mix_rat !> Set dimension of mixing ratio vector for first target: current + recharge + inflow BC [-]
        !this%mixing_ratios_conc%cols(this%mixing_ratios_conc%num_cols)%dim=1+num_rech_mix_rat+num_out_ext_mix_rat !> Set dimension for last target: current + recharge + outflow BC [-]
        !num_Dir=2 !> Commented out: counter for number of Dirichlet BCs (would be 2 for this case)
    else if (this%BCs%labels(1).eq.1 .and. this%BCs%labels(2).eq.2) then !> Check if inflow is Dirichlet (1) and outflow is Neumann (2)
        call Dirichlet_Neumann_BCs_PDE(this) !> Apply Dirichlet-Neumann boundary conditions to PDE matrices
        !this%mixing_ratios_conc%cols(1)%dim=1+num_rech_mix_rat+num_inf_ext_mix_rat !> Set dimension for first target (Dirichlet inflow) [-]
        !this%mixing_ratios_conc%cols(this%mixing_ratios_conc%num_cols)%dim=2+num_rech_mix_rat+num_out_ext_mix_rat !> Set dimension for last target (Neumann outflow needs current + neighbor) [-]
    else if (this%BCs%labels(1).eq.2 .and. this%BCs%labels(2).eq.2) then !> Check if both inflow and outflow are Neumann (2)
        call Neumann_homog_BCs(this) !> Apply Neumann-Neumann boundary conditions to PDE matrices (homogeneous)
        !this%mixing_ratios_conc%cols(1)%dim=2+num_rech_mix_rat+num_inf_ext_mix_rat !> Set dimension for first target (Neumann needs current + neighbor) [-]
        !this%mixing_ratios_conc%cols(this%mixing_ratios_conc%num_cols)%dim=2+num_rech_mix_rat+num_out_ext_mix_rat !> Set dimension for last target (Neumann needs current + neighbor) [-]
    else if (this%BCs%labels(1).eq.3 .and. this%BCs%labels(2).eq.2) then !> Check if inflow is Robin (3) and outflow is Neumann (2)
        call Robin_Neumann_homog_BCs(this) !> Apply Robin-Neumann boundary conditions to PDE matrices (homogeneous)
        !this%mixing_ratios_conc%cols(1)%dim=2+num_rech_mix_rat+num_inf_ext_mix_rat !> Set dimension for first target (Robin needs current + neighbor) [-]
        !this%mixing_ratios_conc%cols(this%mixing_ratios_conc%num_cols)%dim=2+num_rech_mix_rat+num_out_ext_mix_rat !> Set dimension for last target (Neumann needs current + neighbor) [-]
    else if (this%BCs%labels(1).eq.3 .and. this%BCs%labels(2).eq.1) then !> Check if inflow is Robin (3) and outflow is Dirichlet (1)
        call Robin_Dirichlet_BCs(this) !> Apply Robin-Dirichlet boundary conditions to PDE matrices
        !this%mixing_ratios_conc%cols(1)%dim=2+num_rech_mix_rat+num_inf_ext_mix_rat !> Set dimension for first target (Robin needs current + neighbor) [-]
        !this%mixing_ratios_conc%cols(this%mixing_ratios_conc%num_cols)%dim=1+num_rech_mix_rat+num_out_ext_mix_rat !> Set dimension for last target (Dirichlet only needs current) [-]
    else !> Any other combination of boundary conditions
        error stop "Boundary conditions not implemented yet" !> Abort execution if BC combination is not supported
    end if
    if (this%time_discr%int_method.eq.1) then !> Check if time integration method is Euler explicit (method = 1)
        !this%time_discr%theta_t=0d0 !> Set time integration parameter to 0 for fully explicit method [-]
        this%mixing_ratios_conc%cols(1)%dim=2+num_rech_mix_rat+num_inf_ext_mix_rat !> Set dimension for first target (Robin needs current + neighbor) [-]
        this%mixing_ratios_conc%cols(this%mixing_ratios_conc%num_cols)%dim=2+num_rech_mix_rat+num_out_ext_mix_rat !> Set dimension for last target (Neumann needs current + neighbor) [-]
        this%mixing_ratios_R%cols(1)%dim=1 !> Set dimension for kinetic reaction mixing ratios: only current target (explicit) [-]
        !this%mixing_ratios_R_init%cols(1)%dim=1 !> Set dimension for initial kinetic reaction mixing ratios: only current target [-]
        this%mixing_ratios_R%cols(this%mixing_ratios_conc%num_cols)%dim=1 !> Set dimension for kinetic reaction mixing ratios at last target [-]
        !this%mixing_ratios_R_init%cols(this%mixing_ratios_conc%num_cols)%dim=1 !> Set dimension for initial kinetic reaction mixing ratios at last target [-]
        do i=2,this%mixing_ratios_conc%num_cols-1 !> Loop over interior targets (exclude first and last boundaries)
            this%mixing_ratios_conc%cols(i)%dim=3+num_rech_mix_rat !> Set dimension for interior targets: current + upstream + downstream + recharge [-]
            this%mixing_ratios_R%cols(i)%dim=1 !> Set dimension for kinetic reaction mixing ratios: only current target (explicit) [-]
            !this%mixing_ratios_R_init%cols(i)%dim=1 !> Set dimension for initial kinetic reaction mixing ratios: only current target [-]
        end do
        call this%mixing_ratios_conc%allocate_columns() !> Allocate memory for all mixing ratio columns (concentration) [-]
        call this%mixing_ratios_R%allocate_columns() !> Allocate memory for all kinetic reaction mixing ratio columns [-]
        !call this%mixing_ratios_R_init%allocate_columns() !> Allocate memory for initial kinetic reaction mixing ratio columns [-]
        !call this%allocate_mix_conc_indices() !> Commented out: alternative location for allocating mixing water indices
    else !> Time integration method is implicit (Euler implicit or Crank-Nicolson)
        num_mix_rat=num_mix_rat+this%spatial_discr%Num_targets !> Add number of mesh targets to total mixing ratios (implicit couples all targets) [-]
        !call this%mix_conc_indices%allocate_array(this%mixing_ratios_conc%num_cols) !> Commented out: alternative allocation for mixing water indices
        do i=1,this%mixing_ratios_conc%num_cols !> Loop over all targets (including boundaries)
            call this%mixing_ratios_conc%cols(i)%set_dim(num_mix_rat) !> Set dimension for target i: includes all domain targets + boundaries + recharge (implicit coupling) [-]
            call this%mixing_ratios_conc%cols(i)%allocate_vector() !> Allocate memory for mixing ratio vector for target i [-]
            call this%mixing_ratios_R%cols(i)%set_dim(this%mixing_ratios_R%num_cols) !> Set dimension for kinetic reaction mixing ratios: equal to number of targets [-]
            !call this%mixing_ratios_R_init%cols(i)%set_dim(this%mixing_ratios_R%num_cols) !> Set dimension for initial kinetic reaction mixing ratios: equal to number of targets [-]
            call this%mixing_ratios_R%cols(i)%allocate_vector() !> Allocate memory for kinetic reaction mixing ratio vector for target i [-]
            !call this%mixing_ratios_R_init%cols(i)%allocate_vector() !> Allocate memory for initial kinetic reaction mixing ratio vector for target i [-]
            !call this%mix_conc_indices%cols(i)%allocate_vector(this%mixing_ratios_conc%cols(i)%dim+2) !> Commented out: allocate indices with extra space for upstream/downstream counts
        end do
        !if (this%time_discr%int_method.eq.2) then !> Check if time integration method is Euler fully implicit (method = 2)
        !    this%time_discr%theta_t=1d0 !> Set time integration parameter to 1 for fully implicit method [-]
        !else if (this%time_discr%int_method.eq.3) then !> Check if time integration method is Crank-Nicolson (method = 3)
        !    this%time_discr%theta_t=5d-1 !> Set time integration parameter to 0.5 for Crank-Nicolson method [-]
        !else !> Any other time integration method
        !    error stop "Time discretisation not implemented yet" !> Abort execution if integration method is not supported
        !end if
    end if
    call this%allocate_mix_conc_indices() !> Allocate arrays for storing indices of waters that contribute to mixing in each target
!> Step 6: Compute linear system matrices (E, X, A, Y, Z) and source vector (f)
    call this%compute_E_mat_1D(Delta_t,E_mat,E_mat_prev) !> Compute E matrix: evolution matrix for explicit part E = F + θ·Δt·T (tridiagonal) [-]
    call this%compute_X_mat_1D(this%time_discr%theta_t,E_mat_prev) !> Compute X matrix: coefficient matrix for old time level X = F - (1-θ)·Δt·T (tridiagonal) [-]
    call this%compute_A_mat_1D(this%time_discr%theta_t,E_mat) !> Compute A matrix: system matrix for new time level A = F + θ·Δt·T (tridiagonal) [-]
    call this%compute_Y_mat() !> Compute Y matrix: recharge/source term coefficient matrix (diagonal) [-]
    call this%compute_Z_mat() !> Compute Z matrix: boundary condition coefficient matrix (scalar or vector) [-]
    call this%compute_f_vec(this%time_discr%get_Delta_t()) !> Compute f vector: external source term vector [-]
    
!> Step 7: Compute mixing ratios from system matrices
    !if (present(A_mat_lumped)) then !> Check if optional lumped matrix argument is provided
        !call this%compute_lumped_A_mat(A_mat_lumped) !> Compute lumped (mass-lumped) A matrix for mass balance (temporary implementation) [-]
    if (this%time_discr%theta_t>0d0) then !> Check if implicit time integration is used (θ > 0)
!> Step 7a: Allocate arrays for implicit method and invert system matrix A
        call this%allocate_A_mat_inv() !> Allocate memory for inverse of A matrix (A^{-1}) [-]
        call this%allocate_mixing_ratios_mat_conc_mesh() !> Allocate domain mixing ratio matrix for concentrations [-]
        call this%allocate_mixing_ratios_mat_conc_bd() !> Allocate boundary mixing ratio matrix for concentrations [-]
        !allocate(this%mixing_ratios_mat_conc(this%mixing_ratios_conc%num_cols,this%mixing_ratios_conc%num_cols)) !> Commented out: alternative full mixing ratio matrix allocation
        call compute_inverse_tridiag_matrix(this%A_mat,tol_inv,this%A_mat_inv) !> Invert tridiagonal system matrix A using tolerance tol_inv to get A^{-1} [-]
        !print *, this%A_mat_inv(1,:) !> Print inverse of A matrix for debugging (remove in production)
        if (this%time_discr%theta_t.eq.1d0) then !> Check if fully implicit Euler method is used (θ = 1)
!> Step 7b: Fully implicit Euler - mixing ratios from A^{-1} only
            this%mixing_ratios_mat_conc_mesh=transpose(this%A_mat_inv) !> Domain mixing ratios: transpose of A^{-1} (each column is mixing ratios for one target) [-]
            !if (this%spatial_discr%targets_flag==0) then !> Check if targets are cell-centered (flag = 0)
                !> Inflow boundary mixing ratios for cell-centered targets
                this%mixing_ratios_mat_conc_bd(:,1)=this%Z_mat(1)*this%mixing_ratios_mat_conc_mesh(1,:) !> Inflow boundary mixing ratios: Z(1) · (A^{-1})_1 (contribution of inflow BC to all targets) [-]
                !> Outflow boundary mixing ratios for cell-centered targets
                this%mixing_ratios_mat_conc_bd(:,2)=this%Z_mat(2)*this%mixing_ratios_mat_conc_mesh(this%spatial_discr%Num_targets,:) !> Outflow boundary mixing ratios: Z(2) · (A^{-1})_N (contribution of outflow BC to all targets) [-]
            !else if (this%spatial_discr%targets_flag==1) then !> Check if targets are edge-centered (flag = 1)
            !    !> Inflow boundary mixing ratios for edge-centered targets
            !    if (this%Lagr_flag) then !> Check if Lagrangian tracking is enabled (would affect inflow BC treatment)
            !        this%mixing_ratios_mat_conc_bd(:,1)=0d0 !> Set inflow boundary mixing ratios to zero for edge targets with Lagrangian tracking (BCs handled differently) [-] else this%mixing_ratios_mat_conc_bd(:,1)=this%mixing_ratios_mat_conc_mesh(1,:) !> Inflow boundary mixing ratios: Z(1) · (A^{-1})_1 (contribution of inflow BC to all targets) [-] end if !> Outflow boundary mixing ratios for edge-centered targets
            !    else
            !        this%mixing_ratios_mat_conc_bd(:,1)=this%mixing_ratios_mat_conc_mesh(1,:) !> Inflow boundary mixing ratios: Z(1) · (A^{-1})_1 (contribution of inflow BC to all targets) [-]
            !    end if  
            !    !> Outflow boundary mixing ratios for edge-centered targets
            !    this%mixing_ratios_mat_conc_bd(:,2)=this%mixing_ratios_mat_conc_mesh(this%spatial_discr%Num_targets,:) !> Outflow boundary mixing ratios: Z(2) · (A^{-1})_(N+1) (contribution of outflow BC to all targets) [-]
            !end if
        else !> Crank-Nicolson or other partially implicit method (0 < θ < 1)
!> Step 7c: Partially implicit - mixing ratios from X^T · A^{-1}
            call this%X_mat%compute_transpose_tridiag_matrix(X_mat_T) !> Compute transpose of X matrix (X^T) to get coefficient matrix for old time level [-]
            this%mixing_ratios_mat_conc_mesh=X_mat_T%prod_tridiag_mat_mat(transpose(this%A_mat_inv)) !> Domain mixing ratios: X^T · (A^{-1})^T (combined old/new time level contributions) [-]
            !print *, this%mixing_ratios_mat_conc_mesh(:,1) !> Print domain mixing ratios for debugging (remove in production)
        !end if
            !if (this%spatial_discr%targets_flag==0 .or. this%Lagr_flag) then !> Check if targets are cell-centered (flag = 0)
                !> Inflow boundary mixing ratios for cell-centered targets
                this%mixing_ratios_mat_conc_bd(:,1)=&
                    (this%time_discr%theta_t*this%Z_mat(1)+(1d0-this%time_discr%theta_t)*this%Z_mat_prev(1))*this%A_mat_inv(:,1) !> Inflow boundary mixing ratios: Z(1) · (A^{-1})_1 (contribution of inflow BC to all targets) [-]
                !> Outflow boundary mixing ratios for cell-centered targets
                this%mixing_ratios_mat_conc_bd(:,2)=&
                    (this%time_discr%theta_t*this%Z_mat(2)+(1d0-this%time_discr%theta_t)*this%Z_mat_prev(2))*this%A_mat_inv(:,this%spatial_discr%Num_targets) !> Outflow boundary mixing ratios: Z(2) · (A^{-1})_N (contribution of outflow BC to all targets) [-]
            !else if (this%spatial_discr%targets_flag==1) then !> Check if targets are edge-centered (flag = 1)
            !else
            !    !if (this%Lagr_flag) then !> Check if Lagrangian tracking is enabled (would affect inflow BC treatment)
            !    !    this%mixing_ratios_mat_conc_bd(:,1)=0d0 !> Set inflow boundary mixing ratios to zero for edge targets with Lagrangian tracking (BCs handled differently) [-] else this%mixing_ratios_mat_conc_bd(:,1)=this%mixing_ratios_mat_conc_mesh(1,:) !> Inflow boundary mixing ratios: Z(1) · (A^{-1})_1 (contribution of inflow BC to all targets) [-] end if !> Outflow boundary mixing ratios for edge-centered targets
            !    !else
            !        this%mixing_ratios_mat_conc_bd(:,1)=this%mixing_ratios_mat_conc_mesh(1,:) !> Inflow boundary mixing ratios: Z(1) · (A^{-1})_1 (contribution of inflow BC to all targets) [-]
            !    !end if
            !    !> Outflow boundary mixing ratios for edge-centered targets
            !    this%mixing_ratios_mat_conc_bd(:,2)=this%mixing_ratios_mat_conc_mesh(this%spatial_discr%Num_targets,:) !> Outflow boundary mixing ratios: Z(2) · (A^{-1})_(N+1) (contribution of outflow BC to all targets) [-]
            !end if
        end if
        !this%mixing_ratios_mat_Rk=transpose(A_mat_inv) !> Commented out: alternative storage of reaction mixing ratio matrix (transpose of A^{-1})
        !if (num_mix_rat>2*this%spatial_discr%Num_targets .or. num_mix_rat<2*this%spatial_discr%num_targets) then !> Commented out: condition for boundary mixing with no recharge
        !    this%mixing_ratios_conc%cols(1)%col_1=this%bd_mat(1)*this%A_mat_inv(:,1) !> Commented out: inflow boundary mixing ratios (Z(1) · first column of A^{-1})
        !    this%mixing_ratios_conc%cols(this%spatial_discr%Num_targets+2)%col_1=this%bd_mat(2)*this%A_mat_inv(:,this%spatial_discr%num_targets) !> Commented out: outflow boundary mixing ratios (Z(2) · last column of A^{-1})
        !else if (num_mix_rat>this%spatial_discr%Num_targets) then !> Commented out: condition for recharge with no boundary mixing
        !    this%mixing_ratios_mat_Rk=this%A_mat_inv !> Commented out: reaction mixing ratio matrix equals A^{-1} for recharge case
        !end if
        !if (num_inf_ext_mix_rat==1) then !> Commented out: alternative implementation for inflow boundary mixing ratios
        !    this%mixing_ratios_conc%cols(1)%col_1=this%bd_mat(1)*this%A_mat_inv(:,1) !> Commented out: compute inflow boundary mixing ratios
        !    this%mix_conc_indices%cols(1)%col_1(1)=1 !> Commented out: set current water index to 1
        !    this%mix_conc_indices%cols(1)%col_1(this%mix_conc_indices%cols(1)%dim-1)=0 !> Commented out: set upstream water count to 0 (no upstream for inflow BC)
        !    this%mix_conc_indices%cols(1)%col_1(this%mix_conc_indices%cols(1)%dim)=this%spatial_discr%num_targets+num_out_ext_mix_rat !> Commented out: set downstream water count
        !end if
        !if (num_out_ext_mix_rat==1) then !> Commented out: alternative implementation for outflow boundary mixing ratios
        !    this%mixing_ratios_conc%cols(this%spatial_discr%num_targets+num_inf_ext_mix_rat+1)%col_1=this%bd_mat(2)*this%A_mat_inv(:,this%spatial_discr%Num_targets) !> Commented out: compute outflow boundary mixing ratios
        !    this%mix_conc_indices%cols(this%spatial_discr%num_targets+num_inf_ext_mix_rat+1)%col_1(1)=this%spatial_discr%num_targets+num_inf_ext_mix_rat+1 !> Commented out: set current water index
        !    this%mix_conc_indices%cols(this%spatial_discr%num_targets+num_inf_ext_mix_rat+1)%col_1(this%mix_conc_indices%cols(this%spatial_discr%num_targets+num_inf_ext_mix_rat+1)%dim-1)=this%spatial_discr%num_targets+num_inf_ext_mix_rat !> Commented out: set upstream water count
        !    this%mix_conc_indices%cols(this%spatial_discr%num_targets+num_inf_ext_mix_rat+1)%col_1(this%mix_conc_indices%cols(this%spatial_discr%num_targets+num_inf_ext_mix_rat+1)%dim)=0 !> Commented out: set downstream water count to 0 (no downstream for outflow BC)
        !end if
!> Step 7d: Populate mixing ratios for each target (implicit methods)
!         do i=1,this%spatial_discr%num_targets !> Loop over all spatial targets to assign mixing ratios [-]
! !> Assign boundary mixing ratios for target i
!             this%mixing_ratios_conc%cols(i)%col_1(1+num_inf_ext_mix_rat)=this%Z_mat(1)*this%A_mat_inv(i,1) !> Concentration mixing ratio of inflow boundary: contribution of inflow BC to target i [-]
!             this%mixing_ratios_conc%cols(i)%col_1(this%spatial_discr%num_targets+num_inf_ext_mix_rat+num_out_ext_mix_rat)=& !> Index of outflow boundary mixing ratio in target i's vector
!                 this%Z_mat(2)*this%A_mat_inv(i,this%spatial_discr%num_targets) !> Concentration mixing ratio of outflow boundary: contribution of outflow BC to target i [-]
! !> Assign current target mixing ratios
!             this%mixing_ratios_conc%cols(i)%col_1(1)=this%mixing_ratios_mat_conc_mesh(i,i) !> Concentration mixing ratio of current target water (diagonal element of domain mixing ratio matrix) [-]
!             this%mixing_ratios_R%cols(i)%col_1(1)=this%A_mat_inv(i,i) !> Reaction amount mixing ratio of current target (diagonal element of A^{-1}) [-]
!             this%mixing_ratios_R_init%cols(i)%col_1(1)=this%mixing_ratios_R%cols(i)%col_1(1) !> Store initial reaction amount mixing ratio of current target (for reference) [-]
! !> Assign water indices for target i
!             this%mix_conc_indices%cols(i)%col_1(1+num_inf_ext_mix_rat)=num_inf_ext_mix_rat !> Inflow water index (typically 0 or 1) [-]
!             this%mix_conc_indices%cols(i)%col_1(this%spatial_discr%num_targets+num_inf_ext_mix_rat+num_out_ext_mix_rat)=& !> Index position of outflow water index
!                 this%spatial_discr%num_targets+num_inf_ext_mix_rat+num_out_ext_mix_rat !> Outflow water index (last index after all domain and boundary waters) [-]
!             this%mix_conc_indices%cols(i)%col_1(1)=i !> Current water index: target index i [-]
!             this%mix_react_indices%cols(i)%col_1(1)=this%mix_conc_indices%cols(i)%col_1(1) !> Copy current water index to domain indices array [-]
!             this%mix_conc_indices%cols(i)%col_1(this%mix_conc_indices%cols(i)%dim-1)=num_inf_ext_mix_rat+i-1 !> Number of upstream waters (all targets before target i) [-]
!             this%mix_react_indices%cols(i)%col_1(this%mix_react_indices%cols(i)%dim-1)=i-1 !> Number of upstream waters in domain indices (targets before i) [-]
!             this%mix_conc_indices%cols(i)%col_1(this%mix_conc_indices%cols(i)%dim)=this%mixing_ratios_conc%cols(i)%dim-& !> Total dimension of mixing ratio vector minus...
!                 i-num_inf_ext_mix_rat !> ...current index and inflow indices = number of downstream waters [-]
!             this%mix_react_indices%cols(i)%col_1(this%mix_react_indices%cols(i)%dim)=& !> Last position in domain indices array
!                 this%spatial_discr%num_targets-i !> Number of downstream waters in domain (targets after i) [-]
!         end do !> End first loop over all targets for setting up water indices and diagonal mixing ratios
        !this%mixing_ratios_mat_Rk=transpose(A_mat_inv) !> mixing ratios of kinetic reaction amounts
        !if (num_mix_rat>2*this%spatial_discr%Num_targets .or. num_mix_rat<2*this%spatial_discr%num_targets) then !> boundary and no recharge
        !    this%mixing_ratios_conc%cols(1)%col_1=this%bd_mat(1)*this%A_mat_inv(:,1) !> inflow boundary mixing ratios
        !    this%mixing_ratios_conc%cols(this%spatial_discr%Num_targets+2)%col_1=this%bd_mat(2)*this%A_mat_inv(:,this%spatial_discr%num_targets) !> outflow boundary mixing ratios
        !else if (num_mix_rat>this%spatial_discr%Num_targets) then !> recharge with no boundary
        !    this%mixing_ratios_mat_Rk=this%A_mat_inv
        !end if
        !if (num_inf_ext_mix_rat==1) then !> inflow boundary mixing ratios
        !    this%mixing_ratios_conc%cols(1)%col_1=this%bd_mat(1)*this%A_mat_inv(:,1)
        !    this%mix_conc_indices%cols(1)%col_1(1)=1 !> current water
        !    this%mix_conc_indices%cols(1)%col_1(this%mix_conc_indices%cols(1)%dim-1)=0 !> upstream waters
        !    this%mix_conc_indices%cols(1)%col_1(this%mix_conc_indices%cols(1)%dim)=this%spatial_discr%num_targets+num_out_ext_mix_rat !> downstream waters
        !end if
        !if (num_out_ext_mix_rat==1) then !> outflow boundary mixing ratios
        !    this%mixing_ratios_conc%cols(this%spatial_discr%num_targets+num_inf_ext_mix_rat+1)%col_1=this%bd_mat(2)*this%A_mat_inv(:,this%spatial_discr%Num_targets)
        !    this%mix_conc_indices%cols(this%spatial_discr%num_targets+num_inf_ext_mix_rat+1)%col_1(1)=this%spatial_discr%num_targets+num_inf_ext_mix_rat+1 !> current water
        !    this%mix_conc_indices%cols(this%spatial_discr%num_targets+num_inf_ext_mix_rat+1)%col_1(this%mix_conc_indices%cols(this%spatial_discr%num_targets+num_inf_ext_mix_rat+1)%dim-1)=this%spatial_discr%num_targets+num_inf_ext_mix_rat !> upstream waters
        !    this%mix_conc_indices%cols(this%spatial_discr%num_targets+num_inf_ext_mix_rat+1)%col_1(this%mix_conc_indices%cols(this%spatial_discr%num_targets+num_inf_ext_mix_rat+1)%dim)=0 !> downstream waters
        !end if
    !> Step 7e: Populate mixing ratios for each target (implicit methods)
        do i=1,this%spatial_discr%num_targets
            this%mixing_ratios_conc%cols(i)%col_1(1+num_inf_ext_mix_rat)=this%mixing_ratios_mat_conc_bd(i,1) !> concentration mixing ratio of inflow boundary
            this%mixing_ratios_conc%cols(i)%col_1(this%spatial_discr%num_targets+num_inf_ext_mix_rat+num_out_ext_mix_rat)=&
                this%mixing_ratios_mat_conc_bd(i,2) !> concentration mixing ratio of outflow boundary
            this%mixing_ratios_conc%cols(i)%col_1(1)=this%mixing_ratios_mat_conc_mesh(i,i) !> concentration mixing ratio of current water
            this%mixing_ratios_R%cols(i)%col_1(1)=this%A_mat_inv(i,i) !> reaction amount mixing ratio of current water
            !this%mixing_ratios_R_init%cols(i)%col_1(1)=this%mixing_ratios_R%cols(i)%col_1(1) !> initial reaction amount mixing ratio of current water
            this%mix_conc_indices%cols(i)%col_1(1+num_inf_ext_mix_rat)=num_inf_ext_mix_rat !> inflow water index
            this%mix_conc_indices%cols(i)%col_1(this%spatial_discr%num_targets+num_inf_ext_mix_rat+num_out_ext_mix_rat)=&
                this%spatial_discr%num_targets+num_inf_ext_mix_rat+num_out_ext_mix_rat !> outflow water index
            this%mix_conc_indices%cols(i)%col_1(1)=i+num_inf_ext_mix_rat !> current water index
            !this%mix_react_indices%cols(i)%col_1(1)=this%mix_conc_indices%cols(i)%col_1(1) !> current water index
            this%mix_conc_indices%cols(i)%col_1(this%mix_conc_indices%cols(i)%dim-1)=num_inf_ext_mix_rat+i-1 !> upstream waters
            !this%mix_react_indices%cols(i)%col_1(this%mix_react_indices%cols(i)%dim-1)=i-1-&
            !    this%spatial_discr%targets_flag !> upstream waters in domain
            this%mix_conc_indices%cols(i)%col_1(this%mix_conc_indices%cols(i)%dim)=this%mixing_ratios_conc%cols(i)%dim-&
                i-num_inf_ext_mix_rat !> number of downstream waters
            !this%mix_react_indices%cols(i)%col_1(this%mix_react_indices%cols(i)%dim)=&
            !    this%spatial_discr%num_targets-i-this%spatial_discr%targets_flag !> Number of downstream waters in domain (targets after i) [-]
!> Assign mixing ratios for all upstream targets (j < i)
            do j=1,i-1 !> Loop over all targets upstream of target i
                this%mixing_ratios_conc%cols(i)%col_1(1+num_inf_ext_mix_rat+j)=this%mixing_ratios_mat_conc_mesh(j,i) !> Concentration mixing ratio of upstream target j on target i (element (j,i) of domain mixing matrix) [-]
                this%mixing_ratios_R%cols(i)%col_1(1+j)=this%A_mat_inv(i,j) !> Reaction amount mixing ratio of upstream target j on target i (element (i,j) of A^{-1}) [-]
                !this%mixing_ratios_R_init%cols(i)%col_1(1+j)=this%mixing_ratios_R%cols(i)%col_1(1+j) !> Store initial reaction amount mixing ratio of upstream target j (for reference) [-]
                this%mix_conc_indices%cols(i)%col_1(1+num_inf_ext_mix_rat+j)=num_inf_ext_mix_rat+j !> Water index of upstream target j (inflow indices + j) [-]
                ! this%mix_react_indices%cols(i)%col_1(1+j-this%spatial_discr%targets_flag)=&
                !     this%mix_conc_indices%cols(i)%col_1(1+num_inf_ext_mix_rat+j-this%spatial_discr%targets_flag) !> Copy upstream water index to domain indices array [-]
            end do
!> Assign mixing ratios for all downstream targets (j > i)
            do j=i+1,this%spatial_discr%num_targets !> Loop over all targets downstream of target i
                this%mixing_ratios_conc%cols(i)%col_1(num_inf_ext_mix_rat+j)=this%mixing_ratios_mat_conc_mesh(j,i) !> Concentration mixing ratio of downstream target j on target i (element (j,i) of domain mixing matrix) [-]
                this%mixing_ratios_R%cols(i)%col_1(j)=this%A_mat_inv(i,j) !> Reaction amount mixing ratio of downstream target j on target i (element (i,j) of A^{-1}) [-]
                !this%mixing_ratios_R_init%cols(i)%col_1(j)=this%mixing_ratios_R%cols(i)%col_1(j) !> Store initial reaction amount mixing ratio of downstream target j (for reference) [-]
                this%mix_conc_indices%cols(i)%col_1(num_inf_ext_mix_rat+j)=num_inf_ext_mix_rat+j !> Water index of downstream target j (inflow indices + j) [-]
                ! this%mix_react_indices%cols(i)%col_1(j-this%spatial_discr%targets_flag)=&
                !     this%mix_conc_indices%cols(i)%col_1(num_inf_ext_mix_rat+j-this%spatial_discr%targets_flag) !> Copy downstream water index to domain indices array [-]
            end do
!> Assign recharge mixing ratios (if any)
            do j=1,num_rech_mix_rat !> Loop over all recharge waters (currently num_rech_mix_rat=0, so this loop is skipped)
                this%mixing_ratios_conc%cols(i)%col_1(this%spatial_discr%num_targets+num_inf_ext_mix_rat+num_out_ext_mix_rat+j)=& !> Index of recharge water j in target i's mixing ratio vector
                    this%Y_mat%diag(j)*this%A_mat_inv(i,j) !> Mixing ratio of recharge water j: Y(j) · (A^{-1})_{i,j} (contribution of recharge to target i) [-]
                this%mix_conc_indices%cols(i)%col_1(this%spatial_discr%num_targets+num_inf_ext_mix_rat+num_out_ext_mix_rat+j)=& !> Index position of recharge water j
                    this%spatial_discr%num_targets+num_inf_ext_mix_rat+num_out_ext_mix_rat+j !> Water index of recharge water j (after all domain and boundary waters) [-]
            end do
        end do !> End loop over all targets
    else !> Explicit Euler method (θ = 0)
!> Step 8: Assign mixing ratios directly from X matrix for explicit methods
        !> Mixing ratios of reaction amounts is the identity matrix (targets mix only with themselves for reactions)
        this%mix_conc_indices%cols(1)%col_1(1)=1+num_inf_ext_mix_rat !> First element: current water index for first target [-]
        if (this%BCs%labels(1).eq.1) then !> Check if inflow boundary condition is Dirichlet
!> First target with Dirichlet inflow BC
            this%mixing_ratios_conc%cols(1)%col_1(1+num_inf_ext_mix_rat+num_rech_mix_rat)=this%Y_mat%diag(1) !> Mixing ratio from recharge (Y matrix diagonal) [-]
            this%mixing_ratios_conc%cols(1)%col_1(1+num_inf_ext_mix_rat)=this%Z_mat(1) !> Mixing ratio from inflow boundary (Z matrix element 1) [-]
            this%mixing_ratios_conc%cols(1)%col_1(1)=this%X_mat%diag(1) !> Mixing ratio of current target (X matrix diagonal element 1) [-]
            this%mix_conc_indices%cols(1)%col_1(2)=num_inf_ext_mix_rat !> Second element: inflow boundary water index [-]
            this%mix_conc_indices%cols(1)%col_1(3)=this%spatial_discr%num_targets+num_inf_ext_mix_rat+num_out_ext_mix_rat+& !> Third element: total number of waters
                num_rech_mix_rat !> (domain + boundary + recharge) [-]
        else !> Inflow BC is Neumann or Robin (non-Dirichlet)
!> First target with non-Dirichlet inflow BC
            this%mixing_ratios_conc%cols(1)%col_1=[this%X_mat%diag(1),this%X_mat%super(1)] !> Mixing ratios from current target and downstream neighbor (X diagonal and superdiagonal) [-]
            this%mix_conc_indices%cols(1)%col_1(2)=2 !> Second element: downstream water index (target 2) [-]
            this%mix_conc_indices%cols(1)%col_1(3)=0 !> Third element: number of upstream waters (none for first target) [-]
            this%mix_conc_indices%cols(1)%col_1(4)=1 !> Fourth element: number of downstream waters (one neighbor) [-]
        end if
        this%mixing_ratios_R%cols(1)%col_1=1d0 !> Mixing ratio of kinetic reaction rate in first target: identity (1.0) [-]
        !this%mixing_ratios_R_init%cols(1)%col_1=this%mixing_ratios_R%cols(1)%col_1 !> Store initial mixing ratio of kinetic reaction rate in first target [-]
!> Interior targets (explicit Euler)
        do i=2,this%mixing_ratios_conc%num_cols-1 !> Loop over interior targets (exclude first and last)
            this%mixing_ratios_conc%cols(i)%col_1=[this%X_mat%diag(i),this%X_mat%sub(i-1),this%X_mat%super(i)] !> Mixing ratios: diagonal (current), subdiagonal (upstream), superdiagonal (downstream) elements of X [-]
            this%mixing_ratios_R%cols(i)%col_1=1d0 !> Mixing ratio of kinetic reaction rate in target i: identity (1.0) [-]
            !this%mixing_ratios_R_init%cols(i)%col_1=this%mixing_ratios_R%cols(i)%col_1 !> Store initial mixing ratio of kinetic reaction rate in target i [-]
            this%mix_conc_indices%cols(i)%col_1(1)=i !> First element: current water index (target i) [-]
            this%mix_conc_indices%cols(i)%col_1(2)=i-1 !> Second element: upstream water index (target i-1) [-]
            this%mix_conc_indices%cols(i)%col_1(3)=i+1 !> Third element: downstream water index (target i+1) [-]
            this%mix_conc_indices%cols(i)%col_1(4)=1 !> Fourth element: number of upstream waters (one neighbor) [-]
            this%mix_conc_indices%cols(i)%col_1(5)=1 !> Fifth element: number of downstream waters (one neighbor) [-]
        end do
!> Last target (explicit Euler)
        this%mixing_ratios_R%cols(this%mixing_ratios_conc%num_cols)%col_1=1d0 !> Mixing ratio of kinetic reaction rate in last target: identity (1.0) [-]
        !this%mixing_ratios_R_init%cols(this%mixing_ratios_conc%num_cols)%col_1=& !> Store initial mixing ratio for last target
            !this%mixing_ratios_R%cols(this%mixing_ratios_conc%num_cols)%col_1 !> Copy from current mixing ratio [-]
        this%mix_conc_indices%cols(this%mixing_ratios_conc%num_cols)%col_1(1)=1 !> First element: placeholder or flag value [-]
        if (this%BCs%labels(2).eq.1) then !> Check if outflow boundary condition is Dirichlet
!> Last target with Dirichlet outflow BC
            this%mixing_ratios_conc%cols(this%mixing_ratios_conc%num_cols)%col_1(this%mixing_ratios_conc%cols(& !> Access last element of mixing ratio vector
                this%mixing_ratios_conc%num_cols)%dim)=& !> Dimension of last target's mixing ratio vector
                this%X_mat%diag(this%mixing_ratios_conc%num_cols) !> Mixing ratio of current target (X diagonal element for last target) [-]
            this%mix_conc_indices%cols(this%mixing_ratios_conc%num_cols)%col_1(2)=0 !> Second element: number of upstream waters (zero for Dirichlet BC) [-]
            this%mix_conc_indices%cols(1)%col_1(3)=0 !> Third element: number of downstream waters (zero for last target) [-]
        else !> Outflow BC is Neumann or Robin (non-Dirichlet)
!> Last target with non-Dirichlet outflow BC
            this%mixing_ratios_conc%cols(this%mixing_ratios_conc%num_cols)%col_1=[& !> Array constructor for mixing ratios
                this%X_mat%diag(this%mixing_ratios_conc%num_cols),this%X_mat%sub(this%mixing_ratios_conc%num_cols-1)] !> Diagonal (current) and subdiagonal (upstream neighbor) elements of X [-]
            this%mix_conc_indices%cols(this%mixing_ratios_conc%num_cols)%col_1(2)=this%mixing_ratios_conc%num_cols-1 !> Second element: upstream water index (second-to-last target) [-]
            this%mix_conc_indices%cols(this%mixing_ratios_conc%num_cols)%col_1(3)=1 !> Third element: number of upstream waters (one neighbor) [-]
            this%mix_conc_indices%cols(this%mixing_ratios_conc%num_cols)%col_1(4)=0 !> Fourth element: number of downstream waters (none for last target) [-]
        end if
    end if !> End of time integration method selection (implicit vs explicit)
!> Verify that mixing ratios are non-negative and columns sum to 1
    do i=1,this%spatial_discr%num_targets
        do j=1,this%mixing_ratios_conc%cols(i)%dim
            if (this%mixing_ratios_conc%cols(i)%col_1(j) < 0d0) then
                write(*,'(A,I4,A,I4,A,ES15.7)') &
                    'ERROR: negative mixing ratio at target ', i, &
                    ', index ', j, ', value = ', this%mixing_ratios_conc%cols(i)%col_1(j)
                error stop "Negative mixing ratio detected"
            end if
        end do
        mix_sum = sum(this%mixing_ratios_conc%cols(i)%col_1)
        if (abs(mix_sum - 1d0) > 1d-6) then
            write(*,'(A,I4,A,ES14.6)') 'ERROR: mixing ratios column ', i, ' does not sum to 1: ', mix_sum
            error stop "Mixing ratios column does not sum to 1"
        end if
    end do
    call copy_real_array(this%mixing_ratios_R, this%mixing_ratios_R_init)
    call this%set_mix_react_indices() !> Set mixing target water indices from full mixing water indices [-]
end subroutine !> End of compute_mixing_ratios_Delta_t_homog subroutine 