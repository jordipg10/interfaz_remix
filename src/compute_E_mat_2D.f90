!> \file compute_E_mat_2D.f90
!> \brief Computes the E matrix for transient PDE discretization in 2D
!> \details This subroutine computes the E matrix used in the temporal discretization of transient PDEs in 2D.
!> The E matrix combines spatial discretization (transport matrix T), temporal discretization (Δt), 
!> and the inverse of the F matrix (which contains accumulation/storage terms).
!>
!> Mathematical formulation:
!> \f[
!> \mathbf{E} = \Delta t \cdot \mathbf{F}^{-1} \cdot \mathbf{T}
!> \f]
!> where:
!> - \f$ \mathbf{T} \f$ = transport matrix (tridiagonal, contains advection and dispersion terms)
!> - \f$ \Delta t \f$ = time step size [T]
!> - \f$ \mathbf{F} \f$ = diagonal matrix containing accumulation/storage terms [-]
!> - \f$ \mathbf{E} \f$ = combined matrix for time-stepping (tridiagonal, negative semi-definite)
!>
!> Properties:
!> - E is tridiagonal (inherits structure from T)
!> - E is negative semi-definite
!> - Row sums equal 0 when source/sink term r = 0 (conservation property)
!>
!> The computation involves:
!> 1. Copy transport matrix T to E
!> 2. Multiply by time step Δt (homogeneous or heterogeneous)
!> 3. Divide by diagonal elements of F (element-wise inverse)
!>
!> \param[in,out] this Transient PDE object containing spatial discretization, time discretization, transport matrix, and F matrix
!> \param[out] E_mat E matrix: Δt·F⁻¹·T for time-stepping in PDE solver (tridiagonal structure)
!> \param[in] k Time step index for heterogeneous time discretization (optional, used when Δt varies with time) [-]

subroutine compute_E_mat_2D(this,Delta_t,E_mat,E_mat_prev)
!> E=Delta_t(k)*inv(F)*T (tridiagonal, negative semi-definite)
!> rows sum = 0 if r=0
    use PDE_transient_m, only: PDE_2D_transient_c
    use time_discr_m, only: time_discr_homog_c, time_discr_heterog_c
    use arrays_m, only: pentadiag_matrix_c
    use spatial_discr_2D_m, only: mesh_2D_Euler_homog_c
    implicit none !> Enforce explicit variable declarations
    
    class(PDE_2D_transient_c), intent(in) :: this !> Transient PDE object containing spatial discretization, time discretization, transport matrix T, and F matrix [-]
    real(kind=8), intent(in) :: Delta_t !> Time step size (can be homogeneous or heterogeneous) [-]
    type(pentadiag_matrix_c), intent(out) :: E_mat !> E matrix: Δt·F⁻¹·T (pentadiagonal structure with sub-diagonal, diagonal, and super-diagonal)
    type(pentadiag_matrix_c), intent(out) :: E_mat_prev !> E matrix at previous time step (for time integration schemes) [-]
    
    !integer(kind=4), intent(in), optional :: k !> Time step index for heterogeneous time discretization (optional parameter) [-]
    
    integer(kind=4) :: j !> Loop counter for matrix elements (rows/columns) [-]
    integer(kind=4) :: n !> Number of spatial discretization nodes (targets) [-]
    integer(kind=4) :: nx, ny !> Number of cells in x and y directions (for 2D mesh) [-]
    integer(kind=4) :: n_sub !> Number of sub-diagonal elements (x-direction band size) [-]
    integer(kind=4) :: n_sub2 !> Number of second sub-diagonal elements (y-direction band size) [-]
    
    n=this%spatial_discr%Num_targets !> Extract number of spatial nodes from spatial discretization object

    if (this%spatial_discr%adapt_ref.eq.1) then !> Check if adaptive refinement is enabled (spatial mesh adapts during simulation)
        call this%compute_trans_mat_PDE() !> Recompute transport matrix T if mesh has been refined (update advection and dispersion discretization)
        call this%compute_F_mat_PDE() !> Recompute F matrix if mesh has been refined (update accumulation terms)
    end if !> End adaptive refinement check

    call E_mat%allocate_array(n)
    call E_mat_prev%allocate_array(n)

    !> Correct off-diagonal band sizes for 2D structured grid
    select type (mesh=>this%spatial_discr)
    type is (mesh_2D_Euler_homog_c)
        nx     = mesh%Num_cells_x
        ny     = mesh%Num_cells_y
        n_sub  = ny * (nx - 1)      !> x-direction band size
        n_sub2 = (ny - 1) * nx      !> y-direction band size
        !> E_mat
        if (allocated(E_mat%sub))    deallocate(E_mat%sub)
        if (allocated(E_mat%super))  deallocate(E_mat%super)
        if (allocated(E_mat%sub2))   deallocate(E_mat%sub2)
        if (allocated(E_mat%super2)) deallocate(E_mat%super2)
        allocate(E_mat%sub(n_sub),    E_mat%super(n_sub))
        allocate(E_mat%sub2(n_sub2),  E_mat%super2(n_sub2))
        E_mat%sub = 0.0d0 ; E_mat%super  = 0.0d0
        E_mat%sub2= 0.0d0 ; E_mat%super2 = 0.0d0
        !> E_mat_prev
        if (allocated(E_mat_prev%sub))    deallocate(E_mat_prev%sub)
        if (allocated(E_mat_prev%super))  deallocate(E_mat_prev%super)
        if (allocated(E_mat_prev%sub2))   deallocate(E_mat_prev%sub2)
        if (allocated(E_mat_prev%super2)) deallocate(E_mat_prev%super2)
        allocate(E_mat_prev%sub(n_sub),    E_mat_prev%super(n_sub))
        allocate(E_mat_prev%sub2(n_sub2),  E_mat_prev%super2(n_sub2))
        E_mat_prev%sub = 0.0d0 ; E_mat_prev%super  = 0.0d0
        E_mat_prev%sub2= 0.0d0 ; E_mat_prev%super2 = 0.0d0
    end select

    E_mat%sub    = Delta_t*this%trans_mat%sub
    E_mat%diag   = Delta_t*this%trans_mat%diag
    E_mat%super  = Delta_t*this%trans_mat%super
    E_mat%sub2   = Delta_t*this%trans_mat%sub2
    E_mat%super2 = Delta_t*this%trans_mat%super2
    
    E_mat_prev%sub    = Delta_t*this%trans_mat_prev%sub
    E_mat_prev%diag   = Delta_t*this%trans_mat_prev%diag
    E_mat_prev%super  = Delta_t*this%trans_mat_prev%super
    E_mat_prev%sub2   = Delta_t*this%trans_mat_prev%sub2
    E_mat_prev%super2 = Delta_t*this%trans_mat_prev%super2
    
    do j=1,n_sub
        E_mat%super(j) = E_mat%super(j)/this%F_mat_prev%diag(j)
        E_mat%sub(j)   = E_mat%sub(j)  /this%F_mat_prev%diag(j+1)
    end do
    do j=1,n_sub2
        E_mat%super2(j) = E_mat%super2(j)/this%F_mat_prev%diag(j)
        E_mat%sub2(j)   = E_mat%sub2(j)  /this%F_mat_prev%diag(j+1)
    end do
    do j=1,n
        E_mat%diag(j) = E_mat%diag(j)/this%F_mat_prev%diag(j)
    end do

    do j=1,n_sub
        E_mat_prev%super(j) = E_mat_prev%super(j)/this%F_mat%diag(j)
        E_mat_prev%sub(j)   = E_mat_prev%sub(j)  /this%F_mat%diag(j+1)
    end do
    do j=1,n_sub2
        E_mat_prev%super2(j) = E_mat_prev%super2(j)/this%F_mat%diag(j)
        E_mat_prev%sub2(j)   = E_mat_prev%sub2(j)  /this%F_mat%diag(j+1)
    end do
    do j=1,n
        E_mat_prev%diag(j) = E_mat_prev%diag(j)/this%F_mat%diag(j)
    end do
end subroutine 