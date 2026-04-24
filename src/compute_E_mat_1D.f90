!> \file compute_E_mat.f90
!> \brief Computes the E matrix for transient PDE discretization
!> \details This subroutine computes the E matrix used in the temporal discretization of transient PDEs.
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

subroutine compute_E_mat_1D(this,Delta_t,E_mat,E_mat_prev)
!> E=Delta_t(k)*inv(F)*T (tridiagonal, negative semi-definite)
!> rows sum = 0 if r=0
    use PDE_transient_m, only: PDE_1D_transient_c
    use time_discr_m, only: time_discr_homog_c, time_discr_heterog_c
    use arrays_m, only: tridiag_matrix_c
    implicit none !> Enforce explicit variable declarations
    
    class(PDE_1D_transient_c), intent(in) :: this !> Transient PDE object containing spatial discretization, time discretization, transport matrix T, and F matrix [-]
    real(kind=8), intent(in) :: Delta_t !> Time step size (can be homogeneous or heterogeneous) [-]
    type(tridiag_matrix_c), intent(out) :: E_mat !> E matrix: Δt·F⁻¹·T (tridiagonal structure with sub-diagonal, diagonal, and super-diagonal)
    type(tridiag_matrix_c), intent(out) :: E_mat_prev !> E matrix at previous time step (for time integration schemes) [-]
    
    !integer(kind=4), intent(in), optional :: k !> Time step index for heterogeneous time discretization (optional parameter) [-]
    
    integer(kind=4) :: j !> Loop counter for matrix elements (rows/columns) [-]
    integer(kind=4) :: n !> Number of spatial discretization nodes (targets) [-]
    
    n=this%spatial_discr%Num_targets !> Extract number of spatial nodes from spatial discretization object

    if (this%spatial_discr%adapt_ref.eq.1) then !> Check if adaptive refinement is enabled (spatial mesh adapts during simulation)
        call this%compute_trans_mat_PDE() !> Recompute transport matrix T if mesh has been refined (update advection and dispersion discretization)
        call this%compute_F_mat_PDE() !> Recompute F matrix if mesh has been refined (update accumulation terms)
    end if !> End adaptive refinement check
    E_mat%sub=Delta_t*this%trans_mat%sub !> Copy sub-diagonal of transport matrix T to E matrix (lower diagonal elements)
    E_mat%diag=Delta_t*this%trans_mat%diag !> Copy main diagonal of transport matrix T to E matrix (diagonal elements)
    E_mat%super=Delta_t*this%trans_mat%super !> Copy super-diagonal of transport matrix T to E matrix (upper diagonal elements)
    
    E_mat_prev%sub=Delta_t*this%trans_mat_prev%sub !> Copy sub-diagonal of transport matrix T to previous E matrix (lower diagonal elements)
    E_mat_prev%diag=Delta_t*this%trans_mat_prev%diag !> Copy main diagonal of transport matrix T to previous E matrix (diagonal elements)
    E_mat_prev%super=Delta_t*this%trans_mat_prev%super !> Copy super-diagonal of transport matrix T to previous E matrix (upper diagonal elements)
    ! select type (time_discr=>this%time_discr) !> Polymorphic select to determine time discretization type (homogeneous vs heterogeneous Δt)
    ! type is (time_discr_homog_c) !> Case: homogeneous time discretization (constant Δt throughout simulation)
    !     E_mat%sub=E_mat%sub*time_discr%Delta_t !> Multiply sub-diagonal by constant time step: Δt·T_sub
    !     E_mat%diag=E_mat%diag*time_discr%Delta_t !> Multiply main diagonal by constant time step: Δt·T_diag
    !     E_mat%super=E_mat%super*time_discr%Delta_t !> Multiply super-diagonal by constant time step: Δt·T_super
    ! type is (time_discr_heterog_c) !> Case: heterogeneous time discretization (variable Δt, changes with time step k)
    !     E_mat%sub=E_mat%sub*time_discr%Delta_t(k) !> Multiply sub-diagonal by k-th time step: Δt(k)·T_sub
    !     E_mat%diag=E_mat%diag*time_discr%Delta_t(k) !> Multiply main diagonal by k-th time step: Δt(k)·T_diag
    !     E_mat%super=E_mat%super*time_discr%Delta_t(k) !> Multiply super-diagonal by k-th time step: Δt(k)·T_super
    ! end select !> End time discretization type selection
    !print *, size(E_mat%sub),size(this%F_mat%diag) !> Commented: Debug print for checking array sizes
    do j=1,n-1 !> Loop over first n-1 nodes (interior and boundary nodes except last)
        E_mat%super(j)=E_mat%super(j)/this%F_mat_prev%diag(j) !> Divide super-diagonal by F diagonal at row j: (Δt·T_super)ⱼ/Fⱼ (inverse of F applied)
        E_mat%sub(j)=E_mat%sub(j)/this%F_mat_prev%diag(j+1) !> Divide sub-diagonal by F diagonal at row j+1: (Δt·T_sub)ⱼ/Fⱼ₊₁ (shifted index for sub-diagonal structure)
        E_mat%diag(j)=E_mat%diag(j)/this%F_mat_prev%diag(j) !> Divide main diagonal by F diagonal at row j: (Δt·T_diag)ⱼ/Fⱼ (inverse of F applied)
    end do !> End loop over first n-1 nodes
    E_mat%diag(n)=E_mat%diag(n)/this%F_mat_prev%diag(n) !> Divide last diagonal element by F diagonal at row n: (Δt·T_diag)ₙ/Fₙ (boundary node treatment)

    do j=1,n-1 !> Loop over first n-1 nodes for previous E matrix
        E_mat_prev%super(j)=E_mat_prev%super(j)/this%F_mat%diag(j) !> Divide super-diagonal by previous F diagonal at row j
        E_mat_prev%sub(j)=E_mat_prev%sub(j)/this%F_mat%diag(j+1) !> Divide sub-diagonal by previous F diagonal at row j+1
        E_mat_prev%diag(j)=E_mat_prev%diag(j)/this%F_mat%diag(j) !> Divide main diagonal by previous F diagonal at row j
    end do !> End loop over first n-1 nodes for previous E matrix
    E_mat_prev%diag(n)=E_mat_prev%diag(n)/this%F_mat%diag(n) !> Divide last diagonal element of previous E matrix by previous F diagonal at row n
end subroutine 