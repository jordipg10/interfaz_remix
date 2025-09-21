subroutine compute_b_vec_lin_syst(this,theta,conc_old,b_vec,k)
!> A*c^(k+1)=b
!> A=Id-theta*E
!> B=(Id+(1-theta)*E) (tridiagonal, negative semi-definite)
!> b=B*c^k+f
    use PDE_transient_m, only: PDE_1D_transient_c, diag_matrix_c
    implicit none
    
    class(PDE_1D_transient_c), intent(in) :: this
    real(kind=8), intent(in) :: theta
    real(kind=8), intent(in) :: conc_old(:)
    real(kind=8), intent(out) :: b_vec(:) !> must be allocated
    integer(kind=4), intent(in), optional :: k
    
    integer(kind=4) :: i,n
    
    n=size(conc_old)
    if (n/=this%spatial_discr%Num_targets) error stop "Dimension error in subroutine 'compute_b_vec_lin_syst'"
    
    !call X_mat%allocate_array(n)
    !call this%compute_X_mat(theta,k)
    
    b_vec(1)=this%X_mat%diag(1)*conc_old(1)+this%X_mat%super(1)*conc_old(2)
    do i=2,n-1
        b_vec(i)=this%X_mat%sub(i-1)*conc_old(i-1)+this%X_mat%diag(i)*conc_old(i)+this%X_mat%super(i)*conc_old(i+1)
    end do
    b_vec(n)=this%X_mat%sub(n-1)*conc_old(n-1)+this%X_mat%diag(n)*conc_old(n)
    b_vec=b_vec+this%f_vec
end subroutine
    