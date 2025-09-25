subroutine compute_inverse_tridiag_matrix(this,tol,inv_mat)
    use metodos_sist_lin_m, only: Thomas, tridiag_matrix_c
    implicit none
    class(tridiag_matrix_c), intent(in) :: this
    real(kind=8), intent(in) :: tol
    real(kind=8), intent(out) :: inv_mat(:,:) !> must be allocated
    
    integer(kind=4) :: i,j
    real(kind=8), allocatable :: id_col(:)
    
    allocate(id_col(this%num_cols))
    
    do j=1,this%num_cols
        id_col=0d0
        id_col(j)=1d0
        call Thomas(this,id_col,tol,inv_mat(:,j))
    end do
end subroutine