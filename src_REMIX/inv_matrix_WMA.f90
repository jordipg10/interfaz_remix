!> Inverse of square matrix using LU decomposition from first element
subroutine inv_matrix_WMA(A,tol,inv)
    use matrices_m
    use vectors_m
    use metodos_sist_lin_m
    implicit none
    real(kind=8), intent(in) :: A(:,:)
    real(kind=8), intent(in) :: tol
    real(kind=8), intent(out) :: inv(:,:)
    
    integer(kind=4) :: n,j,i
    real(kind=8), parameter :: epsilon=1d-6
    real(kind=8), allocatable :: id(:,:), prod_A_invA(:,:), id_col(:), inv_col(:)
    logical :: nzdiag
    
    if (size(A,1)/=size(A,2)) error stop "Matrix must be square (inv_matrix_WMA)"
    
    !if (det(A)>=tol) error stop "Matrix is not invertible"
    
    call inv_matrix(A,tol,inv) 
end subroutine