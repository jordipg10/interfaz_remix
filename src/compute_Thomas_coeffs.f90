subroutine compute_Thomas_coeffs(A,b,tol,c_tilde,d_tilde)
    !> Computes coefficients of Thomas algorithm
    
    !> A: tridiagonal matrix
    !> c_tilde,d_tilde: coefficients

    use metodos_sist_lin_m
    implicit none
    
    class(tridiag_matrix_c), intent(in) :: A
    real(kind=8), intent(in) :: b(:)
    real(kind=8), intent(in) :: tol
    real(kind=8), intent(out) :: c_tilde(:),d_tilde(:) !> tiene que estar alocatado
    
    integer(kind=4) :: i,n
    real(kind=8) :: denom
    
    c_tilde(1)=A%super(1)/A%diag(1)
    d_tilde(1)=b(1)/A%diag(1)
    do i=2,n-1
        denom=A%diag(i)-A%sub(i-1)*c_tilde(i-1)
        if (abs(denom)<tol) then
            error stop "Singularity in Thomas algorithm"
        end if
        c_tilde(i)=A%super(i)/denom
        d_tilde(i)=(b(i)-A%sub(i-1)*d_tilde(i-1))/denom
    end do
    denom=A%diag(n)-A%sub(n-1)*c_tilde(n-1)
    if (abs(denom)<tol) then
        error stop "Singularity in Thomas algorithm"
    end if
    d_tilde(n)=(b(n)-A%sub(n-1)*d_tilde(n-1))/denom
end subroutine