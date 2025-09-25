subroutine Thomas(A,b,tol,x)
    !> Solves linear system of equations with tridiagonal matrix using Thomas algorithm
    
    !> A: tridiagonal matrix
    !> b: independent term
    !> x: solution of linear system

    use matrices_m, only: tridiag_matrix_c
    use vectors_m, only: inf_norm_vec_real
    implicit none
    class(tridiag_matrix_c), intent(in) :: A
    real(kind=8), intent(in) :: b(:)
    real(kind=8), intent(in) :: tol !> tolerance
    real(kind=8), intent(out) :: x(:) !> tiene que estar alocatado
    
    integer(kind=4) :: i,n
    real(kind=8), parameter :: epsilon=1d-16
    real(kind=8), allocatable :: c_star(:),d_star(:), matrix(:,:)
    
    n=size(b)
    allocate(c_star(n-1),d_star(n),matrix(n,n))
    c_star(1)=A%super(1)/A%diag(1)
    d_star(1)=b(1)/A%diag(1)
    do i=2,n-1
        c_star(i)=A%super(i)/(A%diag(i)-A%sub(i-1)*c_star(i-1))
        d_star(i)=(b(i)-A%sub(i-1)*d_star(i-1))/(A%diag(i)-A%sub(i-1)*c_star(i-1))
    end do
    d_star(n)=(b(n)-A%sub(i-1)*d_star(n-1))/(A%diag(i)-A%sub(i-1)*c_star(i-1))
    x(1)=d_star(n)
    do i=2,n
        x(i)=d_star(n-i+1)-c_star(n-i+1)*x(i-1)
    end do
    x=x(n:1:-1)
    matrix=0d0
    matrix(1,1:2)=[A%diag(1),A%super(1)]
    do i=2,n-1
        matrix(i,i)=A%diag(i)
        matrix(i,i-1)=A%sub(i-1)
        matrix(i,i+1)=A%super(i)
    end do
    matrix(n,(n-1):n)=[A%sub(n-1),A%diag(n)]
    if (inf_norm_vec_real(matmul(matrix,x)-b)>=tol) then
        print *, inf_norm_vec_real(matmul(matrix,x)-b)
        error stop "Thomas solution not accurate enough"
    end if
    deallocate(c_star,d_star,matrix)
end subroutine