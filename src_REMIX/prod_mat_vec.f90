function prod_mat_vec(A,b) result(x) !> Ab=x
    use matrices_m, only: array_c, tridiag_matrix_c
    implicit none
    class(array_c), intent(in) :: A
    real(kind=8), intent(in) :: b(:)
    real(kind=8), allocatable :: x(:) 
    
    integer(kind=4) :: i,j,n,m
    real(kind=8) :: sum
    
    select type (A)
    !type is (tridiag_matrix_vec_c)
    !>    n=size(A%diag)
    !>    if (size(b)/=2*n) error stop "Dimension error in b"
    !>    if (n.eqv.1) then
    !>        x=A%diag(1)*b(1)+A%vector(1)*b(n+1)
    !>    else if (n.eqv.2) then
    !>        x(1)=A%diag(1)*b(1)+A%super(1)*b(2)+A%vector(1)*b(n+1)
    !>        x(2)=A%sub(1)*b(1)+A%diag(2)*b(2)+A%vector(2)*b(n+2)
    !>    else
    !>        x(1)=A%diag(1)*b(1)+A%super(1)*b(2)+A%vector(1)*b(n+1)
    !>        do i=2,n-1
    !>            x(i)=A%sub(i-1)*b(i-1)+A%diag(i)*b(i)+A%super(i)*b(i+1)
    !>            x(i)=x(i)+A%vector(i)*b(n+i)
    !>        end do
    !>        x(n)=A%sub(n-1)*b(n-1)+A%diag(n)*b(n)+A%vector(n)*b(2*n)
    !>    end if
    type is (tridiag_matrix_c)
        n=size(A%diag)
        allocate(x(n))
        if (size(b)/=n) error stop "Dimension error in b"
        x(1)=A%diag(1)*b(1)+A%super(1)*b(2)
        do i=2,n-1
            x(i)=A%sub(i-1)*b(i-1)+A%diag(i)*b(i)+A%super(i)*b(i+1)
        end do
        x(n)=A%sub(n-1)*b(n-1)+A%diag(n)*b(n)
    end select
end function 