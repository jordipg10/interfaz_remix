subroutine Thomas_Toeplitz(a,b,c,d,x)
    !> Solves linear system of equations with tridiagonal Toeplitz matrix using Thomas algorithm
    
    !> a: subdiagonal term
    !> b: diagonal term
    !> c: superdiagoal term
    !> d: independent term
    !> x: solution of linear system
    
    implicit none
    real(kind=8), intent(in) :: a,b,c
    real(kind=8), intent(in) :: d(:)
    real(kind=8), intent(out) :: x(:)
    integer(kind=4) :: i,n
    real(kind=8), allocatable :: c_star(:),d_star(:)
    n=size(d)
    allocate(c_star(n-1),d_star(n))
    c_star(1)=c/b
    d_star=d(1)/b
    do i=2,n-1
        c_star(i)=c/(b-a*c_star(i-1))
        d_star(i)=(d(i)-a*d_star(i-1))/(b-a*c_star(i-1))
    end do
    d_star(n)=(d(n)-a*d_star(n-1))/(b-a*c_star(i-1))
    x(1)=d_star(n)
    do i=2,n
        x(i)=d_star(n-i+1)-c_star(n-i+1)*x(i-1)
    end do
    x=x(n:1:-1)
end subroutine Thomas_Toeplitz
    