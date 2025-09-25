!> This subroutine performs backward substitution in LU decomposition of matrix
subroutine backward_substitution(U,b,x)
    implicit none
    real(kind=8), intent(in) :: U(:,:) !> upper triangular matrix
    real(kind=8), intent(in) :: b(:) !> independent term
    real(kind=8), intent(out) :: x(:) !> solution of linear system
    
    integer(kind=4) :: j,k,n
    real(kind=8) :: sum
    
    n=size(b)
    x(n)=b(n)/U(n,n)
    do j=1,n-1
        sum=0d0
        do k=n-j+1,n
            sum=sum+U(n-j,k)*x(k)
        end do
        x(n-j)=(b(n-j)-sum)/U(n-j,n-j)
    end do
end subroutine