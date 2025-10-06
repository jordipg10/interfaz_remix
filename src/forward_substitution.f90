subroutine forward_substitution(L,b,x)
    implicit none
    real(kind=8), intent(in) :: L(:,:) !> lower triangular matrix
    real(kind=8), intent(in) :: b(:) !> vector
    real(kind=8), intent(out) :: x(:) !> solution of linear system
    integer(kind=4) :: j,k,n
    real(kind=8) :: sum
    n=size(b)
    x(1)=b(1)/L(1,1)
    do j=2,n
        sum=0d0
        do k=1,j-1
            sum=sum+L(j,k)*x(k)
        end do
        x(j)=(b(j)-sum)/L(j,j)
    end do
end subroutine 