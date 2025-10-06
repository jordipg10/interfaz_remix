!> Cholesky decomposition
subroutine Cholesky(A,L)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> positive definite, symmetric matrix
    real(kind=8), intent(out) :: L(:,:) !> lower triangular matrix with positive diagonal terms
    
    real(kind=8), allocatable :: prod_L_Lt(:,:) !> L*L^T
    integer(kind=4) :: i,j,k,n
    real(kind=8) :: sum1, sum2
    n=size(A,1)
    L=0d0
    do j=1,n
        if (j == 1) then
            L(j,j)=sqrt(A(j,j))
        else
            sum1=0d0
            do k=1,j-1
                sum1=sum1+L(j,k)**2
            end do
            L(j,j)=sqrt(A(j,j)-sum1)
        end if
        do i=j+1,n
            sum2=0d0
            do k=1,j-1
                sum2=sum2+L(i,k)*L(j,k)
            end do
            L(i,j)=(A(i,j)-sum2)/L(j,j)
        end do
    end do
    prod_L_Lt=matmul(L,transpose(L))
    if (all(A/=prod_L_Lt)) error stop "Wrong solution"
end subroutine