!> LU decomposition
subroutine LU(A,tol,L,U,error)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> square matrix
    real(kind=8), intent(in) :: tol !> tolerance
    real(kind=8), intent(out) :: L(:,:), U(:,:)
    logical, intent(out) :: error
    
    integer(kind=4) :: n,i,j
    real(kind=8) :: factor
    n=size(A,1)
    U=A
    L=0d0
    error=.false.
    do j=1,n
        if (abs(U(j,j))<tol) then
            error=.true.
            !error stop "Diagonal term is zero in LU decomposition"
        end if
        L(j,j)=1d0
        do i=j+1,n
            factor=U(i,j)/U(j,j)
            L(i,j)=factor
            U(i,1:n)=U(i,1:n)-factor*U(j,1:n)
        end do
    end do

end subroutine LU