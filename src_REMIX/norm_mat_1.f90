!> Computes $L_1$ norm of a matrix
function norm_mat_1(A) result(norm)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> m x n matrix
    real(kind=8) :: norm
            
    integer(kind=4) :: i,j,n,m
    real(kind=8) :: sum
    m=size(A,1)
    n=size(A,2)
    sum=0
    do i=1,m
        sum=sum+abs(A(i,1))
    end do
    norm=sum
    do j=2,n
        sum=0
        do i=1,m
            sum=sum+abs(A(i,j))
        end do
        if (sum>norm) then
            norm=sum
        end if
    end do
end function