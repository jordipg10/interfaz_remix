function inf_norm_mat(A) result(norm)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> m x n matrix
    real(kind=8) :: norm
            
    integer(kind=4) :: i,j,n,m
    real(kind=8) :: sum
    m=size(A,1)
    n=size(A,2)
    sum=0
    do i=1,n
        sum=sum+abs(A(1,i))
    end do
    norm=sum
    do i=2,m
        sum=0
        do j=1,n
            sum=sum+abs(A(i,j))
        end do
        if (sum>norm) then
            norm=sum
        end if
    end do
end function