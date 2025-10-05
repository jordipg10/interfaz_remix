subroutine non_zeros(A,n,ind)
    implicit none
    real(kind=8), intent(in) :: A(:,:)
    integer(kind=4), intent(out) :: n
    integer(kind=4), intent(out) :: ind(:,:)
    
    integer(kind=4) :: i,j
    integer(kind=4), allocatable :: ind_aux(:)
    real(kind=8), parameter :: epsilon=1d-9
    
    n=0
    !ind_aux=[]
    do i=1,size(A,1)
        do j=1,size(A,2)
            if (abs(A(i,j))>=epsilon) then
                n=n+1
                ind_aux=[ind_aux,i,j]
            end if
        end do
    end do
    ind=reshape(ind_aux,[n,2])
end subroutine