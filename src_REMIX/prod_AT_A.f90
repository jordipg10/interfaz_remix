function prod_AT_A(A) result(AT_A) !> AT_A=A^T * A
    implicit none
    real(kind=8), intent(in) :: A(:,:)
    real(kind=8), allocatable :: AT_A(:,:)
    
    integer(kind=4) :: n,m,i,j
    
    n=size(A,1)
    m=size(A,2)
    
    allocate(AT_A(m,m))

    do i=1,m
        do j=1,m
            AT_A(i,j)=dot_product(A(:,i),A(:,j))
        end do
    end do
end function