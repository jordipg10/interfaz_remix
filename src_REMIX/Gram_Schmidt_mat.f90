subroutine Gram_Schmidt_mat(v,u)
    use vectors_m, only : proy_ortog
    implicit none
    real(kind=8), intent(in) :: v(:,:)
    real(kind=8), intent(out) :: u(:,:) !> base ortogonal
    integer(kind=4) :: k,j,n
    real(kind=8), allocatable :: sum_proy_ortog(:)
    n=size(v,1)
    allocate(sum_proy_ortog(n))
    u(1:n,1)=v(1:n,1)
    do k=2,n
        sum_proy_ortog=0d0
        do j=1,k-1
            sum_proy_ortog=sum_proy_ortog+proy_ortog(v(1:n,k),u(1:n,j))
        end do
        u(1:n,k)=v(1:n,k)-sum_proy_ortog
    end do
end subroutine