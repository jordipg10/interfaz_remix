!subroutine Gram_Schmidt_vec(v,u)
!>    use vectores_m, only : proy_ortog
!>    implicit none
!>    real(kind=8), intent(in) :: v(:)
!>    real(kind=8), intent(out) :: u(:) !> vector ortogonal
!>    
!>    integer(kind=4) :: i,,n
!>    real(kind=8) :: sum_proy_ortog
!>    n=size(v)
!>    !allocate(sum_proy_ortog(n))
!>    u(1)=v(1)
!>    sum_proy_ortog=0d0
!>    do i=1,n-1
!>        sum_proy_ortog=sum_proy_ortog+proy_ortog(v(i),u())
!>    end do
!>        u(1:n,k)=v(1:n,k)-sum_proy_ortog
!>    end do
!end subroutine