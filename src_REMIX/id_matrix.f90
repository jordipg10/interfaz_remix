function id_matrix(n)
!> Identity nxn matrix
    implicit none
    integer(kind=4), intent(in) :: n
    integer(kind=4) :: i
    real(kind=8) :: id_matrix(n,n)
    id_matrix=0d0
    forall (i=1:n)
        id_matrix(i,i)=1d0
    end forall
end function