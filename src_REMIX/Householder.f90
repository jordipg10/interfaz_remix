function Householder(x)
    use vectors_m
    use matrices_m, only : id_matrix,compute_det
    implicit none
    real(kind=8), intent(in) :: x(:)
    real(kind=8), allocatable :: Householder(:,:)
    integer(kind=4) :: n
    real(kind=8) :: beta
    real(kind=8), parameter :: epsilon=1d-6
    real(kind=8), allocatable :: u(:)
    n=size(x)
    allocate(Householder(n,n),u(n))
    beta=1d0/(p_norm_vec(x,2)*(abs(x(1))+p_norm_vec(x,2)))
    u(1)=sign(1d0,x(1))*(abs(x(1))+p_norm_vec(x,2))
    u(2:n)=x(2:n)
    Householder=id_matrix(n)-beta*outer_prod_vec(u,u)
    !if (abs(det(Householder))+epsilon<1d0) error stop "Not orthogonal"
end function