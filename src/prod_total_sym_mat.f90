function prod_total_sym_mat(A,y0,b,time) result(y)
!> dy/dt=-Ay+b
    use matrices_m
    implicit none
    class(sq_matrix_c), intent(in) :: A
    real(kind=8), intent(in) :: y0(:)
    real(kind=8), intent(in) :: b(:)
    real(kind=8), intent(in) :: time
    real(kind=8), allocatable :: y(:)
    
    real(kind=8), allocatable :: Q_lambda(:,:),x(:),Pt_x(:)
    integer(kind=4) :: i,n,time_step
    real(kind=8), parameter :: epsilon=1d-6
    
    n=size(A%eigenvalues)
    allocate(Q_lambda(2,n),x(2*n),Pt_x(2*n))
    !Q_lambda(1,:)=1d0
    !if (present(k)) then
    !>    time_step=k
    !else
    !>    time_step=time_discr_obj%Num_time
    !end if
    !select type (time=>time_discr_obj)
    !type is (time_discr_homog_c)
        do i=1,n
            if (abs(A%eigenvalues(i))<epsilon) then
                Q_lambda(2,i)=0d0
            else
                Q_lambda(2,i)=(1d0-exp(-A%eigenvalues(i)*time))/A%eigenvalues(i)
            end if
        end do
        Q_lambda(1,:)=exp(-A%eigenvalues*time)
    !type is (time_discr_heterog_c)
        do i=1,n
            if (abs(A%eigenvalues(i))<epsilon) then
                Q_lambda(2,i)=0d0
            else
                Q_lambda(2,i)=(1d0-exp(-A%eigenvalues(i)*time))/A%eigenvalues(i)
            end if
        end do
        Q_lambda(1,:)=exp(-A%eigenvalues*time)
    !end select
    x(1:n)=y0
    x(n+1:2*n)=b
    Pt_x(1:n)=matmul(transpose(A%eigenvectors),x(1:n))
    Pt_x(n+1:2*n)=matmul(transpose(A%eigenvectors),x(n+1:2*n))
    y=matmul(A%eigenvectors,Q_lambda(1,:)*Pt_x(1:n))+matmul(A%eigenvectors,Q_lambda(2,:)*Pt_x(n+1:2*n))
end function