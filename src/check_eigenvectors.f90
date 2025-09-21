subroutine check_eigenvectors(A,lambda,v,tolerance)
    use matrices_m
    use vectors_m
    implicit none
    real(kind=8), intent(in) :: A(:,:)
    real(kind=8), intent(in) :: lambda(:)
    real(kind=8), intent(in) :: v(:,:)
    real(kind=8), intent(in) :: tolerance
    
    real(kind=8), allocatable :: A_lambda(:,:,:),A_lambda_v(:),inv_v(:,:),prod(:,:)
    integer(kind=4) :: i,j,n
    n=size(lambda)
    allocate(A_lambda(n,n,n),inv_v(n,n))
    do i=1,n
        A_lambda(:,:,i)=A
        do j=1,n
            A_lambda(j,j,i)=A_lambda(j,j,i)-lambda(i)
        end do
    end do
    do i=1,n
        A_lambda_v=matmul(A_lambda(:,:,i),v(:,i))
        if (inf_norm_vec_real(A_lambda_v)>=tolerance) then
            print *, "Error in eigenvector", i, inf_norm_vec_real(A_lambda_v)
            !error stop
        end if
    end do
    !call inv_matrix(v,inv_v)
    !prod=matmul(v,inv_v)
    !do i=1,n
    !>    if (i.eqv.1) then
    !>        if (abs(prod(i,i)-1d0)>=tolerance .or. inf_norm_vec(prod(i+1:n,i))>=tolerance) print *, "Error in first eigenvector"
    !>    else if (i>1 .and. i<n) then
    !>        if (abs(prod(i,i)-1d0)>=tolerance .or. inf_norm_vec(prod(1:i-1,i))>=tolerance .or. inf_norm_vec(prod(i+1:n,i))>=tolerance) then
    !>            print *, "Error in eigenvector", i
    !>            error stop
    !>        end if
    !>    else
    !>        if (abs(prod(i,i)-1d0)>=tolerance .or. inf_norm_vec(prod(1:i-1,i))>=tolerance) print *, "Error in last eigenvector"
    !>    end if
    !end do
end subroutine