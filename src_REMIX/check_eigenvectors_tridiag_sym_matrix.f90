subroutine check_eigenvectors_tridiag_sym_matrix(this,tolerance)
    use matrices_m, only: tridiag_sym_matrix_c
    use vectors_m
    implicit none
    class(tridiag_sym_matrix_c), intent(in) :: this !> tridiagonal symmetric matrix object
    real(kind=8), intent(in) :: tolerance
    
    real(kind=8), allocatable :: A_lambda(:,:,:),A_lambda_P(:),prod(:,:)
    integer(kind=4) :: i,j,n
    
        n=size(this%eigenvalues)
        allocate(A_lambda(n,n,n))
        do i=1,n
            A_lambda(:,:,i)=0d0
            A_lambda(1,1:2,i)=[this%diag(1)-this%eigenvalues(i),this%sub(1)]
            do j=2,n-1
                A_lambda(j,j-1,i)=this%sub(j-1)
                A_lambda(j,j+1,i)=this%sub(j)
                A_lambda(j,j,i)=this%diag(j)-this%eigenvalues(i)
            end do
            A_lambda(n,n-1:n,i)=[this%sub(n-1),this%diag(n)-this%eigenvalues(i)]
        end do
        do i=1,n
            A_lambda_P=matmul(A_lambda(:,:,i),this%eigenvectors(:,i))
            if (inf_norm_vec_real(A_lambda_P)>=tolerance) then
                print *, "Error in eigenvector", i, inf_norm_vec_real(A_lambda_P)
                error stop
            end if
        end do
        prod=matmul(this%eigenvectors,transpose(this%eigenvectors))
        !if ((abs(det(this%eigenvectors))-1d0)>=tolerance) then
        !    print *, "abs(det(P)) is not 1"
        !    error stop
        !end if
end subroutine