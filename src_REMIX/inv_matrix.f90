!> Inverse of square matrix using LU decomposition
subroutine inv_matrix(A,tol,inv)
    use matrices_m, only : id_matrix, compute_det
    use vectors_m, only : inf_norm_vec_real
    use metodos_sist_lin_m, only : LU_lin_syst
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> square matrix 
    real(kind=8), intent(in) :: tol !> tolerance for determinant
    real(kind=8), intent(out) :: inv(:,:) !> inverse matrix of A (must be allocated)
    
    integer(kind=4) :: n,j,i,err
    real(kind=8) :: det
    real(kind=8), parameter :: epsilon=1d-6
    real(kind=8), allocatable :: id(:,:), prod_A_invA(:,:), id_col(:), inv_col(:)
    logical :: nzdiag,error
    
    !print *, size(A,1), size(A,2)
    if (size(A,1)/=size(A,2)) then
        error stop "Matrix must be square (inv_matrix)"
    end if
    if (size(A,1)==1) then
        inv(1,1)=1d0/A(1,1)
        !return
    else
        call compute_det(A,tol,det,error)
        if ((error .eqv. .true.) .or. (abs(det)<tol)) then
            error stop "Matrix is not invertible"
        end if
        n=size(A,1)
        if (n.eq.2) then
            inv(1,1)=A(2,2)
            inv(1,2)=-A(1,2)
            inv(2,1)=-A(2,1)
            inv(2,2)=A(1,1)
            inv=inv/det
        else
            allocate(inv_col(n))
            id=id_matrix(n)
            !nzdiag=.true.
            !do i=1,n
            !>    if (abs(A(i,i))<epsilon) then
            !>        nzdiag=.false.
            !>        exit
            !>    end if
            !end do
            !if (nzdiag.eqv..true.) then
            !>    do j=1,n
            !>        id_col=id(1:n,j)
            !>        call LU_lin_syst(A,id_col,inv_col) !> LU decomposition
            !>        inv(1:n,j)=inv_col
            !>    end do
            !else
            !>    do j=1,n
            !>        id_col=id(1:n,j)
            !>        call Gauss_Jordan(A,id_col,inv_col) !> Gauss-Jordan
            !>        inv(1:n,j)=inv_col
            !>    end do
            !end if
            do j=1,n
                id_col=id(1:n,j)
                !call Gauss_Jordan(A,id_col,tol,inv_col,err) !> Gauss-Jordan
                call LU_lin_syst(A,id_col,tol,inv_col)
                !if (err.eqv.1) then
                !    error stop "Singular equation in Gauss-Jordan"
                !end if
                inv(1:n,j)=inv_col
            end do
            prod_A_invA=matmul(A,inv)
            do i=1,n
                if (inf_norm_vec_real(prod_A_invA(i,:)-id(i,:))>=tol) then
                    error stop "Error in inverse matrix"
                end if
            end do
        end if
    end if
end subroutine