subroutine LU_lin_syst(A,b,tol,x) !> Ax=b
    use vectors_m, only : inf_norm_vec_real
    use matrices_m, only : LU, compute_det
    use metodos_sist_lin_m, only : forward_substitution, backward_substitution
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> square matrix (A=LU)
    real(kind=8), intent(in) :: b(:) !> vector
    real(kind=8), intent(in) :: tol !> tolerance for solution
    real(kind=8), intent(out) :: x(:) !> solution of linear system (must be already allocated)
    
    
    real(kind=8), allocatable :: L(:,:), U(:,:), y(:)
    real(kind=8) :: det
    integer(kind=4) :: n
    logical :: error
    
    n=size(b)
    call compute_det(A,tol,det,error)
    
    if (size(A,1)/=n .or. size(A,2)/=n) then
        error stop "Wrong dimensions in LU_lin_syst"
    else if (ABS(A(1,1))<tol) then
        error stop "A(1,1)=0 in LU_lin_syst"
    else if (abs(det)<tol) then
        error stop "Zero determinant in LU_lin_syst"
    end if
    
    if (n.eq.2) then
        x(2)=(A(1,1)*b(2)-A(2,1)*b(1))/det
        x(1)=(b(1)-A(1,2)*x(2))/A(1,1)
    else
        allocate(L(n,n),U(n,n),y(n))
        call LU(A,tol,L,U,error)
        if (error .eqv. .true.) then
            error stop
        end if
        call forward_substitution(L,b,y)
        call backward_substitution(U,y,x)
        if (inf_norm_vec_real(matmul(A,x)-b) .ge. tol) then
            !print *, "Wrong solution in LU_lin_syst"
            print *, "Residual: ", inf_norm_vec_real(matmul(A,x)-b)
            error stop "Wrong solution in LU_lin_syst"
        end if
    end if

end subroutine