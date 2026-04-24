!> \file metodos_sist_lin_m.f90
!> \brief Linear system solvers module.
!> \details
!> Provides routines for solving dense and tridiagonal linear systems
!> \f$ Ax = b \f$ including matrix inversion via LU decomposition,
!> Thomas algorithm (tridiagonal), Gauss-Jordan elimination, and
!> iterative methods (Jacobi, Gauss-Seidel).
!>
!> \see arrays_m, vectors_m
!> \author Jordi
!> \date Unknown
!> \ingroup algebra

!> \brief Linear system solvers module.
module metodos_sist_lin_m
    use arrays_m, only: id_matrix, tridiag_matrix_c, compute_det, forward_substitution, &
        backward_substitution, LU
    use vectors_m, only: inf_norm_vec_real
    implicit none
    save
    contains
    !> Inverse of square matrix using LU decomposition
subroutine inv_matrix(A,tol,inv)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> square matrix 
    real(kind=8), intent(in) :: tol !> tolerance for determinant
    real(kind=8), intent(out) :: inv(:,:) !> inverse matrix of A (must be allocated)
    
    integer(kind=4) :: n,j,i
    real(kind=8) :: det
    real(kind=8), allocatable :: id(:,:), prod_A_invA(:,:), id_col(:), inv_col(:)
    logical :: error
    
    if (size(A,1)/=size(A,2)) then
        error stop "Matrix must be square (inv_matrix)"
    end if
    if (size(A,1)==1) then
        inv(1,1)=1d0/A(1,1)
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
            do j=1,n
                id_col=id(1:n,j)
                call LU_lin_syst(A,id_col,tol,inv_col)
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

        subroutine Thomas(A,b,tol,x)
        !> Solves linear system of equations with tridiagonal matrix using Thomas algorithm
        
        !> A: tridiagonal matrix
        !> b: independent term
        !> x: solution of linear system
    
        class(tridiag_matrix_c), intent(in) :: A
        real(kind=8), intent(in) :: b(:)
        real(kind=8), intent(in) :: tol !> tolerance
        real(kind=8), intent(out) :: x(:) !> tiene que estar alocatado
        
        integer(kind=4) :: i,n
        real(kind=8), parameter :: epsilon=1d-16
        real(kind=8), allocatable :: c_star(:),d_star(:), matrix(:,:)
        
        n=size(b)
        allocate(c_star(n-1),d_star(n),matrix(n,n))
        c_star(1)=A%super(1)/A%diag(1)
        d_star(1)=b(1)/A%diag(1)
        do i=2,n-1
            c_star(i)=A%super(i)/(A%diag(i)-A%sub(i-1)*c_star(i-1))
            d_star(i)=(b(i)-A%sub(i-1)*d_star(i-1))/(A%diag(i)-A%sub(i-1)*c_star(i-1))
        end do
        d_star(n)=(b(n)-A%sub(i-1)*d_star(n-1))/(A%diag(i)-A%sub(i-1)*c_star(i-1))
        x(1)=d_star(n)
        do i=2,n
            x(i)=d_star(n-i+1)-c_star(n-i+1)*x(i-1)
        end do
        x=x(n:1:-1)
        matrix=0d0
        matrix(1,1:2)=[A%diag(1),A%super(1)]
        do i=2,n-1
            matrix(i,i)=A%diag(i)
            matrix(i,i-1)=A%sub(i-1)
            matrix(i,i+1)=A%super(i)
        end do
        matrix(n,(n-1):n)=[A%sub(n-1),A%diag(n)]
        if (inf_norm_vec_real(matmul(matrix,x)-b)>=tol) then
            print *, inf_norm_vec_real(matmul(matrix,x)-b)
            error stop "Thomas solution not accurate enough"
        end if
        deallocate(c_star,d_star,matrix)
    end subroutine

    subroutine Thomas_Toeplitz(a,b,c,d,x)
        !> Solves linear system of equations with tridiagonal Toeplitz matrix using Thomas algorithm
        
        !> a: subdiagonal term
        !> b: diagonal term
        !> c: superdiagoal term
        !> d: independent term
        !> x: solution of linear system
        
        real(kind=8), intent(in) :: a,b,c
        real(kind=8), intent(in) :: d(:)
        real(kind=8), intent(out) :: x(:)
        integer(kind=4) :: i,n
        real(kind=8), allocatable :: c_star(:),d_star(:)
        n=size(d)
        allocate(c_star(n-1),d_star(n))
        c_star(1)=c/b
        d_star=d(1)/b
        do i=2,n-1
            c_star(i)=c/(b-a*c_star(i-1))
            d_star(i)=(d(i)-a*d_star(i-1))/(b-a*c_star(i-1))
        end do
        d_star(n)=(d(n)-a*d_star(n-1))/(b-a*c_star(i-1))
        x(1)=d_star(n)
        do i=2,n
            x(i)=d_star(n-i+1)-c_star(n-i+1)*x(i-1)
        end do
        x=x(n:1:-1)
    end subroutine Thomas_Toeplitz

    subroutine LU_lin_syst(A,b,tol,x,sing_flag) !> Ax=b
        real(kind=8), intent(in) :: A(:,:) !> square matrix (A=LU)
        real(kind=8), intent(in) :: b(:) !> vector
        real(kind=8), intent(in) :: tol !> tolerance for solution
        real(kind=8), intent(out) :: x(:) !> solution of linear system (must be already allocated)
        logical, intent(out), optional :: sing_flag !> if present, set TRUE on singular matrix instead of error stop
        
        
        real(kind=8), allocatable :: L(:,:), U(:,:), y(:), pb(:)
        real(kind=8) :: det
        integer(kind=4) :: n, i
        integer(kind=4), allocatable :: p(:)
        logical :: error
        
        if (present(sing_flag)) sing_flag = .false.
        
        n=size(b)
        call compute_det(A,tol,det,error)
        
        if (size(A,1)/=n .or. size(A,2)/=n) then
            error stop "Wrong dimensions in LU_lin_syst"
        else if (abs(det)<tol) then
            if (present(sing_flag)) then
                sing_flag = .true.
                return
            end if
            error stop "Zero determinant in LU_lin_syst"
        end if
        
        if (n.eq.2) then
            x(2)=(A(1,1)*b(2)-A(2,1)*b(1))/det
            x(1)=(A(2,2)*b(1)-A(1,2)*b(2))/det
        else
            allocate(L(n,n),U(n,n),y(n),pb(n),p(n))
            call LU(A,tol,L,U,error,perm=p)
            if (error .eqv. .true.) then
                if (present(sing_flag)) then
                    sing_flag = .true.
                    deallocate(L,U,y,pb,p)
                    return
                end if
                error stop "Error in LU decompoisiton"
            end if
            !> Apply row permutation to b
            do i=1,n
                pb(i)=b(p(i))
            end do
            call forward_substitution(L,pb,y)
            call backward_substitution(U,y,x)
            if (inf_norm_vec_real(matmul(A,x)-b) .ge. tol*max(1d0, inf_norm_vec_real(b), maxval(abs(A))*inf_norm_vec_real(x))) then
                if (present(sing_flag)) then
                    sing_flag = .true.
                    deallocate(L,U,y,pb,p)
                    return
                end if
                print *, "DEBUG LU_lin_syst: ||Ax-b||_inf = ", inf_norm_vec_real(matmul(A,x)-b)
                print *, "DEBUG LU_lin_syst: tol          = ", tol
                print *, "DEBUG LU_lin_syst: n            = ", n
                print *, "DEBUG LU_lin_syst: det          = ", det
                print *, "DEBUG LU_lin_syst: ||A||_max    = ", maxval(abs(A))
                print *, "DEBUG LU_lin_syst: ||b||_inf    = ", inf_norm_vec_real(b)
                print *, "DEBUG LU_lin_syst: ||x||_inf    = ", inf_norm_vec_real(x)
                print *, "DEBUG LU_lin_syst: A = "
                do i=1,n
                    print *, "  row", i, ": ", A(i,:)
                end do
                print *, "DEBUG LU_lin_syst: b = ", b
                print *, "DEBUG LU_lin_syst: x = ", x
                print *, "DEBUG LU_lin_syst: Ax-b = ", matmul(A,x)-b
                error stop "Wrong solution in LU_lin_syst"
            end if
        end if
    
    end subroutine

!> Solves linear system Ax=b using Jacobi iterative method
    subroutine Jacobi(A,b,x0,x,niter)
        real(kind=8), intent(in) :: A(:,:)
        real(kind=8), intent(in) :: b(:)
        real(kind=8), intent(inout) :: x0(:)
        real(kind=8), intent(out) :: x(:)
        integer(kind=4), intent(out) :: niter !> number of iterations
        
        real(kind=8), allocatable :: D(:)
        real(kind=8) :: sum
        integer(kind=4) :: i,j,k,n
        real(kind=8), parameter :: tol=1d-12
        n=size(A,1)
        allocate(D(n))
        do i=1,n
            D(i)=A(i,i)
        end do
        niter=0
        do k=1,n
            niter=niter+1 !> we update number of iterations
            do i=1,n
                sum=0
                do j=1,n
                    if (i/=j) then
                        sum=sum+A(i,j)*x0(j)
                    end if
                end do
                x(i)=(1d0/D(i))*(b(i)-sum)
            end do
            if (inf_norm_vec_real(x-x0)<tol) exit
            x0=x
        end do 
    end subroutine Jacobi

    !> Solves linear system Ax=b using Gauss-seidel iterative method
subroutine Gauss_seidel(A,b,x0,x,niter)
    real(kind=8), intent(in) :: A(:,:)
    real(kind=8), intent(in) :: b(:)
    real(kind=8), intent(inout) :: x0(:)
    real(kind=8), intent(out) :: x(:)
    integer(kind=4), intent(out) :: niter !> number of iterations
    
    real(kind=8), allocatable :: D(:)
    real(kind=8) :: sum
    integer(kind=4) :: i,j,k,n
    real(kind=8), parameter :: tol=1d-9
    n=size(A,1)
    allocate(D(n))
    do i=1,n
        D(i)=A(i,i)
    end do
    do k=1,n
        niter=niter+1 !> we update number of iterations
        do i=1,n
            sum=0
            do j=1,n
                if (j<i .and. i>1) then
                    sum=sum+A(i,j)*x(j)
                else if (j>i) then
                    sum=sum+A(i,j)*x0(j)
                else
                    continue
                end if
            end do
            x(i)=(1d0/D(i))*(b(i)-sum)
        end do
        if (inf_norm_vec_real(x-x0)<tol) exit
        x0=x
    end do
end subroutine Gauss_seidel

!> Gauss-Jordan elimination
!> Ax=b
subroutine Gauss_Jordan(A,b,tol,x,error)
    real(kind=8), intent(in) :: A(:,:) !> square matrix
    real(kind=8), intent(in) :: b(:) !> vector
    real(kind=8), intent(in) :: tol !> tolerance for solution
    real(kind=8), intent(out) :: x(:) !> solution of linear system
    integer(kind=4), intent(out), optional :: error !> Error flag: 0->No error, 1->Singular equations
    
    real(kind=8), parameter :: epsilon=1d-16
    real(kind=8) :: aux,factor
    real(kind=8), allocatable :: A_aux(:,:)
    integer(kind=4) :: n,icol,jcol,irow,jrow,ipeak
    n=size(A,1)
    if (n/=size(b)) then
        error stop "Dimension error in Gauss-Jordan"
    end if
    A_aux=A
    x=b
    !> Process n times to get all equations
    main_loop: do icol = 1,n
        !> Find peak pivot for column icol
        irow=icol
        ipeak = 1
        max_pivot: do jrow = 2,n
            if (abs(A_aux(jrow,icol))>abs(A_aux(ipeak,icol))) then
                ipeak = jrow
            end if
        end do max_pivot
        !> Check for singular equations
        if (abs(A_aux(ipeak,icol))<epsilon) then
            if (present(error)) then
                error=1
            end if
            return
        end if
        !> Otherwise, if ipeak /= irow, swap equations irow & ipeak
        if (ipeak/=irow .and. icol<n) then
            do jcol = 1,n
                aux=A_aux(ipeak,jcol)
                A_aux(ipeak,jcol)=A_aux(irow,jcol)
                A_aux(irow,jcol)=aux
            end do
            aux=x(ipeak)
            x(ipeak)=x(irow)
            x(irow)=aux
        end if
        !> Multiply equation irow by -A_aux(jrow,icol)/A_aux(irow,icol) and add it to eqn jrow
        do jrow = 1,n
            if (jrow/=irow) then
                factor=A_aux(jrow,icol)/A_aux(irow,icol)
                do jcol = 1,n
                    A_aux(jrow,jcol)=A_aux(jrow,jcol)-A_aux(irow,jcol)*factor
                end do
                x(jrow)=x(jrow)-x(irow)*factor
            end if
        end do
    end do main_loop
    !> End of main loop over all equations. All off-diagonal terms are now zero.
    !> To get the final answer, we must divide each equation by the coefficient of its diagonal term
    do irow = 1,n
        x(irow)=x(irow)/A_aux(irow,irow)
        A_aux(irow,irow)=1d0
    end do
    if (present(error)) then
        error=0
    end if
    if (inf_norm_vec_real(matmul(A,x)-b)>=tol) then
        print *, inf_norm_vec_real(matmul(A,x)-b)
        error stop "Wrong solution in Gauss-Jordan"
    end if
end subroutine

subroutine compute_Thomas_coeffs(A,b,tol,c_mix,d_tilde)
    !> Computes coefficients of Thomas algorithm
    
    !> A: tridiagonal matrix
    !> c_mix,d_tilde: coefficients

    
    class(tridiag_matrix_c), intent(in) :: A
    real(kind=8), intent(in) :: b(:)
    real(kind=8), intent(in) :: tol
    real(kind=8), intent(out) :: c_mix(:),d_tilde(:) !> tiene que estar alocatado
    
    integer(kind=4) :: i,n
    real(kind=8) :: denom
    
    n=size(b)
    c_mix(1)=A%super(1)/A%diag(1)
    d_tilde(1)=b(1)/A%diag(1)
    do i=2,n-1
        denom=A%diag(i)-A%sub(i-1)*c_mix(i-1)
        if (abs(denom)<tol) then
            error stop "Singularity in Thomas algorithm"
        end if
        c_mix(i)=A%super(i)/denom
        d_tilde(i)=(b(i)-A%sub(i-1)*d_tilde(i-1))/denom
    end do
    denom=A%diag(n)-A%sub(n-1)*c_mix(n-1)
    if (abs(denom)<tol) then
        error stop "Singularity in Thomas algorithm"
    end if
    d_tilde(n)=(b(n)-A%sub(n-1)*d_tilde(n-1))/denom
end subroutine

subroutine compute_inverse_tridiag_matrix(this,tol,inv_mat)
    implicit none
    class(tridiag_matrix_c), intent(in) :: this
    real(kind=8), intent(in) :: tol
    real(kind=8), intent(out) :: inv_mat(:,:) !> must be allocated
    
    integer(kind=4) :: i,j
    real(kind=8), allocatable :: id_col(:)
    
    allocate(id_col(this%num_cols))
    
    do j=1,this%num_cols
        id_col=0d0
        id_col(j)=1d0
        call Thomas(this,id_col,tol,inv_mat(:,j))
    end do
end subroutine


!> Solves Ax=b using LAPACK dgesv (pivoted LU) with row/column equilibration.
!> Scaling improves conditioning when entries of A span many orders of magnitude.
    subroutine scaled_dgesv_lin_syst(A, b, x)
        implicit none
        real(kind=8), intent(in) :: A(:,:)
        real(kind=8), intent(in) :: b(:)
        real(kind=8), intent(out) :: x(:)
        
        integer(kind=4) :: n, info
        integer(kind=4), allocatable :: ipiv(:)
        real(kind=8), allocatable :: A_work(:,:), b_work(:)
        real(kind=8), allocatable :: row_scale(:), col_scale(:)
        real(kind=8) :: row_max, col_max
        integer(kind=4) :: i, j
        
        n = size(b)
        allocate(A_work(n,n), b_work(n), ipiv(n))
        allocate(row_scale(n), col_scale(n))
        
        !> Row equilibration: scale each row so max|A(i,:)| = 1
        do i = 1, n
            row_max = maxval(abs(A(i,:)))
            if (row_max > 0d0) then
                row_scale(i) = 1d0 / row_max
            else
                row_scale(i) = 1d0
            end if
        end do
        
        do i = 1, n
            A_work(i,:) = row_scale(i) * A(i,:)
            b_work(i) = row_scale(i) * b(i)
        end do
        
        !> Column equilibration: scale each column so max|A(:,j)| = 1
        do j = 1, n
            col_max = maxval(abs(A_work(:,j)))
            if (col_max > 0d0) then
                col_scale(j) = 1d0 / col_max
            else
                col_scale(j) = 1d0
            end if
            A_work(:,j) = col_scale(j) * A_work(:,j)
        end do
        
        !> Solve scaled system with LAPACK dgesv (partial pivoting)
        call dgesv(n, 1, A_work, n, ipiv, b_work, n, info)
        
        if (info /= 0) then
            print *, "WARNING scaled_dgesv_lin_syst: dgesv info = ", info
        end if
        
        !> Unscale solution: x_original = col_scale * x_scaled
        x = col_scale * b_work
        
        !> Residual check: ||Ax - b|| / max(1, ||b||, ||A||*||x||)
        if (inf_norm_vec_real(matmul(A,x)-b) > 1d-8 * max(1d0, inf_norm_vec_real(b), maxval(abs(A))*inf_norm_vec_real(x))) then
            print *, "WARNING scaled_dgesv: ||Ax-b||=", inf_norm_vec_real(matmul(A,x)-b), &
                " ||b||=", inf_norm_vec_real(b), " info=", info
        end if
        
        deallocate(A_work, b_work, ipiv, row_scale, col_scale)
    end subroutine


end module 