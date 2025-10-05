!> This module contains methods to solve linear systems of equations
module metodos_sist_lin_m
    use matrices_m
    implicit none
    save
    interface
        ! subroutine Gauss_Jordan(A,b,tol,x,error)
        !     implicit none
        !     real(kind=8), intent(in) :: A(:,:) 
        !     real(kind=8), intent(in) :: b(:)
        !     real(kind=8), intent(in) :: tol !> tolerance for solution
        !     real(kind=8), intent(out) :: x(:)
        !     integer(kind=4), intent(out), optional :: error
        ! end subroutine
        
        ! subroutine forward_substitution(L,b,x)
        !     implicit none
        !     real(kind=8), intent(in) :: L(:,:)
        !     real(kind=8), intent(in) :: b(:)
        !     real(kind=8), intent(out) :: x(:)
        ! end subroutine
        
        ! subroutine backward_substitution(U,b,x)
        !     implicit none
        !     real(kind=8), intent(in) :: U(:,:)
        !     real(kind=8), intent(in) :: b(:)
        !     real(kind=8), intent(out) :: x(:)
        ! end subroutine
        
        ! subroutine LU_lin_syst(A,b,tol,x)
        !     implicit none
        !     real(kind=8), intent(in) :: A(:,:)
        !     real(kind=8), intent(in) :: b(:)
        !     real(kind=8), intent(in) :: tol
        !     real(kind=8), intent(out) :: x(:)
        ! end subroutine
        
        
        ! subroutine Thomas_Toeplitz(a,b,c,d,x)
        !     implicit none
        !     real(kind=8), intent(in) :: a,b,c
        !     real(kind=8), intent(in) :: d(:)
        !     real(kind=8), intent(out) :: x(:)
        ! end subroutine
        
        ! subroutine Jacobi(A,b,x0,x,niter)
        !     implicit none
        !     real(kind=8), intent(in) :: A(:,:)
        !     real(kind=8), intent(in) :: b(:)
        !     real(kind=8), intent(inout) :: x0(:)
        !     real(kind=8), intent(out) :: x(:)
        !     integer(kind=4), intent(out) :: niter !> number of iterations
        ! end subroutine
        
        ! subroutine Gauss_seidel(A,b,x0,x,niter)
        !     implicit none
        !     real(kind=8), intent(in) :: A(:,:)
        !     real(kind=8), intent(in) :: b(:)
        !     real(kind=8), intent(inout) :: x0(:)
        !     real(kind=8), intent(out) :: x(:)
        !     integer(kind=4), intent(out) :: niter !> number of iterations
        ! end subroutine
    end interface

    contains
    !> Inverse of square matrix using LU decomposition
subroutine inv_matrix(A,tol,inv)
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

!> Inverse of square matrix using LU decomposition from first element
subroutine inv_matrix_WMA(A,tol,inv)
    implicit none
    real(kind=8), intent(in) :: A(:,:)
    real(kind=8), intent(in) :: tol
    real(kind=8), intent(out) :: inv(:,:)
    
    integer(kind=4) :: n,j,i
    real(kind=8), parameter :: epsilon=1d-6
    real(kind=8), allocatable :: id(:,:), prod_A_invA(:,:), id_col(:), inv_col(:)
    logical :: nzdiag
    
    if (size(A,1)/=size(A,2)) error stop "Matrix must be square (inv_matrix_WMA)"
    
    !if (det(A)>=tol) error stop "Matrix is not invertible"
    
    call inv_matrix(A,tol,inv) 
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
    
    subroutine Thomas_induction(A,b,i,x_i,tol,x_k)  !> to be developed
        !> A: tridiagonal matrix
        !> b: independent term
        !> x: solution of linear system
    
        class(tridiag_matrix_c), intent(in) :: A
        real(kind=8), intent(in) :: b(:)
        integer(kind=4), intent(in) :: i !> row index
        real(kind=8), intent(in) :: x_i
        real(kind=8), intent(in) :: tol !> tolerance
        real(kind=8), intent(out) :: x_k(:) !> must be allocated
        
        integer(kind=4) :: k
        real(kind=8), allocatable :: c_star(:),d_star(:), matrix(:,:)
        
        !allocate(c_star(this%dim-1),d_star(this%dim))
        !c_star(1)=A%super(1)/A%diag(1)
        !d_star(1)=b(1)/A%diag(1)
        !do i=2,this%dim-1
        !>    c_star(i)=A%super(i)/(A%diag(i)-A%sub(i-1)*c_star(i-1))
        !>    d_star(i)=(b(i)-A%sub(i-1)*d_star(i-1))/(A%diag(i)-A%sub(i-1)*c_star(i-1))
        !end do
        !d_star(this%dim)=(b(this%dim)-A%sub(i-1)*d_star(this%dim-1))/(A%diag(i)-A%sub(i-1)*c_star(i-1))
        !x(this%dim)=d_star(this%dim)
        !x(this%dim-1)=d_star(this%dim-1)-c_star(this%dim-1)*x(this%dim)
        !do k=2,this%dim-1
        !>    sum=0d0
        !>    do j=0,k-1
        !>        prod=1d0
        !>        do l=1,k-j
        !>            prod=prod*c_star(i-j-l)
        !>        end do
        !>        sum=sum+prod*d_star(i-j)
        !>    x(this%dim-k)=d_star(this%dim-k)+sum()
        !end do
        !deallocate(c_star,d_star)
    end subroutine

    subroutine LU_lin_syst(A,b,tol,x) !> Ax=b
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
                error stop "Error in LU decompoisiton"
            end if
            call forward_substitution(L,b,y)
            call backward_substitution(U,y,x)
            if (inf_norm_vec_real(matmul(A,x)-b) .ge. tol) then
                print *, inf_norm_vec_real(matmul(A,x)-b)
                error stop "Wrong solution in LU_lin_syst"
                !print *, "Wrong solution in LU_lin_syst"
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
        
        real(kind=8), allocatable :: L(:,:), R(:,:), D(:), C_mat(:,:), c(:)
        real(kind=8) :: sum
        integer(kind=4) :: i,j,k,n
        real(kind=8), parameter :: tol=1d-12
        n=size(A,1)
        allocate(L(n,n),R(n,n),D(n),C_mat(n,n),c(n))
        forall (i=1:n)
            D(i)=A(i,i)
        end forall
        L=0d0
        R=0d0
        do j=1,n
            L((j+1):n,j)=-A((j+1):n,j)
            R(j,(j+1):n)=-A(j,(j+1):n)
        end do
        C_mat=L+R
        do i=1,n
            C_mat(i,1:n)=(1d0/D(i))*C_mat(i,1:n)
        end do
        forall (i=1:n)
            c(i)=(1d0/D(i))*b(i)
        end forall
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
    
    real(kind=8), allocatable :: L(:,:), R(:,:), D(:), C_mat(:,:), inv(:,:), c(:)
    real(kind=8) :: sum
    integer(kind=4) :: i,j,k,n
    real(kind=8), parameter :: tol=1d-9
    n=size(A,1)
    allocate(L(n,n),R(n,n),D(n),C_mat(n,n),inv(n,n),c(n))
    forall (i=1:n)
        D(i)=A(i,i)
    end forall
    L=0d0
    R=0d0
    do j=1,n
        L((j+1):n,j)=-A((j+1):n,j)
        R(j,(j+1):n)=-A(j,(j+1):n)
    end do
    forall (i=1:n)
        L(i,i)=L(i,i)-D(i)
    end forall
    call inv_matrix(-L,tol,inv)
    C_mat=matmul(inv,R)
    c=matmul(inv,b)
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

subroutine compute_Thomas_coeffs(A,b,tol,c_tilde,d_tilde)
    !> Computes coefficients of Thomas algorithm
    
    !> A: tridiagonal matrix
    !> c_tilde,d_tilde: coefficients

    
    class(tridiag_matrix_c), intent(in) :: A
    real(kind=8), intent(in) :: b(:)
    real(kind=8), intent(in) :: tol
    real(kind=8), intent(out) :: c_tilde(:),d_tilde(:) !> tiene que estar alocatado
    
    integer(kind=4) :: i,n
    real(kind=8) :: denom
    
    c_tilde(1)=A%super(1)/A%diag(1)
    d_tilde(1)=b(1)/A%diag(1)
    do i=2,n-1
        denom=A%diag(i)-A%sub(i-1)*c_tilde(i-1)
        if (abs(denom)<tol) then
            error stop "Singularity in Thomas algorithm"
        end if
        c_tilde(i)=A%super(i)/denom
        d_tilde(i)=(b(i)-A%sub(i-1)*d_tilde(i-1))/denom
    end do
    denom=A%diag(n)-A%sub(n-1)*c_tilde(n-1)
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



end module 