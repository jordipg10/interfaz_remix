module matrices_m
    use vectors_m
    implicit none
    save
    
    type, public, abstract :: array_c !> array superclass
        integer(kind=4) :: num_cols
    contains
        procedure, public :: allocate_array
        procedure, public :: allocate_columns
        procedure, public :: prod_mat_vec
        !procedure, public :: prod_mat_mat
        procedure, public :: get_diag
        procedure, public :: get_sub
        procedure, public :: get_super
        procedure, public :: compute_norm_inf
        procedure, public :: compute_norm_1
    end type
    
    type, public, extends(array_c) :: non_sq_matrix_c
    end type
    
    type, public, extends(non_sq_matrix_c) :: real_array_c
        type(vector_real_c), allocatable :: cols(:) !> columns of the array
    contains
        procedure, public :: is_matrix_full
        procedure, public :: get_can_vec !> gets information of canonical vectors
    end type
    
    type, public, extends(non_sq_matrix_c) :: int_array_c
        type(vector_int_c), allocatable :: cols(:)
    contains
        procedure, public :: get_vector_int
    end type
    
    type, public, extends(array_c) :: sq_matrix_c
        !integer(kind=4) :: dim
        real(kind=8), allocatable :: eigenvalues(:)
        real(kind=8), allocatable :: eigenvectors(:,:)
    contains
        procedure, public :: compute_eigenvalues
        procedure, public :: compute_eigenvectors
    end type
!****************************************************************************************************************************************************
    type, public, extends(sq_matrix_c) :: tridiag_sym_Toeplitz_matrix_c
        real(kind=8) :: sub
        real(kind=8) :: diag
    end type
    
    type, public, extends(tridiag_sym_Toeplitz_matrix_c) :: tridiag_Toeplitz_matrix_c
        real(kind=8) :: super
    contains
        procedure, public :: set_tridiag_Toeplitz_matrix
    end type
!****************************************************************************************************************************************************
    type, public, extends(sq_matrix_c) :: diag_matrix_c
        real(kind=8), allocatable :: diag(:)
    contains
        procedure, public :: set_diag_matrix
        procedure, public :: prod_mat_diag_mat
    end type
    
    type, public, extends(diag_matrix_c) :: tridiag_sym_matrix_c
        real(kind=8), allocatable :: sub(:)
    contains
        procedure, public :: check_eigenvectors_tridiag_sym_matrix
    end type

    type, public, extends(tridiag_sym_matrix_c) :: tridiag_matrix_c
        real(kind=8), allocatable :: super(:)
    contains
        procedure, public :: set_tridiag_matrix
        procedure, public :: compute_transpose_tridiag_matrix
        ! procedure, public :: compute_inverse_tridiag_matrix
        procedure, public :: prod_tridiag_mat_mat
    end type

    !type, public, extends(tridiag_matrix_c) :: tridiag_matrix_vec_c
    !>    real(kind=8), allocatable :: vector(:)
    !end type
!****************************************************************************************************************************************************
    interface
        ! subroutine compute_eigenvalues(this)
        !     import sq_matrix_c
        !     implicit none
        !     class(sq_matrix_c) :: this
        ! end subroutine
        
        ! subroutine compute_eigenvectors(this)
        !     import sq_matrix_c
        !     implicit none
        !     class(sq_matrix_c) :: this
        ! end subroutine
        

        
        
        
        ! function id_matrix(n)
        ! !> Identity nxn matrix
        !     implicit none
        !     integer(kind=4), intent(in) :: n
        !     integer(kind=4) :: i
        !     real(kind=8) :: id_matrix(n,n)
        ! end function
        
        ! function norm_mat_inf(A) result(norm)
        !     implicit none
        !     real(kind=8), intent(in) :: A(:,:)
        !     real(kind=8) :: norm
        ! end function
        
        ! function norm_mat_1(A) result(norm)
        !     implicit none
        !     real(kind=8), intent(in) :: A(:,:)
        !     real(kind=8) :: norm
        ! end function
        
        ! subroutine compute_det(A,tol,det,error)
        ! !> Determinant of square matrix using LU decomposition
        !     implicit none
        !     real(kind=8), intent(in) :: A(:,:)
        !     real(kind=8), intent(in) :: tol
        !     real(kind=8), intent(out) :: det
        !     logical, intent(out) :: error
        ! end subroutine
        
        ! subroutine inv_matrix(A,tol,inv)
        ! !> Inverse of square matrix using Gauss-Jordan
        !     implicit none
        !     real(kind=8), intent(in) :: A(:,:)
        !     real(kind=8), intent(in) :: tol
        !     real(kind=8), intent(out) :: inv(:,:)
        ! end subroutine inv_matrix
        
        ! function prod_mat_vec(this,b) result(x) !> Ab=x
        !     import array_c
        !     implicit none
        !     class(array_c), intent(in) :: this
        !     real(kind=8), intent(in) :: b(:)
        !     real(kind=8), allocatable :: x(:)
        ! end function
        
        ! function prod_tridiag_mat_mat(this,X_mat) result(C_mat) !> AB=C
        !     import tridiag_matrix_c
        !     implicit none
        !     class(tridiag_matrix_c), intent(in) :: this
        !     real(kind=8), intent(in) :: X_mat(:,:)
        !     real(kind=8), allocatable :: C_mat(:,:)
        ! end function 
        
        ! subroutine LU(A,tol,L,U,error)
        !     implicit none
        !     real(kind=8), intent(in) :: A(:,:)
        !     real(kind=8), intent(in) :: tol
        !     real(kind=8), intent(out) :: L(:,:), U(:,:)
        !     logical, intent(out) :: error
        ! end subroutine LU
        
        ! subroutine Cholesky(A,L)
        !     implicit none
        !     real(kind=8), intent(in) :: A(:,:)
        !     real(kind=8), intent(out) :: L(:,:)
        ! end subroutine Cholesky
        
        ! subroutine QR_Householder(A,Q,R)
        !     implicit none
        !     real(kind=8), intent(in) :: A(:,:) !> square matrix
        !     real(kind=8), intent(out) :: Q(:,:), R(:,:)
        ! end subroutine QR_Householder
        
        ! function Householder(x)
        !     implicit none
        !     real(kind=8), intent(in) :: x(:)
        !     real(kind=8), allocatable :: Householder(:,:)
        ! end function Householder
        
        ! subroutine Gram_Schmidt_mat(v,u)
        !     implicit none
        !     real(kind=8), intent(in) :: v(:,:)
        !     real(kind=8), intent(out) :: u(:,:)
        ! end subroutine
        
        ! subroutine non_zeros(A,n,ind)
        !     implicit none
        !     real(kind=8), intent(in) :: A(:,:)
        !     integer(kind=4), intent(out) :: n
        !     integer(kind=4), intent(out) :: ind(:,:)
        ! end subroutine
        
        ! function prod_AT_A(A) result(AT_A)
        !     implicit none
        !     real(kind=8), intent(in) :: A(:,:)
        !     real(kind=8), allocatable :: AT_A(:,:)
        ! end function
        
        ! subroutine check_eigenvectors_tridiag_sym_matrix(this,tolerance)
        !     import tridiag_sym_matrix_c
        !     implicit none
        !     class(tridiag_sym_matrix_c), intent(in) :: this
        !     real(kind=8), intent(in) :: tolerance
        ! end subroutine
        
        ! function prod_diag_tridiag_mat(A,B) result(C)
        !     import diag_matrix_c
        !     import tridiag_matrix_c
        !     class(diag_matrix_c), intent(in) :: A
        !     class(tridiag_matrix_c), intent(in) :: B
        !     type(tridiag_matrix_c) :: C
        ! end function
        
        ! function prod_tridiag_diag_mat(A,B) result(C)
        !     import diag_matrix_c
        !     import tridiag_matrix_c
        !     class(tridiag_matrix_c), intent(in) :: A
        !     class(diag_matrix_c), intent(in) :: B
        !     type(tridiag_matrix_c) :: C
        ! end function
        
        ! subroutine compute_inverse_tridiag_matrix(this,tol,inv_mat)
        !     import tridiag_matrix_c
        !     implicit none
        !     class(tridiag_matrix_c), intent(in) :: this
        !     real(kind=8), intent(in) :: tol
        !     real(kind=8), intent(out) :: inv_mat(:,:)
        ! end subroutine
    end interface
!****************************************************************************************************************************************************
    contains
        subroutine allocate_array(this,n)
            implicit none
            class(array_c) :: this
            integer(kind=4), intent(in), optional :: n
            if (present(n)) then
                this%num_cols=n
            end if
            select type (this)
            class is (diag_matrix_c)
                allocate(this%diag(this%num_cols))
                select type (this)
                class is (tridiag_sym_matrix_c)
                    allocate(this%sub(this%num_cols-1))
                    select type (this)
                    type is (tridiag_matrix_c)
                        allocate(this%super(this%num_cols-1))
                        !select type (this)
                        !type is (tridiag_matrix_vec_c)
                        !>    allocate(this%vector(n))
                        !end select
                    end select
                end select
            type is (real_array_c)
                allocate(this%cols(this%num_cols))
                !do i=1,this%num_cols
                !    call this%cols(i)%allocate_vector()
            type is (int_array_c)
                if (allocated(this%cols)) then
                    deallocate(this%cols)
                end if
                allocate(this%cols(this%num_cols))
            end select
        end subroutine
        
        subroutine allocate_columns(this)
            implicit none
            class(array_c) :: this
            integer(kind=4) :: i
            select type (this)
            type is (real_array_c)
                do i=1,this%num_cols
                    call this%cols(i)%allocate_vector()
                end do
            type is (int_array_c)
                do i=1,this%num_cols
                    call this%cols(i)%allocate_vector()
                end do
            end select
        end subroutine
        
        subroutine set_diag_matrix(this,diag)
            implicit none
            class(diag_matrix_c) :: this
            real(kind=8), intent(in) :: diag(:)
            this%diag=diag
            this%num_cols=size(diag)
        end subroutine
        
        subroutine set_tridiag_Toeplitz_matrix(this,sub,diag,super)
            implicit none
            class(tridiag_Toeplitz_matrix_c) :: this
            real(kind=8), intent(in) :: sub,diag,super
            this%sub=sub
            this%diag=diag
            this%super=super
        end subroutine
        
        subroutine set_tridiag_matrix(this,sub,diag,super)
            implicit none
            class(tridiag_matrix_c) :: this
            real(kind=8), intent(in) :: sub(:),diag(:),super(:)
            if (size(sub)/=size(super) .or. size(diag)/=(size(sub)+1)) error stop "Dimension error in set_tridiag_matrix" 
            this%sub=sub
            this%diag=diag
            this%super=super
        end subroutine
        
        function get_diag(this) result(diag)
            implicit none
            class(array_c), intent(in) :: this
            real(kind=8), allocatable :: diag(:)
            select type (this)
            class is (diag_matrix_c)
                diag=this%diag
            end select
        end function
        
        function get_sub(this) result(sub)
            implicit none
            class(array_c), intent(in) :: this
            real(kind=8), allocatable :: sub(:)
            select type (this)
            class is (tridiag_sym_matrix_c)
                sub=this%sub
            end select
        end function
        
        function get_super(this) result(super)
            implicit none
            class(array_c), intent(in) :: this
            real(kind=8), allocatable :: super(:)
            select type (this)
            class is (tridiag_matrix_c)
                super=this%super
            end select
        end function
        
        function compute_norm_inf(this) result(norm)
            implicit none
            class(array_c), intent(in) :: this
            real(kind=8) :: norm
            
            integer(kind=4) :: i,n
            real(kind=8) :: norm_i
            
            select type (this)
            type is (diag_matrix_c)
                n=size(this%diag)
                norm=abs(this%diag(1))
                do i=2,n
                    norm_i=abs(this%diag(1))
                    if (norm_i>norm) then
                        norm=norm_i
                    end if
                end do
            type is (tridiag_matrix_c)
                n=size(this%diag)
                norm=abs(this%diag(1))+abs(this%super(1))
                do i=2,n-1
                    norm_i=abs(this%sub(i-1))+abs(this%diag(i))+abs(this%super(i))
                    if (norm_i>norm) then
                        norm=norm_i
                    end if
                end do
                norm_i=abs(this%sub(n-1))+abs(this%diag(n))
                if (norm_i>norm) then
                    norm=norm_i
                end if
            end select
        end function
        
        function compute_norm_1(this) result(norm)
            implicit none
            class(array_c), intent(in) :: this
            real(kind=8) :: norm
            
            type(tridiag_matrix_c) :: transpose
            
            select type (this)
            type is (tridiag_matrix_c)
                call this%compute_transpose_tridiag_matrix(transpose)
                norm=transpose%compute_norm_inf()
            end select
        end function
        
        subroutine compute_transpose_tridiag_matrix(this,transpose)
            implicit none
            class(tridiag_matrix_c), intent(in) :: this
            class(tridiag_matrix_c), intent(out) :: transpose
            !> habria que alocatar la traspuesta
            transpose%sub=this%super
            transpose%diag=this%diag
            transpose%super=this%sub
        end subroutine
        
        subroutine get_submatrix(mat,rows,cols,submat)
            implicit none
            real(kind=8), intent(in) :: mat(:,:)
            integer(kind=4), intent(in) :: rows(:),cols(:)
            real(kind=8), intent(out) :: submat(:,:) !> must be allocated
            
            integer(kind=4) :: i,j,n_r,n_c,i_r,i_c
            n_r=size(rows)
            n_c=size(cols)
            do i=1,n_r
                do j=1,n_c
                    submat(i,j)=mat(rows(i),cols(j))
                end do
            end do
        end subroutine
        
        function prod_mat_diag_mat(this,mat) result(prod)
            implicit none
            class(diag_matrix_c), intent(in) :: this
            real(kind=8), intent(in) :: mat(:,:)
            real(kind=8), allocatable :: prod(:,:)
            integer(kind=4) :: j
            if (size(mat,2)/=this%num_cols) then
                error stop "Dimension error in prod_mat_diag_mat"
            end if
            allocate(prod(size(mat,1),this%num_cols))
            do j=1,this%num_cols
                prod(:,j)=mat(:,j)*this%diag(j)
            end do
        end function
        
        function get_vector_int(this) result(vector_int)
            implicit none
            class(int_array_c), intent(in) :: this
            integer(kind=4), allocatable :: vector_int(:)
            
            integer(kind=4) :: i,j,dim_vec
            
            dim_vec=0
            do i=1,this%num_cols
                dim_vec=dim_vec+this%cols(i)%dim
            end do
            allocate(vector_int(dim_vec))
            
            j=1
            do i=1,this%num_cols
                vector_int(j:j+this%cols(i)%dim-1)=this%cols(i)%col_1
                j=j+this%cols(i)%dim
            end do
        end function
        
        subroutine is_matrix_full(this,flag)
            implicit none
            class(real_array_c), intent(in) :: this
            logical, intent(out) :: flag
            
            integer(kind=4) :: i
            
            flag=.true.
            do i=1,this%num_cols
                if (this%cols(i)%dim<this%num_cols) then
                    flag=.false.
                    exit
                end if
            end do
        end subroutine

        !> Spectral radius
function radio_espectral(lambda) result(rho)
    implicit none
    real(kind=8), intent(in) :: lambda(:)
    real(kind=8) :: rho
    rho=max(maxval(lambda),maxval(-lambda))
end function 

function prod_tridiag_mat_mat(this,X_mat) result(C_mat) !> AB=C
    implicit none
    class(tridiag_matrix_c), intent(in) :: this !> nxn
    real(kind=8), intent(in) :: X_mat(:,:) !> nxm
    real(kind=8), allocatable :: C_mat(:,:) !> nxm
    
    integer(kind=4) :: i,j,n,m
    
    m=size(X_mat,2)
    !select type (this)
    !type is (tridiag_matrix_c)
        n=size(this%diag)
        if (n/=size(X_mat,1)) error stop "Dimension error in prod_mat_mat"
        if (.not. allocated(C_mat)) then
            allocate(C_mat(n,m))
        end if
        C_mat(1,:)=this%diag(1)*X_mat(1,:)+this%super(1)*X_mat(2,:)
        do i=2,n-1
            do j=1,m
                C_mat(i,j)=this%sub(i-1)*X_mat(i-1,j)+this%diag(i)*X_mat(i,j)+this%super(i)*X_mat(i+1,j)
            end do
        end do
        C_mat(n,:)=this%sub(n-1)*X_mat(n-1,:)+this%diag(n)*X_mat(n,:)
    !end select
end function 

function prod_tridiag_diag_mat(A,B) result(C) !> A*B=C
    class(tridiag_matrix_c), intent(in) :: A !> tridiagonal matrix
    class(diag_matrix_c), intent(in) :: B !> diagonal matrix
    type(tridiag_matrix_c) :: C

    integer(kind=4) :: i,n
    
    n=size(B%diag)
    
    call C%allocate_array(n)
        
    C%diag(1)=A%diag(1)*B%diag(1)
    
    if (n>1) then
        C%super(1)=A%super(1)*B%diag(2)
        C%sub(n-1)=A%sub(n-1)*B%diag(n-1)
        C%diag(n)=A%diag(n)*B%diag(n)
    end if
    if (n>2) then
        do i=2,n-1
            C%diag(i)=A%diag(i)*B%diag(i)
            C%super(i)=A%super(i)*B%diag(i+1)
            C%sub(i-1)=A%sub(i-1)*B%diag(i-1)
        end do
    end if
end function

function prod_total_sym_mat(A,y0,b,time) result(y)
    !> dy/dt=-Ay+b
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

    function prod_mat_vec(A,b) result(x) !> Ab=x
        implicit none
        class(array_c), intent(in) :: A
        real(kind=8), intent(in) :: b(:)
        real(kind=8), allocatable :: x(:) 
        
        integer(kind=4) :: i,j,n,m
        real(kind=8) :: sum
        
        select type (A)
        !type is (tridiag_matrix_vec_c)
        !>    n=size(A%diag)
        !>    if (size(b)/=2*n) error stop "Dimension error in b"
        !>    if (n.eqv.1) then
        !>        x=A%diag(1)*b(1)+A%vector(1)*b(n+1)
        !>    else if (n.eqv.2) then
        !>        x(1)=A%diag(1)*b(1)+A%super(1)*b(2)+A%vector(1)*b(n+1)
        !>        x(2)=A%sub(1)*b(1)+A%diag(2)*b(2)+A%vector(2)*b(n+2)
        !>    else
        !>        x(1)=A%diag(1)*b(1)+A%super(1)*b(2)+A%vector(1)*b(n+1)
        !>        do i=2,n-1
        !>            x(i)=A%sub(i-1)*b(i-1)+A%diag(i)*b(i)+A%super(i)*b(i+1)
        !>            x(i)=x(i)+A%vector(i)*b(n+i)
        !>        end do
        !>        x(n)=A%sub(n-1)*b(n-1)+A%diag(n)*b(n)+A%vector(n)*b(2*n)
        !>    end if
        type is (tridiag_matrix_c)
            n=size(A%diag)
            allocate(x(n))
            if (size(b)/=n) error stop "Dimension error in b"
            x(1)=A%diag(1)*b(1)+A%super(1)*b(2)
            do i=2,n-1
                x(i)=A%sub(i-1)*b(i-1)+A%diag(i)*b(i)+A%super(i)*b(i+1)
            end do
            x(n)=A%sub(n-1)*b(n-1)+A%diag(n)*b(n)
        end select
    end function 

    function prod_mat_mat(this,X_mat) result(C_mat) !> AB=C
        implicit none
        class(array_c), intent(in) :: this !> nxn
        real(kind=8), intent(in) :: X_mat(:,:) !> nxm
        real(kind=8), allocatable :: C_mat(:,:) !> nxm
        
        integer(kind=4) :: i,j,n,m
        
        m=size(X_mat,2)
        select type (this)
        type is (tridiag_matrix_c)
            n=size(this%diag)
            if (n/=size(X_mat,1)) error stop "Dimension error in prod_mat_mat"
            if (.not. allocated(C_mat)) then
                allocate(C_mat(n,m))
            end if
            C_mat(1,:)=this%diag(1)*X_mat(1,:)+this%super(1)*X_mat(2,:)
            do i=2,n-1
                do j=1,m
                    C_mat(i,j)=this%sub(i-1)*X_mat(i-1,j)+this%diag(i)*X_mat(i,j)+this%super(i)*X_mat(i+1,j)
                end do
            end do
            C_mat(n,:)=this%sub(n-1)*X_mat(n-1,:)+this%diag(n)*X_mat(n,:)
        end select
    end function 

    function prod_diag_tridiag_mat(A,B) result(C) !> A*B=C
        class(diag_matrix_c), intent(in) :: A !> diagonal matrix
        class(tridiag_matrix_c), intent(in) :: B !> tridiagonal matrix  
        type(tridiag_matrix_c) :: C

        integer(kind=4) :: i,n
        
        n=size(A%diag)
        
        call C%allocate_array(n)
        
        C%diag(1)=A%diag(1)*B%diag(1)
        C%super(1)=A%diag(1)*B%super(1)
        
        C%sub(n-1)=A%diag(n)*B%sub(n-1)
        C%diag(n)=A%diag(n)*B%diag(n)
        if (n>2) then
            do i=2,n-1
                C%diag(i)=A%diag(i)*B%diag(i)
                C%super(i)=A%diag(i)*B%super(i)
                C%sub(i-1)=A%diag(i)*B%sub(i-1)
            end do
        end if
        
    end function

    function prod_AT_A(A) result(AT_A) !> AT_A=A^T * A
        implicit none
        real(kind=8), intent(in) :: A(:,:)
        real(kind=8), allocatable :: AT_A(:,:)
        
        integer(kind=4) :: n,m,i,j
        
        n=size(A,1)
        m=size(A,2)
        
        allocate(AT_A(m,m))
    
        do i=1,m
            do j=1,m
                AT_A(i,j)=dot_product(A(:,i),A(:,j))
            end do
        end do
    end function

    subroutine potencia(A,z0,tolerance,rho,niter)
        !> Computes spectral radius
            implicit none
            real(kind=8), intent(in) :: A(:,:)
            real(kind=8), intent(inout) :: z0(:)
            real(kind=8), intent(in) :: tolerance
            real(kind=8), intent(out) :: rho !> spectral radius
            integer(kind=4), intent(out) :: niter !> number of iterations
            
            integer(kind=4) :: n,k
            integer(kind=4), parameter :: iter_max=100
            real(kind=8) :: rho_old,rho_new
            real(kind=8), allocatable :: zk(:)
            n=size(z0)
            allocate(zk(n))
            niter=0
            rho_old=0d0
            
            do
                niter=niter+1 !> we update number of iterations
                zk=matmul(A,z0) 
                rho_new=p_norm_vec(zk,2)/p_norm_vec(z0,2)
                if (abs(rho_new-rho_old)<tolerance .or. niter.eq.iter_max) exit
                z0=zk
                rho_old=rho_new
            end do
            rho=rho_new
        end subroutine

        !> Computes $L_\infty$ norm of a matrix
function norm_mat_inf(A) result(norm)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> m x n matrix
    real(kind=8) :: norm
            
    integer(kind=4) :: i,j,n,m
    real(kind=8) :: sum
    m=size(A,1)
    n=size(A,2)
    sum=0d0
    do i=1,n
        sum=sum+abs(A(1,i))
    end do
    norm=sum
    do i=2,m
        sum=0d0
        do j=1,n
            sum=sum+abs(A(i,j))
        end do
        if (sum>norm) then
            norm=sum
        end if
    end do
end function

!> Computes $L_1$ norm of a matrix
function norm_mat_1(A) result(norm)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> m x n matrix
    real(kind=8) :: norm
            
    integer(kind=4) :: i,j,n,m
    real(kind=8) :: sum
    m=size(A,1)
    n=size(A,2)
    sum=0
    do i=1,m
        sum=sum+abs(A(i,1))
    end do
    norm=sum
    do j=2,n
        sum=0
        do i=1,m
            sum=sum+abs(A(i,j))
        end do
        if (sum>norm) then
            norm=sum
        end if
    end do
end function

!> LU decomposition
subroutine LU(A,tol,L,U,error)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> square matrix
    real(kind=8), intent(in) :: tol !> tolerance
    real(kind=8), intent(out) :: L(:,:), U(:,:)
    logical, intent(out) :: error
    
    integer(kind=4) :: n,i,j
    real(kind=8) :: factor
    n=size(A,1)
    U=A
    L=0d0
    error=.false.
    do j=1,n
        if (abs(U(j,j))<tol) then
            error=.true.
            !error stop "Diagonal term is zero in LU decomposition"
        end if
        L(j,j)=1d0
        do i=j+1,n
            factor=U(i,j)/U(j,j)
            L(i,j)=factor
            U(i,1:n)=U(i,1:n)-factor*U(j,1:n)
        end do
    end do

end subroutine LU





function inf_norm_mat(A) result(norm)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> m x n matrix
    real(kind=8) :: norm
            
    integer(kind=4) :: i,j,n,m
    real(kind=8) :: sum
    m=size(A,1)
    n=size(A,2)
    sum=0
    do i=1,n
        sum=sum+abs(A(1,i))
    end do
    norm=sum
    do i=2,m
        sum=0
        do j=1,n
            sum=sum+abs(A(i,j))
        end do
        if (sum>norm) then
            norm=sum
        end if
    end do
end function

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

    function Householder(x)
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
    
    subroutine Gram_Schmidt_mat(v,u)
        implicit none
        real(kind=8), intent(in) :: v(:,:)
        real(kind=8), intent(out) :: u(:,:) !> base ortogonal
        integer(kind=4) :: k,j,n
        real(kind=8), allocatable :: sum_proy_ortog(:)
        n=size(v,1)
        allocate(sum_proy_ortog(n))
        u(1:n,1)=v(1:n,1)
        do k=2,n
            sum_proy_ortog=0d0
            do j=1,k-1
                sum_proy_ortog=sum_proy_ortog+proy_ortog(v(1:n,k),u(1:n,j))
            end do
            u(1:n,k)=v(1:n,k)-sum_proy_ortog
        end do
    end subroutine
    
    subroutine forward_substitution(L,b,x)
        implicit none
        real(kind=8), intent(in) :: L(:,:) !> lower triangular matrix
        real(kind=8), intent(in) :: b(:) !> vector
        real(kind=8), intent(out) :: x(:) !> solution of linear system
        integer(kind=4) :: j,k,n
        real(kind=8) :: sum
        n=size(b)
        x(1)=b(1)/L(1,1)
        do j=2,n
            sum=0d0
            do k=1,j-1
                sum=sum+L(j,k)*x(k)
            end do
            x(j)=(b(j)-sum)/L(j,j)
        end do
    end subroutine 


    !> Determinant of square matrix using LU decomposition
subroutine compute_det(A,tol,det,error)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> square matrix
    real(kind=8), intent(in) :: tol !> tolerance for determinant
    real(kind=8), intent(out) :: det !> determinant
    logical, intent(out) :: error
    
    real(kind=8) :: det_L, det_U
    real(kind=8), allocatable :: L(:,:), U(:,:)
    integer(kind=4) :: n,i
    
    error=.false.
    n=size(A,1)
    if (n/=size(A,2)) then
        error=.true.
        error stop "Matrix must be square (det)"
    end if
    if (n==1) then
        det=A(1,1)
    else if (n==2) then
        det=A(1,1)*A(2,2)-A(2,1)*A(1,2)
    else
        allocate(L(n,n),U(n,n))
        call LU(A,tol,L,U,error)
        det_L=1d0
        det_U=1d0
        do i=1,n
            det_L=det_L*L(i,i)
            det_U=det_U*U(i,i)
        end do
        det=det_L*det_U
    end if
end subroutine



!> Cholesky decomposition
subroutine Cholesky(A,L)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> positive definite, symmetric matrix
    real(kind=8), intent(out) :: L(:,:) !> lower triangular matrix with positive diagonal terms
    
    real(kind=8), allocatable :: prod_L_Lt(:,:) !> L*L^T
    integer(kind=4) :: i,j,k,n
    real(kind=8) :: sum1, sum2
    n=size(A,1)
    L=0d0
    do j=1,n
        if (j == 1) then
            L(j,j)=sqrt(A(j,j))
        else
            sum1=0d0
            do k=1,j-1
                sum1=sum1+L(j,k)**2
            end do
            L(j,j)=sqrt(A(j,j)-sum1)
        end if
        do i=j+1,n
            sum2=0d0
            do k=1,j-1
                sum2=sum2+L(i,k)*L(j,k)
            end do
            L(i,j)=(A(i,j)-sum2)/L(j,j)
        end do
    end do
    prod_L_Lt=matmul(L,transpose(L))
    if (all(A/=prod_L_Lt)) error stop "Wrong solution"
end subroutine

!> This subroutine performs backward substitution in LU decomposition of matrix
subroutine backward_substitution(U,b,x)
    implicit none
    real(kind=8), intent(in) :: U(:,:) !> upper triangular matrix
    real(kind=8), intent(in) :: b(:) !> independent term
    real(kind=8), intent(out) :: x(:) !> solution of linear system
    
    integer(kind=4) :: j,k,n
    real(kind=8) :: sum
    
    n=size(b)
    x(n)=b(n)/U(n,n)
    do j=1,n-1
        sum=0d0
        do k=n-j+1,n
            sum=sum+U(n-j,k)*x(k)
        end do
        x(n-j)=(b(n-j)-sum)/U(n-j,n-j)
    end do
end subroutine

subroutine eigenvectors_tridiag_toeplitz_matrix(A)
    class(tridiag_Toeplitz_matrix_c) :: A
    
    integer(kind=4) :: i,j
    real(kind=8) :: L2_norm_vj
    real(kind=8), parameter :: pi=4d0*atan(1d0)
    
    if (A%sub*A%super<=0d0) error stop "a*c must be positive"
    allocate(A%eigenvectors(A%num_cols,A%num_cols))
    do j=1,A%num_cols
        do i=1,A%num_cols
            A%eigenvectors(i,j)=(A%sub/A%super)**((i-1d0)/2d0)*sin((A%num_cols-j+1)*pi*i/(A%num_cols+1d0))
        end do
        L2_norm_vj=p_norm_vec(A%eigenvectors(:,j),2)
        A%eigenvectors(:,j)=A%eigenvectors(:,j)/L2_norm_vj
    end do
end subroutine

subroutine eigenvectors_tridiag_sym_matrix(a,b,lambda,v)
    real(kind=8), intent(in) :: a(:) !> diagonal elements
    real(kind=8), intent(in) :: b(:) !> non-diagonal elements
    real(kind=8), intent(in) :: lambda(:) !> eigenvalues
    real(kind=8), intent(out) :: v(:,:) !> eigenvectors
    
    real(kind=8), parameter :: epsilon=1d-12
    real(kind=8) :: L2_norm_vj
    
    integer(kind=4) :: i,j,k,n,m
    integer(kind=4) :: l(2)
    
    n=size(lambda)
    
    do j=1,n
        if (abs(a(1)-lambda(j))>=epsilon .and. abs(a(n)-lambda(j))>=epsilon) then
            v(1,j)=1d0
            v(2,j)=v(1,j)*(lambda(j)-a(1))/b(1)
            do k=3,n-1
                v(k,j)=v(k-1,j)*(lambda(j)-a(k-1))/b(k-1) - v(k-2,j)*b(k-2)/b(k-1)
            end do
            v(n,j)=v(n-1,j)*b(n-1)/(lambda(j)-a(n))
        else if (abs(a(1)-lambda(j))<epsilon) then
            v(2,j)=0d0
            v(1,j)=1d0
            do k=3,n-1
                v(k,j)=v(k-1,j)*(lambda(j)-a(k-1))/b(k-1) - v(k-2,j)*b(k-2)/b(k-1)
            end do
            v(n,j)=v(n-1,j)*b(n-1)/(lambda(j)-a(n))
        else
            v(n-1,j)=0d0
            v(n,j)=1d0
            do k=2,n-2
                v(n-k,j)=v(n-k+1,j)*(lambda(j)-a(n-k+1))/b(n-k) - v(n-k+2,j)*b(n-k+1)/b(n-k)
            end do
            v(1,j)=v(2,j)*b(1)/(lambda(j)-a(1))
        end if
    L2_norm_vj=p_norm_vec(v(:,j),2)
    v(:,j)=v(:,j)/L2_norm_vj
    end do
    do i=1,n
        do j=1,n
            if (abs(lambda(i)-lambda(j))<epsilon .and. i/=j) then
                l=[i,j]
            end if
        end do
    end do
    if (mod(n,2).eq.0) then
        m=n
    else
        m=n-1
    end if
end subroutine

subroutine check_eigenvectors(A,lambda,v,tolerance)
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


subroutine check_eigenvectors_tridiag_sym_matrix(this,tolerance)
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

subroutine compute_eigenvalues(this)
    class(sq_matrix_c) :: this
    
    real(kind=8), allocatable :: d(:),e(:)
    integer(kind=4) :: info
!> Lapack
    select type (this)
    class is (tridiag_sym_matrix_c)
        d=this%diag
        e=this%sub
        !call dsterf(size(d),d,e,info)
        this%eigenvalues=d
    end select
end subroutine

!> USAR LAPACK
subroutine compute_eigenvectors(this)
    class(sq_matrix_c) :: this
end subroutine

subroutine get_can_vec(this,tol,num_can_vec,ind_can_vec,num_non_can_vec,ind_non_can_vec)
!> Computes the number and indices of canonical and non-canonical vectors
implicit none
class(real_array_c), intent(in) :: this !> array object
real(kind=8), intent(in) :: tol !> tolerance
integer(kind=4), intent(out) :: num_can_vec !> number of canonical vectors
integer(kind=4), intent(out), allocatable :: ind_can_vec(:) !> indices of canonical vectors
integer(kind=4), intent(out) :: num_non_can_vec !> number of non-canonical vectors
integer(kind=4), intent(out), allocatable :: ind_non_can_vec(:) !> indices of non-canonical vectors
integer(kind=4) :: i,j
num_can_vec=0 !> we initialize number of canonical vectors
num_non_can_vec=0 !> we initialize number of non-canonical vectors
!> First loop: we count the number of canonical and non-canonical vectors
do i=1,this%num_cols
    if (abs(this%cols(i)%col_1(1)-1d0)<tol .and. &
        inf_norm_vec_real(this%cols(i)%col_1(2:this%cols(i)%dim))<tol) then
        num_can_vec=num_can_vec+1 !> we update number of canonical vectors
    else
        num_non_can_vec=num_non_can_vec+1 !> we update number of non-canonical vectors
    end if
end do
allocate(ind_can_vec(num_can_vec)) !> we allocate the array of indices of canonical vectors
allocate(ind_non_can_vec(num_non_can_vec)) !> we allocate the array of indices of canonical vectors
!> Second loop: we store the indices of canonical and non-canonical vectors
num_can_vec=0 !> we initialize number of canonical vectors
num_non_can_vec=0 !> we initialize number of non-canonical vectors
do i=1,this%num_cols
    if (abs(this%cols(i)%col_1(1)-1d0)<tol .and. &
        inf_norm_vec_real(this%cols(i)%col_1(2:this%cols(i)%dim))<tol) then
        num_can_vec=num_can_vec+1 !> we update number of canonical vectors
        ind_can_vec(num_can_vec)=i !> we store the index of the canonical vector
    else
        num_non_can_vec=num_non_can_vec+1 !> we update number of non-canonical vectors
        ind_non_can_vec(num_non_can_vec)=i !> we store the index of the non-canonical vector
    end if
end do
end subroutine
        
end module 