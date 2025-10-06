module vectors_m
    implicit none
    save
    type, public, abstract :: vector_c
        integer(kind=4) :: dim
    contains
        procedure, public :: allocate_vector
        procedure, public :: set_dim
    end type
    
    type, public, extends(vector_c) :: vector_real_c
        real(kind=8), allocatable :: col_1(:)
    end type
    
    type, public, extends(vector_c) :: vector_int_c
        integer(kind=4), allocatable :: col_1(:)
    end type
    
    contains
        function outer_prod_vec(x,y) !> Computes outer product of two vectors
            implicit none
            real(kind=8), intent(in) :: x(:), y(:) !> vector coefficients
            real(kind=8), allocatable :: outer_prod_vec(:,:)
            integer(kind=4) :: n,m,i,j
            n=size(x)
            m=size(y)
            allocate(outer_prod_vec(n,m))
            do i=1,n
                do j=1,m
                    outer_prod_vec(i,j)=x(i)*y(j)
                end do
            end do
        end function outer_prod_vec
        
        function p_norm_vec(x,p) !> Computes p-norm
            implicit none
            real(kind=8), intent(in) :: x(:) !> vector
            integer(kind=4), intent(in) :: p
            real(kind=8) :: p_norm_vec
            
            integer(kind=4) :: i,n
            real(kind=8) :: sum
            real(kind=8), parameter :: epsilon=1d-9
            if (p<1) error stop "p must be >= 1"
            n=size(x)
            sum=0d0
            do i=1,n
                sum=sum+abs(x(i))**p
            end do
            if (p==2) then
                p_norm_vec=sqrt(sum)
            else
                p_norm_vec=sum**(1d0/p)
            end if
        end function p_norm_vec
        
        function inf_norm_vec_real(x) !> Computes infinite norm
            real(kind=8), intent(in) :: x(:)
            real(kind=8) :: inf_norm_vec_real
            inf_norm_vec_real=maxval(abs(x))
        end function
        
        function inf_norm_vec_int(x) !> Computes infinite norm
            integer(kind=4), intent(in) :: x(:)
            integer(kind=4) :: inf_norm_vec_int
            inf_norm_vec_int=maxval(abs(x))
        end function 
        
        function proy_ortog(u,v) !> Computes orthogonal projection
            real(kind=8), intent(in) :: u(:), v(:)
            real(kind=8), allocatable :: proy_ortog(:)
            integer(kind=4) :: n
            if (size(u)/=size(v)) error stop "u and v must have same dimension"
            n=size(u)
            allocate(proy_ortog(n))
            proy_ortog=(dot_product(u,v)/dot_product(v,v))*v
        end function
        
        function sum_squares(x) !> Computes sum squares elements vector
            implicit none
            real(kind=8), intent(in) :: x(:)
            real(kind=8) :: sum_squares
            integer(kind=4) :: i
            sum_squares=0d0
            do i=1,size(x)
                sum_squares=sum_squares+x(i)**2
            end do
        end function 
        
        subroutine allocate_vector(this,dim)
            implicit none
            class(vector_c) :: this
            integer(kind=4), intent(in), optional :: dim
            if (present(dim)) then
                this%dim=dim
            end if
            select type (this)
            class is (vector_real_c)
                allocate(this%col_1(this%dim))
            class is (vector_int_c)
                allocate(this%col_1(this%dim))
            end select
        end subroutine
        
        subroutine set_dim(this,dim)
            implicit none
            class(vector_c) :: this
            integer(kind=4), intent(in) :: dim
            this%dim=dim
        end subroutine
end module