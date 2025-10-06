!> Solves linear system Ax=b using Jacobi iterative method
subroutine Jacobi(A,b,x0,x,niter)
    use vectors_m
    implicit none
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