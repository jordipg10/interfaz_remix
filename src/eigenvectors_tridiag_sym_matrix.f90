subroutine eigenvectors_tridiag_sym_matrix(a,b,lambda,v)
    use vectors_m, only: p_norm_vec
    implicit none
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