subroutine QR_Householder(A,Q,R)
    use matrices_m, only : id_matrix, Householder !>
    use vectors_m
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> square matrix
    real(kind=8), intent(out) :: Q(:,:), R(:,:)
    !real(kind=8), external :: p_norm_vec, outer_prod_vec
    real(kind=8), allocatable :: Ak(:,:), Ak_tilde(:,:), Hk(:,:), Hk_tilde(:,:), Q_aux(:,:), a1(:)
    real(kind=8) :: beta
    integer(kind=4) :: n,i
    n=size(A,1)
    allocate(Ak(n,n),Hk(n,n),Q_aux(n,n),a1(n))
    Ak=A !A0
    a1=Ak(1:n,1)
    Hk=Householder(a1) !H1
    Q=id_matrix(n)
    do i=1,n-1
        R=matmul(Hk,Ak)
        Ak=R !A1
        Ak_tilde=Ak((i+1):n,(i+1):n) !Atilde_1
        a1=Ak_tilde(1:(n-i),1)
        Hk_tilde=Householder(a1) !Htilde_2
        Hk=0d0
        Hk(1:i,1:i)=id_matrix(i)
        Hk((i+1):n,(i+1):n)=Hk_tilde !H2
        Q_aux=matmul(Q,Hk)
        Q=Q_aux
    end do
end subroutine QR_Householder