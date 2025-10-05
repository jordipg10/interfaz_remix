subroutine potencia(A,z0,tolerance,rho,niter)
!> Computes spectral radius
    use vectors_m
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