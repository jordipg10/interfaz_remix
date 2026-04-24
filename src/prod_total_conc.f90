!> Computes concentrations using eigendecomposition
subroutine prod_total_conc(this,A_mat,time)
    use diffusion_transient_m, only: diffusion_1D_transient_c
    use arrays_m, only: tridiag_matrix_c
    implicit none
    class(diffusion_1D_transient_c) :: this
    class(tridiag_matrix_c), intent(in) :: A_mat
    real(kind=8), intent(in), optional :: time
    
    real(kind=8), allocatable :: y0(:),b(:),c0(:),g(:)
    real(kind=8) :: sumj,sumk1,sumk2,sumk,t,conc_lim_n
    integer(kind=4) :: k,i,j,n
    
    n=this%spatial_discr%Num_targets
    
    y0=sqrt(this%F_mat%diag)*this%conc_init
    b=this%source_term_PDE/sqrt(this%F_mat%diag)
    
    if (.not. allocated(this%conc_obj%conc)) then
        allocate(this%conc_obj%conc(n))
    end if
    
    if (present(time)) then
        do i=1,n
            sumj=0d0
            do j=1,n
                sumk1=0d0
                do k=1,n
                    sumk1=sumk1+A_mat%eigenvectors(k,j)*sqrt(this%F_mat%diag(k))*this%conc_init(k)
                end do
                sumk2=0d0
                do k=1,n
                    sumk2=sumk2+A_mat%eigenvectors(k,j)*this%source_term_PDE(k)/sqrt(this%F_mat%diag(k))
                end do
                sumj=sumj+A_mat%eigenvectors(i,j)*(exp(-A_mat%eigenvalues(j)*time)*sumk1+((1d0-exp(-A_mat%eigenvalues(j)*time))/&
                A_mat%eigenvalues(j))*sumk2)
            end do
            this%conc_obj%conc(i)=sumj/sqrt(this%F_mat%diag(i))
        end do
    else
        do i=1,n
            sumj=0d0
            do j=1,n
                sumk=0d0
                do k=1,n
                    sumk=sumk+A_mat%eigenvectors(k,j)*this%source_term_PDE(k)/sqrt(this%F_mat%diag(k))
                end do
                sumj=sumj+A_mat%eigenvectors(i,j)*sumk/A_mat%eigenvalues(j)
            end do
            this%conc_obj%conc(i)=sumj/sqrt(this%F_mat%diag(i))
        end do
    end if
end subroutine