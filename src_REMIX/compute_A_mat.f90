subroutine compute_A_mat(this,theta,E_mat)
!> A_k*c^(k+1)=b_k
!> A_k=Id-theta*E_k
    use PDE_transient_m, only: PDE_1D_transient_c, tridiag_matrix_c
    implicit none
    
    class(PDE_1D_transient_c) :: this
    real(kind=8), intent(in) :: theta
    class(tridiag_matrix_c), intent(in) :: E_mat
        
    
    this%A_mat%sub=-theta*E_mat%sub
    this%A_mat%diag=1d0-theta*E_mat%diag
    this%A_mat%super=-theta*E_mat%super

end subroutine