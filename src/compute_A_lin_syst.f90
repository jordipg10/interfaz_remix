subroutine compute_A_mat_lin_syst(this,theta,A_mat,k)
!> A*c^(k+1)=b
!> A=Id-theta*E
    use PDE_transient_m
    implicit none
    
    class(PDE_1D_transient_c), intent(in) :: this
    real(kind=8), intent(in) :: theta
    type(tridiag_matrix_c), intent(out) :: A_mat
    integer(kind=4), intent(in), optional :: k
    
    type(tridiag_matrix_c) :: E_mat
    
    call this%compute_E_mat(E_mat,k)
    
    A_mat%sub=-theta*E_mat%sub
    A_mat%diag=1d0-theta*E_mat%diag
    A_mat%super=-theta*E_mat%super

end subroutine