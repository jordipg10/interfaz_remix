subroutine compute_A_mat_ODE(this,A_mat)
!> dy/dt=-Ay+b
!> A=-sqrt(F^(-1))*T*sqrt(F^(-1)): tridiagonal, symmetric, positive semi-definite
    use PDE_transient_m, only: PDE_1D_transient_c, tridiag_matrix_c, diag_matrix_c, prod_tridiag_diag_mat, prod_diag_tridiag_mat
    implicit none
    class(PDE_1D_transient_c), intent(in) :: this
    type(tridiag_matrix_c), intent(out) :: A_mat
    
    integer(kind=4) :: i,j
    type(diag_matrix_c) :: inv_F
    type(tridiag_matrix_c) :: aux_mat
    
    call inv_F%allocate_array(size(this%F_mat%diag))
    inv_F%diag=sqrt(1d0/this%F_mat%diag)
    
    aux_mat=prod_tridiag_diag_mat(this%trans_mat,inv_F)
    A_mat=prod_diag_tridiag_mat(inv_F,aux_mat)
    A_mat%diag=-A_mat%diag
    A_mat%sub=-A_mat%sub
    A_mat%super=-A_mat%super
end subroutine