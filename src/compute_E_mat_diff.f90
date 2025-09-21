!subroutine compute_E_mat_diff(this,E_mat,k)
!!> E=Delta_t(k)*F^(-1)*T (tridiagonal, negative stoich_mat_react_zonemi-definite)
!!> rows sum = 0 if r=0
!>    use diffusion_transient_m
!>    implicit none
!>    
!>    class(diffusion_transient_c) :: this
!>    class(array_c), pointer, intent(out) :: E_mat
!>    integer(kind=4), intent(in), optional :: k
!>    
!>    integer(kind=4) :: n
!>    type(tridiag_sym_matrix_c), target :: E_mat_sym
!>    type(tridiag_matrix_c), target :: E_mat_non_sym
!>    
!>    n=this%spatial_discr%Num_targets
!>    
!>    select type (trans_mat=>this%trans_mat)
!>    type is (tridiag_sym_matrix_c)
!>        E_mat_sym%sub=trans_mat%sub
!>        E_mat_sym%diag=trans_mat%diag
!>        select type (time_discr=>this%time_discr)
!>        type is (time_discr_homog_c)
!>            E_mat_sym%sub=E_mat_sym%sub*time_discr%Delta_t
!>            E_mat_sym%diag=E_mat_sym%diag*time_discr%Delta_t
!>        type is (time_discr_heterog_c)
!>            E_mat_sym%sub=E_mat_sym%sub*time_discr%Delta_t(k)
!>            E_mat_sym%diag=E_mat_sym%diag*time_discr%Delta_t(k)
!>        end select
!>        E_mat_sym%diag=E_mat_sym%diag/this%F_mat%diag
!>    end select
!end subroutine