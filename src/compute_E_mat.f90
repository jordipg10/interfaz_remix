subroutine compute_E_mat(this,E_mat,k)
!> E=Delta_t(k)*inv(F)*T (tridiagonal, negative semi-definite)
!> rows sum = 0 if r=0
    use PDE_transient_m, only: PDE_1D_transient_c, time_discr_homog_c, time_discr_heterog_c, tridiag_matrix_c
    implicit none
    
    class(PDE_1D_transient_c), intent(in) :: this
    type(tridiag_matrix_c), intent(out) :: E_mat
    integer(kind=4), intent(in), optional :: k
    
    integer(kind=4) :: j,n
    
    n=this%spatial_discr%Num_targets

    if (this%spatial_discr%adapt_ref.eq.1) then
        call this%compute_trans_mat_PDE()
        call this%compute_F_mat_PDE()
    end if
    E_mat%sub=this%trans_mat%sub
    E_mat%diag=this%trans_mat%diag
    E_mat%super=this%trans_mat%super
    select type (time_discr=>this%time_discr)
    type is (time_discr_homog_c)
        E_mat%sub=E_mat%sub*time_discr%Delta_t
        E_mat%diag=E_mat%diag*time_discr%Delta_t
        E_mat%super=E_mat%super*time_discr%Delta_t
    type is (time_discr_heterog_c)
        E_mat%sub=E_mat%sub*time_discr%Delta_t(k)
        E_mat%diag=E_mat%diag*time_discr%Delta_t(k)
        E_mat%super=E_mat%super*time_discr%Delta_t(k)
    end select
    do j=1,n-1
        E_mat%super(j)=E_mat%super(j)/this%F_mat%diag(j)
        E_mat%sub(j)=E_mat%sub(j)/this%F_mat%diag(j+1)
        E_mat%diag(j)=E_mat%diag(j)/this%F_mat%diag(j)
    end do
    E_mat%diag(n)=E_mat%diag(n)/this%F_mat%diag(n)
end subroutine 