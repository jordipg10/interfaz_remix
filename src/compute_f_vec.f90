subroutine compute_f_vec(this,k)
!> f=Delta_t(k)*inv(F)*g
    use PDE_transient_m, only: PDE_1D_transient_c
    use time_discr_m, only: time_discr_homog_c, time_discr_heterog_c
    implicit none
    class(PDE_1D_transient_c) :: this
    integer(kind=4), intent(in), optional :: k
    
    this%f_vec=this%source_term_PDE
    select type (time=>this%time_discr)
    type is (time_discr_homog_c)
        this%f_vec=this%f_vec*time%Delta_t/this%F_mat%diag
    type is (time_discr_heterog_c)
        this%f_vec=this%f_vec*time%Delta_t(k)/this%F_mat%diag
    end select
end subroutine