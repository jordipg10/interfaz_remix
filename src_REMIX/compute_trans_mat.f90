!subroutine compute_trans_mat(this,T_sub,T_diag,T_super,k)
!!> T: transition matrix (tridiagonal, negative stoich_mat_react_zonemi-definite)
!!> rows sum = 0 if r=0
!!> F*dc/dt=T*c+g
!>    use transport_transient_1_m
!>    use transport_m
!>    implicit none
!>    
!>    class(diffusion_c) :: this
!>    !integer(kind=4), intent(in) :: opcion
!>    real(kind=8), intent(out), optional :: T_sub(:),T_diag(:),T_super(:)
!>    integer(kind=4), intent(in), optional :: k
!>    
!>    !real(kind=8) :: sign_vel !> sign of velocity parameter
!>    !real(kind=8), allocatable :: sub(:),diag(:),Delta_r(:)
!>    !integer(kind=4) :: i,n
!>    !
!>    !n=this%spatial_discr%Num_targets
!>    
!>    select type (this)
!>    type is (transport_1D_c)
!>        call this%compute_trans_mat_tpt(T_sub,T_diag,T_super,k)
!>    type is (transport_1D_transient_c)
!>        call this%compute_trans_mat_tpt(T_sub,T_diag,T_super,k)
!>    class is (diffusion_c)
!>        call this%compute_trans_mat_diff_bis(T_sub,T_diag,T_super,k)
!>    end select
!end subroutine 