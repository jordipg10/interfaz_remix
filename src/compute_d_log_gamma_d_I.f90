!!> Computes derivative of log_10(activity coefficients) with respect to ionic activity
!subroutine compute_d_log_gamma_d_I(this,ionic_act,params_aq_sol,d_log_gamma_d_I)
!    use aq_phase_m
!    implicit none
!    class(aq_phase_c) :: this
!    real(kind=8), intent(in) :: ionic_act !> ionic activity
!    class(params_aq_sol_t), intent(in) :: params_aq_sol
!    real(kind=8), intent(out) :: d_log_gamma_d_I(:) !> derivative of log_10(activity coefficients) with respect to ionic activity (must be already allocated)
!    
!    integer(kind=4) :: i
!    
!    do i=1,this%num_species
!        d_log_gamma_d_I(i)=-(params_aq_sol%A*this%aq_species(i)%valence**2)/(2d0*sqrt(ionic_act)*(1d0+this%aq_species(i)%params_act_coeff%alpha*sqrt(ionic_act))**2) + this%aq_species(i)%params_act_coeff%beta
!    end do
!end subroutine