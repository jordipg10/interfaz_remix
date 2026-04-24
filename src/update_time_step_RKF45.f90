subroutine update_time_step_RKF45_new(Delta_t_old,tolerance,conc_RK4,conc_RK5,Delta_t_new)
    use PDE_transient_m
    use vectors_m
    implicit none
    real(kind=8), intent(in) :: Delta_t_old
    real(kind=8), intent(in) :: tolerance
    real(kind=8), intent(in) :: conc_RK4(:)
    real(kind=8), intent(in) :: conc_RK5(:)
    real(kind=8), intent(out) :: Delta_t_new
    
    do while(p_norm_vec(conc_RK4-conc_RK5,2)>=tolerance)
        Delta_t_new=9d-1*Delta_t_old*(tolerance/p_norm_vec(conc_RK4-conc_RK5,2))**0.2
    end do
end subroutine update_time_step_RKF45_new