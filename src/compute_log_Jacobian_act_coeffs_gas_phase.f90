!> Computes Jacobian log_10 activity coefficients gases with respect to log_10 concentrations
subroutine compute_log_Jacobian_act_coeffs_gas_phase(this,log_act_coeffs,log_Jacobian_act_coeffs)
    use gas_phase_m
    implicit none
    class(gas_phase_c) :: this
    real(kind=8), intent(in) :: log_act_coeffs(:) !> log_10 activity coefficients gases
    real(kind=8), intent(out) :: log_Jacobian_act_coeffs(:,:) !> (must be allocated)
    
    !error stop "Activity coefficients of gas phase not implemented yet"
    log_Jacobian_act_coeffs=0d0
end subroutine