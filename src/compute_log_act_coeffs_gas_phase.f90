!!> Computes log_10 activity coefficients gases
!!> (To be developed)
subroutine compute_log_act_coeffs_gas_phase(this,log_act_coeffs)
    use gas_phase_m
    implicit none
    class(gas_phase_c) :: this
    real(kind=8), intent(out) :: log_act_coeffs(:) !> log_10 activity coefficients (must be already allocated)
    
    !error stop "Subroutine 'compute_log_act_coeff_gas_phase' not implemented yet"
    log_act_coeffs=0d0
end subroutine