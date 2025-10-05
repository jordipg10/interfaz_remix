!> Computes Jacobian log_10 activity coefficients aqueous species with respect to log_10 concentrations
!! We assume all primary species are aqueous
subroutine compute_log_Jacobian_act_coeffs_aq_phase(this,out_prod,conc,log_Jacobian_act_coeffs)
    use aq_phase_m, only: aq_phase_c
    use matrices_m, only: diag_matrix_c
    implicit none
    class(aq_phase_c) :: this
    real(kind=8), intent(in) :: out_prod(:,:) !> subset outer product between d_log_gamma_d_I and z^2
    real(kind=8), intent(in) :: conc(:) !> subset concentrations aqueous species
    real(kind=8), intent(out) :: log_Jacobian_act_coeffs(:,:) !> (must be allocated)

    type(diag_matrix_c) :: conc_diag
    
    call conc_diag%set_diag_matrix(conc) !> chapuza
        
    log_Jacobian_act_coeffs=5d-1*log(1d1)*conc_diag%prod_mat_diag_mat(out_prod)
     
end subroutine