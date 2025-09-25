subroutine compute_log_Jacobian_act_coeffs_surf_compl(this,ionic_act,log_act_coeffs,conc,log_Jacobian_act_coeffs)
    use surf_compl_m
    implicit none
    class(surface_c) :: this
    real(kind=8), intent(in) :: ionic_act
    real(kind=8), intent(in) :: log_act_coeffs(:)
    real(kind=8), intent(in) :: conc(:) ! concentration of surface complexes in a given target
    real(kind=8), intent(out) :: log_Jacobian_act_coeffs(:,:) ! must be allocated
    
    ! We assume they behave ideally
    log_Jacobian_act_coeffs=0d0
end subroutine