subroutine compute_log_act_coeffs_surf_compl(this,conc,log_acts)
    use surf_compl_m
    implicit none
    class(surface_c) :: this
    real(kind=8), intent(in) :: conc(:)
    real(kind=8), intent(out) :: log_acts(:) ! must be allocated
    
    integer(kind=4) :: i
        
    ! We assume they behave ideally
    log_acts=0d0
end subroutine