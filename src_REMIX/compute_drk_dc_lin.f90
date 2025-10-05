!> Computes linear kinetic reaction rate gradient
subroutine compute_drk_dc_lin(this,drk_dc)
    use lin_kin_reaction_m, only: lin_kin_reaction_c
    implicit none
    class(lin_kin_reaction_c), intent(in) :: this
    real(kind=8), intent(out) :: drk_dc(:) !> must be allocated
    
    drk_dc=0d0
    drk_dc(this%indices_aq_phase)=this%lambda
    
end subroutine