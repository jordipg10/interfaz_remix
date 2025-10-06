!> Davies model to compute log_10 activity coefficient aqueous species
!> Valid for solutions such that \f$0.1 < I \le 0.7\f$
subroutine Davies(this,ionic_act,A,log_act_coeff)
    use aq_species_m, only: aq_species_c
    implicit none
    class(aq_species_c), intent(in) :: this
    real(kind=8), intent(in) :: ionic_act !> ionic activity
    real(kind=8), intent(in) :: A !> parameter
    real(kind=8), intent(out) :: log_act_coeff !> log_10 activity coefficient
    
    log_act_coeff=-A*(this%valence**2)*((sqrt(ionic_act)/(1d0+sqrt(ionic_act)))-3d-1*ionic_act)
    
end subroutine