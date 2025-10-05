!> Computes logarithm activity coefficient of surface complex using Gaines-Thomas convention
!!> Activity of surface complex = equivalent fraction
subroutine compute_log_act_coeff_Gaines_Thomas(this,cation,CEC,log_act_coeff)
    use Gaines_Thomas_m, only: Gaines_Thomas_c, species_c
    implicit none
    class(Gaines_Thomas_c) :: this
    class(species_c), intent(in) :: cation !> adsorbed cation
    real(kind=8), intent(in) :: CEC !> cation exchange capacity
    real(kind=8), intent(out) :: log_act_coeff !> log_10 activity coefficient surface complex
        
    log_act_coeff=log10(cation%valence*1d0)-log10(CEC)
end subroutine