!> Gaines-Thomas convention library
module Gaines_Thomas_m
    use exch_sites_conv_m
    implicit none
    save
    type, public, extends(exch_sites_conv_c) :: Gaines_Thomas_c !> Gaines-Thomas convention class (subclass of "exch_sites_conv_c")
        
    contains
        procedure, public :: compute_log_act_coeff_ads_cat=>compute_log_act_coeff_Gaines_Thomas
    end type
    
    contains

        subroutine compute_log_act_coeff_Gaines_Thomas(this,valence,CEC,log_act_coeff)
        !> Computes logarithm 10 activity coefficient of surface complex using Gaines-Thomas convention
        !! Activity of surface complex = equivalent fraction
            implicit none
            class(Gaines_Thomas_c) :: this
            integer(kind=4), intent(in) :: valence !> exchangeable cation
            real(kind=8), intent(in) :: CEC !> cation exchange capacity
            real(kind=8), intent(out) :: log_act_coeff !> log_10 activity coefficient surface complex
        
            log_act_coeff=log10(valence*1d0)-log10(CEC)
        end subroutine
        

    
    
      
end module