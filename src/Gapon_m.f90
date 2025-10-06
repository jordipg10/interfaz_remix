!> Gapon convention library
module Gapon_m
    use exch_sites_conv_m
    implicit none
    save
    type, public, extends(exch_sites_conv_c) :: Gapon_c !> Gapon class (subclass of "exch_sites_conv_c")
        
    contains
        procedure, public :: compute_log_act_coeff_ads_cat=>compute_log_act_coeff_Gapon
    end type
    
    contains
        subroutine compute_log_act_coeff_Gapon(this,valence,CEC,log_act_coeff)
        !> Activity of surface complex = molar fraction total number exchange sites
            implicit none
            class(Gapon_c) :: this
            integer(kind=4), intent(in) :: valence
            real(kind=8), intent(in) :: CEC
            real(kind=8), intent(out) :: log_act_coeff
            
            log_act_coeff=valence*log10(valence/CEC)

        end subroutine
        
    
    
      
end module