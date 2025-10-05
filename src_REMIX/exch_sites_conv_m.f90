!> Exchange sites convention library
module exch_sites_conv_m
    use species_m, only: species_c
    implicit none
    save
    type, public, abstract :: exch_sites_conv_c !> exchange sites convention abstract superclass
    contains
        procedure(compute_log_act_coeff_ads_cat), public, deferred :: compute_log_act_coeff_ads_cat
    end type
    
    abstract interface
        subroutine compute_log_act_coeff_ads_cat(this,valence,CEC,log_act_coeff)
        !> Computes log_10 activity coefficient of an adsorbed cation
            import exch_sites_conv_c
            import species_c
            implicit none
            class(exch_sites_conv_c) :: this                        !> exchange sites convention
            integer(kind=4), intent(in) :: valence                  !> valence of exchangeable cation
            real(kind=8), intent(in) :: CEC                         !> cation exchange capacity
            real(kind=8), intent(out) :: log_act_coeff              !> log_10 activity coefficient
        end subroutine
        

    end interface
    
    contains
    
      
end module