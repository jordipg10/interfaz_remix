!> \file Vanselow_m.f90
!> \brief Vanselow cation exchange convention module.
!> \details
!> Implements the Vanselow convention for computing log activity coefficients
!> of adsorbed cations. In the Vanselow convention, the activity of an
!> adsorbed species equals its molar fraction:
!> \f[
!>   \log_{10}(\gamma_i) = \log_{10}\!\left(\frac{z_i}{\text{CEC}}\right)
!> \f]
!>
!> \see Gaines_Thomas_m, Gapon_m, exch_sites_conv_m
!> \author Jordi
!> \date Unknown
!> \ingroup chemistry

!> \brief Vanselow convention module.
module Vanselow_m
    use exch_sites_conv_m
    implicit none
    save
    !> \brief Vanselow exchange convention subclass.
    type, public, extends(exch_sites_conv_c) :: Vanselow_c
    contains
        procedure :: compute_log_act_coeff_ads_cat=>compute_log_act_coeff_Vanselow  !< Vanselow activity coefficient
    end type
    
    contains
        !> \brief Compute log activity coefficient using Vanselow convention.
        !> \details Activity = molar fraction.
        !> \f$ \log_{10}(\gamma) = \log_{10}(z / \text{CEC}) \f$
        !> \param[in,out] this     Vanselow object
        !> \param[in]     valence  Cation valence \f$z_i\f$ [-]
        !> \param[in]     CEC      Cation exchange capacity [eq/kg]
        !> \param[out]    log_act_coeff  Log10 activity coefficient [-]
        subroutine compute_log_act_coeff_Vanselow(this,valence,CEC,log_act_coeff)
            implicit none
            class(Vanselow_c) :: this
            integer(kind=4), intent(in) :: valence
            real(kind=8), intent(in) :: CEC
            real(kind=8), intent(out) :: log_act_coeff
        
            log_act_coeff=log10(valence*1d0/CEC)
        end subroutine
    
    
    
      
end module