!> \file Gapon_m.f90
!> \brief Gapon cation exchange convention module.
!> \details
!> Implements the Gapon convention for computing log activity coefficients
!> of adsorbed cations on exchange surfaces. In the Gapon convention:
!> \f[
!>   \log_{10}(\gamma_i) = z_i \, \log_{10}\!\left(\frac{z_i}{\text{CEC}}\right)
!> \f]
!> where \f$z_i\f$ is the cation valence and CEC is the cation exchange
!> capacity [eq/kg].
!>
!> \see Gaines_Thomas_m, Vanselow_m, exch_sites_conv_m
!> \author Jordi
!> \date Unknown
!> \ingroup chemistry

!> \brief Gapon convention module.
module Gapon_m
    use exch_sites_conv_m
    implicit none
    save
    !> \brief Gapon exchange convention subclass.
    type, public, extends(exch_sites_conv_c) :: Gapon_c
    contains
        procedure :: compute_log_act_coeff_ads_cat=>compute_log_act_coeff_Gapon  !< Gapon activity coefficient
    end type
    
    contains
        !> \brief Compute log activity coefficient using Gapon convention.
        !> \details Activity = molar fraction times total exchange sites.
        !> \f$ \log_{10}(\gamma) = z \, \log_{10}(z / \text{CEC}) \f$
        !> \param[in,out] this     Gapon object
        !> \param[in]     valence  Cation valence \f$z_i\f$ [-]
        !> \param[in]     CEC      Cation exchange capacity [eq/kg]
        !> \param[out]    log_act_coeff  Log10 activity coefficient [-]
        subroutine compute_log_act_coeff_Gapon(this,valence,CEC,log_act_coeff)
            implicit none
            class(Gapon_c) :: this
            integer(kind=4), intent(in) :: valence
            real(kind=8), intent(in) :: CEC
            real(kind=8), intent(out) :: log_act_coeff
        
            log_act_coeff=valence*log10(valence*1d0/CEC)
        end subroutine
        
    
    
      
end module