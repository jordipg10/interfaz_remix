!> \file Gaines_Thomas_m.f90
!> \brief Gaines-Thomas convention for cation exchange activity coefficients
!> \details
!>   This module implements the Gaines-Thomas convention for computing
!>   activity coefficients of adsorbed cations in ion exchange reactions.
!>   
!>   **Convention:**
!>   The Gaines-Thomas convention uses equivalent fractions for activities.
!>   Activity of adsorbed cation i: a_i = N_i (equivalent fraction)
!>   where N_i = z_i * c_i / CEC
!>   
!>   **Mathematical formulation:**
!>   For exchange reaction: z_j·X_i + z_i·M_j ⇌ z_j·X_j + z_i·M_i
!>   
!>   Activity coefficient:
!>   γ_i = (CEC / z_i)
!>   log10(γ_i) = log10(z_i) - log10(CEC)
!>   
!>   **Key properties:**
!>   - Rational thermodynamic convention
!>   - Activity coefficient depends on valence and CEC
!>   - Most rigorous from thermodynamic perspective
!>   - Default convention in PHREEQC
!>   
!>   **Applications:**
!>   - Clay mineral cation exchange
!>   - Zeolite ion exchange
!>   - Soil chemistry modeling
!>   
!>   **References:**
!>   - Gaines & Thomas (1953): Soil Sci. Soc. Am. Proc. 17:19-22
!>   - Appelo & Postma (2005): Geochemistry, groundwater and pollution
!>   - Parkhurst & Appelo (2013): PHREEQC Version 3 manual
!>   
!> \author Jordi Petchamé-Guerrero
!> \date October 2025
module Gaines_Thomas_m
    use exch_sites_conv_m !< Import abstract exchange sites convention base class
    implicit none !< Enforce explicit variable declarations
    save !< Preserve module state between calls
    
    !> \brief Gaines-Thomas convention type
    !> \details
    !>   Concrete implementation of the exchange sites convention interface
    !>   for the Gaines-Thomas formulation.
    !>   
    !>   **Extends:** exch_sites_conv_c (abstract base class)
    !>   
    !>   **Implements:**
    !>   - compute_log_act_coeff_ads_cat: Calculates log10(γ) using Gaines-Thomas formula
    !>   
    !>   **Formula:**
    !>   log10(γ_i) = log10(z_i) - log10(CEC)
    !>   
    !>   **No additional attributes:**
    !>   This type contains no member variables - all behavior is in the procedure.
    !>   
    !>   **Usage example:**
    !>   ```fortran
    !>   type(Gaines_Thomas_c) :: GT_conv
    !>   real(kind=8) :: log_gamma
    !>   integer(kind=4) :: z_Na = 1
    !>   real(kind=8) :: CEC = 0.1  ! mol/kg
    !>   call GT_conv%compute_log_act_coeff_ads_cat(z_Na, CEC, log_gamma)
    !>   ! Result: log_gamma = log10(1) - log10(0.1) = 0 - (-1) = 1
    !>   ```
    type, public, extends(exch_sites_conv_c) :: Gaines_Thomas_c !< Gaines-Thomas convention (subclass of exch_sites_conv_c)
        !< No additional attributes needed
    contains
        !> Bind deferred procedure to Gaines-Thomas implementation
        procedure :: compute_log_act_coeff_ads_cat=>compute_log_act_coeff_Gaines_Thomas
    end type
    
    contains

        !> \brief Compute log10 activity coefficient using Gaines-Thomas convention
        !> \details
        !>   Implements the Gaines-Thomas formula for activity coefficients
        !>   of adsorbed cations in cation exchange reactions.
        !>   
        !>   **Mathematical formulation:**
        !>   log10(γ_i) = log10(z_i) - log10(CEC)
        !>   
        !>   where:
        !>   - γ_i = activity coefficient of adsorbed cation i [-]
        !>   - z_i = valence (charge) of cation i [-]
        !>   - CEC = cation exchange capacity [meq/g] or [mol/kg]
        !>   
        !>   **Physical interpretation:**
        !>   Activity a_i = γ_i * N_i where N_i = z_i * c_i / CEC (equivalent fraction)
        !>   The Gaines-Thomas convention gives: a_i = z_i * c_i
        !>   
        !>   **Algorithm:**
        !>   1. Convert valence to real number
        !>   2. Compute log10(valence)
        !>   3. Compute log10(CEC)
        !>   4. Return difference: log10(valence) - log10(CEC)
        !>   
        !>   **Examples:**
        !>   - Na+ (z=1), CEC=0.1: log10(1) - log10(0.1) = 0 - (-1) = 1.0
        !>   - Ca2+ (z=2), CEC=0.1: log10(2) - log10(0.1) = 0.301 - (-1) = 1.301
        !>   - K+ (z=1), CEC=0.05: log10(1) - log10(0.05) = 0 - (-1.301) = 1.301
        !>   
        !>   **Note:**
        !>   CEC must be > 0 to avoid logarithm of zero/negative.
        !> \param[in,out] this Gaines-Thomas convention object
        !> \param[in] valence Valence (charge) of exchangeable cation (e.g., 1 for Na+, 2 for Ca2+) [-]
        !> \param[in] CEC Cation exchange capacity [meq/g] or [mol/kg] (must be > 0)
        !> \param[out] log_act_coeff Log10 activity coefficient [-]
        subroutine compute_log_act_coeff_Gaines_Thomas(this,valence,CEC,log_act_coeff)
            implicit none !< Enforce explicit declarations
            class(Gaines_Thomas_c) :: this !< Gaines-Thomas convention object
            integer(kind=4), intent(in) :: valence !< Valence of exchangeable cation [-]
            real(kind=8), intent(in) :: CEC !< Cation exchange capacity [meq/g] or [mol/kg]
            real(kind=8), intent(out) :: log_act_coeff !< Log10 activity coefficient [-]
        
            !< Compute log10(γ) = log10(valence) - log10(CEC)
            !< Convert valence to real (valence*1d0) before taking logarithm
            log_act_coeff=log10(valence*1d0)-log10(CEC)
        end subroutine

end module !< End of Gaines_Thomas_m module