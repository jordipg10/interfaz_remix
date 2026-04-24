!> \file exch_sites_conv_m.f90
!> \brief Exchange sites convention module for cation exchange activity coefficients
!> \details
!>   This module defines an abstract interface for cation exchange conventions
!>   used to compute activity coefficients of adsorbed cations.
!>   
!>   **Purpose:**
!>   Provides a polymorphic framework for different activity coefficient models.
!>   
!>   **Conventions (subclasses):**
!>   - Gaines-Thomas: Rational activity coefficients (classical thermodynamics)
!>   - Gapon: Activity proportional to equivalent fraction
!>   - Vanselow: Activity proportional to mole fraction
!>   
!>   **Applications:**
!>   - Ion exchange on clay minerals (montmorillonite, illite)
!>   - Cation adsorption on zeolites
!>   - Surface complexation modeling
!>   - Reactive transport simulations
!>   
!>   **Design pattern:**
!>   Abstract base class with deferred procedure for polymorphism.
!>   Each subclass implements convention-specific activity coefficient calculation.
!>   
!>   **References:**
!>   - Appelo & Postma (2005): Geochemistry, groundwater and pollution
!>   - Parkhurst & Appelo (2013): PHREEQC manual
!>   
!> \author Jordi Petchamé-Guerrero
!> \date October 2025
module exch_sites_conv_m
    implicit none !< Enforce explicit variable declarations
    save !< Preserve module state between calls
    private !< Private module scope by default
    !> \brief Abstract base class for exchange sites conventions
    !> \details
    !>   Defines the interface for computing activity coefficients of adsorbed cations
    !>   in cation exchange reactions.
    !>   
    !>   **Polymorphism:**
    !>   Abstract type with deferred procedures allows runtime selection of convention.
    !>   
    !>   **Subclasses must implement:**
    !>   - compute_log_act_coeff_ads_cat: Calculate log10(γ) for adsorbed cation
    !>   
    !>   **Usage pattern:**
    !>   ```fortran
    !>   class(exch_sites_conv_c), pointer :: convention
    !>   type(Gaines_Thomas_c), target :: GT_conv
    !>   convention => GT_conv
    !>   call convention%compute_log_act_coeff_ads_cat(valence, CEC, log_gamma)
    !>   ```
    !>   
    !>   **Theory:**
    !>   Activity coefficient γ relates activity to concentration:
    !>   a_i = γ_i * c_i
    !>   
    !>   Different conventions use different composition variables:
    !>   - Equivalent fraction: N_i = z_i*c_i / CEC
    !>   - Mole fraction: X_i = c_i / Σc_j
    type, public, abstract :: exch_sites_conv_c !< Exchange sites convention (abstract base class)
    contains
        !> Compute log10 activity coefficient for adsorbed cation (deferred to subclasses)
        procedure(compute_log_act_coeff_ads_cat), public, deferred :: compute_log_act_coeff_ads_cat
    end type
    
    !> \brief Abstract interface for activity coefficient calculation
    !> \details
    !>   Defines the signature for the deferred procedure that computes
    !>   log10 activity coefficients of adsorbed cations.
    !>   
    !>   **Function signature:**
    !>   All subclasses must implement this interface with the exact signature.
    !>   
    !>   **Mathematical formulation:**
    !>   Computes: log10(γ_ads) where γ_ads is the activity coefficient
    !>   
    !>   **Convention-specific implementations:**
    !>   - Gaines-Thomas: log10(γ) = 0 (γ = 1, rational convention)
    !>   - Gapon: log10(γ) = (1-z)/2 * log10(CEC)
    !>   - Vanselow: log10(γ) = (z-1) * log10(Σc_j)
    !>   
    !>   where z = valence of cation
    abstract interface
        !> \brief Compute log10 activity coefficient for adsorbed cation
        !> \details
        !>   Interface for convention-specific activity coefficient calculation.
        !>   
        !>   **Input:**
        !>   - valence: Charge of exchangeable cation (e.g., Na+=1, Ca2+=2)
        !>   - CEC: Cation exchange capacity in [meq/g] or [mol/kg]
        !>   
        !>   **Output:**
        !>   - log_act_coeff: Base-10 logarithm of activity coefficient [-]
        !>   
        !>   **Examples:**
        !>   - Na+ (z=1): For Gaines-Thomas, log10(γ) = 0 → γ = 1
        !>   - Ca2+ (z=2): For Gapon, log10(γ) = -0.5*log10(CEC)
        !>   
        !>   **Note:**
        !>   Activity coefficient may depend on:
        !>   - Cation valence (always)
        !>   - CEC (some conventions)
        !>   - Total adsorbed concentration (Vanselow)
        !> \param[in,out] this Exchange sites convention object
        !> \param[in] valence Valence (charge) of exchangeable cation [-]
        !> \param[in] CEC Cation exchange capacity [meq/g] or [mol/kg]
        !> \param[out] log_act_coeff Log10 activity coefficient [-]
        subroutine compute_log_act_coeff_ads_cat(this,valence,CEC,log_act_coeff)
            import exch_sites_conv_c !< Import abstract base class
            implicit none !< Enforce explicit declarations
            class(exch_sites_conv_c) :: this !< Exchange sites convention object
            integer(kind=4), intent(in) :: valence !< Valence of exchangeable cation (e.g., 1 for Na+, 2 for Ca2+) [-]
            real(kind=8), intent(in) :: CEC !< Cation exchange capacity [meq/g] or [mol/kg]
            real(kind=8), intent(out) :: log_act_coeff !< Log10 activity coefficient [-]
        end subroutine

    end interface
    
    contains
    !> \brief Module procedures section
    !> \details
    !>   Currently empty - all procedures are deferred to subclasses.
    !>   
    !>   **Subclass modules:**
    !>   - Gaines_Thomas_m: Implements Gaines-Thomas convention
    !>   - Gapon_m: Implements Gapon convention
    !>   - Vanselow_m: Implements Vanselow convention
      
end module