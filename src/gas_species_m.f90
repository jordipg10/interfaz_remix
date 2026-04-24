!> @file gas_species_m.f90
!> @brief Gas species module for reactive transport modeling
!> @details This module defines the gas_species_c class (subclass of species_c) that represents
!> gaseous species with thermodynamic properties needed for equation of state (EOS) calculations.
!>
!> @par Gas Properties:
!> Three critical properties characterize gas behavior in non-ideal conditions:
!> - **Critical temperature (Tc)**: Temperature above which gas cannot be liquefied [K]
!> - **Critical pressure (Pc)**: Pressure at critical point [atm or Pa]
!> - **Acentric factor (ω)**: Measure of molecular non-sphericity [-]
!>
!> @par Equation of State Applications:
!> These properties are used in cubic equations of state (EOS) like:
!> - **Peng-Robinson EOS**: P = RT/(V-b) - a·α(T)/(V(V+b)+b(V-b))
!> - **Redlich-Kwong EOS**: P = RT/(V-b) - a/(√T·V(V+b))
!>
!> where parameters a, b depend on Tc, Pc, and ω.
!>
!> @par Acentric Factor Definition:
!> ω = -log₁₀(Pᵣˢᵃᵗ|Tᵣ=0.7) - 1.0
!> where Pᵣˢᵃᵗ is reduced vapor pressure at reduced temperature Tᵣ = 0.7
!>
!> @par Typical Values:
!> - Simple molecules (Ar, CH₄): ω ≈ 0.0-0.01
!> - Normal molecules (CO₂, N₂): ω ≈ 0.2-0.4
!> - Complex molecules (hydrocarbons): ω > 0.4
!>
!> @see species_m For base species class
!> @author Jordi Petchamé-Guerrero
!> @date November 2025

module gas_species_m
    use species_m, only: species_c !< Import base species class for inheritance
    implicit none !< Enforce explicit variable declarations for type safety
    save !< Preserve module variables between procedure calls
    private !< Private module scope: internal details hidden from outside modules
    !> @brief Gas species class with thermodynamic properties
    !> @details This type extends the base species class to include critical properties
    !> needed for non-ideal gas behavior calculations using equations of state.
    !>
    !> @par Usage in Reactive Transport:
    !> Gas species participate in:
    !> - Gas-water equilibria (Henry's law, fugacity calculations)
    !> - Gas phase reactions (oxidation, combustion)
    !> - Multiphase flow (gas dissolution/exsolution)
    !>
    !> @par EOS Parameter Calculation:
    !> For Peng-Robinson EOS:
    !> - a = 0.45724·R²·Tc²/Pc
    !> - b = 0.07780·R·Tc/Pc
    !> - α(T) = [1 + κ(1-√(T/Tc))]² where κ = f(ω)
    !>
    !> @extends species_c
    type, public, extends(species_c) :: gas_species_c !< Gas species class with critical properties
        real(kind=8) :: crit_temp !< Critical temperature Tc [K] - temperature above which gas cannot be liquefied
        real(kind=8) :: crit_press !< Critical pressure Pc [atm or Pa] - pressure at critical point
        real(kind=8) :: acentric_fact !< Acentric factor ω [-] - measure of molecular non-sphericity (0 for spherical)
    contains
    !> @name Setter Methods
    !> @{
        procedure :: set_crit_temp !< Set critical temperature
        procedure :: set_crit_press !< Set critical pressure
        procedure :: set_acentric_fact !< Set acentric factor
    !> @}
    end type
    
    contains
        !> @brief Set critical temperature
        !> @details Assigns the critical temperature Tc, the temperature above which the gas
        !> cannot be liquefied regardless of pressure applied. Used in equation of state calculations.
        !>
        !> @par Typical Critical Temperatures:
        !> - He: 5.2 K
        !> - N₂: 126.2 K
        !> - CO₂: 304.1 K
        !> - H₂O: 647.1 K
        !>
        !> @param[inout] this Gas species object
        !> @param[in] crit_temp Critical temperature Tc [K]
        subroutine set_crit_temp(this,crit_temp)
            implicit none !< Enforce explicit variable declarations
            class(gas_species_c), intent(inout) :: this !< Gas species object [-]
            real(kind=8), intent(in) :: crit_temp !< Critical temperature [K]
            this%crit_temp=crit_temp !< Assign critical temperature [K]
        end subroutine
        
        !> @brief Set critical pressure
        !> @details Assigns the critical pressure Pc, the pressure at the critical point where
        !> liquid and gas phases become indistinguishable. Used in equation of state calculations.
        !>
        !> @par Typical Critical Pressures:
        !> - He: 2.27 atm
        !> - N₂: 33.9 atm
        !> - CO₂: 72.8 atm
        !> - H₂O: 217.7 atm
        !>
        !> @param[inout] this Gas species object
        !> @param[in] crit_press Critical pressure Pc [atm or Pa]
        subroutine set_crit_press(this,crit_press)
            implicit none !< Enforce explicit variable declarations
            class(gas_species_c), intent(inout) :: this !< Gas species object [-]
            real(kind=8), intent(in) :: crit_press !< Critical pressure [atm or Pa]
            this%crit_press=crit_press !< Assign critical pressure [atm or Pa]
        end subroutine
        
        !> @brief Set acentric factor
        !> @details Assigns the acentric factor ω, which quantifies molecular non-sphericity
        !> and is used to improve equation of state predictions. Defined as:
        !> ω = -log₁₀(Pᵣˢᵃᵗ|Tᵣ=0.7) - 1.0
        !>
        !> @par Physical Interpretation:
        !> - ω = 0: Perfectly spherical molecule (e.g., noble gases)
        !> - ω > 0: Elongated or polar molecules (most real gases)
        !> - Larger ω: More deviation from ideal spherical behavior
        !>
        !> @par Typical Acentric Factors:
        !> - Ar: 0.000 (reference: spherical)
        !> - N₂: 0.040
        !> - CO₂: 0.225
        !> - H₂O: 0.344
        !> - n-octane: 0.398
        !>
        !> @param[inout] this Gas species object
        !> @param[in] acentric_fact Acentric factor ω [-]
        subroutine set_acentric_fact(this,acentric_fact)
            implicit none !< Enforce explicit variable declarations
            class(gas_species_c), intent(inout) :: this !< Gas species object [-]
            real(kind=8), intent(in) :: acentric_fact !< Acentric factor [-]
            this%acentric_fact=acentric_fact !< Assign acentric factor [-]
        end subroutine
        
end module