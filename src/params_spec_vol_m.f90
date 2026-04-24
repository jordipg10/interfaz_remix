!> @file params_spec_vol_m.f90
!> @brief Specific volume parameters module for aqueous species
!> @details This module contains parameters to compute conventional specific volume (partial molar volume)
!> of solute species in aqueous solution using the Redlich-type equation as implemented in PHREEQC.
!>
!> @par Specific Volume Theory:
!> The specific volume (partial molar volume) V̄ of an ion in solution depends on ionic strength:
!> V̄(I) = V̄⁰ + ΔV̄(I)
!> where V̄⁰ is the standard state value and ΔV̄(I) is the ionic strength correction.
!>
!> @par Redlich Equation:
!> The Redlich-type equation relates specific volume to ionic strength I:
!> V̄ = a₁ + a₂/(ω+I) + a₃/(ω+I)² + a₄/(ω+I)³ + W·I + ion_params
!> where ω is a parameter related to closest approach distance.
!>
!> @par Applications:
!> - Density calculations for aqueous solutions
!> - Volume changes during reactions
!> - Pressure effects on equilibria
!>
!> @see PHREEQC documentation for Redlich equation implementation
!> @author Jordi Petchamé-Guerrero
!> @date November 2025

module params_spec_vol_m
    implicit none !< Enforce explicit variable declarations for type safety
    save !< Preserve module variables between procedure calls
    private !< Default visibility is private; expose only selected entities
    !> @brief Parameters for Redlich-type specific volume equation
    !> @details This type stores the 10 coefficients used in the Redlich equation (PHREEQC convention)
    !> for computing the partial molar volume of aqueous species as a function of ionic strength.
    !>
    !> @par Equation Form:
    !> V̄ = a₁ + a₂/(ω+I) + a₃/(ω+I)² + a₄/(ω+I)³ + W·I + DH_correction(ion_size, i₁, i₂, i₃, i₄)
    !> where I is ionic strength [mol/L]
    !>
    !> @par Parameter Meanings:
    !> - a₁, a₂, a₃, a₄: Polynomial coefficients for 1/(ω+I) expansion [cm³/mol]
    !> - W: Linear ionic strength coefficient [cm³·L/mol²]
    !> - ion_size_param: Ion size parameter å for Debye-Hückel term [Angstroms]
    !> - i₁, i₂, i₃, i₄: Additional ionic strength correction terms [various units]
    type, public :: params_spec_vol_Redlich_c !< Redlich-type specific volume parameters (PHREEQC)
        real(kind=8) :: a1 !< First polynomial coefficient in Redlich equation [cm³/mol]
        real(kind=8) :: a2 !< Second coefficient: numerator for 1/(ω+I) term [cm³/mol]
        real(kind=8) :: a3 !< Third coefficient: numerator for 1/(ω+I)² term [cm³·L/mol²]
        real(kind=8) :: a4 !< Fourth coefficient: numerator for 1/(ω+I)³ term [cm³·L²/mol³]
        real(kind=8) :: W !< Linear ionic strength coefficient (omega parameter) [cm³·L/mol²]
        real(kind=8) :: ion_size_param !< Ion size parameter å for Debye-Hückel correction [Angstroms]
        real(kind=8) :: i1 !< First ionic strength correction parameter [cm³·L/mol²]
        real(kind=8) :: i2 !< Second ionic strength correction parameter [cm³·L²/mol³]
        real(kind=8) :: i3 !< Third ionic strength correction parameter [cm³·L³/mol⁴]
        real(kind=8) :: i4 !< Fourth ionic strength correction parameter [cm³·L⁴/mol⁵]
    contains
        procedure :: set_params !< Assign all 10 Redlich parameters from coefficient array
    end type
    
    contains
        !> @brief Set all Redlich specific volume parameters from coefficient array
        !> @details Assigns the 10 coefficients required for the Redlich equation in the order
        !> used by PHREEQC database files. Performs dimension checking to ensure exactly 10 values.
        !>
        !> @par Coefficient Order (PHREEQC convention):
        !> coeffs(1) = a1, coeffs(2) = a2, coeffs(3) = a3, coeffs(4) = a4, coeffs(5) = W,
        !> coeffs(6) = ion_size_param, coeffs(7) = i1, coeffs(8) = i2, coeffs(9) = i3, coeffs(10) = i4
        !>
        !> @param[inout] this Redlich parameters object to populate
        !> @param[in] coeffs Array of 10 Redlich coefficients in PHREEQC order [various units]
        subroutine set_params(this,coeffs)
            implicit none !< Enforce explicit variable declarations
            class(params_spec_vol_Redlich_c), intent(inout) :: this !< Redlich parameters object to set [-]
            real(kind=8), intent(in) :: coeffs(:) !< Coefficient array, must contain exactly 10 values in PHREEQC order [various units]
            
            !> Check that coefficient array has exactly 10 elements (required for Redlich equation)
            if (size(coeffs)/=10) then
                error stop "Dimension error in coefficients specific volume Redlich" !< Terminate if wrong array size
            else
                !> Assign coefficients in PHREEQC order
                this%a1=coeffs(1) !< First polynomial coefficient [cm³/mol]
                this%a2=coeffs(2) !< Second polynomial coefficient [cm³/mol]
                this%a3=coeffs(3) !< Third polynomial coefficient [cm³·L/mol²]
                this%a4=coeffs(4) !< Fourth polynomial coefficient [cm³·L²/mol³]
                this%W=coeffs(5) !< Omega parameter (linear I term) [cm³·L/mol²]
                this%ion_size_param=coeffs(6) !< Ion size parameter å [Angstroms]
                this%i1=coeffs(7) !< First ionic strength correction [cm³·L/mol²]
                this%i2=coeffs(8) !< Second ionic strength correction [cm³·L²/mol³]
                this%i3=coeffs(9) !< Third ionic strength correction [cm³·L³/mol⁴]
                this%i4=coeffs(10) !< Fourth ionic strength correction [cm³·L⁴/mol⁵]
            end if
        end subroutine
end module