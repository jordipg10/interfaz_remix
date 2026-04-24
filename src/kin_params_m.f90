!> \file kin_params_m.f90
!> \brief Defines abstract base class for kinetic reaction parameters.
!> \details This module provides the abstract kinetic parameters type that serves as a base class for all kinetic reaction
!>          parameter types (e.g., mineral kinetics, Monod kinetics, linear kinetics). Contains the fundamental rate constant
!>          that controls reaction speed.
!>
!> Key concepts:
!> - Rate constant k determines reaction velocity: r = k * f(c) where f(c) is concentration-dependent term
!> - Units of k depend on reaction order: [1/T] for 1st order, [L³/M/T] for 2nd order, etc.
!> - Abstract type allows polymorphism: different kinetic laws extend this base class
!>
!> Physical interpretation:
!> - Larger rate constant → faster reaction (shorter characteristic time)
!> - Temperature dependence typically Arrhenius: k = A * exp(-Ea/RT)
!> - Rate constant connects thermodynamics (equilibrium) to kinetics (time evolution)
!>
!> Derived types:
!> - kin_mineral_params_c: mineral dissolution/precipitation kinetics
!> - Monod_params_c: microbial Monod kinetics
!> - Linear kinetic reaction parameters for simple rate laws
module kin_params_m
    implicit none !< Enforce explicit variable declarations for type safety
    private !< Encapsulate module contents by default
    !> \class kin_params_c
    !> \brief Abstract base class for kinetic reaction parameters
    !> \details Provides fundamental rate constant attribute and setter method for all kinetic reactions.
    !>          Must be extended by concrete kinetic parameter types (cannot be instantiated directly).
    !>
    !> Rate law general form: r = k * f(c, T, ...) where:
    !> - r: reaction rate [M/L³/T]
    !> - k: rate constant [units vary with reaction order]
    !> - f: kinetic law function (e.g., mass action, Monod, TST)
    !> - c: concentrations [M/L³]
    !> - T: temperature [K]
    !>
    !> Common kinetic laws:
    !> - First-order: r = k*c, k [1/T]
    !> - Second-order: r = k*c₁*c₂, k [L³/M/T]
    !> - Monod: r = k*c/(K_s + c), k [M/L³/T]
    !> - TST (minerals): r = k*(1 - Ω), k [mol/m²/s]
    !>
    !> Extends to:
    !> - kin_mineral_params_c: adds surface area, Arrhenius parameters
    !> - Monod_params_c: adds half-saturation constant, inhibition terms
    type, public, abstract :: kin_params_c !< Abstract kinetic parameters superclass (must be extended, cannot instantiate)
        !> \var rate_cst Reaction rate constant
        !> Controls reaction velocity in kinetic rate law
        !> Default: 1.0 [units depend on reaction order]
        !> Physical meaning: intrinsic reaction speed at reference conditions
        !> - Larger value → faster reaction (shorter characteristic time τ = 1/k)
        !> - Temperature dependence: k(T) = k₀ * exp(-Ea/(RT)) (Arrhenius)
        !> - Determined experimentally or from literature
        real(kind=8) :: rate_cst=1d0 !< [variable units] reaction rate constant (default 1.0)
    contains
        procedure :: set_rate_cst !< Setter method to assign reaction rate constant
    end type

contains

    !> \brief Sets the reaction rate constant
    !> \details Assigns a new value to the rate constant for kinetic reactions. The rate constant determines
    !>          the intrinsic speed of the reaction at reference conditions (typically 25°C). Units depend on
    !>          the reaction order and kinetic law type.
    !>
    !> Rate constant units by reaction order:
    !> - Zero-order: [M/L³/T] - constant rate independent of concentration
    !> - First-order: [1/T] - rate proportional to concentration (e.g., radioactive decay)
    !> - Second-order: [L³/M/T] - rate proportional to product of two concentrations
    !> - Monod: [M/L³/T] - maximum rate (V_max)
    !> - Mineral (TST): [mol/m²/s] - intrinsic rate per unit surface area
    !>
    !> Physical considerations:
    !> - Must be positive (reactions proceed in defined direction)
    !> - Temperature dependence handled separately (Arrhenius or other)
    !> - Calibrated from experimental data or literature values
    !> - Uncertainty in k propagates directly to predicted concentrations
    !>
    !> \param[inout] this Kinetic parameters object to modify
    !> \param[in] rate_cst New reaction rate constant value [units vary by reaction order]
    subroutine set_rate_cst(this,rate_cst)
        implicit none !< Enforce explicit variable declarations
        class(kin_params_c), intent(inout) :: this !< Kinetic parameters object whose rate constant will be set
        real(kind=8), intent(in) :: rate_cst !< New reaction rate constant value [units depend on reaction order]
        
        !> Assign input rate constant to object attribute
        !> Direct assignment - no validation (assumes user provides physically reasonable value)
        this%rate_cst=rate_cst !< [variable units] set rate_cst attribute to input value
    end subroutine set_rate_cst
    
end module