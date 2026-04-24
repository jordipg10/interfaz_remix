!> \file eq_reaction_m.f90
!> \brief Defines the equilibrium reaction subclass and related procedures.
!> \details Contains the equilibrium reaction type, public, its interface, and procedures for reading and appending equilibrium reactions,
!>          as well as computing reaction rates for binary systems in equilibrium. Equilibrium reactions reach thermodynamic equilibrium
!>          instantaneously compared to transport timescales, governed by mass action law: K_eq = ∏(a_i^ν_i) where a_i are activities
!>          and ν_i are stoichiometric coefficients.
!>
!> Key equations:
!> - Mass action law: K_eq = ∏(a_i^ν_i) for all species i
!> - Activity: a_i = γ_i * c_i where γ_i is activity coefficient [-] and c_i is concentration [M/L³]
!> - Reaction rate (binary system): r = φ * d²c₂/du² * D * (du/dx)² [M/L³/T]
!>
!> Physical interpretation:
!> - Equilibrium constant K_eq [-] quantifies reaction favorability (large K_eq → products favored)
!> - Binary system: simplest case with 2 primary species and 1 secondary species
!> - Reaction rate couples chemical equilibrium with transport (dispersion gradient effect)
module eq_reaction_m
    use reaction_m, only: reaction_c !> Base reaction type providing abstract interface for all reactions
    use aq_species_m, only: aq_species_c !> Aqueous species type for dissolved species in water
    use phase_m, only: phase_c !> Phase type for solid/gas/liquid phases
    use species_m, only: species_c !> General species type for chemical species
    implicit none !< Enforce explicit variable declarations for type safety
    save !< Preserve module variables between procedure calls
    private !< Default accessibility is private
    !> \class eq_reaction_c
    !> \brief Equilibrium reaction subclass for instantaneous reactions
    !> \details Extends reaction_c base class and contains procedures for reading and handling equilibrium reactions.
    !>          Equilibrium reactions include:
    !>          - Aqueous complexation: e.g., Ca²⁺ + CO₃²⁻ ⇌ CaCO₃⁰
    !>          - Mineral dissolution/precipitation: e.g., CaCO₃(s) ⇌ Ca²⁺ + CO₃²⁻
    !>          - Ion exchange: e.g., Na-X + K⁺ ⇌ K-X + Na⁺
    !>          - Gas dissolution: e.g., CO₂(g) + H₂O ⇌ H₂CO₃*
    !>
    !> Inherited from reaction_c:
    !> - species(:) - array of participating species (reactants and products)
    !> - stoichiometry(:) - stoichiometric coefficients (negative for reactants, positive for products)
    !> - eq_cst - equilibrium constant K_eq [-]
    !> - delta_h - enthalpy change ΔH° [kJ/mol]
    !> - react_type - reaction type index (1=aqueous, 2=mineral, 3=exchange, 6=gas)
    !>
    !> Procedures:
    !> - read_eq_reaction: generic reader for equilibrium reactions from file
    !> - read_dissolution_react_PHREEQC: parse mineral dissolution from PHREEQC database format
    !> - read_association_react_PHREEQC: parse aqueous complexation from PHREEQC database
    !> - read_exchange_react_PHREEQC: parse ion exchange reaction from PHREEQC database
    type, public, extends(reaction_c) :: eq_reaction_c !< Public equilibrium reaction type extending abstract reaction_c base class
    contains
        !procedure :: read_eq_reaction !< Generic reader for equilibrium reactions from file
        procedure :: read_dissolution_react_PHREEQC !< Parse mineral dissolution reaction from PHREEQC database format
        procedure :: read_association_react_PHREEQC !< Parse aqueous complexation reaction from PHREEQC database format
        procedure :: read_exchange_react_PHREEQC !< Parse ion exchange reaction from PHREEQC database format
    end type eq_reaction_c
    
    interface
        !> \brief Reads equilibrium reaction from file
        !> \details Generic interface for reading equilibrium reaction data from input file.
        !>          Parses reaction stoichiometry, species names, equilibrium constant, and enthalpy.
        !> \param[inout] this Equilibrium reaction object to populate
        !> \param[in] species Species object containing available species for validation
        !> \param[in] filename Input filename containing reaction definition
        subroutine read_eq_reaction(this,species,filename)
            import eq_reaction_c !< Import equilibrium reaction type from parent scope
            import species_c !< Import species type from parent scope
            implicit none !< Enforce explicit variable declarations
            class(eq_reaction_c) :: this !< Equilibrium reaction object to populate with file data
            class(species_c), intent(in) :: species !< Species object for validating species names in reaction
            character(len=*), intent(in) :: filename !< Input filename containing reaction definition (e.g., 'reactions.dat')
        end subroutine
        !> \brief Reads dissolution reaction from PHREEQC format
        !> \details Parses mineral dissolution/precipitation reaction from PHREEQC database string format.
        !>          PHREEQC format: "PHASES" block with log_k and delta_h values.
        !>          Example: "Calcite = Ca+2 + CO3-2; log_k = -8.48; -delta_H = -2.297 kcal"
        !> \param[inout] this Equilibrium reaction object to populate
        !> \param[in] string Dissolution reaction string in PHREEQC format
        !> \param[inout] phase Optional phase object to define (e.g., mineral phase)
        subroutine read_dissolution_react_PHREEQC(this,string,phase)
            import eq_reaction_c !< Import equilibrium reaction type from parent scope
            import phase_c !< Import phase type from parent scope
            implicit none !< Enforce explicit variable declarations
            class(eq_reaction_c) :: this !< Equilibrium reaction object to populate with parsed data
            character(len=*), intent(in) :: string !< Dissolution reaction string in PHREEQC format (e.g., "CaCO3 = Ca+2 + CO3-2")
            class(phase_c), intent(inout), optional :: phase !< Optional phase object to define (e.g., calcite mineral phase)
        end subroutine
        !> \brief Reads association reaction from PHREEQC format
        !> \details Parses aqueous complexation reaction from PHREEQC database string format.
        !>          PHREEQC format: "SOLUTION_SPECIES" block with log_k and delta_h values.
        !>          Example: "Ca+2 + CO3-2 = CaCO3; log_k = 3.224; delta_h = -3.545 kcal"
        !>          Primary species appear on left side; secondary (complex) species on right side.
        !> \param[inout] this Equilibrium reaction object to populate
        !> \param[in] string Association reaction string in PHREEQC format
        !> \param[out] prim_flag TRUE if species is primary, FALSE if secondary
        !> \param[out] defined_species Optional aqueous species object to define
        subroutine read_association_react_PHREEQC(this,string,prim_flag,defined_species)
            import eq_reaction_c !< Import equilibrium reaction type from parent scope
            import aq_species_c !< Import aqueous species type from parent scope
            implicit none !< Enforce explicit variable declarations
            class(eq_reaction_c) :: this !< Equilibrium reaction object to populate with parsed data
            character(len=*), intent(in) :: string !< Association reaction string in PHREEQC format (e.g., "Ca+2 + CO3-2 = CaCO3")
            logical, intent(out) :: prim_flag !< TRUE if species is primary (left side), FALSE if secondary (right side)
            type(aq_species_c), intent(out), optional :: defined_species !< Optional aqueous species object to define (the complex species)
        end subroutine
        !> \brief Reads exchange reaction from PHREEQC format
        !> \details Parses ion exchange half-reaction from PHREEQC database string format.
        !>          PHREEQC format: "EXCHANGE_SPECIES" block with log_k values.
        !>          Example: "X- + Na+ = NaX; log_k = 0.0" where X- is exchange site.
        !>          Exchange reactions describe cation/anion sorption to charged mineral surfaces.
        !> \param[inout] this Equilibrium reaction object to populate
        !> \param[in] string Half reaction string in PHREEQC format
        !> \param[out] prim_flag TRUE if species is primary exchange site
        !> \param[out] defined_species Optional species object to define (exchange complex)
        subroutine read_exchange_react_PHREEQC(this,string,prim_flag,defined_species)
            import eq_reaction_c !< Import equilibrium reaction type from parent scope
            import species_c !< Import species type from parent scope
            implicit none !< Enforce explicit variable declarations
            class(eq_reaction_c) :: this !< Equilibrium reaction object to populate with parsed data
            character(len=*), intent(in) :: string !< Half reaction string in PHREEQC format (e.g., "X- + K+ = KX")
            logical, intent(out) :: prim_flag !< TRUE if species is primary exchange site (X-), FALSE if exchange complex (KX)
            type(species_c), intent(out), optional :: defined_species !< Optional species object to define (the exchange complex, e.g., KX)
        end subroutine
    end interface
    
    contains

    !> \brief Computes the reaction rate of a binary system in equilibrium in 1D (De Simoni et al, 2005)
    !> \details Calculates the effective reaction rate for a binary equilibrium system under transport.
    !>          Binary system: 2 primary species (u₁, u₂) and 1 secondary species (c₂) in equilibrium.
    !>          For reaction: A + B ⇌ AB with K_eq = c₂/(c₁*c₂nc)
    !>          The reaction rate couples chemical equilibrium with dispersion:
    !>          r = φ * (∂²c₂/∂u²) * D * (∂u/∂x)²
    !>
    !> Mathematical derivation (De Simoni et al. 2005):
    !> - Component u combines primary and secondary species
    !> - Dispersion creates gradients that perturb local equilibrium
    !> - Second derivative (∂²c₂/∂u²) quantifies nonlinearity of equilibrium
    !> - Result: effective mixing-driven reaction rate [M/L³/T]
    !>
    !> Reference: De Simoni et al. (2005) WRR, "A procedure for the solution of multicomponent
    !>            reactive transport problems"
    !>
    !> \param[in] this Equilibrium reaction object containing K_eq
    !> \param[in] u Component concentration at target location [M/L³]
    !> \param[in] du_dx Spatial gradient of component ∂u/∂x [M/L⁴]
    !> \param[in] D Dispersion coefficient [L²/T]
    !> \param[in] phi Porosity [-]
    !> \param[out] re_mean Effective reaction rate [M/L³/T]
    subroutine react_rate_bin_syst_eq_1D(this,u,du_dx,D,phi,re_mean)
        implicit none !< Enforce explicit variable declarations
        class(eq_reaction_c), intent(in) :: this !< Equilibrium reaction object with K_eq [-]
        real(kind=8), intent(in) :: u !< Component concentration at target location [M/L³]
        real(kind=8), intent(in) :: du_dx !< Spatial gradient of component ∂u/∂x [M/L⁴]
        real(kind=8), intent(in) :: D !< Dispersion coefficient [L²/T]
        real(kind=8), intent(in) :: phi !< Porosity [-]
        real(kind=8), intent(out) :: re_mean !< Effective reaction rate [M/L³/T]
        real(kind=8) :: d2c2v_du2 !< Second derivative ∂²c₂/∂u² [L³/M] quantifying equilibrium nonlinearity
        
        !> Compute second derivative of secondary concentration with respect to component concentration
        !> Formula: ∂²c₂/∂u² = 2*K_eq / (u² + 4*K_eq)^(3/2)
        !> Derived from equilibrium constraint c₂ = K_eq*c₁*c₂nc and mass balance u = c₁ + c₂
        d2c2v_du2=2*this%eq_cst/((u**2+4*this%eq_cst)**1.5) !< [L³/M] second derivative term
        
        !> Compute effective reaction rate: r = φ * (∂²c₂/∂u²) * D * (∂u/∂x)²
        !> Physical meaning: dispersion-induced mixing drives reactions away from equilibrium
        !> Units: [-] * [L³/M] * [L²/T] * [M²/L⁸] = [M/L³/T]
        re_mean=phi*d2c2v_du2*D*du_dx**2 !< [M/L³/T] effective reaction rate from transport
    end subroutine react_rate_bin_syst_eq_1D
end module eq_reaction_m