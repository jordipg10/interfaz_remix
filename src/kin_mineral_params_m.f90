!> \file kin_mineral_params_m.f90
!> \brief Module for kinetic mineral reaction parameters
!>
!> \details
!> This module defines the kinetic parameters class for mineral dissolution/precipitation reactions.
!> It extends the base kinetic parameters class to include mineral-specific parameters such as:
!>
!> **Key Parameters:**
!> - Activation energy for temperature-dependent rate calculation
!> - Multiple parallel reaction mechanisms
!> - Catalytic effects from aqueous species
!> - Experimental constants (p, θ, η) for rate law formulation
!> - Supersaturation threshold for precipitation inhibition
!>
!> **Rate Law Formulation:**
!> The mineral dissolution/precipitation rate typically follows:
!> \f[
!> r = A \sum_{j=1}^{n_{par}} k_j \prod_{i=1}^{n_{cat}} a_i^{p_{ji}} (1 - \Omega^{\theta_j})^{\eta_j}
!> \f]
!> where:
!> - A = reactive surface area [m²]
!> - k_j = rate constant for parallel reaction j [mol/(m²·s)]
!> - a_i = activity of catalyser species i [-]
!> - p_ji = catalytic exponent for species i in reaction j [-]
!> - Ω = saturation index = IAP/K_eq [-]
!> - θ_j, η_j = experimental constants (typically = 1) [-]
!>
!> **Parallel Reactions:**
!> Mineral reactions can proceed via multiple parallel mechanisms (e.g., acid, neutral, base),
!> each with different rate constants and catalytic dependencies.
!>
!> **Temperature Dependence:**
!> The rate constant k varies with temperature via Arrhenius equation:
!> \f[
!> k(T) = k_0 \exp\left[-\frac{E_a}{R}\left(\frac{1}{T} - \frac{1}{T_0}\right)\right]
!> \f]
!> where E_a is the activation energy [J/mol]
!>
!> \note This module supports PHREEQC-style and PFLOTRAN-style mineral kinetics
!> \note Supersaturation threshold prevents unrealistic precipitation rates
!>
!> \see kin_params_m For base kinetic parameters class
!> \see aq_phase_m For aqueous phase and catalyser species definitions
!>
!> \author Jordi Petchamé-Guerrero
!> \date 2025

module kin_mineral_params_m
    use kin_params_m, only: kin_params_c !> Import base kinetic parameters class
    use aq_phase_m, only: aq_phase_c !> Import aqueous phase class (for catalyser species)
    implicit none !> Enforce explicit variable declaration
    save !> Save module variables between calls
    private !> Default visibility is private
!> ================================================================================
!> TYPE DEFINITION: kin_mineral_params_c
!> ================================================================================
    
    type, public, extends(kin_params_c) :: kin_mineral_params_c !> Mineral kinetic parameters subclass
        !> Extends base kin_params_c class with mineral-specific parameters
        
        !> ----------------------------------------------------------------
        !> Temperature dependence
        !> ----------------------------------------------------------------
        real(kind=8) :: act_energy !> Activation energy E_a [J/mol]
            !> Used in Arrhenius equation for temperature-dependent rate constant
            !> Typical range: 40,000 - 100,000 J/mol for mineral reactions
        
        !> ----------------------------------------------------------------
        !> Parallel reactions
        !> ----------------------------------------------------------------
        integer(kind=4) :: num_par_reacts=0 !> Number of parallel reaction mechanisms [-]
            !> e.g., acid mechanism, neutral mechanism, base mechanism
            !> Default: 0 (no parallel reactions defined)
        
        real(kind=8), allocatable :: k(:) !> Rate constants for each parallel reaction [mol/(m²·time)]
            !> Array dimension: (num_par_reacts)
            !> k(j) = rate constant for parallel reaction j at reference temperature
            !> Units depend on time unit convention of the simulation
        
        !> ----------------------------------------------------------------
        !> Catalytic effects
        !> ----------------------------------------------------------------
        integer(kind=4) :: num_cat !> Number of catalyser species [-]
            !> Catalysers are aqueous species that affect reaction rate
            !> Common catalysers: H⁺, OH⁻, Al³⁺, Fe²⁺, organic ligands
        
        integer(kind=4), allocatable :: cat_indices(:) !> Indices of catalyser species in aqueous phase [-]
            !> Array dimension: (num_cat)
            !> Maps catalysers to their position in aqueous phase species array
            !> cat_indices(i) = index of catalyser i in aq_phase%aq_species(:)
        
        real(kind=8), allocatable :: p(:,:) !> Catalytic exponents (power law coefficients) [-]
            !> Array dimension: (num_par_reacts, num_cat)
            !> p(j,i) = exponent for catalyser i in parallel reaction j
            !> Rate contribution: k_j * ∏(a_i^p_ji) where a_i = activity of catalyser i
            !> Positive p: catalysis; Negative p: inhibition
        
        !> ----------------------------------------------------------------
        !> Saturation state dependence
        !> ----------------------------------------------------------------
        real(kind=8), allocatable :: theta(:) !> Saturation state exponent θ [-]
            !> Array dimension: (num_par_reacts)
            !> theta(j) = exponent in saturation term (1 - Ω^θ)^η for reaction j
            !> Typically θ = 1 (linear saturation dependence)
            !> Controls shape of dissolution/precipitation rate vs. saturation curve
        
        real(kind=8), allocatable :: eta(:) !> Saturation state power law exponent η [-]
            !> Array dimension: (num_par_reacts)
            !> eta(j) = power law exponent for saturation term in reaction j
            !> Typically η = 1 (linear rate dependence)
            !> Rate term: (1 - Ω^θ)^η where Ω = IAP/K_eq
        
        real(kind=8) :: supersat_threshold !> Supersaturation threshold Ω_max [-]
            !> Maximum allowed saturation index for precipitation
            !> If Ω > supersat_threshold, precipitation rate may be limited or set to zero
            !> Prevents unrealistic precipitation rates at high supersaturation
            !> Typical value: 10-100 (depends on mineral and nucleation kinetics)
        real(kind=8) :: m                       !< Species travel time through solid portion, used for kinetic rate calculations (size = num_solids)
    contains
        !> Allocation procedures
        procedure :: allocate_constants !> Allocate rate constants and experimental parameters arrays
        procedure :: allocate_cat_indices !> Allocate catalyser indices array
        procedure :: set_m                                    !< Associate with spatial target object
    end type

!> ================================================================================
!> COMMENTED CODE: PFLOTRAN transition state theory formulation
!> ================================================================================
!> The following type definitions are commented out but preserved for reference.
!> They represent an alternative formulation based on PFLOTRAN's transition state
!> theory (TST) approach for mineral kinetics.
!>
!> **Key differences from current implementation:**
!> - Explicit handling of irreversible reactions
!> - Surface area dependence on volume fraction and porosity
!> - Armoring effects (surface passivation by other minerals)
!> - Affinity-based rate limiters
!> - Multiple prefactors with species-specific attenuation
!>
!> **Transition State Rate Law (PFLOTRAN):**
!> \f[
!> r = k \cdot A \cdot \prod_i a_i^{\alpha_i} \cdot f(\Delta G_r)
!> \f]
!> where f(ΔG_r) is an affinity-based function with threshold and limiter

!> PFLOTRAN:
  !type, public :: transition_state_rxn_type !> PFLOTRAN transition state reaction type
  !>  PetscReal :: min_scale_factor !> Minimum scaling factor for rate [-]
  !>  PetscReal :: affinity_factor_sigma !> Affinity factor σ parameter [-]
  !>  PetscReal :: affinity_factor_beta !> Affinity factor β parameter [-]
  !>  PetscReal :: affinity_threshold !> Affinity threshold for rate calculation [J/mol]
  !>  PetscReal :: rate_limiter !> Maximum rate limiter [mol/(m²·s)]
  !>  PetscReal :: surf_area_vol_frac_pwr !> Surface area exponent for volume fraction dependence [-]
  !>  PetscReal :: surf_area_porosity_pwr !> Surface area exponent for porosity dependence [-]
  !>  PetscInt :: irreversible !> Flag: 1 = irreversible reaction (dissolution only or precipitation only) [-]
  !>  PetscReal :: rate !> Base rate constant [mol/(m²·s)]
  !>  PetscReal :: activation_energy !> Activation energy for Arrhenius temperature dependence [J/mol]
  !>  character(len=MAXWORDLENGTH) :: armor_min_name !> Name of armoring mineral (if applicable)
  !>  PetscReal :: armor_pwr !> Armoring power law exponent [-]
  !>  PetscReal :: armor_crit_vol_frac !> Critical volume fraction for armoring effect [-]
  !>  PetscReal :: surf_area_epsilon !> Small epsilon to prevent division by zero in surface area [m²]
  !>  PetscReal :: vol_frac_epsilon !> Small epsilon to prevent division by zero in volume fraction [-]
  !>  type(transition_state_prefactor_type), pointer :: prefactor !> Pointer to prefactor parameters
  !>  type(transition_state_rxn_type), pointer :: next !> Pointer to next reaction (linked list)
  !end type transition_state_rxn_type
  !
  !type, public :: transition_state_prefactor_type !> PFLOTRAN prefactor type
  !>  type(ts_prefactor_species_type), pointer :: species !> Pointer to species-specific parameters
  !>  !> these supercede the those above in transition_state_rxn_type
  !>  PetscReal :: rate !> Prefactor-specific rate constant (overrides base rate) [mol/(m²·s)]
  !>  PetscReal :: activation_energy !> Prefactor-specific activation energy (overrides base E_a) [J/mol]
  !>  type(transition_state_prefactor_type), pointer :: next !> Pointer to next prefactor (linked list)
  !end type transition_state_prefactor_type
  !
  !type, public :: ts_prefactor_species_type !> PFLOTRAN species-specific prefactor parameters
  !>  character(len=MAXWORDLENGTH) :: name !> Species name (e.g., "H+", "Fe++")
  !>  PetscInt :: id !> Species ID/index in chemical system [-]
  !>  PetscReal :: alpha !> Power law exponent α for species activity dependence [-]
  !>  PetscReal :: beta !> Additional exponent β (alternative formulation) [-]
  !>  PetscReal :: attenuation_coef !> Attenuation coefficient for species effect [-]
  !>  type(ts_prefactor_species_type), pointer :: next !> Pointer to next species (linked list)
  !end type ts_prefactor_species_type

!> ================================================================================
!> INTERFACE DECLARATIONS
!> ================================================================================
    interface
        !> Currently no interface procedures defined
        !> Future: could include external procedures for rate calculation
    end interface
    
!> ================================================================================
!> MODULE PROCEDURES
!> ================================================================================
    contains
        
!> --------------------------------------------------------------------------------
!> Subroutine: allocate_constants
!> --------------------------------------------------------------------------------
!> \brief Allocate arrays for kinetic constants and experimental parameters
!>
!> \details
!> Allocates memory for the following arrays based on number of parallel reactions
!> and catalyser species:
!> - k(:) - Rate constants for each parallel reaction
!> - theta(:) - Saturation exponents for each parallel reaction
!> - eta(:) - Power law exponents for each parallel reaction  
!> - p(:,:) - Catalytic exponents matrix (reactions × catalysers)
!>
!> Array dimensions:
!> - k, theta, eta: (num_par_reacts)
!> - p: (num_par_reacts, num_cat)
!>
!> \param[in,out] this Mineral kinetic parameters object
!>
!> \pre num_par_reacts must be set (> 0)
!> \pre num_cat must be set (> 0)
!>
!> \note This subroutine should be called before assigning parameter values
!> \note The commented line for cat_indices allocation is superseded by
!>       allocate_cat_indices subroutine
        
        subroutine allocate_constants(this)
            implicit none !> Enforce explicit variable declaration
            class(kin_mineral_params_c) :: this !> Mineral kinetic parameters object (intent: inout)
            
            !> Allocate kinetic constants arrays
            allocate(this%k(this%num_par_reacts), &         !> Rate constants array [mol/(m²·time)]
                     this%theta(this%num_par_reacts), &      !> Saturation exponents θ [-]
                     this%eta(this%num_par_reacts))          !> Power law exponents η [-]
            
            !> (Commented out - superseded by allocate_cat_indices)
            !allocate(this%cat_indices(this%num_cat)) !> Catalyser indices array [-]
            
            !> Allocate catalytic exponents matrix
            allocate(this%p(this%num_par_reacts,this%num_cat)) !> Catalytic exponents p(j,i) [-]
                !> Dimension: (number of parallel reactions, number of catalysers)
        end subroutine
        
!> --------------------------------------------------------------------------------
!> Subroutine: allocate_cat_indices
!> --------------------------------------------------------------------------------
!> \brief Allocate array for catalyser species indices
!>
!> \details
!> Allocates memory for the cat_indices array which stores the indices of
!> catalyser species in the aqueous phase species array. This allows mapping
!> between catalyser position in the kinetic parameters and the corresponding
!> aqueous species in the chemical system.
!>
!> Array dimension: (num_cat)
!>
!> \param[in,out] this Mineral kinetic parameters object
!>
!> \pre num_cat must be set (> 0)
!>
!> \note This subroutine should be called separately from allocate_constants
!>       to allow flexibility in initialization order
        
        subroutine allocate_cat_indices(this)
            implicit none !> Enforce explicit variable declaration
            class(kin_mineral_params_c) :: this !> Mineral kinetic parameters object (intent: inout)
            
            !> Allocate catalyser indices array
            allocate(this%cat_indices(this%num_cat)) !> Catalyser species indices in aqueous phase [-]
                !> Dimension: (number of catalyser species)
                !> Will store integer indices mapping to aq_phase%aq_species(:)
        end subroutine
        
        subroutine set_m(this, m_val)
    class(kin_mineral_params_c) :: this
    real(kind=8), intent(in) :: m_val
    if (m_val <= 0d0) then
        error stop "m must be positive"
    end if
    this%m = m_val
end subroutine
end module