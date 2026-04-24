!> @file imm_zone_m.f90
!> @brief Immobile zone module for Multi-Rate Mass Transfer (MRMT) modeling
!> @details This module defines the immobile zone type and properties for MRMT reactive
!> transport simulations. Immobile zones represent stagnant or low-permeability regions
!> where solute transport is diffusion-limited, exchanging mass with the mobile zone via
!> first-order mass transfer processes.
!>
!> @par Immobile Zone Characteristics:
!> - No advection occurs in immobile zones (stagnant domains)
!> - Mass transfer with mobile zone governed by exchange rate α [1/T]
!> - Immobile porosity φ_im defines fraction of pore space in immobile domain [-]
!> - Residence time τ = 1/α characterizes exchange timescale [T]
!> - Probability p_i represents relative capacity of each immobile zone [-]
!>
!> @par MRMT Conceptual Model:
!> The Multi-Rate Mass Transfer model partitions total porosity:
!> φ_total = φ_m + Σ φ_im,i
!> where φ_m is mobile porosity and φ_im,i are immobile zone porosities.
!>
!> Mass transfer between mobile and immobile zones follows first-order kinetics:
!> ∂(φ_im,i * C_im,i)/∂t = α_i * φ_m * (C_m - C_im,i)
!>
!> @see MRMT_m Multi-Rate Mass Transfer module
!> @see mob_zone_m Mobile zone module
!> @author Generated documentation
!> @date November 2025

!> @brief Immobile zone module for MRMT transport modeling
module imm_zone_m
    implicit none !< Enforce explicit variable declarations for type safety
    save !< Preserve module variables between procedure calls
    
    !> @brief Immobile zone class for MRMT reactive transport
    !> @details Defines immobile zone properties including porosity, exchange kinetics,
    !> residence time, probability (capacity), flux, and concentration array.
    !> Each immobile zone represents a stagnant domain with first-order mass exchange
    !> with the mobile zone.
    type, public :: imm_zone_c
        real(kind=8) :: imm_por !< Immobile zone porosity φ_im (fraction of pore space in immobile domain, 0 < φ_im < 1) [-]
        real(kind=8) :: exch_rate !< Exchange rate coefficient α (first-order mass transfer rate between mobile and immobile zones) [1/T]
        real(kind=8) :: res_time !< Residence time τ = 1/α (characteristic timescale for mass exchange, larger τ means slower exchange) [T]
        real(kind=8) :: prob !< Probability p_i (relative capacity of this immobile zone, Σ p_i = 1 for all immobile zones) [-]
        real(kind=8) :: flux !< Mass flux between mobile and immobile zones: Φ_i = α_i * φ_m * (C_m - C_im,i) [M/(L³·T)]
        real(kind=8), allocatable :: conc(:) !< Concentration array in immobile zone: conc(species) [M/L³]
    contains
        procedure :: set_imm_por !< Setter method: set immobile porosity with validation
        procedure :: set_exch_rate !< Setter method: set exchange rate coefficient with validation
        procedure :: set_res_time !< Setter method: set residence time with validation
        procedure :: set_prob !< Setter method: set probability with validation
        procedure :: set_flux !< Setter method: set mass flux
        procedure :: set_conc !< Setter method: set concentration array with validation
        procedure :: compute_res_time_from_exch_rate !< Compute residence time from exchange rate
        procedure :: compute_exch_rate_from_res_time !< Compute exchange rate from residence time
    end type
    
contains
    
    !> @brief Set immobile porosity
    !> @details Sets the immobile zone porosity with physical validation.
    !> Immobile porosity φ_im represents the fraction of total pore space in the
    !> immobile (stagnant) domain. Must be positive and typically less than total porosity.
    !>
    !> @param[in,out] this Immobile zone object being configured
    !> @param[in] imm_por Immobile porosity value to assign [-]
    !>
    !> @par Physical Constraints:
    !> - Immobile porosity must be positive: φ_im > 0
    !> - Typically: 0 < φ_im < φ_total (less than total porosity)
    !> - Represents fraction of pore space in immobile domain
    !>
    !> @par MRMT Context:
    !> In Multi-Rate Mass Transfer models, total porosity is partitioned:
    !> φ_total = φ_m + Σ φ_im,i
    !> where φ_m is mobile porosity and φ_im,i are immobile zone porosities
    subroutine set_imm_por(this,imm_por)
        implicit none !< Enforce explicit variable declarations for type safety
        class(imm_zone_c) :: this !< Immobile zone object to be configured [-]
        real(kind=8), intent(in) :: imm_por !< Input immobile porosity value (must be positive) [-]
        
        !> Validate that immobile porosity is physically meaningful (positive value)
        if (imm_por<=0d0) error stop "Immobile porosity must be positive" !< Halt if non-positive porosity [-]
        this%imm_por=imm_por !< Assign validated immobile porosity value [-]
    end subroutine set_imm_por
    
    !> @brief Set exchange rate coefficient
    !> @details Sets the first-order mass transfer exchange rate α between mobile and
    !> immobile zones with physical validation. The exchange rate governs how quickly
    !> solutes equilibrate between mobile and immobile domains.
    !>
    !> @param[in,out] this Immobile zone object being configured
    !> @param[in] exch_rate Exchange rate coefficient to assign [1/T]
    !>
    !> @par Physical Constraints:
    !> - Exchange rate must be positive: α > 0
    !> - Units: [1/T] (inverse time)
    !> - Higher α means faster equilibration between zones
    !> - Lower α means slower mass transfer (more stagnant zone)
    !>
    !> @par Mass Transfer Equation:
    !> ∂(φ_im * C_im)/∂t = α * φ_m * (C_m - C_im)
    !> where C_m is mobile concentration and C_im is immobile concentration
    subroutine set_exch_rate(this,exch_rate)
        implicit none !< Enforce explicit variable declarations for type safety
        class(imm_zone_c) :: this !< Immobile zone object to be configured [-]
        real(kind=8), intent(in) :: exch_rate !< Input exchange rate coefficient (must be positive) [1/T]
        
        !> Validate that exchange rate is physically meaningful (positive value)
        if (exch_rate<=0d0) error stop "Exchange rate must be positive" !< Halt if non-positive exchange rate [-]
        this%exch_rate=exch_rate !< Assign validated exchange rate coefficient [1/T]
    end subroutine set_exch_rate
    
    !> @brief Set residence time
    !> @details Sets the characteristic residence time τ for mass exchange with validation.
    !> Residence time is the inverse of exchange rate (τ = 1/α) and represents the
    !> characteristic timescale for solute equilibration between mobile and immobile zones.
    !>
    !> @param[in,out] this Immobile zone object being configured
    !> @param[in] res_time Residence time to assign [T]
    !>
    !> @par Physical Constraints:
    !> - Residence time must be positive: τ > 0
    !> - Units: [T] (time)
    !> - Larger τ means slower equilibration (more stagnant)
    !> - Smaller τ means faster equilibration (less stagnant)
    !>
    !> @par Relationship to Exchange Rate:
    !> τ = 1/α (residence time is inverse of exchange rate)
    !> A 100-day residence time corresponds to α = 0.01 day⁻¹
    subroutine set_res_time(this,res_time)
        implicit none !< Enforce explicit variable declarations for type safety
        class(imm_zone_c) :: this !< Immobile zone object to be configured [-]
        real(kind=8), intent(in) :: res_time !< Input residence time (must be positive) [T]
        
        !> Validate that residence time is physically meaningful (positive value)
        if (res_time<=0d0) error stop "Residence time must be positive" !< Halt if non-positive residence time [-]
        this%res_time=res_time !< Assign validated residence time [T]
    end subroutine set_res_time
    
    !> @brief Set probability (capacity fraction)
    !> @details Sets the probability p_i representing the relative capacity of this immobile
    !> zone. In MRMT models, probabilities sum to 1 across all immobile zones and represent
    !> the relative contribution of each zone to total immobile storage.
    !>
    !> @param[in,out] this Immobile zone object being configured
    !> @param[in] prob Probability value to assign [-]
    !>
    !> @par Physical Constraints:
    !> - Probability must be between 0 and 1: 0 < p_i ≤ 1
    !> - For all immobile zones: Σ p_i = 1
    !> - Represents relative storage capacity of this zone
    !>
    !> @par MRMT Context:
    !> Probability weights the contribution of each immobile zone:
    !> - p_i = 0.5 means this zone provides 50% of immobile storage
    !> - Multiple zones with different (α_i, p_i) create multi-rate behavior
    subroutine set_prob(this,prob)
        implicit none !< Enforce explicit variable declarations for type safety
        class(imm_zone_c) :: this !< Immobile zone object to be configured [-]
        real(kind=8), intent(in) :: prob !< Input probability value (must be in (0,1]) [-]
        
        !> Validate that probability is physically meaningful (positive and not exceeding 1)
        if (prob<=0d0 .or. prob>1d0) error stop "Probability must be in range (0,1]" !< Halt if invalid probability [-]
        this%prob=prob !< Assign validated probability value [-]
    end subroutine set_prob
    
    !> @brief Set mass flux
    !> @details Sets the mass flux between mobile and immobile zones. Flux is computed
    !> from the first-order mass transfer equation and can be positive (mobile to immobile)
    !> or negative (immobile to mobile).
    !>
    !> @param[in,out] this Immobile zone object being configured
    !> @param[in] flux Mass flux value to assign [M/(L³·T)]
    !>
    !> @par Flux Equation:
    !> Φ_i = α_i * φ_m * (C_m - C_im,i)
    !> where:
    !> - Φ_i > 0: mass transfer from mobile to immobile zone
    !> - Φ_i < 0: mass transfer from immobile to mobile zone
    !> - Φ_i = 0: equilibrium (C_m = C_im,i)
    !>
    !> @note Flux can be positive, negative, or zero depending on concentration gradient
    subroutine set_flux(this,flux)
        implicit none !< Enforce explicit variable declarations for type safety
        class(imm_zone_c) :: this !< Immobile zone object to be configured [-]
        real(kind=8), intent(in) :: flux !< Input mass flux value (can be positive or negative) [M/(L³·T)]
        
        this%flux=flux !< Assign mass flux value (no validation, can be positive or negative) [M/(L³·T)]
    end subroutine set_flux
    
    !> @brief Set concentration array
    !> @details Sets the concentration array in the immobile zone with validation.
    !> Ensures that the provided concentration array has positive dimension and
    !> all concentration values are non-negative (physical requirement).
    !>
    !> @param[in,out] this Immobile zone object being configured
    !> @param[in] conc Concentration array to assign [M/L³]
    !>
    !> @par Array Dimensions:
    !> - Dimension 1: Number of species in immobile zone [-]
    !>
    !> @par Physical Constraints:
    !> - Array dimension must be positive (at least one species)
    !> - All concentration values must be non-negative (physical requirement)
    !> - Concentrations represent mass per volume in immobile pore space
    subroutine set_conc(this,conc)
        implicit none !< Enforce explicit variable declarations for type safety
        class(imm_zone_c) :: this !< Immobile zone object to be configured [-]
        real(kind=8), intent(in) :: conc(:) !< Input concentration array: conc(species) [M/L³]
        
        !> Validate that concentration array has positive dimension (at least one species)
        if (size(conc)<=0) then
            error stop "Concentration array must have positive dimension" !< Halt if array dimension invalid [-]
        !> Validate that all concentration values are physically meaningful (non-negative)
        else if (any(conc<0d0)) then
            error stop "Concentrations must be non-negative" !< Halt if negative concentrations found [-]
        end if
        
        !> Allocate concentration array if not already allocated
        if (.not. allocated(this%conc)) then
            allocate(this%conc(size(conc))) !< Allocate concentration array with same size as input [-]
        !> Check dimension consistency if already allocated
        else if (size(this%conc)/=size(conc)) then
            error stop "Concentration array dimension mismatch" !< Halt if dimensions don't match [-]
        end if
        
        this%conc=conc !< Assign validated concentration array [M/L³]
    end subroutine set_conc
    
    !> @brief Compute residence time from exchange rate
    !> @details Computes and sets residence time from the exchange rate using the
    !> relationship τ = 1/α. This ensures consistency between the two related parameters.
    !>
    !> @param[in,out] this Immobile zone object being configured
    !>
    !> @par Mathematical Relationship:
    !> τ = 1/α
    !> where τ is residence time [T] and α is exchange rate [1/T]
    !>
    !> @pre Exchange rate must be already set and positive
    !> @post Residence time is updated to be consistent with exchange rate
    subroutine compute_res_time_from_exch_rate(this)
        implicit none !< Enforce explicit variable declarations for type safety
        class(imm_zone_c) :: this !< Immobile zone object with exchange rate set [-]
        
        !> Validate that exchange rate has been set and is positive
        if (this%exch_rate<=0d0) error stop "Exchange rate must be positive to compute residence time" !< Halt if invalid exchange rate [-]
        
        this%res_time=1d0/this%exch_rate !< Compute residence time as inverse of exchange rate: τ = 1/α [T]
    end subroutine compute_res_time_from_exch_rate
    
    !> @brief Compute exchange rate from residence time
    !> @details Computes and sets exchange rate from the residence time using the
    !> relationship α = 1/τ. This ensures consistency between the two related parameters.
    !>
    !> @param[in,out] this Immobile zone object being configured
    !>
    !> @par Mathematical Relationship:
    !> α = 1/τ
    !> where α is exchange rate [1/T] and τ is residence time [T]
    !>
    !> @pre Residence time must be already set and positive
    !> @post Exchange rate is updated to be consistent with residence time
    subroutine compute_exch_rate_from_res_time(this)
        implicit none !< Enforce explicit variable declarations for type safety
        class(imm_zone_c) :: this !< Immobile zone object with residence time set [-]
        
        !> Validate that residence time has been set and is positive
        if (this%res_time<=0d0) error stop "Residence time must be positive to compute exchange rate" !< Halt if invalid residence time [-]
        
        this%exch_rate=1d0/this%res_time !< Compute exchange rate as inverse of residence time: α = 1/τ [1/T]
    end subroutine compute_exch_rate_from_res_time
    
    !> @par PFLOTRAN Reference Implementation:
    !> @details The following commented code shows an alternative immobile zone implementation
    !> from PFLOTRAN (Parallel FLOw and TRANsport simulator). This more complex structure
    !> includes decay reactions and multiple immobile species with selective output control.
    !> It is commented out as this module uses a simpler MRMT formulation focused on
    !> first-order mass transfer without explicit decay reactions in immobile zones.
    !>
    !> @note PFLOTRAN uses PETSc (Portable, Extensible Toolkit for Scientific Computation)
    !> data types (PetscInt, PetscBool, PetscReal) for parallel computing capabilities.
    
    !type, public :: immobile_type !< Commented: PFLOTRAN immobile zone type (alternative implementation) [-]
    !
    !>    PetscInt :: nimmobile !< Commented: Number of immobile zones in PFLOTRAN implementation [-]
    !>    PetscBool :: print_all !< Commented: Flag to print all immobile zone output (true/false) [-]
    !
    !>    type(immobile_species_type), pointer :: list !< Commented: Linked list of immobile species objects [-]
    !>    type(immobile_decay_rxn_type), pointer :: decay_rxn_list !< Commented: Linked list of decay reaction objects in immobile zones [-]
    !
    !>    !> immobile species
    !>    character(len=MAXWORDLENGTH), pointer :: names(:) !< Commented: Array of immobile species names (string identifiers) [-]
    !>    PetscBool, pointer :: print_me(:) !< Commented: Array of flags indicating which species to output (true/false per species) [-]
    !
    !>    !> decay rxn
    !>    PetscInt :: ndecay_rxn !< Commented: Number of decay reactions in immobile zones [-]
    !>    PetscInt, pointer :: decayspecid(:) !< Commented: Array of species IDs undergoing decay (index mapping) [-]
    !>    PetscReal, pointer :: decay_rate_constant(:) !< Commented: Array of first-order decay rate constants λ [1/T]
    !end type immobile_type !< Commented: End of PFLOTRAN alternative immobile zone type definition [-]
    
end module !< End of immobile zone module for MRMT modeling