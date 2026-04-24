!> \file kin_mineral_m.f90
!> \brief Defines kinetic mineral reaction subclass for dissolution/precipitation.
!> \details This module extends the kinetic reaction base class to handle mineral dissolution and precipitation reactions
!>          with temperature-dependent kinetics, catalysis effects, and saturation state dependencies.
!>
!> Key equations:
!> - Rate law: r = k * A * f(T) * g(ΔG) * h(catalysts)
!>   where:
!>   * k: intrinsic rate constant [mol/m²/s]
!>   * A: reactive surface area [m²]
!>   * f(T): Arrhenius temperature dependence = exp(-Ea/(RT))
!>   * g(ΔG): saturation state function (e.g., 1-Ω for TST)
!>   * h(catalysts): product of catalyzer activities raised to powers
!> - Saturation index: Ω = IAP/K_eq where IAP is ion activity product
!> - Net rate: r_net = r_diss - r_prec (positive for dissolution)
!>
!> Physical interpretation:
!> - Dissolution: Ω < 1 (undersaturated), mineral dissolves into solution
!> - Precipitation: Ω > 1 (supersaturated), mineral precipitates from solution
!> - Equilibrium: Ω = 1 (saturated), no net reaction
!> - Temperature effect: higher T → faster kinetics (Arrhenius)
!> - Catalysis: H⁺, OH⁻, organic ligands can accelerate reactions
!>
!> Typical applications:
!> - Carbonate weathering: CaCO₃ + H⁺ → Ca²⁺ + HCO₃⁻
!> - Clay formation/dissolution
!> - Secondary mineral precipitation (e.g., gypsum, jarosite)
!> - Reactive transport in porous media
module kin_mineral_m
    use kin_mineral_params_m, only: kin_mineral_params_c !< Mineral kinetic parameters (rate constant, activation energy, surface area)
    use kin_reaction_m, only: kin_reaction_c !< Base kinetic reaction class and write utility
    use mineral_m, only: mineral_c !< Mineral type for phase properties
    use aq_phase_m, only: aq_phase_c !< Aqueous phase and species types
    use species_m, only: species_c
    use reaction_m, only: write_reaction_sup !< Base reaction class for general
    implicit none !< Enforce explicit variable declarations for type safety
    private !< Private module scope: internal details hidden from outside modules
    !> \class kin_mineral_c
    !> \brief Kinetic mineral dissolution/precipitation reaction type
    !> \details Extends kin_reaction_c with mineral-specific parameters and methods for computing
    !>          temperature-dependent, saturation-dependent, and catalyst-affected reaction rates.
    !>
    !> Inherited from kin_reaction_c:
    !> - species(:) - array of reactants/products (mineral + aqueous species)
    !> - stoichiometry(:) - stoichiometric coefficients
    !> - num_species - total number of species
    !> - num_aq_rk - number of aqueous species in reaction
    !> - indices_aq_phase(:) - indices of aqueous species in aqueous phase
    !>
    !> Additional attributes:
    !> - params - mineral kinetic parameters (k, Ea, surface area, catalysts)
    !>
    !> Rate law components:
    !> 1. Intrinsic rate constant k [mol/m²/s] - mineral-specific
    !> 2. Activation energy Ea [kJ/mol] - temperature sensitivity
    !> 3. Reactive surface area A [m²] - contact area for reaction
    !> 4. Saturation state Ω = IAP/K_eq - thermodynamic driving force
    !> 5. Catalyzer activities - H⁺, OH⁻, organic ligands
    !>
    !> Procedures:
    !> - set_mineral_params: assign kinetic parameters
    !> - set_indices_aq_phase_min: link aqueous species to phase
    !> - compute_rk_mineral: calculate reaction rate [mol/L/T]
    !> - compute_drk_dc_mineral: calculate Jacobian ∂r/∂c [1/T]
    !> - write_params: output parameters to file
    !> - write_reaction: output complete reaction info
    type, public, extends(kin_reaction_c) :: kin_mineral_c !< Kinetic mineral reaction extending base kinetic reaction class
        !> \var params Mineral kinetic parameters
        !> Contains:
        !> - rate constant k [mol/m²/s] - intrinsic dissolution/precipitation rate
        !> - activation energy Ea [kJ/mol] - controls temperature dependence
        !> - reactive surface area A [m²/L] - available surface for reaction
        !> - catalyzer parameters - exponents for H⁺, OH⁻, ligands
        !> - parallel reactions - multiple mechanisms (e.g., acid, neutral, base)
        type(kin_mineral_params_c) :: params !< Mineral kinetic parameters for rate computation
    contains
        !> Setter methods
        procedure :: set_mineral_params !< Assign mineral kinetic parameters
        procedure :: set_indices_aq_phase_min !< Map aqueous species indices to aqueous phase
        
        !> Computation methods
        procedure :: compute_rk_mineral !< Calculate mineral reaction rate r [mol/L/T]
        procedure :: compute_drk_dc_mineral !< Calculate reaction rate Jacobian ∂r/∂c [1/T]
        
        !> Output methods
        procedure :: write_params=>write_min_params !< Write kinetic parameters to file
        procedure :: write_reaction=>write_min_reaction !< Write complete reaction information
    end type
        
    interface
        !> \brief Computes mineral dissolution/precipitation reaction rate
        !> \details Calculates the net reaction rate for mineral dissolution or precipitation using:
        !>          r = k * A * exp(-Ea/(RT)) * f(Ω) * ∏(a_cat^n_cat)
        !>          where f(Ω) is saturation state function (e.g., 1-Ω for TST, (1-Ω^m) for nonlinear)
        !>
        !> Rate law components:
        !> - Intrinsic rate: k [mol/m²/s] from params
        !> - Surface area: A [m²/L] reactive surface
        !> - Arrhenius: exp(-Ea/(RT)) temperature dependence
        !> - Saturation: f(Ω) where Ω = IAP/K_eq
        !> - Catalysis: ∏(a_i^n_i) for H⁺, OH⁻, ligands
        !>
        !> Sign convention:
        !> - r > 0: dissolution (mineral → aqueous)
        !> - r < 0: precipitation (aqueous → mineral)
        !> - r = 0: equilibrium (Ω = 1)
        !>
        !> \param[in] this Kinetic mineral reaction object
        !> \param[in] act_cat Activities of catalyzer species (H⁺, OH⁻, etc.) [-]
        !> \param[in] saturation Saturation index Ω = IAP/K_eq [-]
        !> \param[in] react_surf Reactive surface area [m²/L]
        !> \param[in] temp Temperature [K]
        !> \param[out] rk Reaction rate [mol/L/T]
        subroutine compute_rk_mineral(this,act_cat,saturation,react_surf,temp,rk)
            import kin_mineral_c !< Import mineral reaction type from parent scope
            implicit none !< Enforce explicit variable declarations
            class(kin_mineral_c), intent(in) :: this !< Kinetic mineral reaction object with parameters
            real(kind=8), intent(in) :: act_cat(:) !< [-] activities of catalyzer species (H⁺, OH⁻, ligands)
            real(kind=8), intent(in) :: saturation !< [-] saturation index Ω = IAP/K_eq
            real(kind=8), intent(in) :: react_surf !< [m²/L] reactive surface area available for reaction
            real(kind=8), intent(in) :: temp !< [K] absolute temperature
            real(kind=8), intent(out) :: rk !< [mol/L/T] mineral reaction rate (>0 dissolution, <0 precipitation)
        end subroutine
                
        !> \brief Computes Jacobian of mineral reaction rate with respect to concentrations
        !> \details Calculates ∂r/∂c_i for all aqueous species in the reaction, needed for implicit time integration
        !>          and Newton-Raphson methods in reactive transport. Jacobian quantifies how reaction rate changes
        !>          when species concentrations change.
        !>
        !> Mathematical formulation:
        !> ∂r/∂c_i = ∂/∂c_i[k * A * exp(-Ea/(RT)) * f(Ω) * ∏(a_cat^n)]
        !>
        !> Chain rule components:
        !> 1. ∂f(Ω)/∂c_i via ∂Ω/∂c_i (saturation state sensitivity)
        !> 2. ∂(a_i)/∂c_i = γ_i + c_i*(∂γ_i/∂I)*(∂I/∂c_i) (activity coefficient effects)
        !> 3. ∂(∏a_cat^n)/∂c_i for catalyzers (if catalyzer concentrations vary)
        !>
        !> Physical interpretation:
        !> - Large |∂r/∂c_i| → rate very sensitive to species i (stiff system)
        !> - ∂r/∂c_i > 0 → increasing c_i accelerates forward reaction
        !> - ∂r/∂c_i < 0 → increasing c_i slows forward reaction
        !>
        !> Used for:
        !> - Implicit time integration (Backward Euler, BDF)
        !> - Newton-Raphson convergence acceleration
        !> - Stability analysis (eigenvalues of Jacobian)
        !>
        !> \param[in] this Kinetic mineral reaction object
        !> \param[in] conc_sp Concentrations of aqueous species in reaction [M/L³]
        !> \param[in] act_sp Activities of aqueous species in reaction [-]
        !> \param[in] log_act_coeffs_sp Logarithmic activity coefficients log₁₀(γ) [-]
        !> \param[in] act_cat Activities of catalyzer species [-]
        !> \param[in] saturation Saturation index Ω = IAP/K_eq [-]
        !> \param[in] react_surf Reactive surface area [m²/L]
        !> \param[in] temp Temperature [K]
        !> \param[out] drk_dc Jacobian vector ∂r/∂c [1/T] (must be pre-allocated)
        subroutine compute_drk_dc_mineral(this,conc_sp,act_sp,log_act_coeffs_sp,act_cat,saturation,react_surf,temp,drk_dc)
            import kin_mineral_c !< Import mineral reaction type from parent scope
            implicit none !< Enforce explicit variable declarations
            class(kin_mineral_c), intent(in) :: this !< Kinetic mineral reaction object with parameters
            real(kind=8), intent(in) :: conc_sp(:) !< [M/L³] concentrations of aqueous species in reaction
            real(kind=8), intent(in) :: act_sp(:) !< [-] activities of aqueous species in reaction
            real(kind=8), intent(in) :: log_act_coeffs_sp(:) !< [-] log₁₀(γ) activity coefficients of aqueous species
            real(kind=8), intent(in) :: act_cat(:) !< [-] activities of catalyzer species (H⁺, OH⁻, ligands)
            real(kind=8), intent(in) :: saturation !< [-] saturation index Ω = IAP/K_eq
            real(kind=8), intent(in) :: react_surf !< [m²/L] reactive surface area
            real(kind=8), intent(in) :: temp !< [K] absolute temperature
            real(kind=8), intent(out) :: drk_dc(:) !< [1/T] Jacobian ∂r/∂c_i (must be pre-allocated with size num_aq_species)
        end subroutine
                
    end interface
    
contains

    !> \brief Sets mineral kinetic parameters
    !> \details Assigns the mineral-specific kinetic parameters to the reaction object. These parameters
    !>          control the temperature dependence, surface area, catalysis effects, and intrinsic rate.
    !>
    !> Parameters include:
    !> - Rate constant k [mol/m²/s] - intrinsic dissolution/precipitation rate
    !> - Activation energy Ea [kJ/mol] - temperature sensitivity
    !> - Reactive surface area A [m²/L] - available reaction surface
    !> - Catalyzer exponents - powers for H⁺, OH⁻, ligands
    !> - Number of parallel reactions - multiple mechanisms (acid, neutral, base)
    !>
    !> \param[inout] this Kinetic mineral reaction object to configure
    !> \param[in] mineral_params Mineral kinetic parameters to assign
    subroutine set_mineral_params(this,mineral_params)
        implicit none !< Enforce explicit variable declarations
        class(kin_mineral_c) :: this !< Kinetic mineral reaction object to set
        type(kin_mineral_params_c), intent(in) :: mineral_params !< Mineral kinetic parameters (k, Ea, A, catalysts)
        
        !> Direct assignment of all mineral kinetic parameters
        !> Copies rate constant, activation energy, surface area, and catalyzer parameters
        this%params=mineral_params !< Assign complete parameter structure
    end subroutine
                       
    !> \brief Maps aqueous species indices from reaction to aqueous phase
    !> \details Links each aqueous species in the mineral reaction to its corresponding index in the
    !>          global aqueous phase array. Essential for accessing concentrations and activities during
    !>          rate calculations in reactive transport.
    !>
    !> Process:
    !> 1. Allocate index array (size = num_species - 1, excluding mineral itself)
    !> 2. Loop through aqueous species in reaction
    !> 3. Search for each species in aqueous phase
    !> 4. Store phase index for rapid access
    !> 5. Error if species not found (inconsistent system setup)
    !>
    !> Physical meaning:
    !> - Enables efficient lookup: c[indices_aq_phase[i]] gives concentration of species i
    !> - Critical for coupled reactive transport: map local reaction to global phase
    !> - Validates chemical system consistency
    !>
    !> \param[inout] this Kinetic mineral reaction object
    !> \param[in] aq_phase Aqueous phase containing all dissolved species
    subroutine set_indices_aq_phase_min(this,aq_phase,species)
        implicit none !< Enforce explicit variable declarations
        class(kin_mineral_c) :: this !< Kinetic mineral reaction object
        class(aq_phase_c), intent(in) :: aq_phase !< Aqueous phase with all dissolved species
        class(species_c), intent(in) :: species(:) !< Chemical system species array
        
        integer(kind=4) :: i !< Loop counter over species in reaction
        integer(kind=4) :: aq_species_ind !< Index of species in aqueous phase
        logical :: flag !< TRUE if species found in aqueous phase, FALSE otherwise
        integer(kind=4), allocatable :: temp_indices(:) !< Temporary array for aq phase indices
        integer(kind=4), allocatable :: temp_react_sp(:) !< Temporary array for reaction species indices
        
        !> First pass: count variable-activity aqueous species
        allocate(temp_indices(this%num_species), temp_react_sp(this%num_species))
        this%num_aq_rk=0
        do i=1,this%num_species
            call aq_phase%is_species_in_aq_phase(species(this%species_ind(i)),flag,aq_species_ind)
            if (flag .and. (.not. species(this%species_ind(i))%cst_act_flag)) then
                if (aq_species_ind < 1 .or. aq_species_ind > aq_phase%num_species) then
                    write(0,*) "ERROR in set_indices_aq_phase_min:"
                    write(0,*) "  Reaction: ", trim(this%name)
                    write(0,*) "  Species index in chem_syst: ", this%species_ind(i)
                    write(0,*) "  Species name: ", trim(species(this%species_ind(i))%name)
                    write(0,*) "  aq_species_ind=", aq_species_ind
                    write(0,*) "  aq_phase%num_species=", aq_phase%num_species
                    error stop "Error: Species in reaction not found in aqueous phase during index mapping."
                end if
                this%num_aq_rk=this%num_aq_rk+1
                temp_indices(this%num_aq_rk)=aq_species_ind
                temp_react_sp(this%num_aq_rk)=i !< reaction-species position for stoichiometry lookup
            end if
        end do
        
        call this%allocate_indices_aq_phase_kin_react(this%num_aq_rk)
        this%indices_aq_phase(1:this%num_aq_rk) = temp_indices(1:this%num_aq_rk)
        this%indices_react_species(1:this%num_aq_rk) = temp_react_sp(1:this%num_aq_rk)
        deallocate(temp_indices, temp_react_sp)
    end subroutine
       
    !> \brief Writes mineral kinetic parameters to file
    !> \details Outputs the key kinetic parameters to specified file unit for inspection, debugging,
    !>          or restart purposes. Parameters written in order: activation energy, number of parallel
    !>          reactions, and rate constant(s).
    !>
    !> Output format (plain text):
    !> - Line 1: Activation energy Ea [kJ/mol]
    !> - Line 2: Number of parallel reaction mechanisms
    !> - Line 3+: Rate constant(s) k [mol/m²/s] for each mechanism
    !>
    !> Physical interpretation:
    !> - Ea controls temperature sensitivity: larger Ea → stronger T-dependence
    !> - num_par_reacts: multiple mechanisms (e.g., acid, neutral, base dissolution)
    !> - k values: intrinsic rates at reference temperature (usually 25°C)
    !>
    !> \param[in] this Kinetic mineral reaction object
    !> \param[in] unit File unit number for output (must be open)
    subroutine write_min_params(this,unit)
        implicit none !< Enforce explicit variable declarations
        class(kin_mineral_c), intent(in) :: this !< Kinetic mineral reaction object with parameters
        integer(kind=4), intent(in) :: unit !< File unit number for output (pre-opened file)
        
        !> Write activation energy [kJ/mol]
        !> Controls temperature dependence via Arrhenius: k(T) = k₀*exp(-Ea/(RT))
        write(unit,"(10x,A30/)") "Activation energy [kJ/mol]: ", this%params%act_energy
        
        !> Write number of parallel reaction mechanisms
        !> Common: 3 mechanisms (acid, neutral, base) for mineral dissolution
        write(unit,"(10x,A30/)") "Number of parallel reactions: ", this%params%num_par_reacts
        
        !> Write rate constant(s) [mol/m²/s]
        !> If multiple mechanisms, k is array with one value per mechanism
        write(unit,"(10x,A30/)") "Rate constant(s) [mol/m2/s]: ", this%params%k
    end subroutine
       
    !> \brief Writes complete mineral reaction information to file
    !> \details Outputs both the general reaction information (species, stoichiometry) and the
    !>          mineral-specific kinetic parameters. Combines output from base class method
    !>          and mineral parameter writer.
    !>
    !> Output includes:
    !> 1. From write_reaction_sup (base class):
    !>    - Reaction equation with species names and stoichiometry
    !>    - Equilibrium constant
    !>    - Reaction type identifier
    !> 2. From write_params (mineral-specific):
    !>    - Activation energy Ea [kJ/mol]
    !>    - Number of parallel mechanisms
    !>    - Rate constant(s) k [mol/m²/s]
    !>
    !> Use cases:
    !> - Model output for verification
    !> - Debugging parameter values
    !> - Creating restart files
    !> - Documentation of chemical system
    !>
    !> \param[in] this Kinetic mineral reaction object
    !> \param[in] unit File unit number for output (must be open)
    subroutine write_min_reaction(this,unit)
        implicit none !< Enforce explicit variable declarations
        class(kin_mineral_c), intent(in) :: this !< Kinetic mineral reaction object
        integer(kind=4), intent(in) :: unit !< File unit number for output (pre-opened file)
        
        !> Write general reaction information (species, stoichiometry, K_eq)
        !> Inherited method from base kin_reaction_c class
        call write_reaction_sup(this,unit) !< Output reaction equation and equilibrium constant
        
        !> Write mineral-specific kinetic parameters (Ea, k, num_par_reacts)
        !> Calls bound procedure that writes activation energy and rate constants
        call this%write_params(unit) !< Output kinetic parameters (Ea, k values)
    end subroutine
    
end module