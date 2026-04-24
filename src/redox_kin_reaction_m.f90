!> @file redox_kin_reaction_m.f90
!> @brief Redox kinetic reaction module using Monod kinetics for microbial reactions
!> @details This module defines the redox_kin_c class for modeling microbially-mediated redox reactions
!> using Monod (Michaelis-Menten) kinetics. These reactions describe electron transfer processes
!> catalyzed by microorganisms, where the rate depends on substrate availability and inhibitors.
!>
!> @par Monod Rate Law:
!> rk = k_max * X * [Π(c_i/(K_i + c_i))] * [K_inh/(K_inh + c_inh)]
!> where:
!>   - k_max = maximum specific rate constant [1/T]
!>   - X = biomass concentration [M/L³]
!>   - c_i = substrate concentrations (electron donor, acceptor) [M/L³]
!>   - K_i = half-saturation constants [M/L³]
!>   - K_inh, c_inh = inhibition constant and inhibitor concentration [M/L³]
!>
!> @see kin_reaction_m Base kinetic reaction class
!> @see Monod_params_m Monod kinetic parameters (rate constants, half-saturation constants)
!> @author Jordi Petchamé-Guerrero
!> @date November 2025

module redox_kin_reaction_m
    use kin_reaction_m, only: kin_reaction_c !< Import base kinetic reaction class and output subroutine
    use reaction_m, only: write_reaction_sup !< Import write procedure for base reaction information
    use Monod_params_m, only: Monod_params_c !< Import Monod kinetic parameters class
    use aq_phase_m, only : aq_phase_c !< Import aqueous phase and species classes
    implicit none !< Enforce explicit variable declarations for type safety
    save !< Preserve module variables between calls
    private !< Private module scope: internal details hidden from outside modules
    !> @brief Redox kinetic reaction class with Monod kinetics for microbial reactions
    !> @details Extends the base kinetic reaction class to include Monod parameters and biomass yield.
    !> Represents microbially-catalyzed redox reactions where reaction rate depends on:
    !>   - Electron donor concentration (substrate)
    !>   - Electron acceptor concentration (terminal electron acceptor, e.g., O2, NO3-, Fe3+)
    !>   - Inhibitor concentration (e.g., O2 inhibits denitrification)
    !>   - Biomass concentration (microbial catalyst)
    type, public, extends(kin_reaction_c) :: redox_kin_c
        type(Monod_params_c) :: params !< Monod kinetic parameters: rate constant k_max, half-saturation constants K_M, inhibition constant K_inh [-]
        real(kind=8) :: yield=0d0 !< Biomass yield factor Y [M_biomass/M_substrate] for cell synthesis (TODO: consider making subclass) [-]
    contains
    !> @name Setter Methods
    !> @{
        procedure :: set_Monod_params !< Assign Monod kinetic parameters to the reaction
        procedure :: set_yield !< Set biomass yield factor for microbial growth
    !> @}
    
    !> @name Computation Methods
    !> @{
        procedure :: compute_rk_drk_dc_Monod !< Compute reaction rate and Jacobian simultaneously using Monod kinetics
        procedure :: compute_drk_dc_Monod !< Compute Jacobian matrix ∂rk/∂c for Monod kinetics (given rate)
        procedure :: compute_rk_Monod !< Compute reaction rate using Monod kinetics
        procedure :: compute_biomass_growth !< Compute microbial growth rate based on substrate utilization
        procedure :: compute_rk_biomass !< Compute biomass reaction rate from growth rate and yield
    !> @}
    
    !> @name Index Management Methods
    !> @{
        procedure :: rearrange_indices_aq_phase_Monod !< Update aqueous species indices when phase composition changes
    !> @}
    
    !> @name Validation Methods
    !> @{
        !procedure :: check_indices_aq_phase_Monod !< Commented: validate aqueous species indices (not implemented)
    !> @}
    
    !> @name Output Methods
    !> @{
        procedure :: write_params=>write_Monod_params !< Write Monod parameters to file
        procedure :: write_reaction=>write_redox_reaction !< Write complete redox reaction information to file
    !> @}
    end type
    
    !> @note Legacy commented code: Alternative Monod type definition (not currently used)
    !type, public :: monod_type !< Legacy: alternative Monod type structure with linked list
    !>    PetscInt :: id !< Unique identifier for Monod term
    !>    character(len=MAXWORDLENGTH) :: species_name !< Name of species involved in Monod term
    !>    PetscReal :: half_saturation_constant !< Half-saturation constant K_M [M/L³]
    !>    PetscReal :: threshold_concentration !< Threshold concentration for reaction [M/L³]
    !>    type(monod_type), pointer :: next !< Pointer to next Monod term in linked list
    !>  end type monod_type
    
    !> @brief Interface declarations for Monod kinetics subroutines
    !> @details These interfaces define the calling signatures for computing Monod reaction rates
    !> and their derivatives. The actual implementations are in separate source files.
    interface
     
        !> @brief Compute Monod reaction rate
        !> @details Calculates the reaction rate using Monod kinetics with multiple limiting substrates
        !> and inhibition terms. Rate follows multiplicative Monod form.
        !> @param[in] this Redox kinetic reaction object containing Monod parameters
        !> @param[in] conc Concentrations array [M/L³]: [biomass, inhibitor(s), electron acceptor, electron donor]
        !> @param[out] rk Monod reaction rate [M/(L³·T)]
        subroutine compute_rk_Monod(this,conc,rk)
            import redox_kin_c !< Import redox_kin_c type definition
            implicit none !< Enforce explicit declarations
            class(redox_kin_c), intent(in) :: this !< Redox kinetic reaction object [-]
            real(kind=8), intent(in) :: conc(:) !< Concentrations of species involved in Monod rate law [M/L³]
            real(kind=8), intent(out) :: rk !< Computed Monod reaction rate [M/(L³·T)]
        end subroutine

        
        !> @brief Compute Monod reaction rate and Jacobian simultaneously
        !> @details Calculates both the Monod reaction rate and its derivatives with respect to all
        !> species concentrations in a single call (more efficient than separate calls).
        !> The Jacobian ∂rk/∂c is computed analytically using chain rule.
        !> @param[in] this Redox kinetic reaction object containing Monod parameters
        !> @param[in] conc Concentrations array [M/L³]: [biomass, inhibitor(s), electron acceptor, electron donor]
        !> @param[out] rk Monod reaction rate [M/(L³·T)]
        !> @param[out] drk_dc Jacobian array ∂rk/∂c_i [1/T] for all species i
        subroutine compute_rk_drk_dc_Monod(this,conc,rk,drk_dc)
            import redox_kin_c !< Import redox_kin_c type definition
            implicit none !< Enforce explicit declarations
            class(redox_kin_c), intent(in) :: this !< Redox kinetic reaction object [-]
            real(kind=8), intent(in) :: conc(:) !< Concentrations of species involved in Monod rate law [M/L³]
            real(kind=8), intent(out) :: rk !< Computed Monod reaction rate [M/(L³·T)]
            real(kind=8), intent(out) :: drk_dc(:) !< Jacobian derivatives ∂rk/∂c for each species [1/T]
        end subroutine
        
        !> @brief Compute Jacobian for Monod kinetics given reaction rate
        !> @details Calculates the derivatives ∂rk/∂c using the precomputed reaction rate.
        !> This is useful when the rate has already been computed and only the Jacobian is needed.
        !> @param[in] this Redox kinetic reaction object containing Monod parameters
        !> @param[in] conc Concentrations array [M/L³]: [biomass, inhibitor(s), electron acceptor, electron donor]
        !> @param[in] rk Precomputed Monod reaction rate [M/(L³·T)]
        !> @param[out] drk_dc Jacobian array ∂rk/∂c_i [1/T] for all species i
        subroutine compute_drk_dc_Monod(this,conc,rk,drk_dc)
            import redox_kin_c !< Import redox_kin_c type definition
            implicit none !< Enforce explicit declarations
            class(redox_kin_c), intent(in) :: this !< Redox kinetic reaction object [-]
            real(kind=8), intent(in) :: conc(:) !< Concentrations of species involved in Monod rate law [M/L³]
            real(kind=8), intent(in) :: rk !< Precomputed Monod reaction rate [M/(L³·T)]
            real(kind=8), intent(out) :: drk_dc(:) !< Jacobian derivatives ∂rk/∂c for each species [1/T]
        end subroutine
        
      
    end interface !< End of interface declarations
    
    contains !< Begin module subroutine implementations

        !> @brief Set Monod kinetic parameters for redox reaction
        !> @details Assigns the Monod parameters (rate constant, half-saturation constants, inhibition
        !> constants) to the redox kinetic reaction object. These parameters define the functional form
        !> of the rate law.
        !> @param[in,out] this Redox kinetic reaction object to be configured
        !> @param[in] Monod_params Monod parameters object containing k_max, K_M, K_inh values
        subroutine set_Monod_params(this,Monod_params)
            implicit none !< Enforce explicit declarations
            class(redox_kin_c), intent(inout) :: this !< Redox kinetic reaction object to configure [-]
            class(Monod_params_c), intent(in) :: Monod_params !< Monod parameters to assign (rate constant, half-saturation constants) [-]
            this%params=Monod_params !< Copy all Monod parameters to reaction object (structure assignment) [-]
        end subroutine
        
        
        
        !> @brief Rearrange aqueous phase species indices after phase composition change
        !> @details Updates the indices of species involved in Monod kinetics when the aqueous phase
        !> composition changes (e.g., after chemical system reorganization). Maps each species from
        !> the old aqueous phase to its new index in the updated phase.
        !> @param[in,out] this Redox kinetic reaction object with indices to update
        !> @param[in] aq_phase_old Previous aqueous phase configuration
        !> @param[in] aq_phase_new New aqueous phase configuration
        !> @note Stops execution if any species (inhibitor/electron acceptor/electron donor) is not found in new phase
        subroutine rearrange_indices_aq_phase_Monod(this,aq_phase_old,aq_phase_new)
            implicit none !< Enforce explicit declarations
            class(redox_kin_c), intent(inout) :: this !< Redox kinetic reaction object with species indices [-]
            class(aq_phase_c), intent(in) :: aq_phase_old !< Previous aqueous phase configuration (before reorganization) [-]
            class(aq_phase_c), intent(in) :: aq_phase_new !< New aqueous phase configuration (after reorganization) [-]
            
            integer(kind=4) :: i !< Loop counter over Monod terms (inhibitors, acceptors, donors) [-]
            integer(kind=4) :: ind_new !< New index of species in updated aqueous phase [-]
            !type(aq_species_c) :: DOC !< Unused variable: Dissolved Organic Carbon species (legacy, can be removed) [-]
            logical :: flag !< Flag indicating whether species was found in new aqueous phase (true=found, false=not found) [-]
            
            do i=1,this%params%num_terms !< Loop over all Monod terms (each term corresponds to one species: inhibitor, acceptor, or donor)
                !> Search for the i-th species (from old phase) in the new aqueous phase
                call aq_phase_new%is_species_in_aq_phase(aq_phase_old%aq_species(this%indices_aq_phase(i)),flag,ind_new)
                if (flag.eqv..true.) then !< Check if species was found in new phase
                    this%indices_aq_phase(i)=ind_new !< Update index to new position in aqueous phase array [-]
                else !< Species not found in new phase
                    error stop "Inhibitor/electron acceptor/electron donor is not in aqueous phase" !< Stop execution with error message
                end if
            end do !< End loop over Monod terms
        end subroutine
        
        !> @brief Write Monod parameters to file
        !> @details Outputs the inhibition constants to the specified file unit.
        !> Currently only writes k_inh array; could be extended to write all Monod parameters.
        !> @param[in] this Redox kinetic reaction object containing Monod parameters
        !> @param[in] unit File unit number for output (must be already opened)
        subroutine write_Monod_params(this,unit)
            implicit none !< Enforce explicit declarations
            class(redox_kin_c), intent(in) :: this !< Redox kinetic reaction object with parameters to write [-]
            integer(kind=4), intent(in) :: unit !< File unit number for output (opened file handle) [-]
            write(unit,*) this%params%k_inh !< Write inhibition constants K_inh array to file (list-directed format) [M/L³]
        end subroutine
        
        !> @brief Write complete redox reaction information to file
        !> @details Outputs both the base reaction information (stoichiometry, species) and the
        !> Monod-specific parameters to the specified file unit.
        !> @param[in] this Redox kinetic reaction object to write
        !> @param[in] unit File unit number for output (must be already opened)
        subroutine write_redox_reaction(this,unit)
            implicit none !< Enforce explicit declarations
            class(redox_kin_c), intent(in) :: this !< Redox kinetic reaction object to write [-]
            integer(kind=4), intent(in) :: unit !< File unit number for output (opened file handle) [-]
            call write_reaction_sup(this,unit) !< Write base reaction information (stoichiometry, species names, equilibrium constant)
            !call this%write_params(unit) !< Write Monod-specific parameters (k_inh array) to file
        end subroutine
        
        !> @brief Set biomass yield factor for microbial growth
        !> @details Assigns the yield factor Y which relates substrate consumption to biomass production:
        !> Y = M_biomass_produced / M_substrate_consumed [M/M]
        !> The yield must be positive (physical constraint).
        !> @param[in,out] this Redox kinetic reaction object to configure
        !> @param[in] yield Biomass yield factor Y [M_biomass/M_substrate], must be > 0
        !> @note Stops execution if yield is negative
        subroutine set_yield(this,yield)
        class(redox_kin_c), intent(inout) :: this !< Redox kinetic reaction object to configure [-]
        real(kind=8), intent(in) :: yield !< Biomass yield factor Y [M_biomass/M_substrate] for cell synthesis
        if (yield<0d0) error stop "Yield factor must be positive" !< Validate yield is non-negative (physical constraint)
        this%yield=yield !< Assign validated yield factor to reaction object [-]
        end subroutine
        
        !> @brief Compute microbial biomass growth rate
        !> @details Calculates the growth rate of biomass using Monod kinetics with biofilm space limitation.
        !> Follows equation (6) from Wang et al (WRR, 2024).
        !> Rate accounts for:
        !>   1. Available space for biofilm growth: (1 - θ_bf/φ)
        !>   2. Monod kinetics for electron donor, acceptor, and inhibitor
        !>   3. Biomass concentration (catalyst)
        !> @par Mathematical Form:
        !> rk_g = k_max * (1 - θ_bf/φ) * X * [K_inh/(K_inh + c_inh)] * c_acc/(K_acc + c_acc) * c_don/(K_don + c_don)
        !> @param[in] this Redox kinetic reaction object with Monod parameters
        !> @param[in] conc Concentrations array [M/L³]: [X (biomass), c_inh (inhibitor), c_acc (acceptor), c_don (donor)]
        !> @param[in] vol_frac_bf Biofilm volume fraction θ_bf [-] (fraction of pore space occupied by biofilm)
        !> @param[in] porosity Total porosity φ [-] (fraction of total volume that is pore space)
        !> @param[out] rk_g Biomass growth rate [M/(L³·T)]
        !> @see Wang et al., Water Resources Research, 2024, equation (6)
        subroutine compute_biomass_growth(this,conc,vol_frac_bf,porosity,rk_g)
        !< Computes the growth rate of the biomass accounting for biofilm space limitation
        !> We follow equation (6) from the paper of Wang et al (WRR, 2024)
        class(redox_kin_c), intent(in) :: this !< Redox kinetic reaction object containing Monod parameters [-]
        real(kind=8), intent(in) :: conc(:) !< Concentrations [M/L³]: conc=[X, c_inh, c_acceptor, c_donor]
        real(kind=8), intent(in) :: vol_frac_bf !< Biofilm volume fraction θ_bf [-] (0 to φ)
        real(kind=8), intent(in) :: porosity !< Porosity φ [-] (0 to 1)
        real(kind=8), intent(out) :: rk_g !< Biomass growth rate [M/(L³·T)]
        !> Compute growth rate with multiplicative Monod terms and space limitation factor
        !> Space factor: (1 - θ_bf/φ) approaches 0 as biofilm fills pores
        !> Inhibitor term: K_inh/(K_inh + c_inh) approaches 0 when inhibitor is high
        !> Donor term: c_don/(K_don + c_don) approaches 1 when donor is abundant
        !> Acceptor term: c_acc/(K_acc + c_acc) approaches 1 when acceptor is abundant
        rk_g=this%params%rate_cst*(1d0-vol_frac_bf/porosity)*conc(1)*this%params%k_inh(1)*conc(3)*conc(4)/((this%params%k_M(2)+&
            conc(4))*(this%params%k_M(1)+conc(3))*(this%params%k_inh(1)+conc(2))) !< [M/(L³·T)] = [1/T] * [-] * [M/L³] * [M/L³] * [M/L³] * [M/L³] / ([M/L³] * [M/L³] * [M/L³])
        end subroutine
        
        !> @brief Compute substrate consumption rate from biomass growth rate
        !> @details Converts biomass growth rate to substrate reaction rate using the yield factor:
        !> rk = rk_g / Y
        !> where Y = M_biomass/M_substrate.
        !> This gives the rate of substrate consumption required to support the biomass growth.
        !> @param[in] this Redox kinetic reaction object
        !> @param[in] yield Biomass yield factor Y [M_biomass/M_substrate] (TODO: may not be necessary as parameter)
        !> @param[in] rk_g Biomass growth rate [M_biomass/(L³·T)]
        !> @return rk Substrate reaction rate [M_substrate/(L³·T)] (local variable, not output argument)
        !> @note This subroutine computes rk but doesn't return it (missing intent(out)). Likely incomplete.
        subroutine compute_rk_biomass(this,yield,rk_g)
        class(redox_kin_c), intent(in) :: this !< Redox kinetic reaction object [-]
        real(kind=8), intent(in) :: yield !< Biomass yield factor Y [M_biomass/M_substrate] (TODO: may not be necessary as parameter)
        real(kind=8), intent(in) :: rk_g !< Biomass growth rate [M_biomass/(L³·T)]
        real(kind=8) :: rk !< Substrate reaction rate [M_substrate/(L³·T)] (local variable, should be intent(out))
        rk=rk_g/yield !< Convert growth rate to substrate consumption rate via yield: rk = rk_g/Y [M_substrate/(L³·T)]
        end subroutine
        
        !> @note Legacy commented code: Biomass decay subroutine (not currently implemented)
        !> @details Would compute biomass decay rate (death, maintenance, lysis)
        !subroutine compute_biomass_decay(this,yield,rk_g)
        !class(redox_kin_c) :: this !< Redox kinetic reaction object
        !real(kind=8), intent(in) :: yield !< Yield factor (TODO: may not be necessary as parameter)
        !real(kind=8), intent(in) :: rk_g !< Growth rate of the biomass
        !real(kind=8) :: rk !< Reaction rate of the biomass
        !rk=rk_g/yield !< Placeholder computation (decay typically has different formula)
        !end subroutine
        
end module !< End of redox_kin_reaction_m module