!> \file compute_c2_from_c1_Picard.f90
!> \brief Computes secondary species concentrations from primary species using Picard iteration
!> \details This subroutine solves the chemical equilibrium speciation problem using the
!> Picard fixed-point iteration method combined with the mass action law. The method is
!> suitable for systems with both aqueous and solid primary species.
!>
!> **Mathematical Formulation**:
!>
!> The mass action law for secondary species in matrix form:
!> \f[
!> \log c_2 = S_{e,1}^* (\log \gamma_1 + \log c_1) + \log \tilde{K} - \log \gamma_2
!> \f]
!> where:
!> - \f$c_2\f$: Secondary species concentrations [M]
!> - \f$c_1\f$: Primary species concentrations [M]
!> - \f$\gamma_1, \gamma_2\f$: Activity coefficients (primary, secondary) [-]
!> - \f$S_{e,1}^*\f$: Modified stoichiometric matrix [n_sec × n_prim]
!> - \f$\tilde{K}\f$: Modified equilibrium constants [-]
!>
!> **Picard Iteration Algorithm**:
!> 1. Initialize \f$c_2^{(0)}\f$ with initial guess
!> 2. For iteration k:
!>    - Compute ionic strength: \f$I = \frac{1}{2} \sum_i c_i z_i^2\f$
!>    - Update activity coefficients: \f$\gamma_i(I)\f$ via Debye-Hückel or Davies
!>    - Compute free site concentration (for cation exchange)
!>    - Apply mass action law: \f$c_2^{(k+1)} = 10^{S_{e,1}^*(\log\gamma_1 + \log c_1) + \log\tilde{K} - \log\gamma_2}\f$
!>    - Check convergence: \f$\|\frac{c_2^{(k+1)} - c_2^{(k)}}{c_2^{(k)}}\|_\infty < \text{tol}\f$
!> 3. Exit when converged or max iterations reached
!>
!> **Special Handling**:
!> - **Cation Exchange**: Free site concentration computed as \f$X^- = 1 - \sum_{i} X_i^-\f$
!> - **Gas Phase**: If present, includes gas-aqueous equilibria (Henry's law)
!> - **Non-ideal Solutions**: Activity coefficients updated each iteration
!>
!> **Convergence Criteria**:
!> - Relative tolerance on concentration changes
!> - Maximum iteration limit to prevent infinite loops
!>
!> \note This version supports systems with aqueous + solid primary species (minerals, surface sites).
!>       For purely aqueous systems, use compute_c2_from_c1_aq_Picard instead.
!>
!> \warning Convergence is not guaranteed for all chemical systems. Poorly conditioned
!>          systems (extreme pH, high ionic strength) may require smaller tolerances or
!>          alternative methods (Newton-Raphson).
!>
!> \see compute_c2_from_c1_aq_Picard, compute_c2_from_c1_ideal, compute_c2v_from_c1_Picard
!>
!> \author Jordi Petchamé-Guerrero
!> \date 2024
!>
!> \brief Picard iteration for computing secondary species from primary species
!> \details Iteratively solves equilibrium speciation using mass action law with activity coefficient updates.
!>
!> \param[in,out] this Aqueous chemistry object containing chemical system and solution state
!> \param[in,out] c1 Primary species concentrations [M] - updated with free site concentration if cation exchange present
!> \param[in] c2_ig Initial guess for secondary variable activity species concentrations [M]
!> \param[out] c2 Converged secondary variable activity species concentrations [M] (pre-allocated)
!> \param[out] niter Number of Picard iterations performed until convergence or failure [-]
!> \param[out] CV_flag Convergence flag: TRUE if converged, FALSE if max iterations exceeded
!>
subroutine compute_c2_from_c1_Picard(this,c1,c2_ig,c2,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use vectors_m, only: inf_norm_vec_real !< Import chemistry class and linear algebra utilities

    implicit none !< Enforce explicit variable declarations
    
!> \subsection arguments Argument Declarations
    !> \param[in,out] this
    !> \brief Aqueous chemistry object
    !> \details Contains:
    !> - Chemical system definition (species, reactions, equilibrium constants)
    !> - Current solution state (concentrations, activities, ionic strength)
    !> - Speciation algebra matrices (S, U, logK)
    !> - Convergence parameters (tolerances, max iterations)
    class(aqueous_chemistry_c) :: this !< Aqueous chemistry object (polymorphic)
    
    !> \param[in,out] c1
    !> \brief Primary species concentrations [M]
    !> \details Input: Known primary species concentrations
    !> Output: Updated with computed free site concentration (last element) if cation exchange
    real(kind=8), intent(in) :: c1(:) !< Primary concentrations [M]
    
    !> \param[in] c2_ig
    !> \brief Initial guess for secondary species concentrations [M]
    !> \details Starting point for Picard iteration. Good initial guess improves convergence.
    !> Typical choices:
    !> - Previous time step solution
    !> - Ideal solution approximation (γ=1)
    !> - Small positive values (e.g., 1e-10 M)
    real(kind=8), intent(in) :: c2_ig(:) !< Initial guess secondary variable activity concentrations [M]
    
    !> \param[out] c2
    !> \brief Converged secondary variable activity species concentrations [M]
    !> \details Must be pre-allocated to size = num_sec_var_act_species
    !> Contains final concentrations after Picard iteration converges
    real(kind=8), intent(out) :: c2(:) !< Secondary variable activity concentrations [M] (pre-allocated)
    
    !> \param[out] niter
    !> \brief Number of iterations performed [-]
    !> \details Count of Picard iterations executed. Used for:
    !> - Diagnostics and performance monitoring
    !> - Identifying slow-converging cases
    integer(kind=4), intent(out) :: niter !< Number of iterations [-]
    
    !> \param[out] CV_flag
    !> \brief Convergence flag
    !> \details TRUE: Converged within tolerance and max iterations
    !> FALSE: Failed to converge (max iterations exceeded)
    logical, intent(out) :: CV_flag !< TRUE if converges, FALSE otherwise
    
!> \subsection variables Local Variable Declarations
    !> \var n_p
    !> \brief Number of primary species (aqueous + solid) [-]
    !> \details Total count of basis species including minerals and surface sites
    integer(kind=4) :: n_p !< Number of primary species
    
    !> \var n_p_aq
    !> \brief Number of aqueous primary species [-]
    !> \details Count of primary species in aqueous phase only (excludes minerals, surface sites)
    integer(kind=4) :: n_p_aq !< Number of aqueous primary species
    
    !> \var n_2_aq
    !> \brief Number of secondary aqueous species (complexes) [-]
    !> \details Count of aqueous complexes formed from primary species via equilibrium reactions
    integer(kind=4) :: n_2_aq !< Number of secondary aqueous species
    
    !> \var n_v
    !> \brief Number of variable activity species [-]
    !> \details Total count of species with non-constant activity (not used in this version)
    integer(kind=4) :: n_v !< Number of variable activity species
    
    !> \var n_e
    !> \brief Number of equilibrium reactions (total secondary species) [-]
    !> \details Total count including aqueous complexes + gas equilibria + surface complexes
    integer(kind=4) :: n_e !< Number of equilibrium reactions
    
    !> \var n_aq
    !> \brief Total number of aqueous species [-]
    !> \details n_aq = n_p_aq + n_2_aq (primary + secondary aqueous)
    integer(kind=4) :: n_aq !< Total aqueous species count
    
    !> \var i
    !> \brief Loop counter [-]
    !> \details General-purpose index for iterations
    integer(kind=4) :: i !< Loop counter
    
    !> \var k
    !> \brief Damping counter for negative free site concentration [-]
    !> \details Incremented when free site becomes negative, used to reduce c2 guess: c2_ig/(10^k)
    integer(kind=4) :: k !< Damping counter
    
    !> \var log_gamma1_old
    !> \brief Log₁₀ of primary species activity coefficients [-]
    !> \details Stores \f$\log_{10} \gamma_1\f$ from previous iteration for mass action law
    !> Dimension: [n_p]
    real(kind=8), allocatable :: log_gamma1_old(:) !< Log activity coefficients of primary species
    
    !> \var log_gamma2_old
    !> \brief Log₁₀ of secondary species activity coefficients [-]
    !> \details Stores \f$\log_{10} \gamma_2\f$ from previous iteration
    !> Dimension: [n_e] (all equilibrium reactions)
    real(kind=8), allocatable :: log_gamma2_old(:) !< Log activity coefficients of secondary species
    
    !> \var log_c2_new
    !> \brief Log₁₀ of new secondary species concentrations [-]
    !> \details Result of mass action law: \f$\log c_2^{new} = S \log(c_1 \gamma_1) + \log K - \log \gamma_2\f$
    !> Dimension: [n_e]
    real(kind=8), allocatable :: log_c2_new(:) !< Log of new secondary concentrations
    real(kind=8), allocatable :: log_c2_old(:) !< Log of old secondary concentrations
    !> \var c2_old
    !> \brief Secondary species concentrations from previous iteration [M]
    !> \details Used to check convergence: \f$\|c_2^{new} - c_2^{old}\| / \|c_2^{old}\|\f$
    !> Dimension: [n_e]
    real(kind=8), allocatable :: c2_old(:) !< Previous iteration secondary concentrations [M]
    real(kind=8), allocatable :: log_gamma_vec(:) !< Log activity coefficients for all species
    integer(kind=4) :: n_sp !< Number of species
    integer(kind=4) :: n_gas_eq !< Number of gas equilibrium reactions
    integer(kind=4) :: n_gas_eq_cst_act !< Number of gas equilibrium reactions with constant activity
    integer(kind=4) :: n_gas_eq_var_act !< Number of gas equilibrium reactions with variable activity
!****************************************************************************************************************************************************
!> \subsection preprocessing Pre-processing: Initialize dimensions and arrays
    !> Extract problem dimensions from speciation algebra object
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species     !< Total primary species
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species !< Aqueous primary species
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions      !< Equilibrium reactions
    n_2_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species !< Secondary aqueous species
    n_aq=this%aq_phase%num_species             !< Total aqueous species
    n_sp=this%solid_chemistry%reactive_zone%speciation_alg%num_species         !< Total solid species
    n_gas_eq=this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq     !< Number of gas equilibrium reactions
    n_gas_eq_cst_act=this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq     !< Number of gas equilibrium reactions
    n_gas_eq_var_act=this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq     !< Number of gas equilibrium reactions
    !> Allocate working arrays for Picard iteration
    allocate(log_gamma1_old(n_p),log_gamma_vec(n_sp),log_c2_new(n_e),c2_old(n_e),log_c2_old(n_e),log_gamma2_old(n_e))
    
    log_gamma1_old=this%get_log_gamma1()  !< Initialize primary activity coefficients to zero (ideal approximation initially)
    log_gamma2_old=this%get_log_gamma2() !< Get initial secondary activity coefficients from chemistry object
    !log_gamma_vec=this%get_log_gamma() !< Get initial secondary activity coefficients from chemistry object
    c2_old=c2_ig        !< Initialize with user-provided initial guess
    log_c2_old=log10(c2_old) !< Log of initial guess for convergence check
    niter=0             !< Initialize iteration counter to zero
    CV_flag=.false.     !< Initialize convergence flag to FALSE (not converged yet)
    !call this%compute_molalities() !< (Commented) Convert to molalities for activity calculations
    k=0                 !< Initialize damping counter to zero
    
!****************************************************************************************************************************************************
!> \subsection picard_loop Picard Iteration Loop
!> Iteratively solve for secondary species concentrations until convergence or max iterations
    do                  !< Start infinite loop with exit conditions inside
        niter=niter+1   !< Increment iteration counter
        
        !> \subsection update_ionic_strength Update ionic strength and activity coefficients
        !call this%compute_molalities() !< (Commented) Convert molarities to molalities [mol/kg H₂O]
        call this%compute_ionic_strength() !< Compute ionic strength: I = 0.5 * Σ(c_i * z_i²)
        !call this%compute_salinity()      !< (Commented) Compute salinity from TDS
        !call this%compute_molarities()    !< (Commented) Convert back to molarities
        
    !> \subsection compute_activity_coeffs Compute activity coefficients for all species
        !> Compute log₁₀(γ) for aqueous species using Debye-Hückel, Davies, or Pitzer equation
        call this%compute_log_act_coeffs_aq_chem()
        
        !> Update activities for all species: a_i = γ_i * c_i
        call this%compute_activities() !< Compute activities (aqueous + solid phases)
        
        !> Compute water activity coefficient via Raoult's law approximation
        call this%compute_log_act_coeff_wat() !< log γ_w ≈ -0.018 * Σm_i
        
        !> Extract activity coefficients for mass action law
        log_gamma1_old(1:n_p_aq)=this%log_act_coeffs(1:n_p_aq)            !< Aqueous primary γ₁
        log_gamma2_old(1:n_2_aq)=this%log_act_coeffs(n_p_aq+1:n_aq)       !< Aqueous secondary γ₂
        
        !> If gas phase present, compute gas activity coefficients (fugacity coefficients)
        ! do i=1,n_gas_eq_cst_act
        !     !call this%gas_chemistry%compute_log_act_coeffs_gas(i)        !< Compute log φ_i for gases
        !     log_gamma_vec(n_e-n_gas_eq_cst_act+i)=this%gas_chemistry%log_act_coeffs(i) !< Store gas coefficients
        ! end do
        ! do i=1,n_gas_eq_var_act
        !     !call this%gas_chemistry%compute_log_act_coeffs_gas(i)        !< Compute log φ_i for gases
        !     log_gamma_vec(n_e-n_gas_eq+i)=this%gas_chemistry%log_act_coeffs(n_gas_eq_cst_act+i) !< Store gas coefficients
        ! end do

    !> \subsection compute_free_site Compute free site concentration for cation exchange
        !> Free site concentration: X⁻ = 1 - Σ(X-M⁺) where X-M⁺ are occupied sites
        !> Last primary species (c1(n_p)) is the free cation exchange site
        ! c1(n_p)=1d0-sum(this%solid_chemistry%activities(this%solid_chemistry%num_solids-&
        !     this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl+2:this%solid_chemistry%num_solids))
                                                                            !< Free site = 1 - Σ(surface complexes)
        
        !> \subsection handle_negative_free_site Handle negative free site (non-physical situation)
        if (c1(n_p)<0d0) then       !< Free site concentration became negative (unphysical)
            k=k+1                   !< Increment damping counter
            !> Reduce initial guess by factor of 10^k to recover physically valid solution
            call this%set_conc_sec_species(c2_ig/(10**k)) !< Damped secondary concentrations
        else                        !< Free site is positive (physical solution)
            
        !> \subsection mass_action_law Apply mass action law to compute secondary concentrations
            !> Mass action law in matrix form:
            !> log c₂ = Se₁* · (log γ₁ + log c₁) + log K̃ - log γ₂
            !> Where:
            !> - Se₁* is modified stoichiometric matrix [n_e × n_p]
            !> - K̃ are modified equilibrium constants
            log_c2_new=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,log_gamma1_old+log10(c1))+&
                this%solid_chemistry%reactive_zone%speciation_alg%logK_tilde-log_gamma2_old
                                                                            !< Matrix-vector multiply + constants
            
            c2=10**log_c2_new       !< Convert from log space to linear: c₂ = 10^(log c₂)
            
        !> \subsection update_secondary Update secondary species concentrations in chemistry object
            call this%set_conc_sec_species(c2) !< Update concentrations array with new c₂ values
            !if (associated(this%gas_chemistry)) then !< (Commented) Gas phase handling
            !    !call this%gas_chemistry%update_conc_gases(c2(n_2_aq+1:n_e)*this%volume)
            !    call this%gas_chemistry%compute_vol_gas_species_conc()
            !end if
            
        !> \subsection check_convergence Check convergence criterion
             !> Relative convergence test: ||Δc₂/c₂||_∞ < tol
             if (inf_norm_vec_real((c2-c2_old)/c2_old)<this%solid_chemistry%reactive_zone%CV_params%rel_tol .and. &
                 inf_norm_vec_real((log_c2_new-log_c2_old)/log_c2_old)<this%solid_chemistry%reactive_zone%CV_params%log_rel_tol)then
                CV_flag=.true.      !< Set convergence flag to TRUE
                exit                !< Exit Picard loop - converged successfully
                
            !> \subsection check_max_iter Check for maximum iterations exceeded
            else if (niter==this%solid_chemistry%reactive_zone%CV_params%niter_max) then
                print *, inf_norm_vec_real((c2-c2_old)/c2_old) !< Print final error for diagnostics
                print *, inf_norm_vec_real((log_c2_new-log_c2_old)/log_c2_old) !< Print final error for diagnostics
                error stop "Too many Picard iterations in speciation" !< Fatal error - did not converge
                
            !> \subsection continue_iteration Continue to next iteration
            else
                c2_old=c2           !< Update c₂_old for next iteration convergence check
                log_c2_old=log_c2_new !< Update log c₂_old for next iteration
            end if
        end if                      !< End of free site positive branch
    end do                          !< End Picard iteration loop
    
!****************************************************************************************************************************************************
!> \subsection postprocessing Post-processing: Final property calculations
!> After convergence, recompute all solution properties with final concentrations
    
    !call this%compute_molalities() !< (Commented) Convert molarities to molalities for final ionic strength
    
    !> Compute final ionic strength with converged concentrations
    call this%compute_ionic_strength() !< I = 0.5 * Σ(c_i * z_i²) with final c_i
    
    !> Recompute activity coefficients with final ionic strength
    call this%compute_log_act_coeffs_aq_chem()
                                                        !< Final log γ_i using Debye-Hückel at final I
    
    !> Compute final activities for all species
    call this%compute_activities()  !< Final a_i = γ_i * c_i for aqueous and solid phases
    
    !> Legacy code for free site activity coefficient calculation (commented out)
    !> This would compute activity coefficient of free cation exchange site
    !if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
    !    if (ASSOCIATED(this%solid_chemistry%mineral_zone)) then
    !        this%solid_chemistry%activities(this%solid_chemistry%mineral_zone%num_minerals+1)=1d0-SUM(this%solid_chemistry%activities(this%solid_chemistry%mineral_zone%num_minerals+2:this%solid_chemistry%mineral_zone%num_minerals+this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl))
    !        !this%solid_chemistry%log_act_coeffs(this%solid_chemistry%mineral_zone%num_minerals+1)=LOG10(this%solid_chemistry%activities(this%solid_chemistry%mineral_zone%num_minerals+1))-LOG10(this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals+1))
    !    else
    !        this%solid_chemistry%activities(this%solid_chemistry%reactive_zone%num_minerals+1)=1d0-SUM(this%solid_chemistry%activities(this%solid_chemistry%reactive_zone%num_minerals+2:this%solid_chemistry%reactive_zone%num_solids))
    !        !this%solid_chemistry%log_act_coeffs(this%solid_chemistry%reactive_zone%num_minerals+1)=LOG10(this%solid_chemistry%activities(this%solid_chemistry%reactive_zone%num_minerals+1))-LOG10(this%solid_chemistry%concentrations(this%solid_chemistry%reactive_zone%num_minerals+1))
    !    end if
    !end if
    
    !> Compute final water activity coefficient
    call this%compute_log_act_coeff_wat() !< log γ_w via Raoult's law: γ_w ≈ exp(-0.018 * Σm_i)
    
    !> Compute final salinity (total dissolved solids)
    call this%compute_salinity()    !< S = TDS / (1 + TDS) where TDS = Σ(c_i * MW_i)
    
    !call this%compute_molarities() !< (Commented) Convert back from molalities to molarities
    
!> \subsection gas_phase_final Final gas phase calculations (if present)
    if (associated(this%gas_chemistry)) then !< Check if gas chemistry is defined
        !call this%gas_chemistry%update_conc_gases(c2v(n_v2_aq+1:n_e)*this%volume)
        !< (Commented) Update gas moles from concentrations
        
        !> Compute final gas activity coefficients (fugacity coefficients)
        call this%gas_chemistry%compute_log_act_coeffs_gases() !< log φ_i for non-ideal gas mixtures
        
        !> Compute final partial pressures from concentrations
        call this%gas_chemistry%compute_partial_pressures() !< P_i = x_i * P_total via Dalton's law
        
        !> Compute total gas pressure
        call this%gas_chemistry%compute_pressure() !< P_total = Σ P_i
        
        !> Compute total gas volume from concentrations and equation of state
        call this%gas_chemistry%compute_vol_gas_species_conc() !< V_gas via ideal or real gas law
    end if
    
 end subroutine !< End compute_c2_from_c1_Picard