!> \file mixing_iter_comp_ideal.f90
!> \brief Compute component and variable activity concentrations after WMA mixing with ideal solution
!>
!> \details
!> This subroutine implements the Water Mixing Approach (WMA) for equilibrium chemical
!> systems with the following simplifying assumptions:
!> - **Ideal solution behavior:** Activity coefficients γ = 1 for all species
!> - **Aqueous primary species:** All primary species are in the aqueous phase
!> - **No lumping:** Full species set used (no dimensional reduction)
!>
!> **Algorithm Overview:**
!>
!> 1. **Mixing Step:** Compute component concentrations u = U·c̃ where:
!>    - c̃ = mixed variable activity species concentrations
!>    - U = component matrix (basis transformation)
!>    - u = component concentrations
!>
!> 2. **Speciation Step:** Solve for secondary species using Newton-Raphson:
!>    - Input: Component concentrations u
!>    - Output: Variable activity concentrations c_nc (non-conservative species)
!>    - Iterative method with adaptive initialization parameter μ
!>
!> 3. **Convergence Strategy:**
!>    - If Newton fails: Increase initialization parameter μ (0 → 0.25 → 0.5 → 0.75 → 1)
!>    - μ controls interpolation: c₁^(0) = (1-μ)c₁^old + μc₁^current
!>    - If μ = 1 and still no convergence: Error (reduce time step)
!>
!> **Component Formulation:**
!> Components are linear combinations of species that are conserved in equilibrium:
!> \f[
!> u = U \cdot c
!> \f]
!> where U is the component matrix (null space of stoichiometric matrix transpose).
!>
!> **Ideal Solution Assumption:**
!> Under ideal conditions:
!> - Activity a_i = concentration c_i
!> - No activity coefficient corrections needed
!> - Simpler and faster speciation calculation
!>
!> **Difference from Other Variants:**
!> - **vs. mixing_iter_comp_exch_ideal:** No cation exchange reactions
!> - **vs. mixing_iter_comp_ideal_lump:** No lumping (full system)
!> - **vs. mixing_iter_comp:** Non-ideal version (with activity coefficients)
!>
!> \note This subroutine assumes equilibrium reactions only (no kinetics)
!> \note Time step reduction must be done externally if Newton fails
!> \note Parameters rk_tilde and mix_ratio_Rk are currently unused (for future kinetic extension)
!>
!> \see mixing_iter_comp For non-ideal version (with activity coefficients)
!> \see mixing_iter_comp_ideal_lump For lumped version
!> \see compute_c_nc_from_u_aq_Newton_ideal For Newton speciation procedure
!> \see initialise_c1_aq_iterative_method For initialization strategy
!>
!> \author jordi Petchamé-Guerrero
!> \date 2024

!> \brief Compute concentrations after WMA iteration with ideal solution and aqueous primary species
!>
!> \details
!> Water Mixing Algorithm (WMA) iteration for equilibrium chemistry with ideal
!> solution assumption. All primary species are assumed to be aqueous.
!>
!> \param[in,out] this Aqueous chemistry object at current time step
!> \param[in] c1_old Primary species concentrations at beginning of previous time step
!> \param[in] c_mix Variable activity species concentrations after mixing
!> \param[in] rk_tilde Kinetic reaction rate contributions after mixing (currently unused)
!> \param[in,out] mix_ratio_Rk Mixing ratio of kinetic reaction rate in target (currently unused)
!> \param[in] Delta_t Time step size (currently unused but kept for interface consistency)
!> \param[in] theta Time weighting factor (currently unused but kept for interface consistency)
!>                  - θ = 0: Explicit (forward Euler)
!>                  - θ = 0.5: Crank-Nicolson
!>                  - θ = 1: Implicit (backward Euler)
!> \param[out] conc_nc Variable activity (non-conservative) species concentrations (already allocated)
!>
!> \pre this%solid_chemistry%reactive_zone must be initialized
!> \pre Component matrix must be set in speciation_alg
!> \pre c_mix must have correct dimension
!> \pre conc_nc must be pre-allocated with correct dimension
!>
!> \post conc_nc contains equilibrium concentrations of variable activity species
!>
!> \note Parameters rk_tilde, mix_ratio_Rk, Delta_t, theta are reserved for future extensions

subroutine mixing_iter_comp_ideal(this,c1_old,c_hat,mix_ratio_r_old,mix_ratio_r_new,Delta_t,theta,conc_nc,conc_comp)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, initialise_iterative_method
    use vectors_m, only: inf_norm_vec_real
        !> aqueous_chemistry_c: Aqueous chemistry class with speciation methods
        !> inf_norm_vec_real: Infinity norm function for convergence checking
    implicit none !> Enforce explicit variable declaration (no implicit typing)

!> ================================================================================
!> ARGUMENTS
!> ================================================================================
    
    class(aqueous_chemistry_c) :: this !> Aqueous chemistry object (intent: inout)
        !> Contains reactive zone, speciation algorithm, and current chemical state
    
    real(kind=8), intent(in) :: c1_old(:) !> Primary species concentrations at previous time step
        !> Array dimension: (num_prim_species)
        !> Used for Newton initialization: c₁^(0) = (1-μ)c₁^old + μc₁^current
        !> Assumption: All primary species are aqueous
    
    real(kind=8), intent(in) :: c_hat(:) !> Variable activity species concentrations after mixing
        !> Array dimension: (num_species)
        !> Result of spatial mixing before equilibrium speciation
        !> Used to compute component concentrations: u = U·c̃
    
    !real(kind=8), intent(in) :: rk_tilde(:) !> Kinetic reaction rate contributions after mixing (currently unused)
        !> Array dimension: (num_kin_reactions)
        !> Reserved for future extension with kinetic reactions
        !> Would contribute to source/sink terms in reactive transport
    
    real(kind=8), intent(in) :: mix_ratio_r_old !> Mixing ratio of kinetic reaction rate at previous time step (currently unused) [-]
        !> Reserved for future extension with kinetic reactions
        !> Would weight kinetic contributions from different mixing zones
    
    real(kind=8), intent(in) :: mix_ratio_r_new !> Mixing ratio of kinetic reaction rate at new time step (currently unused) [-]
        !> Reserved for future extension with kinetic reactions
        !> Would weight kinetic contributions from different mixing zones
    
    real(kind=8), intent(in) :: Delta_t !> Time step size (currently unused)
        !> Kept for interface consistency with time-dependent versions
        !> Future use: adaptive time stepping based on convergence
    
    real(kind=8), intent(in) :: theta !> Time weighting factor (currently unused) [-]
        !> Kept for interface consistency
        !> θ ∈ [0,1]: 0=explicit, 0.5=Crank-Nicolson, 1=implicit
    
    !real(kind=8), intent(out) :: conc_comp(:) !> (Commented) Component concentrations
        !> Moved to local variable (allocated internally)
        !> No longer exposed as output parameter
    
    real(kind=8), intent(inout) :: conc_nc(:) !> Variable activity species concentrations
        !> Array dimension: (num_species)
        !> Output: equilibrium concentrations after speciation
        !> Must be pre-allocated by caller
    real(kind=8), intent(inout) :: conc_comp(:) !> component concentrations (already allocated)
!> ================================================================================
!> LOCAL VARIABLES
!> ================================================================================
    
    !> ----------------------------------------------------------------
    !> Newton iteration control
    !> ----------------------------------------------------------------
    integer(kind=4) :: niter !> Number of Newton iterations performed [-]
        !> Output from Newton speciation procedure
        !> Used for diagnostics and convergence monitoring
    
    integer(kind=4) :: k_div !> Counter for time step divisions (currently unused) [-]
        !> Reserved for adaptive time stepping
        !> Would track number of time step reductions
    
    integer(kind=4) :: k !> Counter for time steps (currently unused) [-]
        !> Reserved for multi-step integration
        !> Currently incremented but not used
    
    logical :: CV_flag !> Convergence flag from Newton iteration
        !> .true. = Newton converged to tolerance
        !> .false. = Newton failed to converge
    
    !> ----------------------------------------------------------------
    !> Newton initialization and time stepping
    !> ----------------------------------------------------------------
    real(kind=8) :: mu=0d0 !> Newton initialization parameter μ ∈ [0,1] [-]
        !> Controls initial guess: c₁^(0) = (1-μ)c₁^old + μc₁^current
        !> μ = 0: Use previous time step solution
        !> μ = 1: Use current mixed solution
        !> Adaptively increased if Newton fails: 0 → 0.25 → 0.5 → 0.75 → 1
        !> Default: 0.0 (start with previous solution)
    
    real(kind=8) :: Delta_t_bis !> Reduced time step (currently unused)
        !> Reserved for adaptive time stepping
        !> Would store subdivided time step: Δt_bis = Δt/2^k
    
    !> ----------------------------------------------------------------
    !> Concentration arrays
    !> ----------------------------------------------------------------
    real(kind=8), allocatable :: c1(:) !> Primary species concentrations at current time step
        !> Array dimension: (num_prim_species)
        !> Retrieved from chemistry object at current state
        !> All aqueous species
    
    real(kind=8), allocatable :: c1_ig(:) !> Component concentrations u
        !> Array dimension: (num_prim_species)
        !> Computed from mixing: u = U·c̃
        !> Components are conserved quantities in equilibrium

!> ================================================================================
!> PRE-PROCESSING: Setup and allocation
!> ================================================================================
    
    c1=this%get_c1() !> Get primary species concentrations at current time step
        !> Retrieves c₁ from aqueous chemistry object
        !> These are the "current" concentrations for Newton initialization
        !> All primary species are aqueous in this formulation
    
    !allocate(conc_comp(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
        !> Allocate array for component concentrations
        !> Dimension: number of primary species in speciation algorithm
        !> Will store component concentrations u = U·c̃

!> ================================================================================
!> PROCESSING: Mixing and speciation
!> ================================================================================
    
    !> ------------------------------------------------------------------------
    !> Step 1: Compute component concentrations from mixed species
    !> ------------------------------------------------------------------------
    conc_comp=MATMUL(THIS%solid_chemistry%reactive_zone%speciation_alg%comp_mat,c_hat)
        !> Matrix-vector product: u = U·c̃
        !> comp_mat = U: Component matrix (basis transformation)
        !> c_hat = c̃: Mixed variable activity species concentrations
        !> conc_comp = u: Component concentrations (conserved in equilibrium)
    
    !> ------------------------------------------------------------------------
    !> Step 2: Check if mixing changed concentrations significantly
    !> ------------------------------------------------------------------------
    if (inf_norm_vec_real(c_hat(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)-c1)< &
        this%solid_chemistry%reactive_zone%CV_params%abs_tol) then
        !> Compute infinity norm: ||c̃₁ - c₁||∞
        !> c̃₁: Primary species subset of mixed concentrations
        !> c₁: Current primary species concentrations
        !> If ||c̃₁ - c₁||∞ < ε_zero: Concentrations unchanged (within tolerance)
        
        conc_nc=this%get_conc_nc() !> No change: Use current non-conservative species concentrations
            !> Skip speciation since mixing didn't alter primary species
            !> Retrieve existing variable activity concentrations from chemistry object
    
    !> ------------------------------------------------------------------------
    !> Step 3: Significant change detected - perform speciation
    !> ------------------------------------------------------------------------
    else !> Concentrations changed: Need to recompute equilibrium
        
        !> ----------------------------------------------------------------
        !> Initialize Newton's method
        !> ----------------------------------------------------------------
        call initialise_iterative_method(c1_old,c1,mu,c1_ig) !> we compute initial guess for primary concentrations
            !> Initialize primary aqueous concentrations for Newton iteration
            !> Uses current chemistry object state and μ parameter
            !> Sets up initial guess: c₁^(0) = (1-μ)c₁^old + μc₁^current
            !> mu: Initially 0 (use previous solution)
        
        !k=0 !> (Commented) Initialize time step counter
            !> Reserved for multi-step time integration
        
        !k_div=0 !> (Commented) Initialize time step division counter
            !> Reserved for adaptive time stepping
        
        !Delta_t_bis=Delta_t !> (Commented) Initialize reduced time step
            !> Reserved for time step subdivision
        
        !> ----------------------------------------------------------------
        !> Newton iteration loop with adaptive initialization
        !> ----------------------------------------------------------------
        !> (Commented) Outer loop for time step subdivision
        !do !> Outer loop: subdivide time step if needed
            
            !> ============================================================
            !> Inner loop: Newton iteration with adaptive μ
            !> ============================================================
            do !> Infinite loop - exit when converged or error
                
                !> --------------------------------------------------------
                !> Solve for variable activity species using Newton
                !> --------------------------------------------------------
                call this%compute_c_nc_from_u_Newton_ideal(c1_ig,conc_comp,conc_nc,niter,CV_flag)
                    !> Newton-Raphson speciation for ideal aqueous solution
                    !> Input:
                    !>   conc_comp: Component concentrations (target values)
                    !> Output:
                    !>   conc_nc: Variable activity species concentrations
                    !>   niter: Number of iterations performed
                    !>   CV_flag: Convergence flag (.true. if converged)
                    !> Method: Solve for secondary species from component mass balance
                
                !> --------------------------------------------------------
                !> Check convergence and adjust if needed
                !> --------------------------------------------------------
                if (CV_flag .eqv. .false.) then !> Newton did not converge
                    
                    if (mu<1d0) then !> Haven't tried all μ values yet
                        mu=mu+0.25 !> Increase initialization parameter
                            !> Sequence: 0 → 0.25 → 0.5 → 0.75 → 1.0
                            !> Larger μ gives more weight to current mixed state
                        call initialise_iterative_method(c1_old,c1,mu,c1_ig)
                            !> Recompute initial guess with new μ
                            !> Try again with different initialization
                    
                    else !> μ = 1 and still no convergence
                        !mu=0d0 !> (Commented) Reset μ for time step reduction
                            !> Would reset initialization parameter for subdivided time step
                        
                        !k_div=k_div+1 !> (Commented) Increment division counter
                            !> Would track how many times time step was subdivided
                        
                        !Delta_t_bis=Delta_t_bis/2d0 !> (Commented) Halve time step
                            !> Would reduce time step for better convergence
                        
                        error stop "Newton speciation not converging, you must reduce time step"
                            !> Terminate with error message
                            !> User must manually reduce Δt and restart simulation
                            !> Alternative: implement automatic time step reduction
                    end if
                
                else !> Newton converged successfully
                    k=k+1 !> Increment time step counter (currently unused)
                        !> Counts successful speciation calls
                        !> Reserved for future multi-step methods
                    exit !> Exit inner loop - speciation complete
                        !> Convergence achieved, proceed to post-processing
                end if
            
            end do !> End Newton iteration loop
            
            !> (Commented) Check if all subdivided time steps completed
            !if (abs(2d0**(k_div)-k)<this%solid_chemistry%reactive_zone%CV_params%zero) then
                !> Check: Have we completed 2^k_div substeps?
                !> If yes: All subdivisions complete, exit outer loop
            !    exit !> Exit outer loop
            !end if
        !end do !> (Commented) End time step subdivision loop
        
    end if !> End concentration change check

!> ================================================================================
!> POST-PROCESSING: Cleanup and validation
!> ================================================================================
    
    deallocate(c1,c1_ig) !> Deallocate local concentration arrays
        !> Free memory for primary species concentrations c1
        !> Free memory for component concentrations conc_comp
        !> conc_nc remains allocated (output parameter, managed by caller)
    
    !> (Commented) Validation checks
    !call this%check_conc_aq_var_act_species(conc_comp) !> Check component concentrations are valid
        !> Would verify: conc_comp ≥ 0, no NaN/Inf values
        !> Disabled for performance (validation done in debug mode)
    
    !call this%check_act_aq_species() !> Check aqueous species activities are valid
        !> Would verify: activities ≥ 0, consistent with concentrations
        !> Under ideal conditions: a_i = c_i
        !> Disabled for performance (validation done in debug mode)

end subroutine