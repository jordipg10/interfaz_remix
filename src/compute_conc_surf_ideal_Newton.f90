!> \file compute_conc_surf_ideal_Newton.f90
!> \brief Computes surface complex activities from cation concentrations using Newton-Raphson method for ideal ion exchange
!> \details This subroutine solves for adsorbed cation activities on exchange sites using Newton-Raphson iteration
!> in logarithmic form. It implements the ion exchange equilibrium constraint that the sum of site fractions
!> must equal unity, along with mass action equilibrium for each cation exchange half-reaction.
!>
!> The equilibrium condition for each adsorbed cation i is:
!> \f[
!> \log_{10}\left(\frac{a_i^{ads}}{c_i \cdot (1 - \sum_j a_j^{ads})^{z_i}}\right) = \log_{10}(K_i)
!> \f]
!> where:
!> - \f$ a_i^{ads} \f$ = activity (site fraction) of adsorbed cation i on exchange sites [-]
!> - \f$ c_i \f$ = aqueous activity of cation i (equals its aqueous concentration under the ideal \f$ \gamma_i = 1 \f$ assumption) [-]
!> - \f$ K_i \f$ = equilibrium constant for cation exchange half-reaction i [-]
!> - \f$ z_i \f$ = valence of cation i = free-site (x-) stoichiometric coefficient in \f$ M_i^{z_i} + z_i X^- \rightleftharpoons M_iX_{z_i} \f$ [-]
!> - \f$ \sum_j a_j^{ads} < 1 \f$ = constraint that total site occupation is less than unity
!>
!> Residual in logarithmic form:
!> \f[
!> r_i = -\log_{10}(c_i) - z_i\log_{10}(1 - \sum_j a_j^{ads}) + \log_{10}(a_i^{ads}) - \log_{10}(K_i)
!> \f]
!>
!> Jacobian elements:
!> \f[
!> J_{ij} = \frac{z_i\,a_j^{ads}}{1 - \sum_k a_k^{ads}} \quad (i \neq j)
!> \f]
!> \f[
!> J_{ii} = \frac{1 - \sum_k a_k^{ads} + z_i\,a_i^{ads}}{1 - \sum_k a_k^{ads}} = 1 + \frac{z_i\,a_i^{ads}}{1 - \sum_k a_k^{ads}}
!> \f]
!>
!> Newton update:
!> \f[
!> \Delta(\log_{10} a^{ads}) = -\mathbf{J}^{-1} \cdot \mathbf{r}
!> \f]
!> \f[
!> a_i^{ads,(n+1)} = a_i^{ads,(n)} \cdot \left[10^{\Delta \log_{10} a_i^{ads}} - 1\right] + a_i^{ads,(n)}
!> \f]
!>
!> Assumptions:
!> - Ideal activity coefficients for aqueous cations: γᵢ = 1
!> - Ion exchange half-reactions for each adsorbed cation
!> - Initial guess for adsorbed activities provided
!>
!> \param[in,out] this Solid chemistry object containing reactive zone, cation exchange zone, and equilibrium reactions
!> \param[in] conc_cats Aqueous cation activities (equal to concentrations under the ideal γ=1 assumption; dimension = number of cation exchange half-reactions) [-]
!> \param[in] act_ads_cats_ig Initial guess for adsorbed cation activities (site fractions) (must be already allocated) [-]
!> \param[out] niter Number of Newton iterations performed [-]
!> \param[out] CV_flag Convergence flag: TRUE if converged, FALSE otherwise [-]

subroutine compute_conc_surf_ideal_Newton(this,conc_cats,act_ads_cats_ig,niter,CV_flag,logK)
    use solid_chemistry_m, only: solid_chemistry_c !> Import solid chemistry class and utility functions for infinity norm and linear system solver
    use vectors_m, only: inf_norm_vec_real
    use metodos_sist_lin_m, only: LU_lin_syst
    implicit none !> Enforce explicit variable declarations
    class(solid_chemistry_c) :: this !> Solid chemistry object containing reactive zone, cation exchange zone, and chemical system [-]
    real(kind=8), intent(in) :: conc_cats(:) !> Aqueous cation activities (equal to concentrations under the ideal γ=1 assumption; dimension = number of cation exchange half-reactions) [-]
    real(kind=8), intent(in) :: act_ads_cats_ig(:) !> Initial guess for adsorbed cation activities (site fractions on exchange sites) [-]
    integer(kind=4), intent(out) :: niter !> Number of Newton iterations performed [-]
    logical, intent(out) :: CV_flag !> Convergence flag: TRUE if converged within tolerance, FALSE otherwise [-]
    real(kind=8), intent(in) :: logK(:) !> Per-water log10 of the equilibrium constants of the EXCHANGE reactions ONLY (size = num_exch_cats = n_ads_cats), in reactive-zone order (a slice of aqueous_chemistry_c%logK), evaluated at this water's own temperature [-]
    
    integer(kind=4) :: n_ads_cats !> Number of adsorbed cations (number of cation exchange half-reactions, excluding reference cation) [-]
    integer(kind=4) :: n_eq !> Total number of equilibrium reactions in system [-]
    integer(kind=4) :: n_sec_aq !> Number of secondary aqueous species (not used in this routine) [-]
    integer(kind=4) :: n_v_aq !> Number of aqueous non-constant activity species (not used in this routine) [-]
    integer(kind=4) :: i !> Loop counter for adsorbed cations (rows in residual and Jacobian) [-]
    integer(kind=4) :: j !> Loop counter for adsorbed cations (columns in Jacobian) [-]
    integer(kind=4) :: z_cat !> Valence of exchange cation i = free-site (x-) stoichiometric coefficient in its exchange reaction [-]
    real(kind=8), allocatable :: act_ads_cats(:) !> Current adsorbed cation activities (site fractions) at Newton iteration [-]
    real(kind=8), allocatable :: log_act_ads_cats(:) !> Logarithm of adsorbed cation activities (not actively used) [-]
    real(kind=8), allocatable :: Delta_act_ads_cats(:) !> Newton step for adsorbed cation activities in linear form [-]
    real(kind=8), allocatable :: Delta_log_act_ads_cats(:) !> Newton step for logarithm of adsorbed cation activities [-]
    real(kind=8), allocatable :: residual(:) !> Newton residual vector in logarithmic form: r = -log(c) - log(1-Σa) + log(a) - log(K) [-]
    real(kind=8), allocatable :: Jacobian(:,:) !> Jacobian matrix: ∂r/∂(log a) for Newton iteration in logarithmic variables [-]
    real(kind=8) :: lambda_damp !> step-damping factor keeping every activity > 0 and Σa < 1 (free site > 0) [-]
!> Pre-processing
    !n_p_aq=this%reactive_zone%speciation_alg%num_aq_prim_species !> Commented: Number of aqueous primary species (not used)
    n_ads_cats=this%reactive_zone%cat_exch_zone%num_surf_compl-1 !> Compute number of adsorbed cations (excluding reference cation, hence -1)
    allocate(Delta_act_ads_cats(n_ads_cats),residual(n_ads_cats),Jacobian(n_ads_cats,n_ads_cats)) !> Allocate arrays for Newton iteration: activity step, residual, and Jacobian matrix
    allocate(Delta_log_act_ads_cats(n_ads_cats)) !> Allocate array for logarithmic Newton step
!!> Process: COMMENTED OUT - Alternative Newton formulation in linear variables (not currently used)
!    act_ads_cats=act_ads_cats_ig !> Initialize adsorbed cation activities with initial guess
!    niter=0 !> Initialize iteration counter to zero
!    CV_flag=.false. !> Initialize convergence flag to FALSE
!    do !> Begin Newton iteration loop in linear variables
!        niter=niter+1 !> Increment Newton iteration counter
!        if (niter>this%reactive_zone%CV_params%niter_max) then !> Check if maximum iterations exceeded
!            print *, "Residual: ", inf_norm_vec_real(residual) !> Print residual infinity norm for diagnostics
!            print *, "Too many Newton iterations in compute_conc_surf_ideal_Newton" !> Print failure message
!            error stop !> Abort execution due to convergence failure
!        end if
!        do i=1,n_ads_cats !> Loop over adsorbed cations to compute residual and Jacobian
!            residual(i)=act_ads_cats(i)/(conc_cats(i)*(1d0-SUM(act_ads_cats))) - this%reactive_zone%eq_reactions(n_eq-n_ads_cats+i)%eq_cst !> Compute residual: a/(c·(1-Σa)) - K
!            do j=1,n_ads_cats !> Loop over columns of Jacobian
!                Jacobian(i,j)=act_ads_cats(i)/(conc_cats(i)*(1d0-SUM(act_ads_cats))**2) !> Compute off-diagonal Jacobian: ∂rᵢ/∂aⱼ = aᵢ/(c·(1-Σa)²)
!            end do
!            Jacobian(i,i)=(1d0-SUM(act_ads_cats)+act_ads_cats(i))/(conc_cats(i)*(1d0-SUM(act_ads_cats))**2) !> Compute diagonal Jacobian: ∂rᵢ/∂aᵢ = (1-Σa+aᵢ)/(c·(1-Σa)²)
!        end do
!        if (inf_norm_vec_real(residual)<this%reactive_zone%CV_params%abs_tol) then !> Check absolute convergence: ||r||_∞ < abs_tol
!            CV_flag=.true. !> Set convergence flag to TRUE
!            exit !> Exit Newton loop with successful convergence
!        else
!        !> We solve linear system: J·Δa = -r
!            call LU_lin_syst(Jacobian,-residual,this%reactive_zone%CV_params%zero,Delta_act_ads_cats) !> Solve for Newton step in linear variables using LU decomposition
!            if (inf_norm_vec_real(Delta_act_ads_cats/act_ads_cats)<this%reactive_zone%CV_params%abs_tol**2) then !> Check relative convergence: ||Δa/a||_∞ < abs_tol² (temporary workaround)
!                print *, "Relative error: ", inf_norm_vec_real(Delta_act_ads_cats/act_ads_cats) !> Print relative error for diagnostics
!                print *, "Newton speciation not accurate enough in compute_conc_surf_ideal_Newton" !> Print warning message
!                exit !> Exit Newton loop (convergence not achieved to desired accuracy)
!            else
!                call this%update_conc_ads_cats(act_ads_cats,Delta_act_ads_cats) !> Update adsorbed cation activities: a^(n+1) = a^(n) + Δa
!            end if
!        end if
!    end do
!> Process: Newton-Raphson iteration in logarithmic variables for ion exchange equilibrium
    act_ads_cats=act_ads_cats_ig !> Initialize current adsorbed cation activities with initial guess
    !log_act_ads_cats=log(act_ads_cats) !> Commented: Initialize logarithm of activities (not actively used)
    niter=0 !> Initialize Newton iteration counter to zero
    CV_flag=.false. !> Initialize convergence flag to FALSE (not converged yet)
    do !> Newton-Raphson iteration in base-10 log activities for ion-exchange equilibrium
        niter=niter+1 !> Increment Newton iteration counter
        if (niter>this%reactive_zone%CV_params%niter_max) then !> Safety limit on iterations
            print *, "Residual: ", inf_norm_vec_real(residual) !> Print residual infinity norm for diagnostics
            error stop "Too many Newton iterations in compute_conc_surf_ideal_Newton"
        end if
        do i=1,n_ads_cats !> Build residual vector and Jacobian matrix
            z_cat=this%reactive_zone%chem_syst%aq_phase%aq_species(&
                this%reactive_zone%cat_exch_zone%exch_cat_indices(i))%valence !> valence zᵢ = free-site (x-) stoichiometric coefficient
            residual(i)=-LOG10(conc_cats(i)) - z_cat*LOG10(1d0-SUM(act_ads_cats)) + LOG10(act_ads_cats(i)) -& !> log10 residual: -log10(cᵢ) - zᵢ·log10(1-Σaⱼ) + log10(aᵢ) - log10(Kᵢ)
             logK(i)
            do j=1,n_ads_cats !> Loop over columns of Jacobian (all adsorbed cations)
                Jacobian(i,j)=z_cat*act_ads_cats(j)/(1d0-SUM(act_ads_cats)) !> ∂rᵢ/∂(log₁₀ aⱼ) = zᵢ·aⱼ/(1-Σaₖ)
            end do
            Jacobian(i,i)=(1d0-SUM(act_ads_cats)+z_cat*act_ads_cats(i))/(1d0-SUM(act_ads_cats)) !> ∂rᵢ/∂(log₁₀ aᵢ) = 1 + zᵢ·aᵢ/(1-Σaₖ)
        end do
        !> ABSOLUTE convergence: mass-action-law log-residual is negligible (||r||_∞ < abs_tol)
        if (inf_norm_vec_real(residual)<this%reactive_zone%CV_params%abs_tol) then
            CV_flag=.true. !> Converged (residual criterion)
            exit
        end if
        !> Newton step in log10 space (J·Δ(log₁₀ a) = -r), then linear activity increment
        call LU_lin_syst(Jacobian,-residual,this%reactive_zone%CV_params%zero,Delta_log_act_ads_cats)
        Delta_act_ads_cats=act_ads_cats*(1d1**Delta_log_act_ads_cats-1d0) !> Δa = a·(10^(Δ log₁₀ a) - 1); note Δa/a = 10^(Δ log₁₀ a) - 1
        !> RELATIVE convergence: Newton update is negligible relative to the activities
        !> (||Δa/a||_∞ < rel_tol). Well-conditioned (a>0), unlike dividing by log10(a).
        if (inf_norm_vec_real(Delta_act_ads_cats/act_ads_cats)<this%reactive_zone%CV_params%rel_tol) then
            act_ads_cats=act_ads_cats+Delta_act_ads_cats !> apply the (negligible) final step
            CV_flag=.true. !> Converged (relative-step criterion)
            exit
        end if
        !> Damped step: scale Δa so every adsorbed activity stays > 0 and the total adsorbed
        !> site occupancy Σa stays < 1. The free-site (x-) activity is 1 - Σa, so this keeps it
        !> > 0 and enforces the site balance a(x-) + Σ a(adsorbed) = 1 (sum of ALL surface-complex
        !> activities, free site included, equals 1). Also keeps log10(1-Σa) and log10(a) defined.
        lambda_damp=1d0
        do while (.not.(all(act_ads_cats+lambda_damp*Delta_act_ads_cats>0d0) .and. &
                        SUM(act_ads_cats+lambda_damp*Delta_act_ads_cats)<1d0))
            lambda_damp=0.5d0*lambda_damp
            if (lambda_damp<this%reactive_zone%CV_params%abs_tol) exit
        end do
        if (lambda_damp<this%reactive_zone%CV_params%abs_tol) then
            CV_flag=.true. !> step blocked by the Σa<1 constraint (exchanger essentially fully loaded); accept constraint-limited state
            exit
        end if
        act_ads_cats=act_ads_cats+lambda_damp*Delta_act_ads_cats !> a^(n+1) = a^(n) + λ·Δa (damped to satisfy 0<a and Σa<1)
    end do !> End Newton iteration loop
!> Post-processing
    call this%set_act_surf_compl(act_ads_cats) !> Store final adsorbed cation activities in solid chemistry object
    call this%compute_conc_surf_compl() !> Compute surface complex concentrations from activities (convert site fractions to concentrations)
 end subroutine