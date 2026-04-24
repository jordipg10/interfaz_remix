!> \file compute_conc_surf_ideal_Newton.f90
!> \brief Computes surface complex activities from cation concentrations using Newton-Raphson method for ideal ion exchange
!> \details This subroutine solves for adsorbed cation activities on exchange sites using Newton-Raphson iteration
!> in logarithmic form. It implements the ion exchange equilibrium constraint that the sum of site fractions
!> must equal unity, along with mass action equilibrium for each cation exchange half-reaction.
!>
!> The equilibrium condition for each adsorbed cation i is:
!> \f[
!> \log\left(\frac{a_i^{ads}}{c_i \cdot (1 - \sum_j a_j^{ads})}\right) = \log(K_i)
!> \f]
!> where:
!> - \f$ a_i^{ads} \f$ = activity (site fraction) of adsorbed cation i on exchange sites [-]
!> - \f$ c_i \f$ = aqueous concentration of cation i [C]
!> - \f$ K_i \f$ = equilibrium constant for cation exchange half-reaction i [-]
!> - \f$ \sum_j a_j^{ads} < 1 \f$ = constraint that total site occupation is less than unity
!>
!> Residual in logarithmic form:
!> \f[
!> r_i = -\log(c_i) - \log(1 - \sum_j a_j^{ads}) + \log(a_i^{ads}) - \log(K_i)
!> \f]
!>
!> Jacobian elements:
!> \f[
!> J_{ij} = \frac{a_j^{ads}}{1 - \sum_k a_k^{ads}} \quad (i \neq j)
!> \f]
!> \f[
!> J_{ii} = \frac{1 - \sum_k a_k^{ads} + a_i^{ads}}{1 - \sum_k a_k^{ads}} = 1 + \frac{a_i^{ads}}{1 - \sum_k a_k^{ads}}
!> \f]
!>
!> Newton update:
!> \f[
!> \Delta(\log a^{ads}) = -\mathbf{J}^{-1} \cdot \mathbf{r}
!> \f]
!> \f[
!> a_i^{ads,(n+1)} = a_i^{ads,(n)} \cdot \left[\exp(\Delta \log a_i^{ads}) - 1\right] + a_i^{ads,(n)}
!> \f]
!>
!> Assumptions:
!> - Ideal activity coefficients for aqueous cations: γᵢ = 1
!> - Ion exchange half-reactions for each adsorbed cation
!> - Initial guess for adsorbed activities provided
!>
!> \param[in,out] this Solid chemistry object containing reactive zone, cation exchange zone, and equilibrium reactions
!> \param[in] conc_cats Aqueous cation concentrations (dimension = number of cation exchange half-reactions) [C]
!> \param[in] act_ads_cats_ig Initial guess for adsorbed cation activities (site fractions) (must be already allocated) [-]
!> \param[out] niter Number of Newton iterations performed [-]
!> \param[out] CV_flag Convergence flag: TRUE if converged, FALSE otherwise [-]

subroutine compute_conc_surf_ideal_Newton(this,conc_cats,act_ads_cats_ig,niter,CV_flag)
    use solid_chemistry_m, only: solid_chemistry_c !> Import solid chemistry class and utility functions for infinity norm and linear system solver
    use vectors_m, only: inf_norm_vec_real
    use metodos_sist_lin_m, only: LU_lin_syst
    implicit none !> Enforce explicit variable declarations
    class(solid_chemistry_c) :: this !> Solid chemistry object containing reactive zone, cation exchange zone, and chemical system [-]
    real(kind=8), intent(in) :: conc_cats(:) !> Aqueous cation concentrations (dimension = number of cation exchange half-reactions) [C]
    real(kind=8), intent(in) :: act_ads_cats_ig(:) !> Initial guess for adsorbed cation activities (site fractions on exchange sites) [-]
    integer(kind=4), intent(out) :: niter !> Number of Newton iterations performed [-]
    logical, intent(out) :: CV_flag !> Convergence flag: TRUE if converged within tolerance, FALSE otherwise [-]
    
    integer(kind=4) :: n_ads_cats !> Number of adsorbed cations (number of cation exchange half-reactions, excluding reference cation) [-]
    integer(kind=4) :: n_eq !> Total number of equilibrium reactions in system [-]
    integer(kind=4) :: n_sec_aq !> Number of secondary aqueous species (not used in this routine) [-]
    integer(kind=4) :: n_v_aq !> Number of aqueous non-constant activity species (not used in this routine) [-]
    integer(kind=4) :: i !> Loop counter for adsorbed cations (rows in residual and Jacobian) [-]
    integer(kind=4) :: j !> Loop counter for adsorbed cations (columns in Jacobian) [-]
    real(kind=8), allocatable :: act_ads_cats(:) !> Current adsorbed cation activities (site fractions) at Newton iteration [-]
    real(kind=8), allocatable :: log_act_ads_cats(:) !> Logarithm of adsorbed cation activities (not actively used) [-]
    real(kind=8), allocatable :: Delta_act_ads_cats(:) !> Newton step for adsorbed cation activities in linear form [-]
    real(kind=8), allocatable :: Delta_log_act_ads_cats(:) !> Newton step for logarithm of adsorbed cation activities [-]
    real(kind=8), allocatable :: residual(:) !> Newton residual vector in logarithmic form: r = -log(c) - log(1-Σa) + log(a) - log(K) [-]
    real(kind=8), allocatable :: Jacobian(:,:) !> Jacobian matrix: ∂r/∂(log a) for Newton iteration in logarithmic variables [-]
!> Pre-processing
    !n_p_aq=this%reactive_zone%speciation_alg%num_aq_prim_species !> Commented: Number of aqueous primary species (not used)
    n_eq=this%reactive_zone%speciation_alg%num_eq_reactions !> Extract total number of equilibrium reactions from speciation algebra object
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
    do !> Begin infinite Newton iteration loop in logarithmic variables (exits when converged or max iterations exceeded)
        niter=niter+1 !> Increment Newton iteration counter
        if (niter>this%reactive_zone%CV_params%niter_max) then !> Check if maximum iterations exceeded
            print *, "Residual: ", inf_norm_vec_real(residual) !> Print residual infinity norm for diagnostics
            print *, "Too many Newton iterations in compute_conc_surf_ideal_Newton" !> Print failure message
            error stop !> Abort execution due to convergence failure in Newton iteration
        end if
        do i=1,n_ads_cats !> Loop over adsorbed cations to compute residual vector and Jacobian matrix
            residual(i)=-LOG(conc_cats(i)) - LOG(1d0-SUM(act_ads_cats)) + LOG(act_ads_cats(i)) -& !> Compute logarithmic residual part 1: -log(cᵢ) - log(1-Σaⱼ) + log(aᵢ)
             log(this%reactive_zone%chem_syst%eq_reacts(this%reactive_zone%ind_eq_reacts(n_eq-n_ads_cats+i))%eq_cst) !> Subtract log(Kᵢ) from equilibrium reaction constant to complete residual
            do j=1,n_ads_cats !> Loop over columns of Jacobian (all adsorbed cations)
                Jacobian(i,j)=act_ads_cats(j)/(1d0-SUM(act_ads_cats)) !> Compute off-diagonal Jacobian in log variables: ∂rᵢ/∂(log aⱼ) = aⱼ/(1-Σaₖ)
            end do
            Jacobian(i,i)=(1d0-SUM(act_ads_cats)+act_ads_cats(i))/(1d0-SUM(act_ads_cats)) !> Compute diagonal Jacobian: ∂rᵢ/∂(log aᵢ) = (1-Σaₖ+aᵢ)/(1-Σaₖ) = 1 + aᵢ/(1-Σaₖ)
        end do
        if (inf_norm_vec_real(residual)<this%reactive_zone%CV_params%log_rel_tol) then !> Check absolute convergence in logarithmic residual: ||r||_∞ < log_rel_tol
            CV_flag=.true. !> Set convergence flag to TRUE (Newton iteration succeeded)
            exit !> Exit Newton loop with successful convergence
        else
        !> We solve linear system: J·Δ(log a) = -r for Newton step in logarithmic variables
            call LU_lin_syst(Jacobian,-residual,this%reactive_zone%CV_params%zero,Delta_log_act_ads_cats) !> Solve for Newton step in log space using LU decomposition: J·Δ(log a) = -r
            if (inf_norm_vec_real(Delta_log_act_ads_cats/log(act_ads_cats))<this%reactive_zone%CV_params%log_rel_tol**2) then !> Check relative convergence: ||Δ(log a)/log(a)||_∞ < log_rel_tol² (temporary workaround)
                print *, "Relative error: ", inf_norm_vec_real(Delta_log_act_ads_cats/log(act_ads_cats)) !> Print relative error for diagnostics
                print *, "Newton algorithm not accurate enough in compute_conc_surf_ideal_Newton" !> Print warning message indicating convergence criterion not met
                exit !> Exit Newton loop (convergence not achieved to desired relative accuracy)
            else
                Delta_act_ads_cats=act_ads_cats*(EXP(Delta_log_act_ads_cats)-1d0) !> Convert logarithmic step to linear step: Δa = a·(exp(Δ log a) - 1)
                call this%update_act_ads_cats(act_ads_cats,Delta_act_ads_cats) !> Update adsorbed cation activities: a^(n+1) = a^(n) + Δa
                !do while (sum(act_ads_cats)>=1d0) !> Commented: Line search to ensure Σaⱼ < 1 constraint
                !    Delta_act_ads_cats=Delta_act_ads_cats/2d0 !> Commented: Halve step size if constraint violated (temporary workaround)
                !    call this%update_conc_ads_cats(act_ads_cats,Delta_act_ads_cats) !> Commented: Re-update with reduced step
                !    !error stop
                !end do
            end if
        end if
    end do !> End Newton iteration loop
!> Post-processing
    call this%set_act_surf_compl(act_ads_cats) !> Store final adsorbed cation activities in solid chemistry object
    call this%compute_conc_surf_compl() !> Compute surface complex concentrations from activities (convert site fractions to concentrations)
 end subroutine