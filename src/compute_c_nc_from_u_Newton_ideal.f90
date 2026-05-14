!> \file compute_c_nc_from_u_Newton_ideal.f90
!> \brief Computes variable activity species concentrations from component concentrations using Newton-Raphson method
!> \details This subroutine solves the nonlinear speciation problem for ideal solutions
!> using a robust Newton-Raphson algorithm with Levenberg-Marquardt damping, backtracking line search,
!> log-space formulation, and stagnation detection. Given component concentrations and an initial guess
!> for primary species, it computes variable activity species concentrations by iteratively solving:
!> \f[
!> \mathbf{r}(\mathbf{c}_1) = \mathbf{U}_1 \mathbf{c}_1 + \mathbf{U}_2 \cdot \mathbf{c}_{2,v} - \mathbf{u} = \mathbf{0}
!> \f]
!> where:
!> - \f$ \mathbf{c}_1 \f$ = primary species concentrations
!> - \f$ \mathbf{c}_{2,v} \f$ = secondary variable activity species concentrations
!> - \f$ \mathbf{U}_1 \f$ = component matrix for primary species
!> - \f$ \mathbf{U}_2 \f$ = component matrix for secondary variable activity species
!> - \f$ \mathbf{u} \f$ = component concentrations (input)
!>
!> Assumptions:
!> - Ideal solution behavior (ion activity coefficients γ = 1)
!> - Initial guess for primary concentrations is provided as input
!>
!> Robustness features:
!> - Log-space Newton (p = log10(c1)): removes 1/c1 singularity in Jacobian
!> - Levenberg-Marquardt damping: regularizes near-singular Jacobians
!> - Backtracking line search: ensures monotone residual decrease
!> - Stagnation detection: exits early when no progress is being made
!> - Best-solution tracking: returns the best iterate if convergence fails
!> - Perturbation retry recovery: on max-iteration or stagnation exit, restarts once from
!>   the geometric mean (log-space midpoint) of the initial guess and the best iterate:
!>   \f$ c_{1,j}^{\text{retry}} = \sqrt{\max(c_{1,j}^{\text{ig}},\,\texttt{tiny}) \cdot \max(c_{1,j}^{\text{best}},\,\texttt{tiny})} \f$,
!>   resetting all counters and damping parameters for a fresh Newton run
!>
!> Newton iteration in log-space:
!> \f[
!> \mathbf{p}^{(i+1)} = \mathbf{p}^{(i)} + \alpha \Delta\mathbf{p}^{(i)}, \quad \mathbf{p} = \log_{10}(\mathbf{c}_1)
!> \f]
!> where \f$ \Delta\mathbf{p} \f$ solves:
!> \f[
!> \left(\mathbf{J}_p + \lambda \mathbf{I}\right) \Delta\mathbf{p} = -\mathbf{r}
!> \f]
!> with \f$ \mathbf{J}_p(:,j) = \ln(10) \cdot c_{1,j} \cdot \mathbf{J}(:,j) \f$
!>
!> \param[in,out] this Aqueous chemistry class object
!> \param[in] c1_ig Initial guess for primary species concentrations
!> \param[in] conc_comp Component concentrations
!> \param[out] conc_nc Variable activity species concentrations (already allocated)
!> \param[out] niter Number of Newton iterations performed [-]
!> \param[out] CV_flag Convergence flag: TRUE if converged, FALSE otherwise [-]

subroutine compute_c_nc_from_u_Newton_ideal(this,c1_ig,conc_comp,conc_nc,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c  !< Import aqueous chemistry class for speciation methods
    use vectors_m, only: inf_norm_vec_real              !< Import infinity norm utility for convergence checks
    !use metodos_sist_lin_m, only: LU_lin_syst           !< Import LU-based linear system solver for Newton step
    implicit none                                       !< Enforce explicit variable declarations
!> Arguments
    class(aqueous_chemistry_c) :: this                 !< Aqueous chemistry class object [-]
    real(kind=8), intent(in) :: c1_ig(:)               !< Initial guess for primary species concentrations (n_p)
    real(kind=8), intent(in) :: conc_comp(:)           !< Component concentrations (total element concentrations) (n_p)
    real(kind=8), intent(out) :: conc_nc(:)            !< Variable activity species concentrations (already allocated) (n_v)
    integer(kind=4), intent(out) :: niter              !< Number of Newton iterations performed [-]
    logical, intent(out) :: CV_flag                    !< Convergence flag: .true. if converged, .false. otherwise [-]
!> Local variables for speciation algebra dimensions and Jacobian
    real(kind=8), allocatable :: dc2v_dc1(:,:)        !< Jacobian \f$ \partial \mathbf{c}_{2,v}/\partial \mathbf{c}_1 \f$ for ideal solution (n_e x n_p)
    real(kind=8), allocatable :: residual(:)           !< Newton residual \f$ \mathbf{r} = \mathbf{U}_1 \mathbf{c}_1 + \mathbf{U}_2 \mathbf{c}_{2,v} - \mathbf{u} \f$ (n_p)

    real(kind=8), allocatable :: log_gamma_var_act(:)  !< Log10 activity coefficients for variable activity species (all zero for aqueous species) (n_v)
    integer(kind=4) :: i,n_e,n_p,n_v                  !< Loop index; num equilibrium reactions; num primary species; num variable activity species [-]
!> Log-space Newton variables: \f$ \mathbf{p} = \log_{10}(\mathbf{c}_1) \f$
    real(kind=8), allocatable :: log_c1(:)             !< Log10 primary concentrations \f$ p_j = \log_{10}(c_{1,j}) \f$ [-] (n_p)
    real(kind=8), allocatable :: Delta_p(:)            !< Newton step in log-space \f$ \Delta\mathbf{p} \f$ [-] (n_p)
    real(kind=8), allocatable :: mat_lin_syst_p(:,:)   !< Log-space Jacobian \f$ J_p(:,j) = \ln(10)\,c_{1,j}\,J(:,j) \f$ [-] (n_p x n_p)
    real(kind=8), allocatable :: mat_lin_syst_LM(:,:)  !< LM-damped Jacobian \f$ \mathbf{J}_p + \lambda\mathbf{I} \f$ [-] (n_p x n_p)
    real(kind=8), parameter :: ln10 = log(10d0)        !< Natural logarithm of 10: \f$ \ln(10) \approx 2.302585 \f$ [-]
    real(kind=8), parameter :: max_dp = 10d0           !< Maximum allowed log-step magnitude per component (prevents overflow) [-]
    real(kind=8) :: scale_factor                        !< Uniform scaling factor for bounding \f$ \Delta\mathbf{p} \f$ to respect max_dp [-]
!> Levenberg-Marquardt (LM) damping variables
    real(kind=8) :: lambda_LM                          !< LM damping parameter \f$ \lambda \f$; increased on failure, decreased on success [-]
    integer(kind=4) :: i_LM                            !< LM attempt counter within a single Newton step [-]
    integer(kind=4), parameter :: max_LM_tries = 10    !< Maximum number of LM damping attempts per Newton step [-]
    logical :: sing_flag                                !< Flag: .true. if linear system is singular or step contains NaN/Inf [-]
!> Backtracking line search variables
    real(kind=8) :: alpha                              !< Step length \f$ \alpha \in (0,1] \f$, halved on each backtrack [-]
    real(kind=8) :: res_norm                           !< Current residual infinity norm \f$ \|\mathbf{r}\|_\infty \f$
    real(kind=8) :: res_trial_norm                     !< Trial residual infinity norm after candidate step
    real(kind=8), allocatable :: conc_nc_save(:)       !< Saved concentrations for line search rollback (n_v)
    real(kind=8), allocatable :: residual_trial(:)     !< Trial residual vector at candidate step (n_p)

    integer(kind=4) :: i_LS                            !< Line search backtracking counter [-]
    integer(kind=4), parameter :: max_LS_tries = 10    !< Maximum number of half-step backtracks per LM attempt [-]
    logical :: step_accepted                           !< Flag: .true. if line search found a step reducing the residual [-]
!> Stagnation detection and best-solution tracking variables
    real(kind=8), allocatable :: conc_nc_best(:)       !< Best concentrations found across all iterations (n_v)
    real(kind=8) :: best_res_norm                      !< Smallest residual norm encountered: \f$ \min_k \|\mathbf{r}^{(k)}\|_\infty \f$
    real(kind=8) :: res_norm_prev                      !< Previous iteration residual norm for stagnation comparison
    integer(kind=4) :: n_stag                          !< Counter of consecutive stagnation iterations [-]
    integer(kind=4), parameter :: n_stag_max = 20      !< Max allowed stagnation iterations before early exit [-]
    real(kind=8), parameter :: stag_rtol = 1d-3        !< Relative tolerance: stagnation if \f$ \|\mathbf{r}\| \geq (1 - \text{stag\_rtol})\|\mathbf{r}_\text{prev}\| \f$ [-]
    logical :: retry_applied                            !< Whether perturbation retry has been attempted (at most once per call) [-]
    !> Cached loop-invariant parameters (avoid repeated deep struct dereferences per iteration)
    integer(kind=4) :: niter_max                        !< Maximum Newton iterations (cached from CV_params)
    real(kind=8) :: abs_tol                             !< Absolute convergence tolerance (cached)
    real(kind=8) :: rel_tol                             !< Relative convergence tolerance (cached)
    real(kind=8) :: log_abs_tol                         !< Log absolute tolerance (cached)
    real(kind=8) :: log_rel_tol                         !< Log relative tolerance (cached)
    integer(kind=4), allocatable :: ipiv_dgesv(:)       !< Pivot index array for LAPACK dgesv (size n_p)
    integer(kind=4) :: info_dgesv                       !< Return code from LAPACK dgesv (0 = success)
    real(kind=8) :: eps_d                               !< Machine epsilon (cached)
    real(kind=8) :: sqrt_eps_d                           !< eps^(1/2) — practical conditioning noise floor for double precision
    real(kind=8) :: conc_comp_scale                     !< max(||conc_comp||_inf, 1) — loop-invariant
    real(kind=8), allocatable :: comp_mat(:,:)          !< Full component matrix (cached, n_p x n_v) for residual computation
    real(kind=8), allocatable :: U1(:,:)                !< Component matrix slice for primary species comp_mat(:,1:n_p) (cached)
    real(kind=8), allocatable :: U2(:,:)                !< Component matrix slice for secondary species comp_mat(:,n_p+1:n_v) (cached)
    logical :: skip_eval                                !< Skip residual/secondary recomputation after accepted step
    
!> Pre-Process: Initialize variables and allocate arrays
    CV_flag=.false.                                    !< Assume non-convergence until proven otherwise
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions   !< Extract number of equilibrium reactions from speciation algebra
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species   !< Extract number of primary species from speciation algebra
    n_v=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species !< Extract total number of variable activity species
    !> Cache loop-invariant parameters from deep struct dereferences
    niter_max=this%solid_chemistry%reactive_zone%CV_params%niter_max !< Cache max Newton iterations
    abs_tol=this%solid_chemistry%reactive_zone%CV_params%abs_tol !< Cache absolute tolerance
    rel_tol=this%solid_chemistry%reactive_zone%CV_params%rel_tol !< Cache relative tolerance
    log_abs_tol=this%solid_chemistry%reactive_zone%CV_params%log_abs_tol !< Cache log absolute tolerance
    log_rel_tol=this%solid_chemistry%reactive_zone%CV_params%log_rel_tol !< Cache log relative tolerance
    eps_d=epsilon(1d0)                                 !< Cache machine epsilon
    sqrt_eps_d=sqrt(eps_d)                              !< Cache eps^(1/2) ≈ 1.49e-8 (practical noise floor)
    conc_comp_scale=max(inf_norm_vec_real(conc_comp),1d0) !< Cache scaled norm for convergence checks
    !> Cache component matrix slices (avoid repeated deep dereference + slicing each iteration)
    comp_mat=this%solid_chemistry%reactive_zone%speciation_alg%comp_mat !< Cache full component matrix
    U1=comp_mat(:,1:n_p) !< Primary species slice
    U2=comp_mat(:,n_p+1:n_v) !< Secondary species slice
    !> Allocate residual vector (n_p), Jacobian (n_e x n_p), and activity coefficients (n_v)
    allocate(residual(n_p),dc2v_dc1(n_e,n_p),log_gamma_var_act(n_v))
    !> Allocate log-space Newton arrays: log concentrations, step, Jacobian matrices (n_p x n_p), and system matrix
    allocate(log_c1(n_p),Delta_p(n_p),mat_lin_syst_p(n_p,n_p),mat_lin_syst_LM(n_p,n_p),ipiv_dgesv(n_p))
    !> Allocate line search state-save arrays (n_v) and best-solution tracker (n_v)
    allocate(conc_nc_save(n_v),residual_trial(n_p),conc_nc_best(n_v))

    niter=0                                            !< Reset iteration counter
    log_gamma_var_act=this%get_log_gamma_var_act()     !< Retrieve \f$ \log_{10}(\gamma_i) \f$ for all variable activity species (zero for aqueous species)
    lambda_LM=1d-3                                     !< Initial LM damping: small value for near-Newton behavior
    n_stag=0                                           !< Reset stagnation counter
    res_norm_prev=huge(1d0)                            !< Initialize previous residual to machine-maximum so first iterate always improves
    best_res_norm=huge(1d0)                            !< Initialize best residual to machine-maximum
    retry_applied=.false.                              !< Perturbation retry not yet attempted
    skip_eval=.false.                                  !< No accepted step yet; must compute residual

!> Initialize primary species concentrations \f$ \mathbf{c}_1 \f$ with the user-supplied initial guess, clamped to 1d-30 to ensure positivity for log-space operations without underflow
    conc_nc(1:n_p)=max(c1_ig,1d-30)

!> Process: Newton-Raphson iteration loop with LM damping, line search, and log-space formulation
    newton_loop: do                                    !< Begin main Newton iteration loop
        niter=niter+1                                  !< Increment Newton iteration counter
        if (niter>niter_max) then                      !< Exceeded iteration budget?
            if (.not. retry_applied) then               !< Has perturbation retry been tried yet?
                !> Perturbation retry: restart from geometric mean of initial guess and best iterate (log-space midpoint)
                !> \f$ c_{1,j}^{\text{retry}} = \sqrt{\max(c_{1,j}^{\text{ig}},\,10^{-30}) \cdot \max(c_{1,j}^{\text{best}},\,10^{-30})} \f$
                retry_applied = .true.                  !< Mark retry as used (only one attempt allowed per call)
                do i = 1, n_p                           !< Loop over each primary species
                    conc_nc(i) = sqrt(max(c1_ig(i), 1d-30) * max(conc_nc_best(i), 1d-30))  !< Geometric mean of initial guess and best iterate
                end do                                  !< End primary species perturbation loop
                niter = 0                               !< Reset Newton iteration counter for fresh restart
                n_stag = 0                              !< Reset stagnation counter
                res_norm_prev = huge(1d0)               !< Reset previous residual norm to \f$+\infty\f$
                lambda_LM = 1d-3                        !< Reset LM damping parameter to initial value
                best_res_norm = huge(1d0)               !< Reset best residual norm to \f$+\infty\f$
                skip_eval = .false.                     !< Force re-evaluation after retry
                cycle newton_loop                       !< Restart Newton loop with perturbed initial guess
            end if                                      !< End retry guard
        !> Restore best solution and re-speciate for consistency
            conc_nc=conc_nc_best                        !< Restore best-seen concentrations
            call this%set_conc_prim_species(conc_nc(1:n_p))  !< Push best primary concentrations into chemistry object
            call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p),log_gamma_var_act,conc_nc(n_p+1:n_v))  !< Recompute secondary species at best point for consistency
        !> Accept if best residual is within conditioning noise: \f$\|\mathbf{r}_{\text{best}}\|_\infty \le \sqrt{\varepsilon}\cdot\max(\|\mathbf{u}\|,1)\f$
            if (best_res_norm<=sqrt_eps_d*conc_comp_scale) then  !< Within conditioning noise?
                CV_flag=.true.                          !< Accept as converged
            else                                        !< Best residual exceeds noise threshold
                print *, "Too many Newton iterations in speciation, best residual: ", best_res_norm  !< Report failure
            end if                                      !< End conditioning-noise acceptance
            exit newton_loop                            !< Leave Newton loop (converged or failed)
        end if                                          !< End max-iteration guard
    !> Compute secondary species and residual (skipped if reusing accepted line-search results)
        if (skip_eval) then                            !< Reuse residual from accepted line-search step
            skip_eval = .false.                         !< Consume flag; next iteration will recompute unless step accepted again
        else                                            !< Compute fresh secondary species and residual
            call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p),log_gamma_var_act,conc_nc(n_p+1:n_v))  !< Mass-action law
            call this%compute_res_spec(conc_comp,conc_nc,residual)  !< \f$ \mathbf{r} = \mathbf{U}_1 \mathbf{c}_1 + \mathbf{U}_2\mathbf{c}_{2,v} - \mathbf{u} \f$
        end if                                          !< End skip_eval check
        res_norm=inf_norm_vec_real(residual)           !< \f$ \|\mathbf{r}\|_\infty \f$
    !> Track best solution: keep the iterate with the smallest residual norm
        if (res_norm<best_res_norm) then
            best_res_norm=res_norm
            conc_nc_best=conc_nc
        end if
    !> Check convergence: machine-epsilon safeguard, OR (absolute AND relative) arithmetic, OR (absolute AND relative) logarithmic
    !> Machine-epsilon: \f$ \|\mathbf{r}\|_\infty \le \varepsilon_{\text{mach}} \cdot \max(\|\mathbf{u}\|_\infty, 1) \f$ (round-off floor)
        if ((res_norm<=eps_d*conc_comp_scale) .or. &
            (res_norm<abs_tol .and. &
            res_norm<rel_tol*conc_comp_scale) .or. &
            (res_norm>0d0 .and. &
            (log10(res_norm)<log_abs_tol .and. &
            (log10(res_norm)-log10(conc_comp_scale))< &
            log_rel_tol))) then
            CV_flag=.true.                             !< Mark as converged
            exit newton_loop                            !< Leave Newton loop — converged
        end if                                          !< End convergence check
    !> Stagnation detection: check if \f$ \|\mathbf{r}\| < (1-\text{stag\_rtol})\|\mathbf{r}_\text{prev}\| \f$
        if (res_norm<res_norm_prev*(1d0-stag_rtol)) then
            n_stag=0                                   !< Sufficient improvement: reset stagnation counter
            res_norm_prev=res_norm
        else
            n_stag=n_stag+1                            !< No sufficient improvement: increment stagnation counter
            if (n_stag>=n_stag_max) then               !< Stagnation limit reached?
                if (.not. retry_applied) then           !< Has perturbation retry been tried yet?
                    !> Stagnation detected: try perturbation retry before giving up.
                    !> Geometric mean: \f$ c_{1,j}^{\text{retry}} = \sqrt{c_{1,j}^{\text{ig}} \cdot c_{1,j}^{\text{best}}} \f$
                    retry_applied = .true.              !< Mark retry as used (only one attempt allowed per call)
                    do i = 1, n_p                       !< Loop over each primary species
                        conc_nc(i) = sqrt(max(c1_ig(i), 1d-30) * max(conc_nc_best(i), 1d-30))  !< Geometric mean of initial guess and best iterate
                    end do                              !< End primary species perturbation loop
                    niter = 0                           !< Reset Newton iteration counter for fresh restart
                    n_stag = 0                          !< Reset stagnation counter
                    res_norm_prev = huge(1d0)           !< Reset previous residual norm to \f$+\infty\f$
                    lambda_LM = 1d-3                    !< Reset LM damping parameter to initial value
                    best_res_norm = huge(1d0)           !< Reset best residual norm to \f$+\infty\f$
                    skip_eval = .false.                 !< Force re-evaluation after retry
                    cycle newton_loop                   !< Restart Newton loop with perturbed initial guess
                end if                                  !< End stagnation retry guard
                !> Restore best solution and re-speciate for consistency
                conc_nc=conc_nc_best                    !< Restore best-seen concentrations
                call this%set_conc_prim_species(conc_nc(1:n_p))  !< Push best primary concentrations into chemistry object
                call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p),log_gamma_var_act,conc_nc(n_p+1:n_v))  !< Recompute secondary species at best point for consistency
            !> Accept if best residual is within conditioning noise: \f$\|\mathbf{r}_{\text{best}}\|_\infty \le \sqrt{\varepsilon}\cdot\max(\|\mathbf{u}\|,1)\f$
                if (best_res_norm<=sqrt_eps_d*conc_comp_scale) then  !< Within conditioning noise?
                    CV_flag=.true.                      !< Accept as converged
                else                                    !< Best residual exceeds noise threshold
                    !print *, "Newton speciation stagnated after", niter, "iterations, residual: ", best_res_norm  !< Report failure (disabled)
                end if                                  !< End conditioning-noise acceptance at stagnation
                exit newton_loop                        !< Leave Newton loop (converged or stagnated)
            end if                                      !< End stagnation-limit check
        end if                                          !< End stagnation branch
    !> Compute Jacobian \f$ \partial\mathbf{c}_{2,v}/\partial\mathbf{c}_1 \f$ analytically for ideal solution
        call this%compute_dc2v_dc1_ideal(conc_nc(1:n_p),conc_nc(n_p+1:n_v),dc2v_dc1)
    !> Assemble log-space Jacobian directly: \f$ J_p(:,j) = \ln(10)\,c_{1,j}\,(\mathbf{U}_1 + \mathbf{U}_2\,\partial\mathbf{c}_{2,v}/\partial\mathbf{c}_1)_{:,j} \f$
        mat_lin_syst_p=U1+matmul(U2,dc2v_dc1)         !< J = U1 + U2 * dc2v/dc1 using cached slices
        do i=1,n_p
            log_c1(i)=log10(max(conc_nc(i),1d-30)) !< Clamp to 1d-30 to avoid log(0) and underflow in mass-action
            mat_lin_syst_p(:,i)=ln10*conc_nc(i)*mat_lin_syst_p(:,i) !< Transform column j to log-space in-place
        end do
    !> Save current concentrations for line search recovery
        conc_nc_save=conc_nc
    !> Levenberg-Marquardt loop with backtracking line search: try increasing \f$ \lambda \f$ until a good step is found
        step_accepted=.false.                          !< No step accepted yet
        do i_LM=1,max_LM_tries                         !< LM damping loop: try up to max_LM_tries values of \f$\lambda\f$
        !> Apply LM damping: solve \f$ (\mathbf{J}_p + \lambda\mathbf{I})\,\Delta\mathbf{p} = -\mathbf{r} \f$
            mat_lin_syst_LM=mat_lin_syst_p             !< Copy log-space Jacobian
            do i=1,n_p
                mat_lin_syst_LM(i,i)=mat_lin_syst_LM(i,i)+lambda_LM !< Add \f$ \lambda \f$ to diagonal
            end do
            Delta_p = -residual                        !< Set RHS = -r (dgesv overwrites with solution)
            call dgesv(n_p, 1, mat_lin_syst_LM, n_p, ipiv_dgesv, Delta_p, n_p, info_dgesv)  !< LAPACK LU with partial pivoting
        !> Handle singular system: increase damping and retry
            sing_flag = (info_dgesv /= 0)
            if (sing_flag) then
                lambda_LM=min(lambda_LM*10d0,1d8)     !< Increase \f$ \lambda \f$ by factor 10
                cycle
            end if
        !> Check for NaN or huge values in the Newton step \f$ \Delta\mathbf{p} \f$
            sing_flag=.false.
            do i=1,n_p
                if (Delta_p(i)/=Delta_p(i) .or. abs(Delta_p(i))>1d20) then !< NaN check: x/=x is .true. only for NaN
                    sing_flag=.true.
                    exit
                end if
            end do
            if (sing_flag) then
                lambda_LM=min(lambda_LM*10d0,1d8)     !< Treat as singular: increase damping
                cycle
            end if
        !> Bound log-step: uniform scaling \f$ \Delta\mathbf{p} \leftarrow s\,\Delta\mathbf{p} \f$ so \f$ \|\Delta\mathbf{p}\|_\infty \leq \text{max\_dp} \f$
            scale_factor=1d0
            do i=1,n_p
                if (abs(Delta_p(i))>0d0) then
                    scale_factor=min(scale_factor,max_dp/abs(Delta_p(i)))
                end if
            end do
            if (scale_factor<1d0) Delta_p=Delta_p*scale_factor !< Scale only if step exceeds bound
        !> Backtracking line search in log-space: try \f$ \alpha = 1, 0.5, 0.25, \ldots \f$
            alpha=1d0                                  !< Start with full Newton step
            do i_LS=1,max_LS_tries
            !> Update primary concentrations: \f$ c_{1,j} = 10^{p_j^{\text{old}} + \alpha\,\Delta p_j} \f$
                do i=1,n_p
                    conc_nc(i)=10d0**(log_c1(i)+alpha*Delta_p(i))
                end do
            !> Update internal state and recompute secondary species at trial point
                call this%set_conc_prim_species(conc_nc(1:n_p))
                call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p),log_gamma_var_act,conc_nc(n_p+1:n_v))
            !> Evaluate trial residual at the candidate point
                call this%compute_res_spec(conc_comp,conc_nc,residual_trial)
                res_trial_norm=inf_norm_vec_real(residual_trial)
            !> Accept step if trial residual is smaller (monotone decrease criterion)
                if (res_trial_norm<res_norm) then
                    step_accepted=.true.
                    lambda_LM=max(lambda_LM/10d0,1d-12) !< Reward success: decrease damping toward pure Newton
                    exit
                end if
                alpha=alpha*0.5d0                      !< Halve step length for next backtrack attempt
            end do
            if (step_accepted) exit
        !> Line search failed for this \f$ \lambda \f$: restore saved state and increase damping
            conc_nc=conc_nc_save                       !< Restore pre-LM concentrations
            lambda_LM=min(lambda_LM*10d0,1d8)          !< Increase damping by \f$10\times\f$ for next LM attempt
        end do                                          !< End LM damping loop
    !> Cache accepted step results to skip redundant recomputation in next iteration
        if (step_accepted) then                         !< Step was accepted by LM+line-search
            residual = residual_trial                   !< Reuse accepted trial residual (already at current conc_nc)
            skip_eval = .true.                          !< Signal next iteration to skip residual/secondary computation
        end if                                          !< End accepted-step caching
    !> All LM attempts exhausted: restore saved state and reset \f$\lambda\f$ to initial value
        if (.not.step_accepted) then                    !< All max_LM_tries failed?
            conc_nc=conc_nc_save                        !< Fall back to the state before this Newton step
            lambda_LM=1d-3                              !< Reset damping parameter to default
        end if                                          !< End LM-failure recovery
    end do newton_loop                                  !< End main Newton iteration loop
end subroutine                                          !< End compute_c_nc_from_u_Newton_ideal