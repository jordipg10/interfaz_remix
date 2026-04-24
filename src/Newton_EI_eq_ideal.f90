!> \file Newton_EI_eq_ideal.f90
!> \brief Robust Newton method for equilibrium-only reactive mixing (ideal solution).
!>
!> \details
!> This subroutine solves the nonlinear speciation system arising from reactive mixing
!> when only equilibrium reactions are present (no kinetic reactions). It uses:
!> - **Ideal solution** assumption (unit activity coefficients)
!> - **Log-space Newton steps** to eliminate the \f$1/c_1\f$ singularity in the Jacobian
!> - **Levenberg-Marquardt damping** with backtracking line search for global convergence
!> - **Stagnation detection** and **best-solution tracking** for robustness
!> - **Re-speciation fallback** via `compute_c_nc_from_u_Newton_ideal` when
!>   max iterations or stagnation is reached (tried once before giving up)
!> - **Machine-epsilon safeguard** to accept solutions at double-precision limits
!>
!> **Mathematical Formulation:**
!>
!> The nonlinear residual to drive to zero is:
!> \f[
!>   \mathbf{f} = \mathbf{U}\,\mathbf{c}_v - \hat{\mathbf{u}} = \mathbf{0}
!> \f]
!>
!> where:
!> - \f$ \mathbf{U} \f$ = Component matrix mapping species to components
!>   [\f$n_p \times n_v\f$]
!> - \f$ \mathbf{c}_v \f$ = Variable activity species concentrations [C]
!>   (\f$n_v = n_p + n_{eq}\f$)
!> - \f$ \hat{\mathbf{u}} \f$ = Component concentrations after mixing [C]
!>   (size \f$n_p\f$)
!>
!> The secondary variable activity species \f$\mathbf{c}_{2,v}\f$ are determined from the primary species
!> \f$\mathbf{c}_1\f$ via the mass action law (ideal solution), so the system has
!> \f$n_p\f$ unknowns (the primary species concentrations).
!>
!> **Log-Space Jacobian Transformation:**
!>
!> Defining \f$ p_j = \log_{10}(c_{1,j}) \f$, the Newton system is solved for
!> \f$ \Delta\mathbf{p} \f$ instead of \f$ \Delta\mathbf{c}_1 \f$:
!> \f[
!>   \frac{\partial \mathbf{f}}{\partial p_j}
!>     = \ln(10)\,c_{1,j}\,\frac{\partial \mathbf{f}}{\partial c_{1,j}}
!> \f]
!> This cancels the \f$ 1/c_1 \f$ singularity in \f$ \partial c_2/\partial c_1 \f$,
!> producing a well-conditioned system even for trace-level species.
!>
!> The arithmetic Jacobian is:
!> \f[
!>   \frac{\partial \mathbf{f}}{\partial \mathbf{c}_1}
!>     = \mathbf{U}_1 + \mathbf{U}_2\,\frac{\partial \mathbf{c}_{2,v}}{\partial \mathbf{c}_1}
!> \f]
!> where \f$\mathbf{U}_1\f$ and \f$\mathbf{U}_2\f$ are the sub-blocks of
!> \f$\mathbf{U}\f$ corresponding to primary and secondary species, respectively.
!>
!> **Levenberg-Marquardt Damping:**
!>
!> The damped linear system solved at each Newton iteration is:
!> \f[
!>   \left(\frac{\partial \mathbf{f}}{\partial \mathbf{p}}
!>         + \lambda\,\mathbf{I}\right) \Delta\mathbf{p} = -\mathbf{f}
!> \f]
!> where \f$ \lambda \f$ is adapted: decreased by 10\f$\times\f$ on accepted steps,
!> increased by 10\f$\times\f$ on rejected steps, bounded in
!> \f$[10^{-12},\,10^{8}]\f$.
!>
!> **Algorithm:**
!> 1. **Pre-processing:** Allocate arrays, compute ideal activity coefficients.
!> 2. **Newton iteration loop:**
!>    a. Check iteration count; if exceeded, attempt **re-speciation fallback**
!>       (once) by calling `compute_c_nc_from_u_Newton_ideal` with best-seen
!>       \f$\mathbf{c}_1\f$ as initial guess and \f$\hat{\mathbf{u}}\f$ as target,
!>       then restart the Newton loop. If the fallback has already been tried,
!>       accept best if within conditioning noise or exit with failure.
!>    b. Compute secondary species via mass action law (ideal solution).
!>    c. Compute residual \f$ \mathbf{f} = \mathbf{U}\mathbf{c} - \hat{\mathbf{u}} \f$
!>       and its infinity norm.
!>    d. Track best solution (lowest \f$ \|\mathbf{f}\|_\infty \f$).
!>    e. **Convergence check** (triple criteria — see below).
!>    f. If not converged:
!>       - Compute arithmetic Jacobian
!>         \f$ \partial\mathbf{f}/\partial\mathbf{c}_1 \f$.
!>       - Transform to log-space Jacobian
!>         \f$ \partial\mathbf{f}/\partial\mathbf{p} \f$.
!>       - **LM inner loop** (up to `max_LM_tries = 10`):
!>         * Form damped system and solve with LAPACK `dgesv`.
!>         * Check for singular / NaN solutions.
!>         * Uniformly scale \f$ \Delta\mathbf{p} \f$ so
!>           \f$ \|\Delta\mathbf{p}\|_\infty \le \f$ `max_dp` (preserves direction).
!>         * **Backtracking line search** (up to `max_LS_tries = 10`):
!>           halve step length \f$ \alpha \f$ until
!>           \f$ \|\mathbf{f}^{\mathrm{trial}}\|_\infty <
!>              \|\mathbf{f}\|_\infty \f$.
!>         * Accept step or increase \f$ \lambda \f$ and retry.
!>       - If all LM attempts fail, restore saved state and reset \f$ \lambda \f$.
!>    g. **Stagnation detection:** If residual has not improved by factor
!>       \f$ (1 - \texttt{stag\_rtol}) \f$ for `n_stag_max` consecutive iterations,
!>       attempt **re-speciation fallback** (once, same as step 2a) then restart.
!>       If fallback was already tried, accept best if within conditioning noise or exit.
!> 3. **Post-processing:** Deallocate temporaries.
!>
!> **Convergence Criteria (triple test — any branch suffices):**
!>
!> *Machine-epsilon branch:*
!> \f[
!>   \|\mathbf{f}\|_\infty \le \varepsilon_{\mathrm{mach}}
!>     \cdot \max\!\bigl(\|\hat{\mathbf{u}}\|_\infty,\,1\bigr)
!> \f]
!>
!> *Arithmetic branch:*
!> \f[
!>   \|\mathbf{f}\|_\infty < \varepsilon_{\mathrm{abs}}
!>   \quad\text{AND}\quad
!>   \|\mathbf{f}\|_\infty < \varepsilon_{\mathrm{rel}}
!>     \cdot \max\!\bigl(\|\hat{\mathbf{u}}\|_\infty,\,1\bigr)
!> \f]
!>
!> *Logarithmic branch:*
!> \f[
!>   \log_{10}\|\mathbf{f}\|_\infty < \varepsilon_{\mathrm{abs}}^{\log}
!>   \quad\text{AND}\quad
!>   \log_{10}\|\mathbf{f}\|_\infty
!>     - \log_{10}\!\max\!\bigl(\|\hat{\mathbf{u}}\|_\infty,\,1\bigr)
!>     < \varepsilon_{\mathrm{rel}}^{\log}
!> \f]
!>
!> **Re-speciation Fallback:**
!>
!> When the Newton loop exhausts its iteration budget (`niter_max`) or stagnates
!> (`n_stag_max` consecutive iterations without sufficient improvement), a single
!> re-speciation attempt is made before giving up. This calls
!> `compute_c_nc_from_u_Newton_ideal` with the best-seen primary concentrations
!> \f$\mathbf{c}_1^{\mathrm{best}}\f$ as initial guess and the target components
!> \f$\hat{\mathbf{u}}\f$, effectively restarting Newton from a different basin of
!> attraction. All counters (iteration, stagnation, LM damping, best-tracking) are
!> reset. The fallback is tried at most once per call to avoid infinite cycling.
!>
!> **Conditioning Noise Acceptance:**
!>
!> At the max-iterations and stagnation exits, if the best residual satisfies:
!> \f[
!>   \|\mathbf{f}_{\mathrm{best}}\|_\infty
!>     \le \sqrt{\varepsilon_{\mathrm{mach}}}
!>     \cdot \max\!\bigl(\|\hat{\mathbf{u}}\|_\infty,\,1\bigr)
!> \f]
!> the solution is accepted as converged. This threshold (\f$\approx 1.49\times10^{-8}\f$)
!> accounts for accumulated floating-point noise in the residual evaluation.
!>
!> **Assumptions:**
!> - All primary species are aqueous.
!> - Ideal solution behavior (activity coefficients = 1).
!> - No kinetic reactions — only equilibrium speciation.
!>
!> \param[in,out] this     Aqueous chemistry object (class `aqueous_chemistry_c`)
!> \param[in]     u_hat    Component concentrations after mixing
!>                          \f$\hat{\mathbf{u}}\f$ [C] (size \f$n_p\f$)
!> \param[in,out] conc_nc  Variable activity species concentrations [C]
!>                          (size \f$n_v\f$; on entry = initial guess,
!>                          on exit = converged solution if CV_flag is TRUE)
!> \param[out]    niter    Number of Newton iterations performed [-]
!> \param[out]    CV_flag  Convergence flag: `.FALSE.` = did not converge,
!>                          `.TRUE.` = converged
!>
!> \pre  The chemistry object `this` must have valid `solid_chemistry`,
!>       reactive-zone stoichiometry, and speciation algorithm data.
!> \pre  `conc_nc` must be pre-allocated with size \f$n_v\f$ and contain a
!>       physically meaningful initial guess (all positive entries).
!> \post On convergence, `conc_nc` contains the converged species concentrations.
!>
!> \warning Only valid for ideal solutions (unit activity coefficients).
!> \warning LAPACK `dgesv` is called directly for the LM-damped log-space system.
!>
!> \note  The log-space transformation guarantees positivity of concentrations
!>        by construction: \f$ c_{1,j}^{\mathrm{new}} = 10^{p_j + \alpha\Delta p_j} > 0 \f$.
!>
!> \sa Newton_EI_eq_kin_anal_ideal_opt2, compute_c_nc_from_u_Newton_ideal
!> \sa aqueous_chemistry_c::compute_dc2v_dc1_ideal
!>
!> \author Jordi
!> \date   Unknown
!> \ingroup chemistry
    
subroutine Newton_EI_eq_ideal(this,u_hat,conc_nc,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c  !< Import aqueous chemistry derived type
    implicit none                                       !< Require explicit declaration of all variables
!> Arguments
    class(aqueous_chemistry_c) :: this              !< Aqueous chemistry object (polymorphic self-argument)
    real(kind=8), intent(in) :: u_hat(:)            !< [C] Component concentrations after mixing, size(n_p)
    real(kind=8), intent(inout) :: conc_nc(:)       !< [C] Variable activity species concentrations, size(n_v)
    integer(kind=4), intent(out) :: niter           !< [-] Number of Newton iterations performed
    logical, intent(out) :: CV_flag                 !< Convergence flag
!> Local variables — Jacobians, residual, concentrations
    real(kind=8), allocatable :: log_gamma_var_act(:) !< Log activity coefficients for variable activity species (all zero for ideal solution) (n_v) [-]
    real(kind=8), allocatable :: dc2v_dc1(:,:)      !< Jacobian \f$\partial\mathbf{c}_{2,v}/\partial\mathbf{c}_1\f$ (n_eq × n_p) [-]
    real(kind=8), allocatable :: dfk_dc1(:,:)       !< Jacobian of residual w.r.t. primary species \f$\partial\mathbf{f}/\partial\mathbf{c}_1\f$ (n_p × n_p) [-]
    real(kind=8), allocatable :: fk(:)              !< Newton residual vector \f$\mathbf{f}\f$ (n_p) [C]
    real(kind=8), allocatable :: conc_nc_save(:)    !< Saved full concentrations before LM/line-search trial (n_v) [C]
    real(kind=8), allocatable :: conc_nc_best(:)    !< Best concentrations seen during Newton (lowest residual) (n_v) [C]
!> Local variables — log-space Newton system (hoisted out of Newton loop)
    real(kind=8), allocatable :: log_c1(:)          !< Log10 of primary concentrations \f$p_j = \log_{10}(c_{1,j})\f$ (n_p) [-]
    real(kind=8), allocatable :: Delta_p(:)          !< Newton step in log-space \f$\Delta\mathbf{p}\f$ (n_p) [-]
    real(kind=8), allocatable :: dfk_dp(:,:)         !< Log-space Jacobian \f$\partial\mathbf{f}/\partial\mathbf{p}\f$ where \f$(\cdot)_{:,j} = \ln(10)\,c_{1,j}\,(\partial\mathbf{f}/\partial c_{1,j})\f$ (n_p × n_p) [-]
    real(kind=8), allocatable :: dfk_dp_LM(:,:)      !< LM-damped log-space Jacobian \f$\partial\mathbf{f}/\partial\mathbf{p} + \lambda\mathbf{I}\f$ (n_p × n_p) [-]
    real(kind=8), parameter :: ln10 = log(10d0)      !< Natural logarithm of 10, \f$\ln(10) \approx 2.302585\f$ [-]
    real(kind=8), parameter :: max_dp = 10d0         !< Maximum allowed \f$|\Delta p_j|\f$ per iteration (10 orders of magnitude) [-]
    real(kind=8) :: scale_factor                      !< Uniform scaling factor for \f$\Delta\mathbf{p}\f$ to enforce max_dp while preserving direction [-]
!> Local variables — Levenberg-Marquardt and line search
    real(kind=8) :: lambda_LM                       !< Levenberg-Marquardt damping parameter \f$\lambda\f$, adapted in \f$[10^{-12},\,10^{8}]\f$ [-]
    integer(kind=4) :: i_LM                          !< LM inner loop counter [-]
    integer(kind=4), parameter :: max_LM_tries = 10  !< Maximum LM damping attempts per Newton iteration [-]
    real(kind=8) :: alpha                             !< Backtracking line search step length \f$\alpha \in (0,1]\f$ [-]
    integer(kind=4) :: i_LS                           !< Line search iteration counter [-]
    integer(kind=4), parameter :: max_LS_tries = 10   !< Maximum backtracking line search steps per LM trial [-]
    logical :: step_accepted                          !< Whether a step was accepted in the LM + line search procedure
    real(kind=8), allocatable :: fk_trial(:)          !< Trial residual vector for line search evaluation (n_p) [C]
    real(kind=8) :: fk_norm                           !< Current residual infinity norm \f$\|\mathbf{f}\|_\infty\f$ [C]
    real(kind=8) :: fk_trial_norm                     !< Trial residual infinity norm after line search step [C]
!> Local variables — stagnation detection and re-speciation fallback
    integer(kind=4) :: n_stag                       !< Consecutive iterations without sufficient residual improvement [-]
    integer(kind=4), parameter :: n_stag_max = 20    !< Maximum stagnation iterations before early exit or re-speciation fallback [-]
    real(kind=8) :: fk_norm_prev                      !< Previous best residual norm for stagnation comparison [C]
    real(kind=8), parameter :: stag_rtol = 1d-3      !< Relative improvement threshold; stagnation counter resets if \f$\|\mathbf{f}\| < (1-\texttt{stag\_rtol})\,\|\mathbf{f}\|_{\mathrm{prev}}\f$ [-]
    real(kind=8) :: best_fk_norm                      !< Lowest \f$\|\mathbf{f}\|_\infty\f$ seen during the entire Newton loop [C]
    logical :: predictor_applied                      !< Whether the re-speciation fallback has already been tried
    integer(kind=4) :: niter_spec                     !< Number of inner speciation Newton iterations (re-speciation fallback) [-]
    logical :: CV_flag_spec                           !< Convergence flag from inner speciation Newton (re-speciation fallback)
!> Local variables — dimensions and misc
    integer(kind=4) :: n_p                           !< Number of primary species [-]
    integer(kind=4) :: n_v                           !< Number of variable activity species (= n_p + n_eq) [-]
    integer(kind=4) :: n_eq                          !< Number of equilibrium reactions (= n_v - n_p) [-]
    integer(kind=4) :: ii                            !< General-purpose loop index [-]
    logical :: sing_flag                              !< Singularity flag: .TRUE. if dgesv reports failure or NaN/huge values in solution
    integer(kind=4), allocatable :: ipiv_dgesv(:)    !< Pivot index array for LAPACK dgesv (n_p)
    integer(kind=4) :: info_dgesv                    !< Return code from LAPACK dgesv (0 = success)
    logical :: skip_eval                              !< Skip residual/speciation recomputation when accepted step can be reused
    real(kind=8) :: u_norm                            !< Infinity norm of current component concentrations \f$\|\mathbf{U}\mathbf{c}_{v}\|_\infty\f$ [C]
    real(kind=8), allocatable :: Uc(:)                !< Cached product \f$\mathbf{U}\mathbf{c}_{v}\f$ to avoid redundant matmul (n_p) [C]
    real(kind=8), allocatable :: col_scale(:)         !< Column scaling factors for Jacobian equilibration (n_p) [-]
!> Cached loop-invariant parameters (avoid repeated deep struct dereferences per iteration)
    integer(kind=4) :: niter_max                     !< Maximum Newton iterations (cached from CV_params)
    real(kind=8) :: abs_tol                          !< Absolute convergence tolerance (cached)
    real(kind=8) :: rel_tol                          !< Relative convergence tolerance (cached)
    real(kind=8) :: log_abs_tol                      !< Log absolute tolerance (cached)
    real(kind=8) :: log_rel_tol                      !< Log relative tolerance (cached)
    real(kind=8) :: eps_d                            !< Machine epsilon \f$\varepsilon_{\mathrm{mach}}\f$ (cached) [-]
    real(kind=8) :: sqrt_eps_d                       !< \f$\sqrt{\varepsilon_{\mathrm{mach}}}\f$ — practical conditioning noise floor for double precision [-]
    real(kind=8), allocatable :: comp_mat(:,:)       !< Component matrix \f$\mathbf{U}\f$ (cached, n_p × n_v) [-]
    real(kind=8), allocatable :: U1(:,:)             !< Cached primary species block \f$\mathbf{U}_1 = \mathbf{U}(:,1:n_p)\f$ (n_p × n_p) [-]
    real(kind=8), allocatable :: U2(:,:)             !< Cached secondary species block \f$\mathbf{U}_2 = \mathbf{U}(:,n_p+1:n_v)\f$ (n_p × n_eq) [-]
!> -----------------------------------------------------------------------
!> \section newton_eq_preproc Pre-processing: initialisation and allocation
!> -----------------------------------------------------------------------
    niter = 0                       !< Reset iteration counter to zero
    CV_flag = .false.               !< Assume non-convergence until proven otherwise
    n_p = this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species       !< Retrieve number of primary (independent) species
    n_v = this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species    !< Retrieve total number of variable-activity species
    !> Cache loop-invariant parameters from deep struct dereferences
    niter_max = this%solid_chemistry%reactive_zone%CV_params%niter_max   !< Cache max Newton iterations
    abs_tol = this%solid_chemistry%reactive_zone%CV_params%abs_tol       !< Cache absolute tolerance
    rel_tol = this%solid_chemistry%reactive_zone%CV_params%rel_tol       !< Cache relative tolerance
    log_abs_tol = this%solid_chemistry%reactive_zone%CV_params%log_abs_tol !< Cache log absolute tolerance
    log_rel_tol = this%solid_chemistry%reactive_zone%CV_params%log_rel_tol !< Cache log relative tolerance
    eps_d = epsilon(1d0)                                                 !< Cache machine epsilon
    sqrt_eps_d = sqrt(eps_d)                                             !< Cache sqrt(eps) ≈ 1.49e-8 (practical conditioning noise floor)
    comp_mat = this%solid_chemistry%reactive_zone%speciation_alg%comp_mat !< Cache component matrix U
    n_eq = n_v - n_p                 !< Number of equilibrium reactions
    !> Cache component-matrix sub-blocks (avoids creating temporary array slices each iteration)
    U1 = comp_mat(:,1:n_p)          !< Primary species block U1 (n_p × n_p)
    U2 = comp_mat(:,n_p+1:n_v)      !< Secondary species block U2 (n_p × n_eq)
    !> Allocate persistent work arrays used throughout the Newton loop
    allocate(dfk_dc1(n_p,n_p), dc2v_dc1(n_eq, n_p), fk(n_p), conc_nc_best(n_v), conc_nc_save(n_v))  !< Jacobian, residual, best/save concentration arrays
    !> Allocate log-space work arrays once (hoisted out of the Newton loop)
    allocate(dfk_dp(n_p,n_p), dfk_dp_LM(n_p,n_p), Delta_p(n_p), log_c1(n_p), Uc(n_p))  !< Log-Jacobian, LM copy, step, log-conc, cached U*c
    allocate(fk_trial(n_p), ipiv_dgesv(n_p))  !< Trial residual vector and LAPACK pivot array
    allocate(col_scale(n_p))  !< Column scaling factors for Jacobian equilibration
    lambda_LM = 1d-3                !< Set initial LM damping parameter (moderate damping)
    n_stag = 0                      !< Reset stagnation counter
    fk_norm_prev = huge(1d0)        !< Initialise previous residual norm to +∞ (any first residual will be an improvement)
    best_fk_norm = huge(1d0)        !< Initialise best residual norm to +∞
    predictor_applied = .false.     !< Re-speciation fallback has not been tried yet
    skip_eval = .false.             !< No accepted step yet; must compute from scratch
    log_gamma_var_act = this%get_log_gamma_var_act()  !< Retrieve log10 activity coefficients (ideal ⇒ all zeros)
!> -----------------------------------------------------------------------
!> \section newton_eq_process Process: activity coefficients and Newton loop
!> -----------------------------------------------------------------------
    !> Compute ideal activity coefficients (all zeros for ideal solution) — already retrieved above.

    !> -----------------------------------------------------------------------
    !> \subsection newton_eq_loop Newton iteration loop
    !> -----------------------------------------------------------------------
    newton_loop: do                    !< Begin main Newton iteration loop
        niter = niter + 1              !< Increment Newton iteration counter
    !> \subsubsection newton_eq_maxiter Maximum iteration check and re-speciation fallback
    !> If niter exceeds niter_max, the loop has exhausted its budget.
    !> Before giving up, try re-speciation fallback (once), then check conditioning noise.
        if (niter > niter_max) then     !< Check whether iteration budget is exhausted
            if (.not. predictor_applied) then  !< First exhaustion: try re-speciation fallback
                !> Re-speciation fallback: call the inner speciation Newton solver
                !> with the best-seen c1 as initial guess and u_hat as component concentrations.
                !> This effectively restarts from a cleaner speciation state.
                predictor_applied = .true.                  !< Mark fallback as used (only one attempt allowed)
                call this%compute_c_nc_from_u_Newton_ideal(conc_nc_best(1:n_p), u_hat, conc_nc, niter_spec, CV_flag_spec)  !< Re-speciate from best c1 to obtain new initial guess
                niter = 0                                   !< Reset outer Newton counter for fresh restart
                n_stag = 0                                  !< Reset stagnation counter
                fk_norm_prev = huge(1d0)                    !< Reset previous residual norm to +∞
                lambda_LM = 1d-3                            !< Reset LM damping to initial value
                best_fk_norm = huge(1d0)                    !< Reset best residual norm to +∞
                skip_eval = .false.                         !< Force re-evaluation after fallback restart
                cycle newton_loop                           !< Restart Newton loop with new initial guess
            end if                                          !< End re-speciation fallback block
            !> Test conditioning-noise acceptance: \f$\|\mathbf{f}_{\mathrm{best}}\|_\infty \le \sqrt{\varepsilon}\cdot\max(\|\mathbf{U}\mathbf{c}\|,1)\f$
            if (best_fk_norm <= sqrt_eps_d * max(u_norm, 1d0)) then  !< Is best residual within conditioning noise?
                CV_flag = .true.        !< Accept: residual is at machine-precision noise level
                conc_nc = conc_nc_best  !< Restore the best concentrations found during the iteration history
            else                        !< Best residual exceeds conditioning noise — true non-convergence
                !> Report failure: solver exhausted iteration budget without reaching tolerance
                print *, 'Newton_EI_eq_ideal did not converge after', &  !< Diagnostic message to stdout
                    niter_max, 'iterations. Best ||fk|| =', best_fk_norm !< Print iteration count and best residual
            end if                      !< End conditioning-noise acceptance test
            exit newton_loop            !< Leave the Newton loop regardless of acceptance
        end if                          !< End maximum-iteration guard
    !> \subsubsection newton_eq_speciation Speciation: primary → secondary → residual
    !> When the previous line search accepts a step, it already speciates and computes
    !> fk_trial.  Reuse that directly to avoid one O(n_p × n_v) matmul.
        if (skip_eval) then             !< Previous LS accepted a step — reuse its speciation and residual
            skip_eval = .false.         !< Clear the skip flag for subsequent iterations
            fk = fk_trial               !< Reuse residual from accepted line-search step
            Uc = fk + u_hat             !< Recover U*c: O(n_p) add vs O(n_p × n_v) matmul
        else                            !< No cached result — compute speciation and residual from scratch
            call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p), log_gamma_var_act, conc_nc(n_p+1:))  !< Compute \f$\mathbf{c}_{2,v}\f$ from \f$\mathbf{c}_1\f$ via mass-action law (ideal)
            call dgemv('N', n_p, n_v, 1d0, comp_mat, n_p, conc_nc, 1, 0d0, Uc, 1)  !< \f$\mathbf{U}\mathbf{c}_v\f$ via BLAS (no temporary)
            fk = Uc - u_hat             !< \f$\mathbf{f} = \mathbf{U}\mathbf{c}_v - \hat{\mathbf{u}}\f$ — Newton residual
        end if                          !< End speciation/residual evaluation branch
        fk_norm = maxval(abs(fk))       !< \f$\|\mathbf{f}\|_\infty\f$ — current residual infinity norm (inline, no call overhead)
        u_norm = maxval(abs(Uc))        !< \f$\|\mathbf{U}\mathbf{c}\|_\infty\f$ — for convergence normalisation
    !> \subsubsection newton_eq_best Best-solution tracking
    !> Keep a running record of the lowest residual norm and its associated concentrations.
    !> This enables recovery at max-iteration and stagnation exits.
        if (fk_norm < best_fk_norm) then  !< New lowest residual found — update best tracking
            best_fk_norm = fk_norm      !< Update best residual norm
            conc_nc_best = conc_nc      !< Snapshot concentrations at the new best point
        end if                          !< End best-solution tracking update
    !> \subsubsection newton_eq_conv Convergence check (triple arithmetic OR logarithmic criteria)
    !> The Newton iterate is accepted as converged if ANY of the following three
    !> independent criteria is satisfied:
    !>
    !> **Branch 1 — Machine-epsilon safeguard:**
    !> \f$\|\mathbf{f}\|_\infty \le \varepsilon_{\mathrm{mach}}\cdot\max(\|\mathbf{U}\mathbf{c}\|,1)\f$
    !>
    !> **Branch 2 — Arithmetic tolerances (abs AND rel):**
    !> \f$\|\mathbf{f}\|_\infty < \varepsilon_{\mathrm{abs}}\f$ AND
    !> \f$\|\mathbf{f}\|_\infty < \varepsilon_{\mathrm{rel}}\cdot\max(\|\mathbf{U}\mathbf{c}\|,1)\f$
    !>
    !> **Branch 3 — Logarithmic tolerances (log-abs AND log-rel):**
    !> \f$\log_{10}\|\mathbf{f}\|_\infty < \varepsilon_{\mathrm{abs}}^{\log}\f$ AND
    !> \f$\log_{10}\|\mathbf{f}\|_\infty - \log_{10}\max(\|\mathbf{U}\mathbf{c}\|,1) < \varepsilon_{\mathrm{rel}}^{\log}\f$
        if ((fk_norm <= eps_d * max(u_norm, 1d0)) .or. &                               !< Branch 1: machine-epsilon
            (fk_norm < abs_tol .and. &                                                  !< Branch 2a: absolute tolerance
             fk_norm < rel_tol * max(u_norm, 1d0)) .or. &                              !< Branch 2b: relative tolerance
            (fk_norm > 0d0 .and. &                                                      !< Branch 3: guard against log10(0)
             (log10(fk_norm) < log_abs_tol .and. &                                      !< Branch 3a: log-absolute tolerance
              (log10(fk_norm) - log10(max(u_norm, 1d0))) < &                            !< Branch 3b: log-relative tolerance
              log_rel_tol))) then       !< Close Branch 3 and entire convergence test
            CV_flag = .true.            !< Mark as converged
            exit newton_loop            !< Leave the Newton loop — solution accepted
        else                            !< Not converged — proceed with Jacobian and Newton step
        !> \subsubsection newton_eq_jacobian Jacobian assembly and log-space transformation
        !> Compute \f$\partial\mathbf{c}_2/\partial\mathbf{c}_1\f$ via the speciation derivatives,
        !> then form the full arithmetic Jacobian:
        !> \f$\frac{\partial\mathbf{f}}{\partial\mathbf{c}_1} = \mathbf{U}_1 + \mathbf{U}_2\,\frac{\partial\mathbf{c}_2}{\partial\mathbf{c}_1}\f$
            call this%compute_dc2v_dc1_ideal(conc_nc(1:n_p), conc_nc(n_p+1:), dc2v_dc1)  !< Compute \f$\partial\mathbf{c}_{2,v}/\partial\mathbf{c}_1\f$ (ideal speciation derivatives)
            dfk_dc1 = U1              !< Start with primary block \f$\mathbf{U}_1\f$
            call dgemm('N', 'N', n_p, n_p, n_eq, 1d0, U2, n_p, dc2v_dc1, n_eq, 1d0, dfk_dc1, n_p)  !< \f$\mathbf{U}_1 + \mathbf{U}_2\,\partial\mathbf{c}_2/\partial\mathbf{c}_1\f$ in-place (no temporary)
        !> \subsubsection newton_eq_logspace Log-space Jacobian transformation
        !> For each primary species \f$j\f$, compute:
        !>  - \f$p_j = \log_{10}(c_{1,j})\f$ (clamped above tiny to avoid \f$-\infty\f$)
        !>  - \f$\partial f_i/\partial p_j = \ln(10)\,c_{1,j}\,\partial f_i/\partial c_{1,j}\f$
            do ii = 1, n_p              !< Loop over each primary species
                log_c1(ii) = log10(max(conc_nc(ii), tiny(1d0)))   !< \f$p_j = \log_{10}(\max(c_{1,j},\, \texttt{tiny}))\f$
                dfk_dp(:,ii) = ln10 * conc_nc(ii) * dfk_dc1(:,ii) !< Column \f$j\f$ of log-Jacobian: \f$\ln(10)\,c_{1,j}\,(\partial\mathbf{f}/\partial c_{1,j})\f$
            end do                      !< End log-space Jacobian transformation loop
        !> \subsubsection newton_eq_lm Levenberg-Marquardt with backtracking line search in log-space
        !> Save the current state so it can be restored if all LM attempts fail.
        !> Try up to max_LM_tries values of \f$\lambda\f$; for each, solve the damped
        !> log-space system and attempt a backtracking line search.
            conc_nc_save = conc_nc      !< Snapshot concentrations before any trial steps
            step_accepted = .false.     !< No step accepted yet
            do i_LM = 1, max_LM_tries  !< LM damping loop: try different \f$\lambda\f$ values
                !> Column equilibration: normalise each column of the log-Jacobian before
                !> adding LM damping.  This dramatically improves conditioning when species
                !> concentrations span many orders of magnitude.
                dfk_dp_LM = dfk_dp      !< Fresh copy of log-Jacobian (dgesv overwrites it)
                do ii = 1, n_p          !< Compute column scaling factors
                    col_scale(ii) = maxval(abs(dfk_dp_LM(:,ii)))  !< Infinity norm of column
                    if (col_scale(ii) < tiny(1d0)) col_scale(ii) = 1d0  !< Guard against zero columns
                    dfk_dp_LM(:,ii) = dfk_dp_LM(:,ii) / col_scale(ii)  !< Normalise column
                end do                  !< End column equilibration loop
                !> Apply Levenberg damping in the equilibrated space
                do ii = 1, n_p          !< Add damping to each diagonal element
                    dfk_dp_LM(ii, ii) = dfk_dp_LM(ii, ii) + lambda_LM  !< \f$J_{ii} \leftarrow J_{ii} + \lambda\f$
                end do                  !< End diagonal damping loop
                !> Solve the damped linear system with LAPACK dgesv (LU factorisation + back-substitution).
                !> On entry, Delta_p = -fk (RHS); on exit, Delta_p = \f$\Delta\mathbf{p}\f$ (solution).
                Delta_p = -fk           !< Set RHS = \f$-\mathbf{f}\f$
                call dgesv(n_p, 1, dfk_dp_LM, n_p, ipiv_dgesv, Delta_p, n_p, info_dgesv)  !< LAPACK direct solve (ipiv reused each attempt)
                !> Unscale the solution: recover original Delta_p from column-equilibrated system
                if (info_dgesv == 0) then
                    do ii = 1, n_p
                        Delta_p(ii) = Delta_p(ii) / col_scale(ii)
                    end do
                end if
                !> Validate the solution: check for dgesv failure (info /= 0) and NaN / overflow.
                !> If invalid, increase \f$\lambda\f$ and retry with stronger damping.
                sing_flag = (info_dgesv /= 0)  !< dgesv returned an error code
                if (.not. sing_flag) then       !< Only scan if dgesv succeeded
                    sing_flag = any(Delta_p /= Delta_p) .or. maxval(abs(Delta_p)) > 1d20  !< Vectorised NaN (x/=x) + overflow check
                end if                          !< End dgesv-success guard
                if (sing_flag) then     !< Solution is unusable
                    lambda_LM = min(lambda_LM * 10d0, 1d8)  !< Increase damping by 10× (capped at 1e8)
                    cycle               !< Retry with stronger damping
                end if                  !< End singular-solution handler
                !> Uniform step bounding: scale \f$\Delta\mathbf{p}\f$ so that
                !> \f$\|\Delta\mathbf{p}\|_\infty \le \texttt{max\_dp}\f$.
                !> This preserves the Newton direction while preventing excessively large steps.
                scale_factor = maxval(abs(Delta_p))  !< Find largest absolute step component
                if (scale_factor > max_dp) Delta_p = Delta_p * (max_dp / scale_factor)  !< Uniform scaling to enforce max_dp (preserves direction)
                !> \paragraph newton_eq_ls Backtracking line search in log-space
                !> Starting from full step (\f$\alpha=1\f$), halve \f$\alpha\f$ until the trial
                !> residual is strictly smaller than the current residual (sufficient decrease).
                alpha = 1d0             !< Initial step length = full Newton step
                do i_LS = 1, max_LS_tries  !< Backtracking loop
                    !> Update primary species in log-space: \f$c_{1,j}^{\mathrm{trial}} = 10^{p_j + \alpha\,\Delta p_j}\f$
                    !> This guarantees positivity by construction.
                    do ii = 1, n_p      !< Apply log-space update to each primary species
                        conc_nc(ii) = 10d0**(log_c1(ii) + alpha * Delta_p(ii))  !< \f$c_{1,j} = 10^{p_j + \alpha\Delta p_j}\f$
                    end do              !< End log-space primary species update
                    !> Recompute secondary species and trial residual at the candidate point
                    call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p), log_gamma_var_act, conc_nc(n_p+1:))  !< Speciation at trial point
                    fk_trial = -u_hat   !< Initialise with \f$-\hat{\mathbf{u}}\f$
                    call dgemv('N', n_p, n_v, 1d0, comp_mat, n_p, conc_nc, 1, 1d0, fk_trial, 1)  !< Trial residual \f$\mathbf{U}\mathbf{c} - \hat{\mathbf{u}}\f$ via BLAS (no temporary)
                    fk_trial_norm = maxval(abs(fk_trial))  !< \f$\|\mathbf{f}^{\mathrm{trial}}\|_\infty\f$ (inline, no call overhead)
                    !> Accept step if trial residual is strictly smaller than current residual
                    if (fk_trial_norm < fk_norm) then  !< Sufficient decrease achieved — accept step
                        step_accepted = .true.  !< Step is accepted
                        skip_eval = .true.       !< Signal next iteration to skip speciation recomputation
                        lambda_LM = max(lambda_LM / 10d0, 1d-12)  !< Decrease damping by 10× (floored at 1e-12) for next iteration
                        exit            !< Exit line-search loop — proceed to next Newton iteration
                    end if              !< End sufficient-decrease acceptance
                    alpha = alpha * 0.5d0  !< Halve step length and retry
                end do                  !< End backtracking line-search loop
                if (step_accepted) exit !< Exit LM loop if line search found an acceptable step
                !> Line search failed for this \f$\lambda\f$: restore concentrations and increase damping
                conc_nc = conc_nc_save  !< Restore pre-LM concentrations
                lambda_LM = min(lambda_LM * 10d0, 1d8)  !< Increase damping by 10× for next LM attempt
            end do                      !< End Levenberg-Marquardt damping loop
            !> All LM attempts exhausted: restore saved state and reset \f$\lambda\f$ to initial value
            if (.not. step_accepted) then  !< No LM attempt yielded an acceptable step
                conc_nc = conc_nc_save  !< Fall back to the state before this Newton step
                lambda_LM = 1d-3        !< Reset damping parameter to default
            end if                      !< End LM-exhaustion fallback
        !> \subsubsection newton_eq_stag Stagnation detection and re-speciation fallback
        !> If the residual norm has decreased by at least the factor \f$(1 - \texttt{stag\_rtol})\f$
        !> compared to the previous iteration, reset the stagnation counter and update the
        !> reference norm. Otherwise, increment the counter. When it reaches n_stag_max,
        !> the solver concludes the residual has plateaued.
            if (fk_norm < fk_norm_prev * (1d0 - stag_rtol)) then  !< Sufficient improvement detected
                n_stag = 0              !< Reset stagnation counter
                fk_norm_prev = fk_norm  !< Update reference norm for next comparison
            else                        !< Residual did not improve sufficiently
                n_stag = n_stag + 1     !< Increment stagnation counter
                if (n_stag >= n_stag_max) then  !< Stagnation limit reached
                    if (.not. predictor_applied) then  !< Has the re-speciation fallback been tried yet?
                        !> Stagnation detected: try re-speciation from best-seen c1 before giving up.
                        predictor_applied = .true.           !< Mark fallback as used (only one attempt allowed)
                        call this%compute_c_nc_from_u_Newton_ideal(conc_nc_best(1:n_p), u_hat, conc_nc, niter_spec, CV_flag_spec)  !< Re-speciate from best c1 (stagnation restart)
                        niter = 0                            !< Reset outer Newton counter for fresh restart
                        n_stag = 0                           !< Reset stagnation counter
                        fk_norm_prev = huge(1d0)             !< Reset previous residual norm to +∞
                        lambda_LM = 1d-3                     !< Reset LM damping to initial value
                        best_fk_norm = huge(1d0)             !< Reset best residual norm to +∞
                        skip_eval = .false.                  !< Force re-evaluation after fallback restart
                        cycle newton_loop                    !< Restart Newton loop with new initial guess
                    end if              !< End stagnation re-speciation fallback
                    !> Test conditioning-noise acceptance at stagnation exit:
                    !> \f$\|\mathbf{f}_{\mathrm{best}}\|_\infty \le \sqrt{\varepsilon}\cdot\max(\|\mathbf{U}\mathbf{c}\|,1)\f$
                    if (best_fk_norm <= sqrt_eps_d * max(u_norm, 1d0)) then  !< Is best residual within conditioning noise?
                        CV_flag = .true.        !< Accept: residual is within conditioning noise
                        conc_nc = conc_nc_best  !< Restore best-seen concentrations
                    end if              !< End conditioning-noise acceptance at stagnation
                    exit newton_loop    !< Leave Newton loop (converged or stagnated)
                end if                  !< End stagnation-limit handler (n_stag >= n_stag_max)
            end if                      !< End stagnation else-branch
        end if                          !< End convergence else-branch (Jacobian + LM + stagnation)
    end do newton_loop                  !< End of main Newton iteration loop
!> -----------------------------------------------------------------------
!> \section newton_eq_postproc Post-processing: cleanup and deallocation
!> -----------------------------------------------------------------------
    deallocate(fk, dfk_dc1, dc2v_dc1, conc_nc_best, conc_nc_save)  !< Free Jacobian, residual, best/save arrays
    deallocate(dfk_dp, dfk_dp_LM, Delta_p, log_c1, fk_trial, ipiv_dgesv, Uc, col_scale)  !< Free log-space work arrays
    deallocate(comp_mat, U1, U2)        !< Free cached component-matrix blocks
end subroutine Newton_EI_eq_ideal       !< End of Newton_EI_eq_ideal subroutine