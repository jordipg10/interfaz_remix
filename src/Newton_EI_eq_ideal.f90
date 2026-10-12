!> \file Newton_EI_eq_ideal.f90
!> \brief Globalised Newton method (Armijo line search) for equilibrium-only reactive mixing (ideal solution).
!>
!> \details
!> This subroutine solves the nonlinear speciation system arising from reactive mixing
!> when only equilibrium reactions are present (no kinetic reactions). It uses:
!> - **Ideal solution** assumption: the \f$\log_{10}\f$ activity coefficients are read once from the
!>   chemistry object and held fixed (all zero in the ideal case)
!> - **Log-space Newton steps** to eliminate the \f$1/c_1\f$ singularity in the Jacobian
!> - **Armijo backtracking line search** on \f$\phi=\tfrac12\|\mathbf{f}\|_2^2\f$ along the Newton
!>   direction, with steepest descent \f$-\mathbf{J}^T\mathbf{f}\f$ as fallback, for global convergence
!> - **Stagnation detection** and **best-solution tracking** for robustness
!> - **Re-speciation fallback** via `compute_c_nc_from_u_Newton_ideal` when
!>   max iterations or stagnation is reached (tried once before giving up)
!> - **Machine-epsilon safeguard** so that the loop terminates when the residual reaches
!>   double-precision round-off
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
!> \f$\mathbf{c}_1\f$ via the mass action law (ideal solution),
!> \f[
!>   \log_{10}\mathbf{c}_{2,v} = \mathbf{S}^*_{e,v,1}\,\log_{10}\mathbf{c}_1 + \log_{10}\mathbf{K}^*,
!> \f]
!> so the system has \f$n_p\f$ unknowns (the primary species concentrations) and
!> \f$n_{eq} = n_v - n_p\f$ equations of state.
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
!>     = \mathbf{U}_1 + \mathbf{U}_2\,\frac{\partial \mathbf{c}_{2,v}}{\partial \mathbf{c}_1},
!>   \qquad
!>   \frac{\partial c_{2,v,i}}{\partial c_{1,j}} = \frac{c_{2,v,i}\,S^*_{e,v,1}(i,j)}{c_{1,j}}
!> \f]
!> where \f$\mathbf{U}_1\f$ and \f$\mathbf{U}_2\f$ are the sub-blocks of
!> \f$\mathbf{U}\f$ corresponding to primary and secondary species, respectively. The
!> \f$1/c_{1,j}\f$ factor is cancelled by the \f$c_{1,j}\f$ of the log-space transformation.
!>
!> **Globalisation — Armijo backtracking line search:**
!>
!> The (column-equilibrated) log-space Newton system
!> \f[
!>   \mathbf{J}\,\Delta\mathbf{p} = -\mathbf{f}, \qquad
!>   \mathbf{J} = \frac{\partial \mathbf{f}}{\partial \mathbf{p}}
!> \f]
!> is solved with LAPACK `dgesv`. The step length \f$\alpha\f$ along \f$\Delta\mathbf{p}\f$ is chosen by
!> backtracking (\f$\alpha = 1, 1/2, 1/4, \ldots\f$) on the merit function
!> \f$\phi(\alpha)=\tfrac12\|\mathbf{f}(\mathbf{p}+\alpha\Delta\mathbf{p})\|_2^2\f$ until the Armijo condition
!> \f$\phi(\alpha)\le\phi(0)+c\,\alpha\,\phi'(0)\f$, \f$c=10^{-4}\f$, holds. Here
!> \f$\phi'(0)=\left.d\phi/d\alpha\right|_{\alpha=0}=\mathbf{f}^T\mathbf{J}\,\Delta\mathbf{p}\f$ is the
!> directional derivative of \f$\tfrac12\|\mathbf{f}\|_2^2\f$ along \f$\Delta\mathbf{p}\f$ in
!> \f$\mathbf{p}=\log_{10}\mathbf{c}_1\f$ space (before the `max_dp` rescaling below: \f$-\|\mathbf{f}\|_2^2\f$ for the
!> exact Newton direction, \f$-\|\mathbf{J}^T\mathbf{f}\|_2^2\f$ for steepest descent). Before the line search
!> \f$\Delta\mathbf{p}\f$ is scaled uniformly so that \f$\|\Delta\mathbf{p}\|_\infty \le\f$ `max_dp` = 10
!> (ten orders of magnitude) and \f$\phi'(0)\f$ is re-evaluated for the scaled direction; \f$\alpha\f$ is
!> halved at most `max_LS_tries` \f$-1 = 9\f$ times (\f$\alpha_{\min}=2^{-9}\f$).
!> If the solve fails, \f$\Delta\mathbf{p}\f$ is not a descent direction
!> (\f$\phi'(0)\ge 0\f$) or its line search fails, the steepest-descent direction
!> \f$-\mathbf{J}^T\mathbf{f}\f$ is tried once instead. If no step is accepted the iterate is kept and the
!> iteration is treated as stagnated (see step 2g), since repeating it would give the same result.
!>
!> **Algorithm:**
!> 1. **Pre-processing:** Cache dimensions, tolerances and the component-matrix blocks
!>    \f$\mathbf{U}_1,\mathbf{U}_2\f$; allocate work arrays; read the \f$\log_{10}\f$ activity
!>    coefficients of the variable activity species from the object (zero for the ideal solution).
!> 2. **Newton iteration loop:**
!>    a. Check iteration count; if exceeded, attempt **re-speciation fallback**
!>       (once) by calling `compute_c_nc_from_u_Newton_ideal` with best-seen
!>       \f$\mathbf{c}_1\f$ as initial guess and \f$\hat{\mathbf{u}}\f$ as target,
!>       then restart the Newton loop if it converged. Otherwise (or if the fallback
!>       has already been tried) return the best iterate with `CV_flag = .false.`.
!>    b. Compute secondary species via mass action law (ideal solution), unless the
!>       previous line search already did so at the accepted point.
!>    c. Compute residual \f$ \mathbf{f} = \mathbf{U}\mathbf{c} - \hat{\mathbf{u}} \f$
!>       and its infinity norm.
!>    d. Track best solution (lowest \f$ \|\mathbf{f}\|_\infty \f$).
!>    e. **Convergence check** (see below).
!>    f. If not converged:
!>       - Compute arithmetic Jacobian
!>         \f$ \partial\mathbf{f}/\partial\mathbf{c}_1 \f$.
!>       - Transform to log-space Jacobian
!>         \f$ \partial\mathbf{f}/\partial\mathbf{p} \f$.
!>       - Solve the column-equilibrated Newton system with LAPACK `dgesv`; on a
!>         singular / NaN solution switch to steepest descent \f$-\mathbf{J}^T\mathbf{f}\f$.
!>       - Descent check: if \f$\phi'(0)=\mathbf{f}^T\mathbf{J}\Delta\mathbf{p}\ge 0\f$,
!>         switch to steepest descent.
!>       - Uniformly scale \f$ \Delta\mathbf{p} \f$ so
!>         \f$ \|\Delta\mathbf{p}\|_\infty \le \f$ `max_dp` (preserves direction) and
!>         evaluate \f$\phi(0)\f$, \f$\phi'(0)\f$ for the scaled direction.
!>       - **Armijo backtracking line search** (up to `max_LS_tries = 10`):
!>         halve \f$ \alpha \f$ until \f$\phi(\alpha)\le\phi(0)+c\,\alpha\,\phi'(0)\f$.
!>       - If it fails on the Newton direction, retry once with steepest descent;
!>         if that fails too, restore the saved state (no step this iteration).
!>    g. **Stagnation detection:** If residual has not improved by factor
!>       \f$ (1 - \texttt{stag\_rtol}) \f$ for `n_stag_max` consecutive iterations, or no
!>       step could be accepted in step f (the next iterations would be identical),
!>       attempt **re-speciation fallback** (once, same as step 2a) then restart.
!>       If the fallback was already tried or fails, return the best iterate with failure.
!> 3. **Post-processing:** Deallocate temporaries.
!>
!> **Convergence Criteria (machine-epsilon floor, OR arithmetic AND logarithmic):**
!>
!> *Machine-epsilon branch (termination safeguard):*
!> \f[
!>   \|\mathbf{f}\|_\infty \le \varepsilon_{\mathrm{mach}}
!>     \cdot \max\!\bigl(\|\mathbf{U}\mathbf{c}\|_\infty,\,1\bigr)
!> \f]
!>
!> *Arithmetic branch:*
!> \f[
!>   \|\mathbf{f}\|_\infty < \varepsilon_{\mathrm{abs}}
!>   \quad\text{AND}\quad
!>   \|\mathbf{f}\|_\infty < \varepsilon_{\mathrm{rel}}
!>     \cdot \max\!\bigl(\|\mathbf{U}\mathbf{c}\|_\infty,\,\varepsilon_{\mathrm{mach}}\bigr)
!> \f]
!>
!> *Logarithmic branch (required TOGETHER with the arithmetic one):*
!> \f[
!>   \log_{10}\|\mathbf{f}\|_\infty < \varepsilon_{\mathrm{abs}}^{\log}
!>   \quad\text{AND}\quad
!>   \log_{10}\|\mathbf{f}\|_\infty
!>     - \log_{10}\!\max\!\bigl(\|\mathbf{U}\mathbf{c}\|_\infty,\,\varepsilon_{\mathrm{mach}}\bigr)
!>     < \varepsilon_{\mathrm{rel}}^{\log}
!> \f]
!>
!> **Re-speciation Fallback:**
!>
!> When the Newton loop exhausts its iteration budget (`niter_max`) or stagnates
!> (`n_stag_max` consecutive iterations without sufficient improvement, or one iteration in
!> which no step could be accepted), a single
!> re-speciation attempt is made before giving up. This calls
!> `compute_c_nc_from_u_Newton_ideal` with the best-seen primary concentrations
!> \f$\mathbf{c}_1^{\mathrm{best}}\f$ as initial guess and the target components
!> \f$\hat{\mathbf{u}}\f$, effectively restarting Newton from a different basin of
!> attraction. All counters (iteration, stagnation, best-tracking) are
!> reset. The fallback is tried at most once per call to avoid infinite cycling, and only
!> restarts the loop if the inner solve converged; otherwise the best iterate seen so far is
!> returned in `conc_nc` with `CV_flag = .false.`.
!>
!> **Assumptions:**
!> - All primary species are aqueous.
!> - Ideal solution behavior (activity coefficients = 1); the object's \f$\log_{10}\gamma\f$ are
!>   used as fixed values, so a non-ideal object would be treated as fixed-\f$\gamma\f$.
!> - No kinetic reactions — only equilibrium speciation.
!> - Component matrix without constant activity species (De Simoni et al., 2005), so that
!>   \f$n_v - n_p\f$ equals the number of equilibrium reactions.
!>
!> \param[in,out] this     Aqueous chemistry object (class `aqueous_chemistry_c`); on exit its
!>                          secondary variable activity concentrations are those of the last
!>                          speciation performed, and its primary concentrations may have been
!>                          overwritten by the re-speciation fallback (the caller re-synchronises
!>                          the object from `conc_nc`)
!> \param[in]     u_hat    Component concentrations after mixing
!>                          \f$\hat{\mathbf{u}}\f$ [C] (size \f$n_p\f$)
!> \param[in,out] conc_nc  Variable activity species concentrations [C]
!>                          (size \f$n_v\f$; on entry = initial guess,
!>                          on exit = converged solution if CV_flag is TRUE, otherwise the
!>                          iterate with the lowest \f$\|\mathbf{f}\|_\infty\f$ seen since the last (re)start)
!> \param[out]    niter    Number of Newton iterations performed since the last (re)start [-]
!>                          (reset to zero when the re-speciation fallback restarts the loop;
!>                          equals `niter_max` on the max-iteration failure exit)
!> \param[out]    CV_flag  Convergence flag: `.FALSE.` = did not converge,
!>                          `.TRUE.` = converged
!>
!> \pre  The chemistry object `this` must have valid `solid_chemistry`,
!>       reactive-zone stoichiometry, and speciation algorithm data.
!> \pre  `conc_nc` must be pre-allocated with size \f$n_v\f$ and contain a
!>       physically meaningful initial guess (all positive entries).
!> \post On convergence, `conc_nc` contains the converged species concentrations;
!>       otherwise it contains the iterate with the lowest \f$\|\mathbf{f}\|_\infty\f$ seen since
!>       the last (re)start of the loop.
!>
!> \warning Only valid for ideal solutions (unit activity coefficients).
!> \warning LAPACK `dgesv` is called directly for the column-equilibrated log-space Newton system.
!> \warning A diagnostic line is printed to stdout on the max-iteration failure exit
!>          (the stagnation failure exit is silent).
!>
!> \note  The log-space transformation guarantees positivity of concentrations
!>        by construction: \f$ c_{1,j}^{\mathrm{new}} = 10^{p_j + \alpha\Delta p_j} > 0 \f$.
!> \note  The Jacobian formula divides by \f$c_{1,j}\f$, so the initial guess must be
!>        strictly positive (the log-space update then keeps it so).
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
    real(kind=8), allocatable :: conc_nc_save(:)    !< Saved full concentrations before the line-search trials (n_v) [C]
    real(kind=8), allocatable :: conc_nc_best(:)    !< Best concentrations seen during Newton (lowest residual) (n_v) [C]
!> Local variables — log-space Newton system (hoisted out of Newton loop)
    real(kind=8), allocatable :: log_c1(:)          !< Log10 of primary concentrations \f$p_j = \log_{10}(c_{1,j})\f$ (n_p) [-]
    real(kind=8), allocatable :: Delta_p(:)          !< Search direction in log-space \f$\Delta\mathbf{p}\f$: Newton step or steepest descent \f$-\mathbf{J}^T\mathbf{f}\f$ (n_p) [-]
    real(kind=8), allocatable :: dfk_dp(:,:)         !< Log-space Jacobian \f$\partial\mathbf{f}/\partial\mathbf{p}\f$ where \f$(\cdot)_{:,j} = \ln(10)\,c_{1,j}\,(\partial\mathbf{f}/\partial c_{1,j})\f$ (n_p × n_p) [-]
    real(kind=8), allocatable :: dfk_dp_eq(:,:)      !< Column-equilibrated copy of the log-space Jacobian (overwritten by dgesv) (n_p × n_p) [-]
    real(kind=8), parameter :: ln10 = log(10d0)      !< Natural logarithm of 10, \f$\ln(10) \approx 2.302585\f$ [-]
    real(kind=8), parameter :: max_dp = 10d0         !< Maximum allowed \f$|\Delta p_j|\f$ per iteration (10 orders of magnitude) [-]
    real(kind=8) :: scale_factor                      !< \f$\|\Delta\mathbf{p}\|_\infty\f$ before bounding; the step is multiplied by max_dp/scale_factor when it exceeds max_dp [-]
!> Local variables — search direction and Armijo line search
    integer(kind=4) :: i_dir                          !< Direction attempt counter (at most two: Newton then steepest descent, or steepest descent only when the solve fails) [-]
    logical :: use_sd                                 !< Current direction is steepest descent \f$-\mathbf{J}^T\mathbf{f}\f$
    real(kind=8) :: alpha                             !< Backtracking line search step length \f$\alpha \in (0,1]\f$ [-]
    integer(kind=4) :: i_LS                           !< Line search iteration counter [-]
    integer(kind=4), parameter :: max_LS_tries = 10   !< Maximum backtracking line search steps per direction [-]
    real(kind=8), parameter :: c_armijo = 1d-4        !< Armijo sufficient-decrease constant \f$c\f$ [-]
    logical :: step_accepted                          !< Whether the line search accepted a step
    real(kind=8), allocatable :: fk_trial(:)          !< Trial residual vector for line search evaluation (n_p) [C]
    real(kind=8) :: fk_norm                           !< Current residual infinity norm \f$\|\mathbf{f}\|_\infty\f$ [C]
    real(kind=8) :: phi0                              !< Merit function \f$\phi(0)=\tfrac12\|\mathbf{f}\|_2^2\f$ at the current iterate [C²]
    real(kind=8) :: dphi0                             !< \f$\phi'(0)=d\phi/d\alpha|_{\alpha=0}=\mathbf{f}^T\mathbf{J}\Delta\mathbf{p}\f$: directional derivative of \f$\tfrac12\|\mathbf{f}\|_2^2\f$ along \f$\Delta\mathbf{p}\f$ (must be < 0) [C²]
    real(kind=8) :: phi_trial                         !< Merit function \f$\phi(\alpha)\f$ at the trial point [C²]
!> Local variables — stagnation detection and re-speciation fallback
    integer(kind=4) :: n_stag                       !< Consecutive iterations without sufficient residual improvement; set to n_stag_max at once when no step could be accepted [-]
    integer(kind=4), parameter :: n_stag_max = 20    !< Maximum stagnation iterations before early exit or re-speciation fallback [-]
    real(kind=8) :: fk_norm_prev                      !< Reference residual norm for the stagnation test (last value that improved by at least stag_rtol) [C]
    real(kind=8), parameter :: stag_rtol = 1d-3      !< Relative improvement threshold; stagnation counter resets if \f$\|\mathbf{f}\| < (1-\texttt{stag\_rtol})\,\|\mathbf{f}\|_{\mathrm{prev}}\f$ [-]
    real(kind=8) :: best_fk_norm                      !< Lowest \f$\|\mathbf{f}\|_\infty\f$ seen since the last (re)start of the Newton loop [C]
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
    real(kind=8) :: u_scale                          !< Relative-tolerance scale \f$\max(\|\mathbf{U}\mathbf{c}\|_\infty,\,\varepsilon_{\mathrm{mach}})\f$ [C]
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
    comp_mat = this%solid_chemistry%reactive_zone%speciation_alg%comp_mat !< Cache component matrix U
    n_eq = n_v - n_p                 !< Number of equilibrium reactions
    !> Cache component-matrix sub-blocks (avoids creating temporary array slices each iteration)
    U1 = comp_mat(:,1:n_p)          !< Primary species block U1 (n_p × n_p)
    U2 = comp_mat(:,n_p+1:n_v)      !< Secondary species block U2 (n_p × n_eq)
    !> Allocate persistent work arrays used throughout the Newton loop
    allocate(dfk_dc1(n_p,n_p), dc2v_dc1(n_eq, n_p), fk(n_p), conc_nc_best(n_v), conc_nc_save(n_v))  !< Jacobian, residual, best/save concentration arrays
    !> Allocate log-space work arrays once (hoisted out of the Newton loop)
    allocate(dfk_dp(n_p,n_p), dfk_dp_eq(n_p,n_p), Delta_p(n_p), log_c1(n_p), Uc(n_p))  !< Log-Jacobian, equilibrated copy, step, log-conc, cached U*c
    allocate(fk_trial(n_p), ipiv_dgesv(n_p))  !< Trial residual vector and LAPACK pivot array
    allocate(col_scale(n_p))  !< Column scaling factors for Jacobian equilibration
    n_stag = 0                      !< Reset stagnation counter
    fk_norm_prev = huge(1d0)        !< Initialise previous residual norm to +∞ (any first residual will be an improvement)
    best_fk_norm = huge(1d0)        !< Initialise best residual norm to +∞
    predictor_applied = .false.     !< Re-speciation fallback has not been tried yet
    skip_eval = .false.             !< No accepted step yet; must compute from scratch
    log_gamma_var_act = this%get_log_gamma_var_act()  !< Read log10 activity coefficients from the object, held fixed during the iteration (ideal ⇒ all zeros)
!> -----------------------------------------------------------------------
!> \section newton_eq_process Process: Newton loop
!> -----------------------------------------------------------------------

    !> -----------------------------------------------------------------------
    !> \subsection newton_eq_loop Newton iteration loop
    !> -----------------------------------------------------------------------
    newton_loop: do                    !< Begin main Newton iteration loop
        niter = niter + 1              !< Increment Newton iteration counter
    !> \subsubsection newton_eq_maxiter Maximum iteration check and re-speciation fallback
    !> If niter exceeds niter_max, the loop has exhausted its budget.
    !> Before giving up, try re-speciation fallback (once).
        if (niter > niter_max) then     !< Check whether iteration budget is exhausted
            if (.not. predictor_applied) then  !< First exhaustion: try re-speciation fallback
                !> Re-speciation fallback: call the inner speciation Newton solver
                !> with the best-seen c1 as initial guess and u_hat as component concentrations.
                !> This effectively restarts from a cleaner speciation state.
                predictor_applied = .true.                  !< Mark fallback as used (only one attempt allowed)
                call this%compute_c_nc_from_u_Newton_ideal(conc_nc_best(1:n_p), u_hat, conc_nc, niter_spec, CV_flag_spec)  !< Re-speciate from best c1 to obtain new initial guess
                if (CV_flag_spec) then                      !< Restart only from a converged inner solve
                    niter = 0                               !< Reset outer Newton counter for fresh restart
                    n_stag = 0                              !< Reset stagnation counter
                    fk_norm_prev = huge(1d0)                !< Reset previous residual norm to +∞
                    best_fk_norm = huge(1d0)                !< Reset best residual norm to +∞
                    skip_eval = .false.                     !< Force re-evaluation after fallback restart
                    cycle newton_loop                       !< Restart Newton loop with new initial guess
                end if                                      !< End converged-inner-solve restart
            end if                                          !< End re-speciation fallback block
            conc_nc = conc_nc_best      !< Return the best iterate seen
            niter = niter_max           !< This pass only ran the budget check: report the iterations actually performed
            !> Report failure: solver exhausted iteration budget without reaching tolerance
            print *, 'Newton_EI_eq_ideal did not converge after', &  !< Diagnostic message to stdout
                niter_max, 'iterations. Best ||fk|| =', best_fk_norm !< Print iteration count and best residual
            exit newton_loop            !< Leave the Newton loop with CV_flag = .false.
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
    !> The best \f$\mathbf{c}_1\f$ is the initial guess of the re-speciation fallback and the
    !> value returned in conc_nc when the solver fails.
        if (fk_norm < best_fk_norm) then  !< New lowest residual found — update best tracking
            best_fk_norm = fk_norm      !< Update best residual norm
            conc_nc_best = conc_nc      !< Snapshot concentrations at the new best point
        end if                          !< End best-solution tracking update
    !> \subsubsection newton_eq_conv Convergence check (machine-epsilon floor OR arithmetic AND logarithmic)
    !> The Newton iterate is accepted as converged only if
    !>
    !> **Branch 1 — Machine-epsilon safeguard (termination floor):**
    !> \f$\|\mathbf{f}\|_\infty \le \varepsilon_{\mathrm{mach}}\cdot\max(\|\mathbf{U}\mathbf{c}\|,1)\f$
    !>
    !> or BOTH of
    !>
    !> **Branch 2 — Arithmetic tolerances (abs AND rel):**
    !> \f$\|\mathbf{f}\|_\infty < \varepsilon_{\mathrm{abs}}\f$ AND
    !> \f$\|\mathbf{f}\|_\infty < \varepsilon_{\mathrm{rel}}\cdot\max(\|\mathbf{U}\mathbf{c}\|,\varepsilon_{\mathrm{mach}})\f$
    !>
    !> **Branch 3 — Logarithmic tolerances (log-abs AND log-rel):**
    !> \f$\log_{10}\|\mathbf{f}\|_\infty < \varepsilon_{\mathrm{abs}}^{\log}\f$ AND
    !> \f$\log_{10}\|\mathbf{f}\|_\infty - \log_{10}\max(\|\mathbf{U}\mathbf{c}\|,\varepsilon_{\mathrm{mach}}) < \varepsilon_{\mathrm{rel}}^{\log}\f$
    !>
    !> The relative tests scale with the ACTUAL component norm (not max(.,1)), so for
    !> component norms below unity they are tighter than a max(||Uc||,1) scaling would be;
    !> Branch 1 caps the achievable strictness at machine precision so the loop always terminates.
        u_scale = max(u_norm, eps_d)    !< Relative-tolerance scale (never zero)
        if ((fk_norm <= eps_d * max(u_norm, 1d0)) .or. &                               !< Branch 1: machine-epsilon floor
            ((fk_norm < abs_tol .and. &                                                 !< Branch 2a: absolute tolerance
              fk_norm < rel_tol * u_scale) .and. &                                      !< Branch 2b: relative tolerance
             (fk_norm > 0d0 .and. &                                                     !< Branch 3: guard against log10(0)
              log10(fk_norm) < log_abs_tol .and. &                                      !< Branch 3a: log-absolute tolerance
              (log10(fk_norm) - log10(u_scale)) < log_rel_tol))) then                   !< Branch 3b: log-relative tolerance
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
        !> \subsubsection newton_eq_step Newton direction and Armijo backtracking line search in log-space
        !> Solve the column-equilibrated Newton system; if the solve fails, the direction is not a
        !> descent direction of \f$\phi=\tfrac12\|\mathbf{f}\|_2^2\f$, or its line search fails, retry once
        !> with steepest descent \f$-\mathbf{J}^T\mathbf{f}\f$.
            conc_nc_save = conc_nc      !< Snapshot concentrations before any trial steps
            step_accepted = .false.     !< No step accepted yet
            !> Column equilibration of the log-Jacobian (improves conditioning when concentrations span many orders of magnitude)
            dfk_dp_eq = dfk_dp          !< Copy: dgesv overwrites it with the LU factors
            do ii = 1, n_p              !< Compute column scaling factors
                col_scale(ii) = maxval(abs(dfk_dp_eq(:,ii)))  !< Infinity norm of column
                if (col_scale(ii) < tiny(1d0)) col_scale(ii) = 1d0  !< Guard against zero columns
                dfk_dp_eq(:,ii) = dfk_dp_eq(:,ii) / col_scale(ii)  !< Normalise column
            end do                      !< End column equilibration loop
            Delta_p = -fk               !< RHS = \f$-\mathbf{f}\f$ (overwritten with the solution on exit)
            call dgesv(n_p, 1, dfk_dp_eq, n_p, ipiv_dgesv, Delta_p, n_p, info_dgesv)  !< LAPACK LU solve of \f$\mathbf{J}\Delta\mathbf{p} = -\mathbf{f}\f$
            sing_flag = (info_dgesv /= 0)  !< dgesv returned an error code
            if (.not. sing_flag) then   !< Only unscale/scan if dgesv succeeded
                Delta_p = Delta_p / col_scale  !< Recover \f$\Delta\mathbf{p}\f$ of the unscaled system
                sing_flag = any(Delta_p /= Delta_p) .or. maxval(abs(Delta_p)) > 1d20  !< NaN (x/=x) + overflow check
            end if                      !< End dgesv-success guard
            use_sd = sing_flag          !< Unusable Newton step: start directly with steepest descent
            if (use_sd) Delta_p = -matmul(transpose(dfk_dp), fk)  !< \f$-\mathbf{J}^T\mathbf{f}\f$
            dir_loop: do i_dir = 1, 2   !< At most two directions: Newton, then steepest descent
                dphi0 = dot_product(fk, matmul(dfk_dp, Delta_p))  !< \f$\phi'(0)=d\phi/d\alpha|_{\alpha=0}=\mathbf{f}^T\mathbf{J}\Delta\mathbf{p}\f$ (descent iff < 0)
                if (dphi0 >= 0d0) then  !< Not a descent direction: Armijo cannot be satisfied
                    if (use_sd) exit dir_loop  !< \f$\mathbf{J}^T\mathbf{f} = 0\f$ to round-off: stationary point of \f$\phi\f$, no step possible
                    use_sd = .true.     !< Switch to steepest descent
                    Delta_p = -matmul(transpose(dfk_dp), fk)  !< \f$-\mathbf{J}^T\mathbf{f}\f$
                    cycle dir_loop      !< Re-check descent with the new direction
                end if                  !< End descent check
                !> Uniform step bounding: scale \f$\Delta\mathbf{p}\f$ so that \f$\|\Delta\mathbf{p}\|_\infty \le \texttt{max\_dp}\f$ (preserves direction)
                scale_factor = maxval(abs(Delta_p))  !< Find largest absolute step component
                if (scale_factor > max_dp) Delta_p = Delta_p * (max_dp / scale_factor)  !< Uniform scaling to enforce max_dp
                phi0 = 5d-1*dot_product(fk, fk)  !< \f$\phi(0)\f$
                dphi0 = dot_product(fk, matmul(dfk_dp, Delta_p))  !< \f$\phi'(0)=d\phi/d\alpha|_{\alpha=0}\f$ of the (possibly rescaled) direction actually used
                !> \paragraph newton_eq_ls Armijo backtracking line search in log-space
                !> Starting from the full step (\f$\alpha=1\f$), halve \f$\alpha\f$ until
                !> \f$\phi(\alpha) \le \phi(0) + c\,\alpha\,\phi'(0)\f$.
                alpha = 1d0             !< Initial step length = full step
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
                    phi_trial = 5d-1*dot_product(fk_trial, fk_trial)  !< \f$\phi(\alpha)\f$
                    if (phi_trial <= phi0 + c_armijo*alpha*dphi0) then  !< Armijo sufficient decrease
                        step_accepted = .true.  !< Step is accepted
                        skip_eval = .true.       !< Signal next iteration to skip speciation recomputation
                        exit            !< Exit line-search loop — proceed to next Newton iteration
                    end if              !< End Armijo acceptance
                    alpha = alpha * 0.5d0  !< Halve step length and retry
                end do                  !< End backtracking line-search loop
                if (step_accepted) exit dir_loop  !< Done for this Newton iteration
                conc_nc = conc_nc_save  !< Line search failed: restore the pre-step state
                if (use_sd) exit dir_loop  !< Steepest descent failed too: no step this iteration
                use_sd = .true.         !< Retry once with steepest descent
                Delta_p = -matmul(transpose(dfk_dp), fk)  !< \f$-\mathbf{J}^T\mathbf{f}\f$
            end do dir_loop             !< End direction loop
        !> \subsubsection newton_eq_stag Stagnation detection and re-speciation fallback
        !> If the residual norm has decreased by at least the factor \f$(1 - \texttt{stag\_rtol})\f$
        !> compared to the reference, reset the stagnation counter and update the reference norm.
        !> Otherwise increment the counter. A rejected step counts as full stagnation at once: the
        !> state is unchanged, so every following iteration would repeat this one exactly.
            if (.not. step_accepted) then  !< No step from this point
                n_stag = n_stag_max     !< Go straight to the fallback / exit
            else if (fk_norm < fk_norm_prev * (1d0 - stag_rtol)) then  !< Sufficient improvement detected
                n_stag = 0              !< Reset stagnation counter
                fk_norm_prev = fk_norm  !< Update reference norm for next comparison
            else                        !< Residual did not improve sufficiently
                n_stag = n_stag + 1     !< Increment stagnation counter
            end if                      !< End stagnation bookkeeping
            if (n_stag >= n_stag_max) then  !< Stagnation limit reached
                if (.not. predictor_applied) then  !< Has the re-speciation fallback been tried yet?
                    !> Stagnation detected: try re-speciation from best-seen c1 before giving up.
                    predictor_applied = .true.           !< Mark fallback as used (only one attempt allowed)
                    call this%compute_c_nc_from_u_Newton_ideal(conc_nc_best(1:n_p), u_hat, conc_nc, niter_spec, CV_flag_spec)  !< Re-speciate from best c1 (stagnation restart)
                    if (CV_flag_spec) then               !< Restart only from a converged inner solve
                        niter = 0                        !< Reset outer Newton counter for fresh restart
                        n_stag = 0                       !< Reset stagnation counter
                        fk_norm_prev = huge(1d0)         !< Reset previous residual norm to +∞
                        best_fk_norm = huge(1d0)         !< Reset best residual norm to +∞
                        skip_eval = .false.              !< Force re-evaluation after fallback restart
                        cycle newton_loop                !< Restart Newton loop with new initial guess
                    end if                               !< End converged-inner-solve restart
                end if                  !< End stagnation re-speciation fallback
                conc_nc = conc_nc_best  !< Return the best iterate seen
                exit newton_loop        !< Stagnated: leave with CV_flag = .false.
            end if                      !< End stagnation-limit handler (n_stag >= n_stag_max)
        end if                          !< End convergence else-branch (Jacobian + line search + stagnation)
    end do newton_loop                  !< End of main Newton iteration loop
!> -----------------------------------------------------------------------
!> \section newton_eq_postproc Post-processing: cleanup and deallocation
!> -----------------------------------------------------------------------
    deallocate(fk, dfk_dc1, dc2v_dc1, conc_nc_best, conc_nc_save)  !< Free Jacobian, residual, best/save arrays
    deallocate(dfk_dp, dfk_dp_eq, Delta_p, log_c1, fk_trial, ipiv_dgesv, Uc, col_scale)  !< Free log-space work arrays
    deallocate(comp_mat, U1, U2)        !< Free cached component-matrix blocks
end subroutine Newton_EI_eq_ideal       !< End of Newton_EI_eq_ideal subroutine