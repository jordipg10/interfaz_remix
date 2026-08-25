!> \file Newton_EI_eq_kin_anal_ideal_opt2.f90
!> \brief Newton method for WMA reactive mixing iteration with Euler implicit kinetic reactions.
!> We assume the aqueous chemistry has both equilibrium and kinetic reactions, and an ideal solution.
!> Analytical Jacobians are used for efficiency, and option 2 is applied for averaging kinetic rates.
!>
!> \details
!> This subroutine implements a Newton method to solve the nonlinear chemical system arising
!> from a WMA reactive mixing iteration with equilibrium and kinetic reactions. The method uses:
!> - **Euler implicit** time integration for kinetic reactions
!> - **Analytical Jacobians** for efficiency
!> - **Ideal solution** assumption (unit activity coefficients)
!> - **Option 2** for kinetic reaction rates time integration (weighted average of old and new rates)
!> - **Log-space Newton steps** to eliminate the \f$1/c_1\f$ singularity in the Jacobian
!> - **Levenberg-Marquardt damping** with backtracking line search for global convergence
!> - **Stagnation detection** and an **explicit-predictor fallback** for robustness
!>
!> **Mathematical Formulation:**
!>
!> The nonlinear residual to drive to zero is:
!> \f[
!>   \mathbf{f}_k = \mathbf{U}\,\mathbf{c}_{v}
!>                  - \hat{\mathbf{u}}
!>                  - \Delta t \bigl[
!>                      (1-\theta)\,\lambda_{r}^{\mathrm{old}}\,
!>                        (\mathbf{U}\mathbf{S}_k^T)^{\mathrm{old}}\,\mathbf{r}_k^{\mathrm{old}}
!>                    + \theta\,\lambda_{r}^{\mathrm{new}}\,
!>                        (\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}\,\mathbf{r}_k^{\mathrm{new}}
!>                    \bigr]
!>                  = \mathbf{0}
!> \f]
!>
!> where:
!> - \f$ \mathbf{U} \f$ = Component matrix mapping species to components
!>   [\f$n_p \times n_v\f$]
!> - \f$ \mathbf{c}_{v} \f$ = Variable activity species concentrations [C]
!>   (\f$n_v = n_p + n_{eq}\f$)
!> - \f$ \hat{\mathbf{u}} \f$ = Component concentrations after mixing [C]
!>   (size \f$n_p\f$)
!> - \f$ \lambda_{r}^{\mathrm{old}},\;\lambda_{r}^{\mathrm{new}} \f$ = Mixing ratios of
!>   kinetic reaction amounts for old and new time levels [-]
!> - \f$ (\mathbf{U}\mathbf{S}_k^T)^{\mathrm{old/new}} \f$ = Precomputed products of
!>   component matrix and transposed kinetic stoichiometry for old/new reactive zones
!>   [\f$ n_p \times n_{kin} \f$]
!> - \f$ \mathbf{r}_k^{\mathrm{old}} \f$ = Kinetic reaction rates from previous time step
!>   [1/T]
!> - \f$ \mathbf{r}_k^{\mathrm{new}} \f$ = Kinetic reaction rates at current Newton iterate
!>   [1/T]
!> - \f$ \theta \f$ = Time weighting factor (0 = explicit, 1 = fully implicit) [-]
!> - \f$ \Delta t \f$ = Time step size [T]
!>
!> **Log-Space Jacobian Transformation:**
!>
!> Defining \f$ p_j = \log_{10}(c_{1,j}) \f$, the Newton system is solved for
!> \f$ \Delta\mathbf{p} \f$ instead of \f$ \Delta\mathbf{c}_1 \f$:
!> \f[
!>   \frac{\partial \mathbf{f}_k}{\partial p_j}
!>     = \ln(10)\,c_{1,j}\,\frac{\partial \mathbf{f}_k}{\partial c_{1,j}}
!> \f]
!> This cancels the \f$ 1/c_1 \f$ singularity in \f$ \partial c_2/\partial c_1 \f$,
!> producing a well-conditioned system even for trace-level species.
!>
!> **Levenberg-Marquardt Damping:**
!>
!> The damped linear system solved at each Newton iteration is:
!> \f[
!>   \left(\frac{\partial \mathbf{f}_k}{\partial \mathbf{p}}
!>         + \lambda\,\mathbf{I}\right) \Delta\mathbf{p} = -\mathbf{f}_k
!> \f]
!> where \f$ \lambda \f$ is adapted: decreased by 10× on accepted steps, increased
!> by 10× on rejected steps, bounded in \f$[10^{-12},\,10^{8}]\f$.
!>
!> **Algorithm:**
!> 1. **Pre-processing:** Allocate arrays, compute ideal activity coefficients,
!>    retrieve old kinetic rates \f$ \mathbf{r}_k^{\mathrm{old}} \f$, save initial guess.
!> 2. **Newton iteration loop:**
!>    a. Check iteration count; if exceeded, try explicit-predictor fallback (once).
!>    b. Set primary species into the chemistry object.
!>    c. Compute secondary species via mass action law (ideal solution).
!>    d. Set ideal aqueous activities.
!>    e. Compute kinetic rates \f$ \mathbf{r}_k^{\mathrm{new}} \f$ and analytical
!>       Jacobian \f$ \partial\mathbf{r}_k/\partial\mathbf{c} \f$.
!>    f. Compute weighted average rates
!>       \f$ \bar{\mathbf{r}}_k = \theta\,\mathbf{r}_k^{\mathrm{new}}
!>                                + (1-\theta)\,\mathbf{r}_k^{\mathrm{old}} \f$.
!>    g. Compute residual \f$ \mathbf{f}_k \f$ and its infinity norm.
!>    h. Track best solution (lowest \f$ \|\mathbf{f}_k\|_\infty \f$).
!>    i. **Convergence check** (dual criteria — see below).
!>    j. If not converged:
!>       - Compute arithmetic Jacobian \f$ \partial\mathbf{f}_k/\partial\mathbf{c}_1 \f$.
!>       - Transform to log-space Jacobian \f$ \partial\mathbf{f}_k/\partial\mathbf{p} \f$.
!>       - **LM inner loop** (up to `max_LM_tries = 10`):
!>         * Form damped system and solve with LAPACK `dgesv`.
!>         * Check for singular / NaN solutions.
!>         * Uniformly scale \f$ \Delta\mathbf{p} \f$ so
!>           \f$ \|\Delta\mathbf{p}\|_\infty \le \f$ `max_dp` (preserves direction).
!>         * **Backtracking line search** (up to `max_LS_tries = 10`):
!>           halve step length \f$ \alpha \f$ until
!>           \f$ \|\mathbf{f}_k^{\mathrm{trial}}\|_\infty < \|\mathbf{f}_k\|_\infty \f$.
!>         * Accept step or increase \f$ \lambda \f$ and retry.
!>       - If all LM attempts fail, restore saved state and reset \f$ \lambda \f$.
!>    k. **Stagnation detection:** If residual has not improved by factor
!>       \f$ (1 - \texttt{stag\_rtol}) \f$ for `n_stag_max` consecutive iterations,
!>       try the explicit-predictor fallback (once) or exit.
!> 3. **Post-processing:** Increment debug counter, deallocate temporaries.
!>
!> **Convergence Criteria (dual test — either branch suffices):**
!>
!> *Arithmetic branch:*
!> \f[
!>   \|\mathbf{f}_k\|_\infty < \varepsilon_{\mathrm{abs}}
!>   \quad\text{AND}\quad
!>   \|\mathbf{f}_k\|_\infty < \varepsilon_{\mathrm{rel}}
!>     \cdot \max\!\bigl(\|\hat{\mathbf{u}}\|_\infty,\,1\bigr)
!> \f]
!>
!> *Logarithmic branch:*
!> \f[
!>   \log_{10}\|\mathbf{f}_k\|_\infty < \varepsilon_{\mathrm{abs}}^{\log}
!>   \quad\text{AND}\quad
!>   \log_{10}\|\mathbf{f}_k\|_\infty
!>     - \log_{10}\!\max\!\bigl(\|\hat{\mathbf{u}}\|_\infty,\,1\bigr)
!>     < \varepsilon_{\mathrm{rel}}^{\log}
!> \f]
!>
!> **Explicit-Predictor Fallback:**
!>
!> When the Newton loop exceeds `niter_max` or stagnation is detected, a forward-Euler
!> predictor is attempted (once):
!> \f[
!>   \mathbf{u}_{\mathrm{pred}} = \hat{\mathbf{u}}
!>     + \Delta t\bigl[
!>         (1-\theta)\,\lambda_{r}^{\mathrm{old}}\,
!>           (\mathbf{U}\mathbf{S}_k^T)^{\mathrm{old}}\,\mathbf{r}_k^{\mathrm{old}}
!>       + \theta\,\lambda_{r}^{\mathrm{new}}\,
!>           (\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}\,\mathbf{r}_k^{\mathrm{old}}
!>       \bigr]
!> \f]
!> followed by a speciation solve to obtain a new initial guess for the Newton loop.
!>
!> **Assumptions:**
!> - Ideal solution behavior (activity coefficients = 1).
!> - Chemical system associated to aqueous chemistry class object has equilibrium and kinetic reactions.
!>
!> \param[in,out] this   Aqueous chemistry object at current time step
!>                        (class `aqueous_chemistry_c`; provides species, reactions,
!>                        stoichiometry, and convergence parameters)
!> \param[in]     u_hat  Component concentrations after mixing \f$\hat{\mathbf{u}}\f$ [C]
!>                        (assumed-shape array, size \f$n_p\f$)
!> \param[in]     mix_ratio_r_old  Mixing ratio of kinetic reaction amounts from the
!>                        previous time level \f$\lambda_{r}^{\mathrm{old}}\f$ [-]
!> \param[in]     mix_ratio_r_new  Mixing ratio of kinetic reaction amounts at the
!>                        current time level \f$\lambda_{r}^{\mathrm{new}}\f$ [-]
!> \param[in]     Delta_t  Time step size \f$\Delta t\f$ for the (k+1)-th time step [T]
!> \param[in]     theta    Reaction time weighting factor \f$\theta \in [0,1]\f$ [-]
!>                          (0 = explicit, 1 = fully implicit)
!> \param[in,out] conc_nc  Variable activity species concentrations \f$\mathbf{c}_{v}\f$ [C]
!>                          (assumed-shape, size \f$n_v\f$; on entry = initial guess,
!>                          on exit = converged solution if CV_flag is TRUE)
!> \param[out]    niter    Number of Newton iterations performed [-]
!> \param[out]    CV_flag  Convergence flag: `.FALSE.` = did not converge,
!>                          `.TRUE.` = converged
!>
!> \pre  The chemistry object `this` must have valid `solid_chemistry`,
!>       `solid_chemistry_old`, reactive-zone stoichiometry, and speciation
!>       algorithm data initialised before entry.
!> \pre  `conc_nc` must be pre-allocated with size \f$n_v\f$ and contain a
!>       physically meaningful initial guess (all positive entries).
!> \post On convergence, `conc_nc` contains the converged species concentrations
!>       and `this` has updated activities, mean kinetic rates, and kinetic
!>       reaction extents (via `compute_Rk`).
!>
!> \warning Only valid for ideal solutions (unit activity coefficients).
!>          For non-ideal systems, use the corresponding non-ideal variant.
!> \warning The module-level `SAVE` variable `dbg_calls` persists across calls
!>          within the same program execution. Diagnostic output is only printed
!>          on the very first call.
!> \warning LAPACK `dgesv` is called directly (not via `LU_lin_syst`) because the
!>          log-space system requires in-place overwrite of the RHS, and the LM
!>          damping modifies the matrix each inner iteration.
!>
!> \note  The log-space transformation \f$ p_j = \log_{10}(c_{1,j}) \f$ is used
!>        because the mass-action Jacobian \f$ \partial c_2/\partial c_1 \f$
!>        contains \f$ 1/c_1 \f$ terms that become singular for trace species.
!>        Multiplying each Jacobian column by \f$ \ln(10)\,c_{1,j} \f$ absorbs
!>        this singularity, yielding a well-conditioned system independent of
!>        species concentration magnitude.
!> \note  The Newton update in concentration space is performed as
!>        \f$ c_{1,j}^{\mathrm{new}} = 10^{p_j + \alpha\,\Delta p_j} \f$,
!>        which guarantees positivity of concentrations by construction.
!>
!> \sa aqueous_chemistry_m, vectors_m::inf_norm_vec_real, metodos_sist_lin_m::LU_lin_syst
!> \sa aqueous_chemistry_c::compute_rk_Jac_rk_anal, aqueous_chemistry_c::compute_dfk_dc1_EI_ideal
!> \sa aqueous_chemistry_c::compute_c_nc_from_u_Newton_ideal, aqueous_chemistry_c::compute_Rk
!>
!> \author Jordi
!> \date   Unknown
!> \ingroup chemistry

subroutine Newton_EI_eq_kin_anal_ideal_opt2(this,u_hat,mix_ratio_r,Delta_t,theta,conc_nc,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c  !< Import aqueous chemistry class for reactive mixing methods
    use vectors_m, only: inf_norm_vec_real              !< Import infinity norm utility for convergence checks
    !use metodos_sist_lin_m, only: LU_lin_syst           !< Import LU-based linear system solver (unused here; dgesv used instead)
    implicit none                                       !< Enforce explicit variable declarations
!> Arguments
    class(aqueous_chemistry_c) :: this              !< Aqueous chemistry object at current time step [-]
    real(kind=8), intent(in) :: u_hat(:)            !< Component concentrations after mixing (n_p) [C]
    real(kind=8), intent(in) :: mix_ratio_r     !< Mixing ratio of new kinetic reaction amounts \f$\lambda_{r}^{\mathrm{new}}\f$ [-]
    real(kind=8), intent(in) :: Delta_t             !< Time step size \f$\Delta t\f$ for the (k+1)-th time step [T]
    real(kind=8), intent(in) :: theta               !< Reaction time weighting factor \f$\theta \in [0,1]\f$ (0=explicit, 1=fully implicit) [-]
    real(kind=8), intent(inout) :: conc_nc(:)       !< Variable activity species concentrations \f$\mathbf{c}_{v}\f$ (n_v); initial guess on entry, solution on exit [C]
    integer(kind=4), intent(out) :: niter           !< Number of Newton iterations performed [-]
    logical, intent(out) :: CV_flag                 !< Convergence flag: .FALSE. = not converged, .TRUE. = converged
!> Local variables — kinetic rates and Jacobians
    real(kind=8), allocatable :: rk_new(:)          !< New kinetic reaction rates \f$\mathbf{r}_k^{\mathrm{new}}\f$ at current Newton iterate (n_kin) [1/T]
    real(kind=8), allocatable :: drk_dc(:,:)        !< Jacobian \f$\partial\mathbf{r}_k/\partial\mathbf{c}\f$ (n_kin × n_species) [1/(T·C)]
    real(kind=8), allocatable :: dfk_dc1(:,:)       !< Jacobian of residual w.r.t. primary species \f$\partial\mathbf{f}_k/\partial\mathbf{c}_1\f$ (n_p × n_p) [-]
    real(kind=8), allocatable :: rk_old(:)          !< Kinetic reaction rates from previous time step \f$\mathbf{r}_k^{\mathrm{old}}\f$ (n_kin) [1/T]
    real(kind=8), allocatable :: rk_avg(:)          !< Weighted average rates \f$\bar{\mathbf{r}}_k = \theta\,\mathbf{r}_k^{\mathrm{new}} + (1-\theta)\,\mathbf{r}_k^{\mathrm{old}}\f$ (n_kin) [1/T]
!> Local variables — species concentrations and residual
    real(kind=8), allocatable :: c2v(:)             !< Secondary variable activity species concentrations \f$\mathbf{c}_2\f$ (n_eq) [C]
    real(kind=8), allocatable :: Delta_c1(:)        !< Primary species Newton correction \f$\Delta\mathbf{c}_1\f$ (n_p) [C]
    real(kind=8), allocatable :: fk(:)              !< Newton residual vector \f$\mathbf{f}_k\f$ (n_p) [C]
    real(kind=8), allocatable :: conc_nc_save(:)    !< Saved full concentrations before LM/line-search trial (n_v) [C]
    real(kind=8), allocatable :: conc_nc_best(:)    !< Best concentrations seen during Newton (lowest residual) (n_v) [C]
    real(kind=8), allocatable :: log_gamma_var_act(:) !< Log activity coefficients for variable activity species (all zero for ideal solution) [-]
!> Local variables — log-space Newton system
    real(kind=8), allocatable :: log_c1(:)          !< Log10 of primary concentrations \f$p_j = \log_{10}(c_{1,j})\f$ (n_p) [-]
    real(kind=8), allocatable :: Delta_p(:)         !< Newton step in log-space \f$\Delta\mathbf{p}\f$ (n_p) [-]
    real(kind=8), allocatable :: dfk_dp(:,:)        !< Log-space Jacobian \f$\partial\mathbf{f}_k/\partial\mathbf{p}\f$ where \f$(\cdot)_{:,j} = \ln(10)\,c_{1,j}\,(\partial\mathbf{f}_k/\partial c_{1,j})\f$ (n_p × n_p) [-]
    real(kind=8), allocatable :: dfk_dp_LM(:,:)     !< LM-damped log-space Jacobian \f$\partial\mathbf{f}_k/\partial\mathbf{p} + \lambda\mathbf{I}\f$ (n_p × n_p) [-]
    real(kind=8), parameter :: ln10 = log(10d0)     !< Natural logarithm of 10, \f$\ln(10) \approx 2.302585\f$ [-]
    real(kind=8), parameter :: max_dp = 10d0        !< Maximum allowed \f$|\Delta p_j|\f$ per iteration (10 orders of magnitude) [-]
    real(kind=8) :: scale_factor                    !< Uniform scaling factor for \f$\Delta\mathbf{p}\f$ to enforce max_dp while preserving direction [-]
!> Local variables — Levenberg-Marquardt and line search
    real(kind=8) :: lambda_LM                       !< Levenberg-Marquardt damping parameter \f$\lambda\f$, adapted in \f$[10^{-12},\,10^{8}]\f$ [-]
    integer(kind=4) :: i_LM                         !< LM inner loop counter [-]
    integer(kind=4), parameter :: max_LM_tries = 10 !< Maximum LM damping attempts per Newton iteration [-]
    real(kind=8) :: alpha                           !< Backtracking line search step length \f$\alpha \in (0,1]\f$ [-]
    integer(kind=4) :: i_LS                         !< Line search iteration counter [-]
    integer(kind=4), parameter :: max_LS_tries = 10 !< Maximum backtracking line search steps per LM trial [-]
    logical :: step_accepted                        !< Whether a step was accepted in the LM + line search procedure
    real(kind=8), allocatable :: fk_trial(:)        !< Trial residual vector for line search evaluation (n_p) [C]
    real(kind=8), allocatable :: rk_trial(:)        !< Trial kinetic rates for line search evaluation (n_kin) [1/T]
    real(kind=8) :: fk_norm                         !< Current residual infinity norm \f$\|\mathbf{f}_k\|_\infty\f$ [C]
    real(kind=8) :: fk_trial_norm                   !< Trial residual infinity norm after line search step [C]
!> Local variables — stagnation detection and fallback
    integer(kind=4) :: n_stag                       !< Consecutive iterations without sufficient residual improvement [-]
    integer(kind=4), parameter :: n_stag_max = 20   !< Maximum stagnation iterations before early exit or predictor fallback [-]
    real(kind=8) :: fk_norm_prev                    !< Previous best residual norm for stagnation comparison [C]
    real(kind=8), parameter :: stag_rtol = 1d-3     !< Relative improvement threshold; stagnation counter resets if \f$\|\mathbf{f}_k\| < (1-\texttt{stag\_rtol})\,\|\mathbf{f}_k\|_{\mathrm{prev}}\f$ [-]
    real(kind=8) :: best_fk_norm                    !< Lowest \f$\|\mathbf{f}_k\|_\infty\f$ seen during the entire Newton loop [C]
    logical :: predictor_applied                    !< Whether the explicit-predictor fallback has already been tried
    real(kind=8), allocatable :: u_pred(:)          !< Predicted component concentrations from forward-Euler predictor (n_p) [C]
    integer(kind=4) :: niter_spec                   !< Number of inner speciation Newton iterations (predictor fallback) [-]
    logical :: CV_flag_spec                         !< Convergence flag from inner speciation Newton (predictor fallback)
    real(kind=8), save :: predictor_damp = 5d-1      !< Adaptive damping factor for the Euler predictor (persists across calls) [-]
    real(kind=8), parameter :: predictor_damp_min = 1d-2  !< Minimum allowed predictor damping factor [-]
    real(kind=8), parameter :: predictor_damp_max = 1d0   !< Maximum allowed predictor damping factor [-]
    real(kind=8), parameter :: predictor_grow = 1.2d0     !< Growth factor for predictor_damp on success [-]
    real(kind=8), parameter :: predictor_shrink = 0.5d0   !< Shrink factor for predictor_damp on failure [-]
!> Local variables — dimensions and misc
    integer(kind=4) :: n_p                          !< Number of primary (aqueous) species [-]
    integer(kind=4) :: n_v                          !< Number of variable activity species (= n_p + n_eq) [-]
    integer(kind=4) :: i_debug                      !< General-purpose loop index (also used for debug prints)
    logical :: sing_flag                            !< Singularity flag: .TRUE. if dgesv reports failure or NaN/huge values in solution
    integer(kind=4), save :: dbg_calls = 0          !< Persistent call counter limiting diagnostic output to the first invocation
    integer(kind=4), allocatable :: ipiv_dgesv(:)   !< Pivot index array for LAPACK dgesv (n_p)
    integer(kind=4) :: info_dgesv                   !< Return code from LAPACK dgesv (0 = success)
    real(kind=8) :: u_norm                        !< Infinity norm of current component concentrations \f$\|\mathbf{U}\mathbf{c}_{v}\|_\infty\f$ [C]
!> Cached loop-invariant parameters (avoid repeated deep struct dereferences per iteration)
    integer(kind=4) :: niter_max                     !< Maximum Newton iterations (cached from CV_params)
    integer(kind=4) :: n_kin                         !< Number of kinetic reactions (cached)
    real(kind=8) :: abs_tol                          !< Absolute convergence tolerance (cached)
    real(kind=8) :: rel_tol                          !< Relative convergence tolerance (cached)
    real(kind=8) :: log_abs_tol                      !< Log absolute tolerance (cached)
    real(kind=8) :: log_rel_tol                      !< Log relative tolerance (cached)
    real(kind=8) :: eps_d                            !< Machine epsilon (cached)
    real(kind=8) :: sqrt_eps_d                       !< eps^(1/2) — practical conditioning noise floor for double precision
    real(kind=8) :: dt_theta_mr_new                  !< Cached scalar \f$\Delta t\,\theta\,\lambda_{r}^{\mathrm{new}}\f$ [-]
    real(kind=8), allocatable :: comp_mat(:,:)       !< Component matrix U (cached, n_p x n_v)
    real(kind=8), allocatable :: U_SkT_old(:,:)      !< (U*S_k^T)^old product (cached, n_p x n_kin)
    real(kind=8), allocatable :: U_SkT_new(:,:)      !< (U*S_k^T)^new product (cached, n_p x n_kin)
    real(kind=8), allocatable :: source_old(:)       !< Constant kinetic source term: Dt*(1-theta)*mr_old*U_SkT_old*rk_old (cached, n_p)
    real(kind=8), allocatable :: Uc(:)               !< Cached product \f$\mathbf{U}\mathbf{c}_{v}\f$ to avoid redundant matmul (n_p) [C]
    real(kind=8), allocatable :: U_SkT_new_rk_old(:) !< Cached product \f$(\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}\mathbf{r}_k^{\mathrm{old}}\f$ for predictor fallback (n_p) [C]
    logical :: skip_eval                             !< Skip residual/speciation recomputation when accepted step can be reused
    real(kind=8), allocatable :: col_scale(:)        !< Column scaling factors for Jacobian equilibration (n_p) [-]
!> -----------------------------------------------------------------------
!> \section newton_ei_preproc Pre-processing: initialisation and allocation
!> -----------------------------------------------------------------------
    ! print *, "=========================================="
    ! print *, "DEBUG: Entering Newton_EI_eq_kin_anal_ideal_opt2"
    ! print *, "DEBUG: Delta_t           = ", Delta_t
    ! print *, "DEBUG: theta             = ", theta
    ! print *, "DEBUG: mix_ratio_r_old   = ", mix_ratio_r_old
    ! print *, "DEBUG: mix_ratio_r_new   = ", mix_ratio_r_new
    ! print *, "DEBUG: size(u_hat)       = ", size(u_hat)
    ! print *, "DEBUG: u_hat             = ", u_hat
    ! print *, "DEBUG: size(conc_nc)     = ", size(conc_nc)
    ! print *, "DEBUG: conc_nc (input)   = ", conc_nc
    niter=0                                                                  !< Initialise iteration counter
    CV_flag=.false.                                                          !< Initialise convergence flag to not-converged
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species    !< Query number of primary species from speciation algorithm
    n_v=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species !< Query total number of variable activity species (n_p + n_eq)
    !> Cache loop-invariant parameters from deep struct dereferences
    niter_max=this%solid_chemistry%reactive_zone%CV_params%niter_max !< Cache max Newton iterations
    abs_tol=this%solid_chemistry%reactive_zone%CV_params%abs_tol !< Cache absolute tolerance
    rel_tol=this%solid_chemistry%reactive_zone%CV_params%rel_tol !< Cache relative tolerance
    log_abs_tol=this%solid_chemistry%reactive_zone%CV_params%log_abs_tol !< Cache log absolute tolerance
    log_rel_tol=this%solid_chemistry%reactive_zone%CV_params%log_rel_tol !< Cache log relative tolerance
    eps_d=epsilon(1d0)                                           !< Cache machine epsilon
    sqrt_eps_d=sqrt(eps_d)                                      !< Cache eps^(1/2) ≈ 1.49e-8 (practical noise floor)
    n_kin=this%solid_chemistry%mineral_zone%num_minerals_kin+&   !< Cache number of kinetic reactions
        this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
    !> Cache component matrix and U*S_k^T products (used every iteration for residual computation)
    comp_mat=this%solid_chemistry%reactive_zone%speciation_alg%comp_mat !< Cache component matrix U
    U_SkT_old=this%solid_chemistry_old%reactive_zone%U_SkT_prod !< Cache old U*S_k^T product
    U_SkT_new=this%solid_chemistry%reactive_zone%U_SkT_prod     !< Cache new U*S_k^T product
    !> First-call diagnostics: print dimensions and convergence tolerances on the very first invocation
    if (dbg_calls == 0) then                                    !< Guard: only execute on first call
        continue                                                !< Placeholder (debug prints removed)
    end if                                                       !< End first-call diagnostics
    !> Allocate persistent work arrays used throughout the Newton loop
    allocate(dfk_dc1(n_p,n_p),&                                  !< Allocate Jacobian
        drk_dc(n_kin,this%solid_chemistry%reactive_zone%speciation_alg%num_species),&  !< Kinetic-rate Jacobian (n_kin x n_species)
        rk_new(n_kin),fk(n_p),conc_nc_best(n_v))                 !< Rates, residual, best-solution
    !> Allocate log-space work arrays once (hoisted out of the Newton loop)
    allocate(dfk_dp(n_p,n_p), dfk_dp_LM(n_p,n_p), Delta_p(n_p), log_c1(n_p), Uc(n_p))
    allocate(rk_trial(n_kin), fk_trial(n_p))
    allocate(ipiv_dgesv(n_p))                                    !< Allocate pivot array once (reused each LM attempt)
    allocate(col_scale(n_p))                                     !< Allocate column scaling factors for Jacobian equilibration
    drk_dc=0d0         !< Zero-initialise kinetic Jacobian (computed analytically each Newton iteration)
    dfk_dc1=0d0        !< Zero-initialise residual Jacobian
    lambda_LM=1d-3     !< Initial LM damping parameter \f$\lambda_0\f$
    n_stag=0           !< Reset stagnation counter
    fk_norm_prev=huge(1d0) !< Initialise previous-best residual norm to machine-representable maximum
    ! print *, "DEBUG: Arrays allocated successfully"
    ! print *, "DEBUG: size(drk_dc) = ", size(drk_dc,1), " x ", size(drk_dc,2)
    ! print *, "DEBUG: size(dfk_dc1) = ", size(dfk_dc1,1), " x ", size(dfk_dc1,2)
!> -----------------------------------------------------------------------
!> \section newton_ei_process Process: activity coefficients, old rates, and Newton loop
!> -----------------------------------------------------------------------
    !> Compute ideal activity coefficients (all zeros for ideal solution)
        log_gamma_var_act=this%compute_log_act_coeffs_var_act_ideal()
        !print *, "DEBUG: log_gamma_var_act = ", log_gamma_var_act
    !> Retrieve old kinetic reaction rates from the previous time step
        rk_old=this%get_rk_old()
    !> Precompute constant kinetic source term (does not change during Newton iterations)
    !> source_old = Delta_t * (1-theta) * mix_ratio_r * U_SkT_old * rk_old
        source_old = Delta_t*(1d0-theta)*mix_ratio_r*matmul(U_SkT_old, rk_old) !< Constant source term
        dt_theta_mr_new = Delta_t*theta*mix_ratio_r !< Cache scalar product used in residual and predictor computations
        U_SkT_new_rk_old = matmul(U_SkT_new, rk_old)  !< Cache matrix-vector product for predictor fallback
    !> Save original initial guess for recovery via explicit-predictor fallback
        predictor_applied = .false.     !< Explicit-predictor fallback has not been tried yet
        best_fk_norm = huge(1d0)        !< Initialise best residual norm to \f$+\infty\f$
        skip_eval = .false.             !< No accepted step yet; must compute from scratch

    !> -----------------------------------------------------------------------
    !> \subsection newton_ei_loop Newton iteration loop
    !> -----------------------------------------------------------------------
        newton_loop: do
            niter=niter+1               !< Increment Newton iteration counter
            ! if (niter==1 .and. .not. recovery_applied) then
            !     print *, 'CV_params: abs_tol=', this%solid_chemistry%reactive_zone%CV_params%abs_tol, &
            !         ' rel_tol=', this%solid_chemistry%reactive_zone%CV_params%rel_tol, &
            !         ' log_abs_tol=', this%solid_chemistry%reactive_zone%CV_params%log_abs_tol, &
            !         ' log_rel_tol=', this%solid_chemistry%reactive_zone%CV_params%log_rel_tol
            ! end if
        !> \subsubsection newton_ei_maxiter Maximum iteration check and explicit-predictor fallback
            if (niter>niter_max) then  !< Exceeded iteration budget?
                if (.not. predictor_applied) then               !< Has the explicit-predictor fallback been tried yet?
                    !> Explicit predictor fallback: approximate \f$\mathbf{r}_k^{\mathrm{new}} \approx \mathbf{r}_k^{\mathrm{old}}\f$
                    !> (forward Euler) to compute predicted component concentrations, then
                    !> re-speciate and restart the Newton loop from the new initial guess.
                    !> \f$ \mathbf{u}_{\mathrm{pred}} = \hat{\mathbf{u}} + predictor_damp \Delta t\bigl[(1-\theta)\lambda_{r}^{\mathrm{old}}(\mathbf{U}\mathbf{S}_k^T)^{\mathrm{old}} + \theta\lambda_{r}^{\mathrm{new}}(\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}\bigr]\mathbf{r}_k^{\mathrm{old}} \f$
                    predictor_applied = .true.                  !< Mark predictor as used (only one attempt allowed)
                    allocate(u_pred(n_p))                       !< Allocate predicted component vector
                    u_pred = u_hat + predictor_damp * ( &       !< Damped forward-Euler prediction
                        source_old + &                           !< Reuse precomputed constant term
                        dt_theta_mr_new * U_SkT_new_rk_old )     !< Cached \f$\Delta t\theta\lambda_{r}^{\mathrm{new}}(\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}\mathbf{r}_k^{\mathrm{old}}\f$
                    !> Solve speciation to get species concentrations from predicted components
                    call this%compute_c_nc_from_u_Newton_ideal(conc_nc(1:n_p), u_pred, conc_nc, niter_spec, CV_flag_spec)  !< Inner Newton speciation solve
                    deallocate(u_pred)                          !< Free predicted component vector
                    niter = 0                                   !< Reset outer Newton iteration counter for fresh restart
                    n_stag = 0                                  !< Reset stagnation counter
                    fk_norm_prev = huge(1d0)                    !< Reset previous residual norm to \f$+\infty\f$
                    lambda_LM = 1d-3                            !< Reset LM damping parameter to initial value
                    best_fk_norm = huge(1d0)                    !< Reset best residual norm to \f$+\infty\f$
                    skip_eval = .false.                         !< Force re-evaluation after predictor restart
                    cycle newton_loop                           !< Restart Newton loop with new initial guess
                end if                                          !< End predictor fallback block
                !> Accept if best residual is within conditioning noise: \f$\|\mathbf{f}_{\mathrm{best}}\|_\infty \le \sqrt{\varepsilon}\cdot\max(\|\mathbf{U}\mathbf{c}_{v}\|,1)\f$
                if (best_fk_norm<=sqrt_eps_d*max(u_norm,1d0)) then  !< Within conditioning noise?
                    CV_flag=.true.                              !< Accept as converged
                    conc_nc=conc_nc_best                        !< Restore best-seen concentrations
                    call this%set_conc_prim_species(conc_nc(1:n_p))  !< Push best primary concentrations into chemistry object
                    call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p),log_gamma_var_act,conc_nc(n_p+1:))  !< Recompute secondary species at best point
                    call this%set_act_aq_species()               !< Update aqueous activities
                    call this%compute_Rk(theta,Delta_t)          !< Compute kinetic reaction extents at accepted solution
                else                                             !< Best residual exceeds noise threshold
                    print *, 'Newton did not converge after', niter_max, &  !< Report failure
                        'iterations. Best ||fk|| =', best_fk_norm
                end if                                           !< End conditioning-noise acceptance
                exit newton_loop                                 !< Leave Newton loop (converged or failed)
            end if                                               !< End max-iteration guard
        !> \subsubsection newton_ei_speciation Speciation: primary → secondary → activities → rates
        if (skip_eval) then                                     !< Reuse speciation + rates from accepted line-search step
            skip_eval = .false.                                  !< Consume flag; next iteration recomputes unless step accepted again
        else                                                     !< Compute fresh speciation and rates
            call this%set_conc_prim_species(conc_nc(1:n_p))     !< Push current primary concentrations into chemistry object
            call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p),log_gamma_var_act,conc_nc(n_p+1:)) !< Compute secondary species via mass action law (ideal)
            call this%set_act_aq_species()                      !< Set ideal aqueous activities (= concentrations)
            drk_dc=0d0                                          !< Re-zero Jacobian each iteration to clear stale sparse entries
            call this%compute_rk_Jac_rk_anal(rk_new,drk_dc)    !< Compute \f$\mathbf{r}_k^{\mathrm{new}}\f$ and \f$\partial\mathbf{r}_k/\partial\mathbf{c}\f$
        end if                                                   !< End skip_eval check
            ! print *, "DEBUG: rk_new    = ", rk_new
            ! print *, "DEBUG: drk_dc (by rows):"
            ! do i_debug = 1, size(drk_dc, 1)
            !     print *, "  row", i_debug, ":", drk_dc(i_debug, :)
            ! end do
            ! print *, "DEBUG: ||drk_dc||_max = ", maxval(abs(drk_dc))
        !> \subsubsection newton_ei_avg Weighted average of kinetic reaction rates
            rk_avg=theta*rk_new+(1d0-theta)*rk_old               !< \f$\bar{\mathbf{r}}_k = \theta\,\mathbf{r}_k^{\mathrm{new}} + (1-\theta)\,\mathbf{r}_k^{\mathrm{old}}\f$
            ! print *, "DEBUG: rk_avg    = ", rk_avg
            call this%set_rk_mean(rk_avg)                        !< Store weighted average rates in chemistry object
        !> \subsubsection newton_ei_residual Newton residual computation
        !> \f$\mathbf{f}_k = \mathbf{U}\mathbf{c}_{v} - \hat{\mathbf{u}} - \Delta t\bigl[(1-\theta)\lambda_{r}^{\mathrm{old}}(\mathbf{U}\mathbf{S}_k^T)^{\mathrm{old}}\mathbf{r}_k^{\mathrm{old}} + \theta\lambda_{r}^{\mathrm{new}}(\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}\mathbf{r}_k^{\mathrm{new}}\bigr]\f$
            Uc=matmul(comp_mat,conc_nc)                          !< Cache \f$\mathbf{U}\mathbf{c}_{v}\f$ (reused for u_norm)
            fk=Uc-u_hat-&                                        !< \f$\mathbf{U}\mathbf{c}_{v} - \hat{\mathbf{u}}\f$ minus kinetic source terms
                source_old-&                                     !< Precomputed constant: \f$\Delta t(1-\theta)\lambda_{r}^{\mathrm{old}}(\mathbf{U}\mathbf{S}_k^T)^{\mathrm{old}}\mathbf{r}_k^{\mathrm{old}}\f$
                dt_theta_mr_new*matmul(U_SkT_new,rk_new)         !< \f$\Delta t\theta\lambda_{r}^{\mathrm{new}}(\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}\mathbf{r}_k^{\mathrm{new}}\f$
            ! print *, "DEBUG: term1 (comp_mat*conc_nc) = ",matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,conc_nc)
            ! print *, "DEBUG: term2 (u_hat)            = ", u_hat
            ! print *, "DEBUG: term3_old (Dt*(1-th)*mr_old*U_SkT*rk_old) = ", &
            !     Delta_t*(1d0-theta)*mix_ratio_r_old*matmul(this%solid_chemistry_old%reactive_zone%U_SkT_prod,rk_old)
            ! print *, "DEBUG: term3_new (Dt*th*mr_new*U_SkT*rk_new)     = ", &
            !     Delta_t*theta*mix_ratio_r_new*matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,rk_new)
            ! print *, "DEBUG: term3_old + term3_new = ", &
            !     Delta_t*(1d0-theta)*mix_ratio_r_old*matmul(this%solid_chemistry_old%reactive_zone%U_SkT_prod,rk_old) + &
            !     Delta_t*theta*mix_ratio_r_new*matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,rk_new)
            ! print *, "DEBUG: fk (residual)  = ", fk
            ! print *, "DEBUG: ||fk||_inf     = ", inf_norm_vec_real(fk)
        !> \subsubsection newton_ei_conv Convergence check (dual arithmetic OR logarithmic criteria)
            fk_norm = inf_norm_vec_real(fk)                     !< Compute \f$\|\mathbf{f}_k\|_\infty\f$ for convergence and line-search tests
            u_norm = inf_norm_vec_real(Uc)                       !< \f$\|\mathbf{U}\mathbf{c}_{v}\|_\infty\f$ reusing cached product
            !> Per-iteration diagnostics for hard cells (only first 5 + every 10th iteration)
            ! if (niter <= 5 .or. mod(niter,10)==0) then
            !     print *, 'ITER', niter, ' ||fk||=', fk_norm, ' lambda_LM=', lambda_LM, ' c1=', conc_nc(1:n_p)
            ! end if
            !> Track best solution for stagnation recovery
            if (fk_norm < best_fk_norm) then                   !< Is current residual the lowest seen so far?
                best_fk_norm = fk_norm                           !< Update best residual norm
                conc_nc_best = conc_nc                           !< Snapshot concentrations at the new best point
            end if                                               !< End best-solution tracking
            !> Machine-epsilon safeguard: residual at round-off level cannot be reduced further.
            !> Triple convergence check — any branch suffices:
            !>  Branch 1: \f$\|\mathbf{f}_k\|_\infty \le \varepsilon_{\mathrm{mach}}\cdot\max(\|\hat{\mathbf{u}}\|,1)\f$
            !>  Branch 2: \f$\|\mathbf{f}_k\|_\infty < \varepsilon_{\mathrm{abs}}\f$ AND \f$< \varepsilon_{\mathrm{rel}}\cdot\max(\|\hat{\mathbf{u}}\|,1)\f$
            !>  Branch 3: \f$\log_{10}\|\mathbf{f}_k\|_\infty < \varepsilon_{\mathrm{abs}}^{\log}\f$ AND log-relative test
            if ((fk_norm<=eps_d*max(u_norm,1d0)) .or. &          !< Branch 1: machine-epsilon safeguard
                (fk_norm<abs_tol .and. &                          !< Branch 2a: absolute tolerance
                fk_norm<rel_tol*max(u_norm,1d0)) .or. &           !< Branch 2b: relative tolerance
                (fk_norm>0d0 .and. &                              !< Branch 3: guard against log10(0)
                (log10(fk_norm)<log_abs_tol .and. &               !< Branch 3a: log-absolute tolerance
                (log10(fk_norm)-log10(max(u_norm,1d0)))< &        !< Branch 3b: log-relative tolerance
                log_rel_tol))) then
                CV_flag=.true.                                   !< Mark as converged
                ! print *, "DEBUG: ||fk||_inf = ", fk_norm, " abs_tol = ", &
                !     this%solid_chemistry%reactive_zone%CV_params%abs_tol, " rel_tol*||u_hat|| = ", &
                !     this%solid_chemistry%reactive_zone%CV_params%rel_tol*max(inf_norm_vec_real(u_hat),1d0), &
                !     " log_abs_tol = ", this%solid_chemistry%reactive_zone%CV_params%log_abs_tol, &
                !     " log_rel_tol = ", this%solid_chemistry%reactive_zone%CV_params%log_rel_tol
                ! print *, "DEBUG: Final conc_nc = ", conc_nc
                call this%compute_Rk(theta,Delta_t)             !< Compute kinetic reaction extents at converged solution
                exit newton_loop                                 !< Leave Newton loop — converged
            else
            !> \subsubsection newton_ei_jacobian Jacobian assembly and log-space transformation
                !> Compute the arithmetic Jacobian \f$\partial\mathbf{f}_k/\partial\mathbf{c}_1\f$ including kinetic-rate derivatives
                !> Pass the full chem_syst-ordered drk_dc (NOT a slice to 1:n_v): the
                !> callee permutes columns via reactive_zone%ind_var_act_species, whose
                !> values can exceed n_v_zone in multi-zone setups where the zone is a
                !> strict subset of the chem_syst variable-activity species.
                call this%compute_dfk_dc1_EI_ideal(conc_nc(1:n_p),conc_nc(n_p+1:),drk_dc,Delta_t,theta,mix_ratio_r,&  !< Assembly: \f$\mathbf{U}_1 + \mathbf{U}_2\partial\mathbf{c}_2/\partial\mathbf{c}_1 - \Delta t\theta\lambda_{r}^{\mathrm{new}}(\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}\partial\mathbf{r}_k/\partial\mathbf{c}_1\f$
                    dfk_dc1)                                     !< Output: arithmetic Jacobian (n_p \f$\times\f$ n_p)
                !> Transform to log-space Jacobian:
                !> \f$ (\partial\mathbf{f}_k/\partial\mathbf{p})_{:,j} = \ln(10)\,c_{1,j}\,(\partial\mathbf{f}_k/\partial c_{1,j}) \f$
                do i_debug = 1, n_p                              !< Loop over each primary species \f$j = 1, \ldots, n_p\f$
                    log_c1(i_debug) = log10(max(conc_nc(i_debug), tiny(1d0)))  !< \f$p_j = \log_{10}(\max(c_{1,j},\,\texttt{tiny}))\f$ (floor prevents \f$-\infty\f$)
                    dfk_dp(:,i_debug) = ln10 * conc_nc(i_debug) * dfk_dc1(:,i_debug)  !< Column \f$j\f$ of log-Jacobian: \f$\ln(10)\,c_{1,j}\,(\partial\mathbf{f}_k/\partial c_{1,j})\f$
                end do                                           !< End log-space transformation loop
            !> \subsubsection newton_ei_lm Levenberg-Marquardt with backtracking line search in log-space
                conc_nc_save = conc_nc                           !< Snapshot concentrations before any LM trial steps
                step_accepted = .false.                          !< No step accepted yet
                do i_LM = 1, max_LM_tries                       !< LM damping loop: try up to max_LM_tries values of \f$\lambda\f$
                    !> Column equilibration: normalise each column of the log-Jacobian before
                    !> adding LM damping.  This dramatically improves conditioning when species
                    !> concentrations span many orders of magnitude.
                    dfk_dp_LM = dfk_dp                           !< Fresh copy of log-Jacobian (dgesv overwrites it)
                    do i_debug = 1, n_p                          !< Compute column scaling factors
                        col_scale(i_debug) = maxval(abs(dfk_dp_LM(:,i_debug)))  !< Infinity norm of column
                        if (col_scale(i_debug) < tiny(1d0)) col_scale(i_debug) = 1d0  !< Guard against zero columns
                        dfk_dp_LM(:,i_debug) = dfk_dp_LM(:,i_debug) / col_scale(i_debug)  !< Normalise column
                    end do                                       !< End column equilibration loop
                    !> Apply Levenberg damping in the equilibrated space
                    do i_debug = 1, n_p                          !< Add damping to each diagonal element
                        dfk_dp_LM(i_debug, i_debug) = dfk_dp_LM(i_debug, i_debug) + lambda_LM  !< \f$J_{jj} \leftarrow J_{jj} + \lambda\f$
                    end do                                       !< End diagonal damping loop
                    !> Solve \f$(\partial\mathbf{f}_k/\partial\mathbf{p} + \lambda\mathbf{I})\,\Delta\mathbf{p} = -\mathbf{f}_k\f$
                    !> using LAPACK dgesv (LU with partial pivoting).
                    !> \note dgesv is used instead of LU_lin_syst because (a) the RHS is overwritten
                    !>       in-place, (b) the matrix is rebuilt each LM inner iteration, and (c) LM
                    !>       damping already regularises the system so no additional scaling is needed.
                    Delta_p = -fk                               !< Set RHS = \f$-\mathbf{f}_k\f$ (overwritten with \f$\Delta\mathbf{p}\f$ on exit)
                    call dgesv(n_p, 1, dfk_dp_LM, n_p, ipiv_dgesv, Delta_p, n_p, info_dgesv)  !< LAPACK LU direct solve (ipiv reused each attempt)
                    !> Unscale the solution: recover original Delta_p from column-equilibrated system
                    if (info_dgesv == 0) then
                        do i_debug = 1, n_p
                            Delta_p(i_debug) = Delta_p(i_debug) / col_scale(i_debug)
                        end do
                    end if
                    !> Validate the solution: check for dgesv failure (info /= 0) and NaN / overflow
                    sing_flag = (info_dgesv /= 0)                !< dgesv returned an error code
                    do i_debug = 1, n_p                          !< Scan each component of \f$\Delta\mathbf{p}\f$
                        if (Delta_p(i_debug) /= Delta_p(i_debug) .or. abs(Delta_p(i_debug)) > 1d20) then  !< NaN check (x/=x) or overflow guard
                            sing_flag = .true.                   !< Mark solution as unusable
                            exit                                 !< No need to check remaining components
                        end if
                    end do                                       !< End NaN/overflow scan
                    if (sing_flag) then                          !< Solution is singular or contains garbage values
                        lambda_LM = min(lambda_LM * 10d0, 1d8)  !< Increase damping by 10\f$\times\f$ (capped at \f$10^8\f$)
                        cycle                                    !< Retry with stronger damping
                    end if
                    !> Uniform step bounding: scale \f$\Delta\mathbf{p}\f$ so that
                    !> \f$\|\Delta\mathbf{p}\|_\infty \le \texttt{max\_dp}\f$.
                    !> This preserves the Newton direction while preventing excessively large steps.
                    scale_factor = 1d0                           !< Start with full step
                    do i_debug = 1, n_p                          !< Find the tightest bound across all components
                        if (abs(Delta_p(i_debug)) > 0d0) then   !< Skip zero components
                            scale_factor = min(scale_factor, max_dp / abs(Delta_p(i_debug)))  !< Ratio needed to clamp component
                        end if
                    end do                                       !< End step-bound loop
                    if (scale_factor < 1d0) Delta_p = Delta_p * scale_factor  !< Apply uniform scaling (only if needed)
                    !> \paragraph newton_ei_ls Backtracking line search in log-space
                    !> Try \f$\alpha = 1, 0.5, 0.25, \ldots\f$ until
                    !> \f$\|\mathbf{f}_k^{\mathrm{trial}}\|_\infty < \|\mathbf{f}_k\|_\infty\f$.
                    alpha = 1d0                                  !< Initial step length = full Newton step
                    do i_LS = 1, max_LS_tries                    !< Backtracking loop: halve \f$\alpha\f$ each attempt
                        !> Update primary species in log-space: \f$c_{1,j}^{\mathrm{trial}} = 10^{p_j + \alpha\,\Delta p_j}\f$ (guarantees positivity)
                        do i_debug = 1, n_p                      !< Apply log-space update to each primary species
                            conc_nc(i_debug) = 10d0**(log_c1(i_debug) + alpha * Delta_p(i_debug))  !< \f$c_{1,j} = 10^{p_j + \alpha\Delta p_j}\f$
                        end do                                   !< End primary-species update
                        call this%set_conc_prim_species(conc_nc(1:n_p))  !< Push trial primary concentrations into chemistry object
                        call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p), log_gamma_var_act, conc_nc(n_p+1:))  !< Speciation at trial point
                        call this%set_act_aq_species()            !< Update aqueous activities at trial point
                        drk_dc = 0d0                             !< Re-zero kinetic Jacobian before recomputation
                        call this%compute_rk_Jac_rk_anal(rk_trial, drk_dc)  !< Compute trial kinetic rates \f$\mathbf{r}_k^{\mathrm{trial}}\f$ and Jacobian
                        !> Evaluate trial residual using cached matrices and precomputed constant source term
                        fk_trial = matmul(comp_mat, conc_nc) - u_hat - &  !< \f$\mathbf{U}\mathbf{c}_{v}^{\mathrm{trial}} - \hat{\mathbf{u}}\f$
                            source_old - &                       !< Precomputed constant old kinetic source term
                            dt_theta_mr_new*matmul(U_SkT_new, rk_trial)  !< \f$\Delta t\theta\lambda_{r}^{\mathrm{new}}(\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}\mathbf{r}_k^{\mathrm{trial}}\f$
                        fk_trial_norm = inf_norm_vec_real(fk_trial)  !< \f$\|\mathbf{f}_k^{\mathrm{trial}}\|_\infty\f$
                        if (fk_trial_norm < fk_norm) then        !< Sufficient decrease: trial residual < current residual?
                            step_accepted = .true.               !< Step is accepted
                            rk_new = rk_trial                    !< Cache accepted trial rates for next iteration
                            fk = fk_trial                        !< Cache accepted trial residual for next iteration
                            skip_eval = .true.                   !< Signal next iteration to skip speciation/rates recomputation
                            lambda_LM = max(lambda_LM / 10d0, 1d-12)  !< Decrease damping by 10\f$\times\f$ (floored at \f$10^{-12}\f$) for next iteration
                            exit !> Accept: exit line search     !< Leave backtracking loop
                        end if
                        alpha = alpha * 0.5d0                    !< Halve step length and retry
                    end do                                       !< End backtracking line search loop
                    if (step_accepted) exit !> Exit LM loop on accepted step  !< Leave LM loop if line search found an acceptable step
                    !> Line search failed for this \f$\lambda\f$: restore concentrations and increase damping
                    conc_nc = conc_nc_save                       !< Restore pre-LM concentrations
                    lambda_LM = min(lambda_LM * 10d0, 1d8)      !< Increase damping by 10\f$\times\f$ for next LM attempt
                end do                                           !< End LM damping loop
                !> All LM attempts exhausted: restore saved state and reset \f$\lambda\f$ to initial value
                if (.not. step_accepted) then                    !< All max_LM_tries failed?
                    conc_nc = conc_nc_save                       !< Fall back to the state before this Newton step
                    lambda_LM = 1d-3                             !< Reset damping parameter to default
                end if                                           !< End LM-failure recovery
            !> \subsubsection newton_ei_stag Stagnation detection and predictor fallback
            !> If the residual norm decreased by at least \f$(1 - \texttt{stag\_rtol})\f$ compared to the
            !> previous best, reset the stagnation counter. Otherwise, increment it.
                if (fk_norm < fk_norm_prev * (1d0 - stag_rtol)) then  !< Sufficient improvement detected?
                    n_stag = 0                                   !< Reset stagnation counter
                    fk_norm_prev = fk_norm                       !< Update reference norm for next comparison
                else                                             !< Residual did not improve sufficiently
                    n_stag = n_stag + 1                          !< Increment stagnation counter
                    if (n_stag >= n_stag_max) then               !< Stagnation limit reached?
                        if (.not. predictor_applied) then        !< Has the explicit-predictor fallback been tried yet?
                            !> Stagnation detected: try explicit predictor before giving up.
                            !> Same forward-Euler predictor as the max-iteration fallback.
                            predictor_applied = .true.           !< Mark predictor as used (only one attempt allowed)
                            allocate(u_pred(n_p))                !< Allocate predicted component vector
                            u_pred = u_hat + predictor_damp * ( &  !< Damped forward-Euler prediction
                                source_old + &                   !< Reuse precomputed constant term
                                dt_theta_mr_new * U_SkT_new_rk_old )  !< Cached predictor kinetic source term
                            call this%compute_c_nc_from_u_Newton_ideal(conc_nc(1:n_p), u_pred, conc_nc, niter_spec, CV_flag_spec)  !< Inner speciation solve
                            deallocate(u_pred)                   !< Free predicted component vector
                            niter = 0                            !< Reset outer Newton counter for fresh restart
                            n_stag = 0                           !< Reset stagnation counter
                            fk_norm_prev = huge(1d0)             !< Reset previous residual norm to \f$+\infty\f$
                            lambda_LM = 1d-3                     !< Reset LM damping to initial value
                            best_fk_norm = huge(1d0)             !< Reset best residual norm to \f$+\infty\f$
                            skip_eval = .false.                  !< Force re-evaluation after predictor restart
                            cycle newton_loop                    !< Restart Newton loop with new initial guess
                        end if                                   !< End stagnation predictor fallback
                        !> Accept if best residual is within conditioning noise:
                        !> \f$\|\mathbf{f}_{\mathrm{best}}\|_\infty \le \sqrt{\varepsilon}\cdot\max(\|\mathbf{U}\mathbf{c}_{v}\|,1)\f$
                        if (best_fk_norm<=sqrt_eps_d*max(u_norm,1d0)) then  !< Within conditioning noise?
                            CV_flag=.true.                       !< Accept as converged
                            conc_nc=conc_nc_best                 !< Restore best-seen concentrations
                            call this%set_conc_prim_species(conc_nc(1:n_p))  !< Push best primary concentrations into chemistry object
                            call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p),log_gamma_var_act,conc_nc(n_p+1:))  !< Recompute secondary species at best point
                            call this%set_act_aq_species()        !< Update aqueous activities at best point
                            call this%compute_Rk(theta,Delta_t)  !< Compute kinetic reaction extents at accepted solution
                        end if                                   !< End conditioning-noise acceptance at stagnation
                        !print *, 'Newton stagnated after', n_stag, 'iterations without improvement. ||fk|| =', best_fk_norm
                        exit newton_loop                         !< Leave Newton loop (converged or stagnated)
                    end if                                       !< End stagnation-limit check
                end if                                           !< End stagnation branch
            end if                                               !< End non-converged branch (else of convergence check)
        end do newton_loop                                       !< End main Newton iteration loop
!> -----------------------------------------------------------------------
!> \section newton_ei_postproc Post-processing: cleanup and deallocation
!> -----------------------------------------------------------------------
    !print *, "DEBUG: Exiting Newton_EI_eq_kin_anal_ideal_opt2"
    !print *, "DEBUG: niter   = ", niter
    ! print *, "DEBUG: CV_flag = ", CV_flag
    ! print *, "DEBUG: Final conc_nc = ", conc_nc
    ! print *, "=========================================="
    !> Adaptive predictor damping update: adjust predictor_damp based on outcome
    !> If the predictor was invoked and Newton converged, increase damping (toward full step).
    !> If the predictor was invoked but Newton did not converge, decrease damping (more conservative).
    if (predictor_applied) then
        if (CV_flag) then
            predictor_damp = min(predictor_damp * predictor_grow, predictor_damp_max)  !< Success: grow toward full step
        else
            predictor_damp = max(predictor_damp * predictor_shrink, predictor_damp_min)  !< Failure: shrink for safety
        end if
    end if
    dbg_calls = dbg_calls + 1  !< Increment persistent call counter (limits diagnostic output to first call)
    deallocate(fk,dfk_dc1,drk_dc,rk_old,rk_new,rk_avg,conc_nc_best) !< Deallocate persistent work arrays
    deallocate(dfk_dp,dfk_dp_LM,Delta_p,log_c1,rk_trial,fk_trial,ipiv_dgesv,Uc,col_scale) !< Deallocate hoisted log-space work arrays
    deallocate(comp_mat,U_SkT_old,U_SkT_new,source_old,U_SkT_new_rk_old) !< Deallocate cached matrices
end subroutine Newton_EI_eq_kin_anal_ideal_opt2