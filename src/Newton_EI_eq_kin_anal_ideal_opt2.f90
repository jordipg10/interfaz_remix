!> \file Newton_EI_eq_kin_anal_ideal_opt2.f90
!> \brief Globalised Newton method (Armijo line search) for WMA reactive mixing with equilibrium AND
!>        Euler-implicit kinetic reactions (ideal solution, analytical Jacobians, option 2 for the
!>        mean kinetic rate).
!>
!> \details
!> Same algorithm as `Newton_EI_eq_ideal` (equilibrium-only), applied to the kinetic residual:
!> - **Ideal solution** assumption: the \f$\log_{10}\f$ activity coefficients are computed once and
!>   held fixed (all zero in the ideal case)
!> - **Log-space Newton steps** to eliminate the \f$1/c_1\f$ singularity in the Jacobian
!> - **Armijo backtracking line search** on \f$\phi=\tfrac12\|\mathbf{f}_k\|_2^2\f$ along the Newton
!>   direction, with steepest descent \f$-\mathbf{J}^T\mathbf{f}_k\f$ as fallback, for global convergence
!> - **Stagnation detection** and **best-solution tracking** for robustness
!> - **Re-speciation fallback** via `compute_c_nc_from_u_Newton_ideal` when max iterations or
!>   stagnation is reached (tried once before giving up): the kinetic term is frozen at the best
!>   iterate (a Picard step on the kinetic source), which reduces exactly to the equilibrium-only
!>   fallback when there are no kinetic reactions
!> - **Machine-epsilon safeguard** so that the loop terminates when the residual reaches
!>   double-precision round-off
!>
!> **Mathematical Formulation:**
!>
!> The nonlinear residual to drive to zero is (the previous-step \f$(1-\theta)\f$ kinetic
!> contribution is already folded into \f$\hat{\mathbf{u}}\f$ upstream, so only the new-step
!> \f$\theta\f$ term appears here):
!> \f[
!>   \mathbf{f}_k = \mathbf{U}\,\mathbf{c}_{v}
!>                  - \hat{\mathbf{u}}
!>                  - \Delta t\,\theta\,\lambda_{r}\,
!>                        (\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}\,\mathbf{r}_k(\mathbf{c}_v)
!>                  = \mathbf{0}
!> \f]
!>
!> where:
!> - \f$ \mathbf{U} \f$ = Component matrix mapping species to components [\f$n_p \times n_v\f$]
!> - \f$ \mathbf{c}_{v} \f$ = Variable activity species concentrations [C] (\f$n_v = n_p + n_{eq}\f$)
!> - \f$ \hat{\mathbf{u}} \f$ = Component concentrations after mixing [C] (size \f$n_p\f$)
!> - \f$ \lambda_{r} \f$ = Mixing ratio of new-step kinetic reaction amounts [-]
!> - \f$ (\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}} \f$ = Precomputed product of the component matrix and
!>   the transposed kinetic stoichiometry of the current reactive zone [\f$ n_p \times n_{kin} \f$]
!> - \f$ \mathbf{r}_k \f$ = Kinetic reaction rates at the current Newton iterate [1/T]
!> - \f$ \mathbf{r}_k^{\mathrm{old}} \f$ = Kinetic reaction rates from the previous time step (used only
!>   for the mean-rate average \f$\bar{\mathbf{r}}_k\f$ stored at convergence, not in the residual) [1/T]
!> - \f$ \theta \f$ = Time weighting factor (0 = explicit, 1 = fully implicit) [-]
!> - \f$ \Delta t \f$ = Time step size [T]
!>
!> The secondary variable activity species \f$\mathbf{c}_{2,v}\f$ are determined from the primary
!> species \f$\mathbf{c}_1\f$ via the mass action law (ideal solution), so the system has \f$n_p\f$
!> unknowns.
!>
!> **Log-Space Jacobian:**
!>
!> Defining \f$ p_j = \log_{10}(c_{1,j}) \f$, the Newton system is solved for \f$ \Delta\mathbf{p} \f$. Its
!> Jacobian is assembled DIRECTLY in log space by `compute_dfk_dp_EI_ideal`,
!> \f[
!>   \frac{\partial \mathbf{f}_k}{\partial p_j}
!>     = \ln 10\,\Bigl[\mathbf{U}_1(:,j)\,c_{1,j} + \mathbf{U}_2\,(\mathbf{c}_{2,v}\circ\mathbf{S}^*_{e,1}(:,j))\Bigr]
!>       - \Delta t\,\theta\,\lambda_r\,(\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}
!>         \Bigl[\frac{\partial \mathbf{r}_k}{\partial \log_{10} c_{1,j}}
!>         + \frac{\partial \mathbf{r}_k}{\partial \log_{10} \mathbf{c}_{2,v}}\,\mathbf{S}^*_{e,1}(:,j)\Bigr]
!> \f]
!> from the log-form kinetic Jacobian \f$\partial\mathbf{r}_k/\partial\log_{10}\mathbf{c}\f$ of
!> `compute_rk_Jac_rk_anal` (\f$\partial\log_{10}\mathbf{c}_{2,v}/\partial\log_{10}\mathbf{c}_1 =
!> \mathbf{S}^*_{e,1}\f$ for an ideal solution). No \f$1/c\f$ factor is ever formed (they would span as many
!> orders of magnitude as the concentrations), and the update
!> \f$ c_{1,j}^{\mathrm{new}} = 10^{p_j + \alpha\,\Delta p_j} \f$ keeps every concentration positive.
!>
!> **Globalisation — Armijo backtracking line search:**
!>
!> The log-space Newton system \f$\mathbf{J}\,\Delta\mathbf{p} = -\mathbf{f}_k\f$ is equilibrated by rows
!> (infinity norm of each row) and then by columns (infinity norm of each row-scaled column) and solved
!> with LAPACK `dgesv`; the Newton direction is unchanged in exact arithmetic, only the conditioning seen
!> by the LU factorisation improves (species with totals of 1e-1 and 1e-9 pivot alike). The step length
!> \f$\alpha\f$ along \f$\Delta\mathbf{p}\f$ is chosen by backtracking (\f$\alpha = 1, 1/2, 1/4, \ldots\f$) on
!> \f$\phi(\alpha)=\tfrac12\|\mathbf{f}_k(\mathbf{p}+\alpha\Delta\mathbf{p})\|_2^2\f$ until the Armijo
!> condition \f$\phi(\alpha)\le\phi(0)+c\,\alpha\,\phi'(0)\f$, \f$c=10^{-4}\f$, holds, with
!> \f$\phi'(0)=\mathbf{f}_k^T\mathbf{J}\Delta\mathbf{p}\f$. Before the line search \f$\Delta\mathbf{p}\f$ is
!> scaled uniformly so that \f$\|\Delta\mathbf{p}\|_\infty \le\f$ `max_dp` = 2 (two orders of magnitude,
!> i.e. a factor 100 in concentration per iteration, PHREEQC's default step size: a looser cap lets the
!> Jacobian be evaluated at absurd iterates where the kinetic term explodes) and \f$\phi'(0)\f$ is
!> re-evaluated for the scaled direction;
!> \f$\alpha\f$ is halved at most `max_LS_tries` \f$-1 = 9\f$ times. If the solve fails,
!> \f$\Delta\mathbf{p}\f$ is not a descent direction (\f$\phi'(0)\ge 0\f$) or its line search fails,
!> the steepest-descent direction \f$-\mathbf{J}^T\mathbf{f}_k\f$ is tried once instead. If no step is
!> accepted the iterate is kept and the iteration is treated as stagnated, since repeating it would
!> give the same result.
!>
!> **Convergence Criteria (user tolerances, OR round-off floor):**
!> \f[
!>   \bigl(\|\mathbf{f}_k\|_\infty < \varepsilon_{\mathrm{abs}}
!>   \;\wedge\; \|\mathbf{f}_k\|_\infty < \varepsilon_{\mathrm{rel}}\,u_s\bigr)
!>   \;\wedge\;
!>   \bigl(\log_{10}\|\mathbf{f}_k\|_\infty < \varepsilon^{\log}_{\mathrm{abs}}
!>   \;\wedge\; \log_{10}\|\mathbf{f}_k\|_\infty - \log_{10} u_s < \varepsilon^{\log}_{\mathrm{rel}}\bigr),
!>   \qquad u_s=\max(\|\mathbf{U}\mathbf{c}\|_\infty,\varepsilon_{\mathrm{mach}})
!> \f]
!> or, failing those, the round-off floor: every component residual within the forward-error bound of
!> its own evaluation,
!> \f[
!>   |f_{k,i}| \le \rho_i = (n_v+1+n_{kin})\,\varepsilon\,\Bigl(|\hat u_i| + \sum_j |U_{ij}|\,c_j
!>     + \Delta t\,\theta\,\lambda_r \sum_k |(\mathbf{U}\mathbf{S}_k^T)_{ik}|\,|r_k|\Bigr)
!>     + \sum_{j>n_p} |U_{ij}|\,c_j\,\delta_j ,
!> \f]
!> with \f$\delta_j\f$ the mass-action round-off of the secondary species as in
!> `compute_c_nc_from_u_Newton_ideal` (the round-off of the rate evaluation itself is not included).
!> The floor caps the achievable strictness of the user tolerances (an absolute tolerance below one ulp
!> of a large component total cannot be met) and keeps the fallback's convergence consistent with this
!> loop's. A non-converged best iterate is NEVER accepted (no \f$\sqrt{\varepsilon}\f$ acceptance): the
!> caller retries with other initial guesses or stops.
!>
!> **Re-speciation Fallback:**
!>
!> When the Newton loop exhausts `niter_max` or stagnates (`n_stag_max` consecutive iterations without
!> sufficient improvement, or one iteration in which no step could be accepted), a single
!> re-speciation attempt is made: `compute_c_nc_from_u_Newton_ideal` with the best-seen primary
!> concentrations as initial guess and the target components
!> \f$\hat{\mathbf{u}} + \Delta t\,\theta\,\lambda_r (\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}
!> \mathbf{r}_k(\mathbf{c}^{\mathrm{best}})\f$ (kinetic source frozen at the best iterate). All counters
!> are reset and the loop restarts only if the inner solve converged; otherwise the best iterate is
!> returned with `CV_flag = .false.`.
!>
!> **Post-processing at convergence:** the mean kinetic rate \f$\bar{\mathbf{r}}_k = \theta\,\mathbf{r}_k
!> + (1-\theta)\,\mathbf{r}_k^{\mathrm{old}}\f$ is stored (`set_rxn_rate_kin_avg`), and the kinetic
!> reaction amounts and their species contributions are computed (`compute_rxn_amnt_kin`,
!> `compute_contrib_rxn_amnt_kin`).
!>
!> \param[in,out] this   Aqueous chemistry object at current time step (class `aqueous_chemistry_c`);
!>                        on exit its primaries, secondaries, activities and kinetic rates are those of
!>                        the last speciation performed (the converged solution when CV_flag is TRUE)
!> \param[in]     u_hat  Component concentrations after mixing \f$\hat{\mathbf{u}}\f$ [C] (size \f$n_p\f$);
!>                        already includes the previous-step (1-theta) kinetic contribution
!> \param[in]     mix_ratio_r  Mixing ratio of kinetic reaction amounts \f$\lambda_{r}\f$ [-]
!> \param[in]     Delta_t  Time step size \f$\Delta t\f$ [T]
!> \param[in]     theta    Reaction time weighting factor \f$\theta \in [0,1]\f$ [-]
!> \param[in,out] conc_nc  Variable activity species concentrations \f$\mathbf{c}_{v}\f$ [C] (size
!>                          \f$n_v\f$; initial guess on entry, converged solution on exit if CV_flag is
!>                          TRUE, otherwise the iterate with the lowest \f$\|\mathbf{f}_k\|_\infty\f$ seen
!>                          since the last (re)start)
!> \param[out]    niter    Number of Newton iterations performed since the last (re)start [-]
!> \param[out]    CV_flag  Convergence flag: `.FALSE.` = did not converge, `.TRUE.` = converged
!>
!> \pre  The chemistry object `this` must have valid `solid_chemistry`, reactive-zone stoichiometry
!>       (`U_SkT_prod`), kinetic reactions and speciation algorithm data.
!> \pre  `conc_nc` must be pre-allocated with size \f$n_v\f$ and contain a strictly positive initial guess.
!>
!> \warning Only valid for ideal solutions (unit activity coefficients).
!> \warning LAPACK `dgesv` is called directly for the column-equilibrated log-space Newton system.
!> \warning A diagnostic line is printed to stdout on the max-iteration failure exit.
!>
!> \sa Newton_EI_eq_ideal, compute_c_nc_from_u_Newton_ideal
!> \sa aqueous_chemistry_c::compute_rk_Jac_rk_anal, aqueous_chemistry_c::compute_dfk_dp_EI_ideal
!>
!> \author Jordi
!> \ingroup chemistry

subroutine Newton_EI_eq_kin_anal_ideal_opt2(this,u_hat,mix_ratio_r,Delta_t,theta,conc_nc,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c  !< Import aqueous chemistry class for reactive mixing methods
    implicit none                                       !< Enforce explicit variable declarations
!> Arguments
    class(aqueous_chemistry_c) :: this              !< Aqueous chemistry object at current time step [-]
    real(kind=8), intent(in) :: u_hat(:)            !< Component concentrations after mixing (n_p) [C]; already includes the previous-step (1-theta) kinetic contribution
    real(kind=8), intent(in) :: mix_ratio_r         !< Mixing ratio of new kinetic reaction amounts \f$\lambda_{r}\f$ [-]
    real(kind=8), intent(in) :: Delta_t             !< Time step size \f$\Delta t\f$ [T]
    real(kind=8), intent(in) :: theta               !< Reaction time weighting factor \f$\theta \in [0,1]\f$ [-]
    real(kind=8), intent(inout) :: conc_nc(:)       !< Variable activity species concentrations \f$\mathbf{c}_{v}\f$ (n_v); initial guess on entry, solution on exit [C]
    integer(kind=4), intent(out) :: niter           !< Number of Newton iterations performed [-]
    logical, intent(out) :: CV_flag                 !< Convergence flag: .FALSE. = not converged, .TRUE. = converged
!> Local variables — kinetic rates and Jacobians
    real(kind=8), allocatable :: rk_new(:)          !< Kinetic reaction rates \f$\mathbf{r}_k\f$ at the current Newton iterate (n_kin) [1/T]
    real(kind=8), allocatable :: rk_trial(:)        !< Trial kinetic rates for the line search (n_kin) [1/T]
    real(kind=8), allocatable :: rk_best(:)         !< Kinetic rates at the best iterate (frozen kinetic source of the re-speciation fallback) (n_kin) [1/T]
    real(kind=8), allocatable :: rk_old(:)          !< Kinetic reaction rates from the previous time step \f$\mathbf{r}_k^{\mathrm{old}}\f$ (n_kin) [1/T]
    real(kind=8), allocatable :: rk_avg(:)          !< Mean rates \f$\bar{\mathbf{r}}_k = \theta\,\mathbf{r}_k + (1-\theta)\,\mathbf{r}_k^{\mathrm{old}}\f$ stored at convergence (n_kin) [1/T]
    real(kind=8), allocatable :: drk_dlogc(:,:)     !< Log-form Jacobian \f$\partial\mathbf{r}_k/\partial\log_{10}\mathbf{c}\f$ (n_kin x n_species, chem_syst order) [1/T]
!> Local variables — residual and concentrations
    real(kind=8), allocatable :: fk(:)              !< Newton residual \f$\mathbf{f}_k\f$ (n_p) [C]
    real(kind=8), allocatable :: fk_trial(:)        !< Trial residual for the line search (n_p) [C]
    real(kind=8), allocatable :: conc_nc_save(:)    !< Saved concentrations before the line-search trials (n_v) [C]
    real(kind=8), allocatable :: conc_nc_best(:)    !< Best concentrations seen since the last (re)start (lowest residual) (n_v) [C]
    real(kind=8), allocatable :: log_gamma_var_act(:) !< Log10 activity coefficients of the variable activity species (all zero: ideal) [-]
    real(kind=8), allocatable :: u_fix(:)           !< Target components of the re-speciation fallback: u_hat + kinetic source frozen at the best iterate (n_p) [C]
!> Local variables — log-space Newton system
    real(kind=8), allocatable :: log_c1(:)          !< \f$p_j = \log_{10}(c_{1,j})\f$ (n_p) [-]
    real(kind=8), allocatable :: Delta_p(:)         !< Search direction in log-space: Newton step or steepest descent \f$-\mathbf{J}^T\mathbf{f}_k\f$ (n_p) [-]
    real(kind=8), allocatable :: dfk_dp(:,:)        !< Log-space Jacobian \f$\partial\mathbf{f}_k/\partial\mathbf{p}\f$ (n_p x n_p), assembled directly by compute_dfk_dp_EI_ideal [-]
    real(kind=8), allocatable :: dfk_dp_eq(:,:)     !< Row- and column-equilibrated copy of the log-space Jacobian (overwritten by dgesv) (n_p x n_p) [-]
    real(kind=8), allocatable :: row_scale(:)       !< Row scaling factors of the Jacobian equilibration (infinity norm of each row) (n_p) [-]
    real(kind=8), allocatable :: col_scale(:)       !< Column scaling factors of the Jacobian equilibration (infinity norm of each row-scaled column) (n_p) [-]
    integer(kind=4), allocatable :: ipiv_dgesv(:)   !< Pivot indices of LAPACK dgesv (n_p)
    integer(kind=4) :: info_dgesv                   !< Return code of LAPACK dgesv (0 = success)
    real(kind=8), parameter :: ln10 = log(10d0)     !< \f$\ln(10)\f$ [-]
    real(kind=8), parameter :: max_dp = 2d0         !< Maximum allowed \f$|\Delta p_j|\f$ per iteration: 2 orders of magnitude (a factor 100 in concentration per Newton step, PHREEQC's default step size); keeps the Jacobian from being evaluated at absurd iterates [-]
    real(kind=8) :: scale_factor                    !< \f$\|\Delta\mathbf{p}\|_\infty\f$ before bounding [-]
!> Local variables — search direction and line search
    integer(kind=4) :: i_dir                        !< Direction attempt counter (Newton, then steepest descent) [-]
    logical :: use_sd                               !< Current direction is steepest descent
    logical :: sing_flag                            !< dgesv failed or returned NaN/huge values
    real(kind=8) :: alpha                           !< Backtracking step length \f$\alpha \in (0,1]\f$ [-]
    integer(kind=4) :: i_LS                         !< Line search counter [-]
    integer(kind=4), parameter :: max_LS_tries = 10 !< Maximum backtracking steps per direction [-]
    real(kind=8), parameter :: c_armijo = 1d-4      !< Armijo sufficient-decrease constant [-]
    logical :: step_accepted                        !< A step was accepted this iteration
    real(kind=8) :: fk_norm                         !< \f$\|\mathbf{f}_k\|_\infty\f$ [C]
    real(kind=8) :: phi0                            !< \f$\phi(0)=\tfrac12\|\mathbf{f}_k\|_2^2\f$ [C²]
    real(kind=8) :: dphi0                           !< \f$\phi'(0)=\mathbf{f}_k^T\mathbf{J}\Delta\mathbf{p}\f$ (must be < 0) [C²]
    real(kind=8) :: phi_trial                       !< \f$\phi(\alpha)\f$ at the trial point [C²]
!> Local variables — stagnation detection and re-speciation fallback
    integer(kind=4) :: n_stag                       !< Consecutive iterations without sufficient improvement (set to n_stag_max when no step could be accepted) [-]
    integer(kind=4), parameter :: n_stag_max = 20   !< Stagnation limit before the fallback / exit [-]
    real(kind=8) :: fk_norm_prev                    !< Reference residual norm of the stagnation test [C]
    real(kind=8), parameter :: stag_rtol = 1d-3     !< Relative improvement threshold of the stagnation test [-]
    real(kind=8) :: best_fk_norm                    !< Lowest \f$\|\mathbf{f}_k\|_\infty\f$ seen since the last (re)start [C]
    logical :: fallback_applied                     !< The re-speciation fallback has already been tried
    integer(kind=4) :: niter_spec                   !< Inner speciation Newton iterations (fallback) [-]
    logical :: CV_flag_spec                         !< Inner speciation convergence flag (fallback)
!> Local variables — dimensions and cached parameters
    integer(kind=4) :: n_p                          !< Number of primary species [-]
    integer(kind=4) :: n_v                          !< Number of variable activity species (= n_p + n_eq) [-]
    integer(kind=4) :: n_kin                        !< Number of kinetic reactions of the current zone [-]
    integer(kind=4) :: ii                           !< Loop index [-]
    integer(kind=4) :: niter_max                    !< Maximum Newton iterations (CV_params) [-]
    real(kind=8) :: abs_tol                         !< Absolute tolerance (CV_params) [C]
    real(kind=8) :: rel_tol                         !< Relative tolerance (CV_params) [-]
    real(kind=8) :: log_abs_tol                     !< Log absolute tolerance (CV_params) [-]
    real(kind=8) :: log_rel_tol                     !< Log relative tolerance (CV_params) [-]
    real(kind=8) :: eps_d                           !< Machine epsilon [-]
    real(kind=8) :: u_norm                          !< \f$\|\mathbf{U}\mathbf{c}_{v}\|_\infty\f$ [C]
    real(kind=8) :: u_scale                         !< Relative-tolerance scale \f$\max(\|\mathbf{U}\mathbf{c}\|_\infty,\varepsilon_{\mathrm{mach}})\f$ [C]
    real(kind=8) :: dt_theta_mr_new                 !< \f$\Delta t\,\theta\,\lambda_{r}\f$ [T]
    real(kind=8), allocatable :: comp_mat(:,:)      !< Component matrix \f$\mathbf{U}\f$ (n_p x n_v) [-]
    real(kind=8), allocatable :: U_SkT_new(:,:)     !< \f$(\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}\f$ (n_p x n_kin) [-]
    real(kind=8), allocatable :: Uc(:)              !< \f$\mathbf{U}\mathbf{c}_{v}\f$ (n_p) [C]
    logical :: skip_eval                            !< Reuse the speciation, rates and residual of the accepted line-search point
!> Local variables — round-off floor of the residual (termination safeguard; bound of compute_c_nc_from_u_Newton_ideal plus the kinetic term)
    integer(kind=4) :: n_eq                         !< Number of equilibrium reactions (= n_v - n_p) [-]
    real(kind=8) :: zero_tol                        !< Numerically-zero concentration (CV_params%zero), for components with exactly zero \f$\hat u_i\f$ [C]
    real(kind=8), allocatable :: res_floor(:)       !< Round-off bound \f$\rho_i\f$ of each component residual at the current iterate (n_p) [C]
    real(kind=8), allocatable :: abs_comp_mat(:,:)  !< \f$|\mathbf{U}|\f$ (n_p x n_v) [-]
    real(kind=8), allocatable :: abs_U2(:,:)        !< \f$|\mathbf{U}_2|\f$ (n_p x n_eq) [-]
    real(kind=8), allocatable :: abs_U_SkT(:,:)     !< \f$|\mathbf{U}\mathbf{S}_k^T|\f$ (n_p x n_kin) [-]
    real(kind=8), allocatable :: abs_Se_star(:,:)   !< \f$|\mathbf{S}^*_{e,v,1}|\f$ (n_eq x n_p) [-]
    real(kind=8), allocatable :: abs_logK_star(:)   !< \f$|\log_{10}\mathbf{K}^*|\f$ (n_eq) [-]
    real(kind=8), allocatable :: abs_log_gamma_1(:) !< \f$|\log_{10}\gamma|\f$ of the primary species (n_p) [-]
    real(kind=8), allocatable :: abs_log_gamma_2(:) !< \f$|\log_{10}\gamma|\f$ of the secondary variable activity species (n_eq) [-]
    real(kind=8), allocatable :: abs_u_hat(:)       !< \f$|\hat{\mathbf{u}}|\f$ (n_p) [C]
    real(kind=8), allocatable :: row_sum_abs_U(:)   !< \f$\sum_j |U_{ij}|\f$ (n_p) [-]
    real(kind=8), allocatable :: abs_log_c1(:)      !< \f$|\log_{10} c_{1,k}|\f$ at the current iterate (n_p) [-]
    real(kind=8), allocatable :: rel_err_c2(:)      !< Relative round-off bound \f$\delta_j\f$ of each secondary species at the current iterate (n_eq) [-]
!> -----------------------------------------------------------------------
!> \section newton_ei_preproc Pre-processing: initialisation and allocation
!> -----------------------------------------------------------------------
    niter=0
    CV_flag=.false.
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_v=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    niter_max=this%solid_chemistry%reactive_zone%CV_params%niter_max
    abs_tol=this%solid_chemistry%reactive_zone%CV_params%abs_tol
    rel_tol=this%solid_chemistry%reactive_zone%CV_params%rel_tol
    log_abs_tol=this%solid_chemistry%reactive_zone%CV_params%log_abs_tol
    log_rel_tol=this%solid_chemistry%reactive_zone%CV_params%log_rel_tol
    eps_d=epsilon(1d0)
    n_kin=this%solid_chemistry%mineral_zone%num_minerals_kin+&
        this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
    comp_mat=this%solid_chemistry%reactive_zone%speciation_alg%comp_mat
    U_SkT_new=this%solid_chemistry%reactive_zone%U_SkT_prod
    allocate(drk_dlogc(n_kin,this%solid_chemistry%reactive_zone%speciation_alg%num_species))
    allocate(rk_new(n_kin),rk_trial(n_kin),rk_best(n_kin),rk_avg(n_kin))
    allocate(fk(n_p),fk_trial(n_p),conc_nc_best(n_v),conc_nc_save(n_v),u_fix(n_p))
    allocate(dfk_dp(n_p,n_p),dfk_dp_eq(n_p,n_p),Delta_p(n_p),log_c1(n_p),Uc(n_p),row_scale(n_p),col_scale(n_p),ipiv_dgesv(n_p))
    drk_dlogc=0d0
    dfk_dp=0d0
    rk_best=0d0
    n_stag=0
    fk_norm_prev=huge(1d0)
    best_fk_norm=huge(1d0)
    fallback_applied=.false.
    skip_eval=.false.
    log_gamma_var_act=this%compute_log_act_coeffs_var_act_ideal() !> fixed during the iteration (all zeros: ideal)
    rk_old=this%get_rk_old_new_set() !> previous-step rates in the current kinetic set (rows of rk_new): only for the mean rate stored at convergence
    dt_theta_mr_new=Delta_t*theta*mix_ratio_r
    !> Loop-invariant factors of the round-off bound of the residual (see the convergence check)
    n_eq=n_v-n_p
    zero_tol=this%solid_chemistry%reactive_zone%CV_params%zero
    abs_comp_mat=abs(comp_mat)
    abs_U2=abs(comp_mat(:,n_p+1:n_v))
    abs_U_SkT=abs(U_SkT_new)
    abs_u_hat=abs(u_hat)
    row_sum_abs_U=sum(abs_comp_mat,dim=2)
    abs_log_gamma_1=abs(log_gamma_var_act(1:n_p))
    abs_log_gamma_2=abs(log_gamma_var_act(n_p+1:n_v))
    if (n_eq>0) then !> mass-action factors of the secondary species (as used by compute_c2v_from_c1_ideal)
        abs_Se_star=abs(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star)
        abs_logK_star=abs(this%logK_star)
    else
        allocate(abs_Se_star(0,n_p),abs_logK_star(0))
    end if
    allocate(res_floor(n_p),abs_log_c1(n_p),rel_err_c2(n_eq))
!> -----------------------------------------------------------------------
!> \section newton_ei_process Process: Newton loop
!> -----------------------------------------------------------------------
    newton_loop: do
        niter=niter+1
    !> \subsection newton_ei_maxiter Maximum iteration check and re-speciation fallback
        if (niter>niter_max) then
            if (.not. fallback_applied) then
                fallback_applied=.true.
                call respeciate_from_best()
                if (CV_flag_spec) cycle newton_loop !> restart from the re-speciated state (counters reset)
            end if
            conc_nc=conc_nc_best !> return the best iterate seen
            niter=niter_max
            print *, 'Newton_EI_eq_kin_anal_ideal_opt2 did not converge after', niter_max, &
                'iterations. Best ||fk|| =', best_fk_norm
            exit newton_loop
        end if
    !> \subsection newton_ei_eval Speciation, kinetic rates and residual at the current iterate
    !> (reused from the accepted line-search point when available: it already speciated the object,
    !> evaluated the rates and the residual there)
        if (skip_eval) then
            skip_eval=.false.
            fk=fk_trial
            rk_new=rk_trial
            Uc=fk+u_hat+dt_theta_mr_new*matmul(U_SkT_new,rk_new) !> recover U*c of the accepted point
        else
            call this%set_conc_prim_species(conc_nc(1:n_p))
            call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p),log_gamma_var_act,conc_nc(n_p+1:))
            call this%set_act_aq_species()
            drk_dlogc=0d0
            call this%compute_rk_Jac_rk_anal(rk_new,drk_dlogc)
            Uc=matmul(comp_mat,conc_nc)
            fk=Uc-u_hat-dt_theta_mr_new*matmul(U_SkT_new,rk_new)
        end if
        fk_norm=maxval(abs(fk))
        u_norm=maxval(abs(Uc))
    !> \subsection newton_ei_best Best-solution tracking (initial guess of the fallback; returned on failure)
        if (fk_norm<best_fk_norm) then
            best_fk_norm=fk_norm
            conc_nc_best=conc_nc
            rk_best=rk_new
        end if
    !> \subsection newton_ei_conv Convergence check (user tolerances OR round-off floor)
    !> Round-off floor (evaluated only when the user tolerances are not met): every component residual
    !> within the forward-error bound of its own evaluation,
    !> rho_i = (n_v+1+n_kin) eps (|u_hat_i| + sum_j |U_ij| c_j + Dt theta lambda_r sum_k |(U Sk^T)_ik| |r_k|)
    !>         + sum_{j>n_p} |U_ij| c_j delta_j,
    !> i.e. the inner product of n_v+1+n_kin terms plus the propagated mass-action round-off delta_j of
    !> the secondary species (as in compute_c_nc_from_u_Newton_ideal; the round-off of the rate evaluation
    !> itself is not included). It caps the achievable strictness of the user tolerances (an absolute
    !> tolerance below one ulp of a large total cannot be met) and makes the fallback's notion of
    !> convergence consistent with this loop's.
        u_scale=max(u_norm,eps_d)
        CV_flag=(fk_norm<abs_tol .and. fk_norm<rel_tol*u_scale) .and. &
            (fk_norm>0d0 .and. log10(fk_norm)<log_abs_tol .and. (log10(fk_norm)-log10(u_scale))<log_rel_tol)
        if (.not. CV_flag) then
            abs_log_c1=abs(log10(max(conc_nc(1:n_p),tiny(1d0))))
            rel_err_c2=eps_d*(1d0+ln10*dble(n_p+3)*(matmul(abs_Se_star,abs_log_c1+abs_log_gamma_1)+abs_logK_star+abs_log_gamma_2))
            res_floor=dble(n_v+1+n_kin)*eps_d*(abs_u_hat+matmul(abs_comp_mat,abs(conc_nc))+ &
                abs(dt_theta_mr_new)*matmul(abs_U_SkT,abs(rk_new)))+matmul(abs_U2,abs(conc_nc(n_p+1:n_v))*rel_err_c2)
            do ii=1,n_p !> a component with exactly zero total has no root at finite concentrations: its species are accepted as absent below the numerically-zero concentration
                if (u_hat(ii)==0d0) res_floor(ii)=max(res_floor(ii),zero_tol*row_sum_abs_U(ii))
            end do
            CV_flag=all(abs(fk)<=res_floor)
        end if
        if (CV_flag) then
            !> the object holds the speciation and the rates of this iterate; store the mean rate and
            !> the kinetic reaction amounts and contributions of the step
            rk_avg=theta*rk_new+(1d0-theta)*rk_old
            call this%set_rxn_rate_kin_avg(rk_avg)
            call this%compute_rxn_amnt_kin(Delta_t)
            call this%compute_contrib_rxn_amnt_kin(mix_ratio_r)
            exit newton_loop
        end if
    !> \subsection newton_ei_jacobian Jacobian assembly, directly in log-space
    !> drk_dlogc (d r_k / d log10 c, compute_rk_Jac_rk_anal) is passed in full chem_syst order: the callee
    !> permutes its columns through reactive_zone%ind_var_act_species, whose values can exceed n_v in
    !> multi-zone setups, and assembles dfk/dp = ln10 [U_1 c_1 + U_2 dc_2v/dp] - Dt theta lambda_r (U S_k^T)
    !> [d r_k/d log10 c_1 + d r_k/d log10 c_2v Se*_1] without forming any 1/c factor
        call this%compute_dfk_dp_EI_ideal(conc_nc(1:n_p),conc_nc(n_p+1:),drk_dlogc,Delta_t,theta,mix_ratio_r,dfk_dp)
        do ii=1,n_p
            log_c1(ii)=log10(max(conc_nc(ii),tiny(1d0)))
        end do
    !> \subsection newton_ei_step Newton direction and Armijo backtracking line search in log-space
    !> (same scheme as Newton_EI_eq_ideal: Newton step on the equilibrated system, steepest
    !> descent -J^T f when the solve fails, the direction is not a descent direction or its line
    !> search fails; no step accepted -> stagnation)
    !> Equilibration of the linear system J Delta_p = -f: rows first (each row divided by its infinity
    !> norm, so that components with totals of 1e-1 and 1e-9 weigh alike in the pivoting), then columns
    !> (each column of the row-scaled matrix divided by its infinity norm). With R = diag(1/row_scale),
    !> C = diag(1/col_scale): (R J C) y = -R f, Delta_p = C y. The Newton direction is unchanged in exact
    !> arithmetic; only the conditioning seen by dgesv improves. The merit function phi = 0.5 |f|^2 and
    !> its slope are evaluated with the unscaled J and f.
        conc_nc_save=conc_nc
        step_accepted=.false.
        dfk_dp_eq=dfk_dp
        do ii=1,n_p
            row_scale(ii)=maxval(abs(dfk_dp_eq(ii,:)))
            if (row_scale(ii)<tiny(1d0)) row_scale(ii)=1d0
            dfk_dp_eq(ii,:)=dfk_dp_eq(ii,:)/row_scale(ii)
        end do
        do ii=1,n_p
            col_scale(ii)=maxval(abs(dfk_dp_eq(:,ii)))
            if (col_scale(ii)<tiny(1d0)) col_scale(ii)=1d0
            dfk_dp_eq(:,ii)=dfk_dp_eq(:,ii)/col_scale(ii)
        end do
        Delta_p=-fk/row_scale
        call dgesv(n_p,1,dfk_dp_eq,n_p,ipiv_dgesv,Delta_p,n_p,info_dgesv)
        sing_flag=(info_dgesv/=0)
        if (.not. sing_flag) then
            Delta_p=Delta_p/col_scale
            sing_flag=any(Delta_p/=Delta_p) .or. maxval(abs(Delta_p))>1d20
        end if
        use_sd=sing_flag
        if (use_sd) Delta_p=-matmul(transpose(dfk_dp),fk)
        dir_loop: do i_dir=1,2
            dphi0=dot_product(fk,matmul(dfk_dp,Delta_p))
            if (dphi0>=0d0) then !> not a descent direction
                if (use_sd) exit dir_loop !> J^T f = 0 to round-off: stationary point of phi, no step possible
                use_sd=.true.
                Delta_p=-matmul(transpose(dfk_dp),fk)
                cycle dir_loop
            end if
            scale_factor=maxval(abs(Delta_p))
            if (scale_factor>max_dp) Delta_p=Delta_p*(max_dp/scale_factor)
            phi0=5d-1*dot_product(fk,fk)
            dphi0=dot_product(fk,matmul(dfk_dp,Delta_p))
            alpha=1d0
            do i_LS=1,max_LS_tries
                do ii=1,n_p
                    conc_nc(ii)=10d0**(log_c1(ii)+alpha*Delta_p(ii))
                end do
                call this%set_conc_prim_species(conc_nc(1:n_p))
                call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p),log_gamma_var_act,conc_nc(n_p+1:))
                call this%set_act_aq_species()
                drk_dlogc=0d0
                call this%compute_rk_Jac_rk_anal(rk_trial,drk_dlogc)
                fk_trial=matmul(comp_mat,conc_nc)-u_hat-dt_theta_mr_new*matmul(U_SkT_new,rk_trial)
                phi_trial=5d-1*dot_product(fk_trial,fk_trial)
                if (phi_trial<=phi0+c_armijo*alpha*dphi0) then
                    step_accepted=.true.
                    skip_eval=.true.
                    exit
                end if
                alpha=alpha*0.5d0
            end do
            if (step_accepted) exit dir_loop
            conc_nc=conc_nc_save !> line search failed: restore the pre-step state
            if (use_sd) exit dir_loop !> steepest descent failed too: no step this iteration
            use_sd=.true.
            Delta_p=-matmul(transpose(dfk_dp),fk)
        end do dir_loop
        if (.not. step_accepted) then
            !> the object must describe the restored state (the line search speciated trial points)
            call this%set_conc_prim_species(conc_nc(1:n_p))
            call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p),log_gamma_var_act,conc_nc(n_p+1:))
            call this%set_act_aq_species()
            drk_dlogc=0d0
            call this%compute_rk_Jac_rk_anal(rk_new,drk_dlogc)
        end if
    !> \subsection newton_ei_stag Stagnation detection and re-speciation fallback
    !> A rejected step counts as full stagnation at once: the state is unchanged, so every following
    !> iteration would repeat this one exactly
        if (.not. step_accepted) then
            n_stag=n_stag_max
        else if (fk_norm<fk_norm_prev*(1d0-stag_rtol)) then
            n_stag=0
            fk_norm_prev=fk_norm
        else
            n_stag=n_stag+1
        end if
        if (n_stag>=n_stag_max) then
            if (.not. fallback_applied) then
                fallback_applied=.true.
                call respeciate_from_best()
                if (CV_flag_spec) cycle newton_loop
            end if
            conc_nc=conc_nc_best !> return the best iterate seen
            exit newton_loop !> stagnated: leave with CV_flag = .false.
        end if
    end do newton_loop
!> -----------------------------------------------------------------------
!> \section newton_ei_postproc Post-processing: deallocation
!> -----------------------------------------------------------------------
    deallocate(fk,fk_trial,drk_dlogc,rk_old,rk_new,rk_trial,rk_best,rk_avg,conc_nc_best,conc_nc_save,u_fix)
    deallocate(dfk_dp,dfk_dp_eq,Delta_p,log_c1,ipiv_dgesv,Uc,row_scale,col_scale,comp_mat,U_SkT_new,log_gamma_var_act)
    deallocate(res_floor,abs_comp_mat,abs_U2,abs_U_SkT,abs_Se_star,abs_logK_star,abs_log_gamma_1,abs_log_gamma_2, &
        abs_u_hat,row_sum_abs_U,abs_log_c1,rel_err_c2)

contains

    !> Re-speciation fallback: inner speciation Newton (compute_c_nc_from_u_Newton_ideal) from the
    !> best-seen primaries towards the components u_hat + Delta_t*theta*lambda_r*(U*Sk^T)*rk(c_best),
    !> i.e. with the kinetic source frozen at the best iterate (a Picard step on the kinetic term; the
    !> equilibrium-only fallback of Newton_EI_eq_ideal when n_kin = 0). On a converged inner solve the
    !> outer counters are reset so the Newton loop restarts from the new state; otherwise the state of
    !> the best iterate is restored and CV_flag_spec is .false.
    subroutine respeciate_from_best()
        u_fix=u_hat+dt_theta_mr_new*matmul(U_SkT_new,rk_best)
        call this%compute_c_nc_from_u_Newton_ideal(conc_nc_best(1:n_p),u_fix,conc_nc,niter_spec,CV_flag_spec)
        if (CV_flag_spec) then
            niter=0
            n_stag=0
            fk_norm_prev=huge(1d0)
            best_fk_norm=huge(1d0)
            skip_eval=.false.
        else
            conc_nc=conc_nc_best
        end if
    end subroutine respeciate_from_best
end subroutine Newton_EI_eq_kin_anal_ideal_opt2
