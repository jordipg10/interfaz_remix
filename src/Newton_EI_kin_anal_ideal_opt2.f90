!> \file Newton_EI_kin_anal_ideal_opt2.f90
!> \brief Globalised Newton method (Armijo line search) for WMA reactive mixing with ONLY Euler-implicit
!>        kinetic reactions (no equilibrium reactions; ideal solution, analytical Jacobians, option 2 for
!>        the mean kinetic rate).
!>
!> \details
!> Same algorithm as `Newton_EI_eq_kin_anal_ideal_opt2`, specialised to a chemical system WITHOUT
!> equilibrium reactions. The consequences of \f$n_{eq}=0\f$ are:
!> - All variable activity species are primary: \f$n_v = n_p\f$, no secondary species
!> - The component matrix is the identity, \f$\mathbf{U}=\mathbf{I}\f$: the components are the species
!>   themselves (\f$\hat{\mathbf{u}}=\hat{\mathbf{c}}\f$, \f$\mathbf{U}\mathbf{S}_k^T=\mathbf{S}_k^T\f$)
!> - No mass-action law (no `compute_c2v_from_c1_ideal`) and no speciation solve anywhere
!> - The log-space Jacobian is assembled directly here (no `compute_dfk_dp_EI_ideal`)
!> - The re-speciation fallback of the equilibrium version (`compute_c_nc_from_u_Newton_ideal` towards
!>   \f$\hat{\mathbf{u}}\f$ plus the kinetic source frozen at the best iterate) reduces to the explicit
!>   Picard step \f$\mathbf{c}=\hat{\mathbf{c}}+\Delta t\,\theta\,\lambda_r\,\mathbf{S}_k^T\,
!>   \mathbf{r}_k(\mathbf{c}^{\mathrm{best}})\f$, which "converges" iff it is strictly positive
!>
!> Everything else is identical:
!> - **Ideal solution** (unit activity coefficients)
!> - **Log-space Newton steps** (positivity preserved by construction)
!> - **Armijo backtracking line search** on \f$\phi=\tfrac12\|\mathbf{f}_k\|_2^2\f$ along the Newton
!>   direction of the row/column-equilibrated system, with steepest descent \f$-\mathbf{J}^T\mathbf{f}_k\f$
!>   as fallback
!> - **Stagnation detection** and **best-solution tracking**
!> - **Round-off floor** of the residual as termination safeguard; a non-converged best iterate is NEVER
!>   accepted (no \f$\sqrt{\varepsilon}\f$ acceptance)
!>
!> **Mathematical Formulation:**
!>
!> The nonlinear residual to drive to zero is (the previous-step \f$(1-\theta)\f$ kinetic contribution is
!> already folded into \f$\hat{\mathbf{c}}\f$ upstream, so only the new-step \f$\theta\f$ term appears):
!> \f[
!>   \mathbf{f}_k = \mathbf{c} - \hat{\mathbf{c}}
!>                  - \Delta t\,\theta\,\lambda_{r}\,\mathbf{S}_k^{T}\,\mathbf{r}_k(\mathbf{c}) = \mathbf{0}
!> \f]
!>
!> **Log-Space Jacobian:** with \f$p_j=\log_{10}c_j\f$ and the log-form kinetic Jacobian
!> \f$\partial\mathbf{r}_k/\partial\log_{10}\mathbf{c}\f$ of `compute_rk_Jac_rk_anal` (columns permuted
!> from chem_syst order to the zone-local order through `ind_var_act_species`):
!> \f[
!>   \frac{\partial \mathbf{f}_k}{\partial \mathbf{p}}
!>     = \ln 10\,\mathrm{diag}(\mathbf{c})
!>       - \Delta t\,\theta\,\lambda_r\,\mathbf{S}_k^T\,\frac{\partial \mathbf{r}_k}{\partial \log_{10}\mathbf{c}}
!> \f]
!> No \f$1/c\f$ factor is ever formed; the update \f$c_j^{\mathrm{new}}=10^{p_j+\alpha\Delta p_j}\f$
!> keeps every concentration positive.
!>
!> **Convergence Criteria (user tolerances, OR round-off floor):**
!> \f[
!>   \bigl(\|\mathbf{f}_k\|_\infty < \varepsilon_{\mathrm{abs}}
!>   \;\wedge\; \|\mathbf{f}_k\|_\infty < \varepsilon_{\mathrm{rel}}\,c_s\bigr)
!>   \;\wedge\;
!>   \bigl(\log_{10}\|\mathbf{f}_k\|_\infty < \varepsilon^{\log}_{\mathrm{abs}}
!>   \;\wedge\; \log_{10}\|\mathbf{f}_k\|_\infty - \log_{10} c_s < \varepsilon^{\log}_{\mathrm{rel}}\bigr),
!>   \qquad c_s=\max(\|\mathbf{c}\|_\infty,\varepsilon_{\mathrm{mach}})
!> \f]
!> or, failing those, every residual within the forward-error bound of its own evaluation
!> (\f$2+n_{kin}\f$ terms; the round-off of the rate evaluation itself is not included),
!> \f[
!>   |f_{k,i}| \le \rho_i = (2+n_{kin})\,\varepsilon\,\Bigl(|\hat c_i| + c_i
!>     + \Delta t\,\theta\,\lambda_r \sum_k |S^T_{k,ik}|\,|r_k|\Bigr).
!> \f]
!>
!> **Post-processing at convergence:** the mean kinetic rate \f$\bar{\mathbf{r}}_k = \theta\,\mathbf{r}_k
!> + (1-\theta)\,\mathbf{r}_k^{\mathrm{old}}\f$ is stored (`set_rxn_rate_kin_avg`), and the kinetic
!> reaction amounts and their species contributions are computed (`compute_rxn_amnt_kin`,
!> `compute_contrib_rxn_amnt_kin`).
!>
!> \param[in,out] this   Aqueous chemistry object at current time step; on exit its concentrations,
!>                        activities and kinetic rates are those of the last evaluation performed (the
!>                        converged solution when CV_flag is TRUE)
!> \param[in]     c_hat  Species concentrations after mixing \f$\hat{\mathbf{c}}\f$ [C] (size \f$n_p\f$);
!>                        already includes the previous-step (1-theta) kinetic contribution
!> \param[in]     mix_ratio_r  Mixing ratio of kinetic reaction amounts \f$\lambda_{r}\f$ [-]
!> \param[in]     Delta_t  Time step size \f$\Delta t\f$ [T]
!> \param[in]     theta    Reaction time weighting factor \f$\theta \in [0,1]\f$ [-]
!> \param[in,out] conc_nc  Variable activity species concentrations \f$\mathbf{c}\f$ [C] (size \f$n_p\f$;
!>                          initial guess on entry, converged solution on exit if CV_flag is TRUE, otherwise
!>                          the iterate with the lowest \f$\|\mathbf{f}_k\|_\infty\f$ seen since the last
!>                          (re)start)
!> \param[out]    niter    Number of Newton iterations performed since the last (re)start [-]
!> \param[out]    CV_flag  Convergence flag: `.FALSE.` = did not converge, `.TRUE.` = converged
!>
!> \pre  The reactive zone of `this` has NO equilibrium reactions (n_v = n_p); `U_SkT_prod` (= S_k^T) and
!>       the kinetic reactions are set.
!> \pre  `conc_nc` must be pre-allocated with size \f$n_p\f$ and contain a strictly positive initial guess.
!>
!> \warning Only valid for ideal solutions (unit activity coefficients).
!> \warning LAPACK `dgesv` is called directly for the equilibrated log-space Newton system.
!> \warning A diagnostic line is printed to stdout on the max-iteration failure exit.
!>
!> \sa Newton_EI_eq_kin_anal_ideal_opt2, Newton_EI_eq_ideal
!> \sa aqueous_chemistry_c::compute_rk_Jac_rk_anal

subroutine Newton_EI_kin_anal_ideal_opt2(this,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c  !< Import aqueous chemistry class
    implicit none                                       !< Enforce explicit variable declarations
!> Arguments
    class(aqueous_chemistry_c) :: this              !< Aqueous chemistry object at current time step [-]
    real(kind=8), intent(in) :: c_hat(:)            !< Species concentrations after mixing (n_p) [C]; already includes the previous-step (1-theta) kinetic contribution
    real(kind=8), intent(in) :: mix_ratio_r         !< Mixing ratio of new kinetic reaction amounts \f$\lambda_{r}\f$ [-]
    real(kind=8), intent(in) :: Delta_t             !< Time step size \f$\Delta t\f$ [T]
    real(kind=8), intent(in) :: theta               !< Reaction time weighting factor \f$\theta \in [0,1]\f$ [-]
    real(kind=8), intent(inout) :: conc_nc(:)       !< Variable activity species concentrations (n_p); initial guess on entry, solution on exit [C]
    integer(kind=4), intent(out) :: niter           !< Number of Newton iterations performed [-]
    logical, intent(out) :: CV_flag                 !< Convergence flag: .FALSE. = not converged, .TRUE. = converged
!> Local variables — kinetic rates and Jacobians
    real(kind=8), allocatable :: rk_new(:)          !< Kinetic reaction rates at the current Newton iterate (n_kin) [1/T]
    real(kind=8), allocatable :: rk_trial(:)        !< Trial kinetic rates for the line search (n_kin) [1/T]
    real(kind=8), allocatable :: rk_best(:)         !< Kinetic rates at the best iterate (frozen kinetic source of the Picard fallback) (n_kin) [1/T]
    real(kind=8), allocatable :: rk_old(:)          !< Kinetic reaction rates from the previous time step (n_kin) [1/T]
    real(kind=8), allocatable :: rk_avg(:)          !< Mean rates \f$\bar{\mathbf{r}}_k = \theta\,\mathbf{r}_k + (1-\theta)\,\mathbf{r}_k^{\mathrm{old}}\f$ stored at convergence (n_kin) [1/T]
    real(kind=8), allocatable :: drk_dlogc(:,:)     !< Log-form Jacobian \f$\partial\mathbf{r}_k/\partial\log_{10}\mathbf{c}\f$ (n_kin x n_species, chem_syst order) [1/T]
    real(kind=8), allocatable :: drk_dlogc_zone(:,:) !< Columns of drk_dlogc permuted to the zone-local species order (n_kin x n_p) [1/T]
!> Local variables — residual and concentrations
    real(kind=8), allocatable :: fk(:)              !< Newton residual \f$\mathbf{f}_k\f$ (n_p) [C]
    real(kind=8), allocatable :: fk_trial(:)        !< Trial residual for the line search (n_p) [C]
    real(kind=8), allocatable :: conc_nc_save(:)    !< Saved concentrations before the line-search trials (n_p) [C]
    real(kind=8), allocatable :: conc_nc_best(:)    !< Best concentrations seen since the last (re)start (lowest residual) (n_p) [C]
    real(kind=8), allocatable :: c_fix(:)           !< Picard fallback: c_hat + kinetic source frozen at the best iterate (n_p) [C]
!> Local variables — log-space Newton system
    real(kind=8), allocatable :: log_c(:)           !< \f$p_j = \log_{10}(c_j)\f$ (n_p) [-]
    real(kind=8), allocatable :: Delta_p(:)         !< Search direction in log-space: Newton step or steepest descent \f$-\mathbf{J}^T\mathbf{f}_k\f$ (n_p) [-]
    real(kind=8), allocatable :: dfk_dp(:,:)        !< Log-space Jacobian \f$\partial\mathbf{f}_k/\partial\mathbf{p}\f$ (n_p x n_p) [-]
    real(kind=8), allocatable :: dfk_dp_eq(:,:)     !< Row- and column-equilibrated copy of the log-space Jacobian (overwritten by dgesv) (n_p x n_p) [-]
    real(kind=8), allocatable :: row_scale(:)       !< Row scaling factors of the Jacobian equilibration (n_p) [-]
    real(kind=8), allocatable :: col_scale(:)       !< Column scaling factors of the Jacobian equilibration (n_p) [-]
    integer(kind=4), allocatable :: ipiv_dgesv(:)   !< Pivot indices of LAPACK dgesv (n_p)
    integer(kind=4) :: info_dgesv                   !< Return code of LAPACK dgesv (0 = success)
    real(kind=8), parameter :: ln10 = log(10d0)     !< \f$\ln(10)\f$ [-]
    real(kind=8), parameter :: max_dp = 2d0         !< Maximum allowed \f$|\Delta p_j|\f$ per iteration: 2 orders of magnitude (factor 100 in concentration, PHREEQC's default step size) [-]
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
!> Local variables — stagnation detection and Picard fallback
    integer(kind=4) :: n_stag                       !< Consecutive iterations without sufficient improvement (set to n_stag_max when no step could be accepted) [-]
    integer(kind=4), parameter :: n_stag_max = 20   !< Stagnation limit before the fallback / exit [-]
    real(kind=8) :: fk_norm_prev                    !< Reference residual norm of the stagnation test [C]
    real(kind=8), parameter :: stag_rtol = 1d-3     !< Relative improvement threshold of the stagnation test [-]
    real(kind=8) :: best_fk_norm                    !< Lowest \f$\|\mathbf{f}_k\|_\infty\f$ seen since the last (re)start [C]
    logical :: fallback_applied                     !< The Picard fallback has already been tried
    logical :: CV_flag_spec                         !< The Picard fallback produced an admissible (strictly positive) state
!> Local variables — dimensions and cached parameters
    integer(kind=4) :: n_p                          !< Number of primary species (= n_v: kinetic-only) [-]
    integer(kind=4) :: n_kin                        !< Number of kinetic reactions of the current zone [-]
    integer(kind=4) :: ii                           !< Loop index [-]
    integer(kind=4) :: niter_max                    !< Maximum Newton iterations (CV_params) [-]
    real(kind=8) :: abs_tol                         !< Absolute tolerance (CV_params) [C]
    real(kind=8) :: rel_tol                         !< Relative tolerance (CV_params) [-]
    real(kind=8) :: log_abs_tol                     !< Log absolute tolerance (CV_params) [-]
    real(kind=8) :: log_rel_tol                     !< Log relative tolerance (CV_params) [-]
    real(kind=8) :: eps_d                           !< Machine epsilon [-]
    real(kind=8) :: c_norm                          !< \f$\|\mathbf{c}\|_\infty\f$ [C]
    real(kind=8) :: c_scale                         !< Relative-tolerance scale \f$\max(\|\mathbf{c}\|_\infty,\varepsilon_{\mathrm{mach}})\f$ [C]
    real(kind=8) :: dt_theta_mr_new                 !< \f$\Delta t\,\theta\,\lambda_{r}\f$ [T]
    real(kind=8), allocatable :: SkT_new(:,:)       !< \f$\mathbf{S}_k^T\f$ of the current zone (= U_SkT_prod with U = I) (n_p x n_kin) [-]
    integer(kind=4), allocatable :: ind_vas(:)      !< Zone-local -> chem_syst species index map (ind_var_act_species) (n_p) [-]
    logical :: skip_eval                            !< Reuse the rates and residual of the accepted line-search point
!> Local variables — round-off floor of the residual (termination safeguard)
    real(kind=8) :: zero_tol                        !< Numerically-zero concentration (CV_params%zero), for species with exactly zero \f$\hat c_i\f$ [C]
    real(kind=8), allocatable :: res_floor(:)       !< Round-off bound \f$\rho_i\f$ of each residual at the current iterate (n_p) [C]
    real(kind=8), allocatable :: abs_SkT(:,:)       !< \f$|\mathbf{S}_k^T|\f$ (n_p x n_kin) [-]
    real(kind=8), allocatable :: abs_c_hat(:)       !< \f$|\hat{\mathbf{c}}|\f$ (n_p) [C]
!> -----------------------------------------------------------------------
!> \section newton_kin_preproc Pre-processing: initialisation and allocation
!> -----------------------------------------------------------------------
    niter=0
    CV_flag=.false.
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    if (this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species/=n_p) &
        error stop "Newton_EI_kin_anal_ideal_opt2: kinetic-only solver called with equilibrium reactions (n_v /= n_p)"
    niter_max=this%solid_chemistry%reactive_zone%CV_params%niter_max
    abs_tol=this%solid_chemistry%reactive_zone%CV_params%abs_tol
    rel_tol=this%solid_chemistry%reactive_zone%CV_params%rel_tol
    log_abs_tol=this%solid_chemistry%reactive_zone%CV_params%log_abs_tol
    log_rel_tol=this%solid_chemistry%reactive_zone%CV_params%log_rel_tol
    zero_tol=this%solid_chemistry%reactive_zone%CV_params%zero
    eps_d=epsilon(1d0)
    n_kin=this%solid_chemistry%mineral_zone%num_minerals_kin+&
        this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
    SkT_new=this%solid_chemistry%reactive_zone%U_SkT_prod !> U = I for kinetic-only
    ind_vas=this%solid_chemistry%reactive_zone%ind_var_act_species(1:n_p)
    allocate(drk_dlogc(n_kin,this%solid_chemistry%reactive_zone%speciation_alg%num_species),drk_dlogc_zone(n_kin,n_p))
    allocate(rk_new(n_kin),rk_trial(n_kin),rk_best(n_kin),rk_avg(n_kin))
    allocate(fk(n_p),fk_trial(n_p),conc_nc_best(n_p),conc_nc_save(n_p),c_fix(n_p))
    allocate(dfk_dp(n_p,n_p),dfk_dp_eq(n_p,n_p),Delta_p(n_p),log_c(n_p),row_scale(n_p),col_scale(n_p),ipiv_dgesv(n_p))
    allocate(res_floor(n_p))
    drk_dlogc=0d0
    dfk_dp=0d0
    rk_best=0d0
    n_stag=0
    fk_norm_prev=huge(1d0)
    best_fk_norm=huge(1d0)
    fallback_applied=.false.
    skip_eval=.false.
    rk_old=this%get_rk_old_new_set() !> previous-step rates in the current kinetic set: only for the mean rate stored at convergence
    dt_theta_mr_new=Delta_t*theta*mix_ratio_r
    !> Loop-invariant factors of the round-off bound of the residual (see the convergence check)
    abs_SkT=abs(SkT_new)
    abs_c_hat=abs(c_hat)
!> -----------------------------------------------------------------------
!> \section newton_kin_process Process: Newton loop
!> -----------------------------------------------------------------------
    newton_loop: do
        niter=niter+1
    !> \subsection newton_kin_maxiter Maximum iteration check and Picard fallback
        if (niter>niter_max) then
            if (.not. fallback_applied) then
                fallback_applied=.true.
                call picard_from_best()
                if (CV_flag_spec) cycle newton_loop !> restart from the Picard state (counters reset)
            end if
            conc_nc=conc_nc_best !> return the best iterate seen
            niter=niter_max
            print *, 'Newton_EI_kin_anal_ideal_opt2 did not converge after', niter_max, &
                'iterations. Best ||fk|| =', best_fk_norm
            exit newton_loop
        end if
    !> \subsection newton_kin_eval Kinetic rates and residual at the current iterate
    !> (reused from the accepted line-search point when available: it already updated the object,
    !> evaluated the rates and the residual there)
        if (skip_eval) then
            skip_eval=.false.
            fk=fk_trial
            rk_new=rk_trial
        else
            call this%set_conc_prim_species(conc_nc(1:n_p))
            call this%set_act_aq_species()
            drk_dlogc=0d0
            call this%compute_rk_Jac_rk_anal(rk_new,drk_dlogc)
            fk=conc_nc-c_hat-dt_theta_mr_new*matmul(SkT_new,rk_new)
        end if
        fk_norm=maxval(abs(fk))
        c_norm=maxval(abs(conc_nc))
    !> \subsection newton_kin_best Best-solution tracking (initial state of the fallback; returned on failure)
        if (fk_norm<best_fk_norm) then
            best_fk_norm=fk_norm
            conc_nc_best=conc_nc
            rk_best=rk_new
        end if
    !> \subsection newton_kin_conv Convergence check (user tolerances OR round-off floor)
    !> Round-off floor (evaluated only when the user tolerances are not met): every residual within the
    !> forward-error bound of its own evaluation,
    !> rho_i = (2+n_kin) eps (|c_hat_i| + c_i + Dt theta lambda_r sum_k |Sk^T_ik| |r_k|),
    !> i.e. the sum of 2+n_kin terms forming f_k,i (the round-off of the rate evaluation itself is not
    !> included). It caps the achievable strictness of the user tolerances (an absolute tolerance below
    !> one ulp of a large concentration cannot be met).
        c_scale=max(c_norm,eps_d)
        CV_flag=(fk_norm<abs_tol .and. fk_norm<rel_tol*c_scale) .and. &
            (fk_norm>0d0 .and. log10(fk_norm)<log_abs_tol .and. (log10(fk_norm)-log10(c_scale))<log_rel_tol)
        if (.not. CV_flag) then
            res_floor=dble(2+n_kin)*eps_d*(abs_c_hat+abs(conc_nc)+abs(dt_theta_mr_new)*matmul(abs_SkT,abs(rk_new)))
            do ii=1,n_p !> a species with exactly zero total has no root at finite concentrations: accepted as absent below the numerically-zero concentration
                if (c_hat(ii)==0d0) res_floor(ii)=max(res_floor(ii),zero_tol)
            end do
            CV_flag=all(abs(fk)<=res_floor)
        end if
        if (CV_flag) then
            !> the object holds the state and the rates of this iterate; store the mean rate and the
            !> kinetic reaction amounts and contributions of the step
            rk_avg=theta*rk_new+(1d0-theta)*rk_old
            call this%set_rxn_rate_kin_avg(rk_avg)
            call this%compute_rxn_amnt_kin(Delta_t)
            call this%compute_contrib_rxn_amnt_kin(mix_ratio_r)
            exit newton_loop
        end if
    !> \subsection newton_kin_jacobian Jacobian assembly, directly in log-space (U = I, no secondary species)
    !> drk_dlogc (d r_k / d log10 c, compute_rk_Jac_rk_anal) is in chem_syst column order; its columns are
    !> permuted to the zone-local order through ind_var_act_species before the product with S_k^T:
    !> dfk/dp = ln10 diag(c) - Dt theta lambda_r S_k^T (d r_k / d log10 c), without forming any 1/c factor
        do ii=1,n_p
            drk_dlogc_zone(:,ii)=drk_dlogc(:,ind_vas(ii))
        end do
        dfk_dp=-dt_theta_mr_new*matmul(SkT_new,drk_dlogc_zone)
        do ii=1,n_p
            log_c(ii)=log10(max(conc_nc(ii),tiny(1d0)))
            dfk_dp(ii,ii)=dfk_dp(ii,ii)+ln10*conc_nc(ii)
        end do
    !> \subsection newton_kin_step Newton direction and Armijo backtracking line search in log-space
    !> (same scheme as Newton_EI_eq_kin_anal_ideal_opt2: Newton step on the equilibrated system, steepest
    !> descent -J^T f when the solve fails, the direction is not a descent direction or its line search
    !> fails; no step accepted -> stagnation)
    !> Equilibration of J Delta_p = -f: rows first (each row divided by its infinity norm), then columns
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
                    conc_nc(ii)=10d0**(log_c(ii)+alpha*Delta_p(ii))
                end do
                call this%set_conc_prim_species(conc_nc(1:n_p))
                call this%set_act_aq_species()
                drk_dlogc=0d0
                call this%compute_rk_Jac_rk_anal(rk_trial,drk_dlogc)
                fk_trial=conc_nc-c_hat-dt_theta_mr_new*matmul(SkT_new,rk_trial)
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
            !> the object must describe the restored state (the line search evaluated trial points)
            call this%set_conc_prim_species(conc_nc(1:n_p))
            call this%set_act_aq_species()
            drk_dlogc=0d0
            call this%compute_rk_Jac_rk_anal(rk_new,drk_dlogc)
        end if
    !> \subsection newton_kin_stag Stagnation detection and Picard fallback
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
                call picard_from_best()
                if (CV_flag_spec) cycle newton_loop
            end if
            conc_nc=conc_nc_best !> return the best iterate seen
            exit newton_loop !> stagnated: leave with CV_flag = .false.
        end if
    end do newton_loop
!> -----------------------------------------------------------------------
!> \section newton_kin_postproc Post-processing: deallocation
!> -----------------------------------------------------------------------
    deallocate(fk,fk_trial,drk_dlogc,drk_dlogc_zone,rk_old,rk_new,rk_trial,rk_best,rk_avg,conc_nc_best,conc_nc_save,c_fix)
    deallocate(dfk_dp,dfk_dp_eq,Delta_p,log_c,ipiv_dgesv,row_scale,col_scale,SkT_new,ind_vas)
    deallocate(res_floor,abs_SkT,abs_c_hat)

contains

    !> Picard fallback: the kinetic-only counterpart of the re-speciation fallback of
    !> Newton_EI_eq_kin_anal_ideal_opt2. There, compute_c_nc_from_u_Newton_ideal speciates the best-seen
    !> primaries towards u_hat + Delta_t*theta*lambda_r*(U*Sk^T)*rk(c_best) (kinetic source frozen at the
    !> best iterate). With U = I and no secondary species that speciation is the identity, so the state is
    !> simply c = c_hat + Delta_t*theta*lambda_r*Sk^T*rk(c_best); it is admissible iff strictly positive
    !> (a non-positive entry has no root at finite log-concentration, i.e. the inner solve would fail).
    !> On success the outer counters are reset so the Newton loop restarts from the new state; otherwise
    !> the best iterate is restored and CV_flag_spec is .false.
    subroutine picard_from_best()
        c_fix=c_hat+dt_theta_mr_new*matmul(SkT_new,rk_best)
        CV_flag_spec=all(c_fix>0d0)
        if (CV_flag_spec) then
            conc_nc=c_fix
            niter=0
            n_stag=0
            fk_norm_prev=huge(1d0)
            best_fk_norm=huge(1d0)
            skip_eval=.false.
        else
            conc_nc=conc_nc_best
        end if
    end subroutine picard_from_best
end subroutine Newton_EI_kin_anal_ideal_opt2
