!> \file Newton_EI_kin_anal_ideal_opt2.f90
!> \brief Newton method for WMA reactive mixing iteration with Euler implicit kinetic reactions.
!> Designed for kinetic-only chemical systems (no equilibrium reactions), assuming ideal solution.
!> Analytical Jacobians are used for efficiency, and option 2 is applied for averaging kinetic rates.
!>
!> \details
!> This subroutine implements a Newton method to solve the nonlinear chemical system arising
!> from a WMA reactive mixing iteration with ONLY kinetic reactions (no equilibrium). The method uses:
!> - **Euler implicit** time integration for kinetic reactions
!> - **Analytical Jacobians** for efficiency
!> - **Ideal solution** assumption (unit activity coefficients)
!> - **Option 2** for kinetic reaction rates time integration (weighted average of old and new rates)
!> - **Log-space Newton steps** for robust positivity preservation
!> - **Levenberg-Marquardt damping** with backtracking line search for global convergence
!> - **Stagnation detection** and an **explicit-predictor fallback** for robustness
!>
!> **Key difference from Newton_EI_eq_kin_anal_ideal_opt2:**
!> This version is for systems with ONLY kinetic reactions (no equilibrium reactions).
!> As a consequence:
!> - All variable activity species are primary (n_v = n_p, no secondary species)
!> - No mass-action law computation needed (no compute_c2v_from_c1_ideal)
!> - No speciation solve needed
!> - The Jacobian simplifies to df/dc = I - Dt*theta*mr*Sk^T*(drk/dc)
!>
!> **Mathematical Formulation:**
!>
!> The nonlinear residual to drive to zero is (U = I for kinetic-only):
!> \f[
!>   \mathbf{f}_k = \mathbf{c}
!>                  - \hat{\mathbf{c}}
!>                  - \Delta t \bigl[
!>                      (1-\theta)\,\lambda_{r}^{\mathrm{old}}\,
!>                        \mathbf{S}_k^{T,\mathrm{old}}\,\mathbf{r}_k^{\mathrm{old}}
!>                    + \theta\,\lambda_{r}^{\mathrm{new}}\,
!>                        \mathbf{S}_k^{T,\mathrm{new}}\,\mathbf{r}_k^{\mathrm{new}}
!>                    \bigr]
!>                  = \mathbf{0}
!> \f]
!>
!> **Log-Space Jacobian Transformation:**
!>
!> Defining \f$ p_j = \log_{10}(c_j) \f$, the Newton system is solved for
!> \f$ \Delta\mathbf{p} \f$ instead of \f$ \Delta\mathbf{c} \f$:
!> \f[
!>   \frac{\partial \mathbf{f}_k}{\partial p_j}
!>     = \ln(10)\,c_j\,\frac{\partial \mathbf{f}_k}{\partial c_j}
!> \f]
!> which guarantees positivity of concentrations via
!> \f$ c_j^{\mathrm{new}} = 10^{p_j + \alpha\,\Delta p_j} \f$.
!>
!> \param[in,out] this   Aqueous chemistry object at current time step
!> \param[in]     c_hat  Variable activity species concentrations after mixing (n_p) [C]
!> \param[in]     mix_ratio_r_old  Mixing ratio of old kinetic reaction amounts [-]
!> \param[in]     mix_ratio_r_new  Mixing ratio of new kinetic reaction amounts [-]
!> \param[in]     Delta_t  Time step size [T]
!> \param[in]     theta    Reaction time weighting factor [0,1] [-]
!> \param[in,out] conc_nc  Variable activity species concentrations (n_p) [C]
!> \param[out]    niter    Number of Newton iterations performed [-]
!> \param[out]    CV_flag  Convergence flag

subroutine Newton_EI_kin_anal_ideal_opt2(this,c_hat,mix_ratio_r_old,mix_ratio_r_new,Delta_t,theta,conc_nc,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c  !< Import aqueous chemistry class
    use vectors_m, only: inf_norm_vec_real              !< Import infinity norm utility
    implicit none                                       !< Enforce explicit variable declarations
!> Arguments
    class(aqueous_chemistry_c) :: this              !< Aqueous chemistry object at current time step [-]
    real(kind=8), intent(in) :: c_hat(:)            !< Component concentrations after mixing (n_p) [C]
    real(kind=8), intent(in) :: mix_ratio_r_old     !< Mixing ratio of old kinetic reaction amounts [-]
    real(kind=8), intent(in) :: mix_ratio_r_new     !< Mixing ratio of new kinetic reaction amounts [-]
    real(kind=8), intent(in) :: Delta_t             !< Time step size [T]
    real(kind=8), intent(in) :: theta               !< Reaction time weighting factor [0,1] [-]
    real(kind=8), intent(inout) :: conc_nc(:)       !< Variable activity species concentrations (n_p); initial guess on entry, solution on exit [C]
    integer(kind=4), intent(out) :: niter           !< Number of Newton iterations performed [-]
    logical, intent(out) :: CV_flag                 !< Convergence flag: .FALSE. = not converged, .TRUE. = converged
!> Local variables -- kinetic rates and Jacobians
    real(kind=8), allocatable :: rk_new(:)          !< New kinetic reaction rates (n_kin) [1/T]
    real(kind=8), allocatable :: drk_dc(:,:)        !< Jacobian drk/dc (n_kin x n_species) [1/(T*C)]
    real(kind=8), allocatable :: dfk_dc(:,:)        !< Jacobian of residual w.r.t. concentrations (n_p x n_p) [-]
    real(kind=8), allocatable :: rk_old(:)          !< Kinetic reaction rates from previous time step (n_kin) [1/T]
    real(kind=8), allocatable :: rk_avg(:)          !< Weighted average rates (n_kin) [1/T]
!> Local variables -- concentrations and residual
    real(kind=8), allocatable :: fk(:)              !< Newton residual vector (n_p) [C]
    real(kind=8), allocatable :: conc_nc_save(:)    !< Saved concentrations before LM/line-search trial (n_p) [C]
    real(kind=8), allocatable :: conc_nc_best(:)    !< Best concentrations seen during Newton (n_p) [C]
!> Local variables -- log-space Newton system
    real(kind=8), allocatable :: log_c(:)           !< Log10 of concentrations (n_p) [-]
    real(kind=8), allocatable :: Delta_p(:)         !< Newton step in log-space (n_p) [-]
    real(kind=8), allocatable :: dfk_dp(:,:)        !< Log-space Jacobian (n_p x n_p) [-]
    real(kind=8), allocatable :: dfk_dp_LM(:,:)     !< LM-damped log-space Jacobian (n_p x n_p) [-]
    real(kind=8), parameter :: ln10 = log(10d0)     !< ln(10) ~ 2.302585 [-]
    real(kind=8), parameter :: max_dp = 10d0        !< Maximum allowed |Delta_p| per iteration [-]
    real(kind=8) :: scale_factor                    !< Uniform scaling factor for Delta_p [-]
!> Local variables -- Levenberg-Marquardt and line search
    real(kind=8) :: lambda_LM                       !< LM damping parameter [-]
    integer(kind=4) :: i_LM                         !< LM inner loop counter [-]
    integer(kind=4), parameter :: max_LM_tries = 10 !< Maximum LM damping attempts [-]
    real(kind=8) :: alpha                           !< Backtracking line search step length [-]
    integer(kind=4) :: i_LS                         !< Line search iteration counter [-]
    integer(kind=4), parameter :: max_LS_tries = 10 !< Maximum backtracking iterations [-]
    logical :: step_accepted                        !< Whether a step was accepted
    real(kind=8), allocatable :: fk_trial(:)        !< Trial residual vector (n_p) [C]
    real(kind=8), allocatable :: rk_trial(:)        !< Trial kinetic rates (n_kin) [1/T]
    real(kind=8) :: fk_norm                         !< Current residual infinity norm [C]
    real(kind=8) :: fk_trial_norm                   !< Trial residual infinity norm [C]
!> Local variables -- stagnation detection and fallback
    integer(kind=4) :: n_stag                       !< Consecutive stagnation iterations [-]
    integer(kind=4), parameter :: n_stag_max = 20   !< Maximum stagnation before exit [-]
    real(kind=8) :: fk_norm_prev                    !< Previous best residual norm [C]
    real(kind=8), parameter :: stag_rtol = 1d-3     !< Stagnation relative improvement threshold [-]
    real(kind=8) :: best_fk_norm                    !< Lowest ||fk||_inf seen during Newton [C]
    logical :: predictor_applied                    !< Whether explicit-predictor fallback has been tried
    real(kind=8), allocatable :: c_pred(:)          !< Predicted component concentrations from forward-Euler (n_p) [C]
    real(kind=8), save :: predictor_damp = 5d-1     !< Adaptive damping for Euler predictor [-]
    real(kind=8), parameter :: predictor_damp_min = 1d-2  !< Minimum predictor damping [-]
    real(kind=8), parameter :: predictor_damp_max = 1d0   !< Maximum predictor damping [-]
    real(kind=8), parameter :: predictor_grow = 1.2d0     !< Growth factor on success [-]
    real(kind=8), parameter :: predictor_shrink = 0.5d0   !< Shrink factor on failure [-]
!> Local variables -- dimensions and misc
    integer(kind=4) :: n_p                          !< Number of primary species (= n_v for kinetic-only) [-]
    integer(kind=4) :: i                            !< Loop index [-]
    logical :: sing_flag                            !< Singularity flag for dgesv
    integer(kind=4), save :: dbg_calls = 0          !< Persistent call counter for diagnostics
    integer(kind=4), allocatable :: ipiv_dgesv(:)   !< Pivot array for LAPACK dgesv (n_p)
    integer(kind=4) :: info_dgesv                   !< Return code from dgesv
    real(kind=8) :: c_norm                          !< ||c||_inf [C]
    logical :: skip_eval                            !< Skip recomputation when accepted step can be reused
!> Cached loop-invariant parameters
    integer(kind=4) :: niter_max                    !< Max Newton iterations (cached)
    integer(kind=4) :: n_kin                        !< Number of kinetic reactions (cached)
    real(kind=8) :: abs_tol                         !< Absolute tolerance (cached)
    real(kind=8) :: rel_tol                         !< Relative tolerance (cached)
    real(kind=8) :: log_abs_tol                     !< Log absolute tolerance (cached)
    real(kind=8) :: log_rel_tol                     !< Log relative tolerance (cached)
    real(kind=8) :: eps_d                           !< Machine epsilon (cached)
    real(kind=8) :: sqrt_eps_d                      !< sqrt(eps) (cached)
    real(kind=8) :: dt_theta_mr_new                 !< Cached Dt*theta*lambda_r^new [-]
    real(kind=8), allocatable :: SkT_old(:,:)       !< S_k^T at old time (n_p x n_kin) (U=I for kin-only)
    real(kind=8), allocatable :: SkT_new(:,:)       !< S_k^T at new time (n_p x n_kin) (U=I for kin-only)
    real(kind=8), allocatable :: source_old(:)      !< Constant old kinetic source Dt*(1-theta)*mr_old*SkT_old*rk_old (cached, n_p)
    real(kind=8), allocatable :: SkT_new_rk_old(:)  !< Cached S_k^T * rk_old for predictor (n_p) [C]
!> -----------------------------------------------------------------------
!> Pre-processing: initialisation and allocation
!> -----------------------------------------------------------------------
    niter=0
    CV_flag=.false.
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    !> Cache loop-invariant parameters from deep struct dereferences
    niter_max=this%solid_chemistry%reactive_zone%CV_params%niter_max
    abs_tol=this%solid_chemistry%reactive_zone%CV_params%abs_tol
    rel_tol=this%solid_chemistry%reactive_zone%CV_params%rel_tol
    log_abs_tol=this%solid_chemistry%reactive_zone%CV_params%log_abs_tol
    log_rel_tol=this%solid_chemistry%reactive_zone%CV_params%log_rel_tol
    eps_d=epsilon(1d0)
    sqrt_eps_d=sqrt(eps_d)
    n_kin=this%solid_chemistry%mineral_zone%num_minerals_kin+&
        this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
    !> Cache S_k^T products (U=I for kinetic-only, so U*S_k^T = S_k^T)
    SkT_old=this%solid_chemistry_old%reactive_zone%U_SkT_prod
    SkT_new=this%solid_chemistry%reactive_zone%U_SkT_prod
    !> First-call diagnostics
    if (dbg_calls == 0) then
        continue
    end if
    !> Allocate persistent work arrays
    allocate(dfk_dc(n_p,n_p),&
        drk_dc(n_kin,this%solid_chemistry%reactive_zone%speciation_alg%num_species),&
        rk_new(n_kin),fk(n_p),conc_nc_best(n_p))
    allocate(dfk_dp(n_p,n_p), dfk_dp_LM(n_p,n_p), Delta_p(n_p), log_c(n_p))
    allocate(rk_trial(n_kin), fk_trial(n_p))
    allocate(ipiv_dgesv(n_p))
    drk_dc=0d0
    dfk_dc=0d0
    lambda_LM=1d-3
    n_stag=0
    fk_norm_prev=huge(1d0)
!> -----------------------------------------------------------------------
!> Process: retrieve old rates, precompute constants, and Newton loop
!> -----------------------------------------------------------------------
    !> Retrieve old kinetic reaction rates from the previous time step
    rk_old=this%get_rk_old()
    !> Precompute constant kinetic source term (does not change during Newton iterations)
    !> source_old = Delta_t * (1-theta) * mix_ratio_r_old * SkT_old * rk_old  (U=I)
    source_old = Delta_t*(1d0-theta)*mix_ratio_r_old*matmul(SkT_old, rk_old)
    dt_theta_mr_new = Delta_t*theta*mix_ratio_r_new
    SkT_new_rk_old = matmul(SkT_new, rk_old)
    !> Save original initial guess for recovery via explicit-predictor fallback
    predictor_applied = .false.
    best_fk_norm = huge(1d0)
    skip_eval = .false.

    !> -----------------------------------------------------------------------
    !> Newton iteration loop
    !> -----------------------------------------------------------------------
    newton_loop: do
        niter=niter+1
    !> Maximum iteration check and explicit-predictor fallback
        if (niter>niter_max) then
            if (.not. predictor_applied) then
                !> Forward-Euler predictor: approximate r_k^new ~ r_k^old
                predictor_applied = .true.
                allocate(c_pred(n_p))
                c_pred = c_hat + predictor_damp * ( &
                    source_old + &
                    dt_theta_mr_new * SkT_new_rk_old )
                !> For kinetic-only: direct assignment (no speciation solve needed since n_v = n_p)
                do i = 1, n_p
                    conc_nc(i) = max(c_pred(i), tiny(1d0))
                end do
                call this%set_conc_prim_species(conc_nc(1:n_p))
                call this%set_act_aq_species()
                deallocate(c_pred)
                niter = 0
                n_stag = 0
                fk_norm_prev = huge(1d0)
                lambda_LM = 1d-3
                best_fk_norm = huge(1d0)
                skip_eval = .false.
                cycle newton_loop
            end if
            !> Accept if best residual is within conditioning noise
            if (best_fk_norm<=sqrt_eps_d*max(c_norm,1d0)) then
                CV_flag=.true.
                conc_nc=conc_nc_best
                call this%set_conc_prim_species(conc_nc(1:n_p))
                call this%set_act_aq_species()
                call this%compute_Rk(theta,Delta_t)
            else
                print *, 'Newton did not converge after', niter_max, &
                    'iterations. Best ||fk|| =', best_fk_norm
            end if
            exit newton_loop
        end if
    !> Set species and activities (no secondary species for kinetic-only)
        if (skip_eval) then
            skip_eval = .false.
        else
            call this%set_conc_prim_species(conc_nc(1:n_p))
            call this%set_act_aq_species()
            drk_dc=0d0
            call this%compute_rk_Jac_rk_anal(rk_new,drk_dc)
        end if
    !> Weighted average of kinetic reaction rates
        rk_avg=theta*rk_new+(1d0-theta)*rk_old
        call this%set_rk_mean(rk_avg)
    !> Newton residual: f_k = c - c_hat - Dt*[(1-theta)*mr_old*SkT_old*rk_old + theta*mr_new*SkT_new*rk_new]  (U=I)
        fk=conc_nc-c_hat-&
            source_old-&
            dt_theta_mr_new*matmul(SkT_new,rk_new)
    !> Convergence check (dual arithmetic OR logarithmic criteria)
        fk_norm = inf_norm_vec_real(fk)
        c_norm = inf_norm_vec_real(conc_nc)
        !> Track best solution
        if (fk_norm < best_fk_norm) then
            best_fk_norm = fk_norm
            conc_nc_best = conc_nc
        end if
        !> Triple convergence check (same as eq_kin version)
        if ((fk_norm<=eps_d*max(c_norm,1d0)) .or. &
            (fk_norm<abs_tol .and. &
            fk_norm<rel_tol*max(c_norm,1d0)) .or. &
            (fk_norm>0d0 .and. &
            (log10(fk_norm)<log_abs_tol .and. &
            (log10(fk_norm)-log10(max(c_norm,1d0)))< &
            log_rel_tol))) then
            CV_flag=.true.
            call this%compute_Rk(theta,Delta_t)
            exit newton_loop
        else
        !> Jacobian assembly for kinetic-only (U=I, no dc2/dc1 chain rule)
        !> dfk/dc = I - Dt*theta*lambda_r^new * S_k^T * (drk/dc)
            dfk_dc = -dt_theta_mr_new * matmul(SkT_new, drk_dc(:,1:n_p))
            do i = 1, n_p
                dfk_dc(i,i) = 1d0 + dfk_dc(i,i)
            end do
            !> Transform to log-space Jacobian:
            !> (dfk/dp)_{:,j} = ln(10) * c_j * (dfk/dc_j)
            do i = 1, n_p
                log_c(i) = log10(max(conc_nc(i), tiny(1d0)))
                dfk_dp(:,i) = ln10 * conc_nc(i) * dfk_dc(:,i)
            end do
        !> Levenberg-Marquardt with backtracking line search in log-space
            conc_nc_save = conc_nc
            step_accepted = .false.
            do i_LM = 1, max_LM_tries
                !> Apply LM damping: copy undamped Jacobian and add lambda*I to diagonal
                dfk_dp_LM = dfk_dp
                do i = 1, n_p
                    dfk_dp_LM(i, i) = dfk_dp_LM(i, i) + lambda_LM
                end do
                !> Solve (dfk/dp + lambda*I) * Delta_p = -fk via LAPACK dgesv
                Delta_p = -fk
                call dgesv(n_p, 1, dfk_dp_LM, n_p, ipiv_dgesv, Delta_p, n_p, info_dgesv)
                !> Validate solution: check for dgesv failure and NaN/overflow
                sing_flag = (info_dgesv /= 0)
                do i = 1, n_p
                    if (Delta_p(i) /= Delta_p(i) .or. abs(Delta_p(i)) > 1d20) then
                        sing_flag = .true.
                        exit
                    end if
                end do
                if (sing_flag) then
                    lambda_LM = min(lambda_LM * 10d0, 1d8)
                    cycle
                end if
                !> Uniform step bounding: scale Delta_p so ||Delta_p||_inf <= max_dp
                scale_factor = 1d0
                do i = 1, n_p
                    if (abs(Delta_p(i)) > 0d0) then
                        scale_factor = min(scale_factor, max_dp / abs(Delta_p(i)))
                    end if
                end do
                if (scale_factor < 1d0) Delta_p = Delta_p * scale_factor
                !> Backtracking line search in log-space
                alpha = 1d0
                do i_LS = 1, max_LS_tries
                    !> Update species in log-space: c_j^trial = 10^(p_j + alpha*Delta_p_j) (guarantees positivity)
                    do i = 1, n_p
                        conc_nc(i) = 10d0**(log_c(i) + alpha * Delta_p(i))
                    end do
                    call this%set_conc_prim_species(conc_nc(1:n_p))
                    !> No secondary species to compute (kinetic-only)
                    call this%set_act_aq_species()
                    drk_dc = 0d0
                    call this%compute_rk_Jac_rk_anal(rk_trial, drk_dc)
                    !> Evaluate trial residual (U=I)
                    fk_trial = conc_nc - c_hat - &
                        source_old - &
                        dt_theta_mr_new*matmul(SkT_new, rk_trial)
                    fk_trial_norm = inf_norm_vec_real(fk_trial)
                    if (fk_trial_norm < fk_norm) then
                        step_accepted = .true.
                        rk_new = rk_trial
                        fk = fk_trial
                        skip_eval = .true.
                        lambda_LM = max(lambda_LM / 10d0, 1d-12)
                        exit
                    end if
                    alpha = alpha * 0.5d0
                end do
                if (step_accepted) exit
                !> Line search failed: restore concentrations and increase damping
                conc_nc = conc_nc_save
                lambda_LM = min(lambda_LM * 10d0, 1d8)
            end do
            !> All LM attempts exhausted: restore saved state
            if (.not. step_accepted) then
                conc_nc = conc_nc_save
                lambda_LM = 1d-3
            end if
        !> Stagnation detection
            if (fk_norm < fk_norm_prev * (1d0 - stag_rtol)) then
                n_stag = 0
                fk_norm_prev = fk_norm
            else
                n_stag = n_stag + 1
                if (n_stag >= n_stag_max) then
                    if (.not. predictor_applied) then
                        !> Stagnation detected: try explicit predictor before giving up
                        predictor_applied = .true.
                        allocate(c_pred(n_p))
                        c_pred = c_hat + predictor_damp * ( &
                            source_old + &
                            dt_theta_mr_new * SkT_new_rk_old )
                        do i = 1, n_p
                            conc_nc(i) = max(c_pred(i), tiny(1d0))
                        end do
                        call this%set_conc_prim_species(conc_nc(1:n_p))
                        call this%set_act_aq_species()
                        deallocate(c_pred)
                        niter = 0
                        n_stag = 0
                        fk_norm_prev = huge(1d0)
                        lambda_LM = 1d-3
                        best_fk_norm = huge(1d0)
                        skip_eval = .false.
                        cycle newton_loop
                    end if
                    !> Accept if best residual within conditioning noise
                    if (best_fk_norm<=sqrt_eps_d*max(c_norm,1d0)) then
                        CV_flag=.true.
                        conc_nc=conc_nc_best
                        call this%set_conc_prim_species(conc_nc(1:n_p))
                        call this%set_act_aq_species()
                        call this%compute_Rk(theta,Delta_t)
                    end if
                    exit newton_loop
                end if
            end if
        end if
    end do newton_loop
!> -----------------------------------------------------------------------
!> Post-processing: cleanup and deallocation
!> -----------------------------------------------------------------------
    !> Adaptive predictor damping update
    if (predictor_applied) then
        if (CV_flag) then
            predictor_damp = min(predictor_damp * predictor_grow, predictor_damp_max)
        else
            predictor_damp = max(predictor_damp * predictor_shrink, predictor_damp_min)
        end if
    end if
    dbg_calls = dbg_calls + 1
    deallocate(fk,dfk_dc,drk_dc,rk_old,rk_new,rk_avg,conc_nc_best)
    deallocate(dfk_dp,dfk_dp_LM,Delta_p,log_c,rk_trial,fk_trial,ipiv_dgesv)
    deallocate(SkT_old,SkT_new,source_old,SkT_new_rk_old)
end subroutine Newton_EI_kin_anal_ideal_opt2
