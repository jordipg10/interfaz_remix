!> \file RKF45_eq_kin_anal_opt2.f90
!> \brief Adaptive Runge-Kutta-Fehlberg 4(5) integration of the kinetic source for WMA reactive mixing
!>        with equilibrium AND kinetic reactions (option 2 for the mean kinetic rate).
!>
!> \details
!> Explicit counterpart of `Newton_EI_eq_kin_anal_opt2`: instead of solving the Euler-implicit nonlinear
!> system for the end-of-step concentrations with a globalised Newton method, the component ODE
!> \f[
!>   \frac{d\mathbf{u}}{dt} = \mathbf{F}(\mathbf{u}) = \theta\,\lambda_r\,(\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}\,
!>                            \mathbf{r}_k\bigl(\mathbf{c}_v(\mathbf{u})\bigr),
!>   \qquad \mathbf{u}(0)=\hat{\mathbf{u}}
!> \f]
!> is integrated over \f$[0,\Delta t]\f$ with the embedded Runge-Kutta-Fehlberg 4(5) pair and adaptive
!> internal sub-steps \f$h\f$. Each of the six stages requires one equilibrium speciation
!> \f$\mathbf{u}\mapsto\mathbf{c}_v(\mathbf{u})\f$ (`compute_c_nc_from_u_Newton`, initial guess =
!> concentrations of the previous stage) followed by one evaluation of the kinetic rates
!> (`compute_rxn_rate_kin`; no Jacobian is formed).
!>
!> **Butcher tableau (Fehlberg):**
!> \f[
!> \begin{array}{c|cccccc}
!>  0 \\
!>  1/4   & 1/4 \\
!>  3/8   & 3/32 & 9/32 \\
!>  12/13 & 1932/2197 & -7200/2197 & 7296/2197 \\
!>  1     & 439/216 & -8 & 3680/513 & -845/4104 \\
!>  1/2   & -8/27 & 2 & -3544/2565 & 1859/4104 & -11/40 \\ \hline
!>  \mathbf{b}^{(5)} & 16/135 & 0 & 6656/12825 & 28561/56430 & -9/50 & 2/55 \\
!>  \mathbf{b}^{(4)} & 25/216 & 0 & 1408/2565 & 2197/4104 & -1/5 & 0
!> \end{array}
!> \f]
!> After each trial sub-step the local error estimate
!> \f$ \mathrm{err} = \max_i |u^{(5)}_i-u^{(4)}_i| / \mathrm{tol}_i \f$ with
!> \f$ \mathrm{tol}_i = \varepsilon_{\mathrm{abs}} + \varepsilon_{\mathrm{rel}}\max(|u_i|,|u_i^{(5)}|) \f$
!> (`CV_params%abs_tol`, `CV_params%rel_tol`) decides acceptance (\f$\mathrm{err}\le 1\f$). The step is
!> advanced with the 5th-order solution (local extrapolation) and the next sub-step is
!> \f$ h_{\mathrm{new}} = h\,\min\bigl(f_{\max},\max\bigl(f_{\min},\,0.9\,\mathrm{err}^{-1/5}\bigr)\bigr) \f$
!> (\f$f_{\min}=0.1,\;f_{\max}=5\f$), clipped so as not to overshoot \f$\Delta t\f$; a rejected sub-step is
!> repeated with \f$ h\,\max(f_{\min},0.9\,\mathrm{err}^{-1/4}) \f$. A stage whose speciation fails or
!> produces non-finite components is treated as a rejected sub-step with \f$h\f$ halved. The number of
!> accepted plus rejected sub-steps is limited to `max_sub_steps`; exhausting it (or \f$h\f$ falling below
!> \f$\varepsilon_{\mathrm{mach}}\,\Delta t\f$) ends the call with `CV_flag = .false.`.
!>
!> **End-of-step kinetic rate:** by construction
!> \f$ \mathbf{u}^{\mathrm{new}}-\hat{\mathbf{u}} = \theta\,\lambda_r\,(\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}
!> \sum_{\mathrm{sub-steps}} h\sum_{i=1}^{6} b^{(5)}_i\,\mathbf{r}_{k,i} \f$, so the quadrature-weighted mean rate
!> \f$ \bar{\mathbf{r}}_k^{RK} = \frac{1}{\Delta t}\sum_{\mathrm{sub-steps}} h\sum_i b^{(5)}_i\,\mathbf{r}_{k,i} \f$
!> satisfies exactly the identity used downstream by the reactive mixing driver and the mass-balance check
!> (`compute_Re`), \f$\mathbf{U}\mathbf{c}_v-\hat{\mathbf{u}} = \Delta t\,\theta\,\lambda_r(\mathbf{U}\mathbf{S}_k^T)\mathbf{r}_k^{\mathrm{new}}\f$,
!> with \f$\mathbf{r}_k^{\mathrm{new}} = \bar{\mathbf{r}}_k^{RK}\f$ (stored in the object by this routine);
!> \f$\hat{\mathbf{u}}\f$ already carries the \f$(1-\theta)\f$ previous-step kinetic contribution. With
!> \f$\theta=1\f$ the scheme is the plain RKF45 integration of the reaction ODE after mixing. The
!> end-of-step concentrations \f$\mathbf{c}_v\f$ are obtained by a final speciation of \f$\mathbf{u}^{\mathrm{new}}\f$.
!>
!> **Failure handling:** `CV_flag` is `.false.` when the sub-step budget is exhausted, \f$h\f$ underflows,
!> or the final speciation fails. The caller (`reactive_mixing_iter_EI_eq_kin_anal_opt2`) then retries
!> with another initial guess or halves \f$\Delta t\f$, exactly as for the Newton solver.
!>
!> **Post-processing at convergence:** the RK mean rate is stored as the end-of-step rate of the object
!> (layout of `get_rk_new`), the mean kinetic rate \f$\bar{\mathbf{r}}_k = \theta\,\bar{\mathbf{r}}_k^{RK}
!> + (1-\theta)\,\mathbf{r}_k^{\mathrm{old}}\f$ is stored (`set_rxn_rate_kin_avg`), and the kinetic
!> reaction amounts and their species contributions are computed (`compute_rxn_amnt_kin`,
!> `compute_contrib_rxn_amnt_kin`).
!>
!> \param[in,out] this   Aqueous chemistry object at current time step (class `aqueous_chemistry_c`);
!>                        on exit its primaries, secondaries and activities are those of the last
!>                        speciation performed (the end-of-step state when CV_flag is TRUE) and its
!>                        kinetic rates are the RKF45 mean rates
!> \param[in]     u_hat  Component concentrations after mixing \f$\hat{\mathbf{u}}\f$ [C] (size \f$n_p\f$);
!>                        already includes the previous-step (1-theta) kinetic contribution
!> \param[in]     mix_ratio_r  Mixing ratio of kinetic reaction amounts \f$\lambda_{r}\f$ [-]
!> \param[in]     Delta_t  Time step size \f$\Delta t\f$ [T]
!> \param[in]     theta    Reaction time weighting factor \f$\theta \in [0,1]\f$ [-]
!> \param[in,out] conc_nc  Variable activity species concentrations \f$\mathbf{c}_{v}\f$ [C] (size
!>                          \f$n_v\f$; primaries used as initial guess of the first speciation on entry,
!>                          end-of-step concentrations on exit if CV_flag is TRUE, otherwise the last
!>                          successfully speciated state)
!> \param[out]    niter    Number of accepted RKF45 sub-steps [-]
!> \param[out]    CV_flag  Convergence flag: `.FALSE.` = integration failed, `.TRUE.` = step completed
!>
!> \pre  The chemistry object `this` must have valid `solid_chemistry`, reactive-zone stoichiometry
!>       (`U_SkT_prod`), kinetic reactions and speciation algorithm data.
!> \pre  `conc_nc` must be pre-allocated with size \f$n_v\f$ and contain a strictly positive initial guess
!>       of the primary species.
!>
!> \warning Explicit scheme: stiff kinetics force very small internal sub-steps; the sub-step budget
!>          then fails the call and the caller's sub-stepping takes over.
!>
!> \sa Newton_EI_eq_kin_anal_opt2, compute_c_nc_from_u_Newton
!> \sa aqueous_chemistry_c::compute_rxn_rate_kin
!>
!> \author Jordi
!> \ingroup chemistry

subroutine RKF45_eq_kin_anal_opt2(this,u_hat,mix_ratio_r,Delta_t,theta,conc_nc,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c  !< Import aqueous chemistry class for reactive mixing methods
    implicit none                                       !< Enforce explicit variable declarations
!> Arguments
    class(aqueous_chemistry_c) :: this              !< Aqueous chemistry object at current time step [-]
    real(kind=8), intent(in) :: u_hat(:)            !< Component concentrations after mixing (n_p) [C]; already includes the previous-step (1-theta) kinetic contribution
    real(kind=8), intent(in) :: mix_ratio_r         !< Mixing ratio of new kinetic reaction amounts \f$\lambda_{r}\f$ [-]
    real(kind=8), intent(in) :: Delta_t             !< Time step size \f$\Delta t\f$ [T]
    real(kind=8), intent(in) :: theta               !< Reaction time weighting factor \f$\theta \in [0,1]\f$ [-]
    real(kind=8), intent(inout) :: conc_nc(:)       !< Variable activity species concentrations \f$\mathbf{c}_{v}\f$ (n_v); initial guess on entry, solution on exit [C]
    integer(kind=4), intent(out) :: niter           !< Number of accepted RKF45 sub-steps [-]
    logical, intent(out) :: CV_flag                 !< Convergence flag: .FALSE. = integration failed, .TRUE. = step completed
!> Local variables — Runge-Kutta-Fehlberg 4(5) tableau
    integer(kind=4), parameter :: n_stages = 6      !< Number of RK stages [-]
    real(kind=8), parameter :: a21 = 1d0/4d0
    real(kind=8), parameter :: a31 = 3d0/32d0,        a32 = 9d0/32d0
    real(kind=8), parameter :: a41 = 1932d0/2197d0,   a42 = -7200d0/2197d0, a43 = 7296d0/2197d0
    real(kind=8), parameter :: a51 = 439d0/216d0,     a52 = -8d0,           a53 = 3680d0/513d0,    a54 = -845d0/4104d0
    real(kind=8), parameter :: a61 = -8d0/27d0,       a62 = 2d0,            a63 = -3544d0/2565d0,  a64 = 1859d0/4104d0, a65 = -11d0/40d0
    real(kind=8), parameter :: b5(n_stages) = [16d0/135d0, 0d0, 6656d0/12825d0, 28561d0/56430d0, -9d0/50d0, 2d0/55d0] !< 5th-order weights [-]
    real(kind=8), parameter :: b4(n_stages) = [25d0/216d0, 0d0, 1408d0/2565d0,  2197d0/4104d0,   -1d0/5d0,  0d0]      !< 4th-order weights [-]
!> Local variables — step-size control
    real(kind=8), parameter :: safety = 0.9d0       !< Safety factor of the step-size update [-]
    real(kind=8), parameter :: fac_min = 0.1d0      !< Lower bound of the step-size ratio [-]
    real(kind=8), parameter :: fac_max = 5d0        !< Upper bound of the step-size ratio [-]
    integer(kind=4), parameter :: max_sub_steps = 1000 !< Budget of trial (accepted + rejected) sub-steps [-]
    real(kind=8) :: h                               !< Current internal sub-step [T]
    real(kind=8) :: t                               !< Time elapsed within the step [T]
    real(kind=8) :: err                             !< Scaled local error estimate (accept if <= 1) [-]
    real(kind=8) :: fac                             !< Step-size ratio [-]
    integer(kind=4) :: n_trials                     !< Trial sub-steps performed [-]
    logical :: stage_ok                             !< All six stages of the trial sub-step were evaluated
!> Local variables — stage quantities
    real(kind=8), allocatable :: k_stage(:,:)       !< Stage slopes \f$\mathbf{k}_i = \mathbf{F}(\mathbf{u}_i)\f$ (n_p x n_stages) [C/T]
    real(kind=8), allocatable :: rk_stage(:,:)      !< Kinetic rates at each stage \f$\mathbf{r}_{k,i}\f$ (n_kin x n_stages) [1/T]
    real(kind=8), allocatable :: u_cur(:)           !< Components at the start of the current sub-step (n_p) [C]
    real(kind=8), allocatable :: u_stage(:)         !< Components at the current stage (n_p) [C]
    real(kind=8), allocatable :: u_4(:)             !< 4th-order end-of-sub-step components (n_p) [C]
    real(kind=8), allocatable :: u_5(:)             !< 5th-order end-of-sub-step components (n_p) [C]
    real(kind=8), allocatable :: tol_vec(:)         !< Component-wise error tolerance (n_p) [C]
    real(kind=8), allocatable :: conc_stage(:)      !< Speciated concentrations at the current stage; initial guess of the next one (n_v) [C]
    real(kind=8), allocatable :: conc_cur(:)        !< Speciated concentrations at the start of the current sub-step (restored on rejection) (n_v) [C]
    real(kind=8), allocatable :: c2v_ig(:)          !< Initial guess of the secondary variable activity species for the speciation (n_e) [C]
    real(kind=8), allocatable :: log_gamma_var_act_ideal(:) !< Log10 activity coefficients of the variable activity species (ideal: all zero) [-]
!> Local variables — kinetic rates stored at the end of the step
    real(kind=8), allocatable :: rk_old(:)          !< Kinetic reaction rates from the previous time step \f$\mathbf{r}_k^{\mathrm{old}}\f$ (n_kin) [1/T]
    real(kind=8), allocatable :: rk_sum(:)          !< \f$\sum h\sum_i b^{(5)}_i\mathbf{r}_{k,i}\f$ over the accepted sub-steps (n_kin) [-]
    real(kind=8), allocatable :: rk_mean(:)         !< Quadrature-weighted mean of the stage rates \f$\bar{\mathbf{r}}_k^{RK} = \mathrm{rk\_sum}/\Delta t\f$ (n_kin) [1/T]
    real(kind=8), allocatable :: rk_avg(:)          !< Mean rates \f$\theta\,\bar{\mathbf{r}}_k^{RK} + (1-\theta)\,\mathbf{r}_k^{\mathrm{old}}\f$ (n_kin) [1/T]
!> Local variables — dimensions and cached parameters
    integer(kind=4) :: n_p                          !< Number of primary species [-]
    integer(kind=4) :: n_e                          !< Number of equilibrium reactions [-]
    integer(kind=4) :: n_v                          !< Number of variable activity species (= n_p + n_e) [-]
    integer(kind=4) :: n_kin                        !< Number of kinetic reactions of the current zone [-]
    integer(kind=4) :: i_st                         !< Stage index [-]
    integer(kind=4) :: niter_spec                   !< Inner speciation Newton iterations of one stage [-]
    logical :: CV_spec                              !< Inner speciation convergence flag
    real(kind=8) :: abs_tol                         !< Absolute tolerance of the local error (CV_params) [C]
    real(kind=8) :: rel_tol                         !< Relative tolerance of the local error (CV_params) [-]
    real(kind=8) :: theta_mr_new                    !< \f$\theta\,\lambda_{r}\f$ [-]
    real(kind=8), allocatable :: U_SkT_new(:,:)     !< \f$(\mathbf{U}\mathbf{S}_k^T)^{\mathrm{new}}\f$ (n_p x n_kin) [-]
!> -----------------------------------------------------------------------
!> \section rkf45_preproc Pre-processing: initialisation and allocation
!> -----------------------------------------------------------------------
    niter=0
    CV_flag=.false.
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_v=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    abs_tol=this%solid_chemistry%reactive_zone%CV_params%abs_tol
    rel_tol=this%solid_chemistry%reactive_zone%CV_params%rel_tol
    n_kin=this%solid_chemistry%mineral_zone%num_minerals_kin+&
        this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
    U_SkT_new=this%solid_chemistry%reactive_zone%U_SkT_prod
    allocate(k_stage(n_p,n_stages),rk_stage(n_kin,n_stages),rk_sum(n_kin),rk_mean(n_kin),rk_avg(n_kin))
    allocate(u_cur(n_p),u_stage(n_p),u_4(n_p),u_5(n_p),tol_vec(n_p),conc_stage(n_v),conc_cur(n_v),c2v_ig(n_e))
    log_gamma_var_act_ideal=this%compute_log_act_coeffs_var_act_ideal() !> ideal log activity coefficients: initial guess of the secondaries
    rk_old=this%get_rk_old_new_set() !> previous-step rates in the current kinetic set: only for the mean rate stored at the end
    theta_mr_new=theta*mix_ratio_r
    conc_stage=conc_nc !> primaries of the entry guess seed the first speciation
    conc_cur=conc_nc
    u_cur=u_hat
    rk_sum=0d0
    t=0d0
    h=Delta_t !> first trial: one sub-step over the full step
    n_trials=0
!> -----------------------------------------------------------------------
!> \section rkf45_process Process: adaptive RKF45 sub-steps over [0, Delta_t]
!> -----------------------------------------------------------------------
    sub_loop: do while (t<Delta_t)
        n_trials=n_trials+1
        if (n_trials>max_sub_steps .or. h<epsilon(1d0)*Delta_t) then
            print *, 'RKF45_eq_kin_anal_opt2 failed: sub-steps =', n_trials-1, ' h =', h, ' t/Delta_t =', t/Delta_t
            exit sub_loop !> leave with CV_flag = .false.
        end if
        if (t+h>Delta_t) h=Delta_t-t !> do not overshoot the end of the step
    !> \subsection rkf45_stages Six stages: speciation + rate evaluation at each one
        stage_ok=.true.
        stage_loop: do i_st=1,n_stages
            select case (i_st)
            case (1)
                u_stage=u_cur
            case (2)
                u_stage=u_cur+h*a21*k_stage(:,1)
            case (3)
                u_stage=u_cur+h*(a31*k_stage(:,1)+a32*k_stage(:,2))
            case (4)
                u_stage=u_cur+h*(a41*k_stage(:,1)+a42*k_stage(:,2)+a43*k_stage(:,3))
            case (5)
                u_stage=u_cur+h*(a51*k_stage(:,1)+a52*k_stage(:,2)+a53*k_stage(:,3)+a54*k_stage(:,4))
            case (6)
                u_stage=u_cur+h*(a61*k_stage(:,1)+a62*k_stage(:,2)+a63*k_stage(:,3)+a64*k_stage(:,4)+a65*k_stage(:,5))
            end select
            if (.not. is_finite(u_stage)) then
                stage_ok=.false.
                exit stage_loop
            end if
            call speciate(u_stage) !> initial guess: previous stage (start-of-sub-step state for stage 1)
            if (.not. CV_spec) then
                stage_ok=.false.
                exit stage_loop
            end if
            call this%compute_rxn_rate_kin(rk_stage(:,i_st)) !> kinetic rates at the stage state (no Jacobian needed)
            k_stage(:,i_st)=theta_mr_new*matmul(U_SkT_new,rk_stage(:,i_st))
        end do stage_loop
    !> \subsection rkf45_err Embedded error estimate and step-size control
        if (stage_ok) then
            u_5=u_cur+h*matmul(k_stage,b5)
            u_4=u_cur+h*matmul(k_stage,b4)
            stage_ok=is_finite(u_5) .and. is_finite(u_4)
        end if
        if (.not. stage_ok) then !> infeasible stage: reject and halve the sub-step
            conc_stage=conc_cur !> re-seed the speciation from the start-of-sub-step state
            h=0.5d0*h
            cycle sub_loop
        end if
        tol_vec=abs_tol+rel_tol*max(abs(u_cur),abs(u_5))
        err=maxval(abs(u_5-u_4)/tol_vec)
        if (err<=1d0) then !> accept: advance with the 5th-order solution (local extrapolation)
            niter=niter+1
            rk_sum=rk_sum+h*matmul(rk_stage,b5)
            t=t+h
            u_cur=u_5
            !> conc_stage holds the speciation of stage 6 (c = 1/2): a good initial guess of the next sub-step
            conc_cur=conc_stage
            if (err>0d0) then
                fac=min(fac_max,max(fac_min,safety*err**(-0.2d0)))
            else
                fac=fac_max
            end if
            h=h*fac
        else !> reject: retry from the same state with a smaller sub-step
            conc_stage=conc_cur
            h=h*max(fac_min,safety*err**(-0.25d0))
        end if
    end do sub_loop
!> -----------------------------------------------------------------------
!> \section rkf45_final End-of-step speciation and kinetic post-processing
!> -----------------------------------------------------------------------
    if (t>=Delta_t) then
        call speciate(u_cur) !> end-of-step state: exact speciation of u_new, initial guess = last stage
        if (CV_spec) then
            rk_mean=rk_sum/Delta_t !> U c - u_hat = Delta_t*theta*lambda_r*(U Sk^T)*rk_mean exactly (up to the speciation tolerance)
            call set_rk_new(rk_mean)
            rk_avg=theta*rk_mean+(1d0-theta)*rk_old
            call this%set_rxn_rate_kin_avg(rk_avg)
            call this%compute_rxn_amnt_kin(Delta_t)
            call this%compute_contrib_rxn_amnt_kin(mix_ratio_r)
            CV_flag=.true.
        else
            conc_stage=conc_cur
        end if
    else
        conc_stage=conc_cur
    end if
    conc_nc=conc_stage !> end-of-step concentrations on success, last accepted state otherwise
!> -----------------------------------------------------------------------
!> \section rkf45_postproc Post-processing: deallocation
!> -----------------------------------------------------------------------
    deallocate(k_stage,rk_stage,u_cur,u_stage,u_4,u_5,tol_vec,conc_stage,conc_cur,c2v_ig,log_gamma_var_act_ideal)
    deallocate(rk_old,rk_sum,rk_mean,rk_avg,U_SkT_new)

contains

    !> .true. if every entry of v is finite (no NaN, no Inf)
    pure logical function is_finite(v)
        real(kind=8), intent(in) :: v(:)
        is_finite=all(v==v) .and. all(abs(v)<=huge(1d0))
    end function is_finite

    !> Equilibrium speciation of the components u: primaries of conc_stage as initial guess, ideal
    !> secondaries derived from them; conc_stage is overwritten with the result (the object state is
    !> consistent with conc_stage on exit).
    subroutine speciate(u)
        real(kind=8), intent(in) :: u(:)            !< Component concentrations to speciate (n_p) [C]
        real(kind=8), allocatable :: c1_ig(:)       !< Initial guess of the primaries (n_p) [C]
        c1_ig=conc_stage(1:n_p)
        call this%set_conc_prim_species(c1_ig)
        call this%compute_c2v_from_c1_ideal(c1_ig,log_gamma_var_act_ideal,c2v_ig)
        call this%compute_c_nc_from_u_Newton(c1_ig,c2v_ig,u,conc_stage,niter_spec,CV_spec)
        deallocate(c1_ig)
    end subroutine speciate

    !> Stores rk as the end-of-step kinetic rates of the object (layout of get_rk_new: aqueous kinetic
    !> reactions of the chemical system, then the kinetic minerals of the mineral zone)
    subroutine set_rk_new(rk)
        real(kind=8), intent(in) :: rk(:)           !< Kinetic reaction rates in the layout of get_rk_new (n_kin) [1/T]
        integer(kind=4) :: i                        !< Loop index [-]
        integer(kind=4) :: n_aq_kin                 !< Number of aqueous kinetic reactions [-]
        n_aq_kin=this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
        do i=1,n_aq_kin
            this%rxn_rate_kin(i)=rk(i)              !< Aqueous: rates carried by the water
        end do
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            this%solid_chemistry%rxn_rate_kin(i)=rk(n_aq_kin+i) !< Minerals: rates of the cell's solid chemistry
        end do
        !> Faltan gases
    end subroutine set_rk_new
end subroutine RKF45_eq_kin_anal_opt2
