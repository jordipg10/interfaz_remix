!> \file reactive_mixing_iter_RKF45_eq_kin_anal_opt2.f90
!> \brief Reactive mixing iteration with adaptive Runge-Kutta-Fehlberg 4(5) kinetic reactions for
!>        equilibrium-kinetic systems (option 2 for the mean kinetic rate).
!> \details Explicit-integrator counterpart of `reactive_mixing_iter_EI_eq_kin_anal_opt2` with the SAME
!> signature (so both can be targets of the `p_solver` procedure pointer of the reactive transport
!> solvers). It performs a single reactive mixing iteration of the WMA (Water Mixing Approach) with
!> equilibrium and kinetic reactions, solving for the variable activity species concentrations
!> \f$ \mathbf{c}_{v} \f$ and component concentrations \f$ \mathbf{u} \f$ after reactive mixing with
!> `RKF45_eq_kin_anal_opt2` instead of the globalised Newton method `Newton_EI_eq_kin_anal_opt2`.
!>
!> **Reaction step:** the mixed concentrations \f$\hat{\mathbf{c}}_v\f$ already carry every explicit
!> term of the step (the neighbours' \f$\theta\f$-weighted rates and the previous-step \f$(1-\theta)\f$
!> rates, folded in by the mixer for the full \f$\Delta t\f$); only the self term
!> \f$\theta\,\lambda_r(\mathbf{U}\mathbf{S}_k^T)\mathbf{r}_k(\mathbf{c}_v)\f$ is integrated here, as the
!> component ODE
!> \f[
!>   \frac{d\mathbf{u}}{dt} = \theta\,\lambda_r(\mathbf{U}\mathbf{S}_k^T)\,\mathbf{r}_k(\mathbf{c}_v(\mathbf{u})),
!>   \qquad \mathbf{u}(0) = \hat{\mathbf{u}} = \mathbf{U}\hat{\mathbf{c}}_v ,
!> \f]
!> with the embedded RKF45 pair and adaptive internal sub-steps (error control with
!> `CV_params%abs_tol` / `rel_tol`), each stage requiring an equilibrium speciation.
!>
!> **Retry strategy: initial guesses, then time step halving (sub-stepping):**
!>
!> The same outer structure as the Newton driver is kept. The reaction step is advanced from
!> \f$\hat{\mathbf{u}}\f$ in sub-steps of length \f$\Delta t_{bis} = \Delta t/2^{k_{div}}\f$, starting
!> with \f$k_{div} = 0\f$ (one call over the full \f$\Delta t\f$). For every sub-step the integrator is
!> tried with a sequence of initial guesses of the primaries (which only seed the first speciation of
!> the integrator), each at most once:
!> 1. closest downstream target water (if `c1_downstream` is given),
!> 2. linear extrapolation with \f$ \mu \in \{0, 0.25, 0.5, 0.75, 1\} \f$ from the primaries at the
!>    start of the previous and of the current sub-step (initially `c1_old` and \f$\mathbf{c}_1^{(k)}\f$),
!> 3. the mixed concentrations \f$ \hat{\mathbf{c}} \f$.
!> If all of them fail (the integrator exhausts its internal sub-step budget or a speciation does not
!> converge), the sub-step is halved (\f$k_{div} = k_{div}+1\f$) and retried from the same starting
!> state \f$\mathbf{u}^{(j-1)}\f$; the sub-steps already completed are kept. When
!> `CV_params%k_div_max` halvings have been used up the run stops with diagnostics (fatal error).
!>
!> When \f$k_{div} > 0\f$ the end-of-step kinetic quantities are made consistent with the full step:
!> the object's kinetic rates are set to the time-weighted mean of the sub-step mean rates,
!> \f[
!>   \bar{\mathbf{r}}_k = \frac{1}{\Delta t}\sum_j \Delta t_{bis}^{(j)}\,\bar{\mathbf{r}}_k^{(j)}
!>   \qquad\Longrightarrow\qquad
!>   \mathbf{U}\mathbf{c}_v - \hat{\mathbf{u}} = \Delta t\,\theta\,\lambda_r(\mathbf{U}\mathbf{S}_k^T)\bar{\mathbf{r}}_k ,
!> \f]
!> which is what compute_Re and the mass balance assume, and the mean rate
!> \f$\theta\bar{\mathbf{r}}_k + (1-\theta)\mathbf{r}_k^{old}\f$, the kinetic reaction amounts
!> (\f$\Delta t\f$) and their species contributions (\f$\lambda_r\f$) are recomputed for the full step.
!>
!> **Initial guess strategy:** as in the Newton driver, the extrapolated guess
!> \f$ c_{1,j}^{(0)} = (1 + \mu)\,c_{1,j}^{(k)} - \mu\,c_{1,j}^{(k-1)} \f$ is clamped to
!> \f$ c_{1,j}^{(k)} \f$ if negative, and a floor replaces guesses below
!> \f$ \max(10^{-16},\,10^{-8}\hat{c}_j) \f$ (while \f$ \hat{c}_j > 10^{-16} \f$) by \f$ \hat{c}_j \f$.
!>
!> \param[in,out] this Aqueous chemistry object at current time step; on exit its variable activity
!>                     concentrations, aqueous activities, kinetic rates, mean kinetic rates and kinetic
!>                     reaction amounts are those of the solution over the full time step
!> \param[in] c1_old Primary species concentrations at previous time step (n_p)
!> \param[in] c_hat Variable activity species concentrations after mixing (n_v); includes the
!>                  previous-step \f$(1-\theta)\f$ kinetic term
!> \param[in] mix_ratio_r Reaction mixing ratio of this target water \f$ \lambda_{r} \f$ [-]
!> \param[in] Delta_t Full time step \f$ \Delta t \f$ [T]; halved internally into sub-steps only when
!>                    the integrator fails for every initial guess
!> \param[in] theta Reaction time weighting factor \f$ \theta \in [0,1] \f$ [-]
!> \param[in,out] conc_nc Variable activity species concentrations (already allocated, size n_v);
!>                        current concentrations \f$ \mathbf{c}_v^{(k)} \f$ on entry, solution on exit
!> \param[in,out] conc_comp Component concentrations (already allocated, size n_p); overwritten on exit with \f$ \mathbf{U}\mathbf{c}_{v} \f$
!> \param[in] c1_downstream Optional primary concentrations of the closest downstream target water (n_p), used as first initial guess
!>
!> \pre The aqueous chemistry object must have valid speciation algebra, solid chemistry, kinetic
!>      reactions and convergence parameters.
!> \pre `conc_nc` must be pre-allocated with size \f$ n_v \f$ and contain physically meaningful values.
!> \post On exit, `conc_nc` contains the end-of-step species concentrations and `conc_comp` the
!>       corresponding component concentrations.
!>
!> \sa RKF45_eq_kin_anal_opt2, reactive_mixing_iter_EI_eq_kin_anal_opt2, initialise_iterative_method

subroutine reactive_mixing_iter_RKF45_eq_kin_anal_opt2(this,c1_old,c_hat,mix_ratio_r,&
    Delta_t,theta,conc_nc,conc_comp,c1_downstream)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, initialise_iterative_method  !< Import aqueous chemistry class and initial guess helper
    implicit none                                       !< Enforce explicit variable declarations
!> Arguments
    class(aqueous_chemistry_c) :: this                 !< Aqueous chemistry object at current time step [-]
    real(kind=8), intent(in) :: c1_old(:)              !< Primary species concentrations at previous time step \f$ \mathbf{c}_1^{(k-1)} \f$ (n_p)
    real(kind=8), intent(in) :: c_hat(:)               !< Variable activity species concentrations after mixing \f$ \hat{\mathbf{c}}_v \f$ (n_v)
    real(kind=8), intent(in) :: mix_ratio_r            !< Reaction mixing ratio of this target water \f$ \lambda_{r} \f$ [-]
    real(kind=8), intent(in) :: Delta_t                !< Time step \f$ \Delta t \f$ [T]
    real(kind=8), intent(in) :: theta                  !< Reaction time weighting factor \f$ \theta \in [0,1] \f$ [-]
    real(kind=8), intent(inout) :: conc_nc(:)          !< Variable activity species concentrations (already allocated) (n_v)
    real(kind=8), intent(inout) :: conc_comp(:)        !< Component concentrations (already allocated) (n_p)
    real(kind=8), intent(in), optional :: c1_downstream(:) !< Primary concentrations of closest downstream target water (n_p)
!> Local variables
    real(kind=8), allocatable :: c1_sub(:)             !< Primary species concentrations at the start of the current sub-step (initially \f$ \mathbf{c}_1^{(k)} \f$) (n_p)
    real(kind=8), allocatable :: c1_sub_old(:)         !< Primary species concentrations at the start of the previous sub-step (initially `c1_old`) (n_p)
    real(kind=8), allocatable :: u_hat(:)              !< Component concentrations after mixing \f$ \hat{\mathbf{u}} = \mathbf{U}\hat{\mathbf{c}}_v \f$ (n_p)
    real(kind=8), allocatable :: u_sub(:)              !< Component concentrations at the start of the current sub-step \f$ \mathbf{u}^{(j-1)} \f$ (n_p)
    real(kind=8), allocatable :: rk_sum(:)             !< \f$ \sum_j \Delta t_{bis}^{(j)}\,\bar{\mathbf{r}}_k^{(j)} \f$ over the completed sub-steps, layout of get_rk_new (n_kin) [C]
    real(kind=8), allocatable :: rk_mean(:)            !< Time-weighted mean kinetic rate of the step \f$ \bar{\mathbf{r}}_k \f$ (n_kin) [1/T]
    real(kind=8), allocatable :: rk_avg(:)             !< Mean kinetic rate \f$ \theta\bar{\mathbf{r}}_k + (1-\theta)\mathbf{r}_k^{old} \f$ (n_kin) [1/T]
    real(kind=8), allocatable :: u_fix(:)              !< Components of the speciation polish \f$ \hat{\mathbf{u}} + \Delta t\,\theta\,\lambda_r(\mathbf{U}\mathbf{S}_k^T)\bar{\mathbf{r}}_k \f$ (n_p)
    real(kind=8), allocatable :: c1_pol(:)             !< Initial guess of the speciation polish: sub-stepped primaries (n_p)
    real(kind=8), allocatable :: c2v_pol(:)            !< Initial guess of the speciation polish: sub-stepped secondaries (n_v - n_p)
    logical :: CV_spec                                 !< Convergence flag of the speciation polish [-]
    integer(kind=4) :: n_p                             !< Number of primary species [-]
    integer(kind=4) :: niter                           !< Number of accepted RKF45 sub-steps in the last integrator call [-]
    integer(kind=4) :: j_sp                            !< Loop index over primary species for concentration floor [-]
    integer(kind=4) :: k_div                           !< Number of time step halvings performed [-]
    integer(kind=4) :: k_div_max                       !< Maximum number of time step halvings (CV_params%k_div_max) [-]
    integer(kind=8) :: n_sub_done                      !< Completed sub-steps, counted in units of the current \f$ \Delta t_{bis} \f$ [-]
    integer(kind=8) :: n_sub_total                     !< Sub-steps of length \f$ \Delta t_{bis} \f$ in the full step (\f$ 2^{k_{div}} \f$) [-]
    real(kind=8) :: Delta_t_bis                        !< Current sub-step length \f$ \Delta t/2^{k_{div}} \f$ [T]
    real(kind=8) :: mu                                 !< Initialisation parameter \f$ \mu \in [0,1] \f$ [-]
    real(kind=8), parameter :: conc_floor = 1d-16      !< Concentration floor: threshold below which a species is considered effectively zero [-]
    real(kind=8), parameter :: rel_floor = 1d-8        !< Relative floor of the initial guess: a guess below rel_floor * c_hat is replaced by the mixed value [-]
    logical :: CV_flag                                 !< Convergence flag from the integrator: .true. if the step was completed [-]
    logical :: use_downstream                          !< Next attempt uses the downstream-water initial guess [-]
    logical :: c_hat_tried                             !< Whether the c_hat-based initial guess has been attempted as a last resort [-]

    associate(comp_mat => this%solid_chemistry%reactive_zone%speciation_alg%comp_mat)  !< \f$ \mathbf{U} \f$ (n_p x n_v) by reference: no per-call copy
!> -----------------------------------------------------------------------
!> \section rmix_rkf_preproc Pre-processing: dimensions, mixed components, snapshot of current c1, sub-step state
!> -----------------------------------------------------------------------
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species    !< Extract number of primary species from speciation algebra
    k_div_max=this%solid_chemistry%reactive_zone%CV_params%k_div_max         !< Maximum number of time step halvings
    allocate(u_hat(n_p),u_sub(n_p),c1_sub(n_p),c1_sub_old(n_p))              !< Allocate mixed/sub-step component vectors and primary concentration snapshots
    c1_sub=conc_nc(1:n_p)                                                    !< Snapshot \f$ \mathbf{c}_1^{(k)} \f$ before conc_nc is overwritten by the guesses
    c1_sub_old=c1_old                                                        !< Extrapolation partner of the first sub-step: previous time step
    u_hat=matmul(comp_mat,c_hat)                                             !< \f$ \hat{\mathbf{u}} = \mathbf{U}\hat{\mathbf{c}}_v \f$ (the (1-theta) kinetic term is already in c_hat)
    u_sub=u_hat                                                              !< The first sub-step starts from the mixed components
    rk_sum=this%get_rk_new()                                                 !< Size and layout of the current kinetic set (values overwritten)
    rk_sum=0d0                                                               !< No sub-step completed yet
    k_div=0                                                                  !< No time step halving yet
    Delta_t_bis=Delta_t                                                      !< First attempt: one sub-step over the full step
    n_sub_done=0                                                             !< No sub-step completed yet
    n_sub_total=1                                                            !< One sub-step of length Delta_t_bis covers the step
    CV_flag=.false.                                                          !< Nothing solved yet

!> -----------------------------------------------------------------------
!> \section rmix_rkf_process Process: sub-steps, each integrated by RKF45 over a sequence of initial guesses;
!> the sub-step is halved (at most k_div_max times) when every initial guess fails
!> -----------------------------------------------------------------------
    sub_loop: do while (n_sub_done<n_sub_total)                              !< Until the completed sub-steps cover the full time step
        !> Order of the guesses: downstream water (if given), extrapolation with mu = 0, 0.25, ..., 1,
        !> mixed concentrations c_hat. Each is tried once over the current sub-step.
        use_downstream=present(c1_downstream)                                !< Downstream guess goes first when available
        mu=-0.25d0                                                           !< First extrapolation attempt uses mu = 0
        c_hat_tried=.false.                                                  !< c_hat fallback not yet attempted
        try_loop: do                                                         !< One integrator attempt per initial guess
            if (use_downstream) then                                         !< Guess 1: closest downstream target water
                use_downstream=.false.                                       !< Consume the option
                conc_nc(1:n_p)=c1_downstream(1:n_p)                          !< Use downstream target water concentrations as initial guess
            else if (mu<1d0) then                                            !< Guess 2: extrapolation, mu = 0, 0.25, ..., 1 (exact in binary)
                mu=mu+0.25d0                                                 !< Next mu
                call initialise_iterative_method(c1_sub_old,c1_sub,mu,conc_nc(1:n_p)) !< \f$ c_{1,j}^{(0)} = (1+\mu)\,c_{1,j}^{sub} - \mu\,c_{1,j}^{sub,old} \f$
            else if (.not. c_hat_tried) then                                 !< Guess 3: mixed concentrations
                c_hat_tried=.true.                                           !< Consume the option
                conc_nc(1:n_p)=c_hat(1:n_p)                                  !< Use mixed concentrations as initial guess
            else                                                             !< Every initial guess has failed
                exit try_loop                                                !< Leave with CV_flag = .false.
            end if                                                           !< End initial guess selection
            !> Floor near-zero species: replace with mixed values when the guess is negligible
            !> (absolute floor conc_floor and relative floor rel_floor * c_hat)
            do j_sp=1,n_p                                                    !< Loop over each primary species
                if (conc_nc(j_sp) < max(conc_floor, c_hat(j_sp)*rel_floor) .and. & !< Initial guess below absolute or relative floor?
                    c_hat(j_sp) > conc_floor) then                            !< Mixed value is not negligible?
                    conc_nc(j_sp)=c_hat(j_sp)                                 !< Use mixed concentration as initial guess
                end if                                                        !< End floor check
            end do                                                            !< End initial guess floor loop
            call this%RKF45_eq_kin_anal_opt2(u_sub,mix_ratio_r,Delta_t_bis,theta,conc_nc,niter,CV_flag) !< Adaptive RKF45 over the current sub-step
            if (CV_flag) exit try_loop                                       !< Step completed: leave the initial guess loop
        end do try_loop                                                      !< End initial guess loop
        if (CV_flag) then                                                    !< Sub-step completed
            n_sub_done=n_sub_done+1                                          !< One more sub-step of the current length
            rk_sum=rk_sum+Delta_t_bis*this%get_rk_new()                      !< Time-weighted accumulation of the sub-step mean kinetic rates
            u_sub=matmul(comp_mat,conc_nc)                                   !< \f$ \mathbf{u}^{(j)} = \mathbf{U}\mathbf{c}^{(j)} \f$: start of the next sub-step
            c1_sub_old=c1_sub                                                !< Shift the extrapolation partners
            c1_sub=conc_nc(1:n_p)                                            !< Primaries at the start of the next sub-step
        else                                                                 !< Every initial guess failed over this sub-step
            if (k_div>=k_div_max) exit sub_loop                              !< Halving budget exhausted: fatal error below
            k_div=k_div+1                                                    !< One more time step halving
            Delta_t_bis=0.5d0*Delta_t_bis                                    !< Halve the sub-step; retry from the same start state u_sub
            n_sub_done=2*n_sub_done                                          !< Completed time re-expressed in units of the new sub-step
            n_sub_total=2*n_sub_total                                        !< Full step in units of the new sub-step
            print *, "reactive_mixing_iter_RKF45_eq_kin_anal_opt2: water", this%id, ": time step halved, k_div =", k_div, &
                " Delta_t_bis =", Delta_t_bis, " sub-steps done", n_sub_done, " of", n_sub_total
        end if                                                               !< End sub-step outcome
    end do sub_loop                                                          !< End sub-step loop
    if (.not. CV_flag) then                                                  !< Unrecoverable non-convergence
        !> conc_nc holds the last successfully speciated state of the last integrator attempt
        print *, "FATAL: reactive_mixing_iter_RKF45_eq_kin_anal_opt2: RKF45 failed for every initial guess after", &
            k_div, "time step halvings (k_div_max =", k_div_max, ")"
        print *, "  Delta_t=", Delta_t, " Delta_t_bis=", Delta_t_bis, " sub-steps completed=", n_sub_done, " of", n_sub_total
        print *, "  ||U c||_inf=", maxval(abs(matmul(comp_mat,conc_nc))), &
            " abs_tol=", this%solid_chemistry%reactive_zone%CV_params%abs_tol, &
            " rel_tol=", this%solid_chemistry%reactive_zone%CV_params%rel_tol
        print *, "  theta=", theta, " lambda_r=", mix_ratio_r
        print *, "  c1_old=", c1_old                                          !< Print old primary concentrations
        print *, "  c_hat=", c_hat                                            !< Print mixed species concentrations
        print *, "  u_sub=", u_sub                                            !< Print components at the start of the failed sub-step
        print *, "  conc_nc=", conc_nc(1:n_p)                                 !< Print last primary concentrations
        print *, "  rk_old=", this%get_rk_old_new_set()                       !< Print previous-step kinetic rates (current kinetic set, as rk_new)
        print *, "  rk_new=", this%get_rk_new()                               !< Print kinetic reaction rates at the last state
        error stop "Too many time step divisions in reactive_mixing_iter_RKF45_eq_kin_anal_opt2: RKF45 did not complete the step for any initial guess. Check CV parameters."
    end if                                                                    !< End failure handling

!> -----------------------------------------------------------------------
!> \section rmix_rkf_postproc Post-processing: full-step kinetic quantities, update object, component concentrations
!> -----------------------------------------------------------------------
    if (k_div>0) then                                                        !< Sub-stepped: the integrator stored sub-step quantities; make them those of the full step
        rk_mean=rk_sum/Delta_t                                               !< \f$ \bar{\mathbf{r}}_k \f$: U c - u_hat = Delta_t*theta*lambda_r*(U Sk^T)*rk_mean up to the sum of the sub-step speciation residuals
        call set_rk_new(rk_mean)                                             !< End-of-step rates of the object (read by compute_Re, the mass balance and the next step as rk_old)
        !> Speciation polish with the kinetic source frozen at rk_mean: the sub-step speciation residuals
        !> add up to O(n_sub*abs_tol) in U c - u_hat - Delta_t*theta*lambda_r*(U Sk^T)*rk_mean, which the
        !> single-step consistency check of compute_Re rejects. The round-off-level speciation Newton
        !> started from the sub-stepped state removes it in one or two iterations; if it fails, its best
        !> iterate (never worse than the sub-stepped state) is kept.
        u_fix=u_hat+Delta_t*theta*mix_ratio_r*matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,rk_mean)
        c1_pol=conc_nc(1:n_p)                                                !< Sub-stepped primaries as initial guess
        c2v_pol=conc_nc(n_p+1:)                                              !< Sub-stepped secondaries as initial guess
        call this%compute_c_nc_from_u_Newton(c1_pol,c2v_pol,u_fix,conc_nc,niter,CV_spec) !< Exact speciation of u_fix (object state = conc_nc on exit)
        rk_avg=theta*rk_mean+(1d0-theta)*this%get_rk_old_new_set()           !< Mean rate of the full step
        call this%set_rxn_rate_kin_avg(rk_avg)                               !< Store the mean rate
        call this%compute_rxn_amnt_kin(Delta_t)                              !< Kinetic reaction amounts over the full step
        call this%compute_contrib_rxn_amnt_kin(mix_ratio_r)                  !< Their contribution to the species mass balance
    end if                                                                    !< End full-step consistency
    call this%set_conc_var_act_species(conc_nc)                              !< Synchronise the object's variable activity concentrations with the solution
    call this%compute_log_act_coeffs_aq_chem() !> we compute logarithmic activity coefficients for all species
    call this%compute_activities() !> we compute activities
    call this%compute_pH()
    call this%compute_log_act_coeff_wat() !> we compute log_10 activity coefficient of water
    call this%compute_salinity()
    call this%compute_alkalinity() !> we compute the alkalinity of the solution
    call this%compute_saturation_mins() !> we compute the saturation indices of all minerals at the end-of-step state
    conc_comp=matmul(comp_mat,conc_nc)                                       !< Compute final component concentrations \f$ \mathbf{u} = \mathbf{U}\mathbf{c}_v \f$
    end associate                                                            !< End comp_mat association
    deallocate(u_hat,u_sub,c1_sub,c1_sub_old,rk_sum)                         !< Free locally allocated arrays
    if (allocated(rk_mean)) deallocate(rk_mean,rk_avg,u_fix,c1_pol,c2v_pol)  !< Allocated only when the step was subdivided

contains

    !> Stores rk as the end-of-step kinetic rates of the object (counterpart of set_rk_old, layout of
    !> get_rk_new: aqueous kinetic reactions of the chemical system, then the kinetic minerals of the
    !> mineral zone)
    subroutine set_rk_new(rk)
        real(kind=8), intent(in) :: rk(:)                                    !< Kinetic reaction rates in the layout of get_rk_new (n_kin) [1/T]
        integer(kind=4) :: i                                                 !< Loop index [-]
        integer(kind=4) :: n_aq_kin                                          !< Number of aqueous kinetic reactions [-]
        n_aq_kin=this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
        do i=1,n_aq_kin
            this%rxn_rate_kin(i)=rk(i)                                       !< Aqueous: rates carried by the water
        end do
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            this%solid_chemistry%rxn_rate_kin(i)=rk(n_aq_kin+i)             !< Minerals: rates of the cell's solid chemistry
        end do
        !> Faltan gases
    end subroutine set_rk_new
end subroutine                                                               !< End reactive_mixing_iter_RKF45_eq_kin_anal_opt2
