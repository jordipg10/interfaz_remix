!> \file compute_c_nc_from_u_Newton_ideal.f90
!> \brief Computes variable activity species concentrations from component concentrations using Newton-Raphson method
!> \details This subroutine solves the nonlinear speciation problem for ideal solutions with a log-space
!> Newton-Raphson algorithm globalised by an Armijo backtracking line search (the scheme of
!> Newton_EI_eq_ideal), with stagnation detection, best-solution tracking and a round-off-level
!> convergence criterion. Given component concentrations and an initial guess for primary species, it
!> computes variable activity species concentrations by iteratively solving:
!> \f[
!> \mathbf{r}(\mathbf{c}_1) = \mathbf{U}_1 \mathbf{c}_1 + \mathbf{U}_2 \cdot \mathbf{c}_{2,v} - \mathbf{u} = \mathbf{0}
!> \f]
!> where:
!> - \f$ \mathbf{c}_1 \f$ = primary species concentrations
!> - \f$ \mathbf{c}_{2,v} \f$ = secondary variable activity species concentrations (mass action law,
!>   \f$ \log_{10}\mathbf{c}_{2,v} = \mathbf{S}^*_{e,v,1}\log_{10}\mathbf{c}_1 + \log_{10}\mathbf{K}^* \f$)
!> - \f$ \mathbf{U}_1 \f$ = component matrix for primary species
!> - \f$ \mathbf{U}_2 \f$ = component matrix for secondary variable activity species
!> - \f$ \mathbf{u} \f$ = component concentrations (input)
!>
!> Assumptions:
!> - Ideal solution behavior (the object's \f$\log_{10}\gamma\f$ are read once and held fixed: zero)
!> - Initial guess for primary concentrations is provided as input (strictly positive)
!>
!> **Newton iteration in log-space** (\f$ \mathbf{p} = \log_{10}\mathbf{c}_1 \f$, removes the \f$1/c_1\f$
!> singularity of the Jacobian and keeps the concentrations positive by construction):
!> \f[
!> \mathbf{J}_p\,\Delta\mathbf{p} = -\mathbf{r}, \qquad \mathbf{J}_p(:,j) = \ln(10)\,c_{1,j}\,\mathbf{J}(:,j),
!> \qquad \mathbf{J} = \mathbf{U}_1 + \mathbf{U}_2\,\partial\mathbf{c}_{2,v}/\partial\mathbf{c}_1,
!> \qquad \mathbf{p}^{(i+1)} = \mathbf{p}^{(i)} + \alpha\,\Delta\mathbf{p}^{(i)}
!> \f]
!> The column-equilibrated system is solved with LAPACK `dgesv`.
!>
!> **Globalisation — Armijo backtracking line search** (no Levenberg-Marquardt damping): the step length
!> \f$\alpha\f$ along \f$\Delta\mathbf{p}\f$ is chosen by backtracking (\f$\alpha = 1, 1/2, 1/4, \ldots\f$,
!> at most `max_LS_tries` trials) on the merit function
!> \f$\phi(\alpha)=\tfrac12\|\mathbf{W}\mathbf{r}(\mathbf{p}+\alpha\Delta\mathbf{p})\|_2^2\f$,
!> \f$\mathbf{W}=\mathrm{diag}(1/\rho_i)\f$ with \f$\rho_i\f$ the round-off bound of each component residual at
!> the current iterate (frozen during the line search; see the convergence criterion below), until the
!> sufficient-decrease condition \f$\phi(\alpha)\le\phi(0)+c\,\alpha\,\phi'(0)\f$, \f$c=10^{-4}\f$, holds, where
!> \f$\phi'(0)=(\mathbf{W}\mathbf{r})^T\mathbf{W}\mathbf{J}_p\Delta\mathbf{p}\f$ (\f$=-\|\mathbf{W}\mathbf{r}\|_2^2\f$
!> for the exact Newton direction). The weighting makes the merit function measure exactly what the
!> convergence test requires: with the plain \f$\tfrac12\|\mathbf{r}\|_2^2\f$ the residuals of the large
!> components, once at round-off level, are pure noise that hides the progress still needed on a small
!> component (e.g. a 1e-6 total next to 1 M totals), so that only half steps get accepted and the
!> iteration crawls linearly. Before the line search \f$\Delta\mathbf{p}\f$ is scaled uniformly so that
!> \f$\|\Delta\mathbf{p}\|_\infty\le\f$ `max_dp` decades (direction preserved, \f$\phi'(0)\f$ re-evaluated).
!> If the solve fails, the Newton direction is not a descent direction (\f$\phi'(0)\ge0\f$) or its line
!> search fails, the steepest-descent direction \f$-\mathbf{J}_p^T\mathbf{W}^2\mathbf{r}\f$ of \f$\phi\f$ is tried once.
!> If no step is accepted the iterate is kept and the iteration counts as stagnated at once (repeating it
!> would give the same result).
!>
!> **Convergence criterion (the strictest achievable in double precision).** The iterate is accepted when
!> every component residual lies within the rigorous round-off bound of its own evaluation:
!> \f[
!>   |r_i| \le \rho_i = (n_v+1)\,\varepsilon\,\Bigl(|u_i| + \sum_j |U_{ij}|\,c_j\Bigr)
!>                      + \sum_{j>n_p} |U_{ij}|\,c_j\,\delta_j ,
!> \f]
!> where \f$\varepsilon\f$ = `epsilon(1d0)`, the first term is the forward-error bound of the
!> \f$(n_v+1)\f$-term inner product forming \f$r_i\f$ (Higham, Accuracy and Stability of Numerical
!> Algorithms, eq. 3.5; any summation order) and \f$\delta_j\f$ bounds the relative round-off of the
!> mass-action evaluation of secondary species \f$j\f$, \f$c_j = 10^{e_j}\f$,
!> \f$ e_j = \sum_k S^*_{jk}(\log_{10}c_{1,k} + \log_{10}\gamma_k) + \log_{10}K^*_j - \log_{10}\gamma_j \f$:
!> \f[
!>   \delta_j = \varepsilon\Bigl[1 + \ln(10)\,(n_p+3)\Bigl(\sum_k |S^*_{jk}|\,(|\log_{10}c_{1,k}| + |\log_{10}\gamma_k|)
!>              + |\log_{10}K^*_j| + |\log_{10}\gamma_j|\Bigr)\Bigr]
!> \f]
!> (\f$n_p\f$ logarithms, an \f$(n_p+2)\f$-term sum and one exponentiation, each accurate to one ulp;
!> \f$\varepsilon = 2u\f$ makes every bound conservative by a factor 2). A point whose exact residual
!> vanishes has a computed residual within \f$\rho\f$, so no stricter criterion can be met reliably,
!> while any looser one accepts a residual the arithmetic can still resolve. The CV_params tolerances
!> (abs_tol, rel_tol and their logarithmic versions) are therefore not used here. A component whose
!> total is exactly zero has no root at finite concentrations (the species carrying it must vanish);
!> for it the bound is widened to `CV_params%zero`\f$\cdot\sum_j|U_{ij}|\f$, i.e. its species are accepted
!> as absent once they are below the numerically-zero concentration.
!>
!> Robustness features:
!> - Best-solution tracking: the iterate with the smallest \f$\max_i |r_i|/\rho_i\f$ is kept (over the
!>   whole call, perturbation retry included) and returned if convergence fails
!> - Stagnation detection: `n_stag_max` consecutive iterations without a `stag_rtol` relative decrease of
!>   the normalised residual \f$\max_i |r_i|/\rho_i\f$, or one iteration in which no step could be accepted
!> - Perturbation retry: on iteration-budget exhaustion or stagnation, restart once from the geometric
!>   mean (log-space midpoint) of the initial guess and the best iterate,
!>   \f$ c_{1,j}^{\text{retry}} = \sqrt{\max(c_{1,j}^{\text{ig}},\,\texttt{tiny}) \cdot \max(c_{1,j}^{\text{best}},\,\texttt{tiny})} \f$,
!>   resetting the iteration and stagnation counters. A second exhaustion/stagnation returns the best
!>   iterate (re-speciated for consistency) with `CV_flag = .false.`: a residual above the round-off
!>   bound is never accepted as converged.
!>
!> \param[in,out] this Aqueous chemistry class object; on exit its primary and secondary variable
!>                     activity concentrations are those of the returned conc_nc
!> \param[in] c1_ig Initial guess for primary species concentrations
!> \param[in] conc_comp Component concentrations
!> \param[out] conc_nc Variable activity species concentrations (already allocated); the converged
!>                     solution, or the best iterate seen when CV_flag is .false.
!> \param[out] niter Number of Newton iterations performed since the last (re)start [-]
!> \param[out] CV_flag Convergence flag: TRUE if converged, FALSE otherwise [-]
!> \warning A diagnostic line is printed to stdout on the iteration-budget failure exit (the stagnation
!>          failure exit is silent, as in Newton_EI_eq_ideal).
!> \sa Newton_EI_eq_ideal

subroutine compute_c_nc_from_u_Newton_ideal(this,c1_ig,conc_comp,conc_nc,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c  !< Import aqueous chemistry class for speciation methods
    implicit none                                       !< Enforce explicit variable declarations
!> Arguments
    class(aqueous_chemistry_c) :: this                 !< Aqueous chemistry class object [-]
    real(kind=8), intent(in) :: c1_ig(:)               !< Initial guess for primary species concentrations (n_p)
    real(kind=8), intent(in) :: conc_comp(:)           !< Component concentrations (total element concentrations) (n_p)
    real(kind=8), intent(out) :: conc_nc(:)            !< Variable activity species concentrations (already allocated) (n_v)
    integer(kind=4), intent(out) :: niter              !< Number of Newton iterations performed [-]
    logical, intent(out) :: CV_flag                    !< Convergence flag: .true. if converged, .false. otherwise [-]
!> Local variables for speciation algebra dimensions, Jacobian and residual
    integer(kind=4) :: i,n_e,n_p,n_v                  !< Loop index; num equilibrium reactions (= n_v - n_p); num primary species; num variable activity species [-]
    real(kind=8), allocatable :: dc2v_dc1(:,:)        !< Jacobian \f$ \partial \mathbf{c}_{2,v}/\partial \mathbf{c}_1 \f$ for ideal solution (n_e x n_p)
    real(kind=8), allocatable :: residual(:)           !< Newton residual \f$ \mathbf{r} = \mathbf{U}_1 \mathbf{c}_1 + \mathbf{U}_2 \mathbf{c}_{2,v} - \mathbf{u} \f$ (n_p)
    real(kind=8), allocatable :: residual_trial(:)     !< Trial residual vector at the line-search candidate point (n_p)
    real(kind=8), allocatable :: log_gamma_var_act(:)  !< Log10 activity coefficients for variable activity species (all zero for aqueous species) (n_v)
!> Log-space Newton variables: \f$ \mathbf{p} = \log_{10}(\mathbf{c}_1) \f$
    real(kind=8), allocatable :: log_c1(:)             !< Log10 primary concentrations \f$ p_j = \log_{10}(c_{1,j}) \f$ [-] (n_p)
    real(kind=8), allocatable :: Delta_p(:)            !< Search direction in log-space: Newton step or steepest descent \f$-\mathbf{J}_p^T\mathbf{r}\f$ [-] (n_p)
    real(kind=8), allocatable :: mat_lin_syst_p(:,:)   !< Log-space Jacobian \f$ J_p(:,j) = \ln(10)\,c_{1,j}\,J(:,j) \f$ [-] (n_p x n_p)
    real(kind=8), allocatable :: mat_lin_syst_eq(:,:)  !< Column-equilibrated copy of the log-space Jacobian (overwritten by dgesv) [-] (n_p x n_p)
    real(kind=8), allocatable :: col_scale(:)          !< Column scaling factors of the equilibration [-] (n_p)
    integer(kind=4), allocatable :: ipiv_dgesv(:)      !< Pivot index array for LAPACK dgesv (n_p)
    integer(kind=4) :: info_dgesv                      !< Return code from LAPACK dgesv (0 = success)
    real(kind=8), parameter :: ln10 = log(10d0)        !< Natural logarithm of 10: \f$ \ln(10) \approx 2.302585 \f$ [-]
    real(kind=8), parameter :: max_dp = 10d0           !< Maximum allowed log-step magnitude per primary species (ten orders of magnitude) [-]
    real(kind=8) :: scale_factor                        !< \f$\|\Delta\mathbf{p}\|_\infty\f$ before bounding; the step is multiplied by max_dp/scale_factor when it exceeds max_dp [-]
!> Search direction and Armijo line search variables
    integer(kind=4) :: i_dir                           !< Direction attempt counter (at most two: Newton, then steepest descent) [-]
    logical :: use_sd                                   !< Current direction is steepest descent \f$-\mathbf{J}_p^T\mathbf{r}\f$
    logical :: sing_flag                                !< .true. if dgesv failed or the Newton step contains NaN/huge values
    real(kind=8) :: alpha                              !< Step length \f$ \alpha \in (0,1] \f$, halved on each backtrack [-]
    integer(kind=4) :: i_LS                            !< Line search backtracking counter [-]
    integer(kind=4), parameter :: max_LS_tries = 10    !< Maximum number of line-search trials per direction (\f$\alpha_{\min} = 2^{-9}\f$) [-]
    real(kind=8), parameter :: c_armijo = 1d-4         !< Armijo sufficient-decrease constant \f$c\f$ [-]
    logical :: step_accepted                           !< Flag: .true. if the line search accepted a step
    real(kind=8) :: phi0                               !< Merit function \f$\phi(0)=\tfrac12\|\mathbf{W}\mathbf{r}\|_2^2\f$ at the current iterate
    real(kind=8) :: dphi0                              !< \f$\phi'(0)=(\mathbf{W}\mathbf{r})^T\mathbf{W}\mathbf{J}_p\Delta\mathbf{p}\f$: directional derivative of \f$\phi\f$ along the direction actually used (must be < 0)
    real(kind=8) :: phi_trial                          !< Merit function \f$\phi(\alpha)\f$ at the trial point
    real(kind=8), allocatable :: w_res(:)              !< Merit weights \f$1/\rho_i\f$ at the current iterate (n_p)
    real(kind=8), allocatable :: w2r(:)                !< \f$\mathbf{W}^2\mathbf{r}\f$: gradient factor of \f$\phi\f$ (n_p)
    real(kind=8), allocatable :: conc_nc_save(:)       !< Saved concentrations for line-search rollback (n_v)
!> Round-off floor convergence criterion
    real(kind=8) :: eps_d                               !< Machine epsilon epsilon(1d0) = 2u [-]
    real(kind=8) :: zero_tol                            !< Numerically-zero concentration (CV_params%zero), for components with exactly zero total
    real(kind=8), allocatable :: res_floor(:)          !< Round-off bound \f$\rho_i\f$ of each component residual at the current iterate (n_p)
    real(kind=8) :: res_ratio                           !< \f$\max_i |r_i|/\rho_i\f$ at the current iterate: converged iff <= 1 [-]
    real(kind=8), allocatable :: abs_comp_mat(:,:)     !< \f$|\mathbf{U}|\f$ (n_p x n_v)
    real(kind=8), allocatable :: abs_U2(:,:)           !< \f$|\mathbf{U}_2|\f$ (n_p x n_e)
    real(kind=8), allocatable :: abs_Se_star(:,:)      !< \f$|\mathbf{S}^*_{e,v,1}|\f$ (n_e x n_p)
    real(kind=8), allocatable :: abs_logK_star(:)      !< \f$|\log_{10}\mathbf{K}^*|\f$ (n_e)
    real(kind=8), allocatable :: abs_log_gamma_1(:)    !< \f$|\log_{10}\gamma|\f$ of the primary species (n_p)
    real(kind=8), allocatable :: abs_log_gamma_2(:)    !< \f$|\log_{10}\gamma|\f$ of the secondary variable activity species (n_e)
    real(kind=8), allocatable :: row_sum_abs_U(:)      !< \f$\sum_j |U_{ij}|\f$ (n_p)
    logical, allocatable :: zero_comp(:)               !< Component total exactly zero (n_p)
!> Stagnation detection and best-solution tracking variables
    real(kind=8), allocatable :: conc_nc_best(:)       !< Best concentrations found across all iterations (n_v)
    real(kind=8) :: best_res_ratio                     !< Smallest \f$\max_i |r_i|/\rho_i\f$ encountered
    real(kind=8) :: res_ratio_prev                     !< Reference normalised residual for the stagnation test (last value that improved by at least stag_rtol)
    integer(kind=4) :: n_stag                          !< Counter of consecutive stagnation iterations; set to n_stag_max at once when no step could be accepted [-]
    integer(kind=4), parameter :: n_stag_max = 20      !< Max allowed stagnation iterations before the perturbation retry / exit [-]
    real(kind=8), parameter :: stag_rtol = 1d-3        !< Relative tolerance: stagnation if the normalised residual does not fall below \f$(1 - \text{stag\_rtol})\f$ times its reference [-]
    logical :: retry_applied                            !< Whether the perturbation retry has been attempted (at most once per call) [-]
    logical :: skip_eval                                !< Skip residual/secondary recomputation after an accepted step (reuse the trial evaluation)
!> Cached loop-invariant parameters (avoid repeated deep struct dereferences per iteration)
    integer(kind=4) :: niter_max                        !< Maximum Newton iterations (cached from CV_params)
    real(kind=8), allocatable :: comp_mat(:,:)          !< Full component matrix (cached, n_p x n_v) for residual computation
    real(kind=8), allocatable :: U1(:,:)                !< Component matrix slice for primary species comp_mat(:,1:n_p) (cached)
    real(kind=8), allocatable :: U2(:,:)                !< Component matrix slice for secondary species comp_mat(:,n_p+1:n_v) (cached)

!> Pre-Process: Initialize variables and allocate arrays
    CV_flag=.false.                                    !< Assume non-convergence until proven otherwise
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions   !< Extract number of equilibrium reactions from speciation algebra
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species   !< Extract number of primary species from speciation algebra
    n_v=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species !< Extract total number of variable activity species
    if (n_v/=n_p+n_e) error stop "compute_c_nc_from_u_Newton_ideal: num_var_act_species /= num_prim_species + num_eq_reactions"
    !> Cache loop-invariant parameters from deep struct dereferences
    niter_max=this%solid_chemistry%reactive_zone%CV_params%niter_max !< Cache max Newton iterations
    zero_tol=this%solid_chemistry%reactive_zone%CV_params%zero       !< Numerically-zero concentration
    eps_d=epsilon(1d0)                                 !< Cache machine epsilon
    !> Cache component matrix slices (avoid repeated deep dereference + slicing each iteration)
    comp_mat=this%solid_chemistry%reactive_zone%speciation_alg%comp_mat !< Cache full component matrix
    U1=comp_mat(:,1:n_p) !< Primary species slice
    U2=comp_mat(:,n_p+1:n_v) !< Secondary species slice
    !> Loop-invariant factors of the round-off bound
    abs_comp_mat=abs(comp_mat)
    abs_U2=abs(U2)
    allocate(row_sum_abs_U(n_p),zero_comp(n_p))
    do i=1,n_p
        row_sum_abs_U(i)=sum(abs_comp_mat(i,:))
        zero_comp(i)=(conc_comp(i)==0d0)
    end do
    log_gamma_var_act=this%get_log_gamma_var_act()     !< Retrieve \f$ \log_{10}(\gamma_i) \f$ for all variable activity species (zero for aqueous species)
    abs_log_gamma_1=abs(log_gamma_var_act(1:n_p))
    abs_log_gamma_2=abs(log_gamma_var_act(n_p+1:n_v))
    if (n_e>0) then
        if (.not. allocated(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star)) &
            error stop "compute_c_nc_from_u_Newton_ideal: Se_nc_1_star not allocated"
        if (.not. allocated(this%logK_star)) error stop "compute_c_nc_from_u_Newton_ideal: logK_star not allocated"
        if (size(this%logK_star)/=n_e) error stop "compute_c_nc_from_u_Newton_ideal: size(logK_star) /= num_eq_reactions"
        abs_Se_star=abs(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star)
        abs_logK_star=abs(this%logK_star)
    else
        allocate(abs_Se_star(0,n_p),abs_logK_star(0))
    end if
    !> Allocate residual vectors (n_p), Jacobian (n_e x n_p) and round-off bound (n_p)
    allocate(residual(n_p),residual_trial(n_p),dc2v_dc1(n_e,n_p),res_floor(n_p))
    !> Allocate log-space Newton arrays: log concentrations, step, Jacobian matrices (n_p x n_p), equilibration and pivots
    allocate(log_c1(n_p),Delta_p(n_p),mat_lin_syst_p(n_p,n_p),mat_lin_syst_eq(n_p,n_p),col_scale(n_p),ipiv_dgesv(n_p))
    !> Allocate line search state-save array (n_v) and best-solution tracker (n_v)
    allocate(conc_nc_save(n_v),conc_nc_best(n_v))

    niter=0                                            !< Reset iteration counter
    n_stag=0                                           !< Reset stagnation counter
    res_ratio_prev=huge(1d0)                           !< Initialize reference normalised residual to machine-maximum so first iterate always improves
    best_res_ratio=huge(1d0)                           !< Initialize best normalised residual to machine-maximum
    retry_applied=.false.                              !< Perturbation retry not yet attempted
    skip_eval=.false.                                  !< No accepted step yet; must compute residual

!> Initialize primary species concentrations \f$ \mathbf{c}_1 \f$ with the user-supplied initial guess
    conc_nc(1:n_p)=c1_ig
    call this%set_conc_prim_species(conc_nc(1:n_p))    !< Object state consistent with conc_nc even if the guess is already converged

!> Process: Newton-Raphson iteration loop with Armijo line search in log-space
    newton_loop: do                                    !< Begin main Newton iteration loop
        niter=niter+1                                  !< Increment Newton iteration counter
        if (niter>niter_max) then                      !< Exceeded iteration budget?
            if (.not. retry_applied) then               !< Has perturbation retry been tried yet?
                call perturbation_retry()               !< Restart once from the geometric mean of the initial guess and the best iterate
                cycle newton_loop                       !< Restart Newton loop with perturbed initial guess
            end if                                      !< End retry guard
            call finish_with_best()                     !< Return the best iterate (CV_flag stays .false.: above the round-off bound)
            print *, "Too many Newton iterations in speciation, best residual =", best_res_ratio, &
                " times its round-off bound"            !< Report failure
            exit newton_loop                            !< Leave Newton loop
        end if                                          !< End max-iteration guard
    !> Compute secondary species and residual (skipped if reusing the accepted line-search evaluation)
        if (skip_eval) then                            !< Reuse residual from accepted line-search step
            skip_eval = .false.                         !< Consume flag; next iteration will recompute unless step accepted again
            residual = residual_trial                   !< Already evaluated at the current conc_nc
        else                                            !< Compute fresh secondary species and residual
            call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p),log_gamma_var_act,conc_nc(n_p+1:n_v))  !< Mass-action law
            residual=-conc_comp                         !< Initialise with \f$-\mathbf{u}\f$
            call dgemv('N',n_p,n_v,1d0,comp_mat,n_p,conc_nc,1,1d0,residual,1)  !< \f$ \mathbf{r} = \mathbf{U}\mathbf{c}_{nc} - \mathbf{u} \f$ via BLAS (no temporary)
        end if                                          !< End skip_eval check
    !> Round-off bound of each component residual at this iterate, the merit weights and the normalised residual
        call compute_res_floor(conc_nc,res_floor)
        w_res=1d0/max(res_floor,tiny(1d0))             !< \f$1/\rho_i\f$ (frozen during this iteration's line search)
        res_ratio=maxval(abs(residual)*w_res)          !< \f$\max_i |r_i|/\rho_i\f$
    !> Track best solution: keep the iterate with the smallest normalised residual
        if (res_ratio<best_res_ratio) then
            best_res_ratio=res_ratio
            conc_nc_best=conc_nc
        end if
    !> Convergence check: every component residual within the round-off bound of its own evaluation
        if (res_ratio<=1d0) then
            CV_flag=.true.                             !< Mark as converged
            exit newton_loop                            !< Leave Newton loop — converged
        end if                                          !< End convergence check
    !> Compute Jacobian \f$ \partial\mathbf{c}_{2,v}/\partial\mathbf{c}_1 \f$ analytically for ideal solution
        call this%compute_dc2v_dc1_ideal(conc_nc(1:n_p),conc_nc(n_p+1:n_v),dc2v_dc1)
    !> Assemble log-space Jacobian: \f$ J_p(:,j) = \ln(10)\,c_{1,j}\,(\mathbf{U}_1 + \mathbf{U}_2\,\partial\mathbf{c}_{2,v}/\partial\mathbf{c}_1)_{:,j} \f$
        mat_lin_syst_p=U1                               !< Start with the primary block
        if (n_e>0) call dgemm('N','N',n_p,n_p,n_e,1d0,U2,n_p,dc2v_dc1,n_e,1d0,mat_lin_syst_p,n_p)  !< \f$\mathbf{U}_1 + \mathbf{U}_2\,\partial\mathbf{c}_{2,v}/\partial\mathbf{c}_1\f$ in place
        do i=1,n_p
            log_c1(i)=log10(max(conc_nc(i),tiny(1d0))) !< Clamp to tiny to avoid log(0)
            mat_lin_syst_p(:,i)=ln10*conc_nc(i)*mat_lin_syst_p(:,i) !< Transform column j to log-space in-place
        end do
    !> Save current concentrations for line search recovery
        conc_nc_save=conc_nc
        step_accepted=.false.                          !< No step accepted yet
    !> Newton direction from the column-equilibrated system (improves conditioning when the
    !> concentrations span many orders of magnitude); dgesv overwrites the copy with the LU factors
        mat_lin_syst_eq=mat_lin_syst_p
        do i=1,n_p
            col_scale(i)=maxval(abs(mat_lin_syst_eq(:,i)))  !< Infinity norm of column
            if (col_scale(i)<tiny(1d0)) col_scale(i)=1d0    !< Guard against zero columns
            mat_lin_syst_eq(:,i)=mat_lin_syst_eq(:,i)/col_scale(i) !< Normalise column
        end do
        Delta_p=-residual                              !< RHS = \f$-\mathbf{r}\f$ (overwritten with the solution on exit)
        call dgesv(n_p,1,mat_lin_syst_eq,n_p,ipiv_dgesv,Delta_p,n_p,info_dgesv)  !< LAPACK LU with partial pivoting
        sing_flag=(info_dgesv/=0)                      !< dgesv returned an error code
        if (.not. sing_flag) then                      !< Only unscale/scan if dgesv succeeded
            Delta_p=Delta_p/col_scale                  !< Recover \f$\Delta\mathbf{p}\f$ of the unscaled system
            sing_flag=any(Delta_p/=Delta_p) .or. maxval(abs(Delta_p))>1d20  !< NaN (x/=x) + overflow check
        end if                                          !< End dgesv-success guard
        use_sd=sing_flag                               !< Unusable Newton step: start directly with steepest descent
        w2r=w_res**2*residual                          !< \f$\mathbf{W}^2\mathbf{r}\f$
        if (use_sd) Delta_p=-matmul(w2r,mat_lin_syst_p)  !< \f$-\mathbf{J}_p^T\mathbf{W}^2\mathbf{r}\f$
        dir_loop: do i_dir=1,2                         !< At most two directions: Newton, then steepest descent
            dphi0=dot_product(w2r,matmul(mat_lin_syst_p,Delta_p))  !< \f$\phi'(0)=(\mathbf{W}\mathbf{r})^T\mathbf{W}\mathbf{J}_p\Delta\mathbf{p}\f$ (descent iff < 0)
            if (dphi0>=0d0) then                       !< Not a descent direction: Armijo cannot be satisfied
                if (use_sd) exit dir_loop              !< \f$\mathbf{J}_p^T\mathbf{W}^2\mathbf{r}=0\f$ to round-off: stationary point of \f$\phi\f$, no step possible
                use_sd=.true.                          !< Switch to steepest descent
                Delta_p=-matmul(w2r,mat_lin_syst_p)    !< \f$-\mathbf{J}_p^T\mathbf{W}^2\mathbf{r}\f$
                cycle dir_loop                         !< Re-check descent with the new direction
            end if                                      !< End descent check
        !> Uniform step bounding: scale \f$\Delta\mathbf{p}\f$ so that \f$\|\Delta\mathbf{p}\|_\infty \le \text{max\_dp}\f$ (preserves direction)
            scale_factor=maxval(abs(Delta_p))
            if (scale_factor>max_dp) Delta_p=Delta_p*(max_dp/scale_factor)
            phi0=5d-1*sum((w_res*residual)**2)         !< \f$\phi(0)\f$
            dphi0=dot_product(w2r,matmul(mat_lin_syst_p,Delta_p))  !< \f$\phi'(0)\f$ of the (possibly rescaled) direction actually used
        !> Armijo backtracking line search in log-space: \f$ \alpha = 1, 0.5, 0.25, \ldots \f$ until
        !> \f$\phi(\alpha)\le\phi(0)+c\,\alpha\,\phi'(0)\f$
            alpha=1d0                                  !< Start with full step
            do i_LS=1,max_LS_tries
            !> Update primary concentrations: \f$ c_{1,j} = 10^{p_j^{\text{old}} + \alpha\,\Delta p_j} \f$ (positive by construction)
                do i=1,n_p
                    conc_nc(i)=10d0**(log_c1(i)+alpha*Delta_p(i))
                end do
            !> Update internal state and recompute secondary species at trial point
                call this%set_conc_prim_species(conc_nc(1:n_p))
                call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p),log_gamma_var_act,conc_nc(n_p+1:n_v))
            !> Evaluate trial residual and merit function at the candidate point
                residual_trial=-conc_comp
                call dgemv('N',n_p,n_v,1d0,comp_mat,n_p,conc_nc,1,1d0,residual_trial,1)  !< \f$\mathbf{U}\mathbf{c}-\mathbf{u}\f$ via BLAS
                phi_trial=5d-1*sum((w_res*residual_trial)**2)  !< \f$\phi(\alpha)\f$ (NaN fails the test below)
                if (phi_trial<=phi0+c_armijo*alpha*dphi0) then  !< Armijo sufficient decrease
                    step_accepted=.true.
                    skip_eval=.true.                   !< Signal next iteration to reuse this evaluation
                    exit
                end if
                alpha=alpha*0.5d0                      !< Halve step length for next backtrack attempt
            end do
            if (step_accepted) exit dir_loop           !< Done for this Newton iteration
            conc_nc=conc_nc_save                       !< Line search failed: restore the pre-step state
            call this%set_conc_prim_species(conc_nc(1:n_p))
            if (use_sd) exit dir_loop                  !< Steepest descent failed too: no step this iteration
            use_sd=.true.                              !< Retry once with steepest descent
            Delta_p=-matmul(w2r,mat_lin_syst_p)        !< \f$-\mathbf{J}_p^T\mathbf{W}^2\mathbf{r}\f$
        end do dir_loop                                 !< End direction loop
    !> Stagnation detection: a rejected step counts as full stagnation at once (the state is unchanged,
    !> so every following iteration would repeat this one exactly); otherwise check whether the
    !> normalised residual fell below \f$(1-\text{stag\_rtol})\f$ times its reference value
        if (.not. step_accepted) then
            n_stag=n_stag_max                          !< Go straight to the retry / exit
        else if (res_ratio<res_ratio_prev*(1d0-stag_rtol)) then
            n_stag=0                                   !< Sufficient improvement: reset stagnation counter
            res_ratio_prev=res_ratio
        else
            n_stag=n_stag+1                            !< No sufficient improvement: increment stagnation counter
        end if
        if (n_stag>=n_stag_max) then                   !< Stagnation limit reached?
            if (.not. retry_applied) then               !< Has perturbation retry been tried yet?
                call perturbation_retry()               !< Restart once from the geometric mean of the initial guess and the best iterate
                cycle newton_loop                       !< Restart Newton loop with perturbed initial guess
            end if                                      !< End stagnation retry guard
            call finish_with_best()                     !< Return the best iterate (CV_flag stays .false.)
            exit newton_loop                            !< Stagnated: leave the Newton loop
        end if                                          !< End stagnation-limit check
    end do newton_loop                                  !< End main Newton iteration loop

contains
    !> Round-off bound \f$\rho_i\f$ of the computed residual of each component at the point c_nc (see the
    !> header): forward-error bound of the \f$(n_v+1)\f$-term inner product forming \f$r_i\f$, plus the
    !> propagated round-off of the mass-action evaluation of the secondary species. For a component
    !> with exactly zero total the bound is widened to zero_tol*sum_j|U_ij| (its species are accepted
    !> as absent once below the numerically-zero concentration).
    subroutine compute_res_floor(c_nc,floor)
        real(kind=8), intent(in) :: c_nc(:)            !< Variable activity species concentrations (n_v)
        real(kind=8), intent(out) :: floor(:)          !< Round-off bound of each component residual (n_p)
        real(kind=8) :: abs_log_c1(n_p)                !< \f$|\log_{10} c_{1,k}|\f$
        real(kind=8) :: rel_err_c2(n_e)                !< Relative round-off bound \f$\delta_j\f$ of each secondary species
        integer(kind=4) :: k
        abs_log_c1=abs(log10(max(c_nc(1:n_p),tiny(1d0))))
        rel_err_c2=eps_d*(1d0+ln10*dble(n_p+3)*(matmul(abs_Se_star,abs_log_c1+abs_log_gamma_1)+abs_logK_star+abs_log_gamma_2))
        floor=dble(n_v+1)*eps_d*(abs(conc_comp)+matmul(abs_comp_mat,abs(c_nc)))+ &
            matmul(abs_U2,abs(c_nc(n_p+1:n_v))*rel_err_c2)
        do k=1,n_p
            if (zero_comp(k)) floor(k)=max(floor(k),zero_tol*row_sum_abs_U(k))
        end do
    end subroutine compute_res_floor

    !> Perturbation retry: restart from the geometric mean (log-space midpoint) of the initial guess
    !> and the best iterate, resetting the iteration and stagnation counters (the best iterate is kept)
    subroutine perturbation_retry()
        integer(kind=4) :: k
        retry_applied=.true.                           !< Mark retry as used (only one attempt allowed per call)
        do k=1,n_p
            conc_nc(k)=sqrt(max(c1_ig(k),tiny(1d0))*max(conc_nc_best(k),tiny(1d0)))  !< Geometric mean of initial guess and best iterate
        end do
        call this%set_conc_prim_species(conc_nc(1:n_p))
        niter=0                                        !< Reset Newton iteration counter for fresh restart
        n_stag=0                                       !< Reset stagnation counter
        res_ratio_prev=huge(1d0)                       !< Reset reference normalised residual to \f$+\infty\f$
        skip_eval=.false.                              !< Force re-evaluation after retry
    end subroutine perturbation_retry

    !> Restore the best iterate and re-speciate for consistency (object state = returned conc_nc)
    subroutine finish_with_best()
        conc_nc=conc_nc_best
        call this%set_conc_prim_species(conc_nc(1:n_p))
        call this%compute_c2v_from_c1_ideal(conc_nc(1:n_p),log_gamma_var_act,conc_nc(n_p+1:n_v))
    end subroutine finish_with_best
end subroutine                                          !< End compute_c_nc_from_u_Newton_ideal