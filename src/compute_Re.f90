!> \file compute_Re.f90
!> \brief Computes equilibrium reaction amounts from variable activity species concentrations.
!> \details Single routine for both equilibrium-only and equilibrium+kinetic cases (PDF steps
!> 2d.v-vii). The equilibrium contribution to the species mass balance is first stored in the
!> member `contrib_rxn_amnt_eq`:
!> \f[
!>   \rho_{e,v} = c_v^{k+1} - \hat{c}_v - \Delta t\,\lambda_r\,\theta_r\,S_k^T r_K^{k+1}
!> \f]
!> where the theta_r-weighted kinetic term is zero when there are no kinetic reactions and the
!> previous-step (1-theta_r) kinetic contribution is already carried in \f$\hat{c}_v\f$. The
!> equilibrium amounts \f$R_e\f$ are then obtained from \f$\rho_{e,v}/\lambda_r\f$ by a weighted
!> least-squares solve (also with a single equilibrium reaction, where it averages the mass
!> balances of the participating species instead of trusting one of them), and the
!> time-averaged rates are set to \f$\bar{R}_e = R_e/\Delta t\f$. Two tiers are used:
!>
!> **Option A (primary-species based, default):** Least-squares solve using only the
!> primary species difference, which avoids catastrophic cancellation in secondary species:
!> \f[
!>   \min_{R_e} \left\| S_{1,\mathcal{P}}^T \cdot R_e
!>     - \left( \frac{c_1 - \hat{c}_1}{\lambda_r} - \Delta t\,\theta_r\,\left(S_k^T r_K\right)_{1:n_p} \right)_{\mathcal{P}} \right\|^2
!> \f]
!> where \f$S_1 = S_{e,v}(:, 1:n_p)\f$ is the primary-species block of the stoichiometric
!> matrix and \f$\mathcal{P}\f$ the set of primaries with a nonzero column in it. A primary that
!> takes part in no equilibrium reaction (e.g. a conservative tracer) is left out of every
!> least-squares system: its mass balance holds whatever \f$R_e\f$ is, so its \f$\rho_i\f$ is only
!> Newton-residual/round-off noise, which would otherwise enter the noise check, the weight floor
!> and cap, and the residual check. The primary species \f$c_1\f$ are the Newton unknowns and are determined to full
!> double precision, while the mixed primaries \f$\hat{c}_1\f$ come directly from transport.
!> Neither involves the mass action law, so their difference is free from the
!> subtractive cancellation that plagues the secondary species difference.
!> Solved with LAPACK `dgelsd` (SVD-based, rank-revealing).
!>
!> **Option B (fallback):** Overdetermined least-squares using the participating primaries and
!> ALL secondary variable activity species via the full stoichiometric matrix \f$S_{e,v}\f$:
!> \f[
!>   \min_{R_e} \left\| S_{e,v}^T \cdot R_e
!>     - \left( \frac{c_v - \hat{c}_v}{\lambda_r} - \Delta t\,\theta_r\,S_k^T r_K \right) \right\|^2
!> \f]
!> (rows of non-participating primaries removed), solved with LAPACK `dgelsd` (SVD-based,
!> rank-revealing). This fallback triggers when Option A is underdetermined (fewer participating
!> primaries than equilibrium reactions) or fails (LAPACK error, rank deficiency,
!> excessive weighted residual). A noise-level primary signal does not trigger it: \f$R_e\f$ is
!> then set to the placeholder `zero_tol`.
!>
!> **Mass balance check:** unless \f$R_e\f$ is the placeholder, the computed amounts are verified
!> against the species mass balance they must reproduce, \f$\rho_{e,i} = \lambda_r\,(S_{e,v}^T R_e)_i\f$.
!> \f$\rho_{e,v}\f$ carries the Newton component residual \f$f = U\rho_{e,v}\f$ (\f$|f_k| \le\f$ `abs_tol`
!> of the `*_CV_params.dat` file), which lies outside \f$\mathrm{range}(S_{e,v}^T)\f$ and is
!> redistributed over the species by the weighted fit, with cancellation between species shared by
!> several components; so the species residual cannot be compared with `abs_tol` directly. Since the
!> fit is linear, the residual it must leave for the actual \f$f\f$ is computed exactly
!> (\f$r_{pred} = e_f - \lambda_r S_{e,v}^T R_e(e_f)\f$, \f$e_f\f$ the minimum-norm species vector
!> with \f$U e_f = f\f$, \f$R_e(e_f)\f$ the same fit applied to it) and the check is tolerance-free:
!> (1) \f$|f_k| \le\f$ `abs_tol` + roundoff (the speciation converged to its tolerance) and
!> (2) \f$|r_i - r_{pred,i}| \le\f$ roundoff for every species: `mb_eps_factor` ulps of the
!> concentrations forming \f$\rho_i\f$ plus the forward-error bound of the fit,
!> \f$(m+n_{eq})\,\epsilon\,\kappa\,\|R_e\|_\infty\f$ (\f$\kappa\f$ from dgelsd's singular values),
!> through the species' stoichiometric coefficients. A wrong or rcond-truncated \f$R_e\f$ leaves an
!> \f$O(\rho)\f$ residual that \f$r_{pred}\f$ does not contain. The first violation is reported and
!> the run is stopped with `error stop`.
!>
!> Implementation notes: member sizes are cached in local integers; the stoichiometric matrix
!> \f$S_{e,v}\f$ is referenced in place from `reactive_zone%stoich_mat(:,1:n_v)` (no copy); the
!> weighted least-squares matrices are assembled directly in the layout `dgelsd` overwrites; and the
!> `dgelsd` workspace is kept between calls and re-queried only when the problem size changes
!> (the code is serial; these `save` buffers are not thread safe).
!>
!> \param[in,out] this Aqueous chemistry object [-]
!> \param[in] cv_hat Concentrations of ALL variable activity species after mixing [M/L³]
!> \param[in] Delta_t Time step size [T]
!> \param[in] lambda_r Mixing ratio for reactions [-]
!> \param[in] theta_r Reaction time-weighting factor; multiplies ONLY the kinetic reaction term [-]
!> \post Sets `this%contrib_rxn_amnt_eq` and distributes \f$R_e\f$ into the equilibrium reaction-amount
!>       members of the aqueous, solid (minerals, cation exchange) and gas chemistry, updating each
!>       time-averaged rate \f$\bar{R}_e = R_e/\Delta t\f$ via `set_eq_rxn_rates`.

subroutine compute_Re(this,cv_hat,Delta_t,lambda_r,theta_r)
    use aqueous_chemistry_m, only: aqueous_chemistry_c !> Import aqueous chemistry class definition
    implicit none
    class(aqueous_chemistry_c) :: this !> Aqueous chemistry object (polymorphic, passed-object dummy)
    real(kind=8), intent(in) :: cv_hat(:) !> Concentrations of all variable activity species after mixing, size \f$n_v\f$ [M/L³]
    real(kind=8), intent(in) :: Delta_t !> Time step [T]
    real(kind=8), intent(in) :: lambda_r !> Mixing ratio for reactions [-]
    real(kind=8), intent(in) :: theta_r !> Reaction time-weighting factor; multiplies ONLY the kinetic reaction term [-]
!> Local variables
    real(kind=8), allocatable :: Re(:) !> Equilibrium reaction amounts vector, size \f$n_{eq}\f$ [M/L³]
    real(kind=8), allocatable :: c1(:) !> Current primary species concentrations, size \f$n_p\f$ [M/L³]
    real(kind=8), allocatable :: cv(:) !> Current all variable activity species concentrations (primary+secondary), size \f$n_v\f$ [M/L³]
    real(kind=8), allocatable :: A(:,:) !> Weighted least-squares matrix handed to dgelsd (destroyed): \f$W S_1^T\f$ (n_p x n_eq) or \f$W S_{e,v}^T\f$ (n_v x n_eq) [-]
    real(kind=8), allocatable :: rhs(:) !> Weighted least-squares right-hand side (size m); dgelsd returns the solution in rhs(1:n_eq) [M/L³]
    real(kind=8), allocatable :: weights(:) !> Row weights for weighted least-squares, size m [-]
    integer(kind=4), allocatable :: rows(:) !> Species (columns of S_{e,v}) forming the least-squares rows: participating primaries first, then all secondaries [-]
    integer(kind=4) :: n_pp !> Number of primary species taking part in at least one equilibrium reaction [-]
    integer(kind=4) :: n_eq !> Number of equilibrium reactions [-]
    integer(kind=4) :: n_p !> Number of primary species [-]
    integer(kind=4) :: n_v !> Number of variable activity species (primary + secondary) [-]
    integer(kind=4) :: m !> Number of rows of the least-squares system: n_pp (Option A) or n_pp + n_v - n_p (Option B) [-]
    integer(kind=4) :: info !> LAPACK return code: 0 = success, >0 = rank deficiency [-]
    integer(kind=4) :: rank_out !> Effective rank of matrix returned by dgelsd [-]
    integer(kind=4) :: n_min_cst !> Number of constant activity minerals [-]
    integer(kind=4) :: n_min_var !> Number of variable activity minerals [-]
    integer(kind=4) :: n_min !> Total number of minerals (constant + variable activity) [-]
    integer(kind=4) :: n_gas_cst !> Number of constant activity equilibrium gases [-]
    integer(kind=4) :: n_gas_var !> Number of variable activity gas species [-]
    integer(kind=4) :: n_gas_eq !> Total number of equilibrium gases [-]
    integer(kind=4) :: n_gas_eq_var !> Number of variable activity equilibrium gases [-]
    integer(kind=4) :: n_cst_act !> Total number of constant activity species [-]
    integer(kind=4) :: n_aq_eq !> Number of aqueous equilibrium reactions [-]
    integer(kind=4) :: n_exch !> Number of cation exchange reactions [-]
    integer(kind=4) :: wat_flag !> Water flag: 1 if water is included as a species, 0 otherwise [-]
    integer(kind=4) :: i !> Row index (species) [-]
    integer(kind=4) :: j !> Column index (equilibrium reaction) [-]
    real(kind=8) :: zero_tol !> Tolerance below which concentrations are treated as zero [M/L³]
    real(kind=8) :: noise_tol !> Effective tolerance combining absolute and relative checks [M/L³]
    real(kind=8) :: res_norm !> Weighted residual norm after least-squares solve [M/L³]
    real(kind=8) :: rhs_norm !> Weighted RHS norm for the relative residual check [M/L³]
    real(kind=8) :: c_floor !> Floor of the row-weight denominators [M/L³]
    real(kind=8) :: w_cap !> Cap of the row weights (bounds the weight spread) [-]
    real(kind=8) :: mb_res !> Mass balance residual of one species: rho_i - lambda_r*(S_{e,v}^T Re)_i [M/L³]
    real(kind=8) :: mb_tol !> Roundoff allowance of the mass balance check for one species [M/L³]
    real(kind=8) :: f_tol !> Allowance of the Newton component residual check: abs_tol + roundoff of forming U*rho [M/L³]
    real(kind=8) :: kappa_fit !> Condition number sigma_max/sigma_min of the (retained part of the) weighted fit matrix, from dgelsd's singular values [-]
    real(kind=8), allocatable :: f_newton(:) !> Newton component residual U*rho_{e,v} (n_p): the inconsistency the fit must absorb [M/L³]
    real(kind=8), allocatable :: e_f(:) !> Minimum-norm species vector with U*e_f = f_newton (n_v): dgelsd right-hand side / solution buffer [M/L³]
    real(kind=8), allocatable :: U_copy(:,:) !> Copy of the component matrix handed to dgelsd (destroyed) [-]
    real(kind=8), allocatable :: Re_f(:) !> Reaction amounts the accepted weighted fit assigns to e_f (n_eq) [M/L³]
    real(kind=8), allocatable :: r_pred(:) !> Species residual the fit MUST leave for f_newton: e_f - lambda_r*Se^T*Re_f (n_v) [M/L³]
    real(kind=8), allocatable :: c_scale(:) !> Per-species roundoff scale of rho_{e,v}: max(|c_v|,|c_hat_v|) (n_v) [M/L³]
    integer(kind=4) :: k !> Component index [-]
    integer(kind=4) :: rank_f !> dgelsd effective rank in the check solves [-]
    integer(kind=4) :: info_f !> dgelsd return code in the check solves [-]
    logical :: use_fallback !> Flag to trigger Option B (full-system least-squares) fallback [-]
    logical :: re_noise !> True when Re was set to the zero_tol placeholder (no reaction amounts computed) [-]
    real(kind=8), parameter :: w_sec_factor = 1d-1 !> Down-weighting factor for secondary species in Option B [-]
    real(kind=8), parameter :: w_floor_rel = 1d-8 !> Relative floor for row-weight denominators: |c_i| is floored at w_floor_rel*max|c| so the weight spread stays within 1/rcond_tol and near-zero primaries (e.g. the exchange free site) do not swamp the fit [-]
    real(kind=8), parameter :: rel_tol = 1d-12 !> Relative tolerance: differences below rel_tol * |c| are roundoff noise [-]
    real(kind=8), parameter :: rcond_tol = 1d-8 !> RCOND for dgelsd: singular values below rcond_tol*max(sv) are treated as zero, regularizing ill-conditioned (near-collinear) reaction directions [-]
    real(kind=8), parameter :: res_rel_tol = 1d-6 !> Relative residual threshold for solution quality check [-]
    real(kind=8), parameter :: mb_eps_factor = 1d1 !> Roundoff allowance of the mass balance checks, in ulps of the concentrations entering rho_{e,v} (a few ulps: subtractions in different operation order than Newton's residual) [-]
    !> dgelsd workspace kept between calls (serial code: not thread safe); re-queried only when (m, n) changes
    integer(kind=4), save :: ws_m = -1 !> m of the cached workspace query [-]
    integer(kind=4), save :: ws_n = -1 !> n of the cached workspace query [-]
    integer(kind=4), save :: lwork = 0 !> Cached optimal real workspace size for dgelsd [-]
    integer(kind=4), save :: liwork = 0 !> Cached integer workspace size for dgelsd [-]
    real(kind=8), allocatable, save :: work(:) !> LAPACK dgelsd real workspace [-]
    real(kind=8), allocatable, save :: sv(:) !> Singular values from dgelsd, size min(m, n_eq) [-]
    integer(kind=4), allocatable, save :: iwork(:) !> LAPACK dgelsd integer workspace [-]
!> Input validation
    if (Delta_t <= 0d0) error stop 'compute_Re error: Delta_t <= 0'
    if (lambda_r <= 0d0) error stop 'compute_Re error: lambda_r <= 0'
!> Pre-process
    if (lambda_r-this%solid_chemistry%reactive_zone%CV_params%abs_tol>1d0) then !> Check that mixing ratio does not exceed unity (within absolute tolerance)
        write(*,'(A,ES15.7)') 'compute_Re error: lambda_r>1 (physically invalid), value = ', lambda_r !> Report the offending value before stopping
        error stop 'compute_Re error: lambda_r > 1'
    end if
    !> Cache frequently accessed sizes to avoid repeated deep member traversal
    n_eq      = this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions !> Total number of equilibrium reactions in the reactive zone
    n_p       = this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species !> Number of primary (master) species
    n_v       = this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species !> Number of variable activity species (primary + secondary)
    n_cst_act = this%solid_chemistry%reactive_zone%speciation_alg%num_cst_act_species !> Number of constant activity species (minerals at unit activity, fixed-pressure gases)
    n_min_cst = this%solid_chemistry%reactive_zone%num_minerals_cst_act !> Number of pure-phase minerals with constant (unit) activity
    n_min_var = this%solid_chemistry%reactive_zone%num_minerals_var_act !> Number of minerals with variable activity (e.g. solid solutions)
    n_min     = this%solid_chemistry%reactive_zone%num_minerals !> Total mineral count: n_min_cst + n_min_var
    n_gas_cst = this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act !> Number of equilibrium gases at constant activity (fixed partial pressure)
    n_gas_var = this%solid_chemistry%reactive_zone%gas_phase%num_var_act_species !> Number of gas species with variable activity
    n_gas_eq  = this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq !> Total number of equilibrium gas reactions
    n_gas_eq_var = this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act !> Number of equilibrium gases at variable activity
    n_aq_eq   = this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts !> Number of aqueous-phase equilibrium reactions (complexation, dissociation)
    n_exch    = this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats !> Number of cation exchange reactions
    wat_flag  = this%aq_phase%wat_flag !> Water flag: 1 if water is an explicit species in aq_phase, 0 otherwise
    zero_tol  = this%solid_chemistry%reactive_zone%CV_params%zero !> Near-zero threshold from convergence parameters
    if (size(cv_hat) /= n_v) error stop 'compute_Re error: size(cv_hat) /= n_v'
    allocate(Re(n_eq)) !> Allocate equilibrium reaction amounts vector
    Re = 0d0 !> Safe default: zero reaction amounts
    re_noise = .false. !> Set when a noise-level signal makes Re the zero_tol placeholder
    cv = this%get_conc_var_act() !> c_v^{k+1}: used for rho_{e,v} and, if needed, for the Option B weights [M/L³]
    !> PDF step v: contribution of the equilibrium reactions to the species mass balance,
    !> contrib_rxn_amnt_eq (rho_{e,v}) = c_v^{k+1} - c_hat_v - Delta_t*lambda_r*theta_r*Sk_nc^T*r_K^{k+1}.
    !> theta_r weights ONLY the kinetic reactions and the term vanishes when there are none; the
    !> previous-step (1-theta_r) kinetic contribution is already carried in cv_hat.
    if (this%indices_rk%num_cols > 0) then
        this%contrib_rxn_amnt_eq = cv - cv_hat - Delta_t*lambda_r*theta_r*matmul(this%get_rk_new(), this%get_Sk_nc()) !> matmul(vector, matrix) = Sk_nc^T*r_K without forming the transpose
    else
        this%contrib_rxn_amnt_eq = cv - cv_hat
    end if
    associate(Se_nc => this%solid_chemistry%reactive_zone%stoich_mat(:,1:n_v)) !> S_{e,v} (n_eq x n_v) in place: same columns get_Se_nc_react_zone() copies
!> Process
    !> Rows of the least-squares systems: only species with a nonzero stoichiometric column. A primary
    !> in no equilibrium reaction (e.g. a conservative tracer) has a zero row in S_1^T, so it carries no
    !> information on Re while its rho_i (Newton residual / round-off) would still enter the noise
    !> check, the weight floor/cap and the residual check. Secondaries always participate (S_2 is
    !> invertible) and follow the participating primaries in rows(:).
    allocate(rows(n_v))
    n_pp = 0
    do i = 1, n_p
        if (any(Se_nc(:,i) /= 0d0)) then
            n_pp = n_pp + 1
            rows(n_pp) = i
        end if
    end do
    do i = n_p + 1, n_v
        rows(n_pp + i - n_p) = i
    end do
    !> Least squares also for n_eq = 1: averages the weighted primary mass balances instead of reading
    !> Re off one species; a reaction without primaries fails the rank check and goes to Option B.
    !> When n_pp < n_eq, Option A (participating primaries only) is underdetermined:
    !> dgelsd returns the minimum-norm solution, which distributes errors across
    !> all reactions and produces physically meaningless mineral rates.
    !> In that case, skip directly to the weighted full-system solve (Option B)
    !> which adds all secondary variable-activity species and is never underdetermined.
    if (n_pp < n_eq) then
        use_fallback = .true. !> Force Option B for underdetermined primary system
    else
        use_fallback = .false. !> Option A viable: n_pp >= n_eq (determined or overdetermined)
        !> ---- Option A: Primary-species based weighted least-squares ----
        !> Solve min || W*(S_1^T * Re - (c_1 - c_1_hat)/lambda_r) ||^2 over the participating primaries via dgelsd.
        c1 = this%get_c1() !> Retrieve current primary species concentrations after Newton [M/L³]
        m = n_pp
        !> Relative tolerance check: if |c1 - cv_hat| is at roundoff level relative to |c1|,
        !> the difference is pure floating-point noise and cannot be used to recover Re.
        noise_tol = max(zero_tol, rel_tol * maxval(abs(c1(rows(1:m)))) / lambda_r)
        if (maxval(abs(this%contrib_rxn_amnt_eq(rows(1:m)))) / lambda_r < noise_tol) then !> Primary species difference at noise level?
            Re = zero_tol !> Below roundoff: assign negligible placeholder value
            re_noise = .true.
        else
            !> Row weights: w_i = 1/max(|c1_i|,|c1_hat_i|, floor) to equalise relative contribution.
            !> The scale is max(|c_i|,|c_hat_i|) because the roundoff of contrib_i = c_i - c_hat_i is
            !> set by the larger of the two: a species nearly consumed in the step (c_hat_i >> c_i)
            !> would otherwise be over-weighted far beyond its real signal-to-noise.
            !> The floor is RELATIVE to the primary concentration scale (not the absolute
            !> zero_tol) so a near-zero primary — e.g. the cation-exchange free site x-,
            !> whose concentration can be ~1e-16 when the exchanger is nearly saturated —
            !> cannot get a runaway weight (~1/zero_tol) that dwarfs the informative
            !> mg/na/ca rows. With an absolute floor the weight spread exceeds 1/rcond_tol,
            !> so dgelsd (rcond_tol) discards the exchanging-cation rows as "negligible" and
            !> returns Re ~ 0 (in the free-site stoichiometry direction). Capping the spread
            !> at 1/w_floor_rel keeps all rows above the rcond cut.
            allocate(weights(m), A(m, n_eq), rhs(m))
            c_floor = max(zero_tol, w_floor_rel*max(maxval(abs(c1(rows(1:m)))), maxval(abs(cv_hat(rows(1:m))))))
            do i = 1, m
                weights(i) = 1d0 / max(abs(c1(rows(i))), abs(cv_hat(rows(i))), c_floor)
            end do
            !> Cap the weight spread so the weighting alone cannot push the least-squares
            !> conditioning past dgelsd's rcond_tol cut (=1/rcond_tol=1e8), which would drop
            !> informative reaction directions (e.g. the cation-exchange reactions when the
            !> free-site concentration -> 0 gives it a runaway 1/c weight). 1e6 leaves ~100x
            !> margin below the cut; the system is consistent so any full-rank weighting is exact.
            w_cap = minval(weights)/rcond_tol*1d-2
            where (weights > w_cap) weights = w_cap
            !> Weighted system W*S_1^T*Re = W*rho_1/lambda_r assembled directly: row i of S_1^T is Se_nc(:,rows(i))
            do j = 1, n_eq
                do i = 1, m
                    A(i,j) = weights(i)*Se_nc(j,rows(i))
                end do
            end do
            do i = 1, m
                rhs(i) = weights(i)*this%contrib_rxn_amnt_eq(rows(i))/lambda_r !> PDF step vi RHS (primary block) [M/L³]
            end do
            call ls_solve(m, n_eq, A, rhs, rank_out, info) !> SVD-based LS solve (rank-revealing)
            if (info /= 0 .or. rank_out < n_eq) then !> Check dgelsd success and full rank
                use_fallback = .true. !> dgelsd failed or rank-deficient; try Option B
            else
                Re = rhs(1:n_eq) !> Extract solution from first n_eq entries [M/L³]
                !> Residual verification on the WEIGHTED (per-row relative) system, so a species
                !> orders of magnitude below the largest one cannot hide a large relative error:
                !> ||W*(S1_T*Re - rho_1/lambda_r)||_inf / ||W*rho_1/lambda_r||_inf, from the untouched Se_nc and contrib
                rhs_norm = 0d0
                res_norm = 0d0
                do i = 1, m
                    rhs_norm = max(rhs_norm, abs(weights(i)*this%contrib_rxn_amnt_eq(rows(i))/lambda_r))
                    res_norm = max(res_norm, abs(weights(i)*(this%contrib_rxn_amnt_eq(rows(i))/lambda_r - &
                        dot_product(Se_nc(:,rows(i)), Re))))
                end do
                if (rhs_norm > zero_tol .and. res_norm > res_rel_tol * rhs_norm) then
                    use_fallback = .true. !> Poor solution quality; try full system
                end if
            end if
            deallocate(A, rhs) !> weights(1:m) are kept for the mass balance check (rows/weights of the accepted fit)
        end if
        deallocate(c1) !> Deallocate primary species array (no longer needed)
    end if
    !> ---- Option B: weighted overdetermined least-squares with the participating primaries and all secondaries ----
    if (use_fallback) then !> Enter fallback if Option A failed
        m = n_pp + n_v - n_p !> n_v - n_p = n_eq secondaries, so m >= n_eq: never underdetermined
        !> Relative tolerance check for Option B (same rationale as Option A)
        noise_tol = max(zero_tol, rel_tol * maxval(abs(cv(rows(1:m)))) / lambda_r)
        if (maxval(abs(this%contrib_rxn_amnt_eq(rows(1:m)))) / lambda_r < noise_tol) then !> Check if rhs is at noise level
            Re = zero_tol !> Set all reaction amounts to zero
            re_noise = .true.
        else
            !> Compute row weights: primary species get full weight,
            !> secondary species are down-weighted by w_sec_factor.
            !> Scale = max(|c_i|,|c_hat_i|) (roundoff scale of contrib_i, see Option A), floored
            !> RELATIVE to the concentration scale so a near-zero primary such as the exchange
            !> free site x- cannot get a runaway weight that pushes the informative rows below
            !> the dgelsd rcond_tol cut.
            if (allocated(weights)) deallocate(weights) !> Option A's weights (it failed): replaced by Option B's
            allocate(weights(m), A(m, n_eq), rhs(m))
            c_floor = max(zero_tol, w_floor_rel*max(maxval(abs(cv(rows(1:m)))), maxval(abs(cv_hat(rows(1:m))))))
            do i = 1, n_pp
                weights(i) = 1d0 / max(abs(cv(rows(i))), abs(cv_hat(rows(i))), c_floor)
            end do
            do i = n_pp + 1, m
                weights(i) = w_sec_factor / max(abs(cv(rows(i))), abs(cv_hat(rows(i))), c_floor)
            end do
            !> Cap the weight spread (same rationale as Option A)
            w_cap = minval(weights)/rcond_tol*1d-2
            where (weights > w_cap) weights = w_cap
            !> Weighted system W*Se_nc^T*Re = W*rho_v/lambda_r assembled directly
            do j = 1, n_eq
                do i = 1, m
                    A(i,j) = weights(i)*Se_nc(j,rows(i))
                end do
            end do
            do i = 1, m
                rhs(i) = weights(i)*this%contrib_rxn_amnt_eq(rows(i))/lambda_r !> PDF step vi RHS (full system) [M/L³]
            end do
            call ls_solve(m, n_eq, A, rhs, rank_out, info) !> SVD-based LS solve (rank-revealing)
            if (info == 0) then !> Accept dgelsd's (rcond-regularized) minimum-norm least-squares solution
                Re = rhs(1:n_eq) !> Extract solution [M/L³]
            else
                Re = zero_tol !> Last resort: zero only if dgelsd itself failed [M/L³]
            end if
            deallocate(A, rhs) !> weights(1:m) are kept for the mass balance check (rows/weights of the accepted fit)
        end if
    end if
!> Post-process: distribute Re to appropriate reaction categories using cached indices
    !> Constant activity minerals: Re indices [1 : n_min_cst]
    this%solid_chemistry%rxn_amnt_eq(1:n_min_cst) = Re(1:n_min_cst) !> Assign equilibrium amounts for pure-phase minerals (unit activity) [M/L³]
    !> Variable activity minerals: Re indices [n_cst_act - wat_flag + n_aq_eq + 1 : n_eq - n_gas_var - n_exch]
    this%solid_chemistry%rxn_amnt_eq(n_min_cst+1:n_min_cst+n_min_var) = & !> Assign to solid_chemistry%rxn_amnt_eq after constant activity minerals
        Re(n_cst_act-wat_flag+n_aq_eq+1:n_eq-n_gas_var-n_exch) !> Extract from Re skipping constant activity species, water, and aqueous eq reactions [M/L³]
    !> Cation exchange reactions: Re indices [n_eq - n_gas_var - n_exch + 1 : n_eq - n_gas_var]
    this%solid_chemistry%rxn_amnt_eq(n_min+1:n_min+n_exch) = & !> Assign to solid_chemistry%rxn_amnt_eq after all mineral entries
        Re(n_eq-n_gas_var-n_exch+1:n_eq-n_gas_var) !> Extract cation exchange reaction amounts from Re [M/L³]
    !call this%solid_chemistry%set_eq_rxn_rates(Delta_t) !> Update equilibrium reaction rates in solid chemistry: \f$\bar{R}_e = R_e / \Delta t\f$ [M/L³/T]
    !> Gas phase equilibrium reactions (only if gas chemistry is present)
    if (associated(this%gas_chemistry)) then !> Check if gas chemistry object exists (gas phase present in system)
        !> Constant activity gases: Re indices [n_min_cst + 1 : n_cst_act - wat_flag]
        this%gas_chemistry%rxn_amnt_eq(1:n_gas_cst) = Re(n_min_cst+1:n_cst_act-wat_flag) !> Assign equilibrium amounts for fixed-pressure gases [M/L³]
        !> Variable activity gases: Re indices [n_eq - n_gas_eq_var + 1 : n_eq]
        this%gas_chemistry%rxn_amnt_eq(n_gas_cst+1:n_gas_eq) = Re(n_eq-n_gas_eq_var+1:n_eq) !> Assign equilibrium amounts for variable-fugacity gases [M/L³]
        !call this%gas_chemistry%set_eq_rxn_rates(Delta_t) !> Update equilibrium reaction rates in gas chemistry [M/L³/T]
    end if
    !> Aqueous equilibrium reactions: Re indices [n_min_cst + n_gas_cst + 1 : n_min_cst + n_gas_cst + n_aq_eq]
    this%rxn_amnt_eq = Re(n_min_cst+n_gas_cst+1:n_min_cst+n_gas_cst+n_aq_eq) !> Assign aqueous equilibrium reaction amounts (complexation, dissociation) [M/L³]
    call this%set_eq_rxn_rates(Delta_t) !> Update equilibrium reaction rates in aqueous chemistry [M/L³/T]
    !> Mass balance check of the computed reaction amounts (PDF step v read backwards). The species
    !> residual r_i = rho_i - lambda_r*(S_{e,v}^T Re)_i cannot vanish: rho_{e,v} carries the Newton
    !> component residual f = U*rho (|f_k| <= abs_tol), which is not in range(S_{e,v}^T) and is
    !> redistributed over the species by the weighted fit, with cancellation between species shared
    !> by several components (a species can carry several times abs_tol, e.g. hco3- taking
    !> f_hco3 - f_h when the tiny-concentration h+ row is pinned by its large weight). The fit is
    !> linear and, for a full-rank fit, its residual depends on rho only through f, so the residual it
    !> MUST leave for the actual f is computable exactly:  r_pred = e_f - lambda_r*S_{e,v}^T Re(e_f),
    !> with e_f the minimum-norm species vector having U*e_f = f (dgelsd on the underdetermined
    !> system) and Re(e_f) the accepted fit (same rows, weights and rcond regularisation) applied to
    !> e_f. The check is therefore tolerance-free (roundoff only):
    !>  (1) |f_k| <= abs_tol + roundoff for every component: the Newton speciation did converge to
    !>      the tolerance of the *_CV_params.dat file (what the species residual inherits);
    !>  (2) |r_i - r_pred_i| <= roundoff for every variable activity species: Re reproduces rho up
    !>      to that unavoidable inconsistency and nothing else. A wrong or rcond-truncated Re leaves
    !>      an O(rho) residual in a reaction direction that r_pred does not contain.
    !> Checked over ALL n_v species (also when Option A only fitted the primaries). Skipped when Re
    !> is the zero_tol placeholder (no reaction amounts were computed). A violation stops the run.
    if (.not. re_noise) then
        associate(U => this%solid_chemistry%reactive_zone%speciation_alg%comp_mat)
        c_scale = max(abs(cv), abs(cv_hat))
        f_newton = matmul(U, this%contrib_rxn_amnt_eq(1:n_v))
        !> (1) Newton convergence: component residual within abs_tol (+ roundoff of forming U*rho)
        do k = 1, n_p
            f_tol = this%solid_chemistry%reactive_zone%CV_params%abs_tol + &
                mb_eps_factor*epsilon(1d0)*dot_product(abs(U(k,:)), c_scale)
            if (abs(f_newton(k)) > f_tol) then
                write(*,'(A,I0,A,ES12.4,A,ES12.4,A,ES12.4)') &
                    'ERROR: compute_Re: Newton component residual of component ', k, ' = ', f_newton(k), &
                    ' exceeds abs_tol=', this%solid_chemistry%reactive_zone%CV_params%abs_tol, &
                    ' (+ roundoff) =', f_tol
                write(*,'(A,*(ES12.4))') '  U*rho (all components) =', f_newton
                error stop 'compute_Re error: speciation not converged to abs_tol (U*rho > abs_tol)'
            end if
        end do
        !> minimum-norm e_f with U*e_f = f: dgelsd on the underdetermined (n_p x n_v) system; the
        !> buffer holds f on entry (first n_p entries) and e_f (n_v entries) on exit
        allocate(U_copy(n_p, n_v), e_f(n_v))
        U_copy = U
        e_f = 0d0
        e_f(1:n_p) = f_newton
        call ls_solve(n_p, n_v, U_copy, e_f, rank_f, info_f)
        if (info_f /= 0 .or. rank_f < n_p) error stop 'compute_Re error: component matrix rank deficient in the mass balance check'
        !> the accepted fit (rows(1:m), weights(1:m)) applied to e_f
        allocate(A(m, n_eq), rhs(m))
        do j = 1, n_eq
            do i = 1, m
                A(i,j) = weights(i)*Se_nc(j,rows(i))
            end do
        end do
        do i = 1, m
            rhs(i) = weights(i)*e_f(rows(i))/lambda_r
        end do
        call ls_solve(m, n_eq, A, rhs, rank_f, info_f)
        if (info_f /= 0) error stop 'compute_Re error: dgelsd failed in the mass balance check'
        Re_f = rhs(1:n_eq)
        allocate(r_pred(n_v))
        do i = 1, n_v
            r_pred(i) = e_f(i) - lambda_r*dot_product(Se_nc(:,i), Re_f)
        end do
        !> (2) species residual = predicted residual, to roundoff: a few ulps of the concentrations
        !> forming rho_i, plus the forward error of the fit vector Re entering species i through its
        !> stoichiometric coefficients: the standard LAPACK bound p(m,n)*eps*kappa*||Re||, with
        !> kappa = sigma_max/sigma_min of the (retained part of the) weighted fit matrix from dgelsd's
        !> singular values and p ~ m+n_eq. The absolute precision of Re is set by its LARGEST component:
        !> a 1e-11 reaction amount coexisting with a 0.3 one is resolved to ~1e-16 at best.
        kappa_fit = sv(1)/sv(rank_f)
        do i = 1, n_v
            mb_res = this%contrib_rxn_amnt_eq(i) - lambda_r*dot_product(Se_nc(:,i), Re)
            mb_tol = mb_eps_factor*epsilon(1d0)*c_scale(i) + &
                dble(m+n_eq)*epsilon(1d0)*kappa_fit*lambda_r*sum(abs(Se_nc(:,i)))*maxval(abs(Re))
            if (abs(mb_res - r_pred(i)) > mb_tol) then
                write(*,'(A,I0,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,L1)') &
                    'ERROR: compute_Re: mass balance violated for variable activity species ', i, &
                    ': rho_e=', this%contrib_rxn_amnt_eq(i), ' lambda_r*Se^T*Re=', &
                    lambda_r*dot_product(Se_nc(:,i), Re), ' residual=', mb_res, &
                    ' predicted (from the Newton component residual)=', r_pred(i), &
                    ' |difference| > roundoff=', mb_tol, ' use_fallback=', use_fallback
                write(*,'(A,*(ES12.4))') '  U*rho (Newton component residual) =', f_newton
                write(*,'(A,*(ES12.4))') '  rho_e - lambda_r*Se^T*Re (all species) =', &
                    (this%contrib_rxn_amnt_eq(j) - lambda_r*dot_product(Se_nc(:,j), Re), j=1,n_v)
                write(*,'(A,*(ES12.4))') '  predicted residual (all species)      =', r_pred
                write(*,'(A,*(ES12.4))') '  Re =', Re
                write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES12.4)') '  n_p=', n_p, ' n_pp=', n_pp, ' n_eq=', n_eq, ' n_v=', n_v, &
                    ' fit rank=', rank_f, ' kappa=', kappa_fit
                write(*,'(A,*(ES12.4))') '  singular values of the weighted fit matrix =', sv
                write(*,'(A,*(ES12.4))') '  weights =', weights(1:m)
                error stop 'compute_Re error: mass balance violated'
            end if
        end do
        deallocate(U_copy, e_f, A, rhs, Re_f, r_pred, c_scale, f_newton)
        end associate !> U
    end if
    end associate !> Se_nc
    if (allocated(weights)) deallocate(weights)
    deallocate(Re, cv, rows) !> Deallocate local vectors
contains
    !> min ||A x - b||_2 with LAPACK dgelsd (SVD; singular values below rcond_tol*max are zeroed).
    !> A (m x n, lda = m) is destroyed; the solution is returned in b(1:n). size(b) >= max(m, n):
    !> for m < n (underdetermined, minimum-norm solution) b must hold the n-vector solution.
    !> The workspace is kept between calls and re-queried only when (m, n) changes (serial code).
    subroutine ls_solve(m, n, A, b, rank_out, info)
        integer(kind=4), intent(in) :: m !> Number of rows [-]
        integer(kind=4), intent(in) :: n !> Number of columns (unknowns) [-]
        real(kind=8), intent(inout) :: A(:,:) !> Matrix, shape (m, n); overwritten by dgelsd
        real(kind=8), intent(inout) :: b(:) !> Right-hand side, size >= max(m, n); overwritten with the solution in b(1:n)
        integer(kind=4), intent(out) :: rank_out !> Effective rank [-]
        integer(kind=4), intent(out) :: info !> dgelsd return code (0 = success) [-]
        integer(kind=4) :: ldb !> Leading dimension of b: max(m, n)
        real(kind=8) :: wq(1) !> Real workspace-size query result
        integer(kind=4) :: iwq(1) !> Integer workspace-size query result
        ldb = max(m, n)
        if (m /= ws_m .or. n /= ws_n) then !> Problem size changed: re-query the workspace
            if (allocated(sv)) deallocate(sv)
            allocate(sv(min(m, n)))
            call dgelsd(m, n, 1, A, m, b, ldb, sv, rcond_tol, rank_out, wq, -1, iwq, info) !> Workspace query
            if (info /= 0) return
            lwork = int(wq(1))
            liwork = max(1, iwq(1))
            if (allocated(work)) deallocate(work)
            if (allocated(iwork)) deallocate(iwork)
            allocate(work(lwork), iwork(liwork))
            ws_m = m
            ws_n = n
        end if
        call dgelsd(m, n, 1, A, m, b, ldb, sv, rcond_tol, rank_out, work, lwork, iwork, info) !> SVD-based LS solve
    end subroutine ls_solve
end subroutine compute_Re