!> \file compute_Re.f90
!> \brief Computes equilibrium reaction amounts from variable activity species concentrations.
!> \details This subroutine determines the equilibrium reaction amounts using a two-tier approach:
!>
!> **Option A (primary-species based, default):** Least-squares solve using only the
!> primary species difference, which avoids catastrophic cancellation in secondary species:
!> \f[
!>   \min_{R_e} \left\| S_1^T \cdot R_e - \frac{c_1 - \hat{c}_1}{\lambda_r} \right\|^2
!> \f]
!> where \f$S_1 = S_{e,v}(:, 1:n_p)\f$ is the primary-species block of the stoichiometric
!> matrix. The primary species \f$c_1\f$ are the Newton unknowns and are determined to full
!> double precision, while the mixed primaries \f$\hat{c}_1\f$ come directly from transport.
!> Neither involves the mass action law, so their difference is free from the
!> subtractive cancellation that plagues the secondary species difference.
!> Solved with LAPACK `dgelsd` (SVD-based, rank-revealing).
!>
!> **Option B (fallback):** Overdetermined least-squares using ALL variable activity species
!> (primary + secondary) via the full stoichiometric matrix \f$S_{e,v}\f$:
!> \f[
!>   \min_{R_e} \left\| S_{e,v}^T \cdot R_e - \frac{c_v - \hat{c}_v}{\lambda_r} \right\|^2
!> \f]
!> solved with LAPACK `dgelsd` (SVD-based, rank-revealing). This fallback triggers when
!> Option A fails (LAPACK error, rank deficiency, excessive residual, or noise-level primary signal).
!>
!> Frequently accessed member sizes are cached in local integers to avoid repeated
!> deep pointer chain traversal (e.g. `this%solid_chemistry%reactive_zone%...`).
!>
!> \param[in,out] this Aqueous chemistry object [-]
!> \param[in] cv_hat Concentrations of ALL variable activity species after mixing [M/L³]
!> \param[in] Delta_t Time step size [T]
!> \param[in] lambda_r Mixing ratio for reactions [-]

subroutine compute_Re(this,cv_hat,Delta_t,lambda_r)
    use aqueous_chemistry_m, only: aqueous_chemistry_c !> Import aqueous chemistry class definition
    use vectors_m, only: inf_norm_vec_real !> Import infinity-norm utility for real vectors
    implicit none
    class(aqueous_chemistry_c) :: this !> Aqueous chemistry object (polymorphic, passed-object dummy)
    real(kind=8), intent(in) :: cv_hat(:) !> Concentrations of all variable activity species after mixing, size \f$n_v\f$ [M/L³]
    real(kind=8), intent(in) :: Delta_t !> Time step [T]
    real(kind=8), intent(in) :: lambda_r !> Mixing ratio for reactions [-]
!> Local variables
    real(kind=8), allocatable :: Re(:) !> Equilibrium reaction amounts vector, size \f$n_{eq}\f$ [M/L³]
    real(kind=8), allocatable :: c1(:) !> Current primary species concentrations, size \f$n_p\f$ [M/L³]
    real(kind=8), allocatable :: cv(:) !> Current all variable activity species concentrations (primary+secondary), size \f$n_v\f$ [M/L³]
    real(kind=8), allocatable :: S1_T(:,:) !> Transpose of primary stoichiometric submatrix, shape \f$(n_p, n_{eq})\f$ [-]
    real(kind=8), allocatable :: Se_nc(:,:) !> Full stoichiometric matrix for variable activity species, shape \f$(n_{eq}, n_v)\f$ [-]
    real(kind=8), allocatable :: Se_nc_T(:,:) !> Transpose of stoichiometric matrix for variable activity species, shape \f$(n_v, n_{eq})\f$ [-]
    real(kind=8), allocatable :: rhs(:) !> Right-hand side vector for least-squares, size \f$\max(n_p, n_{eq})\f$ or \f$n_v\f$ [M/L³]
    integer(kind=4) :: n_eq !> Number of equilibrium reactions [-]
    integer(kind=4) :: n_p !> Number of primary species [-]
    integer(kind=4) :: n_v !> Number of variable activity species (primary + secondary) [-]
    integer(kind=4) :: info !> LAPACK return code: 0 = success, >0 = rank deficiency [-]
    integer(kind=4) :: lwork !> Optimal workspace size for dgels [-]
    integer(kind=4) :: n_rhs !> Leading dimension for dgels RHS = max(n_p, n_eq) [-]
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
    real(kind=8), allocatable :: work(:) !> LAPACK dgels workspace array [-]
    real(kind=8) :: zero_tol !> Tolerance below which concentrations are treated as zero [M/L³]
    logical :: use_fallback !> Flag to trigger Option B (full-system least-squares) fallback [-]
    real(kind=8), allocatable :: weights(:) !> Row weights for weighted least-squares, size \f$n_p\f$ or \f$n_v\f$ [-]
    integer(kind=4) :: i !> Loop index for applying row weights [-]
    real(kind=8), parameter :: w_sec_factor = 1d-1 !> Down-weighting factor for secondary species in Option B [-]
    real(kind=8), parameter :: rel_tol = 1d-12 !> Relative tolerance: differences below rel_tol * |c| are roundoff noise [-]
    real(kind=8) :: noise_tol !> Effective tolerance combining absolute and relative checks [M/L³]
    real(kind=8), allocatable :: sv(:) !> Singular values from dgelsd, size min(M,N) [-]
    real(kind=8), allocatable :: rhs_save(:) !> Unweighted RHS copy for residual verification [M/L³]
    real(kind=8), allocatable :: S1_T_save(:,:) !> Unweighted S1_T copy for residual verification [-]
    real(kind=8) :: res_norm !> Residual norm after least-squares solve [M/L³]
    real(kind=8) :: rhs_norm !> Norm of unweighted RHS for relative residual [M/L³]
    integer(kind=4), allocatable :: iwork(:) !> Integer workspace for dgelsd [-]
    integer(kind=4) :: rank_out !> Effective rank of matrix returned by dgelsd [-]
    integer(kind=4) :: liwork !> Integer workspace size for dgelsd [-]
    real(kind=8), parameter :: rcond_tol = -1d0 !> RCOND for dgelsd: <0 means use machine precision [-]
    real(kind=8), parameter :: res_rel_tol = 1d-6 !> Relative residual threshold for solution quality check [-]
!> Input validation
    if (Delta_t <= 0d0) error stop 'compute_Re error: Delta_t <= 0'
    if (lambda_r <= 0d0) error stop 'compute_Re error: lambda_r <= 0'
!> Pre-process
    if (lambda_r-this%solid_chemistry%reactive_zone%CV_params%abs_tol>1d0) then !> Check that mixing ratio does not exceed unity (within absolute tolerance)
        error stop 'compute_Re error: lambda_r>1' !> Abort if lambda_r > 1 (physically invalid)
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
    if (size(cv_hat) < n_v) error stop 'compute_Re error: size(cv_hat) < n_v'
    allocate(Re(n_eq)) !> Allocate equilibrium reaction amounts vector
    Re = 0d0 !> Safe default: zero reaction amounts
!> Process
    !> When n_p < n_eq, Option A (primary-species only) is underdetermined:
    !> dgels returns the minimum-norm solution, which distributes errors across
    !> all reactions and produces physically meaningless mineral rates.
    !> In that case, skip directly to the weighted full-system solve (Option B)
    !> which uses all n_v variable-activity species and is overdetermined.
    if (n_p < n_eq) then
        use_fallback = .true. !> Force Option B for underdetermined primary system
    else
        use_fallback = .false. !> Option A viable: n_p >= n_eq (determined or overdetermined)
        !> ---- Option A: Primary-species based weighted least-squares ----
        !> Solve min || W*(S_1^T * Re - (c_1 - c_1_hat)/lambda_r) ||^2 via dgelsd.
        c1 = this%get_c1() !> Retrieve current primary species concentrations after Newton [M/L³]
        n_rhs = max(n_p, n_eq) !> Leading dimension for dgelsd RHS (must be >= max(M,N))
        allocate(rhs(n_rhs)) !> Allocate RHS vector for dgelsd
        rhs = 0d0 !> Zero-initialise (important when n_eq > n_p)
        rhs(1:n_p) = (c1 - cv_hat(1:n_p)) / lambda_r !> Primary species difference scaled by mixing ratio [M/L³]
        !> Relative tolerance check: if |c1 - cv_hat| is at roundoff level relative to |c1|,
        !> the difference is pure floating-point noise and cannot be used to recover Re.
        noise_tol = max(zero_tol, rel_tol * inf_norm_vec_real(c1) / lambda_r)
        if (inf_norm_vec_real(rhs(1:n_p)) < noise_tol) then !> Check if primary species difference is at noise level
            Re = zero_tol !> Below roundoff: assign negligible placeholder value
        else
            !> Save unweighted data for residual verification
            allocate(rhs_save(n_p))
            rhs_save = rhs(1:n_p)
            !> Row weights: w_i = 1/max(|c1_i|, zero_tol) to equalise relative contribution
            allocate(weights(n_p))
            do i = 1, n_p
                weights(i) = 1d0 / max(abs(c1(i)), zero_tol)
            end do
            !> Extract S_1^T from stoichiometric matrix
            Se_nc = this%solid_chemistry%reactive_zone%get_Se_nc_react_zone() !> Full Se_nc, shape (n_eq, n_v)
            allocate(S1_T(n_p, n_eq)) !> Allocate transpose of primary stoichiometric submatrix
            S1_T = transpose(Se_nc(:, 1:n_p)) !> S_1^T, shape (n_p, n_eq)
            deallocate(Se_nc) !> Deallocate full Se_nc (only S1_T needed)
            !> Save unweighted copy for residual check
            allocate(S1_T_save(n_p, n_eq))
            S1_T_save = S1_T
            !> Apply row weights
            do i = 1, n_p
                S1_T(i,:) = weights(i) * S1_T(i,:)
                rhs(i) = weights(i) * rhs(i)
            end do
            deallocate(weights) !> Deallocate weight vector
            !> SVD-based least-squares solve (dgelsd): handles rank deficiency gracefully
            allocate(sv(min(n_p, n_eq)))
            allocate(work(1), iwork(1))
            call dgelsd(n_p, n_eq, 1, S1_T, n_p, rhs, n_rhs, sv, rcond_tol, rank_out, work, -1, iwork, info) !> Workspace query
            lwork = int(work(1)) !> Extract optimal workspace size
            liwork = iwork(1) !> Extract required integer workspace size
            deallocate(work, iwork) !> Deallocate query workspace
            allocate(work(lwork), iwork(max(1, liwork))) !> Allocate optimal workspace
            call dgelsd(n_p, n_eq, 1, S1_T, n_p, rhs, n_rhs, sv, rcond_tol, rank_out, work, lwork, iwork, info) !> SVD-based LS solve
            deallocate(work, iwork, sv, S1_T) !> Deallocate dgelsd workspace and matrix
            if (info /= 0 .or. rank_out < n_eq) then !> Check dgelsd success and full rank
                use_fallback = .true. !> dgelsd failed or rank-deficient; try Option B
            else
                Re = rhs(1:n_eq) !> Extract solution from first n_eq entries [M/L³]
                !> Residual verification: ||S1_T * Re - rhs||_inf / ||rhs||_inf
                rhs_norm = inf_norm_vec_real(rhs_save)
                if (rhs_norm > zero_tol) then
                    do i = 1, n_p
                        rhs_save(i) = rhs_save(i) - dot_product(S1_T_save(i,:), Re)
                    end do
                    res_norm = inf_norm_vec_real(rhs_save)
                    if (res_norm > res_rel_tol * rhs_norm) then
                        use_fallback = .true. !> Poor solution quality; try full system
                    end if
                end if
            end if
            deallocate(S1_T_save, rhs_save) !> Deallocate residual check arrays
        end if
        deallocate(c1) !> Deallocate primary species array (no longer needed)
        deallocate(rhs) !> Deallocate primary-species RHS
    end if
    !> ---- Option B: weighted overdetermined least-squares with all var act species ----
    if (use_fallback) then !> Enter fallback if Option A failed
        cv = this%get_conc_nc() !> Fetch all variable activity concentrations [M/L³]
        Se_nc_T = transpose(this%solid_chemistry%reactive_zone%get_Se_nc_react_zone()) !> Shape (n_v, n_eq)
        allocate(rhs(n_v)) !> Allocate rhs vector of size n_v
        rhs = (cv - cv_hat) / lambda_r !> Compute (c_v - c_v_hat) / lambda_r [M/L³]
        !> Relative tolerance check for Option B (same rationale as Option A)
        noise_tol = max(zero_tol, rel_tol * inf_norm_vec_real(cv) / lambda_r)
        if (inf_norm_vec_real(rhs) < noise_tol) then !> Check if rhs is at noise level
            Re = zero_tol !> Set all reaction amounts to zero
        else
            !> Compute row weights: primary species get full weight,
            !> secondary species are down-weighted by w_sec_factor
            allocate(weights(n_v))
            do i = 1, n_p
                weights(i) = 1d0 / max(abs(cv(i)), zero_tol)
            end do
            do i = n_p + 1, n_v
                weights(i) = w_sec_factor / max(abs(cv(i)), zero_tol)
            end do
            !> Apply row weights: W*Se_nc_T and W*rhs
            do i = 1, n_v
                Se_nc_T(i,:) = weights(i) * Se_nc_T(i,:)
                rhs(i) = weights(i) * rhs(i)
            end do
            deallocate(weights) !> Deallocate weight vector
            !> SVD-based least-squares solve (dgelsd): handles rank deficiency gracefully
            allocate(sv(min(n_v, n_eq)))
            allocate(work(1), iwork(1))
            call dgelsd(n_v, n_eq, 1, Se_nc_T, n_v, rhs, n_v, sv, rcond_tol, rank_out, work, -1, iwork, info) !> Workspace query
            lwork = int(work(1)) !> Extract optimal workspace size
            liwork = iwork(1) !> Extract required integer workspace size
            deallocate(work, iwork) !> Deallocate query workspace
            allocate(work(lwork), iwork(max(1, liwork))) !> Allocate optimal workspace
            call dgelsd(n_v, n_eq, 1, Se_nc_T, n_v, rhs, n_v, sv, rcond_tol, rank_out, work, lwork, iwork, info) !> SVD-based LS solve
            deallocate(work, iwork, sv) !> Deallocate dgelsd workspace
            if (info == 0 .and. rank_out >= n_eq) then !> Check dgelsd success and full rank
                Re = rhs(1:n_eq) !> Extract solution [M/L³]
            else
                Re = zero_tol !> Last resort: zero if dgelsd failed or rank-deficient [M/L³]
            end if
        end if
        deallocate(cv, Se_nc_T, rhs) !> Clean up all fallback arrays
    end if
!> Post-process: distribute Re to appropriate reaction categories using cached indices
    !> Constant activity minerals: Re indices [1 : n_min_cst]
    this%solid_chemistry%Re(1:n_min_cst) = Re(1:n_min_cst) !> Assign equilibrium amounts for pure-phase minerals (unit activity) [M/L³]
    !> Variable activity minerals: Re indices [n_cst_act - wat_flag + n_aq_eq + 1 : n_eq - n_gas_var - n_exch]
    this%solid_chemistry%Re(n_min_cst+1:n_min_cst+n_min_var) = & !> Assign to solid_chemistry%Re after constant activity minerals
        Re(n_cst_act-wat_flag+n_aq_eq+1:n_eq-n_gas_var-n_exch) !> Extract from Re skipping constant activity species, water, and aqueous eq reactions [M/L³]
    !> Cation exchange reactions: Re indices [n_eq - n_gas_var - n_exch + 1 : n_eq - n_gas_var]
    this%solid_chemistry%Re(n_min+1:n_min+n_exch) = & !> Assign to solid_chemistry%Re after all mineral entries
        Re(n_eq-n_gas_var-n_exch+1:n_eq-n_gas_var) !> Extract cation exchange reaction amounts from Re [M/L³]
    call this%solid_chemistry%set_re_mean(Delta_t) !> Update time-averaged equilibrium reaction rates in solid chemistry: \f$\bar{R}_e = R_e / \Delta t\f$ [M/L³/T]
    !> Gas phase equilibrium reactions (only if gas chemistry is present)
    if (associated(this%gas_chemistry)) then !> Check if gas chemistry object exists (gas phase present in system)
        !> Constant activity gases: Re indices [n_min_cst + 1 : n_cst_act - wat_flag]
        this%gas_chemistry%Re(1:n_gas_cst) = Re(n_min_cst+1:n_cst_act-wat_flag) !> Assign equilibrium amounts for fixed-pressure gases [M/L³]
        !> Variable activity gases: Re indices [n_eq - n_gas_eq_var + 1 : n_eq]
        this%gas_chemistry%Re(n_gas_cst+1:n_gas_eq) = Re(n_eq-n_gas_eq_var+1:n_eq) !> Assign equilibrium amounts for variable-fugacity gases [M/L³]
        call this%gas_chemistry%set_re_mean(Delta_t) !> Update time-averaged equilibrium reaction rates in gas chemistry [M/L³/T]
    end if
    !> Aqueous equilibrium reactions: Re indices [n_min_cst + n_gas_cst + 1 : n_min_cst + n_gas_cst + n_aq_eq]
    this%Re = Re(n_min_cst+n_gas_cst+1:n_min_cst+n_gas_cst+n_aq_eq) !> Assign aqueous equilibrium reaction amounts (complexation, dissociation) [M/L³]
    call this%set_re_mean(Delta_t) !> Update time-averaged equilibrium reaction rates in aqueous chemistry [M/L³/T]
    deallocate(Re) !> Deallocate local Re vector
end subroutine