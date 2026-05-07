!> \file compute_Re_kin.f90
!> \brief Computes equilibrium reaction amounts when both equilibrium and kinetic reactions are present.
!> \details Same two-option algorithm as compute_Re, but the RHS includes a kinetic correction term
!> \f$ \lambda_r \, S_{k,2v}^T R_k \f$ that accounts for the kinetic reaction amounts at the current step.
!>
!> **Option A (primary-species based, default):**
!> \f[
!>   \min_{R_e} \left\| S_1^T R_e - \frac{c_1 - \hat{c}_1}{\lambda_r} + S_{k,1}^T R_k \right\|^2
!> \f]
!>
!> **Option B (fallback, all variable activity species):**
!> \f[
!>   \min_{R_e} \left\| S_{e,v}^T R_e - \frac{c_v - \hat{c}_v}{\lambda_r} + S_{k,v}^T R_k \right\|^2
!> \f]
!>
!> \param[in,out] this  Aqueous chemistry object [-]
!> \param[in] cv_hat    Concentrations of ALL variable activity species after mixing [M/L³]
!> \param[in] Delta_t   Time step size [T]
!> \param[in] lambda_r  Mixing ratio for reactions [-]
subroutine compute_Re_kin(this,cv_hat,Delta_t,lambda_r)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use vectors_m, only: inf_norm_vec_real
    implicit none
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at time step k+1
    real(kind=8), intent(in) :: cv_hat(:) !> concentrations of ALL variable activity species (primary + secondary) after mixing [M/L³]
    real(kind=8), intent(in) :: Delta_t !> time step [T]
    real(kind=8), intent(in) :: lambda_r !> reaction mixing ratio [-]
!> Local variables
    real(kind=8), allocatable :: Re(:) !> Equilibrium reaction amounts vector [M/L³]
    real(kind=8), allocatable :: c1(:) !> Current primary species concentrations [M/L³]
    real(kind=8), allocatable :: cv(:) !> Current all variable activity species concentrations [M/L³]
    real(kind=8), allocatable :: Rk(:) !> Kinetic reaction amounts at current step [M/L³]
    real(kind=8), allocatable :: Sk_nc(:,:) !> Full kinetic stoichiometric matrix (n_kin x n_nc) [-]
    real(kind=8), allocatable :: Sk1_T(:,:) !> Transpose of primary block of kinetic stoich (n_p x n_kin) [-]
    real(kind=8), allocatable :: Skv_T(:,:) !> Transpose of full var-act block of kinetic stoich (n_v x n_kin) [-]
    real(kind=8), allocatable :: S1_T(:,:) !> Transpose of primary equilibrium stoich block (n_p x n_eq) [-]
    real(kind=8), allocatable :: Se_nc(:,:) !> Full equilibrium stoich matrix (n_eq x n_v) [-]
    real(kind=8), allocatable :: Se_nc_T(:,:) !> Transpose of full equilibrium stoich matrix (n_v x n_eq) [-]
    real(kind=8), allocatable :: rhs(:) !> Right-hand side vector for least-squares [M/L³]
    real(kind=8), allocatable :: rhs_save(:) !> Unweighted RHS copy for residual verification [M/L³]
    real(kind=8), allocatable :: S1_T_save(:,:) !> Unweighted S1_T copy for residual verification [-]
    real(kind=8), allocatable :: weights(:) !> Row weights for weighted least-squares [-]
    real(kind=8), allocatable :: sv(:) !> Singular values from dgelsd [-]
    real(kind=8), allocatable :: work(:) !> LAPACK workspace [-]
    integer(kind=4), allocatable :: iwork(:) !> LAPACK integer workspace [-]
    integer(kind=4) :: n_eq !> Number of equilibrium reactions [-]
    integer(kind=4) :: n_p !> Number of primary species [-]
    integer(kind=4) :: n_v !> Number of variable activity species [-]
    integer(kind=4) :: n_min_cst !> Number of constant activity minerals [-]
    integer(kind=4) :: n_min_var !> Number of variable activity minerals [-]
    integer(kind=4) :: n_min !> Total number of minerals [-]
    integer(kind=4) :: n_gas_cst !> Number of constant activity equilibrium gases [-]
    integer(kind=4) :: n_gas_var !> Number of variable activity gas species [-]
    integer(kind=4) :: n_gas_eq !> Total number of equilibrium gases [-]
    integer(kind=4) :: n_gas_eq_var !> Number of variable activity equilibrium gases [-]
    integer(kind=4) :: n_cst_act !> Total number of constant activity species [-]
    integer(kind=4) :: n_aq_eq !> Number of aqueous equilibrium reactions [-]
    integer(kind=4) :: n_exch !> Number of cation exchange reactions [-]
    integer(kind=4) :: wat_flag !> Water flag [-]
    integer(kind=4) :: n_rhs !> Leading dimension for dgelsd RHS [-]
    integer(kind=4) :: info !> LAPACK return code [-]
    integer(kind=4) :: lwork !> Optimal workspace size [-]
    integer(kind=4) :: liwork !> Integer workspace size [-]
    integer(kind=4) :: rank_out !> Effective rank from dgelsd [-]
    integer(kind=4) :: i !> Loop index [-]
    real(kind=8) :: zero_tol !> Near-zero threshold [M/L³]
    real(kind=8) :: noise_tol !> Effective noise tolerance [M/L³]
    real(kind=8) :: res_norm !> Residual norm after solve [M/L³]
    real(kind=8) :: rhs_norm !> RHS norm for relative residual [M/L³]
    logical :: use_fallback !> Flag to trigger Option B [-]
    real(kind=8), parameter :: w_sec_factor = 1d-1 !> Down-weighting factor for secondary species [-]
    real(kind=8), parameter :: rel_tol = 1d-12 !> Relative tolerance for noise check [-]
    real(kind=8), parameter :: rcond_tol = -1d0 !> RCOND for dgelsd (<0 = machine precision) [-]
    real(kind=8), parameter :: res_rel_tol = 1d-6 !> Relative residual threshold for quality check [-]
!> Input validation
    if (Delta_t <= 0d0) error stop 'compute_Re_kin error: Delta_t <= 0'
    if (lambda_r <= 0d0) error stop 'compute_Re_kin error: lambda_r <= 0'
    if (lambda_r - this%solid_chemistry%reactive_zone%CV_params%abs_tol > 1d0) &
        error stop 'compute_Re_kin error: lambda_r > 1'
!> Pre-process: cache sizes
    n_eq      = this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    n_p       = this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_v       = this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    n_cst_act = this%solid_chemistry%reactive_zone%speciation_alg%num_cst_act_species
    n_min_cst = this%solid_chemistry%reactive_zone%num_minerals_cst_act
    n_min_var = this%solid_chemistry%reactive_zone%num_minerals_var_act
    n_min     = this%solid_chemistry%reactive_zone%num_minerals
    n_gas_cst = this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
    n_gas_var = this%solid_chemistry%reactive_zone%gas_phase%num_var_act_species
    n_gas_eq  = this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq
    n_gas_eq_var = this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
    n_aq_eq   = this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts
    n_exch    = this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
    wat_flag  = this%aq_phase%wat_flag
    zero_tol  = this%solid_chemistry%reactive_zone%CV_params%zero
    if (size(cv_hat) < n_v) error stop 'compute_Re_kin error: size(cv_hat) < n_v'
    allocate(Re(n_eq))
    Re = 0d0
    !> Fetch kinetic correction: Rk and Sk_nc
    Rk     = this%get_Rk()
    Sk_nc  = this%get_Sk_nc()          !> shape (n_kin, n_nc); n_nc = n_cst_act + n_v
!> Process
    if (n_p < n_eq) then
        use_fallback = .true.
    else
        use_fallback = .false.
        !> ---- Option A: Primary-species based weighted least-squares ----
        c1 = this%get_c1()
        n_rhs = max(n_p, n_eq)
        allocate(rhs(n_rhs))
        rhs = 0d0
        !> Base RHS: (c1 - cv_hat_primary) / lambda_r
        rhs(1:n_p) = (c1 - cv_hat(1:n_p)) / lambda_r
        !> Kinetic correction: subtract S_k1^T * Rk from RHS (primary rows of Sk_nc)
        Sk1_T = transpose(Sk_nc(:, 1:n_p))   !> shape (n_p, n_kin)
        rhs(1:n_p) = rhs(1:n_p) - matmul(Sk1_T, Rk)
        deallocate(Sk1_T)
        noise_tol = max(zero_tol, rel_tol * inf_norm_vec_real(c1) / lambda_r)
        if (inf_norm_vec_real(rhs(1:n_p)) < noise_tol) then
            Re = zero_tol
        else
            allocate(rhs_save(n_p))
            rhs_save = rhs(1:n_p)
            allocate(weights(n_p))
            do i = 1, n_p
                weights(i) = 1d0 / max(abs(c1(i)), zero_tol)
            end do
            Se_nc = this%solid_chemistry%reactive_zone%get_Se_nc_react_zone()  !> shape (n_eq, n_v)
            allocate(S1_T(n_p, n_eq))
            S1_T = transpose(Se_nc(:, 1:n_p))
            deallocate(Se_nc)
            allocate(S1_T_save(n_p, n_eq))
            S1_T_save = S1_T
            do i = 1, n_p
                S1_T(i,:) = weights(i) * S1_T(i,:)
                rhs(i)    = weights(i) * rhs(i)
            end do
            deallocate(weights)
            allocate(sv(min(n_p, n_eq)), work(1), iwork(1))
            call dgelsd(n_p, n_eq, 1, S1_T, n_p, rhs, n_rhs, sv, rcond_tol, rank_out, work, -1, iwork, info)
            lwork = int(work(1)); liwork = iwork(1)
            deallocate(work, iwork)
            allocate(work(lwork), iwork(max(1,liwork)))
            call dgelsd(n_p, n_eq, 1, S1_T, n_p, rhs, n_rhs, sv, rcond_tol, rank_out, work, lwork, iwork, info)
            deallocate(work, iwork, sv, S1_T)
            if (info /= 0 .or. rank_out < n_eq) then
                use_fallback = .true.
            else
                Re = rhs(1:n_eq)
                rhs_norm = inf_norm_vec_real(rhs_save)
                if (rhs_norm > zero_tol) then
                    do i = 1, n_p
                        rhs_save(i) = rhs_save(i) - dot_product(S1_T_save(i,:), Re)
                    end do
                    res_norm = inf_norm_vec_real(rhs_save)
                    if (res_norm > res_rel_tol * rhs_norm) use_fallback = .true.
                end if
            end if
            deallocate(S1_T_save, rhs_save)
        end if
        deallocate(c1, rhs)
    end if
    !> ---- Option B: weighted overdetermined least-squares with all var-act species ----
    if (use_fallback) then
        cv = this%get_conc_nc()
        Se_nc_T = transpose(this%solid_chemistry%reactive_zone%get_Se_nc_react_zone())  !> shape (n_v, n_eq)
        allocate(rhs(n_v))
        rhs = (cv - cv_hat(1:n_v)) / lambda_r
        !> Kinetic correction: subtract Sk_v^T * Rk (var-act rows of Sk_nc)
        Skv_T = transpose(Sk_nc(:, n_cst_act+1:))   !> shape (n_v, n_kin)
        rhs = rhs - matmul(Skv_T, Rk)
        deallocate(Skv_T)
        noise_tol = max(zero_tol, rel_tol * inf_norm_vec_real(cv) / lambda_r)
        if (inf_norm_vec_real(rhs) < noise_tol) then
            Re = zero_tol
        else
            allocate(weights(n_v))
            do i = 1, n_p
                weights(i) = 1d0 / max(abs(cv(i)), zero_tol)
            end do
            do i = n_p+1, n_v
                weights(i) = w_sec_factor / max(abs(cv(i)), zero_tol)
            end do
            do i = 1, n_v
                Se_nc_T(i,:) = weights(i) * Se_nc_T(i,:)
                rhs(i)       = weights(i) * rhs(i)
            end do
            deallocate(weights)
            allocate(sv(min(n_v,n_eq)), work(1), iwork(1))
            call dgelsd(n_v, n_eq, 1, Se_nc_T, n_v, rhs, n_v, sv, rcond_tol, rank_out, work, -1, iwork, info)
            lwork = int(work(1)); liwork = iwork(1)
            deallocate(work, iwork)
            allocate(work(lwork), iwork(max(1,liwork)))
            call dgelsd(n_v, n_eq, 1, Se_nc_T, n_v, rhs, n_v, sv, rcond_tol, rank_out, work, lwork, iwork, info)
            deallocate(work, iwork, sv)
            if (info == 0 .and. rank_out >= n_eq) then
                Re = rhs(1:n_eq)
            else
                Re = zero_tol
            end if
        end if
        deallocate(cv, Se_nc_T, rhs)
    end if
    deallocate(Rk, Sk_nc)
!> Post-process: distribute Re to reaction categories (identical index mapping to compute_Re)
    this%solid_chemistry%Re(1:n_min_cst) = Re(1:n_min_cst)
    this%solid_chemistry%Re(n_min_cst+1:n_min_cst+n_min_var) = &
        Re(n_cst_act-wat_flag+n_aq_eq+1:n_eq-n_gas_var-n_exch)
    this%solid_chemistry%Re(n_min+1:n_min+n_exch) = &
        Re(n_eq-n_gas_var-n_exch+1:n_eq-n_gas_var)
    call this%solid_chemistry%set_re_mean(Delta_t)
    if (associated(this%gas_chemistry)) then
        this%gas_chemistry%Re(1:n_gas_cst) = Re(n_min_cst+1:n_cst_act-wat_flag)
        this%gas_chemistry%Re(n_gas_cst+1:n_gas_eq) = Re(n_eq-n_gas_eq_var+1:n_eq)
        call this%gas_chemistry%set_re_mean(Delta_t)
    end if
    this%Re = Re(n_min_cst+n_gas_cst+1:n_min_cst+n_gas_cst+n_aq_eq)
    call this%set_re_mean(Delta_t)
    deallocate(Re)
end subroutine