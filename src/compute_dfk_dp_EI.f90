!> \file compute_dfk_dp_EI.f90
!> \brief Jacobian of the Newton residual of the reactive mixing with respect to the LOG10 of the primary
!>        concentrations, Euler fully implicit in the chemical reactions, NON-ideal solution (routine compute_dfk_dp_EI)
!> \details Same assembly as compute_dfk_dp_EI_ideal (see that file), with the log-log Jacobian of the
!> secondary species \f$ L = \partial\log_{10}\mathbf{c}_{2,v}/\partial\log_{10}\mathbf{c}_1 \f$ computed by
!> compute_dlogc2v_dlogc1, which accounts for the dependence of the activity coefficients on the ionic
!> strength (extended Debye-Hückel), instead of the constant \f$ \mathbf{S}^*_{e,1} \f$ of the ideal case:
!> \f[
!>   \frac{\partial \mathbf{f}_k}{\partial p_j} = \ln 10\,\Bigl[ \mathbf{U}_1(:,j)\,c_{1,j} + \mathbf{U}_2\,(\mathbf{c}_{2,v}\circ L(:,j)) \Bigr]
!>   - \Delta t\,\theta\,\lambda_r\,(\mathbf{U}\mathbf{S}_k^T)\Bigl[ \frac{\partial \mathbf{r}_k}{\partial \log_{10} c_{1,j}}
!>   + \frac{\partial \mathbf{r}_k}{\partial \log_{10} \mathbf{c}_{2,v}}\, L(:,j) \Bigr]
!> \f]
!> from the log-form kinetic Jacobian \f$ \partial\mathbf{r}_k/\partial\log_{10}\mathbf{c} \f$ of compute_rk_Jac_rk_anal
!> (activity coefficients taken as constant there). No 1/c factor is formed anywhere.
!>
!> We assume the chemical system has equilibrium and kinetic reactions and that all primary species are
!> aqueous.
subroutine compute_dfk_dp_EI(this,c1,c2v,drk_dlogc,Delta_t,theta,mix_ratio_r,dfk_dp)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(in) :: c1(:) !> primary species concentrations (n_p)
    real(kind=8), intent(in) :: c2v(:) !> secondary variable activity species concentrations (n_e)
    real(kind=8), intent(in) :: drk_dlogc(:,:) !> Jacobian of the kinetic reaction rates in LOG form, d r_k / d log10 c (n_kin x n_species, chem_syst column order)
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: mix_ratio_r !> mixing ratio for kinetic reaction rates
    real(kind=8), intent(out) :: dfk_dp(:,:) !> Jacobian of the Newton residual with respect to p = log10 c1 (n_p x n_p, must be already allocated)
!> Variables
    real(kind=8), allocatable :: L(:,:) !> d log10 c2v / d log10 c1 with the ionic strength coupling (n_e x n_p)
    real(kind=8), allocatable :: drk_dlogc1(:,:) !> d r_k / d log10 c1 (n_kin x n_p)
    real(kind=8), allocatable :: drk_dlogc2v(:,:) !> d r_k / d log10 c2v (n_kin x n_e)
    real(kind=8), allocatable :: drk_dlogc_zone(:,:) !> drk_dlogc with columns permuted to zone-local primary-then-secondary order
    real(kind=8), allocatable :: dc2v_dp(:,:) !> d c2v / d p = ln10 * (c2v o L(:,j)) (n_e x n_p)
    integer(kind=4) :: j !> Loop index over the primary species (columns)
    integer(kind=4) :: j_perm !> Loop index for permuting drk_dlogc columns into zone-local order
    integer(kind=4) :: n_p_loc !> Cached number of primary species (zone-local)
    integer(kind=4) :: n_v_loc !> Cached number of variable activity species (zone-local)
    real(kind=8), parameter :: ln10 = log(1d1)
!> Pre-process
    n_p_loc = this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_v_loc = this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    allocate(L(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions, n_p_loc))
    allocate(dc2v_dp(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions, n_p_loc))
!> Process
    !> Log-log Jacobian of the secondary variable activity species (ionic strength coupling included)
    call this%compute_dlogc2v_dlogc1(c1,c2v,L)
    !> drk_dlogc columns are indexed by chem_syst global species index (set in set_indices_rk via
    !> chem_syst%species(j)). Newton uses zone-local primary-then-secondary ordering, so the columns are
    !> permuted into zone-local order before splitting into drk_dlogc1 / drk_dlogc2v (see compute_dfk_dp_EI_ideal).
    allocate(drk_dlogc_zone(size(drk_dlogc,1), n_v_loc))
    do j_perm = 1, n_v_loc
        if (this%solid_chemistry%reactive_zone%ind_var_act_species(j_perm) < 1 .or. &
            this%solid_chemistry%reactive_zone%ind_var_act_species(j_perm) > size(drk_dlogc,2)) then
            print *, "compute_dfk_dp_EI: ind_var_act_species(", j_perm, ") = ", &
                this%solid_chemistry%reactive_zone%ind_var_act_species(j_perm), &
                " out of range [1, ", size(drk_dlogc,2), "]"
            error stop "compute_dfk_dp_EI: ind_var_act_species index out of range for drk_dlogc"
        end if
        drk_dlogc_zone(:, j_perm) = drk_dlogc(:, this%solid_chemistry%reactive_zone%ind_var_act_species(j_perm))
    end do
    !> We separate primary and secondary variable activity species in the log-form Jacobian of the kinetic reaction rates
    drk_dlogc1  = drk_dlogc_zone(:, 1:n_p_loc)
    drk_dlogc2v = drk_dlogc_zone(:, n_p_loc+1:n_v_loc)
    !> d c2v / d p_j = ln10 * c2v * L(:,j)
    do j=1,n_p_loc
        dc2v_dp(:,j)=ln10*c2v*L(:,j)
    end do
    !> Log-space Jacobian of the Newton residual
    associate(sa => this%solid_chemistry%reactive_zone%speciation_alg)
    dfk_dp = matmul(sa%comp_mat(:,n_p_loc+1:n_v_loc), dc2v_dp) - &
        theta*Delta_t*mix_ratio_r*matmul(this%solid_chemistry%reactive_zone%U_SkT_prod, &
        drk_dlogc1 + matmul(drk_dlogc2v, L))
    do j=1,n_p_loc
        dfk_dp(:,j) = dfk_dp(:,j) + ln10*c1(j)*sa%comp_mat(:,j)
    end do
    end associate
!> Post-process
    deallocate(L, drk_dlogc1, drk_dlogc2v, drk_dlogc_zone, dc2v_dp)
end subroutine
