!> Computes Jacobian of Newton residual with respect to primary concentrations in reactive mixing iteration using Euler fully implicit in chemical reactions
!> We assume the chemical system has equilibrium and kinetic reactions
!> We assume all primary species are aqueous
!! We DO NOT apply lumping in this subroutine
subroutine compute_dfk_dc1_EI_ideal(this,c1,c2v,drk_dc,Delta_t,theta,mix_ratio_r,dfk_dc1)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(in) :: c1(:)
    real(kind=8), intent(in) :: c2v(:)
    real(kind=8), intent(in) :: drk_dc(:,:) !> Jacobian of kinetic reaction rates
    !real(kind=8), intent(in) :: porosity !> in solid chemistry associated to this aqueous chemistry
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: mix_ratio_r !> mixing ratio for kinetic reaction rates
    real(kind=8), intent(out) :: dfk_dc1(:,:) !> Jacobian Newton residual - primary concentrations (must be already allocated)
!> Variables
    real(kind=8), allocatable :: dc2v_dc1(:,:) !> Jacobian secondary variable activity concentrations - primary concentrations
    real(kind=8), allocatable :: drk_dc1(:,:) !> Jacobian kinetic reaction rates - primary concentrations
    real(kind=8), allocatable :: drk_dc2v(:,:) !> Jacobian kinetic reaction rates - secondary variable activity concentrations
    real(kind=8), allocatable :: drk_dc_zone(:,:) !> drk_dc with columns permuted to zone-local primary-then-secondary order
    real(kind=8), allocatable :: d_log_gamma_d_I(:) !> derivative of log_10 activity coefficients of variable activity species with respect to ionic activity
    real(kind=8), allocatable :: out_prod(:,:) !> Outer product between d_log_gamma_d_I and squared charges of variable activity species
    integer(kind=4) :: j_perm !> Loop index for permuting drk_dc columns into zone-local order
    integer(kind=4) :: n_p_loc !> Cached number of primary species (zone-local)
    integer(kind=4) :: n_v_loc !> Cached number of variable activity species (zone-local)
!> Pre-process
    n_p_loc = this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_v_loc = this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    allocate(dc2v_dc1(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions, n_p_loc))
!> Process
    !> First we compute d_log_gamma_d_I
        !call this%compute_d_log_gamma_d_I_aq_chem(d_log_gamma_d_I)
    !> Outer product
        !out_prod=outer_prod_vec(d_log_gamma_d_I,this%solid_chemistry%reactive_zone%chem_syst%z2(1:this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
    !> Jacobian of secondary variable activity concentrations
        call this%compute_dc2v_dc1_ideal(c1,c2v,dc2v_dc1)
    !> BUGFIX: drk_dc columns are indexed by chem_syst global species index
    !> (set in set_indices_rk via chem_syst%species(j)). Newton uses zone-local
    !> primary-then-secondary ordering, so we must permute the columns of
    !> drk_dc into zone-local order before splitting into drk_dc1 / drk_dc2v.
    !> Without this permutation the Jacobian carries the wrong sensitivities
    !> whenever zone-local and chem_syst orderings differ (e.g. when a primary
    !> is swapped past a secondary by the equilibrium-elimination algebra), and
    !> Newton stalls with a residual stuck at the noise floor of the wrong
    !> direction.
    allocate(drk_dc_zone(size(drk_dc,1), n_v_loc))
    do j_perm = 1, n_v_loc
        !> Defensive bounds check: ind_var_act_species(j_perm) must index a valid
        !> column of the chem_syst-wide drk_dc. Catches the case where a caller
        !> passes a sliced drk_dc(:,1:n_v) instead of the full chem_syst matrix,
        !> or where a multi-zone setup has zone slot indices exceeding drk_dc cols.
        if (this%solid_chemistry%reactive_zone%ind_var_act_species(j_perm) < 1 .or. &
            this%solid_chemistry%reactive_zone%ind_var_act_species(j_perm) > size(drk_dc,2)) then
            print *, "compute_dfk_dc1_EI_ideal: ind_var_act_species(", j_perm, ") = ", &
                this%solid_chemistry%reactive_zone%ind_var_act_species(j_perm), &
                " out of range [1, ", size(drk_dc,2), "]"
            error stop "compute_dfk_dc1_EI_ideal: ind_var_act_species index out of range for drk_dc"
        end if
        drk_dc_zone(:, j_perm) = drk_dc(:, this%solid_chemistry%reactive_zone%ind_var_act_species(j_perm))
    end do
    !> We separate primary and secondary variable activity species in Jacobian of kinetic reaction rates
        drk_dc1  = drk_dc_zone(:, 1:n_p_loc)
        drk_dc2v = drk_dc_zone(:, n_p_loc+1:n_v_loc)
    !> We compute Jacobian Newton residual
        dfk_dc1=this%solid_chemistry%reactive_zone%speciation_alg%comp_mat(:,1:n_p_loc) + matmul(&
            this%solid_chemistry%reactive_zone%speciation_alg%comp_mat(:,n_p_loc+1:n_v_loc), dc2v_dc1) - &
            theta*Delta_t*mix_ratio_r*matmul(this%solid_chemistry%reactive_zone%U_SkT_prod, &
            drk_dc1 + matmul(drk_dc2v, dc2v_dc1))
!> Post-process
    deallocate(dc2v_dc1, drk_dc1, drk_dc2v, drk_dc_zone)
end subroutine