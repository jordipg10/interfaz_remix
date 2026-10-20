!> Computes Jacobian of Newton residual with respect to primary concentrations in reactive mixing iteration using Euler fully implicit in chemical reactions
!> We assume the chemical system has equilibrium and kinetic reactions
!> We assume all primary species are aqueous
!! We DO NOT apply lumping in this subroutine
!! NOTE: legacy arithmetic form (dfk/dc), not called by the current Newton solvers, which assemble their Jacobian in
!! log-space from the log-form kinetic Jacobian of compute_rk_Jac_rk_anal (compute_dfk_dp_EI(_ideal)). The drk_dc
!! argument here is still the ARITHMETIC Jacobian d r_k/d c; do not feed it the drk_dlogc of compute_rk_Jac_rk_anal.
subroutine compute_dfk_dc_EI_ideal(this,drk_dc,Delta_t,theta,mix_ratio_r,dfk_dc)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use arrays_m, only: id_matrix
    implicit none
!> Arguments
    class(aqueous_chemistry_c), intent(in) :: this !> aqueous chemistry object
    !real(kind=8), intent(in) :: c1(:)
    !real(kind=8), intent(in) :: c2v(:)
    !real(kind=8), intent(in) :: c_v(:) !> Concentrations of variable activity species [C]
    real(kind=8), intent(in) :: drk_dc(:,:) !> Jacobian of kinetic reaction rates
    !real(kind=8), intent(in) :: porosity !> in solid chemistry associated to this aqueous chemistry
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: mix_ratio_r !> mixing ratio for kinetic reaction rates
    real(kind=8), intent(out) :: dfk_dc(:,:) !> Jacobian Newton residual - primary concentrations (must be already allocated)
!> Variables
    real(kind=8), allocatable :: dc2v_dc1(:,:) !> Jacobian secondary variable activity concentrations - primary concentrations
    real(kind=8), allocatable :: drk_dc1(:,:) !> Jacobian kinetic reaction rates - primary concentrations
    real(kind=8), allocatable :: drk_dc2v(:,:) !> Jacobian kinetic reaction rates - secondary variable activity concentrations
    real(kind=8), allocatable :: d_log_gamma_d_I(:) !> derivative of log_10 activity coefficients of variable activity species with respect to ionic activity
    real(kind=8), allocatable :: out_prod(:,:) !> Outer product between d_log_gamma_d_I and squared charges of variable activity species
!> Process
    !> We compute Jacobian of Newton residual
        dfk_dc=id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)-Delta_t*mix_ratio_r*&
            theta*matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,drk_dc)
!> Post-process
end subroutine