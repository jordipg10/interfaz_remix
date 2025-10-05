!> Computes Jacobian of Newton residual with respect to aqueous species concentrations in reactive mixing iteration using Euler fully implicit in chemical reactions
!> We assume the chemical system has only kinetic reactions
!> We assume all species are aqueous
subroutine compute_dfk_dc_aq_EfI(this,drk_dc,porosity,Delta_t,dfk_dc)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use matrices_m, only: id_matrix
    implicit none
!> Arguments
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(in) :: drk_dc(:,:) !> Jacobian of kinetic reaction rates
    real(kind=8), intent(in) :: porosity !> in solid chemistry associated to this aqueous chemistry
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(out) :: dfk_dc(:,:) !> Jacobian Newton residual - aqueous concentrations (must be already allocated)
!> Process
    !> We compute Jacobian Newton residual
        dfk_dc=id_matrix(this%aq_phase%num_species)-(Delta_t/porosity)*matmul(transpose(&
        this%solid_chemistry%reactive_zone%chem_syst%Sk),drk_dc)
end subroutine