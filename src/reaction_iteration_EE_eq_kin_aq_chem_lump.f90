!> Computes chemical part of aqueous component concentrations for a reactive mixing Euler explicit iteration assuming there are equilibrium and kinetic reactions
subroutine reaction_iteration_EE_eq_kin_lump(this,Delta_t,conc_comp_react)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(out) :: conc_comp_react(:) !> reaction part of component concentrations (must be already allocated)
!> Variables
    real(kind=8), allocatable :: rk(:) !> new kinetic reaction rates
!> Process
    allocate(rk(this%indices_rk%num_cols)) !> we initialise kinetic reaction rates
    call this%compute_rk_new(rk) !> we compute new kinetic reaction rates
    conc_comp_react=Delta_t*matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,rk) !> we compute reaction part of new component concentrations
    deallocate(rk) !> we deallocate new kinetic reaction rates
end subroutine