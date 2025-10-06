!> Computes chemical part of aqueous component concentrations for a reactive mixing Euler explicit iteration assuming there are equilibrium and kinetic reactions
subroutine reaction_iteration_EE_eq_kin(this,Delta_t,rk_tilde,conc_comp_react)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    !real(kind=8), intent(in) :: porosity !> porosity
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rate contributions after mixing
    real(kind=8), intent(out) :: conc_comp_react(:) !> reaction part of component concentrations (must be already allocated)
!> Variables
    real(kind=8), allocatable :: rk(:) !> kinetic reaction rates
!> Process
    conc_comp_react=Delta_t*matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,rk_tilde) !> we compute reaction part of new component concentrations
    !> We compute new kinetic reaction rates
    allocate(rk(this%indices_rk%num_cols)) !> we allocate new kinetic reaction rates
    call this%compute_rk(rk) !> we compute new kinetic reaction rates
    !call this%update_rk_old() !> we update old aqueous kinetic reaction rates
    !call this%solid_chemistry%update_rk_old() !> we update old solid kinetic reaction rates
!> Post-process    
    deallocate(rk) !> we deallocate new kinetic reaction rates
end subroutine