!> Computes chemical part of aqueous component concentrations for a reactive mixing Euler explicit iteration assuming there are equilibrium and kinetic reactions
subroutine compute_react_term_EE_eq_kin(this,Delta_t,lambda_r,conc_comp_react)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    !real(kind=8), intent(in) :: R_tilde(:) !> time step
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: lambda_r !> kinetic reaction rate contributions after mixing
    real(kind=8), intent(out) :: conc_comp_react(:) !> reaction part of component concentrations (must be already allocated)
!> Variables
    real(kind=8), allocatable :: rk_old(:) !> kinetic reaction rates
!> Process
    rk_old=this%get_rk_old() !> we get old kinetic reaction rates
    conc_comp_react=Delta_t*lambda_r*matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,rk_old) !> we compute reaction part of new component concentrations
    !> We compute new kinetic reaction rates
    !allocate(rk(this%indices_rk%num_cols)) !> we allocate new kinetic reaction rates
    !call this%update_rk_old() !> we update old aqueous kinetic reaction rates
    !call this%solid_chemistry%update_rk_old() !> we update old solid kinetic reaction rates
!> Post-process    
    deallocate(rk_old) !> we deallocate new kinetic reaction rates
end subroutine