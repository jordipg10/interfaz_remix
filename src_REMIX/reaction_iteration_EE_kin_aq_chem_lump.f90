!> Computes chemical part of variable activity concentrations for a reactive mixing Euler explicit iteration assuming there are only kinetic reactions
subroutine reaction_iteration_EE_kin_lump(this,Delta_t,conc_react)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    !real(kind=8), intent(in) :: porosity !> porosity
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(out) :: conc_react(:) !> reaction part of concentrations (must be already allocated)
!> Variables
    real(kind=8), allocatable :: rk(:) !> kinetic reaction rates
!> Pre-process
    allocate(rk(this%solid_chemistry%mineral_zone%num_minerals_kin+&
        this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts))
!> Process
    call this%compute_rk(rk) !> we compute kinetic reactions rates
    call this%compute_Rk_mean(1d0,Delta_t) !> we update mean reaction amounts
    conc_react=Delta_t*matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,rk)
!> Post-process
    deallocate(rk)
end subroutine