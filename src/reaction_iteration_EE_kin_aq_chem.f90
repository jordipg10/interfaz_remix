!> Computes chemical part of variable activity concentrations for a reactive mixing Euler explicit iteration assuming there are only kinetic reactions
subroutine reaction_iteration_EE_kin(this,Delta_t,rk_tilde,conc_react)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    !real(kind=8), intent(in) :: porosity !> porosity
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rates after mixing
    real(kind=8), intent(out) :: conc_react(:) !> reaction part of variable activity concentrations (must be already allocated)
!> Variables
    real(kind=8), allocatable :: rk(:) !> kinetic reaction rates
! Pre-process
    allocate(rk(this%solid_chemistry%mineral_zone%num_minerals_kin+&
        this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts)) !> we allocate kinetic reaction rates
!> Process
    call this%compute_rk(rk) !> we compute kinetic reactions rates at current time step (CHAPUZA)
    conc_react=Delta_t*matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,rk_tilde) !> we compute reaction part of variable activity concentrations
!> Post-process
    deallocate(rk) !> we deallocate kinetic reaction rates
end subroutine