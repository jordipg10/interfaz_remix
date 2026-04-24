!> Computes chemical part of variable activity concentrations for a reactive mixing Euler explicit iteration assuming there are only kinetic reactions
subroutine compute_react_term_EE_kin(this,Delta_t,lambda_r,conc_react)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: lambda_r !> kinetic reaction rate contributions after mixing
    real(kind=8), intent(out) :: conc_react(:) !> reaction part of variable activity concentrations (must be already allocated)
!> Variables
    real(kind=8), allocatable :: rk_old(:) !> kinetic reaction rates
! Pre-process
    !allocate(rk_old(this%solid_chemistry%mineral_zone%num_minerals_kin+&
    !    this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts)) !> we allocate kinetic reaction rates
!> Process
    rk_old=this%get_rk_old() !> we compute kinetic reactions rates at current time step (CHAPUZA)
    conc_react=Delta_t*lambda_r*matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,rk_old) !> we compute reaction part of variable activity concentrations
!> Post-process
    deallocate(rk_old) !> we deallocate kinetic reaction rates
end subroutine