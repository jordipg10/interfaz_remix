!!> Computes chemical part of aqueous species concentrations for a WMA-EE iteration assuming there are only kinetic reactions
!!> Esta subrutina tiene que ser de la clase reaccion local
!subroutine reaction_iteration_EE_WMA_kin_aq_chem(this,porosity,Delta_t,conc_react)
!    use aqueous_chemistry_m
!    implicit none
!    class(aqueous_chemistry_c) :: this
!    real(kind=8), intent(in) :: porosity
!    real(kind=8), intent(in) :: Delta_t !> time step
!    real(kind=8), intent(out) :: conc_react(:) !> reaction part of concentrations (must be allocated)
!    
!    call this%compute_rk_aq_chem() !> we compute kinetic reaction rates
!    conc_react=Delta_t*matmul(transpose(this%solid_chemistry%reactive_zone%chem_syst%Sk),this%rk)/porosity
!end subroutine