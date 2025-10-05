!!> Computes chemical part of aqueous component concentrations for a WMA-EE iteration assuming there are equilibrium and kinetic reactions
!!> Esta subrutina tiene que ser de la clase reaccion local
!subroutine reaction_iteration_EE_WMA_eq_kin_aq_chem(this,porosity,Delta_t,conc_comp_react)
!    use aqueous_chemistry_m
!    implicit none
!    class(aqueous_chemistry_c) :: this
!    real(kind=8), intent(in) :: porosity !> porosity
!    real(kind=8), intent(in) :: Delta_t !> time step
!    real(kind=8), intent(out) :: conc_comp_react(:) !> reaction part of component concentrations (must be allocated)
!            
!    call this%compute_rk_aq_chem() !> we compute kinetic reaction rates
!    conc_comp_react=Delta_t*matmul(solid_chemistry%reactive_zone%U_SkT_prod,this%rk)/porosity
!end subroutine