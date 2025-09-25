!> This subroutine computes the reaction rate of a binary system in equilibrium in 1D (De Simoni et al, 2005)
subroutine react_rate_bin_syst_eq_1D(this,u,du_dx,D,phi,r_eq)
    use eq_reaction_m, only : eq_reaction_c
    implicit none
    class(eq_reaction_c) :: this
    real(kind=8), intent(in) :: u !> concentration component at this target
    real(kind=8), intent(in) :: du_dx !> gradient of u at this target
    real(kind=8), intent(in) :: D !> dispersion coefficient
    real(kind=8), intent(in) :: phi !> porosity
    real(kind=8), intent(out) :: r_eq !> reaction rate
    
    real(kind=8) :: d2c2nc_du2
    
    d2c2nc_du2=2*this%eq_cst/((u**2+4*this%eq_cst)**1.5) !> second derivative of secondary concentrations with respect to component concentrations
    
    r_eq=phi*d2c2nc_du2*D*du_dx**2
   
end subroutine