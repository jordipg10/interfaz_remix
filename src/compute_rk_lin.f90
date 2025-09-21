!> Computes linear kinetic reaction rate
subroutine compute_rk_lin(this,conc,rk)
    use lin_kin_reaction_m, only: lin_kin_reaction_c
    implicit none
    class(lin_kin_reaction_c) :: this
    real(kind=8), intent(in) :: conc
    real(kind=8), intent(out) :: rk
        
   rk=this%lambda*conc
end subroutine