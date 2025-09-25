!> Computes residual of component concentrations
!!> The component matrix does not include constant activity species
!!> residual: \f$ res = U*c_{nc} - u \f$
subroutine compute_residual(this,conc_comp,c_nc,residual)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
    !> Variables
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(in) :: conc_comp(:)
    real(kind=8), intent(in) :: c_nc(:)
    real(kind=8), intent(out) :: residual(:) !> residual of component concentrations
    
    !> Process
    residual=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,c_nc) - conc_comp
end subroutine