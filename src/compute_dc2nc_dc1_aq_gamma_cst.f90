!> Computes Jacobian of secondary variable activity species concentrations assuming constant activity coefficients
!> \f$\partial c_{2,nc} / \partial c_1\f$ with \f$gamma\f$ constant
!! We assume all primary species are aqueous
subroutine compute_dc2nc_dc1_aq_gamma_cst(this,dc2nc_dc1)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst
    implicit none
!> Arguments
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(out) :: dc2nc_dc1(:,:)
!> Variables
    integer(kind=4) :: i,j    
!> Process
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
            dc2nc_dc1(i,j)=this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i)*&
            this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star(i,j)/this%concentrations(j)
        end do
    end do
end subroutine