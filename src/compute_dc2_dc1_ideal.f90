 !> Computes Jacobian of secondary species concentrations assuming constant activity coefficients
!> \f$\partial c_{2} / \partial c_1\f$ with \f$gamma\f$ constant
!! We assume all primary species are aqueous
subroutine compute_dc2_dc1_ideal(this,c1,c2,dc2_dc1)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(in) :: c1(:) !> chapuza (dim=n_p)
    real(kind=8), intent(in) :: c2(:) !> chapuza (dim=n_eq)
    real(kind=8), intent(out) :: dc2_dc1(:,:)
!> Variables
    integer(kind=4) :: i,j    
!> Process
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
            dc2_dc1(i,j)=c2(i)*this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star(i,j)/c1(j)
        end do
    end do
end subroutine 