 !> Computes Jacobian of secondary species concentrations assuming constant activity coefficients
!> \f$\partial c_{2,nc} / \partial c_1\f$ with \f$gamma\f$ constant
!! We assume primary species are aqueous & solid
subroutine compute_dc2nc_dc1_ideal(this,c1,c2nc,dc2nc_dc1)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(in) :: c1(:) !> chapuza (dim=n_p)
    real(kind=8), intent(in) :: c2nc(:) !> chapuza (dim=n_eq)
    real(kind=8), intent(out) :: dc2nc_dc1(:,:)
!> Variables
    integer(kind=4) :: i,j    
!> Process
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
            dc2nc_dc1(i,j)=c2nc(i)*this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star(i,j)/c1(j)
        end do
    end do
end subroutine