!> Computes Jacobian of ionic activity with respect to primary cocnentrationa
subroutine compute_dI_dc1(this,dc2aq_dc1,dI_dc1)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: dc2aq_dc1(:,:)
    real(kind=8), intent(out) :: dI_dc1(:) !> Jacobian of ionic activity with respect to primary cocnentrationa (must be already allocated)
    
    integer(kind=4) :: i
    
    
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
        dI_dc1(i)=5d-1*this%aq_phase%aq_species(i)%valence**2 + &
        dot_product(this%solid_chemistry%reactive_zone%chem_syst%z2(&
        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%aq_phase%num_species),dc2aq_dc1(:,i))
    end do

end subroutine