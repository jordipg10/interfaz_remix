!> Computes concentration of secondary species from activities of primary species explicitly using mass action law
!! We assume ideal activity coefficients for ions and that water has act=1
subroutine compute_c2_from_a1_ideal(this,a1,log_gamma2,c2)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: a1(:) !> chapuza (dim=n_p)
    real(kind=8), intent(in) :: log_gamma2(:) !> chapuza (dim=n_eq)
    real(kind=8), intent(out) :: c2(:) !> secondary concentrations (must be already allocated)
    
    integer(kind=4) :: n_p,n_e,n_sec_aq,n_nc_aq
    real(kind=8), allocatable :: log_a2(:), log_c2(:)
    
    !n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    !n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    !n_sec_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
    
!> mass action law
    log_a2=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,log10(a1))+&
        this%solid_chemistry%reactive_zone%speciation_alg%logK_tilde
    
    log_c2=log_a2-log_gamma2 !> we compute log_10 secondary species concentrations
    c2=10**log_c2 !> we compute secondary species concentrations
    call this%set_conc_sec_species(c2) !> we set secondary species concentrations
 end subroutine