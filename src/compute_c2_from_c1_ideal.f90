!> Computes concentration of secondary species from concentration of primary species explicitly using mass action law
!! We assume ideal activity coefficients
subroutine compute_c2_from_c1_ideal(this,c1,c2)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst
    implicit none
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: c1(:) !> chapuza (dim=n_p)
    real(kind=8), intent(out) :: c2(:) !> secondary concentrations (must be already allocated)
    
    integer(kind=4) :: n_p,n_e,n_sec_aq,n_nc_aq
    real(kind=8), allocatable :: log_c2(:)
    
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    n_sec_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
            
    log_c2=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,log10(c1))+&
    this%solid_chemistry%reactive_zone%speciation_alg%logK_tilde !> mass action law
    c2=10**log_c2
    call this%set_conc_sec_species(c2)
 end subroutine