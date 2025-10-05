!> Computes concentration of secondary species with variable activity from concentration of primary species explicitly using mass action law
!!  We assume the activity coefficients are ideal
!! This subroutine is meant to be used when primary species are aqueous & solid
subroutine compute_c2nc_from_c1_ideal(this,c1,c2nc)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst
    implicit none
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: c1(:) !> chapuza (dim=n_p)
    real(kind=8), intent(out) :: c2nc(:) !> secondary variable activity concentrations (must be already allocated)
    
    integer(kind=4) :: n_p,n_nc2_aq,n_e
    real(kind=8), allocatable :: log_c2nc(:)
    
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    n_nc2_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species

    log_c2nc=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,log10(c1))+&
    this%solid_chemistry%reactive_zone%speciation_alg%logK_star
    c2nc=10**log_c2nc
    call this%set_conc_sec_var_act_species(c2nc)
 end subroutine