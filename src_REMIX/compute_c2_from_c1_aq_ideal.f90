!> Computes concentration of secondary species from concentration of primary species explicitly using mass action law
!! We assume ideal activity coefficients
!! This subroutine is meant to be used only when all primary species are aqueous
subroutine compute_c2_from_c1_aq_ideal(this,c2)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst
    implicit none
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(out) :: c2(:) !> secondary concentrations (must be already allocated)
    
    integer(kind=4) :: n_p,n_e,n_sec_aq,n_nc_aq
    real(kind=8), allocatable :: log_c2(:)
    
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    n_sec_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
            
    log_c2=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,log10(this%concentrations(1:n_p)))+&
    this%solid_chemistry%reactive_zone%speciation_alg%logK_tilde !> mass action law
    c2=10**log_c2
    call this%set_conc_sec_aq_species(c2(1:n_sec_aq))
        if (associated(this%gas_chemistry)) then
        if (this%gas_chemistry%reactive_zone%gas_phase%num_var_act_species>0) then !> chapuza
            call this%gas_chemistry%update_conc_gases(c2(n_sec_aq+1:n_e)*this%volume) !> we update moles of gases
            !call this%gas_chemistry%compute_vol_gas() !> we compute total volume of gas
        end if
    end if
 end subroutine 