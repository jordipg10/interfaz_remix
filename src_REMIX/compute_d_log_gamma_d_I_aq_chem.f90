!> Computes derivative of log_10(activity coefficients) of aqueous variable activity species with respect to ionic activity
subroutine compute_d_log_gamma_d_I_aq_chem(this,d_log_gamma_d_I)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst
    implicit none
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(inout) :: d_log_gamma_d_I(:) !> derivative of log_10(activity coefficients) with respect to ionic activity (must be already allocated)
    
    integer(kind=4) :: i
    
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
        d_log_gamma_d_I(i)=-this%aq_phase%aq_species(i)%params_act_coeff%alpha*(&
        this%params_aq_sol%A*this%aq_phase%aq_species(i)%valence**2)/(2d0*sqrt(this%ionic_act)*(&
        1d0+this%aq_phase%aq_species(i)%params_act_coeff%beta*sqrt(this%ionic_act))**2) + &
        this%aq_phase%aq_species(i)%params_act_coeff%gamma
    end do
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
        d_log_gamma_d_I(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i)=&
        -this%aq_phase%aq_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)%params_act_coeff%alpha*&
        (this%params_aq_sol%A*this%aq_phase%aq_species(&
        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)%valence**2)/(2d0*sqrt(this%ionic_act)&
        *(1d0+this%aq_phase%aq_species(&
        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)%params_act_coeff%beta*sqrt(this%ionic_act))**2) +&
         this%aq_phase%aq_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)%params_act_coeff%gamma
    end do
end subroutine