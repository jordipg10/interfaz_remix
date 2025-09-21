!> Computes concentration of secondary species with variable activity from concentration of primary species using Picard method and mass action law
!! This subroutine is meant to be used when all primary species are aqueous
subroutine compute_c2nc_from_c1_aq_Picard(this,c2nc_ig,c2nc,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst

    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: c2nc_ig(:) !> initial guess secondary variable activity concentrations
    real(kind=8), intent(out) :: c2nc(:) !> secondary variable activity concentrations (must be already allocated)
    integer(kind=4), intent(out) :: niter !> number of iterations
    logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
!> Variables
    integer(kind=4) :: n_p,n_nc2_aq,n_nc,n_e,n_nc_aq,i
    real(kind=8), allocatable :: log_gamma1_old(:),log_gamma2nc_old(:),log_c2nc_new(:),c2nc_old(:)
!> Pre-process
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    n_nc2_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
    n_nc_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
    
    allocate(log_gamma1_old(n_p),log_gamma2nc_old(n_e),log_c2nc_new(n_e),c2nc_old(n_e))
    
    log_gamma2nc_old=this%get_log_gamma2nc() !> chapuza
    c2nc_old=c2nc_ig
    niter=0 !> we initialise number of iterations
    CV_flag=.false. !> we initialise convergence flag
    !call this%compute_molalities()
    do
        niter=niter+1 !> we update number of iterations
        !call this%compute_molalities() !> we change units to compute ionic activity
        call this%compute_ionic_act() !> we compute ionic activity
        !call this%compute_salinity()
        !call this%compute_molarities()
    !> We compute log_10 of activity coefficients of variable activity species
        call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_act,this%params_aq_sol,this%log_act_coeffs(1:n_nc_aq))
        log_gamma1_old=this%log_act_coeffs(1:n_p)
        log_gamma2nc_old(1:n_nc2_aq)=this%log_act_coeffs(n_p+1:n_nc_aq)
        if (associated(this%gas_chemistry)) then
            if (this%gas_chemistry%reactive_zone%gas_phase%num_var_act_species>0) then !> chapuza
                call this%gas_chemistry%compute_log_act_coeffs_gases()
                log_gamma2nc_old(n_nc2_aq+1:n_e)=this%gas_chemistry%log_act_coeffs
            end if
        end if
    !> We apply mass action law to compute concentration secondary variable activity species
        log_c2nc_new=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,log_gamma1_old+&
        log10(this%concentrations(1:n_p)))+this%solid_chemistry%reactive_zone%speciation_alg%logK_star-log_gamma2nc_old !> mass action law
        c2nc=10**log_c2nc_new
    !> We update secondary aqueous variable activity concentrations in aqueous chemistry object
        call this%set_conc_sec_aq_var_act_species(c2nc(1:n_nc2_aq))
        if (associated(this%gas_chemistry)) then
            if (this%gas_chemistry%reactive_zone%gas_phase%num_var_act_species>0) then !> chapuza
                call this%gas_chemistry%update_conc_gases(c2nc(n_nc2_aq+1:n_e)*this%volume) !> we update moles of gases
                !call this%gas_chemistry%compute_vol_gas() !> we compute total volume of gas
            end if
        end if
    !< We check convergence
         if (inf_norm_vec_real((c2nc-c2nc_old)/c2nc_old)<this%solid_chemistry%reactive_zone%CV_params%rel_tol) then
            CV_flag=.true.
            exit !> CV reached
        else if (niter==this%solid_chemistry%reactive_zone%CV_params%niter_max) then
            print *, inf_norm_vec_real((c2nc-c2nc_old)/c2nc_old)
            error stop "Too many Picard iterations in speciation"
        else
            c2nc_old=c2nc !> we update concentrations for next iteration
        end if
    end do
!> Post-process
    !call this%compute_molalities() !> we change units to compute ionic activity
    call this%compute_ionic_act() !> we compute ionic activity
    call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_act,this%params_aq_sol,this%log_act_coeffs(1:n_nc_aq)) !> we compute log_10 activity coefficients aqueous variable activity species
    call this%compute_activities_aq() !> we compute activities
    call this%compute_log_act_coeff_wat() !> we compute log_10 activity coefficient of water
    call this%compute_salinity()
    !call this%compute_molarities()
    if (associated(this%gas_chemistry)) then !> chapuza
        !call this%gas_chemistry%update_conc_gases(c2nc(n_nc2_aq+1:n_e)*this%volume) !> we update moles of gases
        !call this%gas_chemistry%compute_vol_gas() !> we compute total volume of gas
        call this%gas_chemistry%compute_log_act_coeffs_gases() !> we compute log_10 activity coefficients of gases
        call this%gas_chemistry%compute_partial_pressures() !> we compute activities (ie. partial pressures)
    end if
 end subroutine