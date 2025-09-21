!> Computes concentration of secondary species from concentration of primary species using Picard method and mass action law
!! This subroutine is meant to be used when primary species are aqueous & solid
subroutine compute_c2_from_c1_Picard(this,c1,c2_ig,c2,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst

    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(inout) :: c1(:) !> primary concentrations
    real(kind=8), intent(in) :: c2_ig(:) !> initial guess secondary variable activity concentrations
    real(kind=8), intent(out) :: c2(:) !> secondary variable activity concentrations (must be already allocated)
    integer(kind=4), intent(out) :: niter !> number of iterations
    logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
!> Variables
    integer(kind=4) :: n_p,n_p_aq,n_2_aq,n_nc,n_e,n_aq,i,k
    real(kind=8), allocatable :: log_gamma1_old(:),log_gamma2_old(:),log_c2_new(:),c2_old(:)
!> Pre-process
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    n_2_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
    n_aq=this%aq_phase%num_species
    
    allocate(log_gamma1_old(n_p),log_gamma2_old(n_e),log_c2_new(n_e),c2_old(n_e))
    
    log_gamma1_old=0d0 !> chapuza
    log_gamma2_old=this%get_log_gamma2() !> chapuza
    
    c2_old=c2_ig
    
    niter=0 !> we initialise number of iterations
    CV_flag=.false. !> we initialise convergence flag
    !call this%compute_molalities()
    k=0
    do
        niter=niter+1 !> we update number of iterations
        !call this%compute_molalities() !> we change units to compute ionic activity
        call this%compute_ionic_act() !> we compute ionic activity
        !call this%compute_salinity()
        !call this%compute_molarities()
    !> We compute log_10 of activity coefficients of variable activity species
        call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_act,this%params_aq_sol,this%log_act_coeffs) !> we compute log_10 activity coefficients aqueous variable activity species
        call this%compute_activities() !> we compute activities
        call this%compute_log_act_coeff_wat() !> we compute log_10 activity coefficient of water
        log_gamma1_old(1:n_p_aq)=this%log_act_coeffs(1:n_p_aq)
        log_gamma2_old(1:n_2_aq)=this%log_act_coeffs(n_p_aq+1:n_aq)
        if (associated(this%gas_chemistry)) then !> chapuza
            call this%gas_chemistry%compute_log_act_coeffs_gases()
            log_gamma2_old(n_2_aq+1:n_e)=this%gas_chemistry%log_act_coeffs
        end if
    !> We compute concentration of free site
        c1(n_p)=1d0-sum(this%solid_chemistry%activities(this%solid_chemistry%num_solids-&
        this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl+2:this%solid_chemistry%num_solids))
        if (c1(n_p)<0d0) then
            k=k+1
            call this%set_conc_sec_species(c2_ig/(10**k))
        else
        !> We apply mass action law to compute concentration secondary variable activity species
            log_c2_new=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,log_gamma1_old+log10(c1))+&
            this%solid_chemistry%reactive_zone%speciation_alg%logK_tilde-log_gamma2_old !> mass action law
            c2=10**log_c2_new
        !> We update secondary concentrations
            call this%set_conc_sec_species(c2)
            !if (associated(this%gas_chemistry)) then !> chapuza
            !    !call this%gas_chemistry%update_conc_gases(c2(n_2_aq+1:n_e)*this%volume) !> we update moles of gases
            !    call this%gas_chemistry%compute_vol_gas_conc() !> we compute total volume of gas
            !end if
        !< We check convergence
             if (inf_norm_vec_real((c2-c2_old)/c2_old)<this%solid_chemistry%reactive_zone%CV_params%rel_tol) then
                CV_flag=.true.
                exit !> CV reached
            else if (niter==this%solid_chemistry%reactive_zone%CV_params%niter_max) then
                print *, inf_norm_vec_real((c2-c2_old)/c2_old)
                error stop "Too many Picard iterations in speciation"
            else
                c2_old=c2 !> we update concentrations for next iteration
            end if
        end if
    end do
!> Post-process
    !call this%compute_molalities() !> we change units to compute ionic activity
    call this%compute_ionic_act() !> we compute ionic activity
    call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_act,this%params_aq_sol,this%log_act_coeffs) !> we compute log_10 activity coefficients aqueous variable activity species
    call this%compute_activities() !> we compute activities
    !if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
    !    if (ASSOCIATED(this%solid_chemistry%mineral_zone)) then
    !        this%solid_chemistry%activities(this%solid_chemistry%mineral_zone%num_minerals+1)=1d0-SUM(this%solid_chemistry%activities(this%solid_chemistry%mineral_zone%num_minerals+2:this%solid_chemistry%mineral_zone%num_minerals+this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl))
    !        !this%solid_chemistry%log_act_coeffs(this%solid_chemistry%mineral_zone%num_minerals+1)=LOG10(this%solid_chemistry%activities(this%solid_chemistry%mineral_zone%num_minerals+1))-LOG10(this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals+1))
    !    else
    !        this%solid_chemistry%activities(this%solid_chemistry%reactive_zone%num_minerals+1)=1d0-SUM(this%solid_chemistry%activities(this%solid_chemistry%reactive_zone%num_minerals+2:this%solid_chemistry%reactive_zone%num_solids))
    !        !this%solid_chemistry%log_act_coeffs(this%solid_chemistry%reactive_zone%num_minerals+1)=LOG10(this%solid_chemistry%activities(this%solid_chemistry%reactive_zone%num_minerals+1))-LOG10(this%solid_chemistry%concentrations(this%solid_chemistry%reactive_zone%num_minerals+1))
    !    end if
    !end if
    call this%compute_log_act_coeff_wat() !> we compute log_10 activity coefficient of water
    call this%compute_salinity()
    !call this%compute_molarities()
    if (associated(this%gas_chemistry)) then !> chapuza
        !call this%gas_chemistry%update_conc_gases(c2nc(n_nc2_aq+1:n_e)*this%volume) !> we update moles of gases
        call this%gas_chemistry%compute_log_act_coeffs_gases() !> we compute log_10 activity coefficients of gases
        call this%gas_chemistry%compute_partial_pressures() !> we compute activities (ie. partial pressures)
        call this%gas_chemistry%compute_pressure() !> we compute activities (ie. partial pressures)
        call this%gas_chemistry%compute_vol_gas_conc() !> we compute total volume of gas
    end if
 end subroutine