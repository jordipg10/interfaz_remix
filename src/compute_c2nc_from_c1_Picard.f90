!> Computes concentration of secondary species with variable activity from concentration of primary species using Picard method and mass action law
!! This subroutine is meant to be used when primary species are aqueous & solid
subroutine compute_c2v_from_c1_Picard(this,c1,c2v_ig,c2v,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    use vectors_m, only: inf_norm_vec_real
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: c1(:) !> chapuza (dim=n_p)
    real(kind=8), intent(in) :: c2v_ig(:) !> initial guess secondary variable activity concentrations
    real(kind=8), intent(out) :: c2v(:) !> secondary variable activity concentrations (must be already allocated)
    integer(kind=4), intent(out) :: niter !> number of iterations
    logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
!> Variables
    integer(kind=4) :: n_p,n_p_aq,n_v2_aq,n_v,n_e,n_v_aq,i
    real(kind=8), allocatable :: log_gamma1_old(:),log_gamma2v_old(:),log_c2v_new(:),c2v_old(:)
!> Pre-process
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    n_v2_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
    n_v_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
    
    allocate(log_gamma1_old(n_p),log_c2v_new(n_e),c2v_old(n_e))
    
    log_gamma1_old=this%get_log_gamma1() !> chapuza
    log_gamma2v_old=this%get_log_gamma2v() !> chapuza
    
    !c1=this%get_c1()
    c2v_old=c2v_ig
    
    niter=0 !> we initialise number of iterations
    CV_flag=.false. !> we initialise convergence flag
    !call this%compute_molalities()
    do
        niter=niter+1 !> we update number of iterations
        !call this%compute_molalities() !> we change units to compute ionic activity
        call this%compute_ionic_strength() !> we compute ionic activity
        !call this%compute_salinity()
        !call this%compute_molarities()
    !> We compute log_10 of activity coefficients of aqueous variable activity species
        call this%compute_log_act_coeffs_aq_chem()
        log_gamma1_old(1:n_p_aq)=this%log_act_coeffs(1:n_p_aq)
        log_gamma2v_old(1:n_v2_aq)=this%log_act_coeffs(n_p_aq+1:n_v_aq)
        call this%compute_activities() !> we compute activities
        call this%compute_log_act_coeff_wat() !> we compute log_10 activity coefficient of water
        ! c1(n_p)=1d0-sum(this%solid_chemistry%activities(this%solid_chemistry%num_solids-&
        ! this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl+2:this%solid_chemistry%num_solids)) !> we compute concentration of free site
    !> We apply mass action law to compute concentration secondary variable activity species
        log_c2v_new=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,log_gamma1_old+log10(c1))+&
            this%solid_chemistry%reactive_zone%speciation_alg%logK_star-log_gamma2v_old !> mass action law
        c2v=10**log_c2v_new
    !> We update secondary aqueous variable activity concentrations in aqueous chemistry object
        call this%set_conc_sec_var_act_species(c2v)
        !if (associated(this%gas_chemistry)) then !> chapuza
        !    call this%gas_chemistry%update_conc_gases(c2v(n_v2_aq+1:n_e)*this%volume) !> we update moles of gases
        !    call this%gas_chemistry%compute_vol_gas_species_conc() !> we compute total volume of gas
        !end if
    !< We check convergence
         if (inf_norm_vec_real((c2v-c2v_old)/c2v_old)<this%solid_chemistry%reactive_zone%CV_params%rel_tol) then
            CV_flag=.true.
            exit !> CV reached
        else if (niter==this%solid_chemistry%reactive_zone%CV_params%niter_max) then
            print *, inf_norm_vec_real((c2v-c2v_old)/c2v_old)
            error stop "Too many Picard iterations in speciation"
        else
            c2v_old=c2v !> we update concentrations for next iteration
        end if
    end do
!> Post-process
    !call this%compute_molalities() !> we change units to compute ionic activity
    call this%compute_ionic_strength() !> we compute ionic activity
    call this%compute_log_act_coeffs_aq_chem() !> we compute log_10 activity coefficients aqueous variable activity species
    call this%compute_activities() !> we compute activities
    call this%compute_log_act_coeff_wat() !> we compute log_10 activity coefficient of water
    call this%compute_salinity()
    !call this%compute_molarities()
    ! if (associated(this%gas_chemistry)) then !> chapuza
    !     !call this%gas_chemistry%update_conc_gases(c2v(n_v2_aq+1:n_e)*this%volume) !> we update moles of gases
    !     !call this%gas_chemistry%compute_vol_gas() !> we compute total volume of gas
    !     call this%gas_chemistry%compute_log_act_coeffs_gases() !> we compute log_10 activity coefficients of gases
    !     call this%gas_chemistry%compute_partial_pressures() !> we compute activities (ie. partial pressures)
    !     call this%gas_chemistry%compute_pressure() !> we compute activities (ie. partial pressures)
    ! end if
 end subroutine