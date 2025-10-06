!> Computes concentration of secondary species from concentration of primary species explicitly using mass action law
!! We assume ideal activity coefficients of aqueous species
subroutine compute_conc_surf_ideal(this,conc_cats)
    use solid_chemistry_m, only: solid_chemistry_c
    implicit none
    class(solid_chemistry_c) :: this
    real(kind=8), intent(in) :: conc_cats(:) !> chapuza (dim=n� cation exchange half reactions)
    !real(kind=8), intent(in) :: act_ads_cats_ig(:) !> surface complex concentrations (must be already allocated)
    !integer(kind=4), intent(out) :: niter !> number of iterations
    !logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
    
    integer(kind=4) :: n_ads_cats,n_eq,n_sec_aq,n_gas_eq,i,j
    real(kind=8), allocatable :: act_ads_cats_old(:),log_act_ads_cats_new(:),act_ads_cats_new(:)
!> Pre-processing
    !n_p_aq=this%reactive_zone%speciation_alg%num_aq_prim_species
    n_eq=this%reactive_zone%speciation_alg%num_eq_reactions
    n_gas_eq=this%reactive_zone%gas_phase%num_gases_eq
    !n_ads_cats=this%reactive_zone%cat_exch_zone%num_exch_cats
    !allocate(Delta_act_ads_cats(n_ads_cats),residual(n_ads_cats),Jacobian(n_ads_cats,n_ads_cats))
    !allocate(log_act_ads_cats_new(n_ads_cats))
!> Process
    this%activities(this%num_solids)=this%reactive_zone%chem_syst%eq_reacts(&
        this%reactive_zone%ind_eq_reacts(n_eq-n_gas_eq))%eq_cst*conc_cats(2)/(1d0+&
        this%reactive_zone%chem_syst%eq_reacts(this%reactive_zone%ind_eq_reacts(n_eq-n_gas_eq))%eq_cst*conc_cats(2)+&
        this%reactive_zone%chem_syst%eq_reacts(this%reactive_zone%ind_eq_reacts(n_eq-n_gas_eq-1))%eq_cst*conc_cats(1))
    this%activities(this%num_solids-1)=this%reactive_zone%chem_syst%eq_reacts(&
        this%reactive_zone%ind_eq_reacts(n_eq-n_gas_eq-1))%eq_cst*conc_cats(1)*&
        this%activities(this%num_solids)/(this%reactive_zone%chem_syst%eq_reacts(&
        this%reactive_zone%ind_eq_reacts(n_eq-n_gas_eq))%eq_cst*conc_cats(2))
    this%activities(this%num_solids-2)=1d0-SUM(this%activities(this%num_solids-1:this%num_solids)) !> free site
!> Post-processing
    !call this%set_act_surf_compl(act_ads_cats_new)
    call this%compute_conc_ads_cats()
 end subroutine