!> Computes Jacobian of kinetic reaction rates with respect to aqueous concentrations
subroutine compute_Jacobian_rk_anal(this,drk_dc)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(out) :: drk_dc(:,:) !> Jacobian (must be allocated)
    
    integer(kind=4) :: i,k,n,m,l,rk_ind,cat_ind,p,inh_ind,DOC_ind,num_kin
    integer(kind=4), allocatable :: indices(:)
    real(kind=8), allocatable :: act_cat(:),drk_dc_loc(:)
    real(kind=8) :: saturation

    drk_dc=0d0 !> zero Jacobian by default
    num_kin=0 !> we initialise number of kinetic reactions
!> Linear kinetic reactions
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
        num_kin=num_kin+1
        call this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%compute_drk_dc_lin(drk_dc(num_kin,:))
    end do
    !num_kin=num_kin+this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
!> Mineral kinetic reactions
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
        num_kin=num_kin+1
        !indices=this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(i)%indices_aq_phase
        allocate(drk_dc_loc(this%indices_rk%cols(num_kin)%dim))
        saturation=this%compute_saturation_kin_min(this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))
        call this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%compute_drk_dc_mineral(this%concentrations(&
            this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase)),&
            this%activities(this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase)),&
            this%log_act_coeffs(this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase)),&
            this%activities(this%indices_aq_species(&
            this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(i)%params%cat_indices)),saturation,&
            this%solid_chemistry%react_surfaces(i),this%solid_chemistry%temp,drk_dc_loc)
        drk_dc(num_kin,this%indices_rk%cols(num_kin)%col_1)=drk_dc_loc !> chapuza
        deallocate(drk_dc_loc)
    end do
    !num_kin=num_kin+this%solid_chemistry%mineral_zone%num_minerals_kin
!> Redox kinetic reactions
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        !indices=this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase
        allocate(drk_dc_loc(this%indices_rk%cols(num_kin)%dim))
        call this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%compute_drk_dc_Monod(this%concentrations(&
            this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase)),&
            this%rk(num_kin-this%solid_chemistry%mineral_zone%num_minerals_kin),drk_dc_loc)
        drk_dc(num_kin,this%indices_rk%cols(num_kin)%col_1)=drk_dc_loc !> chapuza
        deallocate(drk_dc_loc)
    end do
end subroutine