!> Computes kinetic reaction rates associated to aqueous chemistry object
subroutine compute_rk(this,rk)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    real(kind=8), intent(out) :: rk(:) !>  (must be allocated)
!> Variables
    integer(kind=4) :: i,n,niter,rk_ind,l,index,num_rk
    integer(kind=4), allocatable :: indices(:)
    real(kind=8), allocatable :: drk_dc_loc(:)
    real(kind=8) :: saturation

    num_rk=0 !> we initialise counter for kinetic reactions
    
!!> Update old kinetic reaction rates
!    call this%update_rk_old() !> we update old values of kinetic reaction rates
!    call this%solid_chemistry%update_rk_old() !> we update old values of kinetic reaction rates
!> We compute linear kinetic reaction rates
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
        num_rk=num_rk+1
        !index=this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%indices_aq_phase(1)
        !allocate(drk_dc_loc(1))
        call this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%compute_rk_lin(&
            this%concentrations(this%indices_aq_species(&
            this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%indices_aq_phase(1))),this%rk(i))
        rk(i)=this%rk(i) !> chapuza
        !allocate(drk_dc_loc(size(this%indices_rk%cols(num_rk)%col_1)))
        !call this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%compute_drk_dc_lin(drk_dc_loc)
        !drk_dc(num_rk,this%indices_rk%cols(num_rk)%col_1) = drk_dc_loc
        !deallocate(drk_dc_loc)
        ! drk_dc(i,this%indices_aq_species(index))=drk_dc_loc(1) !> chapuza
        !deallocate(drk_dc_loc)
    end do
    !num_rk=num_rk+this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
!> We compute mineral kinetic reaction rates
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
        !num_rk=num_rk+1
        !indices=this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(&
        !    this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase
        !allocate(drk_dc_loc(this%indices_rk%cols(num_rk+i)%dim))
        saturation=this%compute_saturation_kin_min(this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))
        call this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%compute_rk_mineral(&
            this%activities(this%indices_aq_species(this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%params%cat_indices)),saturation,&
            this%solid_chemistry%react_surfaces(i),this%solid_chemistry%temp,this%solid_chemistry%rk(i))
        rk(num_rk+i)=this%solid_chemistry%rk(i) !> chapuza
        !call this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
        !    this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%compute_drk_dc_mineral(&
        !        this%concentrations(this%indices_aq_species(&
        !        this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(&
        !        this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase)),&
        !        this%activities(this%indices_aq_species(this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
        !        this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%params%cat_indices)),saturation,&
        !        this%solid_chemistry%react_surfaces(i),this%solid_chemistry%temp,drk_dc_loc)
        !drk_dc(num_rk+i,this%indices_rk%cols(num_rk+this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%col_1)=drk_dc_loc !> chapuza
        !deallocate(drk_dc_loc)
    end do
    num_rk=num_rk+this%solid_chemistry%mineral_zone%num_minerals_kin
!> We compute Monod reaction rates
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        num_rk=num_rk+1
        ! indices=this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase
        !allocate(drk_dc_loc(this%indices_rk%cols(num_rk)%dim))
        call this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%compute_rk_Monod(this%concentrations(&
            this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase)),&
            this%rk(num_rk-this%solid_chemistry%mineral_zone%num_minerals_kin)) !> chapuza
        rk(num_rk)=this%rk(num_rk-this%solid_chemistry%mineral_zone%num_minerals_kin) !> chapuza
        !call this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%compute_drk_dc_Monod(this%concentrations(&
        !    this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase)),&
        !    this%rk(num_rk-this%solid_chemistry%mineral_zone%num_minerals_kin),drk_dc_loc)
        !drk_dc(num_rk,this%indices_rk%cols(num_rk)%col_1)=drk_dc_loc !> chapuza
        !deallocate(drk_dc_loc)
    end do
end subroutine