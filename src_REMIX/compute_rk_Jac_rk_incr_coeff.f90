 !> Computes kinetic reaction rates associated to aqueous chemistry object
subroutine compute_rk_Jac_rk_incr_coeff(this,drk_dc)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, int_array_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    real(kind=8), intent(out) :: drk_dc(:,:) !> Jacobian (must be allocated)
!> Variables
    integer(kind=4) :: i,j,niter,rk_ind,l,index,num_rk
    integer(kind=4), allocatable :: kin_ind(:),species_indices(:)
    real(kind=8), allocatable :: rk_pert(:),c_pert(:)
    real(kind=8) :: saturation,saturation_pert

    type(int_array_c) :: indices_Monod,indices_min !> indices in aqueous phase of species relevant for Monod and mineral reaction rates
    
    allocate(rk_pert(this%solid_chemistry%reactive_zone%chem_syst%num_kin_reacts))
    
    drk_dc=0d0 !> Jacobian is zero by default
    num_rk=0 !> we initialise counter for kinetic reactions
    
!!> Update old kinetic reaction rates
!    call this%update_rk_old() !> we update old values of kinetic reaction rates
!    call this%solid_chemistry%update_rk_old() !> we update old values of kinetic reaction rates
!> We compute linear kinetic reaction rates
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
        num_rk=num_rk+1
        !index=this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%indices_aq_phase(1)
        !allocate(conc_kin(this%solid_chemistry%reactive_zone%chem_syst%num_species),kin_ind(this%solid_chemistry%reactive_zone%chem_syst%num_species))
        call this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%compute_rk_lin(this%concentrations(&
            this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%indices_aq_phase(1))),&
            this%rk(i))
    end do
!> We compute mineral kinetic reaction rates
    call indices_min%allocate_array(this%solid_chemistry%mineral_zone%num_minerals_kin)
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
        num_rk=num_rk+1
        call indices_min%cols(i)%set_dim(this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(i)%num_aq_rk)
        indices_min%cols(i)%col_1=this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(i)%indices_aq_phase !> falta un set
        saturation=this%compute_saturation_kin_min(this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))
        call this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(this%solid_chemistry%mineral_zone%ind_min_chem_syst(i)&
            )%compute_rk_mineral(this%activities(this%indices_aq_species(&
            this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%params%cat_indices)),saturation,&
            this%solid_chemistry%react_surfaces(i),this%solid_chemistry%temp,this%solid_chemistry%rk(i))
    end do
!> We compute Monod reaction rates
    call indices_Monod%allocate_array(this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts)
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        num_rk=num_rk+1
        call indices_Monod%cols(i)%set_dim(this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%num_aq_rk)
        indices_Monod%cols(i)%col_1=this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase !> falta un set
        call this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%compute_rk_Monod(this%concentrations(&
            this%indices_aq_species(indices_Monod%cols(i)%col_1)),&
            this%rk(num_rk-this%solid_chemistry%mineral_zone%num_minerals_kin))
    end do
    num_rk=0 !> we reinitialise counter for kinetic reactions
!> We compute perturbed mineral kinetic reaction rates
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
        num_rk=num_rk+1
        do j=1,indices_min%cols(i)%dim
         !> chapuza
            this%concentrations(this%indices_aq_species(indices_min%cols(i)%col_1(j)))=&
                this%concentrations(this%indices_aq_species(indices_min%cols(i)%col_1(j)))+&
                this%solid_chemistry%reactive_zone%CV_params%eps
            call this%compute_ionic_act()
            call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_act,this%params_aq_sol,this%log_act_coeffs)
            call this%compute_activities()
            saturation_pert=this%compute_saturation_kin_min(this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))
            call this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(i)%compute_rk_mineral(this%activities(&
                this%indices_aq_species(this%solid_chemistry%mineral_zone%chem_syst%min_kin_reacts(&
                this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%params%cat_indices)),saturation_pert,&
                this%solid_chemistry%react_surfaces(i),this%solid_chemistry%temp,rk_pert(num_rk))
            drk_dc(num_rk,this%indices_rk%cols(i)%col_1(j))=(rk_pert(num_rk)-this%rk(num_rk))/&
                this%solid_chemistry%reactive_zone%CV_params%eps
        !> chapuza
            this%concentrations(this%indices_aq_species(indices_min%cols(i)%col_1(j)))=this%concentrations(&
            this%indices_aq_species(indices_min%cols(i)%col_1(j)))-this%solid_chemistry%reactive_zone%CV_params%eps
            call this%compute_ionic_act()
            call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_act,this%params_aq_sol,this%log_act_coeffs)
            call this%compute_activities()
        end do
    end do
!> We compute perturbed Monod reaction rates
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        num_rk=num_rk+1
        do j=1,indices_Monod%cols(i)%dim !> chapuza
            c_pert=this%concentrations(this%indices_aq_species(indices_Monod%cols(i)%col_1))
            c_pert(j)=c_pert(j)+this%solid_chemistry%reactive_zone%CV_params%eps
            call this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%compute_rk_Monod(c_pert,rk_pert(num_rk))
            drk_dc(num_rk,this%indices_rk%cols(i)%col_1(j))=(rk_pert(num_rk)-this%rk(num_rk))/&
                this%solid_chemistry%reactive_zone%CV_params%eps
        end do
    end do
end subroutine