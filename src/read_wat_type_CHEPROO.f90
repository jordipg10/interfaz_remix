!> This subroutine reads a water type from CHEPROO data input
!! We assume that the units of concentrations in input file are molalities
subroutine read_wat_type_CHEPROO(this,n_p_aq,num_cstr,model,Jac_opt,unit,niter,CV_flag,surf_chem)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, solid_chemistry_c, reactive_zone_c, gas_chemistry_c, gas_c, mineral_c, & 
        mineral_zone_c, aq_species_c, species_c, gas_phase_c, cat_exch_c
    implicit none
    class(aqueous_chemistry_c) :: this
    integer(kind=4), intent(in) :: n_p_aq !> number of primary aqueous species
    integer(kind=4), intent(in) :: num_cstr !> number of constrains
    !integer(kind=4), intent(in) :: num_gas_cstr !> number of gas constrains
    integer(kind=4), intent(in) :: model !> activity coefficients model
    integer(kind=4), intent(in) :: Jac_opt !> 0: incremental coeffficinets, 1: analtical
    integer(kind=4), intent(in) :: unit !> file unit
    integer(kind=4), intent(out) :: niter !> number of iterations
    logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
    type(solid_chemistry_c), intent(in), optional :: surf_chem !> surface chemistry (if there are exchange reactions)
    !real(kind=8), intent(in), optional :: c1_surf !> chapuza
    !real(kind=8), intent(in), optional :: CEC !> chapuza
    !real(kind=8), intent(out), optional :: conc_exch(:) !> chapuza
    
    integer(kind=4) :: i,j,k,l,m,niwtype,nbwtype,nrwtype,gas_ind,min_ind,n_comp_aq,ind_cstr,n_gas_constr,n_aq_comp,icon,aq_sp_ind,&
    ind_sp
    integer(kind=4), allocatable :: icons(:),indices_constrains(:,:),gas_indices(:),n_icons(:),prim_indices(:),swap(:)
    real(kind=8), allocatable :: ctots(:),c2_init(:),c2_ig(:),a1(:),c1(:),c1_aq(:),act_ads_cats_ig(:),log_gamma2(:)
    character(len=256) :: prim_sp_name,label,aq_sp_name
    character(len=256), allocatable :: constrains(:)
    real(kind=8) :: temp,conc,ionic_act,guess,ctot
    logical :: flag_gas,flag_min,flag,flag_comp,flag_surf,flag_Se,flag_sp
    integer(kind=4), allocatable :: aux_ind_cstr(:,:),ind_mins(:)
    integer(kind=4), allocatable :: ind_swap(:)
    
    type(reactive_zone_c) :: react_zone
    type(gas_chemistry_c) :: gas_chem
    type(solid_chemistry_c) :: solid_chem
    type(gas_c) :: gas
    type(mineral_c) :: mineral
    type(mineral_zone_c) :: min_zone !> default
    type(aq_species_c) :: aq_species
    type(species_c) :: constrain
    
    !print *, associated(this%gas_chemistry)

    if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
        !call surf_chem%set_reactive_zone(this%reactive_zones(1)) !> chapuza, we assume that the first reactive zone is the one with surface chemistry
        !call this%set_solid_chemistry(surf_chem)
    else
        !call solid_chem%set_reactive_zone(react_zone) !> we set reactive zone by default
        !call solid_chem%set_mineral_zone(min_zone) !> we set mineral zone by default
        !call this%set_solid_chemistry(solid_chem) !> we set solid chemistry by default
        call this%solid_chemistry%reactive_zone%allocate_non_flowing_species(num_cstr)
    end if
    
    !if (this%solid_chemistry%reactive_zone%chem_syst%num_minerals>0) then
        call min_zone%set_chem_syst_min_zone(this%solid_chemistry%reactive_zone%chem_syst)
        call this%solid_chemistry%set_mineral_zone(min_zone)
    !end if
    
    allocate(ind_mins(num_cstr)) !> chapuza
    
    !    read(unit,*) label
    !    if (label=='INITIAL AND BOUNDARY WATER TYPES') then
    !        i=0 !> counter initial water types
    !        read(unit,*) model
    !        read(unit,*) this%num_init_wat_types, this%num_bd_wat_types, this%num_rech_wat_types
        if (n_p_aq<0 .or. n_p_aq>this%aq_phase%num_species) then
            error stop "Number of primary aqueous species not valid"
        else if (num_cstr<0 .or. num_cstr>n_p_aq) then
            error stop "Number of constrains not valid"
        else
            allocate(prim_indices(n_p_aq),icons(n_p_aq),ctots(n_p_aq),indices_constrains(num_cstr,2))
            !if (this%solid_chemistry%reactive_zone%chem_syst%gas_phase%num_species>0) then
            !    call react_zone%set_gas_phase(this%solid_chemistry%reactive_zone%chem_syst%gas_phase)
            !    call gas_chem%set_reactive_zone(react_zone)
            !    call gas_chem%allocate_partial_pressures()
            !    call gas_chem%allocate_conc_gases()
            !    call gas_chem%allocate_log_act_coeffs_gases()
            !    call gas_chem%set_temp(this%temp)
            !    call gas_chem%set_volume(1d0) !> arbitrary
            !    call this%set_gas_chemistry(gas_chem)
            !end if
        end if
        !call this%allocate_conc_comp(n_p_aq)
        !call this%allocate_log_act_coeffs()
        !call this%set_log_act_coeffs() 
        !call this%allocate_activities_aq_species()
        read(unit,*) label
        if (index(label,'guess')/=0) then !> 'icon, guess, ctot, constrain'
            k=1 !> counter primary aqueous species
            l=1 !> counter constrains
            m=1 !> counter gas constrains
            j=1 !> counter minerals in equilibrium
            allocate(n_icons(4))
            n_icons=0 !> number of each icon option
            do 
                read(unit,*) aq_species%name, icon, guess, ctot, constrain%name
                if (aq_species%name=='*') then
                    exit
                else
                    !call this%aq_phase%is_species_in_aq_phase(aq_species,flag,aq_sp_ind)
                    !if (flag.eqv..true.) then
                    !prim_indices(k)=this%indices_aq_species(k)
                        !> Chapuza
                        if (icon==1) then
                            n_icons(1)=n_icons(1)+1
                            if (aq_species%name=='h2o(p)') then
                                guess=1d0/18d-3
                                ctot=guess !> we impose water concentration
                            end if
                        else if (icon==2) then
                            n_icons(2)=n_icons(2)+1
                            n_aq_comp=n_aq_comp+1
                        else if (icon==3) then
                            n_icons(3)=n_icons(3)+1
                            if (aq_species%name=='h+') then
                                !call this%aq_phase%set_ind_prot(this%indices_aq_species(k))
                                call this%set_pH(-log10(ctot))
                            end if
                        else if (icon==4) then
                            n_icons(4)=n_icons(4)+1
                            call this%solid_chemistry%reactive_zone%chem_syst%is_eq_reaction_in_chem_syst(constrain%name,flag,&
                                ind_cstr)
                            if (flag.eqv..true.) then
                                indices_constrains(l,1)=k !> index of species in water type
                                indices_constrains(l,2)=ind_cstr !> index of eq reaction in chem syst
                                !ind_cstr=ind_cstr-this%solid_chemistry%reactive_zone%chem_syst%num_redox_eq_reacts-&
                                !this%solid_chemistry%reactive_zone%chem_syst%aq_phase%num_aq_complexes !> chapuza
                                call this%solid_chemistry%reactive_zone%chem_syst%is_species_in_chem_syst(constrain,flag_sp,ind_sp)
                                if (flag_sp.eqv..true.) then
                                    ind_cstr=ind_cstr-this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts !> index of heterogeneous eq reaction in chem syst
                                    call THIS%solid_chemistry%reactive_zone%non_flowing_species(l)%assign_species(&
                                        this%solid_chemistry%reactive_zone%chem_syst%species(ind_sp))
                                    l=l+1
                                    if (ind_cstr<=this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq_var_act) then
                                        ind_mins(j)=ind_cstr !> chapuza
                                        j=j+1
                                        this%solid_chemistry%reactive_zone%num_minerals_var_act=&
                                            this%solid_chemistry%reactive_zone%num_minerals_var_act+1
                                        this%solid_chemistry%reactive_zone%num_minerals=&
                                            this%solid_chemistry%reactive_zone%num_minerals+1
                                        this%solid_chemistry%mineral_zone%num_minerals_eq_var_act=&
                                            this%solid_chemistry%mineral_zone%num_minerals_eq_var_act+1
                                        this%solid_chemistry%mineral_zone%num_minerals_eq=&
                                            this%solid_chemistry%mineral_zone%num_minerals_eq+1
                                        this%solid_chemistry%mineral_zone%num_minerals=&
                                            this%solid_chemistry%mineral_zone%num_minerals+1
                                    else if (ind_cstr<=this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq) then
                                        ind_mins(j)=ind_cstr !> chapuza repetida ademas
                                        j=j+1
                                        this%solid_chemistry%reactive_zone%num_minerals_cst_act=&
                                            this%solid_chemistry%reactive_zone%num_minerals_cst_act+1
                                        this%solid_chemistry%reactive_zone%num_minerals=&
                                            this%solid_chemistry%reactive_zone%num_minerals+1
                                        this%solid_chemistry%mineral_zone%num_minerals_eq_cst_act=&
                                            this%solid_chemistry%mineral_zone%num_minerals_eq_cst_act+1
                                        this%solid_chemistry%mineral_zone%num_minerals_eq=&
                                            this%solid_chemistry%mineral_zone%num_minerals_eq+1
                                        this%solid_chemistry%mineral_zone%num_minerals=&
                                            this%solid_chemistry%mineral_zone%num_minerals+1
                                    else if (ind_cstr<=this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq+&
                                        this%solid_chemistry%reactive_zone%chem_syst%cat_exch%num_exch_cats) then
                                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats=&
                                            this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats+1
                                    else if (ind_cstr<=this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq+&
                                        this%solid_chemistry%reactive_zone%chem_syst%cat_exch%num_exch_cats+&
                                        this%solid_chemistry%reactive_zone%chem_syst%gas_phase%num_gases_eq_var_act) then
                                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act=&
                                            this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act+1
                                        this%solid_chemistry%reactive_zone%gas_phase%num_var_act_species=&
                                            this%solid_chemistry%reactive_zone%gas_phase%num_var_act_species+1
                                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq=&
                                            this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq+1
                                    else
                                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act=&
                                            this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1
                                        this%solid_chemistry%reactive_zone%gas_phase%num_cst_act_species=&
                                            this%solid_chemistry%reactive_zone%gas_phase%num_cst_act_species+1
                                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq=&
                                            this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq+1
                                    end if
                                else
                                    error stop "Constrain not found in chemical system"
                                end if
                                        
                                !call constrain%is_gas(flag_gas)
                                !if (flag_gas==.true.) then
                                !    call gas%set_name(constrain%name)
                                !    call this%gas_chemistry%reactive_zone%gas_phase%is_gas_in_gas_phase(gas,flag,gas_ind)
                                !    this%gas_chemistry%activities(gas_ind)=ctot
                                !else if (ind_cstr<=this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq) then
                                !    call THIS%solid_chemistry%reactive_zone%non_flowing_species(ind_cstr)%assign_species(this%solid_chemistry%reactive_zone%chem_syst%minerals(this%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin+ind_cstr)%mineral)
                                !else
                                !    call THIS%solid_chemistry%reactive_zone%non_flowing_species(ind_cstr-this%solid_chemistry%reactive_zone%chem_syst%num_minerals_eq-this%solid_chemistry%reactive_zone%chem_syst%gas_phase%num_gases_eq_cst_act)%assign_species(this%solid_chemistry%reactive_zone%chem_syst%cat_exch%surf_compl(this%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin+ind_cstr)%mineral)
                                !end if
                                !if (ind_cstr>this%solid_chemistry%reactive_zone%chem_syst%num_minerals_cst_act .AND. ind_cstr<this%solid_chemistry%reactive_zone%chem_syst%num_redox_eq_reacts) then
                                !else if (ind_cstr>this%solid_chemistry%reactive_zone%chem_syst%num_eq_reacts .AND. ind_cstr<this%solid_chemistry%reactive_zone%chem_syst%num_eq_reacts) then
                                !end if
                                !call this%solid_chemistry%reactive_zone%chem_syst%gas_phase%is_gas_in_gas_phase(gas,flag,gas_ind)
                            else
                                error stop "Constrain not found in chemical system"
                            end if
                        !else if (icon==5) then
                        !    n_icons(6)=n_icons(6)+1
                        else
                            error stop "icon option not implemented yet"
                        end if
                        !this%concentrations(this%indices_aq_species(k))=guess
                        icons(k)=icon
                        !ctots(this%indices_aq_species(k))=ctot
                        this%concentrations(k)=guess
                        !icons(k)=icon
                        ctots(k)=ctot
                        k=k+1 !> aqui hay que verificar dimension
                    !else
                    !    error stop 
                    !end if
                end if
            end do
            call this%solid_chemistry%reactive_zone%cat_exch_zone%compute_num_surf_compl()
            call this%solid_chemistry%reactive_zone%gas_phase%compute_num_species_phase()
            call this%solid_chemistry%reactive_zone%rearrange_non_flowing_species()
            call this%solid_chemistry%reactive_zone%set_num_solids()
            call this%solid_chemistry%reactive_zone%allocate_ind_mins()
            call this%solid_chemistry%reactive_zone%set_ind_mins_chem_syst(&
                ind_mins(1:this%solid_chemistry%reactive_zone%num_minerals))
            call this%solid_chemistry%compute_num_solids_solid_chem()
            if (associated(this%solid_chemistry%mineral_zone)) then
                call this%solid_chemistry%mineral_zone%set_num_mins_kin_min_zone()
            end if
            if (this%solid_chemistry%reactive_zone%gas_phase%num_species==0) then
                nullify(this%gas_chemistry)
            end if
        end if
    do i=1,this%aq_phase%num_species
        call this%aq_phase%aq_species(i)%params_act_coeff%compute_csts(this%aq_phase%aq_species(i)%valence,this%params_aq_sol,model)
    end do
!> We set speciation algebra
    call this%solid_chemistry%reactive_zone%speciation_alg%set_flag_comp(.false.) !> by default
    if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
        flag_surf=.true.
    else
        flag_surf=.false.
        call this%solid_chemistry%allocate_conc_solids()
        call this%solid_chemistry%allocate_activities()
        call this%solid_chemistry%allocate_log_act_coeffs_solid_chem()
    end if
    call this%solid_chemistry%reactive_zone%speciation_alg%set_flag_cat_exch(flag_surf)
    call this%solid_chemistry%reactive_zone%set_speciation_alg_dimensions()
    call this%set_ind_species()
    call this%solid_chemistry%reactive_zone%set_ind_eq_reacts()
    call this%solid_chemistry%reactive_zone%set_stoich_mat_react_zone()
    call this%solid_chemistry%reactive_zone%set_ind_mins_stoich_mat()
    call this%solid_chemistry%reactive_zone%set_ind_gases_stoich_mat()
    allocate(swap(2))
    call this%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(flag_Se,swap)
    if (flag_Se .eqv. .true.) then
        !> We swap indices constrains (chapuza)
        allocate(ind_swap(2))
        aux_ind_cstr=indices_constrains
        !indices_constrains(:,2)=indices_constrains(:,2)-this%solid_chemistry%reactive_zone%chem_syst%num_redox_eq_reacts
        do i=1,num_cstr
            if (indices_constrains(i,2)==swap(1)) then
                ind_swap(1)=i
            else if (indices_constrains(i,2)==swap(2)) then
                ind_swap(2)=i
            else
                continue
            end if
        end do
        indices_constrains(ind_swap(1),:)=aux_ind_cstr(ind_swap(2),:)
        indices_constrains(ind_swap(2),:)=aux_ind_cstr(ind_swap(1),:)
    end if
!> We compute initial concentrations (in molarities)
    !> aqui podrias usar polimorifsmo
    !print *, n_icons([2,4])
    if (sum(n_icons([2,4]))>0) then
        if (flag_surf .eqv. .true.) then
            if (Jac_opt==1) then
                call this%initialise_conc_anal_exch(icons,n_icons,indices_constrains,ctots,niter,CV_flag)
            else
                error stop "Initialisation subroutine not implemented yet"
            end if
        else if (model==0) then
            call this%initialise_conc_anal_ideal(icons,n_icons,indices_constrains,ctots,niter,CV_flag)
        else
            if (Jac_opt==0) then
                call this%initialise_conc_incr_coeff(icons,n_icons,indices_constrains,ctots,niter,CV_flag)
            else if (Jac_opt==1) then
                call this%initialise_conc_anal(icons,n_icons,indices_constrains,ctots,niter,CV_flag)
            else
                error stop "Jac_opt not valid"
            end if
        end if
    else if (model==0) then
        if (flag_surf .eqv. .true.) then
        !> Primary concentrations
            c1_aq=this%get_c1_aq()
            allocate(act_ads_cats_ig(this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats))
            act_ads_cats_ig(1)=1d-6 !> chapuza
            act_ads_cats_ig(2)=1d0-1d-9-act_ads_cats_ig(1) !> chapuza
            !call this%set_conc_sec_species(c2_ig)
            if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats==2) then
                call this%solid_chemistry%compute_conc_surf_ideal(c1_aq(&
                this%solid_chemistry%reactive_zone%cat_exch_zone%exch_cat_indices))
            else
                call this%solid_chemistry%compute_conc_surf_ideal_Newton(c1_aq(&
                this%solid_chemistry%reactive_zone%cat_exch_zone%exch_cat_indices),act_ads_cats_ig,niter,CV_flag)
                !call this%solid_chemistry%compute_conc_surf_ideal_Picard(c1_aq(this%solid_chemistry%reactive_zone%cat_exch_zone%exch_cat_indices),act_ads_cats_ig,niter,CV_flag)
            end if
        else if (this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions>0) then
            call this%set_act_aq_species()
            call this%compute_log_act_coeff_wat()
            a1=this%get_a1()
            log_gamma2=this%get_log_gamma2()
            allocate(c2_init(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
            call this%compute_c2_from_a1_ideal(a1,log_gamma2,c2_init)
        else
            call this%set_act_aq_species()
            call this%compute_log_act_coeff_wat()
            call this%compute_salinity()
            !call this%compute_molarities()
        end if
    else
        c1=this%get_c1()
        allocate(c2_ig(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
        c2_ig=1d-16 !> chapuza
        call this%set_conc_sec_species(c2_ig)
        call this%compute_c2_from_c1_Picard(c1,c2_ig,c2_init,niter,CV_flag)
    end if
end subroutine