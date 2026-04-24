    !> Reads initial target solids and their associated solid zones from a file
    !> We assume file has not already been opened
    subroutine read_tar_sol(this,root,nsrz,ngrz)
    use chemistry_m, only: chemistry_c
    use solid_chemistry_m, only: solid_chemistry_c
    use reactive_zone_m, only: reactive_zone_c
    use mineral_zone_m, only: mineral_zone_c
    implicit none
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: root !> root name of the file
    !type(aqueous_chemistry_c), intent(in) :: water_types(:) !> water types
    !type(solid_chemistry_c), intent(in) :: init_sol_types(:) !> initial solid zones
    !type(gas_chemistry_c), intent(in) :: init_gas_types(:) !> initial gas zones
    integer(kind=4), intent(in) :: nsrz !> number of solid reactive zones
    integer(kind=4), intent(in) :: ngrz !> number of gas reactive zones
    !integer(kind=4), intent(out) :: niter !> number of iterations
    !logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
    
    integer(kind=4) :: first,last,ind_dom,num_tar_sol,ind,i,j,k,m,iszone,num_tar_wat,tar_wat_ind,wtype,nstype,nbwtype,bwtype,&
        mix_wat_ind,ngzns,igzn,num_tar_wat_rech,ind_bar,flag_wat_type,num_tar_wat_bd,ind_bd,unit,ind_tar_sol
    real(kind=8), allocatable :: c_nc(:),u_init(:,:),c1_init(:),c2_init(:),c2_ig(:),gamma_2aq(:)
    integer(kind=4) :: aux_iszone, aux_igzn
    character(len=256) :: label,str,interval
    character(len=:), allocatable :: interval_trim,first_str,last_str
    logical :: flag_comp,flag_surf,flag_aq_phase,flag, flag_Se
    integer(kind=4), allocatable :: swap(:),aux_swap(:)
    type(solid_chemistry_c), target :: solid_chem !> default solid chemistry object
    type(reactive_zone_c), target :: react_zone   !> default reactive zone object
    type(reactive_zone_c), allocatable :: aux_react_zones(:) !> auxiliary reactive zones (chapuza)
    type(mineral_zone_c), target :: min_zone   !> default mineral zone object
    
    !nwtype=this%num_wat_types !> number of water types
    if (nsrz>0) then
        nstype=this%num_mineral_zones*(this%num_init_cat_exch_zones+1)+this%num_init_cat_exch_zones !> number of solid types
    else
        nstype=0 !> chapuza
    end if
    ! if (ngrz>0) then
    !     ngzns=size(init_gas_types)
    ! else
    !     ngzns=0 !> chapuza
    ! end if
    
    !flag_comp=.true. !> by default

    
    if (nstype==0) then !> chapuza
        allocate(aux_react_zones(this%num_reactive_zones))
        do i=1,this%num_reactive_zones
            call aux_react_zones(i)%copy_react_zone(this%reactive_zones(i))
        end do
        if (allocated(this%reactive_zones)) deallocate(this%reactive_zones)
        call this%allocate_reactive_zones(this%num_reactive_zones+1)
        do i=1,size(aux_react_zones)
            call this%reactive_zones(i)%copy_react_zone(aux_react_zones(i))
        end do
        call this%reactive_zones(this%num_reactive_zones)%set_chem_syst_react_zone(this%chem_syst)
        call this%reactive_zones(this%num_reactive_zones)%set_CV_params(this%CV_params)
        deallocate(aux_react_zones)
    end if
        
    call this%allocate_target_solids_dummy() !> we allocate dummy target solids (no associated solids)
    call this%allocate_reactive_zones_dummy() !> we allocate dummy reactive zones (no associated reactions)

    call this%min_zone_dummy%set_chem_syst_min_zone(this%chem_syst) !> we set the chemical system in dummy mineral zone
    
    !call solid_chem%set_reactive_zone(this%reactive_zones(this%num_reactive_zones)) !> we set default reactive zone in default solid chemistry object
    !call solid_chem%set_mineral_zone(this%mineral_zones(1)) !> we set default mineral zone in default solid chemistry object
    
    
    do i=1,this%num_target_solids_dummy !> we set dummy target solids attributes
        !call this%target_solids_dummy(i)%set_chem_syst(this%chem_syst)
        !call this%target_solids_dummy(i)%set_CV_params(this%CV_params)
        call this%target_solids_dummy(i)%set_id(0) !> we set default id
        call this%target_solids_dummy(i)%set_mineral_zone(this%min_zone_dummy) !> we set the dummy mineral zone pointer
        call this%react_zones_dummy(i)%set_CV_params(this%CV_params) !> we set the CV parameters in dummy mineral zone
        call this%react_zones_dummy(i)%set_chem_syst_react_zone(this%chem_syst) !> we set the chemical system in dummy reactive zone
        call this%target_solids_dummy(i)%set_reactive_zone(this%react_zones_dummy(i)) !> we set the default reactive zone pointer
        !call this%target_solids_dummy(i)%reactive_zone%set_speciation_alg_dimensions(.true.) !< esto creo que no es necesario
        !call this%target_solids_dummy(i)%reactive_zone%set_ind_eq_reacts() !> chapuza
        !call this%target_solids_dummy(i)%reactive_zone%set_stoich_mat_react_zone() !> chapuza
        !call this%target_solids_dummy(i)%reactive_zone%set_ind_gases_stoich_mat() !> chapuza
        !call this%target_solids_dummy(i)%reactive_zone%set_ind_mins_stoich_mat() !> chapuza
        !call this%target_solids_dummy(i)%allocate_reaction_rates_solid_chem() !> we allocate reaction rates in solid chemistry object
    end do

    ! if (this%num_mineral_zones==0) then !> chapuza
    !     do i=1,this%num_materials
    !         call this%materials(i)%set_mineral_zone(this%min_zone_dummy) !> we set the dummy mineral zone in each adsorption zone (chapuza)
    !     end do
    ! end if

    unit=56 !> arbitrary unit number
    open(unit, file=root//'_tar_sol.dat', status='old', action='read') !> we open the file with target waters associated with their solids and/or gases
    do
        read(unit,*) label
        if (label=='end') then
            exit
        else if (label=='TARGET SOLIDS') then
            read(unit,*) num_tar_sol
            call this%allocate_target_solids(num_tar_sol)
            do
                read(unit,*) interval, iszone !> interval of target solids and solid zone index
                interval_trim=trim(interval)
                ind_bar=index(interval_trim,'-')
                if (ind_bar>0) then !> we have a range of target solids
                    first_str=interval_trim(1:ind_bar-1) !> first target of interval
                    last_str=interval_trim(ind_bar+1:) !> last target of interval
                    read(first_str,*) first
                    read(last_str,*) last
                    if (first<1 .or. first>last .or. last>num_tar_sol) then
                        error stop "Target solid index out of bounds"
                    else
                        do i=first,last
                            call this%target_solids(i)%copy_solid_chemistry(this%materials(iszone)) !> we assign the initial solid type
                            call this%target_solids(i)%set_id(iszone) !> the id of the target solid is the solid zone index
                            !call this%target_solids(i)%set_reactive_zone(this%reactive_zones(ngrz+iszone)) !> we set the reactive zone of the target solid
                            call this%target_solids(i)%reactive_zone%set_speciation_alg_dimensions(.true.) !< esto creo que no es necesario
                            call this%target_solids(i)%reactive_zone%set_ind_eq_reacts() !> chapuza
                            call this%target_solids(i)%reactive_zone%allocate_ind_var_act_species() !> chapuza
                            call this%target_solids(i)%reactive_zone%set_stoich_mat_react_zone() !> chapuza
                            call this%target_solids(i)%reactive_zone%set_ind_gases_stoich_mat() !> chapuza
                            call this%target_solids(i)%reactive_zone%set_ind_mins_stoich_mat() !> chapuza
                            call this%target_solids(i)%allocate_reaction_rates_solid_chem() !> we allocate reaction rates in solid chemistry object
                        end do
                    end if
                    if (last==num_tar_sol) then
                        exit
                    end if
                else
                    !> we have a single target solid
                    read(interval_trim,*) ind_tar_sol
                    if (ind_tar_sol<1 .or. ind_tar_sol>num_tar_sol) then
                        error stop "Target solid index out of bounds"
                    else
                        call this%target_solids(ind_tar_sol)%copy_solid_chemistry(this%materials(iszone)) !> we set the initial solid type
                        call this%target_solids(ind_tar_sol)%set_id(iszone) !> the id of the target solid is the material index
                        !call this%target_solids(ind_tar_sol)%set_reactive_zone(this%reactive_zones(ngrz+iszone)) !> we set the reactive zone of the target solid
                        call this%target_solids(ind_tar_sol)%reactive_zone%set_speciation_alg_dimensions(.true.) !< esto creo que no es necesario
                        call this%target_solids(ind_tar_sol)%reactive_zone%set_ind_eq_reacts() !> chapuza
                        call this%target_solids(ind_tar_sol)%reactive_zone%allocate_ind_var_act_species() !> chapuza
                        call this%target_solids(ind_tar_sol)%reactive_zone%set_stoich_mat_react_zone() !> chapuza
                        call this%target_solids(ind_tar_sol)%reactive_zone%set_ind_gases_stoich_mat() !> chapuza
                        call this%target_solids(ind_tar_sol)%reactive_zone%set_ind_mins_stoich_mat() !> chapuza
                    end if
                    if (ind_tar_sol==num_tar_sol) then
                        exit
                    end if
                end if
            end do
            !if (nstype>0) then
                !call this%allocate_target_solids(this%num_target_waters) !> we assume bijection with target waters (chapuza)
            !end if
            !if (ngzns>0) then
                !call this%allocate_target_gases(this%num_target_waters) !> we assume bijection with target waters (chapuza)
            !end if
            ! read(unit,*) num_tar_wat_rech
            ! call this%allocate_rech_waters_indices(num_tar_wat_rech)            
            ! read(unit,*) num_tar_wat_bd
            ! call this%allocate_bd_waters_indices(num_tar_wat_bd)
            ! call this%allocate_tar_wat_indices(this%num_target_waters-this%num_rech_waters-this%num_bd_waters)
            !aux_iszone=0
            !aux_igzn=0
            ! if (num_tar_wat_rech>0 .or. num_tar_wat_bd>0) then
            !     do i=1,this%num_target_waters
            !         read(unit,*) tar_wat_ind, wtype, public, iszone, igzn, flag_wat_type
            !         if (tar_wat_ind<1 .or. tar_wat_ind>this%num_target_waters) then
            !             error stop "Target water index out of bounds"
            !         else if (wtype<1 .or. wtype>nwtype) then
            !             error stop "Water type index out of bounds"
            !         else if (iszone<0 .or. iszone>nstype) then
            !             error stop "Solid type index out of bounds"
            !         else if (igzn<0 .or. igzn>ngzns) then
            !             error stop "Gas type index out of bounds"
            !         else if (flag_wat_type==0) then !> boundary water
            !             ind_bd=ind_bd+1 !> counter boundary waters
            !             this%bd_waters_indices(ind_bd)=tar_wat_ind
            !         else if (flag_wat_type==2) then !> external water
            !             ind_ext=ind_ext+1 !> counter external waters
            !             this%ext_waters_indices(ind_ext)=tar_wat_ind
            !         else if (flag_wat_type==1) then !> domain water
            !             ind_dom=ind_dom+1 !> counter domain waters
            !             this%tar_wat_indices(ind_dom)=tar_wat_ind
            !         else
            !             error stop "Water type flag out of bounds"
            !         end if
            !         call this%loop_read_tar_wat_init(flag,this%wat_types,init_sol_types,init_gas_types,nsrz,ngrz,tar_wat_ind,wtype, public,&
            !             iszone,igzn,aux_iszone,aux_igzn,solid_chem)
                    !allocate(rk(this%waters(tar_wat_ind)%indices_rk%num_cols)) !> chapuza
                    !call this%waters(tar_wat_ind)%compute_rk(rk)
                !    !aux_iszone=iszone
                !    this%waters(tar_wat_ind)=water_types(wtype)
                !    ! if (counter_swap==0) then
                !    !     call this%waters(tar_wat_ind)%set_aq_phase(this%chem_syst%aq_phase)
                !    !     !call this%waters(tar_wat_ind)%set_indices_aq_species_aq_chem()
                !    ! end if
                !    if (iszone>0) then
                !        this%target_solids(tar_wat_ind)=init_sol_types(iszone)
                !        !> chapuza intercambio
                !        if (this%reactive_zones(ngzns+nstype*igzn+1)%cat_exch_zone%num_surf_compl>0) then
                !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+1))
                !        else
                !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+iszone))
                !        end if
                !        if (igzn>0) then
                !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+iszone))
                !            this%target_gases(tar_wat_ind)=init_gas_types(igzn)
                !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+iszone))
                !            call this%waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
                !        else
                !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+iszone))
                !        end if
                !        call this%waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(tar_wat_ind))
                !    else
                !        if (igzn>0) then
                !            this%target_gases(tar_wat_ind)=init_gas_types(igzn)
                !            ! call this%reactive_zones(igzn)%set_ind_eq_reacts() !> chapuza
                !            ! call this%reactive_zones(igzn)%set_stoich_mat_react_zone() !> chapuza
                !            ! call this%reactive_zones(igzn)%set_ind_gases_stoich_mat() !> chapuza
                !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(igzn))
                !            call this%waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
                !            call solid_chem%set_reactive_zone(this%reactive_zones(igzn))
                !        else
                !            ! call react_zone%set_ind_eq_reacts() !> chapuza
                !            ! call react_zone%set_stoich_mat_react_zone() !> chapuza
                !            ! call react_zone%set_ind_gases_stoich_mat() !> chapuza
                !            ! call react_zone%set_ind_mins_stoich_mat() !> chapuza
                !            call solid_chem%set_reactive_zone(react_zone)
                !        end if
                !        call solid_chem%set_mineral_zone(min_zone)
                !        call this%waters(tar_wat_ind)%set_solid_chemistry(solid_chem)
                !    end if
                !    !> we check if there is a new reactive zone
                !    if (aux_iszone==0 .or. aux_iszone/=iszone .or. aux_igzn/=igzn) then !> we assume target waters are grouped by their reactive zones
                !        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_speciation_alg_dimensions(flag_comp)
                !        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_ind_eq_reacts() !> chapuza
                !        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_stoich_mat_react_zone() !> chapuza
                !        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_ind_gases_stoich_mat() !> chapuza
                !        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_ind_mins_stoich_mat() !> chapuza
                !        call this%waters(tar_wat_ind)%set_ind_species()
                !        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
                !            flag_Se,swap)
                !        if (flag_Se.eqv..true.) then !> we swap indices of species
                !            aux_swap(1)=this%waters(tar_wat_ind)%ind_var_act_species(swap(1))
                !            aux_swap(2)=this%waters(tar_wat_ind)%ind_var_act_species(swap(2))
                !            !     this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
                !            !print *, this%waters(tar_wat_ind)%ind_var_act_species
                !            !print *, this%waters(tar_wat_ind)%ind_sec_species
                !            this%waters(tar_wat_ind)%ind_var_act_species(swap(1))=aux_swap(2)
                !            this%waters(tar_wat_ind)%ind_var_act_species(swap(2))=aux_swap(1)
                !            ! this%waters(tar_wat_ind)%ind_prim_species(swap(1))=aux_swap(2) !> index of primary species
                !            ! this%waters(tar_wat_ind)%ind_sec_species(swap(2)-&
                !            !     this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=&
                !            !     aux_swap(1) !> index of secondary species
                !        end if
                !    else if (aux_iszone>0 .or. aux_igzn>0) then !> indices remain the same because reactive zone is the same
                !        this%waters(tar_wat_ind)%ind_var_act_species=this%waters(tar_wat_ind-1)%ind_var_act_species
                !        !this%waters(tar_wat_ind)%ind_sec_species=this%waters(tar_wat_ind-1)%ind_sec_species
                !    end if
                !    print *, this%waters(tar_wat_ind)%ind_var_act_species
                !!> Chapuza
                !    !if (associated(this%waters(tar_wat_ind)%solid_chemistry%mineral_zone)) then
                !        if (this%waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin<&
                !            this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin) then
                !            call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod(&
                !                this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                !                this%waters(tar_wat_ind)%solid_chemistry%mineral_zone%ind_min_chem_syst(1:&
                !                this%waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin))
                !        else
                !            call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod()
                !        end if
                !    !end if
                !    call this%waters(tar_wat_ind)%allocate_reaction_rates()
                !    call this%waters(tar_wat_ind)%set_indices_rk()
                !    aux_iszone=iszone
                !    aux_igzn=igzn
    !             end do
    !         else
    !             do i=1,this%num_target_waters
    !                 read(unit,*) tar_wat_ind, wtype, public, iszone, igzn
    !                 if (tar_wat_ind<1 .or. tar_wat_ind>this%num_target_waters) then
    !                     error stop "Target water index out of bounds"
    !                 else if (wtype<1 .or. wtype>nwtype) then
    !                     error stop "Water type index out of bounds"
    !                 else if (iszone<0 .or. iszone>nstype) then
    !                     error stop "Solid type index out of bounds"
    !                 else if (igzn<0 .or. igzn>ngzns) then
    !                     error stop "Gas type index out of bounds"
    !                 else
    !                     ind_dom=ind_dom+1
    !                     this%tar_wat_indices(ind_dom)=tar_wat_ind
    !                 end if
    !                 call this%loop_read_tar_wat_init(.false.,this%wat_types,init_sol_types,init_gas_types,nsrz,ngrz,tar_wat_ind,&
    !                     wtype, public,iszone,igzn,aux_iszone,aux_igzn,solid_chem)
    !             !    !aux_iszone=iszone
    !             !    this%waters(tar_wat_ind)=water_types(wtype)
    !             !    ! if (counter_swap==0) then
    !             !    !     call this%waters(tar_wat_ind)%set_aq_phase(this%chem_syst%aq_phase)
    !             !    !     !call this%waters(tar_wat_ind)%set_indices_aq_species_aq_chem()
    !             !    ! end if
    !             !    if (iszone>0) then
    !             !        this%target_solids(tar_wat_ind)=init_sol_types(iszone)
    !             !        !> chapuza intercambio
    !             !        if (this%reactive_zones(ngzns+nstype*igzn+1)%cat_exch_zone%num_surf_compl>0) then
    !             !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+1))
    !             !        else
    !             !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+iszone))
    !             !        end if
    !             !        if (igzn>0) then
    !             !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+iszone))
    !             !            this%target_gases(tar_wat_ind)=init_gas_types(igzn)
    !             !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+iszone))
    !             !            call this%waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
    !             !        else
    !             !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+iszone))
    !             !        end if
    !             !        call this%waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(tar_wat_ind))
    !             !    else if (igzn>0) then
    !             !        this%target_gases(tar_wat_ind)=init_gas_types(igzn)
    !             !        call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(igzn))
    !             !        call this%waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
    !             !        call solid_chem%set_reactive_zone(this%reactive_zones(igzn))
    !             !        call this%waters(tar_wat_ind)%set_solid_chemistry(solid_chem)
    !             !    else
    !             !        call solid_chem%set_reactive_zone(react_zone)
    !             !        call this%waters(tar_wat_ind)%set_solid_chemistry(solid_chem)
    !             !    end if
    !             !    !> we check if there is a new reactive zone
    !             !    if (aux_iszone==0 .or. aux_iszone/=iszone .or. aux_igzn/=igzn) then !> we assume target waters are grouped by their reactive zones
    !             !        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_speciation_alg_dimensions(flag_comp)
    !             !        call this%waters(tar_wat_ind)%set_ind_species()
    !             !        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
    !             !            flag_Se,swap)
    !             !        if (flag_Se.eqv..true.) then !> we swap indices of species
    !             !            aux_swap(1)=this%waters(tar_wat_ind)%ind_var_act_species(swap(1))
    !             !            aux_swap(2)=this%waters(tar_wat_ind)%ind_var_act_species(swap(2))
    !             !            !     this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
    !             !            !print *, this%waters(tar_wat_ind)%ind_var_act_species
    !             !            !print *, this%waters(tar_wat_ind)%ind_sec_species
    !             !            this%waters(tar_wat_ind)%ind_var_act_species(swap(1))=aux_swap(2)
    !             !            this%waters(tar_wat_ind)%ind_var_act_species(swap(2))=aux_swap(1)
    !             !            ! this%waters(tar_wat_ind)%ind_prim_species(swap(1))=aux_swap(2) !> index of primary species
    !             !            ! this%waters(tar_wat_ind)%ind_sec_species(swap(2)-&
    !             !            !     this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=&
    !             !            !     aux_swap(1) !> index of secondary species
    !             !        end if
    !             !    else if (aux_iszone>0 .or. aux_igzn>0) then !> indices remain the same because reactive zone is the same
    !             !        this%waters(tar_wat_ind)%ind_var_act_species=this%waters(tar_wat_ind-1)%ind_var_act_species
    !             !        !this%waters(tar_wat_ind)%ind_sec_species=this%waters(tar_wat_ind-1)%ind_sec_species
    !             !    end if
    !             !    ! if (flag_Se.eqv..true.) then
    !             !    !     counter_swap=counter_swap+1
    !             !    !     aux_swap=swap
    !             !    !     !print *, this%waters(tar_wat_ind)%indices_aq_species
    !             !    !     do j=tar_wat_ind+1,this%num_target_waters
    !             !    !         call this%waters(j)%set_aq_phase(this%chem_syst%aq_phase)
    !             !    !         call this%waters(j)%solid_chemistry%reactive_zone%set_speciation_alg_dimensions(&
    !             !    !             flag_comp)
    !             !    !         call this%waters(j)%set_ind_prim_sec_species()
    !             !    !         this%waters(j)%ind_prim_species(swap(1))=swap(2)
    !             !    !         this%waters(j)%ind_sec_species(swap(2)-&
    !             !    !             this%waters(j)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=swap(1) !> chapuza
    !             !    !         !call this%waters(j)%set_indices_aq_species_aq_chem()
    !             !    !         ! this%waters(j)%indices_aq_species(swap(2))=swap(1)
    !             !    !         ! this%waters(j)%indices_aq_species(swap(1))=swap(2)
    !             !    !     end do
    !             !    ! end if
    !             !    !print *, this%waters(tar_wat_ind)%ind_var_act_species
    !             !    !print *, this%waters(tar_wat_ind)%ind_sec_species
    !             !    ! if (aux_swap(1)>0 .AND. aux_swap(2)>0) then
    !             !    !     !call this%waters(tar_wat_ind)%set_indices_aq_species_aq_chem()
    !             !    !     this%waters(tar_wat_ind)%ind_prim_species(aux_swap(1))=aux_swap(2)
    !             !    !     this%waters(tar_wat_ind)%ind_sec_species(swap(2)-&
    !             !    !         this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=swap(1) !> chapuza
    !             !    ! end if
    !             !    if (this%waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin<&
    !             !        this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin) then
    !             !        call this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod(&
    !             !            this%waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
    !             !            this%waters(tar_wat_ind)%solid_chemistry%mineral_zone%ind_min_chem_syst(1:&
    !             !            this%waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin))
    !             !    else
    !             !        call this%target_solids%reactive_zone%compute_U_SkT_prod()
    !             !    end if
    !             !    call this%waters(tar_wat_ind)%allocate_reaction_rates()
    !             !    call this%waters(tar_wat_ind)%set_indices_rk()
    !             !    aux_iszone=iszone
    !             !    aux_igzn=igzn
    !             end do
    !         end if
        else
            continue
        end if
    end do
    ! close(unit) !> we close the file with target waters associated with their solids and/or gases
    ! this%waters_init=this%waters
    ! this%target_solids_init=this%target_solids
    ! if (allocated(this%target_gases)) then
    !     this%target_gases_init=this%target_gases
    ! end if
    ! call this%write_aq_comps_init(root) !> write aqueous components of initial target waters
    do i=1,this%num_target_solids
        call this%target_solids_init(i)%copy_solid_chemistry(this%target_solids(i)) !> we assign the initial target solids to the target solids
    end do
end subroutine