    !> Reads initial target gases from a file
    !> We assume file has not already been opened
    subroutine read_tar_gas(this,root,ngrz)
    use chemistry_m, only: chemistry_c, aqueous_chemistry_c, gas_chemistry_c, gas_chemistry_c, reactive_zone_c, &
        mineral_zone_c
    implicit none
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: root !> root name of the file
    !type(aqueous_chemistry_c), intent(in) :: water_types(:) !> water types
    !type(gas_chemistry_c), intent(in) :: init_sol_types(:) !> initial gas zones
    !type(gas_chemistry_c), intent(in) :: init_gas_types(:) !> initial gas zones
    !integer(kind=4), intent(in) :: nsrz !> number of gas reactive zones
    integer(kind=4), intent(in) :: ngrz !> number of gas reactive zones
    !integer(kind=4), intent(out) :: niter !> number of iterations
    !logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
    
    integer(kind=4) :: first,last,ind_dom,num_tar_gas,ind,i,j,k,m,igzone,num_tar_wat,tar_wat_ind,wtype,nstype,nbwtype,bwtype,&
        mix_wat_ind,ngzns,igzn,num_tar_wat_rech,ind_bar,flag_wat_type,num_tar_wat_bd,ind_bd,unit,ind_tar_gas
    real(kind=8), allocatable :: c_nc(:),u_init(:,:),c1_init(:),c2_init(:),c2_ig(:),gamma_2aq(:)
    integer(kind=4) :: aux_igzone, aux_igzn
    character(len=256) :: label,str,interval
    character(len=:), allocatable :: interval_trim,first_str,last_str
    logical :: flag_comp,flag_surf,flag_aq_phase,flag, flag_Se
    integer(kind=4), allocatable :: swap(:),aux_swap(:)
    type(gas_chemistry_c), target :: gas_chem !> default gas chemistry object
    type(reactive_zone_c), target :: react_zone   !> default reactive zone object
    type(reactive_zone_c), allocatable :: aux_react_zones(:) !> auxiliary reactive zones (chapuza)
    type(mineral_zone_c), target :: min_zone   !> default mineral zone object
    
    !nwtype=this%num_wat_types !> number of water types
    !if (ngrz>0) then
    !    nstype=this%num_mineral_zones !> number of gas zones
    !else
    !    nstype=0 !> chapuza
    !end if
    ! if (ngrz>0) then
    !     ngzns=size(init_gas_types)
    ! else
    !     ngzns=0 !> chapuza
    ! end if
    
    !flag_comp=.true. !> by default

    
    !if (nstype==0) then !> chapuza
    !    allocate(aux_react_zones(this%num_reactive_zones))
    !    do i=1,this%num_reactive_zones
    !        call aux_react_zones(i)%assign_react_zone(this%reactive_zones(i))
    !    end do
    !    if (allocated(this%reactive_zones)) deallocate(this%reactive_zones)
    !    call this%allocate_reactive_zones(this%num_reactive_zones+1)
    !    do i=1,size(aux_react_zones)
    !        call this%reactive_zones(i)%assign_react_zone(aux_react_zones(i))
    !    end do
    !    call this%reactive_zones(this%num_reactive_zones)%set_chem_syst_react_zone(this%chem_syst)
    !    call this%reactive_zones(this%num_reactive_zones)%set_CV_params(this%CV_params)
    !    deallocate(aux_react_zones)
    !    call this%allocate_mineral_zones(1) !> chapuza
    !    call this%mineral_zones(1)%set_chem_syst_min_zone(this%chem_syst) !> we set the chemical system in default mineral zone
    !    call gas_chem%set_reactive_zone(this%reactive_zones(this%num_reactive_zones)) !> we set default reactive zone in default gas chemistry object
    !    call gas_chem%set_mineral_zone(this%mineral_zones(1)) !> we set default mineral zone in default gas chemistry object
    !end if
    
    ! ind_Ext=0 !> counter recharge waters
    ! ind_bd=0 !> counter boundary waters
    ! ind_dom=0 !> counter domain waters

    unit=5 !> arbitrary unit number
    open(unit, file=root//'_tar_gas.dat', status='old', action='read') !> we open the file with target waters associated with their gases and/or gases
    do
        read(unit,*) label
        if (label=='end') then
            exit
        else if (label=='TARGET GASES') then
            read(unit,*) num_tar_gas
            call this%allocate_target_gases(num_tar_gas)
            do
                read(unit,*) interval, igzone !> interval of target gases and gas zone index
                interval_trim=trim(interval)
                ind_bar=index(interval_trim,'-')
                if (ind_bar>0) then !> we have a range of target gases
                    first_str=interval_trim(1:ind_bar-1) !> first target of interval
                    last_str=interval_trim(ind_bar+1:) !> last target of interval
                    read(first_str,*) first
                    read(last_str,*) last
                    if (first<1 .or. first>last .or. last>num_tar_gas) then
                        error stop "Target gas index out of bounds"
                    else
                        do i=first,last
                            this%target_gases(i)=this%gas_zones(igzone) !> we set the initial gas type
                            call this%target_gases(i)%set_reactive_zone(this%reactive_zones(igzone)) !> we set the reactive zone of the target gas
                            call this%target_gases(i)%reactive_zone%set_speciation_alg_dimensions(.true.) !< esto creo que no es necesario
                            call this%target_gases(i)%reactive_zone%set_ind_eq_reacts() !> chapuza
                            call this%target_gases(i)%reactive_zone%set_stoich_mat_react_zone() !> chapuza
                            call this%target_gases(i)%reactive_zone%set_ind_gases_stoich_mat() !> chapuza
                            call this%target_gases(i)%reactive_zone%set_ind_mins_stoich_mat() !> chapuza
                            call this%target_gases(i)%allocate_reaction_rates_gas_chem() !> we allocate reaction rates in gas chemistry object
                        end do
                    end if
                    if (last==num_tar_gas) then
                        exit
                    end if
                else
                    !> we have a single target gas
                    read(interval_trim,*) ind_tar_gas
                    if (ind_tar_gas<1 .or. ind_tar_gas>num_tar_gas) then
                        error stop "Target gas index out of bounds"
                    else
                        this%target_gases(ind_tar_gas)=this%gas_zones(igzone) !> we set the initial gas type
                        call this%target_gases(ind_tar_gas)%set_reactive_zone(this%reactive_zones(igzone)) !> we set the reactive zone of the target gas
                        call this%target_gases(ind_tar_gas)%reactive_zone%set_speciation_alg_dimensions(.true.) !< esto creo que no es necesario
                        call this%target_gases(ind_tar_gas)%reactive_zone%set_ind_eq_reacts() !> chapuza
                        call this%target_gases(ind_tar_gas)%reactive_zone%set_stoich_mat_react_zone() !> chapuza
                        call this%target_gases(ind_tar_gas)%reactive_zone%set_ind_gases_stoich_mat() !> chapuza
                        call this%target_gases(ind_tar_gas)%reactive_zone%set_ind_mins_stoich_mat() !> chapuza
                    end if
                    if (ind_tar_gas==num_tar_gas) then
                        exit
                    end if
                end if
            end do
            !if (nstype>0) then
                !call this%allocate_target_gases(this%num_target_waters) !> we assume bijection with target waters (chapuza)
            !end if
            !if (ngzns>0) then
                !call this%allocate_target_gases(this%num_target_waters) !> we assume bijection with target waters (chapuza)
            !end if
            ! read(unit,*) num_tar_wat_rech
            ! call this%allocate_rech_waters_indices(num_tar_wat_rech)            
            ! read(unit,*) num_tar_wat_bd
            ! call this%allocate_bd_waters_indices(num_tar_wat_bd)
            ! call this%allocate_dom_tar_wat_indices(this%num_target_waters-this%num_rech_waters-this%num_bd_waters)
            !aux_igzone=0
            !aux_igzn=0
            ! if (num_tar_wat_rech>0 .or. num_tar_wat_bd>0) then
            !     do i=1,this%num_target_waters
            !         read(unit,*) tar_wat_ind, wtype, igzone, igzn, flag_wat_type
            !         if (tar_wat_ind<1 .or. tar_wat_ind>this%num_target_waters) then
            !             error stop "Target water index out of bounds"
            !         else if (wtype<1 .or. wtype>nwtype) then
            !             error stop "Water type index out of bounds"
            !         else if (igzone<0 .or. igzone>nstype) then
            !             error stop "gas type index out of bounds"
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
            !             this%dom_tar_wat_indices(ind_dom)=tar_wat_ind
            !         else
            !             error stop "Water type flag out of bounds"
            !         end if
            !         call this%loop_read_tar_wat_init(flag,this%wat_types,init_sol_types,init_gas_types,nsrz,ngrz,tar_wat_ind,wtype,&
            !             igzone,igzn,aux_igzone,aux_igzn,gas_chem)
                    !allocate(rk(this%target_waters(tar_wat_ind)%indices_rk%num_cols)) !> chapuza
                    !call this%target_waters(tar_wat_ind)%compute_rk(rk)
                !    !aux_igzone=igzone
                !    this%target_waters(tar_wat_ind)=water_types(wtype)
                !    ! if (counter_swap==0) then
                !    !     call this%target_waters(tar_wat_ind)%set_aq_phase(this%chem_syst%aq_phase)
                !    !     !call this%target_waters(tar_wat_ind)%set_indices_aq_species_aq_chem()
                !    ! end if
                !    if (igzone>0) then
                !        this%target_gases(tar_wat_ind)=init_sol_types(igzone)
                !        !> chapuza intercambio
                !        if (this%reactive_zones(ngzns+nstype*igzn+1)%cat_exch_zone%num_surf_compl>0) then
                !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+1))
                !        else
                !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+igzone))
                !        end if
                !        if (igzn>0) then
                !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+igzone))
                !            this%target_gases(tar_wat_ind)=init_gas_types(igzn)
                !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+igzone))
                !            call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
                !        else
                !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+igzone))
                !        end if
                !        call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
                !    else
                !        if (igzn>0) then
                !            this%target_gases(tar_wat_ind)=init_gas_types(igzn)
                !            ! call this%reactive_zones(igzn)%set_ind_eq_reacts() !> chapuza
                !            ! call this%reactive_zones(igzn)%set_stoich_mat_react_zone() !> chapuza
                !            ! call this%reactive_zones(igzn)%set_ind_gases_stoich_mat() !> chapuza
                !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(igzn))
                !            call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
                !            call gas_chem%set_reactive_zone(this%reactive_zones(igzn))
                !        else
                !            ! call react_zone%set_ind_eq_reacts() !> chapuza
                !            ! call react_zone%set_stoich_mat_react_zone() !> chapuza
                !            ! call react_zone%set_ind_gases_stoich_mat() !> chapuza
                !            ! call react_zone%set_ind_mins_stoich_mat() !> chapuza
                !            call gas_chem%set_reactive_zone(react_zone)
                !        end if
                !        call gas_chem%set_mineral_zone(min_zone)
                !        call this%target_waters(tar_wat_ind)%set_gas_chemistry(gas_chem)
                !    end if
                !    !> we check if there is a new reactive zone
                !    if (aux_igzone==0 .or. aux_igzone/=igzone .or. aux_igzn/=igzn) then !> we assume target waters are grouped by their reactive zones
                !        call this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%set_speciation_alg_dimensions(flag_comp)
                !        call this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%set_ind_eq_reacts() !> chapuza
                !        call this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%set_stoich_mat_react_zone() !> chapuza
                !        call this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%set_ind_gases_stoich_mat() !> chapuza
                !        call this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%set_ind_mins_stoich_mat() !> chapuza
                !        call this%target_waters(tar_wat_ind)%set_ind_species()
                !        call this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%compute_speciation_alg_arrays(&
                !            flag_Se,swap)
                !        if (flag_Se.eqv..true.) then !> we swap indices of species
                !            aux_swap(1)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))
                !            aux_swap(2)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))
                !            !     this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%speciation_alg%num_prim_species)
                !            !print *, this%target_waters(tar_wat_ind)%ind_var_act_species
                !            !print *, this%target_waters(tar_wat_ind)%ind_sec_species
                !            this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))=aux_swap(2)
                !            this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))=aux_swap(1)
                !            ! this%target_waters(tar_wat_ind)%ind_prim_species(swap(1))=aux_swap(2) !> index of primary species
                !            ! this%target_waters(tar_wat_ind)%ind_sec_species(swap(2)-&
                !            !     this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%speciation_alg%num_prim_species)=&
                !            !     aux_swap(1) !> index of secondary species
                !        end if
                !    else if (aux_igzone>0 .or. aux_igzn>0) then !> indices remain the same because reactive zone is the same
                !        this%target_waters(tar_wat_ind)%ind_var_act_species=this%target_waters(tar_wat_ind-1)%ind_var_act_species
                !        !this%target_waters(tar_wat_ind)%ind_sec_species=this%target_waters(tar_wat_ind-1)%ind_sec_species
                !    end if
                !    print *, this%target_waters(tar_wat_ind)%ind_var_act_species
                !!> Chapuza
                !    !if (associated(this%target_waters(tar_wat_ind)%gas_chemistry%mineral_zone)) then
                !        if (this%target_waters(tar_wat_ind)%gas_chemistry%mineral_zone%num_minerals_kin<&
                !            this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%chem_syst%num_minerals_kin) then
                !            call this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%compute_U_SkT_prod(&
                !                this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                !                this%target_waters(tar_wat_ind)%gas_chemistry%mineral_zone%ind_min_chem_syst(1:&
                !                this%target_waters(tar_wat_ind)%gas_chemistry%mineral_zone%num_minerals_kin))
                !        else
                !            call this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%compute_U_SkT_prod()
                !        end if
                !    !end if
                !    call this%target_waters(tar_wat_ind)%allocate_reaction_rates()
                !    call this%target_waters(tar_wat_ind)%set_indices_rk()
                !    aux_igzone=igzone
                !    aux_igzn=igzn
    !             end do
    !         else
    !             do i=1,this%num_target_waters
    !                 read(unit,*) tar_wat_ind, wtype, igzone, igzn
    !                 if (tar_wat_ind<1 .or. tar_wat_ind>this%num_target_waters) then
    !                     error stop "Target water index out of bounds"
    !                 else if (wtype<1 .or. wtype>nwtype) then
    !                     error stop "Water type index out of bounds"
    !                 else if (igzone<0 .or. igzone>nstype) then
    !                     error stop "gas type index out of bounds"
    !                 else if (igzn<0 .or. igzn>ngzns) then
    !                     error stop "Gas type index out of bounds"
    !                 else
    !                     ind_dom=ind_dom+1
    !                     this%dom_tar_wat_indices(ind_dom)=tar_wat_ind
    !                 end if
    !                 call this%loop_read_tar_wat_init(.false.,this%wat_types,init_sol_types,init_gas_types,nsrz,ngrz,tar_wat_ind,&
    !                     wtype,igzone,igzn,aux_igzone,aux_igzn,gas_chem)
    !             !    !aux_igzone=igzone
    !             !    this%target_waters(tar_wat_ind)=water_types(wtype)
    !             !    ! if (counter_swap==0) then
    !             !    !     call this%target_waters(tar_wat_ind)%set_aq_phase(this%chem_syst%aq_phase)
    !             !    !     !call this%target_waters(tar_wat_ind)%set_indices_aq_species_aq_chem()
    !             !    ! end if
    !             !    if (igzone>0) then
    !             !        this%target_gases(tar_wat_ind)=init_sol_types(igzone)
    !             !        !> chapuza intercambio
    !             !        if (this%reactive_zones(ngzns+nstype*igzn+1)%cat_exch_zone%num_surf_compl>0) then
    !             !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+1))
    !             !        else
    !             !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+igzone))
    !             !        end if
    !             !        if (igzn>0) then
    !             !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+igzone))
    !             !            this%target_gases(tar_wat_ind)=init_gas_types(igzn)
    !             !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*igzn+igzone))
    !             !            call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
    !             !        else
    !             !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+igzone))
    !             !        end if
    !             !        call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
    !             !    else if (igzn>0) then
    !             !        this%target_gases(tar_wat_ind)=init_gas_types(igzn)
    !             !        call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(igzn))
    !             !        call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
    !             !        call gas_chem%set_reactive_zone(this%reactive_zones(igzn))
    !             !        call this%target_waters(tar_wat_ind)%set_gas_chemistry(gas_chem)
    !             !    else
    !             !        call gas_chem%set_reactive_zone(react_zone)
    !             !        call this%target_waters(tar_wat_ind)%set_gas_chemistry(gas_chem)
    !             !    end if
    !             !    !> we check if there is a new reactive zone
    !             !    if (aux_igzone==0 .or. aux_igzone/=igzone .or. aux_igzn/=igzn) then !> we assume target waters are grouped by their reactive zones
    !             !        call this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%set_speciation_alg_dimensions(flag_comp)
    !             !        call this%target_waters(tar_wat_ind)%set_ind_species()
    !             !        call this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%compute_speciation_alg_arrays(&
    !             !            flag_Se,swap)
    !             !        if (flag_Se.eqv..true.) then !> we swap indices of species
    !             !            aux_swap(1)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))
    !             !            aux_swap(2)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))
    !             !            !     this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%speciation_alg%num_prim_species)
    !             !            !print *, this%target_waters(tar_wat_ind)%ind_var_act_species
    !             !            !print *, this%target_waters(tar_wat_ind)%ind_sec_species
    !             !            this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))=aux_swap(2)
    !             !            this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))=aux_swap(1)
    !             !            ! this%target_waters(tar_wat_ind)%ind_prim_species(swap(1))=aux_swap(2) !> index of primary species
    !             !            ! this%target_waters(tar_wat_ind)%ind_sec_species(swap(2)-&
    !             !            !     this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%speciation_alg%num_prim_species)=&
    !             !            !     aux_swap(1) !> index of secondary species
    !             !        end if
    !             !    else if (aux_igzone>0 .or. aux_igzn>0) then !> indices remain the same because reactive zone is the same
    !             !        this%target_waters(tar_wat_ind)%ind_var_act_species=this%target_waters(tar_wat_ind-1)%ind_var_act_species
    !             !        !this%target_waters(tar_wat_ind)%ind_sec_species=this%target_waters(tar_wat_ind-1)%ind_sec_species
    !             !    end if
    !             !    ! if (flag_Se.eqv..true.) then
    !             !    !     counter_swap=counter_swap+1
    !             !    !     aux_swap=swap
    !             !    !     !print *, this%target_waters(tar_wat_ind)%indices_aq_species
    !             !    !     do j=tar_wat_ind+1,this%num_target_waters
    !             !    !         call this%target_waters(j)%set_aq_phase(this%chem_syst%aq_phase)
    !             !    !         call this%target_waters(j)%gas_chemistry%reactive_zone%set_speciation_alg_dimensions(&
    !             !    !             flag_comp)
    !             !    !         call this%target_waters(j)%set_ind_prim_sec_species()
    !             !    !         this%target_waters(j)%ind_prim_species(swap(1))=swap(2)
    !             !    !         this%target_waters(j)%ind_sec_species(swap(2)-&
    !             !    !             this%target_waters(j)%gas_chemistry%reactive_zone%speciation_alg%num_prim_species)=swap(1) !> chapuza
    !             !    !         !call this%target_waters(j)%set_indices_aq_species_aq_chem()
    !             !    !         ! this%target_waters(j)%indices_aq_species(swap(2))=swap(1)
    !             !    !         ! this%target_waters(j)%indices_aq_species(swap(1))=swap(2)
    !             !    !     end do
    !             !    ! end if
    !             !    !print *, this%target_waters(tar_wat_ind)%ind_var_act_species
    !             !    !print *, this%target_waters(tar_wat_ind)%ind_sec_species
    !             !    ! if (aux_swap(1)>0 .AND. aux_swap(2)>0) then
    !             !    !     !call this%target_waters(tar_wat_ind)%set_indices_aq_species_aq_chem()
    !             !    !     this%target_waters(tar_wat_ind)%ind_prim_species(aux_swap(1))=aux_swap(2)
    !             !    !     this%target_waters(tar_wat_ind)%ind_sec_species(swap(2)-&
    !             !    !         this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%speciation_alg%num_prim_species)=swap(1) !> chapuza
    !             !    ! end if
    !             !    if (this%target_waters(tar_wat_ind)%gas_chemistry%mineral_zone%num_minerals_kin<&
    !             !        this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%chem_syst%num_minerals_kin) then
    !             !        call this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%compute_U_SkT_prod(&
    !             !            this%target_waters(tar_wat_ind)%gas_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
    !             !            this%target_waters(tar_wat_ind)%gas_chemistry%mineral_zone%ind_min_chem_syst(1:&
    !             !            this%target_waters(tar_wat_ind)%gas_chemistry%mineral_zone%num_minerals_kin))
    !             !    else
    !             !        call this%target_gases%reactive_zone%compute_U_SkT_prod()
    !             !    end if
    !             !    call this%target_waters(tar_wat_ind)%allocate_reaction_rates()
    !             !    call this%target_waters(tar_wat_ind)%set_indices_rk()
    !             !    aux_igzone=igzone
    !             !    aux_igzn=igzn
    !             end do
    !         end if
        else
            continue
        end if
    end do
    ! close(unit) !> we close the file with target waters associated with their gases and/or gases
    ! this%target_waters_init=this%target_waters
    ! this%target_gases_init=this%target_gases
    ! if (allocated(this%target_gases)) then
    !     this%target_gases_init=this%target_gases
    ! end if
    ! call this%write_aq_comps_init(root) !> write aqueous components of initial target waters
    this%target_gases_init=this%target_gases !> we assign the initial target gases to the target gases
end subroutine