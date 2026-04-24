!> Reads initial target waters and their associated target solids and/or gases
!> We assume file has already been opened
subroutine read_target_waters_init(this,root,nsrz,ngrz)
    use chemistry_m, only: chemistry_c, aqueous_chemistry_c, solid_chemistry_c, gas_chemistry_c, reactive_zone_c, &
        mineral_zone_c
    implicit none
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: root !> root name of the file
    integer(kind=4), intent(in) :: nsrz !> number of solid reactive zones
    integer(kind=4), intent(in) :: ngrz !> number of gas reactive zones
    
    integer(kind=4) :: ind_rech,ind_dom,counter_swap,ind,i,j,k,m,nwtype,num_tar_wat,tar_wat_ind,wtype,istype,nstype,nbwtype,bwtype,&
        mix_wat_ind,ngzns,num_tar_wat_rech,aux_col,flag_wat_type,num_tar_wat_bd,ind_bd,unit,ind_bar_wat,&
        num_sol_types,num_gas_types,int_wat_size,iszn,igzn,ind_bar_gas,ind_bar_sol,first_wat,last_wat,first_sol,&
        last_sol,first_gas,last_gas,int_sol_size,int_gas_size,tar_sol_ind,tar_gas_ind
    integer(kind=4), allocatable :: ind_tar_solids(:),ind_tar_gases(:),ind_sol_zones(:),ind_gas_zones(:)
    real(kind=8), allocatable :: c_nc(:),u_init(:,:),c1_init(:),c2_init(:),c2_ig(:),gamma_2aq(:)
    integer(kind=4) :: aux_istype, aux_igzn
    character(len=256) :: label,str,int_wat,int_sol,int_gas
    character(len=:), allocatable :: int_wat_trim,first_str,last_str,int_sol_trim,int_gas_trim
    logical :: flag_comp,flag_surf,flag_aq_phase,flag, flag_Se
    integer(kind=4), allocatable :: swap(:),aux_swap(:)
    type(solid_chemistry_c), target :: solid_chem !> default solid chemistry object
    type(reactive_zone_c), target :: react_zone   !> default reactive zone object
    type(reactive_zone_c), allocatable :: aux_react_zones(:) !> auxiliary reactive zones (chapuza)
    type(mineral_zone_c), target :: min_zone   !> default mineral zone object
    
    nwtype=this%num_wat_types !> number of water types
    if (nsrz>0) then
        nstype=this%num_mineral_zones !> AQUI HAY QUE PONER MATERIALES
    else
        nstype=0 !> chapuza
        !allocate(this%target_solids(1)) !> we allocate a default target solid zone (chapuza)
        !allocate(this%mineral_zones(1)) !> we allocate a default mineral zone (chapuza)
        !call this%target_solids(1)%set_reactive_zone(this%reactive_zones(1)) !> we set the first reactive zone in default solid chemistry object (chapuza)
        !call this%mineral_zones(1)%set_chem_syst_min_zone(this%chem_syst) !> we set the chemical system in default mineral zone (chapuza)
        !call this%target_solids(1)%set_mineral_zone(this%mineral_zones(1)) !> we set the default mineral zone in default solid chemistry object (chapuza)
    end if
    if (ngrz>0) then
        ngzns=this%num_gas_zones !> number of gas zones
    else
        ngzns=0 !> chapuza
    end if
    
    !flag_comp=.true. !> by default

    
    if (nstype==0) then !> chapuza
        allocate(aux_react_zones(this%num_reactive_zones))
        do i=1,this%num_reactive_zones
            call aux_react_zones(i)%assign_react_zone(this%reactive_zones(i))
        end do
        if (allocated(this%reactive_zones)) deallocate(this%reactive_zones)
        call this%allocate_reactive_zones(this%num_reactive_zones+1)
        do i=1,size(aux_react_zones)
            call this%reactive_zones(i)%assign_react_zone(aux_react_zones(i))
        end do
        call this%reactive_zones(this%num_reactive_zones)%set_chem_syst_react_zone(this%chem_syst)
        call this%reactive_zones(this%num_reactive_zones)%set_CV_params(this%CV_params)
        call this%reactive_zones(this%num_reactive_zones)%set_speciation_alg_dimensions(.true.) !> we set the speciation algebra dimensions in the default reactive zone
        call this%reactive_zones(this%num_reactive_zones)%set_ind_eq_reacts() !> chapuza
        call this%reactive_zones(this%num_reactive_zones)%set_stoich_mat_react_zone() !> chapuza
        call this%reactive_zones(this%num_reactive_zones)%set_ind_gases_stoich_mat() !> chapuza
        call this%reactive_zones(this%num_reactive_zones)%set_ind_mins_stoich_mat() !> chapuza
        allocate(swap(2)) !> we allocate the swap array (chapuza)
        call this%reactive_zones(this%num_reactive_zones)%compute_speciation_alg_arrays(flag_Se,swap) !> we compute the speciation algebra arrays in the default reactive zone
        deallocate(aux_react_zones)
        call this%allocate_mineral_zones(1) !> chapuza
        call this%mineral_zones(1)%set_chem_syst_min_zone(this%chem_syst) !> we set the chemical system in default mineral zone
        call solid_chem%set_reactive_zone(this%reactive_zones(this%num_reactive_zones)) !> we set default reactive zone in default solid chemistry object
        call solid_chem%set_mineral_zone(this%mineral_zones(1)) !> we set default mineral zone in default solid chemistry object
        ! call this%allocate_target_solids() !> we allocate default target solids (chapuza)
        ! do i=1,this%num_target_solids
        !     this%target_solids(i)=solid_chem !> we set the default target solids (chapuza)
        ! end do
    end if
    
    ind_rech=0 !> counter recharge waters
    ind_bd=0 !> counter boundary waters
    ind_dom=0 !> counter domain waters

    unit=57 !> arbitrary unit number
    open(unit, file=root//'_tar_wat.dat', status='old', action='read') !> we open the file with target waters associated with their solids and/or gases
    do
        read(unit,*) label
        if (label=='end') then
            exit
        else if (label=='TARGET WATERS') then
            read(unit,*) num_tar_wat !> total number of target waters
            call this%allocate_target_waters(num_tar_wat)
            call this%allocate_target_waters_init(num_tar_wat)
            read(unit,*) num_tar_wat_rech !> number of recharge waters
            call this%allocate_rech_waters_indices(num_tar_wat_rech)            
            read(unit,*) num_tar_wat_bd !> number of boundary waters
            call this%allocate_bd_waters_indices(num_tar_wat_bd)
            call this%allocate_dom_tar_wat_indices(this%num_target_waters-this%num_rech_waters-this%num_bd_waters)
            if (nstype==0) then
                call this%allocate_target_solids(this%num_target_waters) !> we assume bijection with domain target waters (chapuza)
                do i=1,this%num_target_solids
                    this%target_solids(i)=solid_chem !> we set the default solid chemistry object
                end do  
            end if
            if (ngzns==0) then
                call this%allocate_target_gases(this%num_target_waters) !> we assume bijection with domain target waters (chapuza)
            end if
            aux_istype=0 !> auxiliary solid type index
            aux_igzn=0 !> auxiliary gas zone index
            !if (num_tar_wat_rech>0 .or. num_tar_wat_bd>0) then !> we have recharge or boundary waters
                do 
                    read(unit,*) int_wat, wtype, int_sol, int_gas, flag_wat_type
                    !print *, int_wat, wtype, int_sol, int_gas, flag_wat_type
                    flag=.true. !> we set the flag to true (we will use it if we read a domain water)
                    int_wat_trim=trim(int_wat)
                    int_sol_trim=trim(int_sol)
                    int_gas_trim=trim(int_gas)
                    ind_bar_wat=index(int_wat_trim,'-')
                    ind_bar_sol=index(int_sol_trim,'-')
                    ind_bar_gas=index(int_gas_trim,'-')
                    if (wtype<1 .or. wtype>nwtype) then
                        error stop "Water type index out of bounds"
                    ! else if (int_sol<0 .or. int_sol>this%num_materials) then
                    !     error stop "Solid zone index out of bounds"
                    ! else if (int_gas<0 .or. int_gas>this%num_gas_zones) then
                    !     error stop "Gas zone index out of bounds"
                    !else if (flag_wat_type==0) then !> boundary water
                    !    !ind_bd=ind_bd+1 !> counter boundary waters
                    !    !this%bd_waters_indices(ind_bd)=tar_wat_ind
                    !else if (flag_wat_type==2) then !> external water
                    !    !ind_ext=ind_ext+1 !> counter external waters
                    !    !this%ext_waters_indices(ind_ext)=tar_wat_ind
                    !else if (flag_wat_type==1) then !> domain water
                    !    !ind_dom=ind_dom+1 !> counter domain waters
                        !this%dom_tar_wat_indices(ind_dom)=tar_wat_ind
                        !flag=.false.
                    else if (flag_wat_type<0 .or. flag_wat_type>2) then
                        error stop "Water type flag not implemented"
                    end if
                    if (ind_bar_wat>0) then !> we have a range of target waters
                        first_str=int_wat_trim(1:ind_bar_wat-1) !> first target of int_wat
                        last_str=int_wat_trim(ind_bar_wat+1:) !> last target of int_wat
                        read(first_str,*) first_wat
                        read(last_str,*) last_wat
                        int_wat_size=last_wat-first_wat+1 !> size of the waters interval
                        if (first_wat<1 .or. first_wat>last_wat .or. last_wat>num_tar_wat) then
                            error stop "Target water index out of bounds"
                        !else if (tar_wat_ind<1 .or. tar_wat_ind>this%num_target_waters) then
                        !    error stop "Target water index out of bounds"
                        !else if (wtype<1 .or. wtype>nwtype) then
                        !    error stop "Water type index out of bounds"
                        !else if (ind_tar_sol<0 .or. ind_tar_sol>this%num_target_solids) then
                        !    error stop "Target solid index out of bounds"
                        !else if (int_gas<0 .or. int_gas>ngzns) then
                        !    error stop "Gas type index out of bounds"
                        else if (flag_wat_type==0) then !> boundary water
                            do i=1,int_wat_size
                                this%bd_waters_indices(ind_bd+i)=first_wat+i-1 !> we set the indices of boundary waters
                            end do
                            ind_bd=ind_bd+int_wat_size !> counter boundary waters
                        else if (flag_wat_type==2) then !> recharge water
                            !this%ext_waters_indices(ind_ext)=tar_wat_ind
                            do i=1,int_wat_size
                                this%rech_waters_indices(ind_rech+i)=first_wat+i-1 !> we set the indices of recharge waters
                            end do
                            ind_rech=ind_rech+int_wat_size !> counter recharge waters
                        else if (flag_wat_type==1) then !> domain water
                            !this%dom_tar_wat_indices(ind_dom)=tar_wat_ind
                            do i=1,int_wat_size
                                this%tar_wat_indices(ind_dom+i)=first_wat+i-1 !> we set the indices of domain waters
                            end do
                            ind_dom=ind_dom+int_wat_size !> counter domain waters
                            flag=.false.
                        !else
                        !    error stop "Water type flag out of bounds"
                        end if
                        !> target solids interval
                        if (ind_bar_sol>0) then !> we have a range of target solids
                            first_str=int_sol_trim(1:ind_bar_sol-1) !> first target solid of int_sol
                            last_str=int_sol_trim(ind_bar_sol+1:) !> last target solid of int_sol
                            read(first_str,*) first_sol
                            read(last_str,*) last_sol
                            int_sol_size=last_sol-first_sol+1 !> size of the solids interval
                            if (first_sol<0 .or. first_sol>last_sol .or. last_sol>this%num_target_solids) then
                                error stop "Target solid index out of bounds"
                            end if
                            if (int_sol_size/=int_wat_size) then
                                error stop "Dimension error: number of target solids not consistent with target waters interval"
                            end if
                            allocate(ind_tar_solids(int_sol_size), ind_sol_zones(int_sol_size))
                            ind_tar_solids = [(i, i=first_sol, last_sol)]
                            do i=1,int_sol_size
                                ind_sol_zones(i)=this%target_solids(ind_tar_solids(i))%id !> we get the solid zone indices associated to target solids
                            end do
                            ! print *, "ind_tar_solids", ind_tar_solids
                            ! print *, "ind_sol_zones", ind_sol_zones
                        else !> we have a single target solid
                            read(int_sol_trim,*) tar_sol_ind !> target solid index
                            if (tar_sol_ind<0 .or. tar_sol_ind>this%num_target_solids) then
                                error stop "Target solid index out of bounds"
                            else
                                allocate(ind_tar_solids(int_wat_size),ind_sol_zones(int_wat_size))
                                ind_tar_solids=tar_sol_ind !> we set all the target solid indices to the same value
                                if (tar_sol_ind>0) then
                                    ind_sol_zones=this%target_solids(tar_sol_ind)%id !> we get the solid zone index associated to target solid
                                else 
                                    ind_sol_zones=0 !> we set to zero if no solids involved
                                end if
                            end if
                        end if
                        !> target gases interval
                        if (ind_bar_gas>0) then !> we have a range of target gases
                            first_str=int_gas_trim(1:ind_bar_gas-1) !> first target gas of int_gas
                            last_str=int_gas_trim(ind_bar_gas+1:) !> last target gas of int_gas
                            read(first_str,*) first_gas
                            read(last_str,*) last_gas
                            int_gas_size=last_gas-first_gas+1 !> size of the gases interval
                            if (first_gas<0 .or. first_gas>last_gas .or. last_gas>this%num_target_gases) then
                                error stop "Target gas index out of bounds"
                            end if
                            if (int_gas_size/=int_wat_size) then
                                error stop "Dimension error: number of target gases not consistent with target waters interval"
                            end if
                            allocate(ind_tar_gases(int_gas_size), ind_gas_zones(int_gas_size))
                            ind_tar_gases = [(i, i=first_gas, last_gas)]
                            do i=1,int_gas_size
                                ind_gas_zones(i)=this%target_gases(ind_tar_gases(i))%id !> we get the gas zone indices associated to target gases
                            end do
                            ! print *, "ind_tar_gases", ind_tar_gases
                            ! print *, "ind_gas_zones", ind_gas_zones
                        else !> we have a single target gas
                            read(int_gas_trim,*) tar_gas_ind !> target gas index
                            if (tar_gas_ind<0 .or. tar_gas_ind>this%num_target_gases) then
                                error stop "Target gas index out of bounds"
                            else
                                allocate(ind_tar_gases(int_wat_size),ind_gas_zones(int_wat_size))
                                ind_tar_gases=tar_gas_ind !> we set all the target gas indices to the same value
                                if (tar_gas_ind>0) then
                                    ind_gas_zones=this%target_gases(tar_gas_ind)%id !> we get the gas zone index associated to target gas
                                else
                                    ind_gas_zones=0 !> we set to zero if no gases involved
                                end if
                            end if
                        end if
                        ! if (int_sol>0) then
                        !     ind_tar_sol=this%get_tar_sol_ind(int_sol) !> we get the indices of the target solids associated to solid zone int_sol
                        !     if (size(ind_tar_sol)/=int_size) then
                        !         error stop "Dimension error: number of target solids associated to solid zone not consistent with &
                        !             int_wat size"
                        !     end if
                        ! else
                        !     allocate(ind_tar_sol(int_size))
                        !     ind_tar_sol=0 !> we set to zero the target solid indices
                        ! end if
                        ! if (int_gas>0) then
                        !     ind_tar_gases=this%get_tar_gas_ind(int_gas) !> we get the indices of the target gases associated to gas zone int_gas
                        !     if (size(ind_tar_gases)/=int_size) then
                        !         error stop "Dimension error: number of target gases associated to gas zone not consistent with & 
                        !             int_wat size"
                        !     end if
                        ! else
                        !     allocate(ind_tar_gases(int_size))
                        !     ind_tar_gases=0 !> we set to zero the target gas indices
                        ! end if
                        do i=1,int_wat_size
                            call this%read_tar_wat_line(flag,ind_sol_zones(i),ind_gas_zones(i),first_wat+i-1,wtype,&
                                ind_tar_solids(i),ind_tar_gases(i),aux_istype,aux_igzn)                        
                        end do
                        if (last_wat==num_tar_wat) then
                            exit
                        end if
                        deallocate(ind_tar_solids,ind_tar_gases,ind_sol_zones,ind_gas_zones)
                    else
                        !> we have a single target water
                        read(int_wat_trim,*) tar_wat_ind !> target water index
                        read(int_sol_trim,*) tar_sol_ind !> target solid index
                        read(int_gas_trim,*) tar_gas_ind !> target gas index
                        if (tar_wat_ind<1 .or. tar_wat_ind>num_tar_wat) then
                            error stop "Target water index out of bounds"
                        else if (ind_bar_sol>0) then
                            error stop "Dimension error: a single target water cannot point to a range of target solids"
                        else if (ind_bar_gas>0) then
                            error stop "Dimension error: a single target water cannot point to a range of target gases"
                        else if (tar_sol_ind<0 .or. tar_sol_ind>this%num_target_solids) then
                            error stop "Target solid index out of bounds"
                        else if (tar_gas_ind<0 .or. tar_gas_ind>this%num_target_gases) then
                            error stop "Target gas index out of bounds"
                        else if (flag_wat_type==0) then !> boundary water
                            ind_bd=ind_bd+1 !> counter boundary waters
                            this%bd_waters_indices(ind_bd)=tar_wat_ind
                        else if (flag_wat_type==2) then !> recharge water
                            ind_rech=ind_rech+1 !> counter recharge waters
                            this%rech_waters_indices(ind_rech)=tar_wat_ind
                        else if (flag_wat_type==1) then !> domain water
                            ind_dom=ind_dom+1 !> counter domain waters
                            this%tar_wat_indices(ind_dom)=tar_wat_ind
                            !call this%loop_read_tar_wat_init(flag,this%wat_types,init_sol_types,init_gas_types,nsrz,ngrz,tar_wat_ind,wtype,&
                            !    istype,int_gas,aux_istype,aux_int_gas,solid_chem)
                        end if
                        if (tar_sol_ind>0) then !> we have a solid zone
                            iszn=this%target_solids(tar_sol_ind)%id !> we get the indices of the target solids associated to solid zone int_sol
                            ! if (size(ind_sol_zn)/=1) then
                            !     error stop "Dimension error: number of target solids associated to solid zone must be 1"
                            ! end if
                        else
                            !allocate(ind_tar_solids(1))
                            iszn=0 !> we set to zero the solid zone index
                        end if
                        if (tar_gas_ind>0) then !> we have a gas zone
                            igzn=this%target_gases(tar_gas_ind)%id !> we get the index of the gas zone associated to target gas
                            ! if (size(ind_tar_gases)/=1) then
                            !     error stop "Dimension error: number of target gases associated to gas zone must be 1"
                            ! end if
                        else
                            !allocate(ind_tar_gases(1))
                            igzn=0 !> we set to zero the gas zone index
                        end if
                        call this%read_tar_wat_line(flag,iszn,igzn,tar_wat_ind,wtype,&
                                tar_sol_ind,tar_gas_ind,aux_istype,aux_igzn)
                            !this%target_waters_init(ind_tar_solids)=init_sol_types(sol_zone) !> we set the initial solid type
                            !call this%target_solids_init(ind_dom)%set_target(ind_dom) !> we set the target index
                        if (tar_wat_ind==num_tar_wat) then
                            exit
                        end if
                        !deallocate(ind_tar_solids,ind_tar_gases)
                    end if
                    !allocate(rk(this%target_waters(tar_wat_ind)%indices_rk%num_cols)) !> chapuza
                    !call this%target_waters(tar_wat_ind)%compute_rk(rk)
                !    !aux_istype=istype
                !    this%target_waters(tar_wat_ind)=water_types(wtype)
                !    ! if (counter_swap==0) then
                !    !     call this%target_waters(tar_wat_ind)%set_aq_phase(this%chem_syst%aq_phase)
                !    !     !call this%target_waters(tar_wat_ind)%set_indices_aq_species_aq_chem()
                !    ! end if
                !    if (istype>0) then
                !        this%target_solids(tar_wat_ind)=init_sol_types(istype)
                !        !> chapuza intercambio
                !        if (this%reactive_zones(ngzns+nstype*int_gas+1)%cat_exch_zone%num_surf_compl>0) then
                !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+1))
                !        else
                !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+istype))
                !        end if
                !        if (int_gas>0) then
                !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+istype))
                !            this%target_gases(tar_wat_ind)=init_gas_types(int_gas)
                !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+istype))
                !            call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
                !        else
                !            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+istype))
                !        end if
                !        call this%target_waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(tar_wat_ind))
                !    else
                !        if (int_gas>0) then
                !            this%target_gases(tar_wat_ind)=init_gas_types(int_gas)
                !            ! call this%reactive_zones(int_gas)%set_ind_eq_reacts() !> chapuza
                !            ! call this%reactive_zones(int_gas)%set_stoich_mat_react_zone() !> chapuza
                !            ! call this%reactive_zones(int_gas)%set_ind_gases_stoich_mat() !> chapuza
                !            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(int_gas))
                !            call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
                !            call solid_chem%set_reactive_zone(this%reactive_zones(int_gas))
                !        else
                !            ! call react_zone%set_ind_eq_reacts() !> chapuza
                !            ! call react_zone%set_stoich_mat_react_zone() !> chapuza
                !            ! call react_zone%set_ind_gases_stoich_mat() !> chapuza
                !            ! call react_zone%set_ind_mins_stoich_mat() !> chapuza
                !            call solid_chem%set_reactive_zone(react_zone)
                !        end if
                !        call solid_chem%set_mineral_zone(min_zone)
                !        call this%target_waters(tar_wat_ind)%set_solid_chemistry(solid_chem)
                !    end if
                !    !> we check if there is a new reactive zone
                !    if (aux_istype==0 .or. aux_istype/=istype .or. aux_int_gas/=int_gas) then !> we assume target waters are grouped by their reactive zones
                !        call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_speciation_alg_dimensions(flag_comp)
                !        call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_ind_eq_reacts() !> chapuza
                !        call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_stoich_mat_react_zone() !> chapuza
                !        call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_ind_gases_stoich_mat() !> chapuza
                !        call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_ind_mins_stoich_mat() !> chapuza
                !        call this%target_waters(tar_wat_ind)%set_ind_species()
                !        call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
                !            flag_Se,swap)
                !        if (flag_Se.eqv..true.) then !> we swap indices of species
                !            aux_swap(1)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))
                !            aux_swap(2)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))
                !            !     this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
                !            !print *, this%target_waters(tar_wat_ind)%ind_var_act_species
                !            !print *, this%target_waters(tar_wat_ind)%ind_sec_species
                !            this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))=aux_swap(2)
                !            this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))=aux_swap(1)
                !            ! this%target_waters(tar_wat_ind)%ind_prim_species(swap(1))=aux_swap(2) !> index of primary species
                !            ! this%target_waters(tar_wat_ind)%ind_sec_species(swap(2)-&
                !            !     this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=&
                !            !     aux_swap(1) !> index of secondary species
                !        end if
                !    else if (aux_istype>0 .or. aux_int_gas>0) then !> indices remain the same because reactive zone is the same
                !        this%target_waters(tar_wat_ind)%ind_var_act_species=this%target_waters(tar_wat_ind-1)%ind_var_act_species
                !        !this%target_waters(tar_wat_ind)%ind_sec_species=this%target_waters(tar_wat_ind-1)%ind_sec_species
                !    end if
                !    print *, this%target_waters(tar_wat_ind)%ind_var_act_species
                !!> Chapuza
                !    !if (associated(this%target_waters(tar_wat_ind)%solid_chemistry%mineral_zone)) then
                !        if (this%target_waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin<&
                !            this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin) then
                !            call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod(&
                !                this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                !                this%target_waters(tar_wat_ind)%solid_chemistry%mineral_zone%ind_min_chem_syst(1:&
                !                this%target_waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin))
                !        else
                !            call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod()
                !        end if
                !    !end if
                !    call this%target_waters(tar_wat_ind)%allocate_reaction_rates()
                !    call this%target_waters(tar_wat_ind)%set_indices_rk()
                !    aux_istype=istype
                !    aux_int_gas=int_gas
                end do
            !else !> we have no recharge or boundary waters
            !    do! i=1,this%num_target_waters
            !        read(unit,*) int_wat, wtype, int_sol, int_gas
            !        int_wat_trim=trim(int_wat)
            !        ind_bar=index(int_wat_trim,'-')
            !        if (wtype<1 .or. wtype>nwtype) then
            !            error stop "Water type index out of bounds"
            !        else if (ind_tar_solids<0 .or. ind_tar_solids>this%num_target_solids) then
            !            error stop "Target solid index out of bounds"
            !        else if (ind_tar_gases<0 .or. ind_tar_gases>this%num_target_gases) then
            !            error stop "Target gas index out of bounds"
            !        end if
            !        if (ind_bar>0) then !> we have a range of target waters
            !            first_str=int_wat_trim(1:ind_bar-1) !> first target of int_wat
            !            last_str=int_wat_trim(ind_bar+1:) !> last target of int_wat
            !            read(first_str,*) first
            !            read(last_str,*) last
            !            if (first<1 .or. first>last .or. last>num_tar_wat) then
            !                error stop "Target water index out of bounds"
            !            !else if (tar_wat_ind<1 .or. tar_wat_ind>this%num_target_waters) then
            !            !    error stop "Target water index out of bounds"
            !            !else if (wtype<1 .or. wtype>nwtype) then
            !            !    error stop "Water type index out of bounds"
            !            !else if (ind_tar_solids<0 .or. ind_tar_solids>this%num_target_solids) then
            !            !    error stop "Target solid index out of bounds"
            !            !else if (int_gas<0 .or. int_gas>ngzns) then
            !            !    error stop "Gas type index out of bounds"
            !            !else if (flag_wat_type==0) then !> boundary water
            !            !    ind_bd=ind_bd+1 !> counter boundary waters
            !            !    this%bd_waters_indices(ind_bd)=tar_wat_ind
            !            !else if (flag_wat_type==2) then !> external water
            !            !    ind_ext=ind_ext+1 !> counter external waters
            !            !    this%ext_waters_indices(ind_ext)=tar_wat_ind
            !            !else if (flag_wat_type==1) then !> domain water
            !            !    ind_dom=ind_dom+1 !> counter domain waters
            !            !    this%dom_tar_wat_indices(ind_dom)=tar_wat_ind
            !            !else
            !            !    error stop "Water type flag out of bounds"
            !            end if
            !            do i=first,last
            !                this%dom_tar_wat_indices(i)=i
            !                call this%read_tar_wat_line(.false.,nsrz,ngrz,i,wtype,&
            !                    int_sol,int_gas,aux_istype,aux_int_gas)
            !            end do
            !            if (last==num_tar_wat) then
            !                exit
            !            end if
            !        else
            !            !> we have a single target water
            !            read(int_wat_trim,*) tar_wat_ind !> target water index
            !            if (tar_wat_ind<1 .or. tar_wat_ind>num_tar_wat) then
            !                error stop "Target water index out of bounds"
            !            else
            !                this%dom_tar_wat_indices(tar_wat_ind)=tar_wat_ind
            !                call this%read_tar_wat_line(.false.,nsrz,ngrz,tar_wat_ind,wtype,&
            !                    int_sol,int_gas,aux_istype,aux_int_gas)
            !                !call this%loop_read_tar_wat_init(flag,this%wat_types,init_sol_types,init_gas_types,nsrz,ngrz,tar_wat_ind,wtype,&
            !                !    istype,int_gas,aux_istype,aux_int_gas,solid_chem)
            !                !this%target_waters_init(ind_tar_solids)=init_sol_types(sol_zone) !> we set the initial solid type
            !                !call this%target_solids_init(ind_dom)%set_target(ind_dom) !> we set the target index
            !            end if
            !            if (tar_wat_ind==num_tar_wat) then
            !                exit
            !            end if
            !        end if
            !    end do
            !    !do i=1,this%num_target_waters
            !    !    read(unit,*) tar_wat_ind, wtype, istype, int_gas
            !    !    if (tar_wat_ind<1 .or. tar_wat_ind>this%num_target_waters) then
            !    !        error stop "Target water index out of bounds"
            !    !    else if (wtype<1 .or. wtype>nwtype) then
            !    !        error stop "Water type index out of bounds"
            !    !    else if (istype<0 .or. istype>nstype) then
            !    !        error stop "Solid type index out of bounds"
            !    !    else if (int_gas<0 .or. int_gas>ngzns) then
            !    !        error stop "Gas type index out of bounds"
            !    !    else
            !    !        ind_dom=ind_dom+1
            !    !        this%dom_tar_wat_indices(ind_dom)=tar_wat_ind
            !    !    end if
            !    !    call this%loop_read_tar_wat_init(.false.,this%wat_types,init_sol_types,init_gas_types,nsrz,ngrz,tar_wat_ind,&
            !    !        wtype,istype,int_gas,aux_istype,aux_int_gas,solid_chem)
            !    !!    !aux_istype=istype
            !    !!    this%target_waters(tar_wat_ind)=water_types(wtype)
            !    !!    ! if (counter_swap==0) then
            !    !!    !     call this%target_waters(tar_wat_ind)%set_aq_phase(this%chem_syst%aq_phase)
            !    !!    !     !call this%target_waters(tar_wat_ind)%set_indices_aq_species_aq_chem()
            !    !!    ! end if
            !    !!    if (istype>0) then
            !    !!        this%target_solids(tar_wat_ind)=init_sol_types(istype)
            !    !!        !> chapuza intercambio
            !    !!        if (this%reactive_zones(ngzns+nstype*int_gas+1)%cat_exch_zone%num_surf_compl>0) then
            !    !!            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+1))
            !    !!        else
            !    !!            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+istype))
            !    !!        end if
            !    !!        if (int_gas>0) then
            !    !!            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+istype))
            !    !!            this%target_gases(tar_wat_ind)=init_gas_types(int_gas)
            !    !!            call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+nstype*int_gas+istype))
            !    !!            call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
            !    !!        else
            !    !!            call this%target_solids(tar_wat_ind)%set_reactive_zone(this%reactive_zones(ngzns+istype))
            !    !!        end if
            !    !!        call this%target_waters(tar_wat_ind)%set_solid_chemistry(this%target_solids(tar_wat_ind))
            !    !!    else if (int_gas>0) then
            !    !!        this%target_gases(tar_wat_ind)=init_gas_types(int_gas)
            !    !!        call this%target_gases(tar_wat_ind)%set_reactive_zone(this%reactive_zones(int_gas))
            !    !!        call this%target_waters(tar_wat_ind)%set_gas_chemistry(this%target_gases(tar_wat_ind))
            !    !!        call solid_chem%set_reactive_zone(this%reactive_zones(int_gas))
            !    !!        call this%target_waters(tar_wat_ind)%set_solid_chemistry(solid_chem)
            !    !!    else
            !    !!        call solid_chem%set_reactive_zone(react_zone)
            !    !!        call this%target_waters(tar_wat_ind)%set_solid_chemistry(solid_chem)
            !    !!    end if
            !    !!    !> we check if there is a new reactive zone
            !    !!    if (aux_istype==0 .or. aux_istype/=istype .or. aux_int_gas/=int_gas) then !> we assume target waters are grouped by their reactive zones
            !    !!        call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%set_speciation_alg_dimensions(flag_comp)
            !    !!        call this%target_waters(tar_wat_ind)%set_ind_species()
            !    !!        call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(&
            !    !!            flag_Se,swap)
            !    !!        if (flag_Se.eqv..true.) then !> we swap indices of species
            !    !!            aux_swap(1)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))
            !    !!            aux_swap(2)=this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))
            !    !!            !     this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
            !    !!            !print *, this%target_waters(tar_wat_ind)%ind_var_act_species
            !    !!            !print *, this%target_waters(tar_wat_ind)%ind_sec_species
            !    !!            this%target_waters(tar_wat_ind)%ind_var_act_species(swap(1))=aux_swap(2)
            !    !!            this%target_waters(tar_wat_ind)%ind_var_act_species(swap(2))=aux_swap(1)
            !    !!            ! this%target_waters(tar_wat_ind)%ind_prim_species(swap(1))=aux_swap(2) !> index of primary species
            !    !!            ! this%target_waters(tar_wat_ind)%ind_sec_species(swap(2)-&
            !    !!            !     this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=&
            !    !!            !     aux_swap(1) !> index of secondary species
            !    !!        end if
            !    !!    else if (aux_istype>0 .or. aux_int_gas>0) then !> indices remain the same because reactive zone is the same
            !    !!        this%target_waters(tar_wat_ind)%ind_var_act_species=this%target_waters(tar_wat_ind-1)%ind_var_act_species
            !    !!        !this%target_waters(tar_wat_ind)%ind_sec_species=this%target_waters(tar_wat_ind-1)%ind_sec_species
            !    !!    end if
            !    !!    ! if (flag_Se.eqv..true.) then
            !    !!    !     counter_swap=counter_swap+1
            !    !!    !     aux_swap=swap
            !    !!    !     !print *, this%target_waters(tar_wat_ind)%indices_aq_species
            !    !!    !     do j=tar_wat_ind+1,this%num_target_waters
            !    !!    !         call this%target_waters(j)%set_aq_phase(this%chem_syst%aq_phase)
            !    !!    !         call this%target_waters(j)%solid_chemistry%reactive_zone%set_speciation_alg_dimensions(&
            !    !!    !             flag_comp)
            !    !!    !         call this%target_waters(j)%set_ind_prim_sec_species()
            !    !!    !         this%target_waters(j)%ind_prim_species(swap(1))=swap(2)
            !    !!    !         this%target_waters(j)%ind_sec_species(swap(2)-&
            !    !!    !             this%target_waters(j)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=swap(1) !> chapuza
            !    !!    !         !call this%target_waters(j)%set_indices_aq_species_aq_chem()
            !    !!    !         ! this%target_waters(j)%indices_aq_species(swap(2))=swap(1)
            !    !!    !         ! this%target_waters(j)%indices_aq_species(swap(1))=swap(2)
            !    !!    !     end do
            !    !!    ! end if
            !    !!    !print *, this%target_waters(tar_wat_ind)%ind_var_act_species
            !    !!    !print *, this%target_waters(tar_wat_ind)%ind_sec_species
            !    !!    ! if (aux_swap(1)>0 .AND. aux_swap(2)>0) then
            !    !!    !     !call this%target_waters(tar_wat_ind)%set_indices_aq_species_aq_chem()
            !    !!    !     this%target_waters(tar_wat_ind)%ind_prim_species(aux_swap(1))=aux_swap(2)
            !    !!    !     this%target_waters(tar_wat_ind)%ind_sec_species(swap(2)-&
            !    !!    !         this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=swap(1) !> chapuza
            !    !!    ! end if
            !    !!    if (this%target_waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin<&
            !    !!        this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin) then
            !    !!        call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod(&
            !    !!            this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
            !    !!            this%target_waters(tar_wat_ind)%solid_chemistry%mineral_zone%ind_min_chem_syst(1:&
            !    !!            this%target_waters(tar_wat_ind)%solid_chemistry%mineral_zone%num_minerals_kin))
            !    !!    else
            !    !!        call this%target_waters(tar_wat_ind)%solid_chemistry%reactive_zone%compute_U_SkT_prod()
            !    !!    end if
            !    !!    call this%target_waters(tar_wat_ind)%allocate_reaction_rates()
            !    !!    call this%target_waters(tar_wat_ind)%set_indices_rk()
            !    !!    aux_istype=istype
            !    !!    aux_int_gas=int_gas
            !    !end do
            !end if
        else
            continue
        end if
    end do
    close(unit) !> we close the file with target waters associated with their solids and/or gases
    this%target_waters_init=this%target_waters
    !this%target_solids_init=this%target_solids
    !if (allocated(this%target_gases)) then
    !    this%target_gases_init=this%target_gases
    !end if
    !call this%write_aq_comps_init(root) !> write aqueous components of initial target waters
end subroutine
