subroutine read_init_min_zones_CHEPROO(this,unit,num_rz,surf_chem)
!> We assume all minerals are pure phases
    use chemistry_m, only: chemistry_c
    use aqueous_chemistry_m, only: solid_chemistry_c, reactive_zone_c, mineral_zone_c, mineral_c
    implicit none
    class(chemistry_c) :: this !> chemistry object
    integer(kind=4), intent(in) :: unit !> file
    !type(solid_chemistry_c), intent(out), allocatable :: materials(:) !> initial mineral zones
    integer(kind=4), intent(out) :: num_rz !> number of mineral reactive zones ( we assume each mineral zone has different minerals )
    type(solid_chemistry_c), intent(in), optional :: surf_chem
    
    integer(kind=4) :: num_mins_cst_kin,num_mins_var_kin,num_mins_loc_kin,num_surf_rz,num_gas_rz,num_rz_old,i,j,k,nmtype,nrwtype,&
    icon,num_min_zones,num_mins_rz,num_mins_glob,num_mins_loc,num_mins_loc_eq,ind_rz,imtype,min_ind,num_mins_var_eq,&
    num_mins_cst_eq
    integer(kind=4), allocatable :: init_min_indices(:),min_indices(:,:)
    character(len=256) :: str,constrain,label,min_name,name
    real(kind=8) :: guess,c_tot,temp,vol_frac,react_surf,conc
    logical :: min_flag,flag_comp
    type(mineral_c) :: mineral
    type(reactive_zone_c) :: aux_react_zone
    type(mineral_zone_c), allocatable :: min_zones(:)
    type(reactive_zone_c), allocatable :: react_zones(:),aux_react_zones(:)
    
    flag_comp=.true. !> flag to compute speciation algebra dimensions (by default: TRUE)

    read(unit,*) nmtype
    if (nmtype<0) then
        error stop "Number of mineral zones must be positive"
    else if (nmtype==0) then
        num_rz=0
        return
    else
        !nmrz=1 !> counter mineral reactive zones
        allocate(react_zones(nmtype)) !> chapuza
        call this%allocate_mineral_zones(nmtype)
        allocate(min_indices(this%chem_syst%num_minerals,nmtype))
        num_rz=nmtype !> number of mineral reactive zones
            allocate(this%materials(nmtype))
            num_min_zones=0 !> counter mineral zones
            do
                read(unit,*) imtype, temp !> index of mineral zone, temperature (C)
                if (imtype<1 .or. imtype>nmtype) error stop "Mineral zone index out of range"
                read(unit,*) name !> name of mineral zone
                read(unit,*) str !> headings
                    if (str=='*') then
                        if (num_min_zones<nmtype) then
                            continue
                        else
                            exit
                        end if
                    else if (index(str,'mineral')/=0) then
                        num_min_zones=num_min_zones+1
                        num_mins_loc=0 !> counter minerals in this zone
                        num_mins_loc_eq=0 !> counter minerals in eq in this zone
                        num_mins_loc_kin=0 !> counter minerals kinetic in this zone
                        num_mins_cst_eq=0
                        num_mins_cst_kin=0
                        num_mins_var_eq=0
                        num_mins_var_kin=0
                        do
                            read(unit,*) mineral%name
                            if (mineral%name=='*') then
                                call this%mineral_zones(imtype)%set_chem_syst_min_zone(this%chem_syst)
                                call this%mineral_zones(imtype)%set_id(imtype)
                                call this%materials(imtype)%set_mineral_zone(this%mineral_zones(imtype)) !> chapuza
                                call react_zones(imtype)%set_chem_syst_react_zone(this%chem_syst)
                                call this%materials(imtype)%mineral_zone%set_chem_syst_min_zone(this%chem_syst)
                                call react_zones(imtype)%set_CV_params(this%CV_params)
                                call react_zones(imtype)%allocate_ind_mins(num_mins_loc_eq)
                                call this%materials(imtype)%mineral_zone%allocate_ind_chem_syst_min_zone(num_mins_loc)
                                call this%materials(imtype)%mineral_zone%allocate_ind_min_Sk(num_mins_loc_kin)
                                call this%materials(imtype)%mineral_zone%set_num_mins_eq_min_zone(num_mins_loc_eq)
                                call this%materials(imtype)%mineral_zone%set_num_mins_eq_cst_act_min_zone(num_mins_cst_eq)
                                call this%materials(imtype)%mineral_zone%set_num_mins_kin_cst_act_min_zone(num_mins_cst_kin)
                                call this%materials(imtype)%mineral_zone%set_num_mins_eq_var_act_min_zone(num_mins_var_eq)
                                !call react_zones(imtype)%allocate_non_flowing_species()
                                if (present(surf_chem)) then
                                    call react_zones(imtype)%set_cat_exch_zone(surf_chem%reactive_zone%cat_exch_zone)
                                    call this%materials(imtype)%set_CEC(surf_chem%CEC)
                                else
                                    !call react_zones(imtype)%set_surf_chem(this%surf_chem)
                                    !call this%materials(imtype)%set_surf_chem(this%surf_chem)
                                end if
                                !call react_zones(imtype)%set_num_mins(num_mins_loc)
                                call react_zones(imtype)%set_num_solids()
                                call react_zones(imtype)%set_num_mins_cst_act(num_mins_cst_eq)
                                !call this%materials(imtype)%mineral_zone%set_num_mins_cst_act_min_zone(num_mins_cst_eq+num_mins_cst_kin)
                                call react_zones(imtype)%set_num_mins_var_act(num_mins_var_eq)
                                !call this%materials(imtype)%mineral_zone%set_num_mins_var_act_min_zone(num_mins_var_eq+num_mins_var_kin)
                                call react_zones(imtype)%set_num_non_flowing_species() !> sets number of non-flowing species in reactive zone
                                call react_zones(imtype)%set_speciation_alg_dimensions(flag_comp)
                                !> aqui habria que comprobar si se estan repitiendo zonas reactivas
                                call this%materials(imtype)%set_reactive_zone(react_zones(imtype))
                                call this%materials(imtype)%compute_num_solids_solid_chem()
                                call this%materials(imtype)%allocate_vol_fracts()
                                call this%materials(imtype)%allocate_react_surfaces()
                                call this%materials(imtype)%allocate_conc_solids()
                                call this%materials(imtype)%allocate_log_act_coeffs_solid_chem()
                                call this%materials(imtype)%allocate_activities()
                                call this%materials(imtype)%set_temp(temp+273.18) !> Kelvin
                                call this%materials(imtype)%allocate_var_act_species_indices(num_mins_var_eq+num_mins_var_kin)
                                call this%materials(imtype)%allocate_cst_act_species_indices(num_mins_cst_eq+num_mins_cst_kin)
                                call this%materials(imtype)%set_indices_solids()
                                exit
                            else
                                !> We check if mineral exists in the chemical system
                                call this%chem_syst%is_mineral_in_chem_syst(mineral,min_flag,min_ind)
                                if (min_flag .eqv. .true.) then
                                    num_mins_loc=num_mins_loc+1
                                    min_indices(num_mins_loc,num_min_zones)=min_ind !> we save index of mineral in chemical system
                                    if (this%chem_syst%minerals(min_ind)%mineral%cst_act_flag.eqv..true. .and. &
                                    min_ind>this%chem_syst%num_minerals_kin) then
                                        num_mins_cst_eq=num_mins_cst_eq+1
                                        num_mins_loc_eq=num_mins_loc_eq+1
                                    else if (this%chem_syst%minerals(min_ind)%mineral%cst_act_flag.eqv..false. .and. &
                                        min_ind>this%chem_syst%num_minerals_kin) then
                                        num_mins_var_eq=num_mins_var_eq+1
                                        num_mins_loc_eq=num_mins_loc_eq+1
                                    else if (this%chem_syst%minerals(min_ind)%mineral%cst_act_flag.eqv..true. .and. &
                                        min_ind<=this%chem_syst%num_minerals_kin) then
                                        num_mins_cst_kin=num_mins_cst_kin+1
                                        num_mins_loc_kin=num_mins_loc_kin+1
                                    else
                                        num_mins_var_kin=num_mins_var_kin+1
                                        num_mins_loc_kin=num_mins_loc_kin+1
                                    end if
                                else
                                    error stop "Mineral not found in chemical system"
                                end if
                            end if
                        end do
                    else
                        exit
                    end if
                !end do
                if (num_min_zones==nmtype) exit
            end do
            rewind(unit)
        !> Second iteration
            num_min_zones=0 !> counter mineral zones
            do
                read(unit,*) label
                if (label=='INITIAL MINERAL ZONES') then
                    read(unit,*) nmtype !> number of mineral zones
                    do
                        read(unit,*) imtype, temp !> index of mineral zone, temperature (C)
                        if (imtype<1 .or. imtype>nmtype) error stop "Mineral zone index out of range"
                        read(unit,*) name !> name of mineral zone
                        read(unit,*) str !> headings
                        if (index(str,'mineral')/=0) then
                            num_min_zones=num_min_zones+1
                            num_mins_loc=0 !> counter minerals in this zone
                            num_mins_loc_kin=0 !> counter kinetic minerals in this zone
                            num_mins_loc_eq=0 !> counter equilibrium minerals in this zone
                            do
                                read(unit,*) mineral%name, vol_frac, react_surf !> react_surf: [m^2 min./m^3 rock]
                                if (mineral%name=='*') then
                                    exit
                                else
                                    num_mins_loc=num_mins_loc+1
                                    if (min_indices(num_mins_loc,num_min_zones)>this%chem_syst%num_minerals_kin) then !> mineral in equilibrium
                                        num_mins_loc_eq=num_mins_loc_eq+1 !> counter equilibrium minerals in this zone
                                        this%materials(imtype)%mineral_zone%ind_min_chem_syst(&
                                            this%materials(imtype)%mineral_zone%num_minerals_kin+num_mins_loc_eq)=&
                                            min_indices(num_mins_loc,num_min_zones) !> we set index of mineral in equilibrium in mineral zone
                                        this%materials(imtype)%reactive_zone%ind_mins_chem_syst(num_mins_loc_eq)=&
                                            min_indices(num_mins_loc,num_min_zones) !> we set index of mineral in equilibrium in reactive zone
                                        this%materials(imtype)%vol_fracts(&
                                            this%materials(imtype)%mineral_zone%num_minerals_kin+num_mins_loc_eq)=vol_frac !> we set volumetric fraction
                                        this%materials(imtype)%react_surfaces(&
                                            this%materials(imtype)%mineral_zone%num_minerals_kin+num_mins_loc_eq)=react_surf !> we set reactive surface
                                        this%materials(imtype)%concentrations(&
                                            this%materials(imtype)%mineral_zone%num_minerals_kin+num_mins_loc_eq)=1d0 !> we assume minerals are pure phases
                                        this%materials(imtype)%activities(&
                                        this%materials(imtype)%mineral_zone%num_minerals_kin+num_mins_loc_eq)=1d0 !> we assume minerals are pure phases
                                    else !> mineral kinetic
                                        num_mins_loc_kin=num_mins_loc_kin+1 !> counter kinetic minerals in this zone
                                        this%materials(imtype)%mineral_zone%ind_min_chem_syst(num_mins_loc_kin)=&
                                            min_indices(num_mins_loc,num_min_zones) !> we set index of kinetic mineral in mineral zone
                                        this%materials(imtype)%vol_fracts(num_mins_loc_kin)=vol_frac !> we set volumetric fraction
                                        this%materials(imtype)%react_surfaces(num_mins_loc_kin)=react_surf !> we set reactive surface
                                        this%materials(imtype)%concentrations(num_mins_loc_kin)=1d0 !> we assume minerals are pure phases
                                        this%materials(imtype)%activities(num_mins_loc_kin)=1d0 !> we assume minerals are pure phases
                                    end if
                                end if
                            end do
                            call this%materials(imtype)%reactive_zone%set_non_flowing_species()
                            call this%materials(imtype)%reactive_zone%set_ind_eq_reacts()
                            call this%materials(imtype)%reactive_zone%set_stoich_mat_react_zone()
                            call this%materials(imtype)%reactive_zone%set_ind_mins_stoich_mat()
                            call this%materials(imtype)%reactive_zone%set_ind_gases_stoich_mat()
                            call this%materials(imtype)%mineral_zone%set_ind_min_Sk()
                            if (num_min_zones==nmtype) then
                                exit
                            end if
                        else
                            exit
                        end if
                    end do
                    exit
                else
                    continue
                end if
            end do
                if (allocated(this%reactive_zones)) then
                    num_rz_old=size(this%reactive_zones)
                    allocate(aux_react_zones(num_rz_old))
                    num_gas_rz=0
                    do i=1,num_rz_old
                        call aux_react_zones(i)%assign_react_zone(this%reactive_zones(i))
                        if (this%reactive_zones(i)%gas_phase%num_gases_eq>0 .and. this%reactive_zones(i)%num_solids==0) then
                            num_gas_rz=num_gas_rz+1
                        end if
                    end do
                    !num_surf_rz=(num_rz_old-num_gas_rz)/(1+num_gas_rz)
                    num_rz=num_gas_rz+nmtype*(1+num_gas_rz)
                    call this%allocate_reactive_zones(num_rz)
                    do i=1,num_gas_rz
                        call this%reactive_zones(i)%assign_react_zone(aux_react_zones(i))
                    end do
                    do i=1,nmtype
                        call this%reactive_zones(num_gas_rz+i)%assign_react_zone(this%materials(i)%reactive_zone)
                    end do
                    do i=1,num_gas_rz
                        !do j=1,num_surf_rz
                        !    call this%reactive_zones(num_gas_rz+nmtype+i*num_surf_rz+j)%assign_react_zone(aux_react_zones(&
                        !    num_gas_rz+num_surf_rz+(i-1)*num_surf_rz+j))
                        !end do
                        do j=1,nmtype
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%assign_react_zone(&
                                this%materials(j)%reactive_zone)
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_gas_phase(&
                                this%reactive_zones(i)%gas_phase)
                            !call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_gas_phase(&
                            !    this%reactive_zones(i)%gas_phase)
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_non_flowing_species()
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_speciation_alg_dimensions()
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_ind_eq_reacts()
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_stoich_mat_react_zone()
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_ind_mins_stoich_mat()
                            call this%reactive_zones(num_rz-num_gas_rz*nmtype+(i-1)*nmtype+j)%set_ind_gases_stoich_mat()
                        end do
                    end do
                else
                    call this%allocate_reactive_zones(nmtype)
                    do i=1,nmtype
                        call this%reactive_zones(i)%assign_react_zone(this%materials(i)%reactive_zone)
                    end do
                end if
    end if
end subroutine read_init_min_zones_CHEPROO