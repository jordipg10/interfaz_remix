subroutine read_init_gas_zones_CHEPROO(this,unit,ngrz)
    use chemistry_m, only: chemistry_c, gas_chemistry_c, reactive_zone_c, gas_c, gas_phase_c
    use gas_phase_m, only: are_gas_phases_equal
    implicit none
    class(chemistry_c) :: this
    integer(kind=4), intent(in) :: unit !> file
    !type(gas_chemistry_c), intent(out), allocatable :: this%gas_zones(:)
    integer(kind=4), intent(out) :: ngrz !> number of gas reactive zones
    !type(reactive_zone_c), intent(inout), allocatable, optional :: reactive_zones(:)
    
    integer(kind=4) :: num_rz,i,j,k,ngtype,igtype,nrwtype,gas_ind,num_gases_loc,num_gases_var,num_gases_cst,n_gas_eq,&
    n_gas_kin,num_surf_rz,n_gas_eq_cst_act,n_gas_eq_var_act
    integer(kind=4), allocatable :: ind_gases(:) !> indices of gases in gas phase of chemical system
    integer(kind=4), allocatable :: ind_gas_zones(:) !> indices of gases in gas phase of chemical system
    character(len=256) :: str,constrain,label,name
    real(kind=8) :: guess,conc,temp,part_press,vol
    type(gas_c) :: gas
    type(gas_phase_c) :: gas_phase
    type(reactive_zone_c) :: react_zone
    type(reactive_zone_c), allocatable :: aux_react_zones(:),react_zones(:)
    logical :: flag_comp,flag

    read(unit,*) ngtype !> number of gas zones
    if (ngtype<0) then
        error stop "Number of gas zones cannot be negative"
    else if (ngtype==0) then
        ngrz=0
        return
    end if
    
    !allocate(this%gas_zones(ngtype))
    call this%allocate_gas_zones(ngtype)
    allocate(ind_gas_zones(ngtype))
    allocate(ind_gases(this%chem_syst%gas_phase%num_species))
    allocate(react_zones(ngtype))

    this%num_gas_zones=0 !> counter gas zones
    
    ngrz=1 !> counter gas reactive zones

    flag_comp=.true. !> flag to compute speciation algebra dimensions (by default: TRUE)
    
    do
        read(unit,*) igtype, temp, vol !> index, temperature and volume of gas zone
        if (igtype<1 .or. igtype>ngtype) error stop "Index of gas zone out of range"
        read(unit,*) name !> name of gas zone
        call this%gas_zones(igtype)%set_name(trim(name)) !> we set name of gas zone
        read(unit,*) str !> headings
        if (str=='*') then
            if (this%num_gas_zones<ngtype) then
                continue
            else
                exit
            end if
        else if (index(str,'gas')/=0) then
            this%num_gas_zones=this%num_gas_zones+1 !> counter gas zones
            ind_gas_zones(this%num_gas_zones)=igtype !> we store index of gas zone
            num_gases_loc=0 !> counter gases in this zone
            num_gases_var=0
            num_gases_cst=0
            n_gas_eq=0
            n_gas_eq_cst_act=0
            n_gas_eq_var_act=0
            n_gas_kin=0
            do
                read(unit,*) gas%name, part_press
                if (gas%name=='*') then
                    call react_zones(igtype)%set_chem_syst_react_zone(this%chem_syst)
                    call react_zones(igtype)%set_CV_params(this%CV_params)
                    call react_zones(igtype)%gas_phase%allocate_gases(num_gases_loc)
                    call react_zones(igtype)%gas_phase%set_num_var_act_species_phase(num_gases_var)
                    call react_zones(igtype)%gas_phase%set_num_cst_act_species_phase(num_gases_cst)
                    call react_zones(igtype)%gas_phase%compute_num_species_phase()
                    call react_zones(igtype)%gas_phase%set_num_gases_eq(n_gas_eq) !> 
                    call react_zones(igtype)%gas_phase%set_num_gases_eq_cst_act(n_gas_eq_cst_act) !> 
                    call react_zones(igtype)%gas_phase%set_num_gases_eq_var_act(n_gas_eq_var_act) !> 
                    call react_zones(igtype)%gas_phase%set_num_gases_kin(n_gas_kin) !> 
                    call react_zones(igtype)%set_num_solids() !> sets number of solids in reactive zone
                    call react_zones(igtype)%set_num_non_flowing_species() !> sets number of non-flowing species in reactive zone
                    call react_zones(igtype)%set_speciation_alg_dimensions(flag_comp)
                    !if (this%num_gas_zones>1) then
                    !    call are_gas_phases_equal(react_zones(igtype)%gas_phase,react_zones(ind_this%gas_zones(this%num_gas_zones-1))%gas_phase,flag)
                    !    if (flag==.true.) then
                    !        call this%gas_zones(igtype)%set_reactive_zone(react_zones(ind_this%gas_zones(this%num_gas_zones-1)))
                    !    else
                    !        call this%gas_zones(igtype)%set_reactive_zone(react_zones(igtype))
                    !    end if
                    !end if
                    call this%gas_zones(igtype)%set_reactive_zone(react_zones(igtype)) !> we set reactive zone (chapuza)
                    call this%gas_zones(igtype)%set_temp(temp+273.15) !> Kelvin
                    call this%gas_zones(igtype)%set_volume(vol)
                    call this%gas_zones(igtype)%allocate_partial_pressures()
                    call this%gas_zones(igtype)%allocate_conc_gases()
                    call this%gas_zones(igtype)%allocate_log_act_coeffs_gases()
                    call this%gas_zones(igtype)%allocate_var_act_species_indices(num_gases_var)
                    call this%gas_zones(igtype)%allocate_cst_act_species_indices(num_gases_cst)
                    call this%gas_zones(igtype)%allocate_ind_gases_eq_cst_act()
                    exit
                else
                    call this%chem_syst%gas_phase%is_gas_in_gas_phase(gas,flag,gas_ind)
                    if (flag .eqv. .true.) then
                        num_gases_loc=num_gases_loc+1
                        ind_gases(num_gases_loc)=gas_ind !> we store index of gas in gas phase of chemical system
                        if (this%chem_syst%gas_phase%gases(gas_ind)%cst_act_flag.eqv. .true. .and. gas_ind<=&
                        this%chem_syst%gas_phase%num_gases_eq) then
                            num_gases_cst=num_gases_cst+1
                            n_gas_eq=n_gas_eq+1
                            n_gas_eq_cst_act=n_gas_eq_cst_act+1
                        !> we modify equilibrium constant
                            !this%chem_syst%eq_reacts(this%chem_syst%aq_phase%num_aq_complexes+&
                            !this%chem_syst%num_minerals_eq+gas_ind)%eq_cst=this%chem_syst%eq_reacts(&
                            !this%chem_syst%aq_phase%num_aq_complexes+this%chem_syst%num_minerals_eq+gas_ind)%eq_cst/part_press
                        else if (this%chem_syst%gas_phase%gases(gas_ind)%cst_act_flag.eqv. .false. .and. gas_ind<=&
                            this%chem_syst%gas_phase%num_gases_eq) then
                            num_gases_var=num_gases_var+1
                            n_gas_eq=n_gas_eq+1
                            n_gas_eq_var_act=n_gas_eq_var_act+1
                        else if (this%chem_syst%gas_phase%gases(gas_ind)%cst_act_flag.eqv. .true. .and. gas_ind>&
                            this%chem_syst%gas_phase%num_gases_eq) then
                            num_gases_cst=num_gases_cst+1
                            n_gas_kin=n_gas_kin+1
                        else
                            num_gases_var=num_gases_var+1
                            n_gas_kin=n_gas_kin+1
                        end if
                        !if (gas_ind<=this%chem_syst%gas_phase%num_gases_eq) then
                        !    n_gas_eq=n_gas_eq+1
                        !else
                        !    n_gas_kin=n_gas_kin+1
                        !end if
                    else
                        error stop "Gas not found in gas phase"
                    end if
                end if
            end do
            if (this%num_gas_zones==ngtype) exit
        else
            exit
        end if
    end do
    rewind(unit)
    this%num_gas_zones=0 !> counter gas zones
    !num_gases_glob=1 !> counter gases in chemical system
    do
        read(unit,*) label
        if (label=='INITIAL AND BOUNDARY GAS ZONES') then
            read(unit,*) ngtype !> number of gas zones
            do i=1,ngtype
                read(unit,*) igtype, temp, vol !> index, temperature and volume of gas zone
                read(unit,*) name !> name of gas zone
                read(unit,*) str !> headings
                this%num_gas_zones=this%num_gas_zones+1 !> counter gas zones
                k=0 !> counter gases in equilibrium with constant activity in this zone
                do j=1,this%gas_zones(igtype)%reactive_zone%gas_phase%num_species
                    read(unit,*) gas%name, part_press
                    call this%gas_zones(igtype)%reactive_zone%gas_phase%gases(j)%assign_species(this%chem_syst%gas_phase%gases(&
                        ind_gases(j))) !> we set gas
                    this%gas_zones(igtype)%activities(j)=part_press !> activities are partial pressures
                    if (ind_gases(j) <= this%chem_syst%gas_phase%num_gases_eq_cst_act) then
                        k=k+1
                        this%gas_zones(igtype)%ind_gases_eq_cst_act(k)=ind_gases(j) !> we store index of gas in gas phase of chemical system
                    end if
                end do
                call this%gas_zones(igtype)%compute_conc_gases_ideal() !> we assume gases are ideal
                call this%gas_zones(igtype)%compute_pressure() !> we compute pressure
                call this%gas_zones(igtype)%compute_log_act_coeffs_gases() !> we compute log10 activity coefficients
                call this%gas_zones(igtype)%set_indices_gases() !> we set indices of gases in gas phase of chemical system
                if (this%num_gas_zones>1) then
                    call are_gas_phases_equal(this%gas_zones(igtype)%reactive_zone%gas_phase,this%gas_zones(&
                        ind_gas_zones(this%num_gas_zones-1))%reactive_zone%gas_phase,flag)
                    if (flag .eqv. .true.) then
                        call this%gas_zones(igtype)%set_reactive_zone(this%gas_zones(ind_gas_zones(&
                            this%num_gas_zones-1))%reactive_zone)
                    else
                        !call this%gas_zones(igtype)%set_reactive_zone(react_zones(igtype))
                        ngrz=ngrz+1 !> counter gas reactive zones
                        call this%gas_zones(igtype)%reactive_zone%set_non_flowing_species() !> sets non-flowing species in reactive zone
                        call this%gas_zones(igtype)%reactive_zone%set_ind_eq_reacts() !> sets indices of equilibrium reactions in reactive zone
                        call this%gas_zones(igtype)%reactive_zone%set_stoich_mat_react_zone() !> sets stoichiometric matrix of reactive zone
                        call this%gas_zones(igtype)%reactive_zone%set_ind_gases_stoich_mat() !> sets indices of gases in stoichiometric matrix of reactive zone
                    end if
                else
                    call this%gas_zones(igtype)%reactive_zone%set_non_flowing_species() !> sets non-flowing species in reactive zone
                    call this%gas_zones(igtype)%reactive_zone%set_ind_eq_reacts() !> sets indices of equilibrium reactions in reactive zone
                    call this%gas_zones(igtype)%reactive_zone%set_stoich_mat_react_zone() !> sets stoichiometric matrix of reactive zone
                    call this%gas_zones(igtype)%reactive_zone%set_ind_gases_stoich_mat() !> sets indices of gases in stoichiometric matrix of reactive zone
                end if
                read(unit,*) str
            end do
        else if (label=='end') then
            backspace(unit)
            exit
        else
            continue
        end if
    end do
    !if (present(reactive_zones)) then
        if (allocated(this%reactive_zones)) then
            !> First gases, then surfaces
            num_surf_rz=size(this%reactive_zones)
            allocate(aux_react_zones(num_surf_rz))
            do i=1,num_surf_rz
                call aux_react_zones(i)%assign_react_zone(this%reactive_zones(i))
            end do
            num_rz=num_surf_rz+ngtype*(1+num_surf_rz)
            !deallocate(this%reactive_zones)
            !allocate(this%reactive_zones(num_rz))
            call this%allocate_reactive_zones(num_rz)
            do i=1,num_surf_rz
                call this%reactive_zones(ngrz+i)%assign_react_zone(aux_react_zones(i))
            end do
            do i=1,ngrz
                call this%reactive_zones(i)%set_chem_syst_react_zone(this%gas_zones(i)%reactive_zone%chem_syst)
                call this%reactive_zones(i)%set_gas_phase(this%gas_zones(i)%reactive_zone%gas_phase)
                !call this%gas_zones(i)%set_reactive_zone(this%reactive_zones(i))
            end do
            do i=1,ngrz
                do j=1,num_surf_rz
                    call this%reactive_zones(ngrz+i*num_surf_rz+j)%set_chem_syst_react_zone(&
                        this%gas_zones(I)%reactive_zone%chem_syst)
                    call this%reactive_zones(ngrz+i*num_surf_rz+j)%set_gas_phase(this%gas_zones(i)%reactive_zone%gas_phase)
                    call this%reactive_zones(ngrz+i*num_surf_rz+j)%set_cat_exch_zone(aux_react_zones(i)%cat_exch_zone)
                end do
            end do
        else
            !allocate(this%reactive_zones(ngrz))
            call this%allocate_reactive_zones(ngrz)
            do i=1,ngrz
                call this%reactive_zones(i)%assign_react_zone(this%gas_zones(i)%reactive_zone)
                !call this%gas_zones(i)%set_reactive_zone(this%reactive_zones(i)) !> chapuza
            end do
        end if
    !end if
end subroutine