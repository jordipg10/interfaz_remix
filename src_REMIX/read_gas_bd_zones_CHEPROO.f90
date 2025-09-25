subroutine read_gas_bd_zones_CHEPROO(this,unit,gas_bd_zones,reactive_zones)
    use chemistry_m,   only: chemistry_c, gas_chemistry_c, reactive_zone_c, gas_c, gas_phase_c
    implicit none
    class(chemistry_c) :: this
    integer(kind=4), intent(in) :: unit !> file
    type(gas_chemistry_c), intent(out), allocatable :: gas_bd_zones(:)
    type(reactive_zone_c), intent(inout), allocatable, optional :: reactive_zones(:)
    
    integer(kind=4) :: i,j,k,ngtype,igtype,nrwtype,gas_ind,num_gas_zones,num_gases_loc,num_gases_var,num_gases_cst
    integer(kind=4), allocatable :: ind_gases(:)
    character(len=256) :: str,constrain,label
    real(kind=8) :: guess,c_tot,temp,part_press
    type(gas_c) :: gas
    type(gas_phase_c) :: gas_phase
    type(reactive_zone_c) :: react_zone
    logical :: flag

    read(unit,*) ngtype !> number of gas surface boundary zones
    
    allocate(gas_bd_zones(ngtype))
    allocate(ind_gases(this%chem_syst%gas_phase%num_species)) !> chapuza

    num_gas_zones=0 !> counter gas zones
    !num_gases_glob=1 !> counter gases in chemical system
    
    do
        read(unit,*) igtype
        if (igtype<1 .or. igtype>ngtype) error stop
        read(unit,*) str
        if (str=='*') then
            if (num_gas_zones<ngtype) then
                continue
            else
                exit
            end if
        else if (index(str,'gas')/=0) then
            num_gas_zones=num_gas_zones+1
            num_gases_loc=0 !> counter gases in this zone
            num_gases_var=0
            num_gases_cst=0
            do
                read(unit,*) gas%name
                if (gas%name=='*') then
                    call react_zone%set_chem_syst_react_zone(this%chem_syst)
                    call react_zone%gas_phase%allocate_gases(num_gases_loc)
                    call react_zone%gas_phase%set_num_var_act_species_phase(num_gases_var)
                    call react_zone%gas_phase%set_num_cst_act_species_phase(num_gases_cst)
                    call gas_bd_zones(igtype)%set_reactive_zone(react_zone)
                    call gas_bd_zones(igtype)%set_temp()
                    call gas_bd_zones(igtype)%allocate_partial_pressures()
                    call gas_bd_zones(igtype)%allocate_conc_gases()
                    call gas_bd_zones(igtype)%allocate_log_act_coeffs_gases()
                    call gas_bd_zones(igtype)%allocate_var_act_species_indices(num_gases_var)
                    call gas_bd_zones(igtype)%allocate_cst_act_species_indices(num_gases_cst)
                    call gas_bd_zones(igtype)%set_indices_gases()
                    exit
                else
                    call this%chem_syst%gas_phase%is_gas_in_gas_phase(gas,flag,gas_ind)
                    if (flag .eqv. .true.) then
                        num_gases_loc=num_gases_loc+1
                        ind_gases(num_gases_loc)=gas_ind
                        if (this%chem_syst%gas_phase%gases(gas_ind)%cst_act_flag .eqv. .true.) then
                            num_gases_cst=num_gases_cst+1
                        else
                            num_gases_var=num_gases_var+1
                        end if
                    else
                        error stop "Gas not found in gas phase"
                    end if
                end if
            end do
            if (num_gas_zones==ngtype) exit
        else
            exit
        end if
    end do
    rewind(unit)
    num_gas_zones=0 !> counter gas zones
    !num_gases_glob=1 !> counter gases in chemical system
    do
        read(unit,*) label
        if (label=='INITIAL GAS BOUNDARY ZONES') then
            read(unit,*) ngtype !> number of gas zones
            do i=1,ngtype
                read(unit,*) igtype
                read(unit,*) str
                !if (str=='*') then
                !    if (num_gas_zones<ngtype) then
                !        continue
                !    else
                !        exit
                !    end if
                !else if (index(str,'gas')/=0) then
                    num_gas_zones=num_gas_zones+1
                    !num_gases_loc=0 !> counter gases in this zone
                    do j=1,gas_bd_zones(igtype)%reactive_zone%gas_phase%num_species
                        read(unit,*) gas%name, part_press
                        !if (gas%name=='*') then
                        !    call gas_bd_zones(igtype)%compute_conc_gases_ideal() !> we assume gases are ideal
                        !    exit
                        !else if (gas%name/=this%solid_chemistry%reactive_zone%chem_syst%gas_phase%gases(num_gases_glob)%name) then
                        !    error stop
                        !else
                            !num_gases_loc=num_gases_loc+1
                            gas_bd_zones(igtype)%reactive_zone%gas_phase%gases(j)=this%chem_syst%gas_phase%gases(ind_gases(j)) !> we set gas
                            gas_bd_zones(igtype)%activities(j)=part_press !> activities are partial pressures
                        !end if
                    end do
                    call gas_bd_zones(igtype)%compute_conc_gases_ideal() !> we assume gases are ideal
                    call gas_bd_zones(igtype)%compute_log_act_coeffs_gases()
                    call gas_bd_zones(igtype)%compute_pressure()
                !else
                !    exit
                !end if
            end do
        else if (label=='end') then
            backspace(unit)
            exit
        else
            continue
        end if
    end do
    if (present(reactive_zones)) then
        allocate(reactive_zones(ngtype))
        do i=1,ngtype
            reactive_zones(i)=gas_bd_zones(i)%reactive_zone
        end do
    end if
end subroutine