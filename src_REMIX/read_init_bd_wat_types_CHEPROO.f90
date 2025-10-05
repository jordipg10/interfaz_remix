subroutine read_init_bd_wat_types_CHEPROO(this,unit,init_cat_exch_zones,&
        gas_chem)
    use chemistry_m, only: chemistry_c
    use aqueous_chemistry_m, only: aqueous_chemistry_c, aq_species_c, aq_phase_c, gas_chemistry_c, & 
        reactive_zone_c, mineral_zone_c, solid_chemistry_c, mineral_c, chem_system_c
    implicit none
    class(chemistry_c) :: this !> chemistry object
    integer(kind=4), intent(in) :: unit !> file unit
    type(solid_chemistry_c), intent(inout), allocatable :: init_cat_exch_zones(:)
    type(gas_chemistry_c), intent(in), optional :: gas_chem !> chapuza
    
    integer(kind=4) :: i,j,k,l,nwtype,icon,n_p_aq,gas_ind,min_ind,model,niter,sp_ind
    integer(kind=4), allocatable :: cols(:)
    integer(kind=4), allocatable :: num_aq_prim_array(:)
    integer(kind=4), allocatable :: num_cstr_array(:)
    character(len=256) :: prim_sp_name,constrain,label,name
    real(kind=8) :: guess,c_tot,temp,conc
    real(kind=8), allocatable :: eq_csts(:)
    logical :: CV_flag,flag,flag_surf,flag_comp,flag_Se
    
    type(reactive_zone_c) :: react_zone !> default object
    type(reactive_zone_c), allocatable :: react_zones(:) !> default objects
    type(mineral_zone_c), allocatable :: min_zones(:) !> default objects
    type(solid_chemistry_c) :: solid_chem !> default object
    !type(solid_chemistry_c), allocatable :: solid_chems(:) !> default objects
    type(solid_chemistry_c), allocatable :: cat_exch_zones(:) !> 
    type(aq_species_c) :: aq_species
    type(mineral_c) :: mineral
    type(aq_phase_c) :: old_aq_phase
    

    
    read(unit,*) this%act_coeffs_model
    
    read(unit,*) nwtype
    
    if (size(init_cat_exch_zones)==1) then
        allocate(cat_exch_zones(nwtype)) !> we allocate cat exch zones
        do i=1,nwtype
            call cat_exch_zones(i)%assign_solid_chemistry(init_cat_exch_zones(1)) !> chapuza
        end do
        deallocate(init_cat_exch_zones)
        allocate(init_cat_exch_zones(nwtype))
    end if
        
    call this%allocate_wat_types(nwtype) !> we allocate water types
    call this%allocate_sol_types() !> we allocate solid types associated to water types (chapuza)
    call this%allocate_react_zones_wat_types() !> we allocate reactive zones associated to water types (chapuza)
    call this%allocate_gas_zones_wat_types() !> we allocate gas zones associated to water types (chapuza)
    allocate(react_zones(nwtype),min_zones(nwtype))
    allocate(cols(2))
    allocate(num_aq_prim_array(nwtype),num_cstr_array(nwtype))
    
    num_aq_prim_array=0
    num_cstr_array=0
    
     do i=1,nwtype
        read(unit,*) j, temp !> we read index water type and temperature (in Celsius)
        if (j<1 .or. j>nwtype) error stop
        !ind_wat_type(i)=j
        call this%wat_types(j)%set_aq_phase(this%chem_syst%aq_phase)
        call this%wat_types(j)%set_indices_aq_species_aq_chem() !> we set default indices
        !call this%wat_types(j)%set_ind_diss_solids_aq_chem() !> we set default indices dissolved solids
        call this%wat_types(j)%set_temp(temp+273.18) !> Kelvin
        call this%wat_types(j)%set_density()
        !call this%wat_types(j)%set_solid_chemistry(solid_chem)
        call this%wat_types(j)%allocate_conc_aq_species()
        call this%wat_types(j)%allocate_log_act_coeffs_aq_chem()
        call this%wat_types(j)%allocate_activities_aq_species()
        read(unit,*) name !> we read name of water type
        call this%wat_types(j)%set_name(trim(name))
        read(unit,*) label
        if (index(label,'icon')/=0) then !> 'icon, guess, ctot, constrain'
            k=0 !> counter primary aqueous species
            l=0 !> counter dissolved solids
            do
                read(unit,*) aq_species%name, icon!, guess, ctot, constrain%name
                if (aq_species%name=='*') then
                    exit
                else
                    call this%wat_types(j)%aq_phase%is_species_in_aq_phase(aq_species,flag,sp_ind)
                    if (flag .eqv. .true.) then
                        k=k+1
                        if (sp_ind/=k) then
                            this%wat_types(j)%indices_aq_species(sp_ind)=k
                            this%wat_types(j)%indices_aq_phase(k)=sp_ind
                        end if
                        ! if (aq_species%name=='h2o') then
                        !     !call this%wat_types(j)%set_ind_wat_aq_chem(k)
                        ! else
                        !     l=l+1
                        !     !this%wat_types(j)%ind_diss_solids(l)=k
                        !     ! if (aq_species%name=='h+') then
                        !     !     call this%wat_types(j)%set_ind_prot_aq_chem(k)
                        !     ! !> aqui faltan especies
                        !     ! end if
                        ! end if
                        num_aq_prim_array(j)=num_aq_prim_array(j)+1
                        if (icon==4) then
                            num_cstr_array(j)=num_cstr_array(j)+1
                        !if (icon==1) then
                        !    n_icons(1)=n_icons(1)+1
                        !else if (icon==2) then
                        !    n_icons(2)=n_icons(2)+1
                        !    !n_aq_comp=n_aq_comp+1
                        !else if (icon==3) then
                        !    n_icons(3)=n_icons(3)+1
                        !    if (aq_species%name=='h+') then
                        !        call this%aq_phase%set_ind_proton(sp_ind)
                        !        call this%set_pH(-log10(ctot))
                        !    end if
                        !else if (icon==4) then
                        !    n_icons(4)=n_icons(4)+1 !> number of constrains
                        !    !num_cstr=num_cstr+1
                        !    !call this%chem_syst%is_eq_reaction_in_chem_syst(constrain,flag,ind_cstr)
                        !    !if (flag==.true.) then
                        !    !    indices_constrains(l)=ind_cstr
                        !    !    l=l+1
                        !    !else
                        !    !    error stop
                        !    !end if
                        !!else if (icon==5) then
                        !!    n_icons(6)=n_icons(6)+1
                        else if (icon<1 .or. icon>4) then
                            error stop "icon option not implemented yet"
                        end if
                        !this%concentrations(sp_ind)=guess
                        !icons(sp_ind)=icon
                        !ctots(sp_ind)=ctot
                        !constrains(sp_ind)=constrain
                        !k=k+1 !> aqui hay que verificar dimension
                    else
                        error stop "Aqueous species not present in aqueous phase"
                    end if
                end if
            end do
        else
            error stop "Error reading water type"
        end if
    end do
    rewind(unit)
    do
        read(unit,*) label
        if (label=='INITIAL AND BOUNDARY WATER TYPES') then
            read(unit,*) model
            read(unit,*) nwtype
            do i=1,nwtype
                read(unit,*) j, temp !> we read index water type and temperature (in Celsius)
                read(unit,*) name
                call this%react_zones_wat_types(j)%set_CV_params(this%CV_params)
                call this%react_zones_wat_types(j)%set_chem_syst_react_zone(this%chem_syst)
                call min_zones(j)%set_chem_syst_min_zone(this%chem_syst)
                call this%sol_types(j)%set_reactive_zone(this%react_zones_wat_types(j)) !> we set default reactive zone
                call this%sol_types(j)%set_mineral_zone(min_zones(j)) !> we set default mineral zone
                ! if (present(gas_chem)) then
                !    call this%wat_types(j)%set_gas_chemistry(gas_chem) !> chapuza
                ! end if
                if (size(init_cat_exch_zones)>1) then
                    call this%wat_types(j)%set_solid_chemistry(cat_exch_zones(j)) !> chapuza
                    call this%wat_types(j)%read_wat_type_CHEPROO(num_aq_prim_array(j),num_cstr_array(j),this%num_gas_zones,this%act_coeffs_model,&
                        this%Jac_opt,unit,niter,CV_flag)
                    call init_cat_exch_zones(j)%assign_solid_chemistry(this%wat_types(j)%solid_chemistry) !> chapuza
                else if (present(gas_chem)) then
                    call this%wat_types(j)%set_solid_chemistry(this%sol_types(j)) !> chapuza
                    call this%wat_types(j)%read_wat_type_CHEPROO(num_aq_prim_array(j),num_cstr_array(j),this%num_gas_zones,this%act_coeffs_model,&
                        this%Jac_opt,unit,niter,CV_flag,gas_chem)
                end if
            end do
            exit
         else
            continue
         end if
    end do
    !> we eliminate constant activity species from component matrix of chemical system
    call this%chem_syst%change_spec_alg_chem_syst(this%CV_params%zero)
    !> we eliminate constant activity species from component matrix of water types
    do i=1,nwtype
        call this%wat_types(i)%change_spec_alg_aq_chem()
    end do
!> Chapuza
    ! do i=1,nwtype
    !     !> rearrange indices
    !     ! print *, this%wat_types(i)%indices_aq_species
    !     ! print *, this%wat_types(i)%ind_prim_species
    !     ! print *, this%wat_types(i)%ind_diss_solids
    !     ! call this%wat_types(i)%set_ind_prim_sec_species()
    !     !call this%wat_types(i)%set_indices_aq_species_aq_chem()
    !     !call this%wat_types(i)%set_ind_diss_solids_aq_chem()
    ! end do
    !do i=1,this%num_bd_this%wat_types
    !    call this%bd_this%wat_types(i)%rearrange_state_vars(old_aq_phase)
    !end do
!> Chapuza
    do i=1,this%chem_syst%num_minerals_kin
        !> indices reactants
        call this%chem_syst%min_kin_reacts(i)%set_indices_aq_phase_min(this%chem_syst%aq_phase)
    end do
!> Chapuza
    do i=1,this%chem_syst%num_lin_kin_reacts
        !> indices reactants
        call this%chem_syst%lin_kin_reacts(i)%set_index_aq_phase_lin(this%chem_syst%aq_phase)
    end do
end subroutine