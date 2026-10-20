!> Lectura quimica CHEPROO
subroutine read_chemistry_CHEPROO(this,root,path_pb,path_DB,unit_chem_syst_file,unit_loc_chem_file)
    use chemistry_m, only: chemistry_c
    use solid_chemistry_m, only: solid_chemistry_c
    implicit none
    class(chemistry_c) :: this
    character(len=*), intent(in) :: root
    character(len=*), intent(in) :: path_pb
    character(len=*), intent(in) :: path_DB
    integer(kind=4), intent(in) :: unit_chem_syst_file
    integer(kind=4), intent(in) :: unit_loc_chem_file
    
    integer(kind=4) :: i,ngrz,nmrz,nsrz,ndrz
    integer(kind=4), allocatable :: ind_wat_type(:),num_aq_prim_array(:),num_cstr_array(:)
    character(len=256) :: label
    
    !type(aqueous_chemistry_c), allocatable :: wat_types(:),bd_wat_types(:),rech_wat_types(:)
    !type(solid_chemistry_c), allocatable :: init_min_zones(:),init_cat_exch_zones(:),init_sol_zones(:)
    !type(gas_chemistry_c), allocatable :: init_gas_zones(:)
    
    !type(species_c) :: species,water,constrain
    !type(species_c), allocatable :: aq_species(:),cst_act_species(:),prim_species(:)
    !type(mineral_c) :: mineral
    !class(mineral_c), allocatable :: minerals(:)
    !type(reactive_zone_c), allocatable :: reactive_zones(:),rz_mins(:),rz_surf(:)
    !class(kin_params_c), pointer :: p_kin_params=>null()
    !type(Monod_params_c), allocatable, target :: Monod_params_array(:)
    !type(Monod_params_c), target :: Monod_params
    !class(kin_reaction_c), pointer :: p_kin_react=>null()
    !type(kin_reaction_poly_c) :: kin_react_ptr
    !class(kin_reaction_poly_c), allocatable :: kin_reacts(:)
    !type(redox_kin_c), target :: Monod
    !type(eq_reaction_c) :: eq_react
    !type(eq_reaction_c), allocatable :: eq_reacts(:)
    logical :: flag_comp
    !type(kin_mineral_c), target :: kin_mineral
    
    !class(chem_system_c), pointer :: p_chem_syst=>null()
    
    !type(aq_phase_c), target :: aq_phase_new
    
    !type(surface_c), allocatable :: init_cat_exch_zones(:)
    
    flag_comp=.true.
    ngrz=0 !> number of gas reactive zones by default
    nmrz=0 !> number of mineral reactive zones by default
    ndrz=0 !> number of cation exchange reactive zones by default
!> Chemical system
    open(unit_chem_syst_file,file=path_pb//root//'_sist_quim.dat',status='old',action='read')
    call this%chem_syst%read_chem_system_CHEPROO(path_DB,unit_chem_syst_file)
    close(unit_chem_syst_file)
!> Chapuza
    call this%min_zone_dummy%set_chem_syst_min_zone(this%chem_syst) !> we set the chemical system for the dummy mineral zone
!> Local chemistry
    open(unit_loc_chem_file,file=path_pb//root//'_quim_loc.dat',status='old',action='read')
!> Read the whole _quim_loc.dat file in a single pass, in the order the sections
!> appear in the file: water types first, then mineral zones, then surface
!> adsorption zones, then gas zones. Water types no longer depend on the gas or
!> solid zones (they build their own reactive zone from the icon=4 constraints),
!> so the previous two-pass (rewind) scheme is no longer needed.
    do
        read(unit_loc_chem_file,*) label
        if (label=='end') then
            exit
        else if (label=='INITIAL AND BOUNDARY WATER TYPES') then
            call this%read_init_bd_wat_types_CHEPROO(unit_loc_chem_file)
        else if (label=='INITIAL MINERAL ZONES') then
            call this%read_init_min_zones_CHEPROO(unit_loc_chem_file,nmrz)
            write(0,'(A,I0)') ' [RCC-DBG] min_zones done nmrz=',nmrz
            block; real(kind=8),allocatable::hpq(:); integer::ihpq
                do ihpq=1,200; allocate(hpq(1024*ihpq)); hpq=1d0; deallocate(hpq); end do
                write(0,'(A)') ' [HPC] post read_init_min_zones'
            end block
        else if (label=='INITIAL SURFACE ADSORPTION ZONES') then
            call this%read_init_cat_exch_zones_CHEPROO(unit_loc_chem_file,ndrz)
            write(0,'(A,I0)') ' [RCC-DBG] cat_exch done ndrz=',ndrz
            block; real(kind=8),allocatable::hpq(:); integer::ihpq
                do ihpq=1,200; allocate(hpq(1024*ihpq)); hpq=1d0; deallocate(hpq); end do
                write(0,'(A)') ' [HPC] post read_init_cat_exch'
            end block
        else if (label=='INITIAL AND BOUNDARY GAS ZONES') then
            call this%read_gas_zones_CHEPROO(unit_loc_chem_file,ngrz)
        else
            continue
        end if
    end do
    close(unit_loc_chem_file)
    write(0,'(A)') ' [RCC-DBG] read loop done (pre materials block)'
    !do i=1,size(init_cat_exch_zones)
    !    call init_cat_exch_zones(i)%set_reactive_zone(this%reactive_zones(size(init_gas_zones)+1))
    !end do
!!> Chapuza
!    if (size(init_cat_exch_zones)>0) then
!        allocate(init_cat_exch_zones(this%num_init_wat_types))
!        if (size(init_cat_exch_zones)==1) then
!            do i=1,this%num_init_wat_types
!                init_cat_exch_zones(ind_wat_type(i))=init_cat_exch_zones(1)
!                call this%init_wat_types(ind_wat_type(i))%read_wat_type_CHEPROO(num_aq_prim_array(ind_wat_type(i)),num_cstr_array(ind_wat_type(i)),this%model,this%Jac_opt,unit,niter,CV_flag,init_cat_exch_zones(ind_wat_type(i))%solid_chem)
!                !init_cat_exch_zones(ind_wat_type(i))%solid_chem%concentrations(1)=conc_exch(num_aq_prim_array(ind_wat_type(i))+1)
!                !init_cat_exch_zones(ind_wat_type(i))%solid_chem%concentrations(2:init_cat_exch_zones(ind_wat_type(i))%solid_chem%reactive_zone%cat_exch_zone%num_surf_compl)=conc_exch(this%init_wat_types(ind_wat_type(i))%aq_chem%aq_phase%num_species+2:this%init_wat_types(ind_wat_type(i))%aq_chem%solid_chemistry%reactive_zone%speciation_alg%num_species)
!            end do
!            do i=1,this%num_bd_wat_types
!                init_cat_exch_zones(ind_wat_type(this%num_init_wat_types+i))=init_cat_exch_zones(1)
!                call this%bd_wat_types(ind_wat_type(this%num_init_wat_types+i))%read_wat_type_CHEPROO(num_aq_prim_array(ind_wat_type(this%num_init_wat_types+i)),num_cstr_array(ind_wat_type(this%num_init_wat_types+i)),this%model,this%Jac_opt,unit,niter,CV_flag,init_cat_exch_zones(ind_wat_type( i))%solid_chem)
!                !init_cat_exch_zones(ind_wat_type(this%num_init_wat_types+i))%solid_chem%concentrations(1)=conc_exch(num_aq_prim_array(ind_wat_type(i))+1)
!                !init_cat_exch_zones(ind_wat_type(i))%solid_chem%concentrations(2:init_cat_exch_zones(ind_wat_type(i))%solid_chem%reactive_zone%cat_exch_zone%num_surf_compl)=conc_exch(this%init_wat_types(ind_wat_type(i))%aq_chem%aq_phase%num_species+2:this%init_wat_types(ind_wat_type(i))%aq_chem%solid_chemistry%reactive_zone%speciation_alg%num_species)
!            end do
!        end if
!    else
!        do i=1,this%num_init_wat_types
!            call this%init_wat_types(ind_wat_type(i))%read_wat_type_CHEPROO(num_aq_prim_array(ind_wat_type(i)),num_cstr_array(ind_wat_type(i)),this%model,this%Jac_opt,unit,niter,CV_flag)
!        end do
!        do i=1,this%num_bd_wat_types
!            call this%bd_wat_types(ind_wat_type(this%num_init_wat_types+i))%read_wat_type_CHEPROO(num_aq_prim_array(ind_wat_type(this%num_init_wat_types+i)),num_cstr_array(ind_wat_type(this%num_init_wat_types+i)),this%model,this%Jac_opt,unit,niter,CV_flag)
!        end do
!    end if
!> Chapuza 
    if (ndrz>0 .and. nmrz==0) then
    !    allocate(init_sol_zones(size(init_min_zones)+size(init_cat_exch_zones)))
    !    do i=1,size(init_min_zones)
    !        init_sol_zones(size(init_cat_exch_zones)+i)=init_min_zones(i)
    !    end do
        call this%set_num_materials(this%num_init_cat_exch_zones)
        call this%allocate_materials()
       do i=1,this%num_init_cat_exch_zones
           call this%materials(i)%copy_solid_chemistry(this%init_cat_exch_zones(i))
       end do
    else if (ndrz>0 .and. nmrz>0) then
        !> Solid zones are ordered [mineral-only zones (1..nmrz), combined mineral+exchange
        !> zones (nmrz+1..nmrz+ndrz)], matching the reactive-zone ordering built in
        !> read_init_cat_exch_zones_CHEPROO. Preserve the mineral-only materials and APPEND the
        !> combined solid chemistries (this%init_cat_exch_zones, built merged there), so
        !> <root>_tar_sol.dat can select the calcite-only zone (index<=nmrz) or the combined
        !> calcite+exchange zone (index>nmrz). The combined chemistries' mineral_zone pointers
        !> target the persistent this%mineral_zones array (not this%materials), so the realloc
        !> is safe; read_tar_sol repoints each target solid's reactive_zone to the live zone.
        block
            type(solid_chemistry_c), allocatable :: aux_materials(:)
            integer(kind=4) :: im, nmat0
            nmat0 = this%num_materials
            allocate(aux_materials(nmat0))
            do im=1,nmat0
                call aux_materials(im)%copy_solid_chemistry(this%materials(im))
            end do
            deallocate(this%materials)
            call this%set_num_materials(nmat0 + this%num_init_cat_exch_zones)
            call this%allocate_materials()
            do im=1,nmat0
                call this%materials(im)%copy_solid_chemistry(aux_materials(im))
            end do
            do im=1,this%num_init_cat_exch_zones
                call this%materials(nmat0+im)%copy_solid_chemistry(this%init_cat_exch_zones(im))
            end do
            deallocate(aux_materials)
        end block
    end if
    ! if (allocated(init_cat_exch_zones) .and. (.not. allocated(init_min_zones))) then
    !    allocate(init_sol_zones(size(init_cat_exch_zones)))
    !    do i=1,size(init_cat_exch_zones)
    !        init_sol_zones(i)=init_cat_exch_zones(i)
    !    end do
    !    !nsrz=size(init_cat_exch_zones)
    ! else if (nmrz>0 .and. (.not. allocated(init_cat_exch_zones))) then
    !     !allocate(init_sol_zones(size(init_min_zones)))
    !     !do i=1,nmrz
    !     !    init_sol_zones(i)=this%materials(i)
    !     !end do
    !     nsrz=nmrz
    ! else if (allocated(init_cat_exch_zones) .and. allocated(init_min_zones)) then
    !     allocate(init_sol_zones(size(init_min_zones)))
    !     do i=1,size(init_min_zones)
    !         init_sol_zones(i)=init_min_zones(i)
    !     end do
    !     nsrz=nmrz
    ! else
    !     allocate(init_sol_zones(0)) !> Chapuza
    !     nsrz=0
    ! end if  

!> total number of solid reactive zones
    nsrz=nmrz+ndrz
    write(0,'(A)') ' [RCC-DBG] materials block done'
    block; real(kind=8),allocatable::hpq(:); integer::ihpq
        do ihpq=1,200; allocate(hpq(1024*ihpq)); hpq=1d0; deallocate(hpq); end do
        write(0,'(A)') ' [HPC] post materials block'
    end block
    write(0,'(A,I0,A,I0,A,I0)') ' MARK A: sections read. nmrz=',nmrz,' ndrz=',ndrz,' nsrz=',nsrz
!!> Chapuza
!    do i=1,size(reactive_zones)
!        call reactive_zones(i)%set_num_solids()
!        call reactive_zones(i)%set_ind_non_flow_species()
!        call reactive_zones(i)%set_speciation_alg_dimensions(flag_comp)
!        call reactive_zones(i)%set_eq_reactions()
!        call reactive_zones(i)%set_stoich_mat_react_zone()
!        call reactive_zones(i)%set_stoich_mat_sol_rz()
!        call reactive_zones(i)%set_CV_params(this%CV_params) !> chapuza
!    end do
!> We verify pointers (chapuza)
    !if (size(reactive_zones)>0) then
    !    do i=1,size(init_sol_zones)
    !        if (.not. associated(init_sol_zones(i)%reactive_zone%chem_syst)) then
    !            call init_sol_zones(i)%reactive_zone%set_chem_syst_react_zone(this%chem_syst)
    !            !print *, associated(init_sol_zones(i)%solid_chem%reactive_zone%chem_syst)
    !        end if
    !        call init_sol_zones(i)%set_reactive_zone(reactive_zones(i))
    !    end do
    !    do i=1,size(init_gas_zones)
    !        if (.not. associated(init_gas_zones(i)%reactive_zone%chem_syst)) then
    !            call init_gas_zones(i)%reactive_zone%set_chem_syst_react_zone(this%chem_syst)
    !            !print *, associated(init_sol_zones(i)%solid_chem%reactive_zone%chem_syst)
    !        end if
    !        call init_gas_zones(i)%set_reactive_zone(reactive_zones(i))
    !    end do
    !end if

!> Target solids
    if (nsrz>0) then
        call this%read_tar_sol(path_pb//root,nsrz,ngrz)
    else
        !> When no solid reactive zones, allocate dummy solids/reactive zones
        !> needed by read_waters_init for waters with gas but no solid phase
        call this%allocate_target_solids_dummy()
        call this%allocate_reactive_zones_dummy()
        do i=1,this%num_target_solids_dummy
            call this%target_solids_dummy(i)%set_id(0)
            call this%target_solids_dummy(i)%set_mineral_zone(this%min_zone_dummy)
            call this%react_zones_dummy(i)%set_CV_params(this%CV_params)
            call this%react_zones_dummy(i)%set_chem_syst_react_zone(this%chem_syst)
            call this%react_zones_dummy(i)%set_speciation_alg_dimensions(.true.)
            call this%react_zones_dummy(i)%set_ind_eq_reacts()
            call this%react_zones_dummy(i)%allocate_ind_var_act_species()
            call this%react_zones_dummy(i)%set_stoich_mat_react_zone()
            call this%react_zones_dummy(i)%set_ind_gases_stoich_mat()
            call this%react_zones_dummy(i)%set_ind_mins_stoich_mat()
            call this%target_solids_dummy(i)%set_reactive_zone(this%react_zones_dummy(i))
        end do
    end if
    write(0,'(A,I0)') ' MARK B: target solids done. ngrz=',ngrz
    block; real(kind=8),allocatable::hpq(:); integer::ihpq
        do ihpq=1,200; allocate(hpq(1024*ihpq)); hpq=1d0; deallocate(hpq); end do
        write(0,'(A)') ' [HPC] post read_tar_sol'
    end block
!> Target gases
    if (ngrz>0) then
        call this%read_tar_gas(path_pb//root,ngrz)
    end if
    write(0,'(A)') ' MARK C: target gases done. calling read_waters_init'
!> Target waters
    call this%read_waters_init(path_pb//root,nsrz,ngrz)
    write(0,'(A)') ' MARK D: read_waters_init done (read_chemistry_CHEPROO returning)'
    block; real(kind=8),allocatable::hpq(:); integer::ihpq
        do ihpq=1,200; allocate(hpq(1024*ihpq)); hpq=1d0; deallocate(hpq); end do
        write(0,'(A)') ' [HPC] post read_waters_init'
    end block
!> The zone waters and zone solids (one per target water and reactive zone) are allocated and linked by
!> set_zone_waters, called from set_target_waters (RT_m): num_target_waters is only final there (with
!> edge-centred targets and Dirichlet BCs the boundary waters are added to the target waters), and the
!> solvers index zone_waters(idx,:) by the target water index.
end subroutine