!> subroutine read_aq_chem_init(this,filename,solid_chem_init,num_tar_sol,line,num_tar_aq)
!>    use aqueous_chemistry_m
!>    implicit none
!>    class(aqueous_chemistry_c) :: this
!>    character(len=*), intent(in) :: filename
!>    class(solid_chemistry_c), intent(in) :: solid_chem_init(:)
!>    integer(kind=4), intent(in) :: num_tar_sol(:) !> number of targets associated to each initial solid chemistry
!>    integer(kind=4), intent(inout) :: line
!>    integer(kind=4), intent(out) :: num_tar_aq !> number of targets associated to aqueous chemistry
!>    
!>    !type(aqueous_chemistry_c), allocatable :: initial_water_types(:)
!>    integer(kind=4) :: j,solid_chem_ind
!>    type(aq_species_c) :: species
!>    type(species_c) :: constrain
!>    character(len=256) :: str
!>    real(kind=8) :: aux,conc,temp
!>    logical :: flag
!>    !if (size(eq_reacts)>0 .and. size(kin_reacts)>0) then
!>                    !p_chem_syst=>my_chem_syst_eq_kin
!>                    !> chapuza
!>                    !my_chem_eq_kin%aqueous_chem=this%aqueous_chem
!>                    !my_chem_eq_kin%solid_chem=this%solid_chem
!>                    !my_chem_eq_kin%water_sources=this%water_sources
!>                    !p_chem=>my_chem_eq_kin
!>                !end if
!>                !p_chem_syst=>my_chem_syst_eq_kin
!>                !select type (chem_syst=>this%solid_chemistry%reactive_zone%chem_syst)
!>                !type is (chem_system_eq_c)
!>                !>    p_chem_syst=>my_chem_syst_eq
!>                !type is (chem_system_kin_c)
!>                !>    p_chem_syst=>my_chem_syst_kin
!>                !>    prim_species=aq_species
!>                !type is (chem_system_eq_kin_c)
!>                !>    p_chem_syst=>my_chem_syst_eq_kin
!>                !end select
!>                !call p_chem_syst%set_aq_species(aq_species)
!>                !call p_chem_syst%set_num_species()
!>                !call p_chem_syst%set_prim_species(prim_species)
!>                !call p_chem_syst%compute_num_prim_species()
!>                !call p_chem_syst%set_cst_act_species(cst_act_species)
!>                !call p_chem_syst%set_num_cst_act_species()
!>                !call p_chem_syst%set_solid_species(solids)
!>                !call p_chem_syst%set_num_solid_species()
!>                !call p_chem_syst%set_num_var_act_species()
!>                !call p_chem_syst%set_var_act_species()
!>                !call p_chem_syst%set_species()
!>                !call p_chem_syst%set_mineral_zone(mineral_zone)
!>                !call p_chem_syst%compute_num_eq_reactions()
!>                !call p_chem_syst%set_eq_reactions(eq_reacts)
!>                !call p_chem_syst%set_kin_reactions(kin_reacts)
!>                !!print *, chem_syst%num_eq_reactions, chem_syst%num_kin_reactions
!>                !call p_chem_syst%compute_num_kin_reactions()
!>                !call p_chem_syst%compute_num_reactions()
!>                !!print *, chem_syst%get_num_eq_reactions(), chem_syst%get_num_kin_reactions()
!>                !!print *, chem_syst%mineral_zone%num_minerals, chem_syst%mineral_zone%num_minerals_eq
!>                !call p_chem_syst%construct_stoich_mat_chem_system()
!>                !!call p_chem_syst%allocate_eq_csts()
!>                !call p_chem_syst%set_eq_csts()
!>                !call p_chem_syst%compute_num_prim_species()
!>                !call p_chem_syst%compute_num_cons_comps()
!>                !call p_chem_syst%compute_num_mob_cons_comps()
!>                !call p_chem_syst%compute_stoich_submats()
!>                !call p_chem_syst%compute_logK_star()
!>                !call p_chem_syst%compute_comp_mat()
!>                !call p_chem_syst%compute_comp_submats()
!>                !call p_chem_syst%compute_elim_mat()
!>                !call p_chem_syst%compute_mob_comp_mat()
!>                !call p_chem_syst%compute_mob_cons_comp_mat()
!>                !call this%set_chem_syst(p_chem_syst)
!>    !open(unit=5,file=filename,status='old',action='read',access='direct')
!>    !read(5,rec=line) str
!>    !if (str=='INITIAL AND BOUNDARY WATER TYPES') then
!>        !do
!>                !read(5,*) niwtype, nbwtype, nrwtype
!>                !allocate(initial_water_types(niwtype))
!>                !allocate(n_tar(niwtype))
!>                !allocate(gamma_1(size(prim_species)),gamma_2(p_chem_syst%get_num_eq_reactions()))
!>                !gamma_1=1d0
!>                !gamma_2=1d0
!>            
!>                !do i=1,niwtype
!>                    !if (i>niwtype) exit
!>                    read(5,*) solid_chem_ind, temp
!>                    if (solid_chem_ind<1 .or. solid_chem_ind>size(solid_chem_init)) error stop
!>                    !if (sum(n_tar(1:i))>size(this%aqueous_chem)) error stop
!>                    this%temp=temp
!>                    !call this%set_solid_chemistry(
!>                    call this%set_solid_chemistry(solid_chem_init(solid_chem_ind))
!>                    num_tar_aq=num_tar_sol(solid_chem_ind)
!>                    call this%allocate_concentrations()
!>                    select type (chem_syst=>this%solid_chemistry%reactive_zone%chem_syst)
!>                    class is (chem_system_eq_c)
!>                        j=1
!>                        do
!>                            line=line+1
!>                            read(5,*) species%name
!>                            if (species%name=='*') then
!>                                !!> Chapuza
!>                                !select type (p_chem_syst)
!>                                !class is (chem_system_eq_c)
!>                                    call this%compute_c2nc_from_c1_pure_phase()
!>                                    !i=i+1
!>                                    exit
!>                                !end select
!>                            end if
!>                            call belongs_to(species,chem_syst%prim_species,flag)
!>                            if (flag==.true.) then
!>                                backspace(5)
!>                                read(5,*) species%name, conc, constrain%name
!>                                !print *, str4, conc, constrain%name
!>                                call belongs_to(constrain,chem_syst%cst_act_species,flag)
!>                                if (flag==.true. .or. constrain%name=='') then
!>                                    this%concentrations(j)=conc
!>                                    j=j+1
!>                                    if (j>this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1) error stop
!>                                else
!>                                    error stop "This constrain is wrong"
!>                                end if
!>                            else if (species%name=='') then
!>                                continue
!>                            else
!>                                error stop "This is not a primary species"
!>                            end if
!>                        end do
!>                    end select
!>                    !call belongs_to(water,cst_act_species,flag,index)
!>                    !if (flag==.true.) then
!>                    !>    initial_water_types(i)%concentrations(p_chem_syst%num_species)=1d0
!>                    !else
!>                    !>    continue
!>                    !end if
!>                    !call this%get_chem_syst(p_chem_syst)
!>                    !call initial_water_types(i)%compute_conc_comp_mob(p_chem_syst)
!>                    !call compute_c2aq_from_c1_explicit(p_chem_syst,initial_water_types(i)%concentrations(1:p_chem_syst%num_prim_species),gamma_1,gamma_2,initial_water_types(i)%concentrations(chem_syst%num_prim_species+1:chem_syst%num_species))
!>               !> end do
!>    !else
!>    !>    error stop "Wrong record in read initial water types"
!>    !end if
!end subroutine