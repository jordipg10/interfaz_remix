!> Lectura sistema quimico
subroutine read_chem_system(this,filename)
    use chem_system_m, only: chem_system_c
    implicit none
    class(chem_system_c) :: this
    character(len=*), intent(in) :: filename
    
    real(kind=8), allocatable :: Sk(:,:),logK(:),gamma_1(:),gamma_2(:)
    integer(kind=4) :: i,j,num_sp,num_prim_sp,num_cst_act_sp,k,num_aq_sp,num_sec_aq_sp,wat_ind, nbwtype, nrwtype,index, & 
        kin_react_type,n_r,num_mins,num_species,num_surf_compl,num_exch_cats
    integer(kind=4), allocatable :: n_tar(:)
    real(kind=8) :: aux,conc,temp
    character(len=256) :: str,str1,str2,str3,str4,str5,path_name,Monod_name,file_kin_params
    logical :: flag,flag_wat
    
    character(len=256), allocatable :: aq_species_str(:),prim_species_str(:),cst_act_species_str(:),minerals_str(:), & 
        solid_species_str(:),kin_react_names(:)
    ! type(species_c) :: species
    ! type(species_c), allocatable :: surf_compl(:)
    ! type(aq_species_c) :: water
    ! type(aq_species_c), allocatable :: aq_species(:),exch_cats(:),prim_species(:)
    ! type(mineral_c) :: mineral
    ! type(mineral_c), allocatable :: mins(:)
    ! !type(cat_exch_c) :: cat_exch_obj
    ! type(surface_c) :: surf_ads_obj
    ! !type(reactive_zone_c) :: reactive_zone
    ! !type(solid_c) :: solid
    ! !class(solid_c), allocatable :: solids(:)
    ! type(gas_c) :: gas
    ! type(gas_c), allocatable :: gases(:)
    ! class(kin_params_c), pointer :: p_kin_params=>null()
    ! !type(Monod_params_c), allocatable, target :: Monod_params_array(:)
    ! !type(Monod_params_c), target :: Monod_params
    ! class(kin_reaction_c), pointer :: p_kin_react=>null()
    ! type(kin_reaction_poly_c) :: kin_react_ptr
    ! class(kin_reaction_poly_c), allocatable :: kin_reacts(:)
    ! !type(redox_kin_c), target :: Monod
    ! type(eq_reaction_c) :: eq_react
    ! type(eq_reaction_c), allocatable :: eq_reacts(:)
    ! !type(kin_lin_c), target :: linear
    !type(kin_mineral_c), target :: kin_mineral
    !type(glob_stoich_c), allocatable :: kin_stoich(:)
    !type(loc_stoich_c), allocatable :: eq_stoich(:,:)
    
    !class(chem_system_c), pointer :: p_chem_syst=>null()
    !type(chem_system_eq_kin_c), target :: my_chem_syst_eq_kin
    !type(chem_system_eq_c), target :: my_chem_syst_eq
    !type(chem_system_kin_c), target :: my_chem_syst_kin
    
    !type(chemistry_eq_kin_c), target :: my_chem_eq_kin
!>    
!>    open(unit=55,file=filename,status='old',action='read')
!>    num_sp=0
!>    num_cst_act_sp=0
!>    n_r=0
!>    do
!>        read(55,*) str
!>        if (str=='end') exit
!>        if (str=='PRIMARY AQUEOUS SPECIES') then
!>            num_prim_sp=0
!>            do
!>                read(55,*) str1
!>                !print *, len_trim(species%name)
!>                print *, trim(str1)
!>                !print *, species%name(1:1)
!>                if (str1=='*') then
!>                    backspace(55)
!>                    allocate(prim_species(num_prim_sp))
!>                    do i=1,num_prim_sp
!>                        backspace(55)
!>                        read(55,*) prim_species(num_prim_sp-i+1)%name, prim_species(num_prim_sp-i+1)%valence, prim_species(num_prim_sp-i+1)%params_act_coeff%ion_size_param, prim_species(num_prim_sp-i+1)%params_act_coeff%a_TJ, prim_species(num_prim_sp-i+1)%params_act_coeff%b_TJ
!>                        print *, prim_species(num_prim_sp-i+1)%name
!>                        prim_species(num_prim_sp-i+1)%cst_act_flag=.false.
!>                        backspace(55)
!>                    end do
!>                    !call this%aq_phase%set_aq_species(aq_species)
!>                    exit
!>                else
!>                    num_sp=num_sp+1
!>                    num_prim_sp=num_prim_sp+1
!>                end if
!>            end do
!>        else if (str=='AQUEOUS COMPLEXES') then
!>            num_aq_sp=num_prim_sp
!>            do
!>                read(55,*) species%name
!>                if (species%name=='*') then
!>                    backspace(55)
!>                    allocate(aq_species(num_aq_sp))
!>                    aq_species(1:num_prim_sp)=prim_species
!>                    num_sec_aq_sp=num_aq_sp-num_prim_sp
!>                    do i=1,num_sec_aq_sp
!>                        backspace(55)
!>                        read(55,*) aq_species(num_aq_sp-i+1)%name, aq_species(num_aq_sp-i+1)%valence, aq_species(num_aq_sp-i+1)%params_act_coeff%ion_size_param, aq_species(num_prim_sp-i+1)%params_act_coeff%a_TJ, aq_species(num_prim_sp-i+1)%params_act_coeff%b_TJ
!>                        if (aq_species(num_aq_sp-i+1)%name=='h2o') then
!>                            n_r=n_r+1
!>                            num_cst_act_sp=num_cst_act_sp+1
!>                            aq_species(num_aq_sp-i+1)%cst_act_flag=.true.
!>                            aq_species(num_aq_sp-i+1)%valence=0
!>                        else
!>                            aq_species(num_aq_sp-i+1)%cst_act_flag=.false.
!>                        end if
!>                        backspace(55)
!>                    end do
!>                    call this%aq_phase%set_aq_species(aq_species)
!>                    call this%aq_phase%rearrange_aq_species() !> sorts aqueous species into variable activity & constant activity
!>                    exit
!>                else
!>                    num_sp=num_sp+1
!>                    num_aq_sp=num_aq_sp+1
!>                end if
!>            end do               
!>        else if (str=='MINERALS') then
!>            num_mins=0
!>            do
!>                read(55,*) species%name
!>                if (species%name=='*') then
!>                    backspace(55)
!>                    allocate(mins(num_mins))
!>                    do i=1,num_mins
!>                        backspace(55)
!>                        read(55,*) mins(i)%name
!>                        backspace(55)
!>                    end do
!>                    call this%set_minerals(mins)
!>                    exit
!>                else
!>                    n_r=n_r+1
!>                    num_cst_act_sp=num_cst_act_sp+1
!>                    num_sp=num_sp+1
!>                    num_mins=num_mins+1
!>                end if
!>            end do
!>        else if (str=='GAseS') then
!>            num_species=0
!>            do
!>                read(55,*) species%name
!>                if (species%name=='*') then
!>                    backspace(55)
!>                    allocate(gases(num_species))
!>                    do i=1,num_species
!>                        backspace(55)
!>                        read(55,*) gases(i)%name
!>                    !> falta a�adir la cst_act_flag
!>                        !gases(i)%cst_act_flag==.false.
!>                        !if (aq_species(i)%name=='h2o') then
!>                        !>    aq_species(i)%cst_act_flag==.true.
!>                        !else
!>                        !>    aq_species(i)%cst_act_flag==.false.
!>                        !end if
!>                        backspace(55)
!>                    end do
!>                    !call this%set_gases(gases)
!>                    exit
!>                else
!>                    n_r=n_r+1
!>                    num_sp=num_sp+1
!>                    num_species=num_species+1
!>                end if
!>            end do
!>        else if (str=='SURFACE') then
!>            num_surf_compl=0
!>            do
!>                read(55,*) species%name
!>                if (species%name=='*') then
!>                    backspace(55)
!>                    call surf_compl_obj%allocate_surf_compl(num_surf_compl)
!>                    do i=1,surf_compl_obj%num_surf_compl
!>                        backspace(55)
!>                        read(55,*) surf_compl_obj%surf_compl(i)%name
!>                        backspace(55)
!>                        surf_compl_obj%surf_compl(i)%cst_act_flag=.true.
!>                    end do
!>                    call this%set_surf_compl_obj(surf_compl_obj)
!>                    exit
!>                else
!>                    n_r=n_r+1
!>                    num_cst_act_sp=num_cst_act_sp+1
!>                    num_sp=num_sp+1
!>                    num_surf_compl=num_surf_compl+1
!>                end if
!>            end do
!>        else if (str=='EXCHANGE') then
!>            num_exch_cats=0
!>            do
!>                read(55,*) species%name
!>                if (species%name=='*') then
!>                    backspace(55)
!>                    allocate(exch_cats(num_exch_cats))
!>                    do i=1,num_exch_cats
!>                        backspace(55)
!>                        read(55,*) exch_cats(i)%name
!>                        !mins(i)%cst_act_flag==.true.
!>                        !if (aq_species(i)%name=='h2o') then
!>                        !>    aq_species(i)%cst_act_flag==.true.
!>                        !else
!>                        !>    aq_species(i)%cst_act_flag==.false.
!>                        !end if
!>                        backspace(55)
!>                    end do
!>                    call cat_exch_obj%set_exch_cats(exch_cats)
!>                    exit
!>                else
!>                    n_r=n_r+1
!>                    num_sp=num_sp+1
!>                    num_exch_cats=num_exch_cats+1
!>                end if
!>            end do
!>            call this%set_cat_exch_obj(cat_exch_obj)
!>        else 
!>            continue
!>        end if
!>    end do
!>    close(55)
!>    call this%set_num_species(num_sp)
!>    call this%allocate_species()
!>    call this%set_num_reacts(n_r)
!>    call this%allocate_reacts()
!>    call this%set_num_cst_act_species(num_cst_act_sp)
!>    call this%allocate_cst_act_species()
!>    call this%set_species()
!>    !call this%read_master25()
!>    call water%set_name('h2o')
!>    call this%is_aq_species_in_chem_syst(water,flag_wat,wat_ind)
!>    if (flag_wat==.true.) then
!>        call eq_react%read_eq_reaction(water,'reacciones_disociacion.dat')
!>        call eq_react%append_eq_reaction(this%eq_reacts)
!>        this%num_eq_reacts=this%num_eq_reacts+1
!>    end if
!>    call this%read_kinetics_DB()
!>    !select type (kin_react=>this%kin_reacts(1)%kin_reaction)
!>    !type is (lin_kin_reaction_c)
!>    !>    print *, kin_react%lambda
!>    !end select
end subroutine