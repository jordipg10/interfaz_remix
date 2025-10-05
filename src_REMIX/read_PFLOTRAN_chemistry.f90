!> Reads PFLOTRAN CHEMISTRY block
subroutine read_chem_system_PFLOTRAN(this,path,unit)
    use chem_system_m, only: chem_system_c, species_c, aq_species_c, mineral_c, gas_c, surface_c, kin_params_c, kin_reaction_c,&
        kin_reaction_poly_c, eq_reaction_c
    use array_ops_m, only: append_int_1D_array
    implicit none
    class(chem_system_c) :: this
    character(len=*), intent(in) :: path
    integer(kind=4), intent(in) :: unit
    
    real(kind=8), allocatable :: Sk(:,:),logK(:),gamma_1(:),gamma_2(:),eq_csts(:)
    integer(kind=4) :: n_eq_homog,unit_DB,unit_kinetics,unit_redox,i,j,num_sp,ind_var_act_sp,num_var_act_sp,num_cst_act_sp,k,&
    num_aq_sp,num_sec_aq_sp,exch_cat_ind, n_min_kin, n_gas_kin,index,kin_react_type,n_r,num_mins,num_species,num_surf_compl,&
    num_exch_cats,n_eq,n_k,num_mins_eq_indices,n_mr,n_lin_kin
    integer(kind=4), allocatable :: n_tar(:),mins_eq_indices(:),gases_eq_indices(:),cols(:)
    real(kind=8) :: aux,conc,temp,SI,half_sat_cst,thr_conc
    real(kind=8), parameter :: tol=1d-12
    character(len=256) :: str,str1,str2,str3,str4,str5,Monod_name,file_kin_params,label,file_DB
    character(len=:), allocatable :: str_block_trim,str_trim,valence_str,exch_cat_val,exch_cat_name,label_trim
    logical :: flag,eq_label,exch_cat_flag,flag_comp,flag_surf,flag_Se
    
    character(len=256), allocatable :: aq_species_str(:),prim_species_str(:),cst_act_species_str(:),minerals_str(:),&
    solid_species_str(:),kin_react_names(:)
    type(species_c) :: species
    type(species_c), allocatable :: surf_compl(:)
    type(aq_species_c) :: exch_cat
    type(aq_species_c), allocatable :: aq_species(:),exch_cats(:),prim_species(:)
    type(mineral_c) :: mineral
    type(mineral_c), allocatable :: mins(:)
    type(surface_c) :: surf_ads_obj
    !type(surface_c) :: surf_compl_obj
    !type(reactive_zone_c) :: reactive_zone
    !type(solid_c) :: solid
    !class(solid_c), allocatable :: solids(:)
    type(gas_c) :: gas
    type(gas_c), allocatable :: gases(:)
    class(kin_params_c), pointer :: p_kin_params=>null()
    !type(Monod_params_c), allocatable, target :: Monod_params_array(:)
    !type(Monod_params_c), target :: Monod_params
    class(kin_reaction_c), pointer :: p_kin_react=>null()
    type(kin_reaction_poly_c) :: kin_react_ptr
    class(kin_reaction_poly_c), allocatable :: kin_reacts(:)
    !type(redox_kin_c), target :: Monod
    type(eq_reaction_c) :: eq_react
    type(eq_reaction_c), allocatable :: eq_reacts(:)
    
    !type(reaction_c), allocatable :: reacts(:)
    !type(kin_lin_c), target :: linear
    !type(kin_mineral_c), target :: kin_mineral
    !type(glob_stoich_c), allocatable :: kin_stoich(:)
    !type(loc_stoich_c), allocatable :: eq_stoich(:,:)
    
    !class(chem_system_c), pointer :: p_chem_syst=>null()
    !type(chem_system_eq_kin_c), target :: my_chem_syst_eq_kin
    !type(chem_system_eq_c), target :: my_chem_syst_eq
    !type(chem_system_kin_c), target :: my_chem_syst_kin
    
    !type(chemistry_eq_kin_c), target :: my_chem_eq_kin
    
    !open(unit,file=filename,status='old',action='read')
    
    num_sp=0
    num_aq_sp=0
    num_cst_act_sp=0
    num_var_act_sp=0
    num_surf_compl=0
    num_exch_cats=0
    num_species=0
    num_mins=0
    n_eq=0
    n_eq_homog=0
    n_k=0
    n_min_kin=0
    n_gas_kin=0
    n_lin_kin=0
    n_mr=0 !> n� microbial reactions
        
    do
        read(unit,"(A100)") label
        label_trim=trim(label)
        if (index(label_trim,'END')/=0) then
            rewind(unit)
            exit
        else if (index(label_trim,'PRIMARY_SPECIES')/=0) then
        !else if (index(label_trim=='AQUEOUS SPECIES') then
            !this%aq_phase%wat_flag=0
            do
                read(unit,*) str
                !if (str=='*') then
                !>    exit
                !else if (str=='h2o') then
                !>    error stop "Water cannot be a primary species"
                !else
                !>    continue
                !end if
                if (str=='/') then
                    exit
                else if (str=='H2O') then
                    num_cst_act_sp=num_cst_act_sp+1
                    this%aq_phase%wat_flag=1
                else
                    num_var_act_sp=num_var_act_sp+1
                end if
                num_aq_sp=num_aq_sp+1
                !num_var_act_sp=num_var_act_sp+1
                num_sp=num_sp+1
            end do
        else if (index(label_trim,'SECONDARY_SPECIES')/=0) then
            !num_aq_sp=num_prim_sp+num_sec_aq_sp
            do
                read(unit,*) str
                if (str=='/') then
                    exit
                else if (str=='H2O') then
                    num_cst_act_sp=num_cst_act_sp+1
                    this%aq_phase%wat_flag=1
                else
                    num_var_act_sp=num_var_act_sp+1
                end if
                num_aq_sp=num_aq_sp+1
                num_sp=num_sp+1
                n_eq=n_eq+1
                n_eq_homog=n_eq_homog+1
            end do            
        else if (index(label_trim,'MINERALS')/=0) then
            do
                read(unit,*) str
                if (str=='/') exit
                num_mins=num_mins+1
                num_cst_act_sp=num_cst_act_sp+1
                num_sp=num_sp+1
                !if (eq_label==.true.) then
                !>    n_eq=n_eq+1
                !else
                !>    n_k=n_k+1
                !>    n_min_kin=n_min_kin+1
                !end if
            end do 
        else if (index(label_trim,'GAS_SPECIES')/=0) then
            do
                read(unit,*) str, eq_label
                if (str=='/') exit
                num_species=num_species+1
                num_var_act_sp=num_var_act_sp+1
                num_sp=num_sp+1
                if (eq_label.eqv..true.) then
                    n_eq=n_eq+1
                else
                    n_k=n_k+1
                    n_gas_kin=n_gas_kin+1
                end if
            end do 
        else if (index(label_trim,'SORPTION')/=0) then
            do
                read(unit,*) str
                if (str=='/') exit
                !str_trim=trim(str)
                !num_surf_compl=num_surf_compl+1
                !num_var_act_sp=num_var_act_sp+1
                !num_sp=num_sp+1
                !exch_cat_ind=index(str_trim,'-')
                !if (exch_cat_ind>0 .and. exch_cat_ind<len(str_trim)-1) then
                !>    num_exch_cats=num_exch_cats+1
                !>    n_eq=n_eq+1
                !end if
            end do
        else if (index(label_trim,'MICROBIAL_REACTION')/=0) then
            do 
                read(unit,"(A100)") str
                str_trim=trim(str)
                if (index(str_trim,'REACTION')/=0) then
                    !str_react=str_trim(9:len(str_trim))
                    n_r=n_r+1
                    n_mr=n_mr+1
                    !read(unit,"(A100)") sub_str
                    !sub_str_trim=trim(sub_str)
                    !if (index(sub_str_trim,'RATE_CONSTANT')/=0) then
                        !read(sub_str_trim(14:len(sub_str_trim)),*) Monod%params%rate_cst
                else if (index(str_trim,'CONCENTRATION_UNITS')/=0) then
                    if (index(str_trim,'ACTIVITY')/=0) then
                    end if
                else if (index(str_trim,'MONOD')/=0 ) then
                    read(unit,*) label, species%name
                    read(unit,*) label, half_sat_cst
                    read(unit,*) label, thr_conc
                else
                    continue
                end if
            end do
        else if (index(label_trim,'DATABASE')/=0) then
            write(file_DB,*) label_trim(10:len(label))
        else 
            !error stop "This label is not implemented yet"
            continue
        end if
    end do
    
    call this%aq_phase%allocate_aq_species(num_aq_sp)
    call this%allocate_cst_act_sp_indices(num_cst_act_sp)
    call this%allocate_var_act_sp_indices(num_var_act_sp)
    call this%compute_num_species()
    !call this%allocate_species(num_sp)
    call this%allocate_minerals(num_mins)
    call this%gas_phase%allocate_gases(num_species)
    call this%cat_exch%allocate_surf_compl(num_surf_compl)
    call this%cat_exch%set_num_exch_cats(num_exch_cats)
    !this%num_eq_reacts_homog=n_eq_homog !> falta un set aqui
    call this%allocate_reacts(n_eq,n_k)
    !call this%set_num_minerals_kin(n_min_kin)
    call this%allocate_min_kin_reacts(n_min_kin)
    call this%allocate_redox_kin_reacts(n_mr)
    call this%allocate_lin_kin_reacts(n_lin_kin)
    call this%compute_num_reacts()
    
    !call this%solid_chemistry%reactive_zone%speciation_alg%set_dimensions(this%num_species,this%num_eq_reacts,this%num_cst_act_species,this%num_eq_reacts-this%cat_exch%num_surf_compl)
    
    ind_var_act_sp=0
    do
        read(unit,*) label
        if (label=='END') then
            exit
        else if (label=='PRIMARY_SPECIES') then
            i=0 !> counter primary aqueous species
            do
                read(unit,*) str
                if (str=='*') exit
                i=i+1
                str_trim=trim(str)
                call this%aq_phase%aq_species(i)%set_name(str_trim)
                if (str_trim=='h2o') then
                    call this%aq_phase%aq_species(this%aq_phase%num_species)%set_cst_act_flag(.true.)
                    !call this%cst_act_species(1)%set_name(str_trim)
                    !num_cst_act_sp=num_cst_act_sp+1
                else
                    ind_var_act_sp=ind_var_act_sp+1
                    call this%aq_phase%aq_species(i)%set_cst_act_flag(.false.)
                    !call this%var_act_species(ind_var_act_sp)%set_name(str_trim)
                    !num_var_act_sp=num_var_act_sp+1
                end if
                !call this%aq_phase%aq_species(i)%read_species(str_trim)
                !call this%aq_phase%aq_species(i)%set_cst_act_flag(.false.)
            end do
        !else if (index(label=='SECONDARY AQUEOUS SPECIES') then
        !>    do
        !>        read(unit,*) str
        !>        if (str=='*') exit
        !>        i=i+1
        !>        str_trim=trim(str)
        !>        call this%aq_phase%aq_species(i)%set_name(str_trim)
        !>        !if (this%aq_phase%aq_species(i)%name=='h2o') then
        !>        !>    call this%aq_phase%aq_species(i)%set_cst_act_flag(.true.)
        !>        !else
        !>        !>    call this%aq_phase%aq_species(i)%set_cst_act_flag(.false.)
        !>        !end if
        !>    end do
        else if (label=='SECONDARY_SPECIES') then
            do
                read(unit,*) str
                if (str=='*') exit
                i=i+1
                str_trim=trim(str)
                if (str_trim=='h2o') then
                    call this%aq_phase%aq_species(i)%set_cst_act_flag(.true.)
                    call this%aq_phase%aq_species(i)%set_name(str_trim)
                    !call this%cst_act_species(1)%set_name(str_trim)
                else
                    ind_var_act_sp=ind_var_act_sp+1
                    call this%aq_phase%aq_species(i)%set_cst_act_flag(.false.)
                    call this%aq_phase%aq_species(i)%set_name(str_trim)
                    !call this%var_act_species(ind_var_act_sp)%set_name(str_trim)
                end if
                !call this%aq_phase%aq_species(i)%set_cst_act_flag(.false.)
            end do
            !num_aq_sp=num_prim_sp
            !do
            !>    read(unit,*) species%name
            !>    if (species%name=='*') then
            !>        backspace(unit)
            !>        allocate(aq_species(num_aq_sp))
            !>        aq_species(1:num_prim_sp)=prim_species
            !>        num_sec_aq_sp=num_aq_sp-num_prim_sp
            !>        do i=1,num_sec_aq_sp
            !>            backspace(unit)
            !>            read(unit,*) aq_species(num_aq_sp-i+1)%name, aq_species(num_aq_sp-i+1)%valence, aq_species(num_aq_sp-i+1)%params_act_coeff%ion_size_param, aq_species(num_prim_sp-i+1)%params_act_coeff%a_TJ, aq_species(num_prim_sp-i+1)%params_act_coeff%b_TJ
            !>            if (aq_species(num_aq_sp-i+1)%name=='h2o') then
            !>                n_r=n_r+1
            !>                num_cst_act_sp=num_cst_act_sp+1
            !>                aq_species(num_aq_sp-i+1)%cst_act_flag=.true.
            !>                aq_species(num_aq_sp-i+1)%valence=0
            !>            else
            !>                aq_species(num_aq_sp-i+1)%cst_act_flag=.false.
            !>            end if
            !>            backspace(unit)
            !>        end do
            !>        call this%aq_phase%set_aq_species(aq_species)
            !>        call this%aq_phase%rearrange_aq_species() !> sorts aqueous species into variable activity & constant activity
            !>        exit
            !>    else
            !>        num_sp=num_sp+1
            !>        num_aq_sp=num_aq_sp+1
            !>    end if
            !end do               
        else if (label=='MINERALS') then
            allocate(mins_eq_indices(0))
            i=0 !> minerals
            do
                read(unit,*) str, eq_label
                if (str=='*') exit
                i=i+1
                !j=j+1
                str_trim=trim(str)
                call this%minerals(i)%set_phase_name(str_trim)
                !call this%cst_act_species(j)%set_name(str_trim)
                !call this%cst_act_species(i)%set_cst_act_flag(.true.)
                if (eq_label.eqv..true.) then
                    call append_int_1D_array(mins_eq_indices,i)
                end if
            end do
        else if (label=='GAS_SPECIES') then
            allocate(gases_eq_indices(0))
            i=0 !> gases
            do
                read(unit,*) str, eq_label
                if (str=='*') exit
                i=i+1
                str_trim=trim(str)
                call this%gas_phase%gases(i)%set_name(str_trim)
                !call this%var_act_species(i)%set_name(str_trim)
                !call this%var_act_species(i)%set_cst_act_flag(.false.)
                if (eq_label.eqv..true.) then
                    call append_int_1D_array(gases_eq_indices,i)
                end if
            end do
        else if (label=='SORPTION') then
            i=0 !> surface complexes
            j=0 !> exchangeable cations
            do
                read(unit,*) str
                if (str=='*') exit
                i=i+1
                str_trim=trim(str)
                call this%cat_exch%surf_compl(i)%set_name(str_trim)
                !call this%var_act_species(i)%set_name(str_trim)
                !call this%var_act_species(i)%set_cst_act_flag(.false.)
                exch_cat_ind=index(str_trim,'-')
                if (exch_cat_ind>0 .and. exch_cat_ind<len(str_trim)-1) then
                    j=j+1
                    if (exch_cat_ind>2) then
                        exch_cat_val=trim(str(2:exch_cat_ind-1))
                        allocate(character(len=(len(str_trim)-exch_cat_ind+1+len(exch_cat_val))) :: exch_cat_name)
                        write(exch_cat_name(1:len(str_trim)-exch_cat_ind),"(A)") str_trim(exch_cat_ind+1:len(str_trim))
                        write(exch_cat_name(len(str_trim)-exch_cat_ind+1:len(str_trim)-exch_cat_ind+2),"(A1)") '+'
                        write(exch_cat_name(len(exch_cat_name)-len(exch_cat_val)+1:len(exch_cat_name)),"(A)") exch_cat_val
                    end if
                    call exch_cat%set_name(exch_cat_name)
                    !call this%is_aq_species_in_chem_syst(exch_cat,exch_cat_flag,exch_cat_ind)
                    if (exch_cat_flag.eqv..true.) then
                        !call this%cat_exch%exch_cats(j)%assign_species(this%aq_phase%aq_species(exch_cat_ind))
                    else
                        error stop "Exchangeable cation is not in the chemical system"
                    end if
                    deallocate(exch_cat_name)
                end if
                
            end do
        !else if (str=='INITIAL AND BOUNDARY WATER TYPES') then
            !call this%
        !else if (str=='EXCHANGE') then
            !num_exch_cats=0
            !do
            !>    read(unit,*) species%name
            !>    if (species%name=='*') then
            !>        backspace(unit)
            !>        allocate(exch_cats(num_exch_cats))
            !>        do i=1,num_exch_cats
            !>            backspace(unit)
            !>            read(unit,*) exch_cats(i)%name
            !>            !mins(i)%cst_act_flag==.true.
            !>            !if (aq_species(i)%name=='h2o') then
            !>            !>    aq_species(i)%cst_act_flag==.true.
            !>            !else
            !>            !>    aq_species(i)%cst_act_flag==.false.
            !>            !end if
            !>            backspace(unit)
            !>        end do
            !>        call cat_exch_obj%set_exch_cats(exch_cats)
            !>        exit
            !>    else
            !>        n_r=n_r+1
            !>        num_sp=num_sp+1
            !>        num_exch_cats=num_exch_cats+1
            !>    end if
            !end do
            !call this%set_cat_exch_obj(cat_exch_obj)
        else if (label=='MICROBIAL_REACTION') then
            do 
                read(unit,*) str
                str_trim=trim(str)
                if (str_trim=='REACTION') then
                    !str_react=str_trim(9:len(str_trim))
                    n_r=n_r+1
                    n_mr=n_mr+1
                    !read(unit,"(A100)") sub_str
                    !sub_str_trim=trim(sub_str)
                    !if (index(sub_str_trim,'RATE_CONSTANT')/=0) then
                        !read(sub_str_trim(14:len(sub_str_trim)),*) Monod%params%rate_cst
                else if (str_trim=='CONCENTRATION_UNITS') then
                    !if (index(str_trim,'ACTIVITY')/=0) then
                    !end if
                else if (str_trim=='MONOD') then
                    read(unit,*) label, species%name
                    read(unit,*) label, half_sat_cst
                    read(unit,*) label, thr_conc
                else
                    continue
                end if
            end do
        else 
            continue
        end if
    end do
    !close(unit)
    !call this%set_num_reacts()
    !allocate(reacts(this%num_reacts))
    !unit_master_25=2
    !if (this%num_eq_reacts>0) then
    !call this%read_PFLOTRAN_DB(unit_DB,file_DB)
    !end if
    !i=1
    !j=1
    !num_mins_eq_indices=size(mins_eq_indices)
    !do
    !>    if (mins_eq_indices(i)==j) then
    !>        call this%eq_reacts(i)%copy_attributes(reacts(j))
    !>    else if (j<this%num_reacts) then
    !>        j=j+1
    !>    else if (i<num_mins_eq_indices) then
    !>        i=i+1
    !>        j=1
    !>    else
    !>        exit
    !>    end if
    !end do
    !do i=1,this%cat_exch%num_exch_cats
    !>    call this%eq_reacts(i+num_mins_eq_indices)%copy_attributes(reacts(this%num_minerals+i))
    !end do
    !unit_kinetics=3
    !if (this%num_minerals_kin>0) then
    !>    call this%read_kinetics_DB(unit_kinetics)
    !end if
    unit_redox=4
    if (this%num_redox_kin_reacts>0) then
        call this%read_Monod_reacts(path,unit_redox)
    end if
    !call this%set_species()
    call this%set_stoich_mat()
    
    eq_csts=this%get_eq_csts()
    call this%speciation_alg%compute_arrays(this%Se,eq_csts,tol,flag_Se,cols)
end subroutine