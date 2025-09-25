!> Reads PFLOTRAN database (eg 'hanford.dat')
subroutine read_PFLOTRAN_DB(this,unit,filename)
    use chem_system_m, only: chem_system_c
!> Arguments
    implicit none
    class(chem_system_c) :: this
    integer(kind=4), intent(in) :: unit
    character(len=*), intent(in) :: filename
!>    
!>    integer(kind=4) :: num_temp_pts,wat_ind,i,j,int_var,num_reactants,ind_pr_sp,ind_eq_reacts,ind_mins,ind_gases,ind_exch_cats,ind_cst_act_sp,ind_var_act_sp,num_temp_data,ind_aq_compl,num_pr_aq_sp,num_aq_sp,num_aq_compl,num_cst_act_sp,num_var_act_sp,num_eq_reacts,num_mins,num_species,num_exch_cats,counter_aq_Compl,counter_mins,counter_surf_compl,counter_gases,aq_sp_ind,min_ind,gas_ind,surf_compl_ind
!>    integer(kind=4), allocatable :: eq_react_indices(:),indices_exch_cats(:),aux_ind(:),indices_aq_sp(:)
!>    logical :: aq_sp_flag,min_flag,gas_flag,surf_compl_flag
!>    real(kind=8) :: unk,valence,mol_weight,diff_vol,mol_vol
!>    real(kind=8), allocatable :: global_stoich_coeffs(:),aux_coeffs(:),temp_pts(:),log_K(:)
!>    character(len=256) :: str,name,convention
!>    character(len=256), allocatable :: str1(:),str2(:),species_str(:),global_stoich_names(:)
!>        
!>    type(aq_species_c) :: aq_species
!>    type(mineral_c) :: mineral
!>    type(gas_c) :: gas
!>    type(species_c) :: surf_compl
!
!>    open(unit,file=filename,status='old',action='read')
!>    
!>    allocate(eq_react_indices(this%num_eq_reacts),indices_aq_sp(this%aq_phase%num_species),indices_exch_cats(this%cat_exch%num_exch_cats))
!!> First iteration of database
!>    num_aq_sp=0 !> number of aqueous species
!>    num_pr_aq_sp=0 !> number of primary aqueous species (CHEPROO)
!>    ind_aq_compl=0 !> index aqueous complexes
!>    num_aq_compl=0 !> number of aqueous complexes
!>    num_cst_act_sp=0 !> number of constant activity species
!>    num_var_act_sp=0 !> number of variable activity species
!>    num_eq_reacts=0 !> number of equilibrium reactions
!>    ind_mins=0 !> index minerals
!>    num_mins=0 !> number of minerals
!>    ind_gases=0 !> index gases
!>    num_species=0 !> number of gases
!>    ind_exch_cats=0 !> index exchangeable cations
!>    num_exch_cats=0 !> number of exchangeable cations
!>    wat_ind=0
!!> Temperature
!>    read(unit,*) str, num_temp_pts
!>    allocate(temp_pts(num_temp_pts),log_K(num_temp_pts))
!!> Primary species and colloid fields
!>    do
!>        read(unit,*,iostat=int_var) aq_species%name, aq_species%params_act_coeff%ion_size_param, valence, molecular_weight
!>        !if (int_var<0) then
!>        !>    exit
!>        if (aq_species%name=='null') then
!>            exit
!>        !else if (aq_species%name=='h2o') then
!>        !>    error stop "Water cannot be a primary species"
!>        else!> if (aq_species%name==this%aq_phase%aq_species(num_aq_sp+1)%name) then
!>            call this%is_aq_species_in_chem_syst(aq_species,aq_sp_flag,aq_sp_ind)
!>            if (aq_sp_flag==.true.) then 
!>                if (aq_species%name=='h2o') then
!>                    !error stop "Water is not a primary species"
!>                    !wat_ind=aq_sp_ind
!>                    call aq_species%set_cst_act_flag(.true.)
!>                    call this%cst_act_species(1)%assign_species(aq_species)
!>                    num_cst_act_sp=num_cst_act_sp+1
!>                else
!>                    call aq_species%set_cst_act_flag(.false.)
!>                    call this%var_act_species(aq_sp_ind)%assign_species(aq_species)
!>                    num_var_act_sp=num_var_act_sp+1
!>                end if
!>                indices_aq_sp(num_aq_sp+1)=aq_sp_ind
!>                this%aq_phase%aq_species(aq_sp_ind)=aq_species
!>                num_aq_sp=num_aq_sp+1
!>                num_pr_aq_sp=num_pr_aq_sp+1
!>            else
!>                continue
!>            end if
!>        !else
!>        !>    continue
!>        end if
!>    end do
!!> Secondary species
!>    do
!>        read(unit,*,iostat=int_var) aq_species%name, num_reactants
!>        if (aq_species%name=='null') then
!>            exit
!>        else!> if (aq_species%name==this%aq_phase%aq_species(num_aq_sp+1)%name) then
!>            call this%is_aq_species_in_chem_syst(aq_species,aq_sp_flag,aq_sp_ind)
!>            if (aq_sp_flag==.true.) then
!>                if (aq_species%name=='h2o') then
!>                    wat_ind=aq_sp_ind
!>                    call aq_species%set_cst_act_flag(.true.)
!>                    !call this%cst_act_species(1)%assign_species(aq_species)
!>                    num_cst_act_sp=num_cst_act_sp+1
!>                else
!>                    !if (wat_ind>0 .and. aq_sp_ind>wat_ind) then
!>                    !>    call this%var_act_species(aq_sp_ind-1)%assign_species(aq_species)
!>                    !else
!>                    !>    call this%var_act_species(aq_sp_ind)%assign_species(aq_species)
!>                    !end if    
!>                    call aq_species%set_cst_act_flag(.false.)
!>                    num_var_act_sp=num_var_act_sp+1
!>                end if
!>            !> We assume all primary species are aqueous
!>                call this%eq_reacts(aq_sp_ind-num_pr_aq_sp)%allocate_reaction(num_reactants+1)
!>                call this%eq_reacts(aq_sp_ind-num_pr_aq_sp)%set_react_type(1)
!>                call this%eq_reacts(aq_sp_ind-num_pr_aq_sp)%set_react_name(aq_species%name)
!>                eq_react_indices(aq_sp_ind-num_pr_aq_sp)=ind_aq_compl+1
!>                indices_aq_sp(num_aq_sp+1)=aq_sp_ind
!>                num_aq_compl=num_aq_compl+1
!>                num_aq_sp=num_aq_sp+1
!>            else
!>                continue
!>            end if
!>            !if (ind_pr_sp+ind_aq_compl<this%aq_phase%num_species) then
!>            !>    ind_aq_compl=ind_aq_compl+1
!>            !>    ind_eq_reacts=ind_eq_reacts+1
!>            !else
!>            !>    exit
!>            !end if
!>        !else
!>        !>    continue
!>        end if
!>        ind_aq_compl=ind_aq_compl+1
!>    end do
!>    num_eq_reacts=num_eq_reacts+num_aq_compl
!!> Gases
!>    do
!>        read(unit,*,iostat=int_var) gas%name, gas%mol_vol, diff_vol, num_reactants
!>        if (gas%name=='null') then
!>            exit
!>        else!> if (name==this%gas_phase%gases(num_species+1)%name) then
!>            call this%gas_phase%is_gas_in_gas_phase(gas,gas_flag,gas_ind)
!>            if (gas_flag==.true.) then
!>                call this%gas_phase%gases(gas_ind)%set_molecular_weight(mol_weight)
!>                !> falta set diff_vol
!>                num_species=num_species+1
!>                call this%var_act_species(num_var_act_sp+gas_ind)%set_name(gas%name)
!>                call this%var_act_species(num_var_act_sp+gas_ind)%set_cst_act_flag(.false.)
!>                call this%var_act_species(num_var_act_sp+gas_ind)%set_valence(0)
!>                !num_var_act_sp=num_var_act_sp+1
!>                call this%eq_reacts(num_eq_reacts+gas_ind)%allocate_reaction(num_reactants+1)
!>                call this%eq_reacts(num_eq_reacts+gas_ind)%set_react_type(6)
!>                eq_react_indices(num_eq_reacts+gas_ind)=ind_gases+1
!>                !num_eq_reacts=num_eq_reacts+1
!>                !if (num_species<this%gas_phase%num_species) then
!>                !>    num_species=num_species+1
!>                !>    
!>                !else
!>                !>    exit
!>                !end 
!>            else
!>                continue
!>            end if
!>        end if
!>        ind_gases=ind_gases+1
!>    end do
!>    num_var_act_sp=num_var_act_sp+num_species
!>    num_eq_reacts=num_eq_reacts+num_species
!>    !> Minerals
!>    do
!>        read(unit,*,iostat=int_var) mineral%name, mineral%mol_vol, num_reactants
!>        if (mineral%name=='null') then
!>            exit
!>        else!> if (mineral%name==this%minerals(num_mins+1)%name) then
!>            call this%is_mineral_in_chem_syst(mineral,min_flag,min_ind)
!>            if (min_flag==.true.) then
!>                this%minerals(min_ind)=mineral
!>                call this%cst_act_species(num_cst_act_sp+min_ind)%set_species_phase(mineral)
!>                !call this%cst_act_species(ind_cst_act_sp)%set_name(mineral%name)
!>                !call this%cst_act_species(ind_cst_act_sp)%set_cst_act_flag(.true.)
!>                !call this%cst_act_species(ind_cst_act_sp)%set_valence(0)
!>                call this%eq_reacts(num_eq_reacts+min_ind)%allocate_reaction(num_reactants+1)
!>                call this%eq_reacts(num_eq_reacts+min_ind)%set_react_type(2)
!>                eq_react_indices(num_eq_reacts+min_ind)=ind_mins+1
!>                num_mins=num_mins+1
!>                !num_cst_act_sp=num_cst_act_sp+1
!>                !num_eq_reacts=num_eq_reacts+1
!>            else
!>                continue
!>            end if
!>            !if (ind_mins<this%num_minerals) then
!>            !>    ind_mins=ind_mins+1
!>            !>    ind_cst_act_sp=ind_cst_act_sp+1
!>            !>    ind_eq_reacts=ind_eq_reacts+1
!>            !else
!>            !>    exit
!>            !end if
!>        !else
!>        !>    continue
!>        end if
!>        ind_mins=ind_mins+1
!>    end do
!>    num_cst_act_sp=num_cst_act_sp+num_mins
!>    num_eq_reacts=num_eq_reacts+num_mins
!>    !> Surface complexes
!>    if (this%cat_exch%num_surf_compl>0) then
!>        call this%var_act_species(num_var_act_sp+1)%set_name('x-')
!>        call this%var_act_species(num_var_act_sp+1)%set_cst_act_flag(.false.)
!>        call this%var_act_species(num_var_act_sp+1)%set_valence(-1)
!>    end if
!>    do
!>        read(unit,*,iostat=int_var) surf_compl%name, num_reactants
!>        if (surf_compl%name=='null') then
!>            exit
!>        else!> if (name==this%cat_exch%surf_compl(num_exch_cats+1)%name) then
!>            call this%cat_exch%is_surf_compl_in(surf_compl,surf_compl_flag,surf_compl_ind)
!>            if (surf_compl_flag==.true.) then
!>                num_exch_cats=num_exch_cats+1
!>                call this%var_act_species(num_var_act_sp+surf_compl_ind)%set_name(surf_compl%name)
!>                call this%var_act_species(num_var_act_sp+surf_compl_ind)%set_cst_act_flag(.false.)
!>                !call this%var_act_species(num_var_act_sp+surf_compl_ind)%set_valence(0)
!>                call this%eq_reacts(num_eq_reacts+surf_compl_ind-1)%allocate_reaction(num_reactants+1)
!>                call this%eq_reacts(num_eq_reacts+surf_compl_ind-1)%set_react_type(3)
!>                eq_react_indices(num_eq_reacts+surf_compl_ind-1)=ind_exch_cats+1
!>                indices_exch_cats(num_exch_cats)=surf_compl_ind
!>                !num_exch_cats=num_exch_cats+1
!>                !num_eq_reacts=num_eq_reacts+1
!>                !num_var_act_sp=num_var_act_sp+1
!>            else
!>                continue
!>            end if
!>        end if
!>        ind_exch_cats=ind_exch_cats+1
!>    end do
!>    num_var_act_sp=num_var_act_sp+this%cat_exch%num_surf_compl
!>    num_eq_reacts=num_eq_reacts+num_exch_cats
!>    !allocate(surf_compl_indices(num_exch_cats))
!>    !do
!>    !>    read(unit,*,iostat=int_var) name, num_reactants
!>    !>    if (name=='null') then
!>    !>        exit
!>    !>    else if (str==this%cat_exch%exch_cats(l)%name) then
!>    !>        call this%eq_reacts(ind_eq_reacts)%allocate_reaction(num_reactants+1)
!>    !>        ind_surf_compl=ind_surf_compl+1
!>    !>        if (ind_surf_compl>this%cat_exch%num_exch_cats) exit
!>    !>        ind_eq_reacts=ind_eq_reacts+1
!>    !>        n_eq=n_eq+1
!>    !>        backspace(15)
!>    !>        read(15,*) str, num_reactants
!>    !>        print *, str, num_reactants
!>    !>        backspace(15)
!>    !>        call react%allocate_reaction(num_reactants+1)
!>    !>        read(15,*) str, num_reactants, ((react%stoichiometry(j), react%species(j)%name), j=1,num_reactants), logK
!>    !>        call react%set_eq_cst(10**logK)
!>    !>        call this%cat_exch%exch_cats(l)%set_valence(int(valence))
!>    !>        if (l<this%cat_exch%num_exch_cats) then
!>    !>            l=l+1
!>    !>        else
!>    !>            continue
!>    !>        end if
!>    !>        reacts_old=reacts_new
!>    !>        deallocate(reacts_new)
!>    !>        allocate(reacts_new(n_eq))
!>    !>        do j=1,n_eq-1
!>    !>            reacts_new(j)=reacts_old(j)
!>    !>        end do
!>    !>        reacts_new(n_eq)=react
!>    !>        deallocate(reacts_old,react%species,react%stoichiometry)  
!>    !>    else
!>    !>        !error stop "This species is not present in the master25 database"
!>    !>        continue
!>    !>    end if
!>    !end do
!>    
!>    rewind(unit)
!>    
!!> Second iteration of database
!>    counter_aq_compl=1 !> counter aqueous complexes
!>    ind_aq_compl=0 !> index aqueous complexes 
!>    ind_cst_act_sp=0 !> index constant activity species
!>    ind_var_act_sp=0 !> index variable activity species
!>    !ind_eq_reacts=1 !> index equilibrium reactions
!>    ind_mins=0 !> index minerals
!>    counter_mins=1 !> counter minerals
!>    ind_gases=0 !> index gases
!>    counter_gases=1 !> counter gases
!>    ind_exch_cats=0 !> index exchangeable cations
!>    counter_surf_compl=1 !> counter surface complexes
!>    !> Temperature
!>    read(unit,*) str, num_temp_data
!>    !> Primary aqueous species
!>    do
!>        read(unit,*,iostat=int_var) name
!>        if (name=='null') exit
!>    end do
!>    ind_var_act_sp=num_pr_aq_sp
!>    !> Aqueous complexes
!>    if (num_aq_compl>0) then
!>        do
!>            !read(unit,*,iostat=int_var) aq_species%name, num_reactants, ((this%eq_reacts(ind_eq_reacts)%stoichiometry(j), this%eq_reacts(ind_eq_reacts)%species(j)%name), j=1,num_reactants), log_K, aq_species%ionic_radius, valence!, aq_species%params_act_coeff%b_TJ, aq_species%molecular_weight
!>            !read(unit,*,iostat=int_var) aq_species%name, num_reactants, (eq_react%stoichiometry(j), eq_react%species(j), j=1,num_reactants), log_K, aq_species%ionic_radius, valence!, aq_species%params_act_coeff%b_TJ, aq_species%molecular_weight
!>            !if (counter_aq_compl==eq_react_indices(ind_eq_reacts)) then
!>            call is_int_in_1D_array(counter_aq_compl,eq_react_indices(1:num_aq_compl),aq_sp_flag,ind_eq_reacts)
!>            if (aq_sp_flag==.true.) then
!>                ind_aq_compl=ind_aq_compl+1
!>                read(unit,*,iostat=int_var) aq_species%name, num_reactants, ((this%eq_reacts(ind_eq_reacts)%stoichiometry(j), this%eq_reacts(ind_eq_reacts)%species(j)%name), j=1,num_reactants), log_K, aq_species%params_act_coeff%ion_size_param, valence, aq_species%molecular_weight
!>                !if (aq_species%name=='null') exit
!>                if (aq_species%name=='h2o') then
!>                    ind_cst_act_sp=ind_cst_act_sp+1
!>                    call aq_species%set_cst_act_flag(.true.)
!>                    call this%cst_act_species(1)%assign_species(aq_species)
!>                else
!>                    ind_var_act_sp=ind_var_act_sp+1
!>                    call aq_species%set_cst_act_flag(.false.)
!>                    if (indices_aq_sp(num_pr_aq_sp+ind_aq_compl)>wat_ind .and. wat_ind>0) then
!>                        call this%var_act_species(indices_aq_sp(num_pr_aq_sp+ind_aq_compl)-1)%assign_species(aq_species)
!>                    else
!>                        call this%var_act_species(indices_aq_sp(num_pr_aq_sp+ind_aq_compl))%assign_species(aq_species)
!>                    end if
!>                end if
!>                this%aq_phase%aq_species(indices_aq_sp(num_pr_aq_sp+ind_aq_compl))=aq_species
!>                !do i=1,num_temp_pts
!>                !>    call this%eq_reacts(ind_eq_reacts)%set_eq_cst(10**(-log_K(i)))
!>                !end do
!>                call this%eq_reacts(ind_eq_reacts)%set_single_species(aq_species,this%eq_reacts(ind_eq_reacts)%num_species)
!>                this%eq_reacts(ind_eq_reacts)%stoichiometry(this%eq_reacts(ind_eq_reacts)%num_species)=-1d0
!>                call this%eq_reacts(ind_eq_reacts)%change_sign_stoichiometry()
!>                !if (ind_eq_reacts<this%num_eq_reacts) then
!>                !>    ind_eq_reacts=ind_eq_reacts+1
!>                !else
!>                !>    exit
!>                !end if
!>                !ind_aq_compl=ind_aq_compl+1
!>                !if (ind_pr_sp+ind_aq_compl<this%aq_phase%num_species) then
!>                !>    ind_aq_compl=ind_aq_compl+1
!>                !>    ind_eq_reacts=ind_eq_reacts+1
!>                !else
!>                !>    exit
!>                !end if
!>                !if (ind_aq_compl==num_aq_compl) exit
!>            else
!>                read(unit,*,iostat=int_var) name
!>                if (name=='null') exit
!>            end if
!>            counter_aq_compl=counter_aq_compl+1
!>        end do
!>        !ind_var_act_sp=ind_var_act_sp+num_aq_compl
!>    else
!>        do
!>            read(unit,*,iostat=int_var) name
!>            if (name=='null') exit
!>        end do
!>    end if
!!> Gases
!>    if (num_species>0) then
!>        do
!>            !read(unit,*,iostat=int_var) name, mol_weight, diff_vol, num_reactants, ((this%eq_reacts(ind_eq_reacts)%stoichiometry(j), this%eq_reacts(ind_eq_reacts)%species(j)%name), j=1,num_reactants), log_K
!>            !if (counter_gases==eq_react_indices(ind_eq_reacts)) then
!>            call is_int_in_1D_array(counter_gases,eq_react_indices(num_aq_compl+num_mins+1:num_aq_compl+num_mins+num_species),gas_flag,ind_eq_reacts)
!>            if (gas_flag==.true.) then
!>                ind_eq_reacts=ind_eq_reacts+num_aq_compl+num_mins
!>                read(unit,*,iostat=int_var) name, mol_vol, num_reactants, ((this%eq_reacts(ind_eq_reacts)%stoichiometry(j), this%eq_reacts(ind_eq_reacts)%species(j)%name), j=1,num_reactants), log_K, mol_weight
!>                !if (name=='null') then
!>                !>    exit
!>                !else
!>                    !call this%eq_reacts(ind_eq_reacts)%set_eq_cst(10**(-log_K))
!>                    call this%eq_reacts(ind_eq_reacts)%species(this%eq_reacts(ind_eq_reacts)%num_species)%set_name(name)
!>                    !call this%eq_reacts(ind_eq_reacts)%species(this%eq_reacts(ind_eq_reacts)%num_species)%set_mol_vol(mol_vol)
!>                    call this%eq_reacts(ind_eq_reacts)%species(this%eq_reacts(ind_eq_reacts)%num_species)%set_molecular_weight(mol_weight)
!>                    call this%eq_reacts(ind_eq_reacts)%species(this%eq_reacts(ind_eq_reacts)%num_species)%set_valence(0)
!>                    call this%eq_reacts(ind_eq_reacts)%species(this%eq_reacts(ind_eq_reacts)%num_species)%set_cst_act_flag(.false.) !> "true" solo si gas tiene presion constante
!>                    this%eq_reacts(ind_eq_reacts)%stoichiometry(this%eq_reacts(ind_eq_reacts)%num_species)=-1d0
!>                    call this%eq_reacts(ind_eq_reacts)%change_sign_stoichiometry()
!>                    ind_var_act_sp=ind_var_act_sp+1
!>                    ind_gases=ind_gases+1
!>                    !if (ind_eq_reacts<this%num_eq_reacts) then
!>                    !>    ind_eq_reacts=ind_eq_reacts+1
!>                    !else
!>                    !>    exit
!>                    !end if
!>                    !ind_gases=ind_gases+1
!>                    !if (ind_gases==num_species) exit
!>                !end if
!>            else
!>                read(unit,*,iostat=int_var) name
!>                if (name=='null') exit
!>            end if
!>            counter_gases=counter_gases+1
!>        end do
!>    else
!>        do
!>            read(unit,*,iostat=int_var) name
!>            if (name=='null') exit
!>        end do
!>    end if
!!> Minerals
!>    if (num_mins>0) then
!>        do
!>            !read(unit,*,iostat=int_var) name, mol_vol, num_reactants, ((this%eq_reacts(ind_eq_reacts)%stoichiometry(j), this%eq_reacts(ind_eq_reacts)%species(j)%name), j=1,num_reactants), log_K
!>            !if (counter_mins==eq_react_indices(ind_eq_reacts)) then
!>            call is_int_in_1D_array(counter_mins,eq_react_indices(num_aq_compl+1:num_aq_compl+num_mins),min_flag,ind_eq_reacts)
!>            if (min_flag==.true.) then
!>                ind_eq_reacts=ind_eq_reacts+num_aq_compl
!>                read(unit,*,iostat=int_var) name, mol_vol, num_reactants, ((this%eq_reacts(ind_eq_reacts)%stoichiometry(j), this%eq_reacts(ind_eq_reacts)%species(j)%name), j=1,num_reactants), log_K
!>                !if (name=='null') then
!>                !>    exit
!>                !else
!>                    !call this%eq_reacts(ind_eq_reacts)%set_eq_cst(10**(-log_K))
!>                    call this%eq_reacts(ind_eq_reacts)%species(this%eq_reacts(ind_eq_reacts)%num_species)%set_name(name)
!>                    call this%eq_reacts(ind_eq_reacts)%species(this%eq_reacts(ind_eq_reacts)%num_species)%set_valence(0)
!>                    call this%eq_reacts(ind_eq_reacts)%species(this%eq_reacts(ind_eq_reacts)%num_species)%set_cst_act_flag(.true.)
!>                    this%eq_reacts(ind_eq_reacts)%stoichiometry(this%eq_reacts(ind_eq_reacts)%num_species)=-1d0
!>                    call this%eq_reacts(ind_eq_reacts)%change_sign_stoichiometry()
!>                    ind_mins=ind_mins+1
!>                    ind_cst_act_sp=ind_cst_act_sp+1
!>                    !if (ind_eq_reacts<this%num_eq_reacts) then
!>                    !>    ind_eq_reacts=ind_eq_reacts+1
!>                    !else
!>                    !>    exit
!>                    !end if
!>                    !ind_mins=ind_mins+1
!>                    !if (ind_mins==num_mins) exit
!>                !end if
!>            else
!>                read(unit,*,iostat=int_var) name
!>                if (name=='null') exit
!>            end if
!>            counter_mins=counter_mins+1
!>        end do
!>    else
!>        do
!>            read(unit,*,iostat=int_var) name
!>            if (name=='null') exit
!>        end do
!>    end if
!!> Surface complexes
!>    if (num_exch_cats>0) then
!>        do
!>            !read(unit,*,iostat=int_var) name, num_reactants, ((this%eq_reacts(ind_eq_reacts)%stoichiometry(j), this%eq_reacts(ind_eq_reacts)%species(j)%name), j=1,num_reactants), log_K, valence
!>            !if (counter_surf_compl==eq_react_indices(ind_eq_reacts)) then
!>            call is_int_in_1D_array(counter_surf_compl,eq_react_indices(num_aq_compl+num_mins+num_species+1:num_aq_compl+num_mins+num_species+num_exch_cats),surf_compl_flag,ind_eq_reacts)
!>            if (surf_compl_flag==.true.) then
!>                ind_exch_cats=ind_exch_cats+1
!>                ind_eq_reacts=ind_eq_reacts+num_aq_compl+num_mins+num_species
!>                read(unit,*,iostat=int_var) surf_compl%name, num_reactants, ((this%eq_reacts(ind_eq_reacts)%stoichiometry(j), this%eq_reacts(ind_eq_reacts)%species(j)%name), j=1,num_reactants), log_K, surf_compl%valence
!>                this%cat_exch%surf_compl(indices_exch_cats(ind_exch_cats))=surf_compl
!>                !call this%cat_exch%surf_compl(ind_surf_compl)%set_valence(nint(valence))
!>                !call this%cat_exch%surf_compl(ind_surf_compl)%set_cst_act_flag(.false.)
!>                call this%var_act_species(ind_var_act_sp+indices_exch_cats(ind_exch_cats))%assign_species(surf_compl)
!>                !call this%eq_reacts(ind_eq_reacts)%set_eq_cst(10**(-log_K))
!>                call this%eq_reacts(ind_eq_reacts)%set_single_species(surf_compl,this%eq_reacts(ind_eq_reacts)%num_species)
!>                this%eq_reacts(ind_eq_reacts)%stoichiometry(this%eq_reacts(ind_eq_reacts)%num_species)=-1d0
!>                call this%eq_reacts(ind_eq_reacts)%change_sign_stoichiometry()
!>            else
!>                read(unit,*,iostat=int_var) name
!>                if (name=='null') exit
!>            end if
!>            counter_surf_compl=counter_surf_compl+1
!>        end do
!>        ind_var_act_sp=ind_var_act_sp+this%cat_exch%num_surf_compl
!>    else
!>        do
!>            read(unit,*,iostat=int_var) name
!>            if (name=='null') exit
!>        end do
!>    end if
!>    !do
!>    !>    read(unit,*,iostat=int_var) name, num_reactants
!>    !>    if (name=='null') then
!>    !>        exit
!>    !>    else if (str==this%cat_exch%exch_cats(l)%name) then
!>    !>        call this%eq_reacts(ind_eq_reacts)%allocate_reaction(num_reactants+1)
!>    !>        ind_surf_compl=ind_surf_compl+1
!>    !>        if (ind_surf_compl>this%cat_exch%num_exch_cats) exit
!>    !>        ind_eq_reacts=ind_eq_reacts+1
!>    !>        n_eq=n_eq+1
!>    !>        backspace(15)
!>    !>        read(15,*) str, num_reactants
!>    !>        print *, str, num_reactants
!>    !>        backspace(15)
!>    !>        call react%allocate_reaction(num_reactants+1)
!>    !>        read(15,*) str, num_reactants, ((react%stoichiometry(j), react%species(j)%name), j=1,num_reactants), logK
!>    !>        call react%set_eq_cst(10**logK)
!>    !>        call this%cat_exch%exch_cats(l)%set_valence(int(valence))
!>    !>        if (l<this%cat_exch%num_exch_cats) then
!>    !>            l=l+1
!>    !>        else
!>    !>            continue
!>    !>        end if
!>    !>        reacts_old=reacts_new
!>    !>        deallocate(reacts_new)
!>    !>        allocate(reacts_new(n_eq))
!>    !>        do j=1,n_eq-1
!>    !>            reacts_new(j)=reacts_old(j)
!>    !>        end do
!>    !>        reacts_new(n_eq)=react
!>    !>        deallocate(reacts_old,react%species,react%stoichiometry)  
!>    !>    else
!>    !>        !error stop "This species is not present in the master25 database"
!>    !>        continue
!>    !>    end if
!>    !end do
!>    close(unit)
!>    call this%aq_phase%rearrange_aq_species()
end subroutine