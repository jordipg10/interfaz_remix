!> Lectura sistema quimico basada en CHEPROO
!> Suponemos que el archivo ya ha sido abierto
subroutine read_chem_system_CHEPROO(this,path_DB,unit)
    use chem_system_m, only: chem_system_c, species_c, aq_species_c, mineral_c, surface_c, gas_c, kin_params_c, kin_reaction_c, &
    kin_reaction_poly_c, eq_reaction_c
    implicit none
    class(chem_system_c) :: this !> chemical system object
    character(len=*), intent(in) :: path_DB !> path to databases
    integer(kind=4), intent(in) :: unit !> file unit
    
    real(kind=8), allocatable :: Sk(:,:),logK(:),gamma_1(:),gamma_2(:)
    integer(kind=4) :: n_gas_eq_cst_act,n_gas_eq_var_act,n_gas_kin,exch_cat_valence,n_eq_homog,unit_master_25,unit_kinetics,&
        unit_redox,i,j,num_sp,ind_var_act_sp,num_var_act_sp,num_cst_act_sp,k,num_aq_sp,num_sec_aq_sp,exch_cat_ind, n_min_kin,&
        n_gas_eq,index,kin_react_type,n_r,num_mins,num_gases,num_surf_compl,num_exch_cats,n_eq,n_k,num_mins_eq_indices,n_redox,&
        n_lin_kin,n_lin_eq,num_aq_compl,num_mins_eq,num_cst_act_gases,num_var_act_gases,n_gas_kin_cst_act,n_gas_kin_var_act,l, &
        ind_diss_solids,num_eq_cst_act_mins,num_var_act_mins,n_redox_eq,n_redox_kin,num_kin_cst_act_mins,num_eq_var_act_mins,&
        num_kin_var_act_mins
    integer(kind=4), allocatable :: n_tar(:),mins_eq_indices(:),gases_eq_indices(:),indices_lin_reacts(:,:)
    real(kind=8) :: aux,conc,temp,SI,lambda,yield
    character(len=256) :: str,str1,str2,str3,str4,str5,Monod_name,file_kin_params,label
    character(len=:), allocatable :: str_block_trim,str_trim,valence_str,exch_cat_val,exch_cat_name,str_trim_1,str_trim_2
    logical :: flag,eq_label,exch_cat_flag,cst_act_label,flag_1,flag_2
    
    character(len=256), allocatable :: aq_species_str(:),prim_species_str(:),cst_act_species_str(:),minerals_str(:),&
        solid_species_str(:),kin_react_names(:)
    type(species_c) :: species
    type(species_c), allocatable :: surf_compl(:)
    type(aq_species_c) :: exch_cat,aq_sp_1,aq_sp_2
    type(aq_species_c), allocatable :: aq_species(:),exch_cats(:),prim_species(:)
    type(mineral_c) :: mineral
    type(mineral_c), allocatable :: mins(:)
    type(surface_c) :: cat_exch_obj
    type(gas_c) :: gas
    type(gas_c), allocatable :: gases(:)
    class(kin_params_c), pointer :: p_kin_params=>null()
    class(kin_reaction_c), pointer :: p_kin_react=>null()
    type(kin_reaction_poly_c) :: kin_react_ptr
    class(kin_reaction_poly_c), allocatable :: kin_reacts(:)
    type(eq_reaction_c) :: eq_react
    type(eq_reaction_c), allocatable :: eq_reacts(:)
    
    
    logical :: flag_comp,flag_surf
!> We initialise counters
    num_sp=0
    num_aq_sp=0
    num_aq_compl=0
    num_cst_act_sp=0
    num_cst_act_gases=0
    num_var_act_gases=0
    num_eq_cst_act_mins=0
    num_eq_var_act_mins=0
    num_kin_cst_act_mins=0
    num_kin_var_act_mins=0
    num_var_act_sp=0
    num_surf_compl=0
    num_exch_cats=0
    num_gases=0
    num_mins=0
    num_mins_eq=0
    n_eq=0
    n_k=0
    n_min_kin=0
    n_gas_eq=0
    n_gas_eq_cst_act=0
    n_gas_eq_var_act=0
    n_gas_kin_cst_act=0
    n_gas_kin_var_act=0
    n_gas_kin=0
    n_lin_kin=0
    n_lin_eq=0
    n_redox_eq=0
    n_redox_kin=0
!> First iteration of chemical system file
    do
        read(unit,*) label
        if (label=='end') then
            rewind(unit)
            exit
        else if (label=='PRIMARY AQUEOUS SPECIES') then !> suponemos ordenadas en primarias y secundarias
            do
                read(unit,*) str
                if (str=='*') then
                    exit
                else if (str=='h2o(p)') then
                !> we assume water has activity = 1
                    num_cst_act_sp=num_cst_act_sp+1
                    this%aq_phase%wat_flag=1
                    this%aq_phase%ind_wat=num_aq_sp+1
                else
                    num_var_act_sp=num_var_act_sp+1
                end if
                num_aq_sp=num_aq_sp+1
                num_sp=num_sp+1
            end do
        else if (label=='AQUEOUS COMPLEXES') then
            do
                read(unit,*) str
                if (str=='*') then
                    exit
                else if (str=='h2o') then
                !> we assume water has constant activity (approx 1)
                    num_cst_act_sp=num_cst_act_sp+1
                    this%aq_phase%wat_flag=1
                    this%aq_phase%ind_wat=num_aq_sp+1
                else
                    num_var_act_sp=num_var_act_sp+1
                end if
                num_aq_sp=num_aq_sp+1
                num_aq_compl=num_aq_compl+1
                num_sp=num_sp+1
                n_eq=n_eq+1
            end do            
        else if (label=='MINERALS') then
            do
                read(unit,*) str, eq_label, cst_act_label
                if (str=='*') exit
                num_mins=num_mins+1
                num_sp=num_sp+1
                if (cst_act_label.eqv..true. .AND. eq_label.eqv..true.) then
                    n_eq=n_eq+1
                    num_mins_eq=num_mins_eq+1
                    num_eq_cst_act_mins=num_eq_cst_act_mins+1
                    num_cst_act_sp=num_cst_act_sp+1
                else if (cst_act_label.eqv..false. .AND. eq_label.eqv..true.) then
                    n_eq=n_eq+1
                    num_mins_eq=num_mins_eq+1
                    num_eq_var_act_mins=num_eq_var_act_mins+1
                    num_var_act_sp=num_var_act_sp+1
                else if (cst_act_label.eqv..true. .AND. eq_label.eqv..false.) then
                    num_cst_act_sp=num_cst_act_sp+1
                    num_kin_cst_act_mins=num_kin_cst_act_mins+1
                    n_k=n_k+1
                    n_min_kin=n_min_kin+1
                else
                    num_kin_var_act_mins=num_kin_var_act_mins+1
                    num_var_act_sp=num_var_act_sp+1
                    n_k=n_k+1
                    n_min_kin=n_min_kin+1
                end if
            end do 
        else if (label=='GASES') then
            do
                read(unit,*) str, eq_label, cst_act_label
                if (str=='*') exit
                num_gases=num_gases+1
                num_sp=num_sp+1
                if (cst_act_label.eqv..true. .AND. eq_label.eqv..true.) then
                    n_eq=n_eq+1
                    n_gas_eq=n_gas_eq+1
                    n_gas_eq_cst_act=n_gas_eq_cst_act+1
                    num_cst_act_gases=num_cst_act_gases+1
                    num_cst_act_sp=num_cst_act_sp+1
                else if (cst_act_label.eqv..false. .AND. eq_label.eqv..true.) then
                    n_eq=n_eq+1
                    n_gas_eq=n_gas_eq+1
                    n_gas_eq_var_act=n_gas_eq_var_act+1
                    num_var_act_gases=num_var_act_gases+1
                    num_var_act_sp=num_var_act_sp+1
                else if (cst_act_label.eqv..true. .AND. eq_label.eqv..false.) then 
                    n_gas_kin=n_gas_kin+1
                    n_gas_kin_cst_act=n_gas_kin_cst_act+1
                    num_cst_act_gases=num_cst_act_gases+1
                    num_cst_act_sp=num_cst_act_sp+1
                else
                    n_gas_kin=n_gas_kin+1
                    n_gas_kin_var_act=n_gas_kin_var_act+1
                    num_var_act_gases=num_var_act_gases+1
                    num_var_act_sp=num_var_act_sp+1
                end if
            end do 
        else if (label=='SURFACE COMPLEXES') then
            do
                read(unit,*) str
                if (str=='*') exit
                str_trim=trim(str)
                num_surf_compl=num_surf_compl+1
                num_var_act_sp=num_var_act_sp+1
                num_sp=num_sp+1
                exch_cat_ind=index(str_trim,'-')
                if (exch_cat_ind>0 .and. exch_cat_ind<len(str_trim)) then
                    num_exch_cats=num_exch_cats+1
                    n_eq=n_eq+1
                end if
            end do
        else if (label=='REDOX REACTIONS') then
            do
                read(unit,*) str, eq_label
                if (str=='*') then
                    exit
                !else if (str=='aerobic degradation DOC (review)') then
                !    n_redox_eq=n_redox_eq+1
                !    n_eq=n_eq+1
                !else
                !    n_redox_kin=n_redox_kin+1
                !    n_k=n_k+1
                end if
                if (eq_label.eqv..true.) then
                    n_redox_eq=n_redox_eq+1
                    n_eq=n_eq+1
                else
                    n_redox_kin=n_redox_kin+1
                    n_k=n_k+1
                end if
            end do
        else if (label=='LINEAR REACTIONS') then
            do
                read(unit,*) str1, str2, lambda, eq_label
                if (str1=='*') exit
                if (eq_label.eqv..true.) then
                    n_lin_eq=n_lin_eq+1
                    n_eq=n_eq+1
                else
                    n_lin_kin=n_lin_kin+1
                    n_k=n_k+1
                end if
            end do
        else 
            continue
        end if
    end do
!> We set & allocate attributes
    !> Aqueous phase
    call this%aq_phase%allocate_aq_species(num_aq_sp)
    call this%aq_phase%allocate_ind_diss_solids()
    call this%aq_phase%set_num_aq_complexes(num_aq_compl)
    call this%aq_phase%set_num_var_act_species_phase(num_aq_sp-this%aq_phase%wat_flag)
    call this%aq_phase%set_num_cst_act_species_phase(this%aq_phase%wat_flag)
    !> Gas phase
    call this%gas_phase%allocate_gases(num_gases)
    call this%gas_phase%set_num_gases_eq(n_gas_eq) !> 
    call this%gas_phase%set_num_gases_eq_cst_act(n_gas_eq_cst_act) !> 
    call this%gas_phase%set_num_gases_eq_var_act(n_gas_eq_cst_act) !> 
    call this%gas_phase%set_num_gases_kin_cst_act(n_gas_kin_cst_act) !> 
    call this%gas_phase%set_num_gases_eq_var_act(n_gas_kin_var_act) !> 
    call this%gas_phase%set_num_gases_kin(n_gas_kin) !> 
    call this%gas_phase%set_num_var_act_species_phase(num_var_act_gases)
    call this%gas_phase%set_num_cst_act_species_phase(num_cst_act_gases)
    !> Cation exchange
    call this%cat_exch%allocate_surf_compl(num_surf_compl)
    call this%cat_exch%allocate_exch_cat_indices(num_exch_cats)
    call this%cat_exch%set_num_exch_cats(num_exch_cats)
    !> Chemical system
    call this%allocate_cst_act_sp_indices(num_cst_act_sp)
    call this%allocate_var_act_sp_indices(num_var_act_sp)
    call this%allocate_species()
    call this%allocate_minerals(num_mins)
    call this%set_num_minerals_eq(num_mins_eq)
    call this%set_num_minerals_eq_cst_act(num_eq_cst_act_mins)
    call this%set_num_minerals_eq_var_act(num_eq_var_act_mins)
    call this%allocate_reacts(n_eq,n_k)
    call this%set_num_redox_eq_reacts(n_redox_eq)
    call this%set_num_aq_eq_reacts()
    call this%allocate_min_kin_reacts(n_min_kin)
    call this%set_num_minerals_kin_cst_act(num_kin_cst_act_mins)
    call this%set_num_minerals_kin_var_act(num_kin_var_act_mins)
    call this%allocate_redox_kin_reacts(n_redox_kin)
    call this%allocate_lin_kin_reacts(n_lin_kin)
    call this%set_num_aq_kin_reacts()
    call this%compute_num_reacts()
    call this%compute_num_solids_chem_syst()
    !> Speciatrion algebra    
    call this%speciation_alg%set_flag_comp(.false.) !> chapuza
    if (num_surf_compl>0) then
        flag_surf=.true.
    else
        flag_surf=.false.
    end if
    call this%speciation_alg%set_flag_cat_exch(flag_surf)
    call this%speciation_alg%set_dimensions(this%num_species,this%num_eq_reacts,this%num_cst_act_species,this%aq_phase%num_species,&
    this%aq_phase%num_species-this%aq_phase%wat_flag,this%num_minerals_kin,num_gases-n_gas_eq)
!> Second iteration of chemical system file
    ind_var_act_sp=0 !< counter variable activity species
    ind_diss_solids=0 !> counter dissolved solids
    i=0 !> counter aqueous species
    do
        read(unit,*) label
        if (label=='end') then
            exit
        else if (label=='PRIMARY AQUEOUS SPECIES' .OR. label=='AQUEOUS COMPLEXES') then !> suponemos ordenadas en primarias y secundarias
            do
                read(unit,*) str
                if (str=='*') exit
                i=i+1
                str_trim=trim(str)
                call this%aq_phase%aq_species(i)%set_name(str_trim)
                if (str=='h2o(p)' .or. str=='h2o') then
                    call this%aq_phase%aq_species(i)%set_cst_act_flag(.true.)
                    this%cst_act_sp_indices(1)=i
                else
                    if (str=='h+') then
                        call this%aq_phase%set_ind_prot(i)
                    end if
                    ind_var_act_sp=ind_var_act_sp+1
                    ind_diss_solids=ind_diss_solids+1
                    call this%aq_phase%aq_species(i)%set_cst_act_flag(.false.)
                    this%var_act_sp_indices(ind_var_act_sp)=i
                    this%aq_phase%ind_diss_solids(ind_diss_solids)=i
                end if
                call this%species(i)%assign_species(this%aq_phase%aq_species(i))
            end do
        else if (label=='MINERALS') then
            i=0 !> counter minerals equilibrium var activity
            j=0 !> counter minerals kinetic variable activity
            k=0 !> counter minerals equilibrium constant activity
            l=0 !> counter minerals kinetic constant activity
            do
                read(unit,*) str, eq_label, cst_act_label
                if (str=='*') exit
                str_trim=trim(str)
                if ((eq_label .eqv. .true.) .and. (cst_act_label .eqv. .false.)) then !> equilibrium variable activity
                    i=i+1
                    call this%minerals(this%num_minerals_kin+i)%set_phase_name(str_trim)
                    call this%minerals(this%num_minerals_kin+i)%mineral%set_name(str_trim)
                    call this%minerals(this%num_minerals_kin+i)%mineral%set_cst_act_flag(cst_act_label)
                    call this%minerals(this%num_minerals_kin+i)%mineral%set_valence(0)
                else if (eq_label .eqv. .true. .and. cst_act_label .eqv. .true.) then !> equilibrium constant activity
                    k=k+1
                    call this%minerals(this%num_minerals-this%num_minerals_eq_cst_act+k)%set_phase_name(str_trim)
                    call this%minerals(this%num_minerals-this%num_minerals_eq_cst_act+k)%mineral%set_name(str_trim)
                    call this%minerals(this%num_minerals-this%num_minerals_eq_cst_act+k)%mineral%set_cst_act_flag(cst_act_label)
                    call this%minerals(this%num_minerals-this%num_minerals_eq_cst_act+k)%mineral%set_valence(0)
                else if (eq_label .eqv. .false. .and. cst_act_label .eqv. .false.) then !> kinetic variable activity
                    j=j+1
                    call this%minerals(j)%set_phase_name(str_trim)
                    call this%minerals(j)%mineral%set_name(str_trim)
                    call this%minerals(j)%mineral%set_cst_act_flag(cst_act_label)
                    call this%minerals(j)%mineral%set_valence(0)
                else !> kinetic constant activity
                    l=l+1
                    call this%minerals(this%num_minerals_kin_var_act+l)%set_phase_name(str_trim)
                    call this%minerals(this%num_minerals_kin_var_act+l)%mineral%set_name(str_trim)
                    call this%minerals(this%num_minerals_kin_var_act+l)%mineral%set_cst_act_flag(cst_act_label)
                    call this%minerals(this%num_minerals_kin_var_act+l)%mineral%set_valence(0)
                end if
            end do
        else if (label=='GASES') then
            i=0 !> counter gases equilibrium with constant activity
            j=0 !> counter gases equilibrium with variable activity
            k=0 !> counter gases kinetic with constant activity
            l=0 !> counter gases equilibrium with variable activity
            !> first gases in equilibrium with constant activity, then equilibrium with variable activity, & 
            !! then kinetic with constant activity, then kinetic with variable activity
            do
                read(unit,*) str, eq_label, cst_act_label
                if (str=='*') exit
                str_trim=trim(str)
                if ((eq_label .eqv. .true.) .and. (cst_act_label .eqv. .true.)) then
                    i=i+1
                    call this%gas_phase%gases(i)%set_name(str_trim)
                    call this%gas_phase%gases(i)%set_cst_act_flag(cst_act_label)
                    call this%gas_phase%gases(i)%set_valence(0)
                else if ((eq_label .eqv. .true.) .and. (cst_act_label .eqv. .false.)) then
                    j=j+1
                    call this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+j)%set_name(str_trim)
                    call this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+j)%set_cst_act_flag(cst_act_label)
                    call this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+j)%set_valence(0)
                else if ((eq_label .eqv. .false.) .and. (cst_act_label .eqv. .true.)) then
                    k=k+1
                    call this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+this%gas_phase%num_gases_eq_var_act+k)%set_name(&
                        str_trim)
                    call this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+&
                        this%gas_phase%num_gases_eq_var_act+k)%set_cst_act_flag(cst_act_label)
                    call this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+this%gas_phase%num_gases_eq_var_act+&
                        k)%set_valence(0)
                else
                    l=l+1
                    call this%gas_phase%gases(this%gas_phase%num_species-this%gas_phase%num_gases_kin_var_act+l)%set_name(str_trim)
                    call this%gas_phase%gases(this%gas_phase%num_species-this%gas_phase%num_gases_kin_var_act+l)%set_cst_act_flag(&
                        cst_act_label)
                    call this%gas_phase%gases(this%gas_phase%num_species-this%gas_phase%num_gases_kin_var_act+l)%set_valence(0)
                end if
            end do
        else if (label=='SURFACE COMPLEXES') then
            i=0 !> counter surface complexes
            j=0 !> counter exchangeable cations
            do
                read(unit,*) str
                if (str=='*') exit
                i=i+1
                str_trim=trim(str)
                call this%cat_exch%surf_compl(i)%set_name(str_trim)
                call this%cat_exch%surf_compl(i)%set_cst_act_flag(.false.)
                exch_cat_ind=index(str_trim,'-')
                if (exch_cat_ind>0 .and. exch_cat_ind<len(str_trim)) then
                    j=j+1
                    if (exch_cat_ind>1) then
                        if (exch_cat_ind==2) then
                            exch_cat_valence=1
                            allocate(character(len=(len(str_trim)-exch_cat_ind+1)) :: exch_cat_name)
                            write(exch_cat_name(1:len(str_trim)-exch_cat_ind),"(A)") str_trim(exch_cat_ind+1:len(str_trim))
                            write(exch_cat_name(len(str_trim)-exch_cat_ind+1:len(str_trim)-exch_cat_ind+1),"(A1)") '+'
                        else
                            exch_cat_val=trim(str(2:exch_cat_ind-1))
                            read(exch_cat_val,*) exch_cat_valence
                            allocate(character(len=(len(str_trim)-exch_cat_ind+1+len(exch_cat_val))) :: exch_cat_name)
                            write(exch_cat_name(1:len(str_trim)-exch_cat_ind),"(A)") str_trim(exch_cat_ind+1:len(str_trim))
                            write(exch_cat_name(len(str_trim)-exch_cat_ind+1:len(str_trim)-exch_cat_ind+1),"(A1)") '+'
                            write(exch_cat_name(len(exch_cat_name)-len(exch_cat_val)+1:len(exch_cat_name)),"(A)") exch_cat_val
                        end if
                    end if
                    call exch_cat%set_name(exch_cat_name)
                    call exch_cat%set_valence(exch_cat_valence)
                    call this%aq_phase%is_species_in_aq_phase(exch_cat,exch_cat_flag,exch_cat_ind)
                    if (exch_cat_flag.eqv..true.) then
                        this%cat_exch%exch_cat_indices(j)=exch_cat_ind
                    else
                        error stop "Exchangeable cation is not in the chemical system"
                    end if
                    deallocate(exch_cat_name)
                end if
                
            end do
        else if (label=='REDOX REACTIONS') then
            i=0 !> counter redox kinetic reactions
            j=0 !> counter redox equilibrium reactions
            do
                read(unit,*) str, eq_label!, yield !> reads name of reaction, equilibrium flag, and yield
                if (str=='*') then
                    exit
                end if
                str_trim=trim(str)
                if (eq_label.eqv..false.) then !> redox kinetic reactions
                    i=i+1
                    call this%redox_kin_reacts(i)%set_react_name(str_trim)
                    call this%redox_kin_reacts(i)%set_react_type(4)
                    call this%kin_reacts(this%num_lin_kin_reacts+i)%set_kin_reaction(this%redox_kin_reacts(i))
                    if (str=='aerobic degradation DOC (review)') then !> chapuza
                        call this%redox_kin_reacts(i)%set_yield(yield)
                    end if
                else !> redox eq reacts (they are the first in equilibrium reactions attribute of chemical system)
                    j=j+1
                    call this%eq_reacts(j)%set_react_name(str_trim)
                    call this%eq_reacts(j)%set_react_type(4)
                end if
            end do
        else if (label=='LINEAR REACTIONS') then
            i=0 !> counter linear kinetic reactions
            j=0 !> counter linear equilibrium reactions
            allocate(indices_lin_reacts(n_lin_kin+n_lin_eq,2)) !> first kinetic, then equilibrium
            do
                read(unit,*) str1, str2, lambda, eq_label
                if (str1=='*') exit
                str_trim_1=trim(str1)
                str_trim_2=trim(str2)
                call aq_sp_1%set_name(str_trim_1)
                call aq_sp_2%set_name(str_trim_2)
                if (eq_label.eqv..false.) then
                    i=i+1
                    call this%aq_phase%is_species_in_aq_phase(aq_sp_1,flag_1,indices_lin_reacts(i,1))
                    call this%aq_phase%is_species_in_aq_phase(aq_sp_2,flag_2,indices_lin_reacts(i,2))
                    if (flag_1.eqv..false. .or. flag_2.eqv..false.) then
                        error stop "Species in linear kinetic reaction not found in chemical system"
                    else
                        call this%lin_kin_reacts(i)%set_react_name(str_trim_1//'-->'//str_trim_2)
                        call this%lin_kin_reacts(i)%set_react_type(5)
                        call this%lin_kin_reacts(i)%allocate_reaction(2)
                        call this%lin_kin_reacts(i)%set_lambda(lambda)
                        call this%kin_reacts(i)%set_kin_reaction(this%lin_kin_reacts(i))
                    end if
                else
                    !j=j+1
                    !call this%eq_reacts(j)%set_react_name(str_trim_1//str_trim_2)
                    !call this%eq_reacts(j)%set_react_type(5)
                end if
            end do
        else 
            continue
        end if
    end do
!> We read databases
    !> Master25 (equilibrium reactions)
    unit_master_25=2
    call this%read_master25(path_DB,unit_master_25)
    !> Kinetics (minerals)
    unit_kinetics=3
    if (this%num_minerals_kin>0) then
        call this%read_kinetics_DB(path_DB,unit_kinetics)
    end if
    !> Monod reactions
    unit_redox=4
    if (this%num_redox_kin_reacts>0 .or. this%num_redox_eq_reacts>0) then
        call this%read_Monod_reacts(path_DB,unit_redox)
    end if
    !> Linear kinetcir eactions
    if (this%num_lin_kin_reacts>0) then
        do i=1,this%num_lin_kin_reacts
            call this%lin_kin_reacts(i)%set_all_species(this%aq_phase%aq_species(indices_lin_reacts(i,:)))
            call this%lin_kin_reacts(i)%set_stoichiometry([-1d0,1d0])
        end do
    end if
!> We rearrange species and equilibrium reactions
    call this%rearrange_species()
    call this%compute_z2()
    if (this%speciation_alg%flag_comp .eqv. .true.) then
        call this%rearrange_eq_reacts()
    end if
!> We set stoichiometric matrices
    call this%set_stoich_mat()
    call this%set_stoich_mat_gas()
    call this%set_stoich_mat_sol()
end subroutine