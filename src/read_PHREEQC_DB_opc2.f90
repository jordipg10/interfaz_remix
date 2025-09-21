!> Reads PHREEQC database assuming we already know the species and phases of our chemical system
subroutine read_PHREEQC_DB_opc2(this,unit,filename)
    use chem_system_m, only: chem_system_c, species_c, aq_species_c, gas_phase_c, gas_c, mineral_c, phase_c, kin_params_c, &
        kin_reaction_c, kin_reaction_poly_c, eq_reaction_c, params_spec_vol_Redlich_c
    use array_ops_m, only: append_int_1D_array
    implicit none
    class(chem_system_c) :: this
    integer(kind=4), intent(in) :: unit
    character(len=*), intent(in) :: filename
    
    real(kind=8), allocatable :: Sk(:,:),logK(:),gamma_1(:),gamma_2(:),coeffs_Redlich(:)
    integer(kind=4) :: aq_counter,min_counter,gas_counter,eq_ind,valence,line_count,i,j,k,l,num_surf_compl,num_exch_sp,num_mins,&
    num_mins_eq,num_nf_sp,sp_ind,end_of_file, nbwtype, nrwtype,min_ind,eq_react_ind,num_gas,num_eq_ph,chem_syst_type,n_eq,n_k,&
    num_react_zones,aq_sp_ind,solid_chem_ind,niter,eq_phase_ind
    integer(kind=4), allocatable :: eq_react_indices(:),min_eq_ind(:),n_tar_sol(:),min_indices(:),nf_sp_indices(:),aq_sp_indices(:)
    real(kind=8) :: gfw_elem,aux,conc_init,temp,value,target_SI,logK_25,delta_h_25,alk_contrib
    character(len=256) :: gfw_form,phase_str,str,str_block,str3,str4,str5,path_name,Monod_name,file_kin_params,chem_syst_file,&
    sol_str,eq_phase_str
    logical :: flag_sp,flag_min,flag_gas,flag_react,prim_flag
    character(len=:), allocatable :: str_block_trim,str_trim,valence_str
    
    character(len=:), allocatable :: aq_species_str(:),prim_species_str(:),cst_act_species_str(:),minerals_str(:),&
    solid_species_str(:),kin_react_names(:)
    
    type(species_c) :: species,water,constrain
    type(species_c), allocatable :: aq_species_array(:),cst_act_species(:),prim_species(:)
    type(aq_species_c) :: aq_species,def_species,elem
    type(gas_phase_c), target :: gas_phase
    type(gas_c) :: gas
    type(gas_c), allocatable :: gases(:)
    type(mineral_c), target :: mineral
    class(phase_c), pointer :: phase=>null()
    class(kin_params_c), pointer :: p_kin_params=>null()
    class(kin_reaction_c), pointer :: p_kin_react=>null()
    type(kin_reaction_poly_c) :: kin_react_ptr
    class(kin_reaction_poly_c), allocatable :: kin_reacts(:)
    type(eq_reaction_c) :: eq_react
    type(eq_reaction_c), allocatable :: eq_reacts(:)
    
    type(chem_system_c) :: my_chem_syst
    
    type(params_spec_vol_Redlich_c) :: params_spec_vol_Redlich
    
    open(unit=5,file=filename,status='old',action='read')
    line_count=1
    !allocate(str_block(30))
    allocate(this%eq_reacts(0),this%aq_phase%aq_species(0))
    allocate(coeffs_Redlich(10))
    call eq_react%allocate_reaction(2)
    aq_counter=1 !> aqueous species counter
    min_counter=0 !> mineral counter
    gas_counter=0 !> gas counter
    do
        read(5,"(A100)",iostat=end_of_file) str_block
        if (end_of_file<0) exit
        !if (str_block(1:4)=='-end') then
            !if (end_of_file==0) then
            !>    end_of_file=end_of_file+1
            !>    rewind(5)
            !else
            !>    exit
            !end if
        if (str_block(1:len_trim(str_block))=='SOLUTION_MASTER_SPECIES') then
            !read(5,"(A30)") str
            !call this%read_solution_species_PHREEQC(str)
            
            do
                read(5,"(A256)") str
                str_trim=trim(str)
                if (index(str_trim,'#')==1 .or. index(str_trim,'#')==2) then
                    cycle
                !else if (str_trim(1:1)/='') then
                    !call this%read_sol_master_species_PHREEQC(str_trim,aq_species,prim_flag)
                    !call aq_species%append_aq_species(this%aq_phase%aq_species)
                !else if (index(str_trim,'log_k')/=0) then !> logK a 25 �C
                !>    backspace(5)
                !>    read(5,*) str, logK_25
                !>    call eq_react%set_eq_cst(logK_25)
                !else if (index(str,'delta_h')/=0) then
                !>    backspace(5)
                !>    read(5,*) str, delta_h_25
                !>    call eq_react%set_delta_h(delta_h_25)
                !else if (index(str,'-analytic')/=0) then
                !>    call eq_react%allocate_coeffs_logK_T()
                !>    backspace(5)
                !>    read(5,*) str, (eq_react%coeffs_logK_T(j), j=1,6)
                !else if (index(str,'-gamma')/=0) then
                !>    backspace(5)
                !>    read(5,*) str, aq_species%params_act_coeff%a_TJ, aq_species%params_act_coeff%b_TJ
                else if (index(str_trim,this%aq_phase%aq_species(aq_counter)%name)/=0) then 
                    backspace(5)
                    read(5,*) elem%name, aq_species%name, alk_contrib, gfw_form, gfw_elem
                    if (aq_counter<this%aq_phase%num_species) then
                        aq_counter=aq_counter+1
                    end if           
                else if (index(str_trim,'Alkalinity')/=0) then
                    !> por hacer
                else if (str=='') then
                    exit                    
                end if
            end do
        else if (str_block(1:len_trim(str_block))=='SOLUTION_SPECIES') then
            !read(5,"(A30)") str
            !call this%read_solution_species_PHREEQC(str)
            do
                read(5,"(A256)") str
                str_trim=trim(str)
                !print *, index(str_trim,'H')
                if (index(str_trim,'#')==1 .or. index(str_trim,'#')==2) then
                    !print *, str_trim
                    cycle
                else if (index(str_trim,this%aq_phase%aq_species(aq_counter)%name)/=0) then
                    if (aq_counter<this%aq_phase%num_species) then
                        aq_counter=aq_counter+1
                    end if
                    if (eq_react%num_species>1) then
                        !call this%append_reaction(eq_react)
                        call eq_react%deallocate_reaction()
                    end if
                    !call this%aq_phase%append_aq_species(aq_species)
                    call eq_react%read_association_react_PHREEQC(str_trim,prim_flag)
                else if (index(str_trim,'log_k')==3) then !> logK a 25 �C
                    backspace(5)
                    read(5,*) str, logK_25
                    call eq_react%set_eq_cst(logK_25)
                else if (index(str_trim,'delta_h')==3) then
                    backspace(5)
                    read(5,*) str, delta_h_25
                    call eq_react%set_delta_h(delta_h_25)
                else if (index(str_trim,'analytic')==3) then
                    !call eq_react%allocate_coeffs_logK_T()
                    !backspace(5)
                    !read(5,*) str, (eq_react%coeffs_logK_T(j), j=1,6)
                    !read(5,*) str, eq_react%coeffs_logK_T
                else if (index(str_trim,'gamma')==3) then
                    backspace(5)
                    read(5,*) str, this%aq_phase%aq_species(aq_counter-1)%params_act_coeff%a_TJ, &
                    this%aq_phase%aq_species(aq_counter-1)%params_act_coeff%b_TJ
                else if (index(str_trim,'Vm')==3) then
                    !backspace(5)
                    !read(5,*) str, (coeffs_Redlich(j), j=1,10)
                    !read(5,*) str, coeffs_Redlich
                    call params_spec_vol_Redlich%set_params(coeffs_Redlich)
                    call this%aq_phase%aq_species(aq_counter-1)%set_params_spec_vol_Redlich(params_spec_vol_Redlich)
                else if (str=='') then
                    exit                    
                end if
            end do
        else if (str_block(1:len_trim(str_block))=='PHAseS') then
            do
                read(5,"(A50)") phase_str
                !print *, phase_str
                if (index(str_trim,'#')==1 .or. index(str_trim,'#')==2) then
                    !print *, str_trim
                    cycle
                else if (index(phase_str,'(g)')==0) then
                    phase=>mineral
                    if (min_counter<this%num_minerals) then
                        min_counter=min_counter+1
                    end if
                else
                    phase=>gas_phase
                    if (gas_counter<this%gas_phase%num_species) then
                        gas_counter=gas_counter+1
                    end if
                    exit
                end if
                do !while (str_trim(1:1)/='')
                    read(5,"(A256)") str
                    str_trim=trim(str)
                    !print *, str_trim
                    if (index(str_trim,'#')==1 .or. index(str_trim,'#')==2) then
                        !print *, str_trim
                        cycle
                    else if (index(str_trim,this%minerals(min_counter)%name)/=0 .or. &
                        index(str_trim,this%gas_phase%gases(gas_counter)%name)/=0) then
                        if (eq_react%num_species>1) then
                            !call this%append_reaction(eq_react)
                            call eq_react%deallocate_reaction()
                        end if
                        !call this%append_phase(phase)
                        call eq_react%read_dissolution_react_PHREEQC(str_trim)
                        !call aq_species%append_aq_species(this%aq_phase%aq_species)
                    else if (index(str_trim,'log_k')/=0) then !> logK a 25 �C
                        backspace(5)
                        read(5,*) str, logK_25
                        call eq_react%set_eq_cst(logK_25)
                    else if (index(str,'delta_h')/=0) then
                        backspace(5)
                        read(5,*) str, delta_h_25
                        call eq_react%set_delta_h(delta_h_25)
                    else if (index(str,'analytic')/=0) then
                        !> call eq_react%allocate_coeffs_logK_T()
                        !backspace(5)
                        !read(5,*) str, (eq_react%coeffs_logK_T(j), j=1,6)
                    else if (index(str_trim,'Vm')==3) then
                        backspace(5)
                        read(5,*) str, this%minerals(min_counter)%mineral%mol_vol
                    else if (index(str_trim,'T_c')==3) then
                        backspace(5)
                        read(5,*) str, this%gas_phase%gases(gas_counter)%crit_temp
                    else if (index(str_trim,'P_c')==3) then
                        backspace(5)
                        read(5,*) str, this%gas_phase%gases(gas_counter)%crit_press
                    else if (index(str_trim,'Omega')==3) then
                        backspace(5)
                        read(5,*) str, this%gas_phase%gases(gas_counter)%acentric_fact
                    else if (str_trim(1:1)/='') then
                        backspace(5)
                        exit
                    else
                        continue
                    end if
                    !call this%append_phase(phase)
                    !call this%append_reaction(eq_react)
                end do
            end do
    !>    else if (str=='EQUILIBRIUM_PHAseS') then
    !>        num_mins=0
    !>        num_mins_eq=0
    !>        num_gas=0
    !>        num_eq_ph=0
    !>        allocate(min_indices(0),nf_sp_indices(0))
    !>        do
    !>            read(5,*) species%name, target_SI, conc_init
    !>            if (species%name=='*') then
    !>                call reactive_zone%set_chem_syst_react_zone(this%solid_chemistry%reactive_zone%chem_syst)
    !>                if (num_mins>0 .and. num_gas==0) then
    !>                    backspace(5)
    !>                    call reactive_zone%allocate_non_flowing_species(num_mins)
    !>                    call solid_chem_init%set_reactive_zone(reactive_zone)
    !>                    call solid_chem_init%allocate_concentrations()
    !>                    do i=1,num_mins
    !>                        backspace(5)
    !>                        read(5,*) reactive_zone%non_flowing_species(i)%name, target_SI, solid_chem_init%concentrations(i)
    !>                        reactive_zone%non_flowing_species(i)%cst_act_flag=.true.
    !>                        backspace(5)
    !>                    end do
    !>                end if
    !>                call reactive_zone%set_eq_reactions()
    !>                call reactive_zone%set_kin_reactions()
    !>                call reactive_zone%set_min_kin_reactions()
    !>                call reactive_zone%compute_num_reactions()
    !>                call reactive_zone%set_solid_chemistry%reactive_zone%speciation_alg()
    !>                call reactive_zone%construct_se()
    !>                call reactive_zone%compute_solid_chemistry%reactive_zone%speciation_alg_arrays()
    !>                !else if (num_species>0 .and. num_mins==0) then
    !>                    !backspace(555)
    !>                    !call gas_phase%allocate_gases(num_species)
    !>                    !do i=1,num_species
    !>                    !>    backspace(555)
    !>                    !>    read(5,*) gas_phase%gases(i)%name
    !>                    !>    backspace(555)
    !>                    !end do
    !>                    !call this%set_gas_phase(mineral_zone)
    !>                exit
    !>            else
    !>                num_eq_ph=num_eq_ph+1
    !>                call mineral%set_name(species%name)
    !>                call this%solid_chemistry%reactive_zone%chem_syst%is_mineral_in_chem_syst(mineral,flag_min,min_ind)
    !>                
    !>                if (flag_min==.true.) then
    !>                    num_nf_sp=num_nf_sp+1
    !>                    num_mins=num_mins+1
    !>                    call append_int_1D_array(min_indices,min_ind)
    !>                    call append_int_1D_array(nf_sp_indices,sp_ind)
    !>                    if (target_SI==0.0) then !> target saturation index
    !>                        num_mins_eq=num_mins_eq+1
    !>                    end if
    !>                else
    !>                    call gas%set_name(eq_phase_str)
    !>                    call this%solid_chemistry%reactive_zone%chem_syst%is_gas_in_chem_syst(gas,flag_gas,eq_phase_ind)
    !>                    if (flag_gas==.true.) then
    !>                        num_gas=num_gas+1
    !>                    else
    !>                        print *, "Pure phase ", eq_phase_str, " does not belong to chemical system"
    !>                        error stop
    !>                    end if
    !>                end if
    !>            end if
    !>        end do
    !>        !call this%set_target_solids([solid_chem_init])
    !>        call this%set_reactive_zones([reactive_zone]) !> chapuza
    !>    else if (str=='SURFACE') then
    !>        num_surf_compl=0
    !>        do 
    !>            read(5,*) species%name
    !>            if (eq_phase_str=='*') then
    !>                exit
    !>            else
    !>                call this%solid_chemistry%reactive_zone%chem_syst%is_eq_reaction_in_chem_syst(species,flag_react,eq_react_ind)
    !>                if (flag_react==.true.) then
    !>                    num_surf_compl=num_surf_compl+1
    !>                    num_nf_sp=num_nf_sp+1
    !>                    n_eq=n_eq+1
    !>                    call append_int_1D_array(eq_react_indices,eq_react_ind)
    !>                end if
    !>            end if
    !>        end do
    !>    else if (str=='EXCHANGE') then
    !>        num_exch_sp=0
    !>        do
    !>            read(5,*) species%name
    !>            if (eq_phase_str=='*') then
    !>                n_eq=compute_binomial_coeff(num_exch_sp,2) !> we assume all exchangeable cations exchange with each other
    !>                do i=1,num_exch_sp-1
    !>                    do j=i+1,num_exch_sp
    !>                        call this%solid_chemistry%reactive_zone%chem_syst%is_eq_reaction_in_chem_syst(this%solid_chemistry%reactive_zone%chem_syst%species(nf_sp_indices(i)),flag_react,eq_react_ind,this%solid_chemistry%reactive_zone%chem_syst%species(nf_sp_indices(j)))
    !>                        if (flag_react==.true.) then
    !>                            n_eq=n_eq+1
    !>                            call append_int_1D_array(eq_react_indices,eq_react_ind)
    !>                        end if
    !>                    end do
    !>                end do
    !>                exit
    !>            else
    !>                call this%solid_chemistry%reactive_zone%chem_syst%is_species_in_chem_syst(species,flag_sp,sp_ind)
    !>                if (flag_sp==.true.) then
    !>                    num_exch_sp=num_exch_sp+1
    !>                    num_nf_sp=num_nf_sp+1
    !>                    call append_int_1D_array(nf_sp_indices,sp_ind)
    !>                end if
    !>            end if
    !>        end do 
    !>    !else if (str=='REACTIVE_ZONES') then !> aqui se leen las zonas reactivas (cada una se lee a ella misma)
    !>    !>    line_count=line_count+2
    !>    !>    call reactive_zone%set_chem_syst_react_zone(this%solid_chemistry%reactive_zone%chem_syst)
    !>    !>    call reactive_zone%read_reactive_zone_Lagr(filename,line_count)
        else 
            continue
        end if
    end do
    close(5)
    !call this%initialise_target_solids(10) !> autentica chapuza ad hoc
    !!print *, this%target_solids(1)%reactive_zone%num_non_flowing_species
    !!print *, this%target_solids(1)%reactive_zone%non_flowing_species(1)%name
    !call this%initialise_target_waters_init(initial_waters)
    !call this%initialise_ext_waters()
    !!print *, associated(this%solid_chemistry%reactive_zone%chem_syst%kin_reacts(1)%kin_reaction)
    !!print *, associated(this%reactive_zones(1)%kin_reactions(1)%kin_reaction)
    !!select type (kin_react=>this%solid_chemistry%reactive_zone%chem_syst%kin_reacts(1)%kin_reaction)
    !!type is (lin_kin_reaction_c)
    !!>    !print *, kin_react%lambda, kin_react%species(1)%name
    !!end select
end subroutine