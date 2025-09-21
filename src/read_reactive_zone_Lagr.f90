!!> Reads reactive zone: "non-flowing" species
!subroutine read_reactive_zone_Lagr(this,filename,irec)
!    use reactive_zone_Lagr_m, only: reactive_zone_c, species_c, mineral_c, gas_c, eq_reaction_c, kin_params_c, kin_reaction_c,&
!     kin_reaction_poly_c, chem_system_c, speciation_algebra_c
!    use special_fcts_m
!    use array_ops_m
!    implicit none
!    class(reactive_zone_c) :: this
!    character(len=*), intent(in) :: filename !> name of input file
!    integer(kind=4), intent(in) :: irec !> record of filename where reactive zone information begins
!    
!    real(kind=8), allocatable :: Sk(:,:),logK(:),gamma_1(:),gamma_2(:),c2nc_init_guess(:),concentrations_init(:)
!    integer(kind=4) :: i,j,k,l,num_surf_compl,num_exch_sp,num_mins,num_mins_eq,num_nf_sp,sp_ind,niwtype, nbwtype, nrwtype,min_ind,&
!    eq_react_ind,num_gas,num_eq_ph,chem_syst_type,n_eq,n_k,num_react_zones,nischem,solid_chem_ind,niter,eq_phase_ind
!    integer(kind=4), allocatable :: eq_react_indices(:),min_eq_ind(:),n_tar_sol(:),min_indices(:),nf_sp_indices(:)
!    real(kind=8) :: aux,conc,temp,target_SI,max_moles
!    real(kind=8), parameter :: tol=1d-9
!    character(len=256) :: file_eq_react,eq_react_str,str,str1,str2,str3,str4,str5,path_name,Monod_name,file_kin_params,&
!    chem_syst_file,eq_phase_str
!    logical :: flag_min,flag_gas,flag_sp,flag_react
!    
!    character(len=256), allocatable :: aq_species_str(:),prim_species_str(:),cst_act_species_str(:),minerals_str(:),&
!    solid_species_str(:),kin_react_names(:)
!    type(species_c) :: species,water,constrain
!    type(species_c), allocatable :: aq_species(:),cst_act_species(:),prim_species(:),non_flowing_species(:),species_array(:)
!    type(species_c), allocatable, target :: aq_species_tgt(:)
!    type(mineral_c) :: mineral
!    type(mineral_c), allocatable :: minerals_eq(:),minerals(:)
!    type(gas_c) :: gas
!    type(gas_c), allocatable :: gases(:)
!    class(kin_params_c), pointer :: p_kin_params=>null()
!    class(kin_reaction_c), pointer :: p_kin_react=>null()
!    type(kin_reaction_poly_c) :: kin_react_ptr
!    class(kin_reaction_poly_c), allocatable :: kin_reacts(:)
!    type(eq_reaction_c) :: eq_react
!    type(eq_reaction_c), allocatable :: eq_reacts(:)
!    
!    class(chem_system_c), pointer :: p_chem_syst=>null()
!    type(chem_system_c) :: my_chem_syst
!    
!    type(speciation_algebra_c) :: spec_alg
!   
!        
!    open(unit=555,file=filename,status='old',action='read',form='formatted')
!   
!    num_nf_sp=0
!    do
!        read(555,*) str
!        !print *, str
!        if (str=='end') exit
!        !if (str=='EQUILIBRIUM_PHAseS') then
!        !>    num_mins=0
!        !>    num_mins_eq=0
!        !>    num_gas=0
!        !>    num_eq_ph=0
!        !>    do
!        !>            read(5,*) species%name, target_SI
!        !>            if (species%name=='*') then
!        !>                if (num_mins>0 .and. num_gas==0) then
!        !>                    backspace(555)
!        !>                    call mineral_zone%allocate_minerals(num_mins)
!        !>                    do i=1,num_mins
!        !>                        backspace(555)
!        !>                        read(5,*) mineral_zone%minerals(i)%name, target_SI
!        !>                        backspace(555)
!        !>                    end do
!        !>                    call this%set_mineral_zone(mineral_zone)
!        !>                end if
!        !>                !else if (num_species>0 .and. num_mins==0) then
!        !>                    !backspace(555)
!        !>                    !call gas_phase%allocate_gases(num_species)
!        !>                    !do i=1,num_species
!        !>                    !>    backspace(555)
!        !>                    !>    read(5,*) gas_phase%gases(i)%name
!        !>                    !>    backspace(555)
!        !>                    !end do
!        !>                    !call this%set_gas_phase(mineral_zone)
!        !>                exit
!        !>            else
!        !>                num_eq_ph=num_eq_ph+1
!        !>                call mineral%set_name(species%name)
!        !>                call this%solid_chemistry%reactive_zone%chem_syst%is_mineral_in_chem_syst(mineral,flag_min,min_ind)
!        !>                !call this%solid_chemistry%reactive_zone%chem_syst%is_species_in_chem_syst(species,flag_sp,sp_ind)
!        !>                call this%solid_chemistry%reactive_zone%chem_syst%is_eq_reaction_in_chem_syst(species,flag_react,eq_react_ind)
!        !>                if (flag_react==.true. .and. flag_min==.true.) then
!        !>                    num_nf_sp=num_nf_sp+1
!        !>                    num_mins=num_mins+1
!        !>                    call append_int_1D_array(min_indices,min_ind)
!        !>                    call append_int_1D_array(nf_sp_indices,sp_ind)
!        !>                    call append_int_1D_array(eq_react_indices,eq_react_ind)
!        !>                    !call nf_species%set_name(mineral%name)
!        !>                    !call nf_species%set_cst_act_flag(.true.)
!        !>                    !concentrations_init=[concentrations_init,max_moles]
!        !>                    !call mineral%set_default_values()
!        !>                    !call mineral%set_name(eq_phase_str)
!        !>                    !minerals=[minerals,mineral] !> autentica chapuza
!        !>                    !call min_zone%set_mineral_eq(mineral,min_ind)
!        !>                    if (target_SI==0.0) then !> target saturation index
!        !>                        num_mins_eq=num_mins_eq+1
!        !>                    end if
!        !>                else
!        !>                    call gas%set_name(eq_phase_str)
!        !>                    call this%solid_chemistry%reactive_zone%chem_syst%is_gas_in_chem_syst(gas,flag_gas,eq_phase_ind)
!        !>                    if (flag_gas==.true.) then
!        !>                        num_gas=num_gas+1
!        !>                        !call gas%set_name(eq_phase_str)
!        !>                        !gases=[gases,gas] !> autentica chapuza
!        !>                    else
!        !>                        print *, "Pure phase ", eq_phase_str, " does not belong to chemical system"
!        !>                        error stop
!        !>                    end if
!        !>                end if
!        !>            end if
!        !>        end do
!        !else if (str=='SURFACE') then
!        !>    num_surf_compl=0
!        !>    do 
!        !>        read(555,*) species%name
!        !>        if (eq_phase_str=='*') then
!        !>            !if (num_surf_compl>0) then
!        !>            !>    
!        !>            !end if
!        !>            exit
!        !>        else
!        !>            call this%solid_chemistry%reactive_zone%chem_syst%is_eq_reaction_in_chem_syst(species,flag_react,eq_react_ind)
!        !>            if (flag_react==.true.) then
!        !>                num_surf_compl=num_surf_compl+1
!        !>                num_nf_sp=num_nf_sp+1
!        !>                n_eq=n_eq+1
!        !>                call append_int_1D_array(eq_react_indices,eq_react_ind)
!        !>            end if
!        !>        end if
!        !>    end do
!        !>    !deallocate(eq_react_indices)
!        !else if (str=='EXCHANGE') then
!        !>    num_exch_sp=0
!        !>    do
!        !>        read(555,*) species%name
!        !>        if (eq_phase_str=='*') then
!        !>            n_eq=compute_binomial_coeff(num_exch_sp,2) !> we assume all exchangeable cations exchange with each other
!        !>            do i=1,num_exch_sp-1
!        !>                do j=i+1,num_exch_sp
!        !>                    call this%solid_chemistry%reactive_zone%chem_syst%is_eq_reaction_in_chem_syst(this%solid_chemistry%reactive_zone%chem_syst%species(nf_sp_indices(i)),flag_react,eq_react_ind,this%solid_chemistry%reactive_zone%chem_syst%species(nf_sp_indices(j)))
!        !>                    if (flag_react==.true.) then
!        !>                        n_eq=n_eq+1
!        !>                        call append_int_1D_array(eq_react_indices,eq_react_ind)
!        !>                    end if
!        !>                end do
!        !>            end do
!        !>            !n_eq=n_eq+ceiling(num_exch_sp/2d0)
!        !>            !if (num_surf_compl>0) then
!        !>            !>    
!        !>            !end if
!        !>            exit
!        !>        else
!        !>            call this%solid_chemistry%reactive_zone%chem_syst%is_species_in_chem_syst(species,flag_sp,sp_ind)
!        !>            if (flag_sp==.true.) then
!        !>                num_exch_sp=num_exch_sp+1
!        !>                num_nf_sp=num_nf_sp+1
!        !>                call append_int_1D_array(nf_sp_indices,sp_ind)
!        !>            end if
!        !>        end if
!        !>    end do 
!        !else
!        !>    continue
!        !end if
!    end do
!    close(555)
!    !if (n_e>0 .and. n_k>0) then
!    !>    !this%reactive_zones(l)=>reactive_zone_eq_kin
!    !>    !this%solid_chemistry%reactive_zone%chem_syst=>my_chem_syst_eq_kin
!    !>    this%solid_chemistry%reactive_zone%chem_syst=>my_chem_syst_eq_kin
!    !else if (n_e>0) then
!    !>    !this%reactive_zones(l)=>reactive_zone_eq
!    !>    !this%solid_chemistry%reactive_zone%chem_syst=>my_chem_syst_eq
!    !>    this%solid_chemistry%reactive_zone%chem_syst=>my_chem_syst_eq
!    !else if (n_k>0) then
!    !>    !this%reactive_zones(l)=>reactive_zone_kin
!    !>    !this%solid_chemistry%reactive_zone%chem_syst=>my_chem_syst_kin
!    !>    this%solid_chemistry%reactive_zone%chem_syst=>my_chem_syst_kin
!    !else
!    !>    error stop "No reactions involved"
!    !end if
!    !call this%solid_chemistry%reactive_zone%chem_syst%set_aq_species(aq_species)
!    !call this%aq_phase%set_num_species()
!    !call this%solid_chemistry%reactive_zone%chem_syst%set_prim_species(prim_species)
!    !call this%solid_chemistry%reactive_zone%chem_syst%compute_num_prim_species()
!    !call this%solid_chemistry%reactive_zone%chem_syst%set_cst_act_species(cst_act_species)
!    !call this%solid_chemistry%reactive_zone%chem_syst%set_num_cst_act_species()
!    !call this%solid_chemistry%reactive_zone%chem_syst%set_solid_species(solids)
!    !call this%solid_chemistry%reactive_zone%chem_syst%set_num_solid_species()
!    !call this%solid_chemistry%reactive_zone%chem_syst%set_num_var_act_species()
!    !call this%solid_chemistry%reactive_zone%chem_syst%set_var_act_species()
!    !call this%solid_chemistry%reactive_zone%chem_syst%set_species()
!    
!    
!    !call this%set_mineral_zone(mineral_zone)
!    call this%set_non_flowing_species_from_chem_syst(nf_sp_indices)
!    
!    
!    !call this%set_num_eq_reactions(n_eq)
!    call this%set_eq_reactions(eq_react_indices)
!    !call this%set_kin_reactions(kin_reacts)
!                    !!!print *, chem_syst%num_eq_reactions, chem_syst%num_kin_reactions
!                    !!call my_chem_syst%compute_num_kin_reactions()
!    !call this%compute_num_reactions()
!    
!    
!                    !!!print *, chem_syst%get_num_eq_reactions(), chem_syst%get_num_kin_reactions()
!                    !!!print *, chem_syst%mineral_zone%num_minerals, chem_syst%mineral_zone%num_minerals_eq
!    call this%set_stoich_mat_react_zone()
!    
!    !call spec_alg%set_dimensions(size(prim_species)+n_eq+size(cst_act_species),n_eq,size(cst_act_species))
!    !call spec_alg%compute_arrays(this%se,this%get_eq_csts(),this%Sk)
!    !call this%set_solid_chemistry%reactive_zone%speciation_alg(spec_alg)             
!    !call aq_spec%set_solid_chemistry%reactive_zone%speciation_alg(spec_alg)
!    !call aq_spec%set_aq_species(aq_species)
!    !call aq_spec%allocate_act_coeffs()
!                    
!    !call this%set_aq_speciation(aq_spec)
!    !call this%compute_aq_solid_chemistry%reactive_zone%speciation_alg()
!    !call this%initialise_solid_chem(init_solid_chems,n_tar_sol)
!    !call this%initialise_aqueous_chem(initial_water_types,n_tar)
!    !call this%initialise_water_sources()
!end subroutine