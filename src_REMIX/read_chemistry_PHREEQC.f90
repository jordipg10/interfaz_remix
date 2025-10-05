! !> Based on PHREEQC
! subroutine read_chemistry_PHREEQC(this,path,filename)
!     use chemistry_m, only: chemistry_c
!     use aqueous_chemistry_m, only: aqueous_chemistry_c, solid_chemistry_c, reactive_zone_c
!     use chem_system_m, only: eq_reaction_c, kin_params_c, kin_reaction_c, kin_reaction_poly_c, eq_reaction_c, species_c, & 
!         aq_species_c, mineral_c, chem_system_c, gas_c, aq_phase_c, aq_phase_c
!     use special_fcts_m, only: compute_binomial_coeff
!     use array_ops_m, only: append_int_1D_array
!     implicit none
!     class(chemistry_c) :: this
!     character(len=*), intent(in) :: path
!     character(len=*), intent(in) :: filename
    
!     real(kind=8), allocatable :: Sk(:,:),logK(:),gamma_1(:),gamma_2(:)
!     integer(kind=4) :: line_count,i,j,k,l,num_surf_compl,num_exch_sp,num_mins,num_mins_eq,num_nf_sp,sp_ind,niwtype, nbwtype, &
!     nrwtype,min_ind,eq_react_ind,num_gas,num_eq_ph,chem_syst_type,n_eq,n_k,num_react_zones,aq_sp_ind,solid_chem_ind,niter,&
!     eq_phase_ind
!     integer(kind=4), allocatable :: eq_react_indices(:),min_eq_ind(:),n_tar_sol(:),min_indices(:),nf_sp_indices(:),&
!     aq_sp_indices(:)
!     real(kind=8) :: aux,conc_init,temp,value,target_SI
!     character(len=256) :: str,str1,str2,str3,str4,str5,path_name,Monod_name,file_kin_params,chem_syst_file,sol_str,eq_phase_str
!     logical :: flag_sp,flag_min,flag_gas,flag_react
    
!     type(aqueous_chemistry_c) :: initial_water
!     type(aqueous_chemistry_c), allocatable :: initial_waters(:)
!     type(solid_chemistry_c) :: solid_chem_init
    
!     character(len=256), allocatable :: aq_species_str(:),prim_species_str(:),cst_act_species_str(:),minerals_str(:),&
!     solid_species_str(:),kin_react_names(:)
!     type(species_c) :: species,water,constrain
!     type(species_c), allocatable :: aq_species_array(:),cst_act_species(:),prim_species(:)
!     type(aq_species_c) :: aq_species
!     type(mineral_c) :: mineral
!     type(mineral_c), allocatable :: minerals(:),minerals_eq(:)
!     type(reactive_zone_c) :: reactive_zone
!     type(gas_c) :: gas
!     type(gas_c), allocatable :: gases(:)
!     class(kin_params_c), pointer :: p_kin_params=>null()
!     class(kin_reaction_c), pointer :: p_kin_react=>null()
!     type(kin_reaction_poly_c) :: kin_react_ptr
!     class(kin_reaction_poly_c), allocatable :: kin_reacts(:)
!     type(eq_reaction_c) :: eq_react
!     type(eq_reaction_c), allocatable :: eq_reacts(:)
    
!     type(chem_system_c) :: my_chem_syst
    
!     open(unit=5,file=filename,status='old',action='read')
!     line_count=1
!     do
!         read(5,*) str
!         print *, str
!         if (str=='end') exit
!         if (str=='CHEMICAL_SYSTEM') then
!             line_count=line_count+2
!             read(5,*) chem_syst_file
!             !call my_chem_syst%read_chem_system(path,chem_syst_file) !> chemical system associated to chemistry class
!             call this%set_chem_syst(my_chem_syst)
!             !print *, this%chem_syst%lin_kin_reacts(1)%lambda
!             !print *, "hola 1"
!             !select type (kin_react=>this%chem_syst%kin_reacts(1)%kin_reaction)
!             !type is (lin_kin_reaction_c)
!             !>    call kin_react%set_lambda(1d0)
!             !>    print *, kin_react%lambda
!             !>    print *, "hola 2"
!             !end select
!         else if (str=='SOLUTION') then
!             line_count=line_count+3
!             call initial_water%set_aq_phase(this%chem_syst%aq_phase)
!             call initial_water%allocate_conc_aq_species()
!             call initial_water%allocate_log_act_coeffs_aq_chem()
!         !> The loop below should belong to aqueous chemistry class
!             do
!                 line_count=line_count+1
!                 read(5,*) sol_str, value
!                 if (sol_str=='*') then
!                     exit
!                 else if (sol_str=='temp') then
!                     call initial_water%set_temp(value)
!                 else if (sol_str=='pressure') then
!                     call initial_water%set_pressure(value)
!                 else if (sol_str=='pH') then
!                     call initial_water%set_pH(value)
!                 else if (sol_str=='pe') then
!                     call initial_water%set_pe(value)
!                 else if (sol_str=='density') then
!                     call initial_water%set_density(value)
!                 else
!                     call aq_species%set_name(sol_str)
!                     call aq_species%set_cst_act_flag(.false.)
!                     !call this%chem_syst%is_aq_species_in_chem_syst(aq_species,flag_sp,aq_sp_ind)
!                     if (flag_sp.eqv..true.) then
!                         call initial_water%set_conc_single_species(value,aq_sp_ind)
!                     else
!                         error stop "Aqueous species does not belong to chemical system"
!                     end if
!                 end if
!             end do
!             initial_waters=[initial_waters,initial_water] !> autentica chapuza
!             deallocate(initial_water%concentrations,initial_water%log_act_coeffs,initial_water%log_Jacobian_act_coeffs)
!             !call this%initialise_target_waters(initial_waters)
!         else if (str=='EQUILIBRIUM_PHAseS') then
!             num_mins=0
!             num_mins_eq=0
!             num_gas=0
!             num_eq_ph=0
!             allocate(min_indices(0),nf_sp_indices(0))
!             do
!                 read(5,*) species%name, target_SI, conc_init
!                 if (species%name=='*') then
!                     call reactive_zone%set_chem_syst_react_zone(this%chem_syst)
!                     if (num_mins>0 .and. num_gas==0) then
!                         backspace(5)
!                         call reactive_zone%allocate_non_flowing_species(num_mins)
!                         call solid_chem_init%set_reactive_zone(reactive_zone)
!                         call solid_chem_init%allocate_conc_solids()
!                         do i=1,num_mins
!                             backspace(5)
!                             read(5,*) reactive_zone%non_flowing_species(i)%name, target_SI, solid_chem_init%concentrations(i)
!                             reactive_zone%non_flowing_species(i)%cst_act_flag=.true.
!                             backspace(5)
!                         end do
!                     end if
!                     call reactive_zone%set_eq_reactions()
!                     !call reactive_zone%set_kin_reactions()
!                     !call reactive_zone%set_min_kin_reactions()
!                     !call reactive_zone%compute_num_reactions()
!                     !call reactive_zone%set_speciation_alg_dimensions()
!                     !call reactive_zone%set_se()
!                     !call reactive_zone%compute_speciation_alg_arrays()
!                     !else if (num_species>0 .and. num_mins==0) then
!                         !backspace(555)
!                         !call gas_phase%allocate_gases(num_species)
!                         !do i=1,num_species
!                         !>    backspace(555)
!                         !>    read(5,*) gas_phase%gases(i)%name
!                         !>    backspace(555)
!                         !end do
!                         !call this%set_gas_phase(mineral_zone)
!                     exit
!                 else
!                     num_eq_ph=num_eq_ph+1
!                     call mineral%set_phase_name(species%name)
!                     call this%chem_syst%is_mineral_in_chem_syst(mineral,flag_min,min_ind)
                    
!                     if (flag_min.eqv..true.) then
!                         num_nf_sp=num_nf_sp+1
!                         num_mins=num_mins+1
!                         call append_int_1D_array(min_indices,min_ind)
!                         call append_int_1D_array(nf_sp_indices,sp_ind)
!                         if (target_SI==0.0) then !> target saturation index
!                             num_mins_eq=num_mins_eq+1
!                         end if
!                     else
!                         call gas%set_name(eq_phase_str)
!                         call this%chem_syst%gas_phase%is_gas_in_gas_phase(gas,flag_gas,eq_phase_ind)
!                         if (flag_gas.eqv..true.) then
!                             num_gas=num_gas+1
!                         else
!                             print *, "Pure phase ", eq_phase_str, " does not belong to chemical system"
!                             error stop
!                         end if
!                     end if
!                 end if
!             end do
!             !call this%set_target_solids([solid_chem_init])
!             call this%set_reactive_zones([reactive_zone]) !> chapuza
!         else if (str=='SURFACE') then
!             num_surf_compl=0
!             do 
!                 read(5,*) species%name
!                 if (eq_phase_str=='*') then
!                     exit
!                 else
!                     call this%chem_syst%is_eq_reaction_in_chem_syst(species%name,flag_react,eq_react_ind)
!                     if (flag_react.eqv..true.) then
!                         num_surf_compl=num_surf_compl+1
!                         num_nf_sp=num_nf_sp+1
!                         n_eq=n_eq+1
!                         call append_int_1D_array(eq_react_indices,eq_react_ind)
!                     end if
!                 end if
!             end do
!         else if (str=='EXCHANGE') then
!             num_exch_sp=0
!             do
!                 read(5,*) species%name
!                 if (eq_phase_str=='*') then
!                     n_eq=compute_binomial_coeff(num_exch_sp,2) !> we assume all exchangeable cations exchange with each other
!                     do i=1,num_exch_sp-1
!                         do j=i+1,num_exch_sp
!                             !call this%chem_syst%is_eq_reaction_in_chem_syst(this%chem_syst%species(nf_sp_indices(i)),flag_react,eq_react_ind,this%chem_syst%species(nf_sp_indices(j)))
!                             if (flag_react.eqv..true.) then
!                                 n_eq=n_eq+1
!                                 call append_int_1D_array(eq_react_indices,eq_react_ind)
!                             end if
!                         end do
!                     end do
!                     exit
!                 else
!                     call this%chem_syst%is_species_in_chem_syst(species,flag_sp,sp_ind)
!                     if (flag_sp.eqv..true.) then
!                         num_exch_sp=num_exch_sp+1
!                         num_nf_sp=num_nf_sp+1
!                         call append_int_1D_array(nf_sp_indices,sp_ind)
!                     end if
!                 end if
!             end do 
!         !else if (str=='REACTIVE_ZONES') then !> aqui se leen las zonas reactivas (cada una se lee a ella misma)
!         !>    line_count=line_count+2
!         !>    call reactive_zone%set_chem_syst_react_zone(this%chem_syst)
!         !>    call reactive_zone%read_reactive_zone_Lagr(filename,line_count)
!         else 
!             continue
!         end if
!     end do
!     close(5)
!     !print *, this%target_solids(1)%reactive_zone%num_non_flowing_species
!     !print *, this%target_solids(1)%reactive_zone%non_flowing_species(1)%name
!     !print *, associated(this%chem_syst%kin_reacts(1)%kin_reaction)
!     !print *, associated(this%reactive_zones(1)%kin_reactions(1)%kin_reaction)
!     !select type (kin_react=>this%chem_syst%kin_reacts(1)%kin_reaction)
!     !type is (lin_kin_reaction_c)
!     !>    !print *, kin_react%lambda, kin_react%species(1)%name
!     !end select
! end subroutine