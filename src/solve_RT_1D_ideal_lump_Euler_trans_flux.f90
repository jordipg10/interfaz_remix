!> Solves reactive mixing problem assuming ideal conditions
!> Applies lumping to the mixing ratios of reaction amounts
!> Computes concentrations, activities, activity coefficients, reaction rates and volumetric fractions of minerals (if present)
subroutine solve_RT_ideal_lump_Euler_trans_flux_1D(this,dir,root)
    use RT_m, only: RT_1D_transient_c, RT_1D_stat_c, RT_c
    use aqueous_chemistry_m, only: aqueous_chemistry_c, &
        reactive_mixing_iter_EE_eq_kin_ideal, reactive_mixing_iter_EE_kin, &
        reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2, &
        reactive_mixing_iter_EI_kin_anal_ideal_opt2, &
        compute_r_tilde_impl_opt1, mixing_iter_comp_ideal, compute_c_mix
    implicit none
!> Arguments
    class(RT_1D_transient_c) :: this !> reactive transport object
    character(len=*), intent(in) :: dir !> Directory for output files
    character(len=*), intent(in) :: root !> root path
    !integer(kind=4), intent(in) :: unit !> file unit
    !class(real_array_c), intent(in) :: mixing_ratios !> mixing ratios matrix for concentrations
    !class(real_array_c), intent(in) :: mixing_ratios_R_init !> initial mixing ratios matrix for reaction amounts
    !class(int_array_c), intent(in) :: this%transport%mix_conc_indices !> matrix that contains indices of target waters that mix with each target water
    !class(int_array_c), intent(in) :: this%transport%mix_react_indices !> matrix that contains indices of domain target waters that mix with each target water
    !real(kind=8), intent(in) :: F_mat(:) !> storage matrix (diagonal)
    !class(time_discr_c), intent(in) :: time_discr_tpt !> time discretisation object (used to solve transport)
    !real(kind=8), intent(in) :: theta_r !> integration method for chemical reactions
    !class(real_array_c), intent(inout) :: mixing_ratios_R !> final mixing ratios matrix for reaction amounts
!> Variables
    integer(kind=4) :: i !> counter target waters
    integer(kind=4) :: j !> counter target solids 
    integer(kind=4) :: l !> counter reactive zones
    integer(kind=4) :: k !> counter time steps
    integer(kind=4) :: kk !> counter time steps chem_out_options
    integer(kind=4) :: ii !> counter target waters chem_out_options
    integer(kind=4) :: num_tar_wat !> number of target waters
    integer(kind=4) :: num_tar_sol !> number of target solids
    integer(kind=4) :: n_p !> number of primary species
    integer(kind=4) :: n_v !> number of variable activity species
    integer(kind=4) :: n_v_aq !> number of aqueous variable activity species
    integer(kind=4) :: n_v_aq_2 !> number of aqueous secondary variable activity species
    integer(kind=4) :: mix_ind !> mixing waters starting index
    integer(kind=4) :: num_can_vec !> number of canonical vectors
    integer(kind=4) :: num_non_can_vec !> number of canonical vectors
    integer(kind=4) :: unit !> unit for output file
    integer(kind=4) :: k_flux !> counter flux change
    integer(kind=4), allocatable :: tar_gas_indices(:) !> indices target gases in each reactive zone
    integer(kind=4), allocatable :: tar_sol_indices(:) !> indices target solids in each reactive zone
    integer(kind=4), allocatable :: tar_wat_indices(:) !> indices target waters in each reactive zone
    integer(kind=4), allocatable :: perm(:) !> permutation vector aqueous concentrations
    integer(kind=4), allocatable :: ind_can_vec(:) !> indices of canonical vectors in mixing ratios
    integer(kind=4), allocatable :: ind_non_can_vec(:) !> indices of non-canonical vectors in mixing ratios
    REAL(KIND=8) :: time !> time
    REAL(KIND=8) :: Delta_t !> time step
    REAL(KIND=8) :: theta !> time weighting factor
    REAL(KIND=8) :: y !> sum of upstream mixing ratios
    REAL(KIND=8), allocatable :: c_hat(:) !> variable activity mobile species concentrations after mixing
    real(kind=8), allocatable :: rk_tilde_up(:) !> kinetic reaction rates after mixing
    real(kind=8), allocatable :: rk_tilde_down(:) !> kinetic reaction rates after mixing
    real(kind=8), allocatable :: rk_tilde(:) !> kinetic reaction rates
    REAL(KIND=8), allocatable :: conc_old(:,:) !> concentrations before mixing
    REAL(KIND=8), allocatable :: conc_nc(:) !> concentrations variable activity species
    REAL(KIND=8), allocatable :: conc_comp(:) !> concentrations components
    type(aqueous_chemistry_c), allocatable :: waters_new(:) !> target waters time step k+1
    type(aqueous_chemistry_c), allocatable :: waters_old(:) !> target waters time step k
    type(aqueous_chemistry_c), allocatable :: waters_old_old(:) !> target waters time step k-1
    type(aqueous_chemistry_c), allocatable :: mix_waters(:) !> mixing waters
    integer(kind=4) :: num_lump !> number of lumpings
    integer(kind=4) :: j_mix !> loop index for mixing waters
    integer(kind=4) :: num_mix_loc !> number of mixing waters
    real(kind=8), allocatable :: conc_old_mix(:,:) !> conc_old from each mixing water
    integer(kind=4), allocatable :: ind_aq_sp_mix(:,:) !> indices_aq_species from each mixing water
!> Procedure pointers
    !> reactive mixing subroutines
    procedure(reactive_mixing_iter_EE_eq_kin_ideal), pointer :: p_solver=>null()
    procedure(compute_r_tilde_impl_opt1), pointer :: compute_r_tilde=>null()
!> We initialise target waters
    waters_old=this%chemistry%waters
    waters_old_old=waters_old
    waters_new=waters_old

select type (this)
type is (RT_1D_transient_c)
    time=0d0
    k_flux=1 !> counter flux change
    kk=2 !> counter time steps chem_out_options
    ii=1 !> counter target waters chem_out_options
    !print *, this%chemistry%chem_out_options%time_steps
    this%chemistry%chem_out_options%time_steps(this%chemistry%chem_out_options%num_time_steps)=this%transport%time_discr%Num_time !> chapuza
!> Chapuza
    ! if (this%chemistry%chem_syst%num_kin_reacts==this%chemistry%chem_syst%num_aq_kin_reacts) then
    !     allocate(rk_tilde(this%chemistry%chem_syst%num_kin_reacts))
    !     allocate(rk(this%chemistry%chem_syst%num_kin_reacts))
    ! else
    !     continue
    ! end if
    unit=7 !> arbitrary unit for output file
    open(unit,file=dir//root//'.output',form="formatted",access="sequential",status="unknown")
!> We select reactive mixing subroutine depending on the nature of the chemical system and the methods to compute Jacobians and integrate in time
    if (this%chemistry%chem_syst%num_kin_reacts>0 .and. this%chemistry%chem_syst%speciation_alg%num_eq_reactions>0) then !> equilibrium and kinetic reactions
        if (this%int_method_chem_reacts==1) then !> Euler explicit
            theta=0d0 !> Euler explicit
            !p_solver=>reactive_mixing_iter_EE_eq_kin_ideal_lump
            compute_r_tilde=>compute_r_tilde_impl_opt1
            mix_ind=1 !> we mix all waters
            !mix_waters_indices=>this%transport%mix_conc_indices
        else if (this%chemistry%Jac_opt==1) then !> analytical Jacobian
            mix_ind=2 !> we mix all waters except one
            if (this%int_method_chem_reacts==2) then !> Euler fully implicit
                theta=1d0 !> Euler fully implicit
            else if (this%int_method_chem_reacts==3) then !> Crank-Nicolson
                theta=5d-1 !> Crank-Nicolson parameter
            else
                error stop "Integration method for chemical reactions not implemented yet"
            end if
            if (this%chemistry%rk_avg_opt==1) then
                !p_solver=>reactive_mixing_iter_EI_eq_kin_anal_ideal_opt1
            else if (this%chemistry%rk_avg_opt==2) then
                !p_solver=>reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2
            else
                error stop "rk average option not implemented yet"
            end if
            !if (this%chemistry%cons_opt==1) then !> rk_tilde explicit
            !    !p_solver=>reactive_mixing_iter_CN_eq_kin_ideal_opt1
            !    compute_rk_tilde=>compute_rk_tilde_expl
            !else if (this%chemistry%cons_opt==2) then !> rk_tilde implicit
            !    if (this%chemistry%r_down_opt==1) then
            !        compute_rk_tilde=>compute_rk_tilde_impl_opt1
            !    else if (this%chemistry%r_down_opt==2) then
            !        compute_rk_tilde=>compute_rk_tilde_impl_opt2
            !    else if (this%chemistry%r_down_opt==3) then
            !        compute_rk_tilde=>compute_rk_tilde_impl_opt3
            !    else
            !        error stop "rk down option not implemented yet"
            !    end if
            !    !p_solver=>reactive_mixing_iter_CN_eq_kin_ideal_opt2
            !    !compute_rk_tilde=>compute_rk_tilde_impl_opt2
            !    !compute_rk_tilde=>compute_rk_tilde_impl_opt3
            !else
            !    error stop "WMA consistent option not implemented yet"
            !end if
        ! else if (this%int_method_chem_reacts==3 .and. this%chemistry%Jac_opt==1) then !> Crank-Nicolson, analytical Jacobian
        !     theta=5d-1 !> Crank-Nicolson parameter
        !     mix_ind=2 !> we mix all waters except one
        !     if (this%chemistry%rk_avg_opt==1) then
        !         !p_solver=>reactive_mixing_iter_EI_eq_kin_anal_ideal_opt1
        !     else if (this%chemistry%rk_avg_opt==2) then
        !         p_solver=>reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2
        !     else
        !         error stop "rk option not implemented yet"
        !     end if
        !     if (this%chemistry%cons_opt==1) then !> rk_tilde explicit
        !         !p_solver=>reactive_mixing_iter_CN_eq_kin_ideal_opt1
        !         compute_rk_tilde=>compute_rk_tilde_expl
        !     else if (this%chemistry%cons_opt==2) then !> rk_tilde implicit
        !         !p_solver=>reactive_mixing_iter_CN_eq_kin_ideal_opt2
        !         if (this%chemistry%r_down_opt==1) then
        !             compute_rk_tilde=>compute_rk_tilde_impl_opt1
        !         else if (this%chemistry%r_down_opt==2) then
        !             compute_rk_tilde=>compute_rk_tilde_impl_opt2
        !         else if (this%chemistry%r_down_opt==3) then
        !             compute_rk_tilde=>compute_rk_tilde_impl_opt3
        !         else
        !             error stop "rk down option not implemented yet"
        !         end if
            ! else
            !     error stop "WMA consistent option not implemented yet"
            ! end if
        else
            error stop "Integration method for chemical reactions not implemented yet"
        end if
    else if (this%chemistry%chem_syst%speciation_alg%num_eq_reactions>0) then !> only equilibrium reactions
        !p_solver=>mixing_iter_comp_ideal !> only equilibrium reactions
    else !> only kinetic reactions
        if (this%int_method_chem_reacts==1) then !> Euler explicit
            theta=0d0 !> Euler explicit
            p_solver=>reactive_mixing_iter_EE_kin
            !compute_rk_tilde=>compute_rk_tilde_expl
            mix_ind=1 !> we mix all waters
        else if (this%chemistry%Jac_opt==1) then !> analytical Jacobian
            if (this%int_method_chem_reacts==2) then !> Euler fully implicit
                theta=1d0
            else if (this%int_method_chem_reacts==3) then !> Crank-Nicolson
                theta=5d-1 !> Crank-Nicolson parameter
            else
                error stop "Integration method for chemical reactions not implemented yet"
            end if
            if (this%chemistry%rk_avg_opt==1) then
                !p_solver=>reactive_mixing_iter_EI_kin_anal_ideal_opt1
            else if (this%chemistry%rk_avg_opt==2) then
                !p_solver=>reactive_mixing_iter_EI_kin_anal_ideal_opt2
            else
                error stop "rk average option not implemented yet"
            end if
            mix_ind=2 !> we mix all waters except one
            !if (this%chemistry%cons_opt==1) then !> mix kinetic reaction rates from previous time step
            !    !p_solver=>reactive_mixing_iter_EfI_eq_kin_anal_ideal_opt1
            !    compute_rk_tilde=>compute_rk_tilde_expl
            !else if (this%chemistry%cons_opt==2) then !> solve upstream to downstream
            !    !p_solver=>reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2
            !    if (this%chemistry%r_down_opt==1) then
            !        compute_rk_tilde=>compute_rk_tilde_impl_opt1
            !    else if (this%chemistry%r_down_opt==2) then
            !       compute_rk_tilde=>compute_rk_tilde_impl_opt2
            !    else if (this%chemistry%r_down_opt==3) then
            !       compute_rk_tilde=>compute_rk_tilde_impl_opt3
            !    else if (this%chemistry%r_down_opt==4) then
            !        !compute_rk_tilde=>compute_rk_tilde_impl_opt4
            !    else
            !        error stop "rk down option not implemented yet"
            !    end if
            !else
            !    error stop "WMA consistent option not implemented yet"
            !end if
        ! else if (this%int_method_chem_reacts==3 .and. this%chemistry%Jac_opt==1) then !> Crank-Nicolson, analytical Jacobian
        !     theta=5d-1 !> Crank-Nicolson parameter
        !     mix_ind=2 !> we mix all waters except one
        !     if (this%chemistry%rk_avg_opt==1) then
        !         !p_solver=>reactive_mixing_iter_EI_eq_kin_anal_ideal_opt1
        !     else if (this%chemistry%rk_avg_opt==2) then
        !         p_solver=>reactive_mixing_iter_EI_kin_anal_ideal_opt2
        !     else
        !         error stop "rk option not implemented yet"
        !     end if
        !     if (this%chemistry%cons_opt==1) then !> mix kinetic reaction rates from previous time step
        !         !p_solver=>reactive_mixing_iter_CN_eq_kin_ideal_opt1
        !         !compute_rk_tilde=>compute_rk_tilde_expl
        !     else if (this%chemistry%cons_opt==2) then !> solve upstream to downstream
        !         !p_solver=>reactive_mixing_iter_CN_eq_kin_ideal_opt2
        !         compute_rk_tilde=>compute_rk_tilde_impl_opt1
        !         !compute_rk_tilde=>compute_rk_tilde_impl_opt2
        !         !compute_rk_tilde=>compute_rk_tilde_impl_opt3
        !     else
        !         error stop "WMA option not implemented yet"
        !     end if
        else
            error stop "Integration method for chemical reactions not implemented yet"
        end if
    end if
    !> We get the canonical vectors in the mixing ratios array
    call this%transport%mixing_ratios_conc%get_can_vec(this%chemistry%CV_params%abs_tol,num_can_vec,ind_can_vec,num_non_can_vec,&
        ind_non_can_vec)
    !> We set old concentrations in Dirichlet boundary target waters (chapuza)
    do i=1,this%chemistry%num_waters
        call waters_new(i)%set_conc_old() !> 
        call waters_new(i)%set_conc_old_old() !> 
        call waters_new(i)%solid_chemistry%set_conc_old() !> 
        call waters_new(i)%solid_chemistry%set_conc_old_old() !> 
        !call initialise_iterative_method()
        !call waters_new(ind_non_can_vec(i))%update_rk_old()
        !call waters_new(ind_can_vec(i))%solid_chemistry%update_rk_old()
    end do
    !do i=1,num_non_can_vec
    !    call waters_new(ind_can_vec(i))%set_conc_old() !>
    !    call waters_new(ind_can_vec(i))%solid_chemistry%set_conc_old() !>
    !    call waters_new(ind_can_vec(i))%solid_chemistry%set_conc_old_old() !>
    !end do
        !> Time loop
        do k=1,this%transport%time_discr%Num_time
            Delta_t=this%transport%time_discr%get_Delta_t(k)
            time=time+Delta_t
            if (k==this%chemistry%chem_out_options%time_steps(kk)) then
               write(unit,"(2x,'t = ',*(ES15.5))") time
               write(unit,"(A10,*(A15))") 'Water', (this%chemistry%chem_syst%aq_phase%aq_species( &
                    this%chemistry%chem_out_options%ind_aq_species(j))%name, &
                    j=1,this%chemistry%chem_out_options%num_aq_species)
            end if
            !> we check if flux has changed
            if (this%transport%tpt_props_heterog%flux_trans%time_ind(k_flux)==k) then
                call this%transport%tpt_props_heterog%update_flux_int_trans(k_flux) !> we update flux
                call this%transport%stab_params_tpt%compute_stab_params_tpt_1D(this%transport%tpt_props_heterog,&
                    this%transport%spatial_discr,Delta_t)
                call this%transport%update_mixing_ratios_Delta_t_homog(Delta_t) !> we update mixing ratios
                if (k_flux<this%transport%tpt_props_heterog%flux_trans%ntime) then
                    k_flux=k_flux+1 !> we update counter flux change
                end if
            end if
            !> We update concentrations and kinetic reaction rates in domain target waters (chapuza)
            !print *, this%chemistry%tar_wat_indices
            do i=1,this%chemistry%num_target_waters
                call waters_new(this%chemistry%tar_wat_indices(i))%update_conc_old() !> 
                call waters_new(this%chemistry%tar_wat_indices(i))%solid_chemistry%update_conc_old() !> 
                call waters_new(this%chemistry%tar_wat_indices(i))%update_rk_old()
                call waters_new(this%chemistry%tar_wat_indices(i))%solid_chemistry%update_rk_old()
            end do
            !> Domain target waters loop
            do i=1,num_non_can_vec
                n_p=this%chemistry%waters(&
                    this%chemistry%tar_wat_indices(&
                    i))%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
                n_v=this%chemistry%waters(this%chemistry%tar_wat_indices(i)&
                    )%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
                n_v_aq=this%chemistry%waters(this%chemistry%tar_wat_indices(i)&
                    )%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
                n_v_aq_2=this%chemistry%waters(this%chemistry%tar_wat_indices(i)&
                    )%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                !print *, this%chemistry%waters(this%chemistry%tar_wat_indices(i))%solid_chemistry%reactive_zone%num_solids
                
                !if (num_can_vec>0 .and. ind_non_can_vec(i)==ind_can_vec(cntr_can_vec)) then
                !    cntr_can_vec=cntr_can_vec+1 !> we update counter of canonical vectors
                !    !call waters_new(ind_non_can_vec(i))%compute_rk() !> we compute kinetic reaction rates
                !    conc_nc=waters_new(ind_non_can_vec(i))%get_conc_nc()
                !    !> chapuza
                !    !waters_new(ind_non_can_vec(i))%Rk_est=0d0
                !    !waters_new(ind_non_can_vec(i))%solid_chemistry%Rk_est=0d0
                !    continue
                !else
                    if (this%chemistry%waters(this%chemistry%tar_wat_indices(&
                    i))%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl &
                        >0) then !> variable activity species are aqueous and solid
                        !p_solver=>mixing_iter_comp_exch_ideal !> only equilibrium reactions
                    end if
                    allocate(conc_nc(n_v))
                    allocate(conc_comp(n_p))
                    !allocate(rk_tilde_up(n_v),rk_tilde_down(n_v))
                    !allocate(rk_tilde(n_v))
                    !allocate(conc_old(n_v,mixing_ratios_conc%cols(i)%dim)) !> chapuza
                    !print *, this%chemistry%tar_wat_indices
                    if (this%transport%mix_conc_indices%cols(ind_non_can_vec(i))%col_1(1)/=this%chemistry%tar_wat_indices(&
                        i)) then
                        print *, this%transport%mix_conc_indices%cols(ind_non_can_vec(i))%col_1(1),&
                            this%chemistry%tar_wat_indices(i)
                        error stop "Target waters not in the right order"
                    end if
                    !> We solve mixing caused by transport
                    !c_hat=waters_old(this%chemistry%tar_wat_indices(i))%get_conc_nc()
                    num_mix_loc=this%transport%mix_conc_indices%cols(ind_non_can_vec(i))%dim-3
                    allocate(conc_old_mix(size(waters_old(this%transport%mix_conc_indices%cols(ind_non_can_vec(i))%col_1(1))%conc_old), num_mix_loc))
                    allocate(ind_aq_sp_mix(size(waters_old(this%transport%mix_conc_indices%cols(ind_non_can_vec(i))%col_1(1))%indices_aq_species), num_mix_loc))
                    do j_mix=1,num_mix_loc
                        conc_old_mix(:,j_mix)=waters_old(this%transport%mix_conc_indices%cols(ind_non_can_vec(i))%col_1(j_mix+1))%conc_old
                        ind_aq_sp_mix(:,j_mix)=waters_old(this%transport%mix_conc_indices%cols(ind_non_can_vec(i))%col_1(j_mix+1))%indices_aq_species
                    end do
                    call compute_c_mix(waters_old(this%transport%mix_conc_indices%cols(ind_non_can_vec(i))%col_1(1)),&
                        conc_old_mix,ind_aq_sp_mix,&
                        this%transport%mixing_ratios_conc%cols(ind_non_can_vec(i))%col_1,c_hat)
                    deallocate(conc_old_mix,ind_aq_sp_mix)
                    !> chapuza
                    !allocate(mix_waters(this%chemistry%num_target_waters-num_can_vec-mix_ind+1))
                    !do j=1,this%chemistry%num_target_waters-mix_ind+1
                    !    if (j==
                    !    mix_waters(j)=waters_new(this%transport%mix_conc_indices%cols(i)%col_1(&
                    !    mix_ind+j-1))
                    !mix_waters=waters_new(this%transport%mix_react_indices%cols(ind_non_can_vec(i))%col_1(mix_ind:this%transport%mix_react_indices%cols(ind_non_can_vec(i))%dim-2))
                    !call compute_rk_tilde(mix_waters,mixing_ratios_R_init%cols(ind_non_can_vec(i))%col_1(mix_ind:),&
                    !    this%transport%mix_react_indices%cols(ind_non_can_vec(i))%col_1(this%transport%mix_react_indices%cols(ind_non_can_vec(i))%dim-1),&
                    !    this%transport%mix_react_indices%cols(ind_non_can_vec(i))%col_1(this%transport%mix_react_indices%cols(ind_non_can_vec(i))%dim),theta,Delta_t,rk_tilde)
                    !!> chapuza
                    !waters_new(this%transport%mix_react_indices%cols(ind_non_can_vec(i))%col_1(mix_ind:this%transport%mix_react_indices%cols(ind_non_can_vec(i))%dim-2))=&
                    !    mix_waters
                    !call waters_new(this%chemistry%tar_wat_indices(ind_non_can_vec(i)))%solid_chemistry%modify_mix_ratios_rk(&
                    !    mixing_ratios_R_init%cols(ind_non_can_vec(i))%col_1(1),c_mix,Delta_t,rk_tilde,&
                    !    mixing_ratios_R%cols(ind_non_can_vec(i))%col_1(1),num_lump)
                    !!> chapuza
                    !this%chemistry%num_lump=this%chemistry%num_lump+num_lump !> we update number of lumpings
                    !y=sum(mixing_ratios_R%cols(i)%col_1(2:1+this%transport%mix_conc_indices%cols(i)%col_1(this%transport%mix_conc_indices%cols(i)%dim-1)))
                    !print *, y*Delta_t
                    !do j=1,n_v_aq
                    !    if (waters_old(ind_non_can_vec(i))%ind_var_act_species(j) /= &
                    !        waters_old(ind_non_can_vec(i))%indices_aq_species(&
                    !        this%chemistry%chem_syst%aq_phase%ind_diss_solids(j))) then
                    !        c_mix(waters_old(ind_non_can_vec(i))%ind_var_act_species(j))=c_mix_aux(&
                    !            waters_old(ind_non_can_vec(i))%indices_aq_species(&
                    !            this%chemistry%chem_syst%aq_phase%ind_diss_solids(j)))
                    !    end if
                    !end do
                    !> We solve reactive mixing iteration
                     call p_solver(waters_new(this%chemistry%tar_wat_indices(i)),waters_old_old(&
                        this%chemistry%tar_wat_indices(i))%get_c1(),c_hat,this%transport%mixing_ratios_R%cols(ind_non_can_vec(i))%col_1(1),&
                        Delta_t,theta,conc_nc,conc_comp)
                    !> Compute pH, salinity, and ionic strength from updated concentrations
                    call waters_new(this%chemistry%tar_wat_indices(i))%compute_pH()
                    call waters_new(this%chemistry%tar_wat_indices(i))%compute_salinity()
                    call waters_new(this%chemistry%tar_wat_indices(i))%compute_ionic_strength()
                    !> We compute equilibrium reaction rates from mass balance equation
                    if (this%chemistry%waters(this%chemistry%tar_wat_indices(i)&
                        )%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions>0 .and. &
                        this%chemistry%waters(this%chemistry%tar_wat_indices(i))%indices_rk%num_cols>0)&
                            then
                        call waters_new(this%chemistry%tar_wat_indices(i))%compute_Re_kin(&
                            c_hat(n_p+1:n_v),Delta_t,&
                            this%transport%mixing_ratios_R%cols(ind_non_can_vec(i))%col_1(1)) !> chapuza
                    else if (this%chemistry%waters(this%chemistry%tar_wat_indices(i)&
                        )%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions>0) then
                        call waters_new(this%chemistry%tar_wat_indices(i))%compute_Re(&
                            c_hat(1:n_v),Delta_t,&
                            this%transport%mixing_ratios_R%cols(ind_non_can_vec(i))%col_1(1)) !> chapuza
                    end if
                    ! call waters_new(ind_non_can_vec(i))%compute_re_mean(c_mix(n_p+1:n_p+n_v_aq_2),Delta_t,&
                    !     rk_tilde)
                    !> We compute solid chemistry state variables
                    if (associated(waters_new(this%chemistry%tar_wat_indices(&
                        i))%solid_chemistry%mineral_zone)) &
                        then
                    !> We compute mass volumetric fractions of minerals from mass balance equation
                        call waters_new(this%chemistry%tar_wat_indices(&
                            i))%solid_chemistry%compute_mass_bal_mins(&
                            Delta_t)
                    !> We compute concentrations of minerals
                        call waters_new(this%chemistry%tar_wat_indices(&
                            i))%solid_chemistry%compute_conc_minerals_iter(Delta_t)
                    end if
                    !> We compute gas chemistry state variables
                    if (associated(waters_new(this%chemistry%tar_wat_indices(i))%gas_chemistry)) then
                    !> We compute concentrations of gases
                        call waters_new(this%chemistry%tar_wat_indices(&
                            i))%gas_chemistry%compute_conc_gases_iter(&
                            Delta_t,&
                            waters_new(i)%volume,[waters_new(i)%re_mean,waters_new(i)%rk])
                    !> We compute volume of gas
                        call waters_new(this%chemistry%tar_wat_indices(&
                            i))%gas_chemistry%compute_vol_gas_species_conc()
                    !> We compute activity coefficients of gases
                        call waters_new(this%chemistry%tar_wat_indices(&
                            i))%gas_chemistry%compute_log_act_coeffs_gases()
                    end if
                !> Accumulate aqueous reaction amounts
                    waters_new(this%chemistry%tar_wat_indices(i))%Rk_accum=&
                        waters_new(this%chemistry%tar_wat_indices(i))%Rk_accum+&
                        waters_new(this%chemistry%tar_wat_indices(i))%Rk
                !> Accumulate solid reaction amounts
                    waters_new(this%chemistry%tar_wat_indices(i))%solid_chemistry%Rk_accum=&
                        waters_new(this%chemistry%tar_wat_indices(i))%solid_chemistry%Rk_accum+&
                        waters_new(this%chemistry%tar_wat_indices(i))%solid_chemistry%Rk
                !> Deallocate
                    deallocate(c_hat)
                !end if
                !> We write output
                if ((k==this%chemistry%chem_out_options%time_steps(kk)) .and. (ind_non_can_vec(i)==&
                    this%chemistry%chem_out_options%ind_waters(ii))) then
                    !do j=1,n_v
                        write(unit,"(I10,*(ES15.5))") this%chemistry%chem_out_options%ind_waters(ii), (conc_nc(j), j=1,n_v)
                    !end do
                    if (ii<this%chemistry%chem_out_options%num_waters) then
                        ii=ii+1
                    else if (kk<this%chemistry%chem_out_options%num_time_steps) then
                        kk=kk+1
                        ii=1
                    else
                        exit
                    end if
                end if
                deallocate(conc_nc)
            end do
            !> We update target waters
            waters_old_old=waters_old
            waters_old=waters_new
        end do
!> We set the new target waters to the chemistry object
    this%chemistry%waters=waters_new
    close(unit)
end select
end subroutine