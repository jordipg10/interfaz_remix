!> Solves reactive mixing problem assuming ideal conditions
!> Computes concentrations, activities, activity coefficients, reaction rates and volumetric fractions of minerals (if present)
subroutine solve_reactive_mixing_cons(this,root,mixing_ratios_conc,mixing_ratios_Rk_init,mixing_waters_indices,&
        mixing_waters_indices_dom,time_discr_tpt,int_method_chem_reacts,mixing_ratios_Rk)
    use chemistry_m, only: chemistry_c
    use aqueous_chemistry_m, only: aqueous_chemistry_c, real_array_c, int_array_c, mixing_iter_comp_ideal, & 
        water_mixing_iter_EE_eq_kin_ideal, water_mixing_iter_EfI_eq_kin_anal_ideal, water_mixing_iter_EfI_eq_kin_anal_ideal_opt2, &
        mixing_iter_comp_exch_ideal, inf_norm_vec_real, water_mixing_iter_EE_kin_ideal, water_mixing_iter_EfI_kin_anal_ideal, &
        compute_c_tilde, compute_rk_tilde_expl, compute_rk_tilde_impl_opt1, water_mixing_iter_EI_eq_kin_anal_ideal_opt2, &
        water_mixing_iter_EI_kin_anal_ideal_opt2, compute_rk_tilde_impl_opt2, compute_rk_tilde_impl_opt3, compute_rk_tilde_impl_opt4
    use time_discr_m, only: time_discr_c
    implicit none
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: root !> root path
    !integer(kind=4), intent(in) :: unit !> file unit
    class(real_array_c), intent(in) :: mixing_ratios_conc !> mixing ratios matrix for concentrations
    class(real_array_c), intent(in) :: mixing_ratios_Rk_init !> initial mixing ratios matrix for reaction amounts
    class(int_array_c), intent(in) :: mixing_waters_indices !> matrix that contains indices of target waters that mix with each target water
    class(int_array_c), intent(in) :: mixing_waters_indices_dom !> matrix that contains indices of domain target waters that mix with each target water
    !real(kind=8), intent(in) :: F_mat(:) !> storage matrix (diagonal)
    class(time_discr_c), intent(in) :: time_discr_tpt !> time discretisation object (used to solve transport)
    integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method for chemical reactions
    class(real_array_c), intent(inout) :: mixing_ratios_Rk !> final mixing ratios matrix for reaction amounts
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
    integer(kind=4) :: n_nc !> number of variable activity species
    integer(kind=4) :: n_nc_aq !> number of aqueous variable activity species
    integer(kind=4) :: n_nc_aq_2 !> number of aqueous secondary variable activity species
    integer(kind=4) :: mix_ind !> mixing waters starting index
    integer(kind=4) :: num_can_vec !> number of canonical vectors
    integer(kind=4) :: num_non_can_vec !> number of canonical vectors
    integer(kind=4) :: unit !> unit for output file
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
    REAL(KIND=8), allocatable :: c_tilde(:) !> variable activity mobile species concentrations after mixing
    real(kind=8), allocatable :: rk_tilde_up(:) !> kinetic reaction rates after mixing
    real(kind=8), allocatable :: rk_tilde_down(:) !> kinetic reaction rates after mixing
    real(kind=8), allocatable :: rk_tilde(:) !> kinetic reaction rates
    REAL(KIND=8), allocatable :: conc_old(:,:) !> concentrations before mixing
    REAL(KIND=8), allocatable :: conc_nc(:) !> concentrations variable activity species
    !REAL(KIND=8), allocatable :: conc_comp(:) !> concentrations components
    type(aqueous_chemistry_c), allocatable :: target_waters_new(:) !> target waters time step k+1
    type(aqueous_chemistry_c), allocatable :: target_waters_old(:) !> target waters time step k
    type(aqueous_chemistry_c), allocatable :: target_waters_old_old(:) !> target waters time step k-1
    type(aqueous_chemistry_c), allocatable :: mix_waters(:) !> mixing waters
    integer(kind=4) :: num_lump !> number of lumpings
!> Procedure pointers
    !> reactive mixing subroutines
    procedure(water_mixing_iter_EI_kin_anal_ideal_opt2), pointer :: p_solver=>null()
    procedure(compute_rk_tilde_expl), pointer :: compute_rk_tilde=>null()
!> We initialise target waters
    target_waters_old=this%target_waters
    target_waters_old_old=target_waters_old
    target_waters_new=target_waters_old

    time=0d0
    kk=2 !> counter time steps chem_out_options
    ii=1 !> counter target waters chem_out_options
    this%chem_out_options%time_steps(this%chem_out_options%num_time_steps)=time_discr_tpt%Num_time !> chapuza
!> Chapuza
    ! if (this%chem_syst%num_kin_reacts==this%chem_syst%num_aq_kin_reacts) then
    !     allocate(rk_tilde(this%chem_syst%num_kin_reacts))
    !     allocate(rk(this%chem_syst%num_kin_reacts))
    ! else
    !     continue
    ! end if
    unit=7 !> arbitrary unit for output file
    open(unit,file=root//'.output',form="formatted",access="sequential",status="unknown")
!> We select reactive mixing subroutine depending on the nature of the chemical system and the methods to compute Jacobians and integrate in time
    if (this%chem_syst%num_kin_reacts>0 .and. this%chem_syst%num_eq_reacts>0) then !> equilibrium and kinetic reactions
        if (int_method_chem_reacts==1) then !> Euler explicit
            theta=0d0 !> Euler explicit
            !p_solver=>water_mixing_iter_EE_eq_kin_ideal
            compute_rk_tilde=>compute_rk_tilde_expl
            mix_ind=1 !> we mix all waters
            !mix_waters_indices=>mixing_waters_indices
        else if (this%Jac_opt==1) then !> analytical Jacobian
            mix_ind=2 !> we mix all waters except one
            if (int_method_chem_reacts==2) then !> Euler fully implicit
                theta=1d0 !> Euler fully implicit
            else if (int_method_chem_reacts==3) then !> Crank-Nicolson
                theta=5d-1 !> Crank-Nicolson parameter
            else
                error stop "Integration method for chemical reactions not implemented yet"
            end if
            if (this%rk_avg_opt==1) then
                !p_solver=>water_mixing_iter_EI_eq_kin_anal_ideal_opt1
            else if (this%rk_avg_opt==2) then
                p_solver=>water_mixing_iter_EI_eq_kin_anal_ideal_opt2
            else
                error stop "rk average option not implemented yet"
            end if
            if (this%cons_opt==1) then !> rk_tilde explicit
                !p_solver=>water_mixing_iter_CN_eq_kin_ideal_opt1
                compute_rk_tilde=>compute_rk_tilde_expl
            else if (this%cons_opt==2) then !> rk_tilde implicit
                if (this%rk_down_opt==1) then
                    compute_rk_tilde=>compute_rk_tilde_impl_opt1
                else if (this%rk_down_opt==2) then
                    compute_rk_tilde=>compute_rk_tilde_impl_opt2
                else if (this%rk_down_opt==3) then
                    compute_rk_tilde=>compute_rk_tilde_impl_opt3
                else
                    error stop "rk down option not implemented yet"
                end if
                !p_solver=>water_mixing_iter_CN_eq_kin_ideal_opt2
                !compute_rk_tilde=>compute_rk_tilde_impl_opt2
                !compute_rk_tilde=>compute_rk_tilde_impl_opt3
            else
                error stop "WMA consistent option not implemented yet"
            end if
        ! else if (int_method_chem_reacts==3 .and. this%Jac_opt==1) then !> Crank-Nicolson, analytical Jacobian
        !     theta=5d-1 !> Crank-Nicolson parameter
        !     mix_ind=2 !> we mix all waters except one
        !     if (this%rk_avg_opt==1) then
        !         !p_solver=>water_mixing_iter_EI_eq_kin_anal_ideal_opt1
        !     else if (this%rk_avg_opt==2) then
        !         p_solver=>water_mixing_iter_EI_eq_kin_anal_ideal_opt2
        !     else
        !         error stop "rk option not implemented yet"
        !     end if
        !     if (this%cons_opt==1) then !> rk_tilde explicit
        !         !p_solver=>water_mixing_iter_CN_eq_kin_ideal_opt1
        !         compute_rk_tilde=>compute_rk_tilde_expl
        !     else if (this%cons_opt==2) then !> rk_tilde implicit
        !         !p_solver=>water_mixing_iter_CN_eq_kin_ideal_opt2
        !         if (this%rk_down_opt==1) then
        !             compute_rk_tilde=>compute_rk_tilde_impl_opt1
        !         else if (this%rk_down_opt==2) then
        !             compute_rk_tilde=>compute_rk_tilde_impl_opt2
        !         else if (this%rk_down_opt==3) then
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
    else if (this%chem_syst%num_eq_reacts>0) then !> only equilibrium reactions
        p_solver=>mixing_iter_comp_ideal !> only equilibrium reactions
    else !> only kinetic reactions
        if (int_method_chem_reacts==1) then !> Euler explicit
            theta=0d0 !> Euler explicit
            !p_solver=>water_mixing_iter_EE_kin_ideal
            compute_rk_tilde=>compute_rk_tilde_expl
            mix_ind=1 !> we mix all waters
        else if (this%Jac_opt==1) then !> analytical Jacobian
            if (int_method_chem_reacts==2) then !> Euler fully implicit
                theta=1d0
            else if (int_method_chem_reacts==3) then !> Crank-Nicolson
                theta=5d-1 !> Crank-Nicolson parameter
            else
                error stop "Integration method for chemical reactions not implemented yet"
            end if
            if (this%rk_avg_opt==1) then
                !p_solver=>water_mixing_iter_EI_kin_anal_ideal_opt1
            else if (this%rk_avg_opt==2) then
                p_solver=>water_mixing_iter_EI_kin_anal_ideal_opt2
            else
                error stop "rk average option not implemented yet"
            end if
            mix_ind=2 !> we mix all waters except one
            if (this%cons_opt==1) then !> mix kinetic reaction rates from previous time step
                !p_solver=>water_mixing_iter_EfI_eq_kin_anal_ideal_opt1
                compute_rk_tilde=>compute_rk_tilde_expl
            else if (this%cons_opt==2) then !> solve upstream to downstream
                !p_solver=>water_mixing_iter_EI_eq_kin_anal_ideal_opt2
                if (this%rk_down_opt==1) then
                    compute_rk_tilde=>compute_rk_tilde_impl_opt1
                else if (this%rk_down_opt==2) then
                   compute_rk_tilde=>compute_rk_tilde_impl_opt2
                else if (this%rk_down_opt==3) then
                   compute_rk_tilde=>compute_rk_tilde_impl_opt3
                else if (this%rk_down_opt==4) then
                    !compute_rk_tilde=>compute_rk_tilde_impl_opt4
                else
                    error stop "rk down option not implemented yet"
                end if
            else
                error stop "WMA consistent option not implemented yet"
            end if
        ! else if (int_method_chem_reacts==3 .and. this%Jac_opt==1) then !> Crank-Nicolson, analytical Jacobian
        !     theta=5d-1 !> Crank-Nicolson parameter
        !     mix_ind=2 !> we mix all waters except one
        !     if (this%rk_avg_opt==1) then
        !         !p_solver=>water_mixing_iter_EI_eq_kin_anal_ideal_opt1
        !     else if (this%rk_avg_opt==2) then
        !         p_solver=>water_mixing_iter_EI_kin_anal_ideal_opt2
        !     else
        !         error stop "rk option not implemented yet"
        !     end if
        !     if (this%cons_opt==1) then !> mix kinetic reaction rates from previous time step
        !         !p_solver=>water_mixing_iter_CN_eq_kin_ideal_opt1
        !         !compute_rk_tilde=>compute_rk_tilde_expl
        !     else if (this%cons_opt==2) then !> solve upstream to downstream
        !         !p_solver=>water_mixing_iter_CN_eq_kin_ideal_opt2
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
    call mixing_ratios_conc%get_can_vec(this%CV_params%abs_tol,num_can_vec,ind_can_vec,num_non_can_vec,ind_non_can_vec)
    !> We set old concentrations in Dirichlet boundary target waters (chapuza)
    do i=1,this%num_target_waters
        call target_waters_new(i)%set_conc_old() !> 
        call target_waters_new(i)%set_conc_old_old() !> 
        call target_waters_new(i)%solid_chemistry%set_conc_old() !> 
        call target_waters_new(i)%solid_chemistry%set_conc_old_old() !> 
        !call initialise_iterative_method()
        !call target_waters_new(ind_non_can_vec(i))%update_rk_old()
        !call target_waters_new(ind_can_vec(i))%solid_chemistry%update_rk_old()
    end do
    !do i=1,num_non_can_vec
    !    call target_waters_new(ind_can_vec(i))%set_conc_old() !>
    !    call target_waters_new(ind_can_vec(i))%solid_chemistry%set_conc_old() !>
    !    call target_waters_new(ind_can_vec(i))%solid_chemistry%set_conc_old_old() !>
    !end do
    !> Time loop
    !cntr_can_vec=1 !> counter canonical vectors
        do k=1,time_discr_tpt%Num_time
            Delta_t=time_discr_tpt%get_Delta_t(k)
            time=time+Delta_t
            if (k==this%chem_out_options%time_steps(kk)) then
               write(unit,"(2x,'t = ',*(ES15.5))") time
               write(unit,"(20x,*(A15))") (this%chem_syst%aq_phase%aq_species(this%chem_out_options%ind_aq_species(j))%name, &
                j=1,this%chem_out_options%num_aq_species)
            end if
            !> We update concentrations and kinetic reaction rates in domain target waters (chapuza)
            do i=1,num_non_can_vec
                call target_waters_new(ind_non_can_vec(i))%update_conc_old() !> 
                call target_waters_new(ind_non_can_vec(i))%solid_chemistry%update_conc_old() !> 
                !call initialise_iterative_method()
                call target_waters_new(ind_non_can_vec(i))%update_rk_old()
                call target_waters_new(ind_non_can_vec(i))%solid_chemistry%update_rk_old()
            end do
            !> Domain target waters loop
            do i=1,num_non_can_vec
                n_p=this%target_waters(&
                    this%dom_tar_wat_indices(ind_non_can_vec(i)))%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
                n_nc=this%target_waters(&
                    this%dom_tar_wat_indices(ind_non_can_vec(i)))%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
                n_nc_aq=this%target_waters(&
                   this%dom_tar_wat_indices(ind_non_can_vec(i)))%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
                n_nc_aq_2=this%target_waters(&
                    this%dom_tar_wat_indices(&
                    ind_non_can_vec(i)))%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                    if (this%target_waters(&
                        this%dom_tar_wat_indices(ind_non_can_vec(i)))%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl &
                        >0) then !> variable activity species are aqueous and solid
                        !p_solver=>mixing_iter_comp_exch_ideal !> only equilibrium reactions
                    end if
                    allocate(conc_nc(n_nc))
                    allocate(rk_tilde(n_nc))
                    if (mixing_waters_indices%cols(ind_non_can_vec(i))%col_1(1)/=this%dom_tar_wat_indices(ind_non_can_vec(i))) then
                        print *, mixing_waters_indices%cols(ind_non_can_vec(i))%col_1(1),&
                            this%dom_tar_wat_indices(ind_non_can_vec(i))
                        error stop "Target waters not in the right order"
                    end if
                    !conc_old(:,1)=target_waters_old(ind_non_can_vec(i))%get_conc_nc()
                    !do j=2,mixing_waters_indices%cols(i)%dim
                    !    !print *, mixing_waters_indices%cols(i)%col_1(j)
                    !    conc_old(:,j)=target_waters_old(mixing_waters_indices%cols(i)%col_1(j))%get_conc_nc()
                    !    !allocate(perm(n_nc_aq))
                    !    !perm = [target_waters_old(ind_non_can_vec(i))%ind_prim_species, &
                    !    !        target_waters_old(ind_non_can_vec(i))%ind_sec_species(1:n_nc_aq_2)]
                    !    !print *, ind_non_can_vec(i)
                    !    conc_old(1:n_nc_aq,j) = target_waters_old(mixing_waters_indices%cols(i)%col_1(j))%get_conc_nc_aq_perm(&
                    !        target_waters_old(ind_non_can_vec(i))%ind_var_act_species)
                    !    !deallocate(perm)
                    !    ! do l=1,n_nc_aq
                    !    !     conc_old(l,j)=target_waters_old(mixing_waters_indices%cols(i)%col_1(j))%concentrations(&
                    !    !         target_waters_old(ind_non_can_vec(i))%indices_aq_species(l))
                    !    ! end do
                    !end do
                    !> We solve mixing caused by transport
                    c_tilde=target_waters_old(this%dom_tar_wat_indices(ind_non_can_vec(i)))%get_conc_nc()
                    !c_tilde(1:n_nc_aq)=compute_c_tilde(target_waters_old(this%dom_tar_wat_indices(ind_non_can_vec(i)),mixing_ratios_conc%cols(i)%col_1,&
                    !    conc_old(1:n_nc_aq,:))
                    !print *, target_waters_old(this%dom_tar_wat_indices(ind_non_can_vec(i))%ind_var_act_species
                    !print *, target_waters_old(this%dom_tar_wat_indices(ind_non_can_vec(i))%indices_aq_species(&
                    !    this%chem_syst%aq_phase%ind_diss_solids)
                    !call target_waters_old(this%dom_tar_wat_indices(ind_non_can_vec(i))%compute_rk(rk)
                    !call target_waters_old(this%dom_tar_wat_indices(ind_non_can_vec(i))%update_rk_old()
                    !call target_waters_old(this%dom_tar_wat_indices(ind_non_can_vec(i))%solid_chemistry%update_rk_old()
                    call compute_c_tilde(target_waters_old(mixing_waters_indices%cols(ind_non_can_vec(i))%col_1(&
                        1:mixing_waters_indices%cols(ind_non_can_vec(i))%dim-2)),&
                        mixing_ratios_conc%cols(ind_non_can_vec(i))%col_1,c_tilde(1:n_nc_aq)) !> we assume mixing ratios are the same for conc and rk (chapuza)
                    !> chapuza
                    !allocate(mix_waters(this%num_target_waters_dom-num_can_vec-mix_ind+1))
                    !do j=1,this%num_target_waters_dom-mix_ind+1
                    !    if (j==
                    !    mix_waters(j)=target_waters_new(mixing_waters_indices%cols(i)%col_1(&
                    !    mix_ind+j-1))
                    if (target_waters_new(this%dom_tar_wat_indices(ind_non_can_vec(i)))%indices_rk%num_cols>0) then
                        mix_waters=target_waters_new(mixing_waters_indices_dom%cols(ind_non_can_vec(i))%col_1(&
                            mix_ind:mixing_waters_indices_dom%cols(ind_non_can_vec(i))%dim-2))
                        call compute_rk_tilde(mix_waters,mixing_ratios_Rk_init%cols(ind_non_can_vec(i))%col_1(mix_ind:),&
                            mixing_waters_indices_dom%cols(ind_non_can_vec(i))%col_1(mixing_waters_indices_dom%cols(&
                            ind_non_can_vec(i))%dim-1),&
                            mixing_waters_indices_dom%cols(ind_non_can_vec(i))%col_1(mixing_waters_indices_dom%cols(&
                            ind_non_can_vec(i))%dim),theta,Delta_t,rk_tilde)
                        !> chapuza
                        target_waters_new(mixing_waters_indices_dom%cols(ind_non_can_vec(i))%col_1(mix_ind:&
                            mixing_waters_indices_dom%cols(ind_non_can_vec(i))%dim-2))=mix_waters
                        call target_waters_new(this%dom_tar_wat_indices(ind_non_can_vec(i)))%solid_chemistry%modify_mix_ratios_rk(&
                            mixing_ratios_Rk_init%cols(ind_non_can_vec(i))%col_1(1),c_tilde,Delta_t,rk_tilde,&
                            mixing_ratios_Rk%cols(ind_non_can_vec(i))%col_1(1),num_lump)
                        !> chapuza
                        this%num_lump=this%num_lump+num_lump !> we update number of lumpings
                    end if
                    !y=sum(mixing_ratios_Rk%cols(i)%col_1(2:1+mixing_waters_indices%cols(i)%col_1(mixing_waters_indices%cols(i)%dim-1)))
                    !print *, y*Delta_t
                    !do j=1,n_nc_aq
                    !    if (target_waters_old(ind_non_can_vec(i))%ind_var_act_species(j) /= &
                    !        target_waters_old(ind_non_can_vec(i))%indices_aq_species(&
                    !        this%chem_syst%aq_phase%ind_diss_solids(j))) then
                    !        c_tilde(target_waters_old(ind_non_can_vec(i))%ind_var_act_species(j))=c_tilde_aux(&
                    !            target_waters_old(ind_non_can_vec(i))%indices_aq_species(&
                    !            this%chem_syst%aq_phase%ind_diss_solids(j)))
                    !    end if
                    !end do
                    !> We solve reactive mixing iteration
                     call p_solver(target_waters_new(this%dom_tar_wat_indices(ind_non_can_vec(i))),target_waters_old_old(&
                        this%dom_tar_wat_indices(ind_non_can_vec(i)))%get_c1(),c_tilde,rk_tilde,mixing_ratios_Rk%cols(&
                        ind_non_can_vec(i))%col_1(1),Delta_t,theta,conc_nc)
                    !> We compute equilibrium reaction rates from mass balance equation
                    if (this%target_waters(this%dom_tar_wat_indices(ind_non_can_vec(i))&
                        )%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions>0 .and. &
                        this%target_waters(this%dom_tar_wat_indices(ind_non_can_vec(i)))%indices_rk%num_cols>0) then
                        call target_waters_new(this%dom_tar_wat_indices(ind_non_can_vec(i)))%compute_Re_mean_rk(&
                            c_tilde(n_p+1:n_p+n_nc_aq_2),Delta_t,&
                            theta,rk_tilde(n_p+1:n_p+n_nc_aq_2)) !> chapuza
                    else if (this%target_waters(this%dom_tar_wat_indices(ind_non_can_vec(i))&
                        )%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions>0) then
                        call target_waters_new(this%dom_tar_wat_indices(ind_non_can_vec(i)))%compute_Re_mean(&
                            c_tilde(n_p+1:n_p+n_nc_aq_2),Delta_t,&
                            theta) !> chapuza
                    end if
                    ! call target_waters_new(ind_non_can_vec(i))%compute_r_eq(c_tilde(n_p+1:n_p+n_nc_aq_2),Delta_t,&
                    !     rk_tilde)
                    !> We compute solid chemistry state variables
                    if (associated(target_waters_new(this%dom_tar_wat_indices(ind_non_can_vec(i)))%solid_chemistry%mineral_zone)) &
                        then
                    !> We compute mass volumetric fractions of minerals from mass balance equation
                        call target_waters_new(this%dom_tar_wat_indices(ind_non_can_vec(i)))%solid_chemistry%compute_mass_bal_mins(&
                            Delta_t)
                    !> We compute concentrations of minerals
                        call target_waters_new(this%dom_tar_wat_indices(&
                            ind_non_can_vec(i)))%solid_chemistry%compute_conc_minerals_iter_EE(Delta_t)
                    end if
                    !> We compute gas chemistry state variables
                    if (associated(target_waters_new(this%dom_tar_wat_indices(ind_non_can_vec(i)))%gas_chemistry)) then
                    !> We compute concentrations of gases
                        call target_waters_new(this%dom_tar_wat_indices(ind_non_can_vec(i)))%gas_chemistry%compute_conc_gases_iter(&
                            Delta_t,&
                            target_waters_new(i)%volume,[target_waters_new(i)%r_eq,target_waters_new(i)%rk])
                    !> We compute volume of gas
                        call target_waters_new(this%dom_tar_wat_indices(ind_non_can_vec(i)))%gas_chemistry%compute_vol_gas_conc()
                    !> We compute activity coefficients of gases
                        call target_waters_new(this%dom_tar_wat_indices(&
                            ind_non_can_vec(i)))%gas_chemistry%compute_log_act_coeffs_gases()
                    end if
                !> Deallocate
                    deallocate(c_tilde,rk_tilde)
                !end if
                !> We write output
                if ((k==this%chem_out_options%time_steps(kk)) .and. (ind_non_can_vec(i)==&
                    this%chem_out_options%ind_target_waters(ii))) then
                    !do j=1,n_nc
                        write(unit,"(I10,*(ES15.5))") this%chem_out_options%ind_target_waters(ii), (conc_nc(j), j=1,n_nc)
                    !end do
                    if (ii<this%chem_out_options%num_target_waters) then
                        ii=ii+1
                    else if (kk<this%chem_out_options%num_time_steps) then
                        kk=kk+1
                        ii=1
                    else
                        exit
                    end if
                end if
                deallocate(conc_nc)
            end do
        !> We update target waters
            target_waters_old_old=target_waters_old
            target_waters_old=target_waters_new
        !> Chapuza
            !compute_rk_tilde=>compute_rk_tilde_impl_opt2
        end do
!> We set the new target waters to the chemistry object
    this%target_waters=target_waters_new
    close(unit)
    !write(*,*) "Number of lumpings: ", this%num_lump
 end subroutine