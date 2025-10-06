!> Solves reactive mixing problem using lumping in the mixing ratios of reaction rates
!> Computes concentrations, activities, activity coefficients, reaction rates and volumetric fractions of minerals (if present)
subroutine solve_reactive_mixing_lump(this,root,mixing_ratios,mixing_waters_indices,time_discr_tpt,&
    int_method_chem_reacts)
    use chemistry_m, only: chemistry_c
    use aqueous_chemistry_m, only: aqueous_chemistry_c, real_array_c, int_array_c, mixing_iter_comp, & 
        water_mixing_iter_EE_eq_kin, water_mixing_iter_EfI_eq_kin_anal_lump, mixing_iter_comp_exch, &
        compute_c_tilde, inf_norm_vec_real, mixing_iter_comp_lump, mixing_iter_comp_exch_lump
    use time_discr_m, only: time_discr_c
    implicit none
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: root
    !integer(kind=4), intent(in) :: unit
    class(real_array_c), intent(in) :: mixing_ratios !> mixing ratios matrix
    class(int_array_c), intent(in) :: mixing_waters_indices !> matrix that contains indices of target waters that mix with each target water
    !real(kind=8), intent(in) :: F_mat(:) !> storage matrix (diagonal)
    class(time_discr_c), intent(in) :: time_discr_tpt !> time discretisation object (used to solve transport)
    integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method for chemical reactions
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
    integer(kind=4) :: unit !> unit of output file
    integer(kind=4), allocatable :: tar_gas_indices(:) !> indices target gases in each reactive zone
    integer(kind=4), allocatable :: tar_sol_indices(:) !> indices target solids in each reactive zone
    integer(kind=4), allocatable :: tar_wat_indices(:) !> indices target waters in each reactive zone
    integer(kind=4), allocatable :: perm(:) !> permutation vector aqueous concentrations
    REAL(KIND=8) :: time !> time
    REAL(KIND=8) :: Delta_t !> time step
    real(kind=8) :: theta !> time weighting factor
    REAL(KIND=8), allocatable :: c_tilde(:) !> variable activity mobile species concentrations after mixing
    real(kind=8), allocatable :: c_tilde_aux(:) !> chapuza
    REAL(KIND=8), allocatable :: conc_old(:,:) !> concentrations before mixing
    REAL(KIND=8), allocatable :: conc_nc(:) !> concentrations variable activity species
    !REAL(KIND=8), allocatable :: conc_comp(:) !> concentrations components
    type(aqueous_chemistry_c), allocatable :: target_waters_new(:) !> target waters time step k+1
    type(aqueous_chemistry_c), allocatable :: target_waters_old(:) !> target waters time step k
    type(aqueous_chemistry_c), allocatable :: target_waters_old_old(:) !> target waters time step k-1
    !type(reactive_zone_c) :: react_zone
    !type(aq_phase_c), target :: aq_phase
!> Procedure pointers
    !> reactive mixing subroutines
    procedure(mixing_iter_comp_lump), pointer :: p_solver=>null()
    ! procedure(mixing_iter_comp_ideal), pointer :: p_solver_bis=>null()
!> We initialise target waters
    target_waters_old=this%target_waters
    target_waters_old_old=target_waters_old
    target_waters_new=target_waters_old

    time=0d0
    kk=2 !> counter time steps chem_out_options
    ii=1 !> counter target waters chem_out_options
    this%chem_out_options%time_steps(this%chem_out_options%num_time_steps)=time_discr_tpt%Num_time !> chapuza
            
    unit=7 !> unit of output file
    open(unit,file=root//'.out',form="formatted",access="sequential",status="unknown")
!> We select reactive mixing subroutine depending on the nature of the chemical system and the methods to compute Jacobians and integrate in time
    if (this%chem_syst%num_kin_reacts>0 .and. this%chem_syst%num_eq_reacts>0) then !> equilibrium and kinetic reactions
        if (int_method_chem_reacts==1) then !> Euler explicit
            theta=0d0 !> Euler explicit
            !p_solver=>water_mixing_iter_EE_eq_kin_ideal
            !compute_rk_tilde=>compute_rk_tilde_expl
            !mix_ind=1 !> we mix all waters
            !mix_waters_indices=>mixing_waters_indices
        else if (this%Jac_opt==1) then !> analytical Jacobian
            !mix_ind=2 !> we mix all waters except one
            if (int_method_chem_reacts==2) then !> Euler fully implicit
                theta=1d0 !> Euler fully implicit
            else if (int_method_chem_reacts==3) then !> Crank-Nicolson
                theta=5d-1 !> Crank-Nicolson parameter
            else
                error stop "Integration method for chemical reactions not implemented yet"
            end if
            if (this%rk_avg_opt==1) then
                !p_solver=>water_mixing_iter_EI_eq_kin_anal_opt1
            else if (this%rk_avg_opt==2) then
                !p_solver=>water_mixing_iter_EI_eq_kin_anal_opt2_lump
            else
                error stop "rk average option not implemented yet"
            end if
            ! if (this%cons_opt==1) then !> rk_tilde explicit
            !     !p_solver=>water_mixing_iter_CN_eq_kin_opt1
            !     compute_rk_tilde=>compute_rk_tilde_expl
            ! else if (this%cons_opt==2) then !> rk_tilde implicit
            !     if (this%rk_down_opt==1) then
            !         compute_rk_tilde=>compute_rk_tilde_impl_opt1
            !     else if (this%rk_down_opt==2) then
            !         compute_rk_tilde=>compute_rk_tilde_impl_opt2
            !     else if (this%rk_down_opt==3) then
            !         compute_rk_tilde=>compute_rk_tilde_impl_opt3
            !     else
            !         error stop "rk down option not implemented yet"
            !     end if
            !     !p_solver=>water_mixing_iter_CN_eq_kin_ideal_opt2
            !     !compute_rk_tilde=>compute_rk_tilde_impl_opt2
            !     !compute_rk_tilde=>compute_rk_tilde_impl_opt3
            ! else
            !     error stop "WMA consistent option not implemented yet"
            ! end if
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
        p_solver=>mixing_iter_comp_lump !> only equilibrium reactions
    else !> only kinetic reactions
        if (int_method_chem_reacts==1) then !> Euler explicit
            theta=0d0 !> Euler explicit
            !p_solver=>water_mixing_iter_EE_kin_ideal
            !compute_rk_tilde=>compute_rk_tilde_expl
            !mix_ind=1 !> we mix all waters
        else if (this%Jac_opt==1) then !> analytical Jacobian
            if (int_method_chem_reacts==2) then !> Euler fully implicit
                theta=1d0
            else if (int_method_chem_reacts==3) then !> Crank-Nicolson
                theta=5d-1 !> Crank-Nicolson parameter
            else
                error stop "Integration method for chemical reactions not implemented yet"
            end if
            !mix_ind=2 !> we mix all waters except one
            ! if (this%cons_opt==1) then !> mix kinetic reaction rates from previous time step
            !     !p_solver=>water_mixing_iter_EfI_eq_kin_anal_ideal_opt1
            !     compute_rk_tilde=>compute_rk_tilde_expl
            ! else if (this%cons_opt==2) then !> solve upstream to downstream
            !     !p_solver=>water_mixing_iter_EfI_eq_kin_anal_ideal_opt2
            !     if (this%rk_down_opt==1) then
            !         compute_rk_tilde=>compute_rk_tilde_impl_opt1
            !     else if (this%rk_down_opt==2) then
            !         compute_rk_tilde=>compute_rk_tilde_impl_opt2
            !     else if (this%rk_down_opt==3) then
            !         compute_rk_tilde=>compute_rk_tilde_impl_opt3
            !     else
            !         error stop "rk down option not implemented yet"
            !     end if
            ! else
            !     error stop "WMA consistent option not implemented yet"
            ! end if
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
!> Time loop
        do k=1,time_discr_tpt%Num_time
            Delta_t=time_discr_tpt%get_Delta_t(k)
            time=time+Delta_t
            if (k==this%chem_out_options%time_steps(kk)) then
                write(unit,"(2x,'t = ',*(ES15.5))") time
                write(unit,"(20x,*(A15))") (this%chem_syst%aq_phase%aq_species(this%chem_out_options%ind_aq_species(j))%name, &
                j=1,this%chem_out_options%num_aq_species)
            end if
        !> Target waters loop
            do i=1,this%num_target_waters_dom
                n_p=this%target_waters(this%dom_tar_wat_indices(i))%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
                n_nc=this%target_waters(this%dom_tar_wat_indices(i)&
                    )%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
                n_nc_aq=this%target_waters(this%dom_tar_wat_indices(i)&
                    )%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
                n_nc_aq_2=this%target_waters(this%dom_tar_wat_indices(i)&
                    )%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                if (abs(mixing_ratios%cols(i)%col_1(1)-1d0)<&
                this%target_waters(this%dom_tar_wat_indices(i))%solid_chemistry%reactive_zone%CV_params%abs_tol .and. &
                    inf_norm_vec_real(mixing_ratios%cols(i)%col_1(2:mixing_ratios%cols(i)%dim))<&
                this%target_waters(this%dom_tar_wat_indices(i))%solid_chemistry%reactive_zone%CV_params%abs_tol) then
                    continue
                else
                    if (this%target_waters(this%dom_tar_wat_indices(i))%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl&
                        >0) then !> variable activity species are aqueous and solid
                        p_solver=>mixing_iter_comp_exch_lump !> only equilibrium reactions
                    end if
                    allocate(conc_nc(n_nc))
                    !allocate(conc_comp(n_p))
                    !allocate(conc_old(n_nc,mixing_ratios%cols(i)%dim)) !> chapuza
                    if (mixing_waters_indices%cols(i)%col_1(1)/=this%dom_tar_wat_indices(i)) then
                        print *, mixing_waters_indices%cols(i)%col_1(1), this%dom_tar_wat_indices(i)
                        error stop "Target waters not in the right order"
                    end if
                    !conc_old(:,1)=target_waters_old(this%dom_tar_wat_indices(i))%get_conc_nc()
                    !do j=2,mixing_waters_indices%cols(i)%dim
                    !    !print *, mixing_waters_indices%cols(i)%col_1(j)
                    !    conc_old(:,j)=target_waters_old(mixing_waters_indices%cols(i)%col_1(j))%get_conc_nc()
                    !    !allocate(perm(n_nc_aq))
                    !    !perm = [target_waters_old(this%dom_tar_wat_indices(i))%ind_prim_species, &
                    !    !        target_waters_old(this%dom_tar_wat_indices(i))%ind_sec_species(1:n_nc_aq_2)]
                    !    !print *, this%dom_tar_wat_indices(i)
                    !    conc_old(1:n_nc_aq,j) = target_waters_old(mixing_waters_indices%cols(i)%col_1(j))%get_conc_nc_aq_perm(&
                    !        target_waters_old(this%dom_tar_wat_indices(i))%ind_var_act_species)
                    !    !deallocate(perm)
                    !    ! do l=1,n_nc_aq
                    !    !     conc_old(l,j)=target_waters_old(mixing_waters_indices%cols(i)%col_1(j))%concentrations(&
                    !    !         target_waters_old(this%dom_tar_wat_indices(i))%indices_aq_species(l))
                    !    ! end do
                    !end do
                !> We solve mixing caused by transport
                    c_tilde=target_waters_old(this%dom_tar_wat_indices(i))%get_conc_nc()
                    !c_tilde(1:n_nc_aq)=compute_c_tilde(target_waters_old(this%dom_tar_wat_indices(i)),mixing_ratios%cols(i)%col_1,&
                    !    conc_old(1:n_nc_aq,:))
                    !print *, target_waters_old(this%dom_tar_wat_indices(i))%ind_var_act_species
                    !print *, target_waters_old(this%dom_tar_wat_indices(i))%indices_aq_species(&
                    !    this%chem_syst%aq_phase%ind_diss_solids)
                    !c_tilde(1:n_nc_aq)=compute_c_tilde_bis(target_waters_old(mixing_waters_indices%cols(i)%col_1),&
                    !    mixing_ratios%cols(i)%col_1)
                    !c_tilde_aux=c_tilde(1:n_nc_aq) !> chapuza
                    !do j=1,n_nc_aq
                    !    if (target_waters_old(this%dom_tar_wat_indices(i))%ind_var_act_species(j) /= &
                    !        target_waters_old(this%dom_tar_wat_indices(i))%indices_aq_species(&
                    !        this%chem_syst%aq_phase%ind_diss_solids(j))) then
                    !        c_tilde(target_waters_old(this%dom_tar_wat_indices(i))%ind_var_act_species(j))=c_tilde_aux(&
                    !            target_waters_old(this%dom_tar_wat_indices(i))%indices_aq_species(&
                    !            this%chem_syst%aq_phase%ind_diss_solids(j)))
                    !    end if
                    !end do
                !> We solve reactive mixing iteration
                    call p_solver(target_waters_new(this%dom_tar_wat_indices(i)),target_waters_old_old(&
                        this%dom_tar_wat_indices(i))%get_c1(),target_waters_old(&
                        this%dom_tar_wat_indices(i))%get_c2nc(),c_tilde,Delta_t,conc_nc)
                !> We compute equilibrium reaction rates from mass balance equation
                    call target_waters_new(this%dom_tar_wat_indices(i))%compute_Re_mean_lump(c_tilde(n_p+1:n_p+n_nc_aq_2),Delta_t,&
                        theta)
            !> We compute solid chemistry state variables
                    if (associated(target_waters_new(this%dom_tar_wat_indices(i))%solid_chemistry%mineral_zone)) then
                    !> We compute mass volumetric fractions of minerals from mass balance equation
                        call target_waters_new(this%dom_tar_wat_indices(i))%solid_chemistry%compute_mass_bal_mins(Delta_t)
                    !> We compute concentrations of minerals
                        call target_waters_new(this%dom_tar_wat_indices(i))%solid_chemistry%compute_conc_minerals_iter_EE(Delta_t)
                    end if
                !> We compute gas chemistry state variables
                    if (associated(target_waters_new(this%dom_tar_wat_indices(i))%gas_chemistry)) then
                    !> We compute concentrations of gases
                        call target_waters_new(this%dom_tar_wat_indices(i))%gas_chemistry%compute_conc_gases_iter(Delta_t,&
                            target_waters_new(i)%volume,[target_waters_new(i)%r_eq,target_waters_new(i)%rk])
                    !> We compute volume of gas   
                        call target_waters_new(this%dom_tar_wat_indices(i))%gas_chemistry%compute_vol_gas_conc()
                    !> We compute activity coefficients of gases    
                        call target_waters_new(this%dom_tar_wat_indices(i))%gas_chemistry%compute_log_act_coeffs_gases()
                    end if
                !> We write output
                    if (k==this%chem_out_options%time_steps(kk) .and. i==this%chem_out_options%ind_target_waters(ii)) then
                        !do j=1,this%target_waters(i)%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
                            write(unit,"(I10,*(ES15.5))") i, (conc_nc(j), j=1,this%chem_out_options%num_aq_species)
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
                !> Deallocate
                    deallocate(c_tilde,conc_nc)
                end if
            end do
        !> We update target waters
            target_waters_old_old=target_waters_old
            target_waters_old=target_waters_new
        end do
!> We set the new target waters to the chemistry object
    this%target_waters=target_waters_new
    close(unit)
end subroutine