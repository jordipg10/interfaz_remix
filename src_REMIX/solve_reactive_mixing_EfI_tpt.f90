!!> Solves reactive mixing for Euler implicit in transport
!subroutine solve_reactive_mixing_EI_tpt(this,mixing_ratios,mixing_waters_indices,F_mat,time_discr_tpt,int_method_chem_reacts)
!    use chemistry_m
!    use BCs_m
!    use metodos_sist_lin_m
!    implicit none
!    class(chemistry_c) :: this
!    class(real_array_c), intent(in) :: mixing_ratios !> mixing ratios
!    class(int_array_c), intent(in) :: mixing_waters_indices !> matrix that contains indices of target waters that mix with each target water
!    class(diag_matrix_c), intent(in) :: F_mat !> storage matrix
!    class(time_discr_c), intent(in) :: time_discr_tpt !> time discretisation transport
!    integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method chemical reactions
!    
!    integer(kind=4) :: i,j,n,l,k,niter,CV_flag,num_tar_wat,num_tar_sol
!    integer(kind=4), allocatable :: tar_sol_indices(:),tar_wat_indices(:),mix_wat_indices(:)
!    
!    real(kind=8), parameter :: mu=1d0, tolerance=1d-16
!    real(kind=8), allocatable :: mixing_ratios_vec(:),r_eq_tilde(:,:),c1_old_old(:),c1_old(:),A_mat_inv(:,:),r_eq_mat(:,:),rk_mat(:,:)
!    type(aqueous_chemistry_c), allocatable :: target_waters_old(:),target_waters_new(:),target_waters_old_old(:)
!
!    !procedure(water_mixing_iter_EI_tpt_EE_rk_eq_kin), pointer :: p_solver=>null()
!    procedure(transport_iter_comp_EE_aq_chem), pointer :: p_solver=>null()
!    
!    target_waters_old=this%target_waters_init
!    target_waters_old_old=target_waters_old
!    target_waters_new=target_waters_old
!    
!    !if (time_discr_tpt%int_method/=2) then
!    !    error stop
!    
!    !if (size(mixing_ratios,1)/=size(mixing_ratios,2)) then
!    !    error stop
!    !end if
!    
!!> Chapuza
!    !allocate(mix_wat_indices(n-1,n))
!    !do j=1,n !> por columnas
!    !    do i=1,j-1
!    !        mix_wat_indices(i,j)=i
!    !    end do
!    !    do i=j+1,n
!    !        mix_wat_indices(i-1,j)=i
!    !    end do
!    !end do
!    
!    n=mixing_ratios%num_cols
!    allocate(rk_mat(this%chem_syst%num_kin_reacts,n))
!    rk_mat=0d0 !> autentica chapuza
!    if (this%num_reactive_zones>0) then
!        if (this%chem_syst%num_kin_reacts>0) then !> eq and kin reactions
!            if (int_method_chem_reacts==1) then
!                !p_solver=>water_mixing_iter_EI_tpt_EE_rk_eq_kin
!                p_solver=>water_mixing_iter_EE_eq_kin
!            else if (int_method_chem_reacts==2) then
!                !p_solver=>water_mixing_iter_EI_tpt_EfI_rk_eq_kin
!                p_solver=>water_mixing_iter_EfI_eq_kin
!            else
!                error stop "Integration method for chemical reactions not implemented yet"
!            end if
!        else
!            !p_solver=>transport_iter_EI_aq_chem !> only eq reactions
!            p_solver=>transport_iter_comp_EE_aq_chem !> only eq reactions
!        end if
!        do l=1,this%num_reactive_zones
!            call this%link_target_solids_reactive_zone(l,tar_sol_indices)
!            call this%link_target_waters_target_solids(tar_sol_indices,tar_wat_indices)
!            num_tar_wat=size(tar_wat_indices)
!            num_tar_sol=size(tar_sol_indices) !> we assume = num_tar_wat
!            allocate(r_eq_tilde(this%reactive_zones(l)%num_eq_reactions,num_tar_sol))
!            allocate(r_eq_mat(this%reactive_zones(l)%num_eq_reactions,num_tar_sol))
!            
!        !> Time loop
!            do k=1,time_discr_tpt%Num_time
!            !!> We compute rk_mat
!            !    do i=1,num_tar_wat
!            !        call target_waters_new(tar_wat_indices(i))%compute_rk_aq_chem()
!            !        rk_mat(:,i)=target_waters_new(tar_wat_indices(i))%rk
!            !    end do
!            !> Target waters loop
!                do i=1,num_tar_wat
!                !> Esto debe ir fuera del time loop
!                    mix_wat_indices=mixing_waters_indices%cols(tar_wat_indices(i))%col_1
!                    mixing_ratios_vec=mixing_ratios%cols(tar_wat_indices(i))%col_1
!                !> We initialise iterative method
!                    call initialise_iterative_method(target_waters_old_old(tar_wat_indices(i))%concentrations(1:target_waters_old_old(tar_wat_indices(i))%speciation_alg%num_prim_species),target_waters_old(tar_wat_indices(i))%concentrations(1:target_waters_old(tar_wat_indices(i))%speciation_alg%num_prim_species),mu,target_waters_new(tar_wat_indices(i))%concentrations(1:target_waters_new(tar_wat_indices(i))%speciation_alg%num_prim_species))
!                !> We solve reactive mixing iteration and speciate
!                    !call p_solver(target_waters_new(tar_wat_indices(i)),tar_wat_indices(i),mixing_ratios(:,i:i),target_waters_old(mix_wat_indices(:,i)),rk_mat,F_mat%diag(tar_sol_indices(i)),time_discr_tpt%get_Delta_t(k))
!                    call p_solver(target_waters_new(tar_wat_indices(i)),mixing_ratios_vec,target_waters_old(mix_wat_indices),F_mat%diag(tar_sol_indices(i)),time_discr_tpt%get_Delta_t(k))
!                    !do j=1,n-1
!                    !    call target_waters_old(mix_wat_indices(j,i))%compute_molarities() !> in order to compute reaction rates
!                    !end do
!                    !call target_waters_old(tar_wat_indices(i))%compute_molarities() !> in order to compute reaction rates
!                    !call target_waters_new(tar_wat_indices(i))%compute_molarities() !> in order to compute reaction rates
!                !> We compùte equilibrium reaction rates from mass balance equation
!                    call target_waters_new(tar_wat_indices(i))%compute_r_eq_aq_chem(target_waters_old(tar_wat_indices(i)),mixing_ratios_vec,target_waters_old(mix_wat_indices),time_discr_tpt%get_Delta_t(k),F_mat%diag(tar_sol_indices(i)))
!                    !call target_waters_new(tar_wat_indices(i))%compute_r_eq_tilde_aq_chem(tar_wat_indices(i),target_waters_old(tar_wat_indices(i)),mixing_ratios(:,i),target_waters_old(mix_wat_indices(:,i)),time_discr_tpt%get_Delta_t(k),F_mat%diag(tar_sol_indices(i)),r_eq_tilde(:,tar_sol_indices(i)))
!                    !call target_waters_new(tar_wat_indices(i))%solid_chemistry%compute_mass_bal_mins(target_waters_new(tar_wat_indices(i))%r_eq,time_discr_tpt%get_Delta_t(k))
!                    !call target_waters_new(tar_wat_indices(i))%solid_chemistry%compute_conc_solids_iter_EE(target_waters_new(tar_wat_indices(i))%r_eq,time_discr_tpt%get_Delta_t(k))
!                    !do j=1,n-1
!                    !    call target_waters_old(mix_wat_indices(j,i))%compute_molalities() !> in order to speciate
!                    !end do
!                    !call target_waters_old(tar_wat_indices(i))%compute_molalities() !> in order to speciate
!                    !call target_waters_new(tar_wat_indices(i))%compute_molalities() !> in order to speciate
!                    !deallocate(mixing_ratios_vec,mix_wat_indices)
!                end do
!                !do j=1,this%reactive_zones(l)%num_eq_reactions
!                !    call LU_lin_syst(transpose(mixing_ratios),r_eq_tilde(j,:),this%target_waters(1)%CV_params%zero,r_eq_mat(j,:)) !> chapuza en la tolerancia
!                !end do
!                !do j=1,this%reactive_zones(l)%num_eq_reactions
!                !    call LU_lin_syst(transpose(mixing_ratios),r_eq_tilde(j,:),this%target_waters(1)%CV_params%zero,r_eq_mat(j,:)) !> chapuza en la tolerancia
!                !end do
!                !> chapuza
!                do i=1,num_tar_wat
!                    !target_waters_new(tar_wat_indices(i))%r_eq=r_eq_mat(:,i) !> falta un set
!                    call target_waters_new(tar_wat_indices(i))%solid_chemistry%compute_mass_bal_mins(target_waters_new(tar_wat_indices(i))%r_eq,time_discr_tpt%get_Delta_t(k))
!                    call target_waters_new(tar_wat_indices(i))%solid_chemistry%compute_conc_solids_iter_EE(target_waters_new(tar_wat_indices(i))%r_eq,time_discr_tpt%get_Delta_t(k))
!                    !do j=1,n-1
!                    !    call target_waters_old(mix_wat_indices(j,i))%compute_molalities() !> in order to speciate
!                    !end do
!                    !call target_waters_old(tar_wat_indices(i))%compute_molalities() !> in order to speciate
!                    !call target_waters_new(tar_wat_indices(i))%compute_molalities() !> in order to speciate
!                end do
!                target_waters_old_old=target_waters_old
!                target_waters_old=target_waters_new
!            end do
!            deallocate(tar_sol_indices,tar_wat_indices)
!        end do
!    else
!        if (this%chem_syst%num_eq_reacts>0) then
!            if (int_method_chem_reacts==1) then
!                !p_solver=>water_mixing_iter_EE_eq_kin
!            else if (int_method_chem_reacts==2) then
!                !p_solver=>water_mixing_iter_EfI_eq_kin
!            else
!                error stop "Integration method for chemical reactions not implemented yet"
!            end if
!        else if (this%chem_syst%num_kin_reacts>0) then
!            if (int_method_chem_reacts==1) then
!                !p_solver=>water_mixing_iter_EE_kin 
!            else if (int_method_chem_reacts==2) then
!                !p_solver=>water_mixing_iter_EfI_kin
!            else
!                error stop "Integration method for chemical reactions not implemented yet"
!            end if
!        else
!            p_solver=>transport_iter_species_EE_aq_chem
!        end if
!        allocate(tar_wat_indices(this%num_target_waters))
!        do i=1,this%num_target_waters
!            tar_wat_indices(i)=i
!        end do
!        do k=1,time_discr_tpt%Num_time
!            do i=1,this%num_target_waters
!            !> Esto va fuera del time loop
!                mix_wat_indices=mixing_waters_indices%cols(tar_wat_indices(i))%col_1
!                mixing_ratios_vec=mixing_ratios%cols(tar_wat_indices(i))%col_1
!                call p_solver(target_waters_new(tar_wat_indices(i)),mixing_ratios_vec,target_waters_old(mix_wat_indices),F_mat%diag(tar_wat_indices(i)),time_discr_tpt%get_Delta_t(k))
!                !deallocate(mix_wat_indices,mixing_ratios_vec)
!            end do
!            target_waters_old=target_waters_new
!        end do
!        deallocate(tar_wat_indices)
!    end if
!    this%target_waters=target_waters_new
!end subroutine