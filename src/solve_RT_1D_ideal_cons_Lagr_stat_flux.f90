!> \file solve_RT_1D_ideal_cons_Lagr_stat_flux.f90
!> \brief Solve 1D reactive transport with ideal speciation, consistent WMA, and Lagrangian particle tracking
!> \details
!>   Solves the 1D reactive transport problem using:
!>   - Ideal speciation (gamma = 1)
!>   - Consistent Water Mixing Approach (WMA)
!>   - Lagrangian particle tracking with stationary flux
!>
!>   Operator splitting decouples transport (advection-dispersion via particle mixing)
!>   from reactions (equilibrium/kinetic chemistry at each target).
!>
!> \param[inout] this Reactive transport object (transient)
!> \param[in] dir  Directory for output files
!> \param[in] root Root name for output file
!> \author Jordi Petchame-Guerrero
!> \date November 2025
subroutine solve_RT_ideal_cons_Lagr_stat_flux_1D(this, dir, root)
    use RT_m, only: RT_1D_transient_c, move_particles_stat_flux_1D, move_particles_stat_flux_EC_1D
    use aqueous_chemistry_m, only: &
        reactive_mixing_iter_EI_eq_anal_ideal, &
        reactive_mixing_iter_EI_kin_anal_ideal_opt2, &
        compute_r_tilde_impl_opt1, compute_r_tilde_impl_opt2, compute_c_mix, &
        compute_r_tilde_impl_opt3, compute_r_tilde_impl_opt4, &
        reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2
    implicit none

    ! Arguments
    class(RT_1D_transient_c) :: this
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: root

    ! Loop counters
    integer(kind=4) :: i, j, k, kk, ii, j_mix

    ! Problem dimensions
    integer(kind=4) :: n_p, n_v, n_eq
    integer(kind=4) :: num_can_vec, num_non_can_vec
    integer(kind=4) :: num_mix_loc, num_lump
    integer(kind=4) :: out_unit

    ! Index arrays
    integer(kind=4), allocatable :: ind_can_vec(:), ind_non_can_vec(:)

    ! Time and mixing
    real(kind=8) :: time, Delta_t, lambda_R

    ! Concentration arrays
    real(kind=8), allocatable :: c_hat(:), c_mix(:), conc_comp(:)
    real(kind=8), allocatable :: conc_nc(:), r_tilde(:)
    real(kind=8), allocatable :: conc_old_mix(:,:)
    integer(kind=4), allocatable :: ind_aq_sp_mix(:,:)

    ! Procedure pointers
    procedure(reactive_mixing_iter_EI_kin_anal_ideal_opt2), pointer :: p_solver => null()
    procedure(compute_r_tilde_impl_opt1), pointer :: p_r_tilde => null()
    procedure(move_particles_stat_flux_1D), pointer :: p_displace_particles => null()

    ! Cached sizes for allocation optimisation
    integer(kind=4) :: prev_n_v, prev_n_p

    ! Count equilibrium reactions across all target waters
    n_eq = 0
    do i = 1, this%chemistry%num_target_waters
        n_eq = n_eq + this%target_waters(i)%aq_chem%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    end do

    select type (this)
    type is (RT_1D_transient_c)
        time = 0d0
        kk = 2
        ii = 1
        this%chemistry%chem_out_options%time_steps(this%chemistry%chem_out_options%num_time_steps) = &
            this%transport%time_discr%Num_time
        open(newunit=out_unit, file=dir//root//'.output', form="formatted", access="sequential", status="unknown")

        !--- Solver selection based on chemical system type ---
        if (this%chemistry%chem_syst%num_kin_reacts > 0 .and. n_eq > 0) then
            ! Equilibrium + kinetic reactions
            if (this%chemistry%rk_avg_opt == 1) then
                error stop "rk average option 1 not implemented yet"
            else if (this%chemistry%rk_avg_opt == 2) then
                p_solver => reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2
            else
                error stop "rk average option not implemented yet"
            end if
        else if (n_eq > 0) then
            ! Equilibrium-only reactions
            if (this%int_method_chem_reacts == 2) then
                p_solver => reactive_mixing_iter_EI_eq_anal_ideal
            else
                error stop "Integration method for equilibrium-only reactions not implemented yet"
            end if
        else
            ! Kinetic-only reactions
            if (this%chemistry%rk_avg_opt == 1) then
                error stop "rk average option 1 not implemented yet"
            else if (this%chemistry%rk_avg_opt == 2) then
                p_solver => reactive_mixing_iter_EI_kin_anal_ideal_opt2
            else
                error stop "rk average option not implemented yet"
            end if
        end if

        ! Downstream reaction rate estimation (common to all chemical system types)
        select case (this%chemistry%r_down_opt)
        case (1)
            p_r_tilde => compute_r_tilde_impl_opt1
        case (2)
            p_r_tilde => compute_r_tilde_impl_opt2
        case (3)
            p_r_tilde => compute_r_tilde_impl_opt3
        case (4)
            p_r_tilde => compute_r_tilde_impl_opt4
        case default
            error stop "Downstream reaction rate estimation option not implemented yet"
        end select

        ! Particle displacement method
        if (this%transport%spatial_discr%targets_flag == 0) then
            p_displace_particles => move_particles_stat_flux_1D
        else if (this%transport%spatial_discr%targets_flag == 1) then
            p_displace_particles => move_particles_stat_flux_EC_1D
        else
            error stop "Stationary flux Lagrangian scheme not implemented yet"
        end if

        ! Extract canonical vectors from mixing ratios matrix
        call this%transport%mixing_ratios_conc%get_can_vec( &
            this%chemistry%CV_params%abs_tol, num_can_vec, ind_can_vec, num_non_can_vec, ind_non_can_vec)

        ! Initialise concentration history
        do i = 1, this%chemistry%num_waters
            call this%chemistry%waters(i)%set_conc_old()
            call this%chemistry%waters(i)%set_conc_old_old()
            call this%chemistry%waters(i)%solid_chemistry%set_conc_old()
            call this%chemistry%waters(i)%set_solid_chemistry_old()
        end do

        ! Write output file header
        write(out_unit, "(2x,'Species :',5x,*(A15))") &
            (this%chemistry%chem_syst%aq_phase%aq_species( &
            this%chemistry%chem_out_options%ind_aq_species(j))%name, &
            j=1, this%chemistry%chem_out_options%num_aq_species)

        ! Initialise allocation tracking for inner loop
        prev_n_v = -1
        prev_n_p = -1

        !====================================================================
        ! Main time loop
        !====================================================================
        do k = 1, this%transport%time_discr%Num_time
            Delta_t = this%transport%time_discr%get_Delta_t(k)
            time = time + Delta_t

            if (k == this%chemistry%chem_out_options%time_steps(kk)) then
                write(out_unit, "(/,2x,'t = ',*(ES15.5),/)") time
            end if

            ! Shift time levels: current -> old, old -> old_old
            do i = 1, this%chemistry%num_target_waters
                call this%target_waters(i)%update_old_attributes()
            end do

            call p_displace_particles(this, k)
            call this%introduce_particle(k)

            !----------------------------------------------------------------
            ! Target waters loop (non-canonical vectors only)
            !----------------------------------------------------------------
            do i = 1, num_non_can_vec
                associate( &
                    idx   => ind_non_can_vec(i), &
                    aq    => this%target_waters(ind_non_can_vec(i))%aq_chem, &
                    sc    => this%target_waters(ind_non_can_vec(i))%aq_chem%solid_chemistry, &
                    rz    => this%target_waters(ind_non_can_vec(i))%aq_chem%solid_chemistry%reactive_zone, &
                    sa    => this%target_waters(ind_non_can_vec(i))%aq_chem%solid_chemistry%reactive_zone%speciation_alg, &
                    mc    => this%transport%mix_conc_indices%cols(ind_non_can_vec(i)), &
                    mr    => this%transport%mix_react_indices%cols(ind_non_can_vec(i)), &
                    mr_c  => this%transport%mixing_ratios_conc%cols(ind_non_can_vec(i)), &
                    mr_Ri => this%transport%mixing_ratios_R_init%cols(ind_non_can_vec(i)), &
                    mr_R  => this%transport%mixing_ratios_R%cols(ind_non_can_vec(i)) &
                )

                n_p = sa%num_prim_species
                n_v = sa%num_var_act_species

                ! Reallocate only when sizes change
                if (n_v /= prev_n_v .or. n_p /= prev_n_p) then
                    if (allocated(c_hat)) deallocate(c_hat, conc_comp, r_tilde)
                    allocate(c_hat(n_v), conc_comp(n_p), r_tilde(n_v))
                    prev_n_v = n_v
                    prev_n_p = n_p
                end if

                ! Validate target water ordering
                if (mc%col_1(1) /= this%chemistry%tar_wat_indices(idx)) then
                    print *, mc%col_1(1), this%chemistry%tar_wat_indices(idx)
                    print *, "Problematic target water", idx
                    error stop "Target waters not in the right order in mixing waters indices"
                end if

                ! Transport mixing: compute mixed concentrations
                num_mix_loc = mc%dim - 3
                allocate(conc_old_mix( &
                    size(this%chemistry%waters(mc%col_1(1))%conc_old), num_mix_loc))
                allocate(ind_aq_sp_mix( &
                    size(this%chemistry%waters(mc%col_1(1))%indices_aq_species), num_mix_loc))
                do j_mix = 1, num_mix_loc
                    conc_old_mix(:, j_mix) = this%chemistry%waters(mc%col_1(j_mix + 1))%conc_old
                    ind_aq_sp_mix(:, j_mix) = this%chemistry%waters(mc%col_1(j_mix + 1))%indices_aq_species
                end do
                call compute_c_mix(this%chemistry%waters(mc%col_1(1)), &
                    conc_old_mix, ind_aq_sp_mix, mr_c%col_1, c_mix)
                deallocate(conc_old_mix, ind_aq_sp_mix)

                ! Compute mixed reaction rates
                call p_r_tilde( &
                    this%chemistry%waters(mr%col_1(2:mr%dim - 2)), &
                    rz%ind_var_act_species, sa%comp_mat, &
                    mr_Ri%col_1(2:), &
                    mr%col_1(mr%dim - 1), mr%col_1(mr%dim), &
                    this%transport%time_discr%theta_r, Delta_t, r_tilde)

                ! Lumping of reaction mixing ratios
                call aq%modify_mix_ratios_reacts( &
                    mr_Ri%col_1, c_mix, Delta_t, r_tilde, mr_R%col_1, num_lump)
                this%chemistry%num_lump = this%chemistry%num_lump + num_lump

                ! Update concentrations with reaction contributions
                c_hat = c_mix + Delta_t * r_tilde
                conc_nc = aq%get_conc_nc()

                ! Compute reaction mixing ratio
                if (this%chemistry%r_down_opt == 4) then
                    lambda_R = this%transport%compute_mix_ratio_R_opt4(idx)
                else
                    lambda_R = mr_R%col_1(1)
                end if

                ! Reactive mixing iteration
                call p_solver(aq, aq%get_c1_old_old(), c_hat, &
                    mr_R%col_1(1), lambda_R, Delta_t, &
                    this%transport%time_discr%theta_r, conc_nc, conc_comp)

                ! Update derived quantities
                call aq%compute_pH()
                call aq%compute_salinity()
                call aq%compute_ionic_strength()

                ! Accumulate reaction amounts
                aq%Rk_accum = aq%Rk_accum + aq%Rk
                sc%Rk_accum = sc%Rk_accum + sc%Rk

                ! Compute equilibrium reaction rates
                if (sa%num_eq_reactions > 0 .and. aq%indices_rk%num_cols > 0) then
                    call aq%compute_Re_kin(c_hat(n_p + 1:n_v), Delta_t, lambda_R)
                else if (sa%num_eq_reactions > 0) then
                    call aq%compute_Re(c_hat(1:n_v), Delta_t, lambda_R)
                end if

                ! Update mineral state
                if (associated(sc%mineral_zone)) then
                    call sc%compute_mass_bal_mins(Delta_t)
                    call sc%compute_conc_minerals_iter(Delta_t)
                end if

                ! Update gas state
                if (associated(aq%gas_chemistry)) then
                    call aq%gas_chemistry%compute_conc_gases_iter( &
                        Delta_t, aq%volume, [aq%re_mean, aq%rk_mean])
                    call aq%gas_chemistry%compute_vol_gas_species_conc()
                    call aq%gas_chemistry%compute_log_act_coeffs_gases()
                end if

                ! Write output at specified time steps and locations
                if ((k == this%chemistry%chem_out_options%time_steps(kk)) .and. &
                    (idx == this%chemistry%chem_out_options%ind_waters(ii))) then
                    write(out_unit, "(I10,*(ES15.5))") &
                        this%chemistry%chem_out_options%ind_waters(ii), &
                        (aq%concentrations(aq%indices_aq_species( &
                        this%chemistry%chem_out_options%ind_aq_species(j))), &
                        j=1, this%chemistry%chem_out_options%num_aq_species)
                    if (ii < this%chemistry%chem_out_options%num_waters) then
                        ii = ii + 1
                    else if (kk < this%chemistry%chem_out_options%num_time_steps) then
                        kk = kk + 1
                        ii = 1
                    end if
                end if

                deallocate(c_mix, conc_nc)
                end associate
            end do ! target waters
        end do ! time loop

        if (allocated(c_hat)) deallocate(c_hat, conc_comp, r_tilde)
        close(out_unit)
    end select

    print *, "Number of lumpings performed: ", this%chemistry%num_lump
end subroutine
