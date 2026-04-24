!> \file solve_RT_1D_ideal_lump_Lagr_trans_flux.f90
!> \brief Solves 1D reactive transport using Lagrangian formulation with lumped WMA and transient flux
!> \details This subroutine implements a coupled transport-chemistry solver for 1D reactive transport
!> problems using:
!> - **Lagrangian formulation**: Particles move with the flow
!> - **Lumped Water Mixing Algorithm (WMA)**: Reaction amounts are lumped to reduce stiffness
!> - **Transient flux**: Flow field can change with time
!> - **Ideal mixing**: Activity coefficients assumed constant during mixing
!>
!> The algorithm follows this sequence for each time step:
!> 1. Update transport properties (if flux changes)
!> 2. Compute mixing ratios based on Lagrangian transport
!> 3. Mix conservative components using mixing ratios
!> 4. Solve chemical reactions (equilibrium + kinetic)
!> 5. Update mineral/gas phase compositions
!> 6. Advance particles to new positions
!> 7. Introduce new particles at boundaries
!>
!> The lumped WMA reduces the number of canonical vectors (independent waters) by lumping
!> similar mixing ratios, which improves computational efficiency for large systems.
!>
!> \param[in,out] this Reactive transport object containing all simulation data
!> \param[in] root Root filename prefix for output files (e.g., "simulation" → "simulation.output")

subroutine solve_RT_ideal_lump_Lagr_trans_flux_1D(this,dir,root)
    use RT_m, only: RT_1D_transient_c, RT_1D_stat_c, RT_c !> Import RT_1D module types for reactive transport
    use aqueous_chemistry_m, only: aqueous_chemistry_c, compute_r_tilde_impl_opt1, compute_c_mix, &
        reactive_mixing_iter_EE_eq_kin_ideal, reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2, &
        mixing_iter_comp_ideal, reactive_mixing_iter_EE_kin, reactive_mixing_iter_EI_kin_anal_ideal_opt2 !> Import aqueous chemistry types and functions
    implicit none !> Enforce explicit variable declarations
    class(RT_1D_transient_c) :: this !> Reactive transport object containing transport, chemistry, and discretization data
    character(len=*), intent(in) :: dir !> Directory for output files
    character(len=*), intent(in) :: root !> Root filename prefix for output files (e.g., "problem" → "problem.output")
    !integer(kind=4), intent(in) :: unit !> file unit (commented out - defined internally instead)
    !class(real_array_c), intent(in) :: mixing_ratios !> mixing ratios matrix for concentrations (commented out)
    !class(real_array_c), intent(in) :: mixing_ratios_R_init !> initial mixing ratios matrix for reaction amounts (commented out)
    !class(int_array_c), intent(in) :: this%transport%mix_conc_indices !> matrix of target water indices that mix (commented out)
    !class(int_array_c), intent(in) :: this%transport%mix_react_indices !> matrix of domain water indices that mix (commented out)
    !real(kind=8), intent(in) :: F_mat(:) !> storage matrix diagonal (commented out)
    !class(time_discr_c), intent(in) :: time_discr_tpt !> time discretisation for transport (commented out)
    !real(kind=8), intent(in) :: theta_r !> integration method for chemical reactions (commented out)
    !class(real_array_c), intent(inout) :: mixing_ratios_R !> final mixing ratios matrix for reaction amounts (commented out)
    integer(kind=4) :: i !> Loop counter for target waters: i = 1, 2, ..., num_waters [-]
    integer(kind=4) :: j !> Loop counter for target solids, species, or reactions [-]
    integer(kind=4) :: l !> Loop counter for reactive zones: l = 1, 2, ..., num_reactive_zones [-]
    integer(kind=4) :: k !> Loop counter for time steps: k = 1, 2, ..., Num_time [-]
    integer(kind=4) :: kk !> Counter for chemistry output time steps (indexed into chem_out_options%time_steps) [-]
    integer(kind=4) :: ii !> Counter for chemistry output target waters (indexed into chem_out_options%ind_waters) [-]
    integer(kind=4) :: num_tar_wat !> Number of target waters in the domain [-]
    integer(kind=4) :: num_tar_sol !> Number of target solids in the domain [-]
    integer(kind=4) :: n_p !> Number of primary species in speciation algorithm [-]
    integer(kind=4) :: n_v !> Total number of variable activity species (aqueous + solid) [-]
    integer(kind=4) :: n_v_aq !> Number of aqueous variable activity species [-]
    integer(kind=4) :: n_v_aq_2 !> Number of aqueous secondary variable activity species (complexes) [-]
    integer(kind=4) :: mix_ind !> Starting index for mixing waters (1 for explicit, 2 for implicit) [-]
    integer(kind=4) :: num_can_vec !> Number of canonical vectors in mixing ratios (independent waters) [-]
    integer(kind=4) :: num_non_can_vec !> Number of non-canonical vectors (dependent waters requiring mixing) [-]
    integer(kind=4) :: unit !> File unit number for output file (typically 7) [-]
    integer(kind=4) :: k_flux !> Counter for flux changes in transient flow field [-]
    integer(kind=4), allocatable :: tar_gas_indices(:) !> Indices of target gases in each reactive zone [-]
    integer(kind=4), allocatable :: tar_sol_indices(:) !> Indices of target solids in each reactive zone [-]
    integer(kind=4), allocatable :: tar_wat_indices(:) !> Indices of target waters in each reactive zone [-]
    integer(kind=4), allocatable :: perm(:) !> Permutation vector for reordering aqueous concentrations [-]
    integer(kind=4), allocatable :: ind_can_vec(:) !> Indices of canonical vectors in mixing ratios array [-]
    integer(kind=4), allocatable :: ind_non_can_vec(:) !> Indices of non-canonical vectors requiring mixing [-]
    REAL(KIND=8) :: time !> Current simulation time [T]
    REAL(KIND=8) :: Delta_t !> Time step size for current iteration [T]
    REAL(KIND=8) :: theta !> Time integration weighting factor: 0=explicit Euler, 1=implicit Euler, 0.5=Crank-Nicolson [-]
    REAL(KIND=8) :: y !> Sum of upstream mixing ratios (for verification/debugging) [-]
    REAL(KIND=8), allocatable :: c_hat(:) !> Variable activity species concentrations after mixing but before reaction [M/L³]
    real(kind=8), allocatable :: rk_tilde_up(:) !> Kinetic reaction rates after mixing from upstream waters [1/T]
    real(kind=8), allocatable :: rk_tilde_down(:) !> Kinetic reaction rates after mixing from downstream waters [1/T]
    real(kind=8), allocatable :: rk_tilde(:) !> Effective kinetic reaction rates after mixing [1/T]
    REAL(KIND=8), allocatable :: conc_comp(:) !> Concentrations before mixing for all mixing waters (n_v × num_mix_waters) [M/L³]
    REAL(KIND=8), allocatable :: conc_nc(:) !> Concentrations of variable activity species after reaction [M/L³]
    type(aqueous_chemistry_c), allocatable :: waters_new(:) !> Target waters at time step k+1 (after reaction)
    type(aqueous_chemistry_c), allocatable :: waters_old(:) !> Target waters at time step k (current state)
    type(aqueous_chemistry_c), allocatable :: waters_old_old(:) !> Target waters at time step k-1 (previous state)
    type(aqueous_chemistry_c), allocatable :: mix_waters(:) !> Subset of waters participating in mixing at current cell
    integer(kind=4) :: num_lump !> Number of lumpings performed in current reactive mixing iteration [-]
    integer(kind=4) :: j_mix !> loop index for mixing waters
    integer(kind=4) :: num_mix_loc !> number of mixing waters
    real(kind=8), allocatable :: conc_old_mix(:,:) !> conc_old from each mixing water
    integer(kind=4), allocatable :: ind_aq_sp_mix(:,:) !> indices_aq_species from each mixing water
    procedure(reactive_mixing_iter_EE_eq_kin_ideal), pointer :: p_solver=>null() !> Procedure pointer to reactive mixing solver (selected based on problem type)
    procedure(compute_r_tilde_impl_opt1), pointer :: compute_r_tilde=>null() !> Procedure pointer for computing kinetic reaction rates after mixing
!**************************************************************************************************
!> \name Initialization
!> \brief Initialize target waters for time integration
!**************************************************************************************************
    waters_old=this%chemistry%waters !> Copy current target waters to old state (time step k)
    waters_old_old=waters_old !> Copy old state to old-old state (time step k-1)

!**************************************************************************************************
!> \name Type Selection and Initialization
!> \brief Select RT_1D object type and initialize simulation parameters
!**************************************************************************************************
select type (this) !> Polymorphic type selection for reactive transport object
type is (RT_1D_transient_c) !> Transient reactive transport case
    time=0d0 !> Initialize simulation time to zero [T]
    k_flux=1 !> Initialize flux change counter to 1 (first flux state) [-]
    kk=2 !> Initialize chemistry output time step counter to 2 (first output at kk=1 is initial condition) [-]
    ii=1 !> Initialize chemistry output target water counter to 1 [-]
    this%chemistry%chem_out_options%time_steps(this%chemistry%chem_out_options%num_time_steps)=this%transport%time_discr%Num_time !> Set last output time step to final simulation time step [-]
!> Commented out chapuza (workaround) for allocating reaction rate arrays
    ! if (this%chemistry%chem_syst%num_kin_reacts==this%chemistry%chem_syst%num_aq_kin_reacts) then
    !     allocate(rk_tilde(this%chemistry%chem_syst%num_kin_reacts))
    !     allocate(rk(this%chemistry%chem_syst%num_kin_reacts))
    ! else
    !     continue
    ! end if
    unit=7 !> Set output file unit number to 7 [-]
    open(unit,file=dir//root//'.output',form="formatted",access="sequential",status="unknown") !> Open formatted output file for writing results

!**************************************************************************************************
!> \name Solver Selection
!> \brief Select reactive mixing subroutine based on chemical system and integration method
!> \details The solver is chosen based on:
!> - Presence of equilibrium and/or kinetic reactions
!> - Time integration method (explicit Euler, implicit Euler, Crank-Nicolson)
!> - Jacobian computation method (analytical vs. incremental coefficients)
!> - Reaction rate averaging option
!**************************************************************************************************
    if (this%chemistry%chem_syst%num_kin_reacts>0 .and. this%chemistry%chem_syst%speciation_alg%num_eq_reactions>0) then !> Case: Both equilibrium and kinetic reactions present
        if (this%int_method_chem_reacts==1) then !> Explicit Euler method
            theta=0d0 !> Set time weighting to 0 (explicit) [-]
            !p_solver=>reactive_mixing_iter_EE_eq_kin_ideal_lump !> Point to explicit Euler solver with equilibrium and kinetic reactions
            compute_r_tilde=>compute_r_tilde_impl_opt1 !> Point to explicit kinetic rate computation
            mix_ind=1 !> Start mixing from index 1 (include all waters) [-]
            !mix_waters_indices=>this%transport%mix_conc_indices
    else if (this%chemistry%Jac_opt==1) then !> Implicit/Crank-Nicolson with analytical Jacobian
            mix_ind=2 !> Start mixing from index 2 (exclude one water for implicit schemes) [-]
            if (this%int_method_chem_reacts==2) then !> Implicit Euler method
                theta=1d0 !> Set time weighting to 1 (fully implicit) [-]
            else if (this%int_method_chem_reacts==3) then !> Crank-Nicolson method
                theta=5d-1 !> Set time weighting to 0.5 (Crank-Nicolson) [-]
            else !> Unsupported integration method
                error stop "Integration method for chemical reactions not implemented yet"
            end if
            if (this%chemistry%rk_avg_opt==1) then !> Reaction rate averaging option 1
                !p_solver=>reactive_mixing_iter_EI_eq_kin_anal_ideal_opt1 !> Option 1 not implemented
            else if (this%chemistry%rk_avg_opt==2) then !> Reaction rate averaging option 2
                !p_solver=>reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2 !> Point to implicit/CN solver with option 2 averaging
            else !> Unsupported averaging option
                error stop "rk average option not implemented yet"
            end if
            !> Commented out options for different WMA consistency approaches
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
        !> Commented out Crank-Nicolson implementation (alternative approach)
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
        else !> Unsupported integration method
            error stop "Integration method for chemical reactions not implemented yet"
        end if
    else if (this%chemistry%chem_syst%speciation_alg%num_eq_reactions>0) then !> Case: Only equilibrium reactions (no kinetics)
        !p_solver=>mixing_iter_comp_ideal !> Point to equilibrium-only solver with lumped WMA
    else !> Case: Only kinetic reactions (no equilibrium)
        if (this%int_method_chem_reacts==1) then !> Explicit Euler method for kinetic reactions only
            theta=0d0 !> Set time weighting to 0 (explicit) [-]
            p_solver=>reactive_mixing_iter_EE_kin !> Point to explicit Euler solver for kinetic reactions only
            mix_ind=1 !> Start mixing from index 1 (include all waters) [-]
    else if (this%chemistry%Jac_opt==1) then !> Implicit/Crank-Nicolson with analytical Jacobian for kinetic reactions
            if (this%int_method_chem_reacts==2) then !> Implicit Euler method
                theta=1d0 !> Set time weighting to 1 (fully implicit) [-]
            else if (this%int_method_chem_reacts==3) then !> Crank-Nicolson method
                theta=5d-1 !> Set time weighting to 0.5 (Crank-Nicolson) [-]
            else !> Unsupported integration method
                error stop "Integration method for chemical reactions not implemented yet"
            end if
            if (this%chemistry%rk_avg_opt==1) then !> Reaction rate averaging option 1
                !p_solver=>reactive_mixing_iter_EI_kin_anal_ideal_opt1 !> Option 1 not implemented
            else if (this%chemistry%rk_avg_opt==2) then !> Reaction rate averaging option 2
                !p_solver=>reactive_mixing_iter_EI_kin_anal_ideal_opt2 !> Point to implicit/CN solver for kinetics with option 2 averaging
            else !> Unsupported averaging option
                error stop "rk average option not implemented yet"
            end if
            mix_ind=2 !> Start mixing from index 2 (exclude one water for implicit schemes) [-]
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
        else !> Unsupported integration method for kinetic reactions
            error stop "Integration method for chemical reactions not implemented yet"
        end if
    end if
    
!**************************************************************************************************
!> \name Vector Classification
!> \brief Identify canonical and non-canonical vectors in mixing ratios
!> \details Canonical vectors are independent (representing distinct waters) while non-canonical
!> vectors are linear combinations of canonical ones. Lumping reduces the number of canonical
!> vectors, improving computational efficiency.
!**************************************************************************************************
    call this%transport%mixing_ratios_conc%get_can_vec(this%chemistry%CV_params%abs_tol,num_can_vec,ind_can_vec,num_non_can_vec,&
        ind_non_can_vec) !> Classify vectors based on coefficient of variation tolerance
    
!**************************************************************************************************
!> \name Concentration Storage Initialization
!> \brief Store old concentrations for time integration
!**************************************************************************************************
    do i=1,this%chemistry%num_waters !> Loop over all target waters
        call this%chemistry%waters(i)%set_conc_old() !> Store current aqueous concentrations as "old" (time step k)
        call this%chemistry%waters(i)%set_conc_old_old() !> Store previous concentrations as "old_old" (time step k-1)
        call this%chemistry%waters(i)%solid_chemistry%set_conc_old() !> Store current solid concentrations as "old"
        call this%chemistry%waters(i)%solid_chemistry%set_conc_old_old() !> Store previous solid concentrations as "old_old"
    end do
    !> Commented out alternative approach for non-canonical vectors only
    !do i=1,num_non_can_vec
    !    call this%chemistry%waters(ind_can_vec(i))%set_conc_old() !>
    !    call this%chemistry%waters(ind_can_vec(i))%solid_chemistry%set_conc_old() !>
    !    call this%chemistry%waters(ind_can_vec(i))%solid_chemistry%set_conc_old_old() !>
    !end do
    
!**************************************************************************************************
!> \name Main Time Loop
!> \brief Advance solution through time for transient transport and chemistry
!> \details For each time step:
!> 1. Update transport properties (if flux changes)
!> 2. Mix conservative components
!> 3. Solve chemical reactions
!> 4. Update mineral/gas phases
!> 5. Advance particles
!> 6. Write output
!**************************************************************************************************
        do k=1,this%transport%time_discr%Num_time !> Loop over all time steps from k=1 to final time
            Delta_t=this%transport%time_discr%get_Delta_t(k) !> Get time step size for current iteration [T]
            time=time+Delta_t !> Advance current time by Delta_t [T]
            !> Write output header for current time step if this is an output time
            if (k==this%chemistry%chem_out_options%time_steps(kk)) then !> Check if current time step is scheduled for output
               write(unit,"(2x,'t = ',*(ES15.5))") time !> Write current simulation time [T]
               write(unit,"(A10,*(A15))") 'Water', (this%chemistry%chem_syst%aq_phase%aq_species( &
                this%chemistry%chem_out_options%ind_aq_species(j))%name, &
                j=1,this%chemistry%chem_out_options%num_aq_species) !> Write column headers with selected aqueous species names
            end if
            !> Update flux if time-dependent flux changes at current time step
            if (this%transport%tpt_props_heterog%flux_trans%time_ind(k_flux)==k) then !> Check if flux changes at time step k
                call this%transport%tpt_props_heterog%update_flux_int_trans(k_flux) !> Update interpolated transient flux values
                call this%transport%stab_params_tpt%compute_stab_params_tpt_1D(this%transport%tpt_props_heterog,&
                    this%transport%spatial_discr,Delta_t) !> Recompute stability parameters with new flux and time step
                call this%transport%update_mixing_ratios_Delta_t_homog(Delta_t) !> Update mixing ratios for new flux conditions
                if (k_flux<this%transport%tpt_props_heterog%flux_trans%ntime) then !> Check if more flux changes remain
                    k_flux=k_flux+1 !> Increment flux change counter to prepare for next change [-]
                end if
            end if
            !> Update concentrations and reaction rates from previous time step to current
            do i=1,this%chemistry%num_target_waters !> Loop over domain target waters
                call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%update_conc_old() !> Shift conc_old_old ← conc_old ← conc for aqueous phase
                call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%solid_chemistry%update_conc_old() !> Shift conc_old_old ← conc_old ← conc for solid phase
                call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%update_rk_old() !> Shift rk_old ← rk for aqueous kinetic rates [1/T]
                call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%solid_chemistry%update_rk_old() !> Shift rk_old ← rk for solid kinetic rates [1/T]
            end do
            
!**************************************************************************************************
!> \name Reactive Mixing Loop
!> \brief Solve reactive mixing for each domain target water
!**************************************************************************************************
            do i=1,this%chemistry%num_target_waters !> Loop over all domain target waters
                !> Extract speciation algorithm dimensions for current target water's reactive zone
                n_p=this%chemistry%waters(&
                    this%chemistry%tar_wat_indices(&
                    i))%solid_chemistry%reactive_zone%speciation_alg%num_prim_species !> Number of primary species [-]
                n_v=this%chemistry%waters(this%chemistry%tar_wat_indices(i)&
                    )%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species !> Total number of variable activity species (aq + solid) [-]
                n_v_aq=this%chemistry%waters(this%chemistry%tar_wat_indices(i)&
                    )%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species !> Number of aqueous variable activity species [-]
                n_v_aq_2=this%chemistry%waters(this%chemistry%tar_wat_indices(i)&
                    )%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species !> Number of secondary aqueous variable activity species [-]
                
                !> Commented out canonical vector handling logic
                !if (num_can_vec>0 .and. ind_non_can_vec(i)==ind_can_vec(cntr_can_vec)) then
                !    cntr_can_vec=cntr_can_vec+1 !> we update counter of canonical vectors
                !    !call this%chemistry%waters(ind_non_can_vec(i))%compute_rk() !> we compute kinetic reaction rates
                !    conc_nc=this%chemistry%waters(ind_non_can_vec(i))%get_conc_nc()
                !    !> chapuza
                !    !this%chemistry%waters(ind_non_can_vec(i))%Rk_est=0d0
                !    !this%chemistry%waters(ind_non_can_vec(i))%solid_chemistry%Rk_est=0d0
                !    continue
                !else
                    if (this%chemistry%waters(this%chemistry%tar_wat_indices(&
                        i))%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl &
                        >0) then !> Variable activity species include both aqueous and solid (cation exchange present)
                        !p_solver=>mixing_iter_comp_exch_ideal !> Would use cation exchange solver (only equilibrium reactions)
                    end if
                    allocate(conc_nc(n_v)) !> Allocate array for concentrations of variable activity species [M/L³]
                    allocate(conc_comp(n_p)) !> Commented: allocate upstream/downstream kinetic rates
                    !allocate(rk_tilde(n_v)) !> Commented: allocate mixed kinetic rates
                    !allocate(conc_old(n_v,mixing_ratios_conc%cols(i)%dim)) !> Commented: allocate old concentrations matrix (chapuza)
                    !> Verify that mixing water indices are correctly ordered
                    if (this%transport%mix_conc_indices%cols(ind_non_can_vec(i))%col_1(1)/=this%chemistry%tar_wat_indices(&
                        i)) then !> Check if first mixing water index matches current domain water index
                        print *, this%transport%mix_conc_indices%cols(ind_non_can_vec(i))%col_1(1),&
                            this%chemistry%tar_wat_indices(i) !> Print mismatched indices for debugging
                        error stop "Target waters not in the right order" !> Halt execution if indices don't match
                    end if
                    
!> \subsection Transport Mixing Step
!> \brief Mix conservative components based on Lagrangian transport
                    !c_hat=waters_old(this%chemistry%tar_wat_indices(i))%get_conc_nc() !> Get variable activity species concentrations from previous time step [M/L³]
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
                    !> Commented out complex lumping logic for kinetic reaction rates (chapuza = workaround)
                    !allocate(mix_waters(this%chemistry%num_target_waters-num_can_vec-mix_ind+1))
                    !do j=1,this%chemistry%num_target_waters-mix_ind+1
                    !    if (j==
                    !    mix_waters(j)=this%chemistry%waters(this%transport%mix_conc_indices%cols(i)%col_1(&
                    !    mix_ind+j-1))
                    !mix_waters=this%chemistry%waters(this%transport%mix_react_indices%cols(ind_non_can_vec(i))%col_1(mix_ind:this%transport%mix_react_indices%cols(ind_non_can_vec(i))%dim-2))
                    !call compute_rk_tilde(mix_waters,mixing_ratios_R_init%cols(ind_non_can_vec(i))%col_1(mix_ind:),&
                    !    this%transport%mix_react_indices%cols(ind_non_can_vec(i))%col_1(this%transport%mix_react_indices%cols(ind_non_can_vec(i))%dim-1),&
                    !    this%transport%mix_react_indices%cols(ind_non_can_vec(i))%col_1(this%transport%mix_react_indices%cols(ind_non_can_vec(i))%dim),theta,Delta_t,rk_tilde)
                    !!> chapuza
                    !this%chemistry%waters(this%transport%mix_react_indices%cols(ind_non_can_vec(i))%col_1(mix_ind:this%transport%mix_react_indices%cols(ind_non_can_vec(i))%dim-2))=&
                    !    mix_waters
                    !call this%chemistry%waters(this%chemistry%tar_wat_indices(ind_non_can_vec(i)))%solid_chemistry%modify_mix_ratios_rk(&
                    !    mixing_ratios_R_init%cols(ind_non_can_vec(i))%col_1(1),c_hat,Delta_t,rk_tilde,&
                    !    mixing_ratios_R%cols(ind_non_can_vec(i))%col_1(1),num_lump)
                    !!> chapuza
                    !this%chemistry%num_lump=this%chemistry%num_lump+num_lump !> we update number of lumpings
                    !y=sum(mixing_ratios_R%cols(i)%col_1(2:1+this%transport%mix_conc_indices%cols(i)%col_1(this%transport%mix_conc_indices%cols(i)%dim-1)))
                    !> Commented out alternative concentration indexing logic
                    !do j=1,n_v_aq
                    !    if (waters_old(ind_non_can_vec(i))%ind_var_act_species(j) /= &
                    !        waters_old(ind_non_can_vec(i))%indices_aq_species(&
                    !        this%chemistry%chem_syst%aq_phase%ind_diss_solids(j))) then
                    !        c_hat(waters_old(ind_non_can_vec(i))%ind_var_act_species(j))=c_hat_aux(&
                    !            waters_old(ind_non_can_vec(i))%indices_aq_species(&
                    !            this%chemistry%chem_syst%aq_phase%ind_diss_solids(j)))
                    !    end if
                    !end do
                    
!> \subsection Chemical Reaction Step
!> \brief Solve chemical reactions for mixed water
                          call p_solver(this%chemistry%waters(this%chemistry%tar_wat_indices(i)),waters_old_old(&
                                this%chemistry%tar_wat_indices(i))%get_c1(),c_hat,this%transport%mixing_ratios_R%cols(ind_non_can_vec(i))%col_1(1),&
                                Delta_t,theta,conc_nc,conc_comp) !> Solve reactive mixing iteration: equilibrium speciation + kinetic reactions → new concentrations [M/L³]
                                
!> \subsection Update Derived Quantities
!> \brief Compute pH, salinity, and ionic strength from updated concentrations
                    call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%compute_pH()
                    call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%compute_salinity()
                    call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%compute_ionic_strength()
                    
!> \subsection Equilibrium Reaction Rates
!> \brief Compute equilibrium reaction rates from mass balance
                    if (this%chemistry%waters(this%chemistry%tar_wat_indices(i)&
                        )%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions>0 .and. &
                        this%chemistry%waters(this%chemistry%tar_wat_indices(i))%indices_rk%num_cols>0) then !> Case: equilibrium reactions AND kinetic reactions present
                        call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%compute_Re_kin(&
                            c_hat(n_p+1:n_v),Delta_t,&
                            this%transport%mixing_ratios_R%cols(ind_non_can_vec(i))%col_1(1)) !> Compute equilibrium reaction rates using lumped WMA [mol/T]
                    else if (this%chemistry%waters(this%chemistry%tar_wat_indices(i)&
                        )%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions>0) then !> Case: only equilibrium reactions (no kinetics)
                        call this%chemistry%waters(this%chemistry%tar_wat_indices(i))%compute_Re(&
                            c_hat(1:n_v),Delta_t,&
                            this%transport%mixing_ratios_R%cols(ind_non_can_vec(i))%col_1(1)) !> Compute equilibrium reaction rates using lumped WMA [mol/T]
                    end if
                    !> Commented out alternative equilibrium rate computation
                    ! call this%chemistry%waters(ind_non_can_vec(i))%compute_re_mean(c_hat(n_p+1:n_p+n_v_aq_2),Delta_t,&
                    !     rk_tilde)
                    
!> \subsection Mineral Phase Update
!> \brief Update mineral concentrations and volumetric fractions
                    if (associated(this%chemistry%waters(this%chemistry%tar_wat_indices(&
                        i))%solid_chemistry%mineral_zone)) then !> Check if mineral zone is associated with current target water
                        call this%chemistry%waters(this%chemistry%tar_wat_indices(&
                            i))%solid_chemistry%compute_mass_bal_mins(&
                            Delta_t) !> Compute mass balance for minerals (dissolution/precipitation) over time step [mol]
                        call this%chemistry%waters(this%chemistry%tar_wat_indices(&
                            i))%solid_chemistry%compute_conc_minerals_iter(Delta_t) !> Update mineral concentrations from mass balance [M/L³]
                    end if
                    
!> \subsection Gas Phase Update
!> \brief Update gas concentrations, volumes, and activity coefficients
                    if (associated(this%chemistry%waters(this%chemistry%tar_wat_indices(i))%gas_chemistry)) then !> Check if gas chemistry is associated
                        call this%chemistry%waters(this%chemistry%tar_wat_indices(&
                            i))%gas_chemistry%compute_conc_gases_iter(&
                            Delta_t,&
                            this%chemistry%waters(i)%volume,[&
                            this%chemistry%waters(this%chemistry%tar_wat_indices(i))%re_mean,&
                            this%chemistry%waters(this%chemistry%tar_wat_indices(i))%rk]) !> Update gas concentrations from equilibrium and kinetic reactions [M/L³]
                        call this%chemistry%waters(this%chemistry%tar_wat_indices(&
                            i))%gas_chemistry%compute_vol_gas_species_conc() !> Compute gas volume from concentrations using ideal gas law [L³]
                        call this%chemistry%waters(this%chemistry%tar_wat_indices(&
                            i))%gas_chemistry%compute_log_act_coeffs_gases() !> Compute logarithm of gas activity coefficients (usually log(γ)=0 for ideal gases) [-]
                    end if
                !> Accumulate aqueous reaction amounts
                    this%chemistry%waters(this%chemistry%tar_wat_indices(i))%Rk_accum=&
                        this%chemistry%waters(this%chemistry%tar_wat_indices(i))%Rk_accum+&
                        this%chemistry%waters(this%chemistry%tar_wat_indices(i))%Rk
                !> Accumulate solid reaction amounts
                    this%chemistry%waters(this%chemistry%tar_wat_indices(i))%solid_chemistry%Rk_accum=&
                        this%chemistry%waters(this%chemistry%tar_wat_indices(i))%solid_chemistry%Rk_accum+&
                        this%chemistry%waters(this%chemistry%tar_wat_indices(i))%solid_chemistry%Rk
                    deallocate(c_hat) !> Deallocate mixed concentration array
                    
!> \subsection Output Writing
!> \brief Write concentrations to output file for selected waters and times
                    if ((k==this%chemistry%chem_out_options%time_steps(kk)) .and. (ind_non_can_vec(i)==&
                        this%chemistry%chem_out_options%ind_waters(ii))) then !> Check if current time and water are scheduled for output
                        write(unit,"(I10,*(ES15.5))") this%chemistry%chem_out_options%ind_waters(ii), (conc_nc(j), j=1,n_v) !> Write water index and all variable activity species concentrations
                        if (ii<this%chemistry%chem_out_options%num_waters) then !> Check if more target waters need output at this time
                            ii=ii+1 !> Increment target water counter [-]
                        else if (kk<this%chemistry%chem_out_options%num_time_steps) then !> Check if more output times remain
                            kk=kk+1 !> Increment output time step counter [-]
                            ii=1 !> Reset target water counter to first water [-]
                        else !> All outputs complete
                            exit !> Exit reactive mixing loop early
                        end if
                    end if
                    deallocate(conc_nc) !> Deallocate concentration array
            end do !> End loop over domain target waters
            
!**************************************************************************************************
!> \name Time Step Finalization
!> \brief Update state for next time step and introduce new Lagrangian particles
!**************************************************************************************************
            waters_old_old=waters_old !> Shift time history: old_old ← old
            waters_old=this%chemistry%waters !> Shift time history: old ← current
            call this%introduce_particle(k) !> Introduce new Lagrangian particles at boundaries for time step k+1
        end do !> End main time loop
        
!> Commented out redundant assignment (target waters already updated in-place)
    !this%chemistry%waters=this%chemistry%waters
    close(unit) !> Close output file
end select !> End type selection
end subroutine !> End solve_RT_ideal_lump_Lagr_trans_flux_1D