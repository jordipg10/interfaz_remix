!> @file solve_reactive_mixing_lump.f90
!> @brief Solves reactive mixing problem using lumping optimization for mixing ratios and reaction rates
!> @details Computes concentrations, activities, activity coefficients, reaction rates and volumetric 
!> fractions of minerals (if present) using a lumping technique to reduce computational cost.
!>
!> @author Jordi Petchamé-Guerrero
!> @date November 2025

subroutine solve_reactive_mixing_lump(this,dir,root,mixing_ratios_conc,mixing_ratios_R,mix_conc_indices,time_discr,&
    int_method_chem_reacts)
    !> Import main chemistry class
    use chemistry_m, only: chemistry_c
    !> Import aqueous chemistry classes, arrays, and reactive mixing subroutines
    use aqueous_chemistry_m, only: aqueous_chemistry_c, mixing_iter_comp, reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2, &
        reactive_mixing_iter_EI_kin_anal_ideal_opt2, reactive_mixing_iter_EI_eq_anal_ideal, &
        mixing_iter_comp_ideal, compute_c_mix
    !> Import time discretisation class
    use time_discr_m, only: time_discr_c
    use arrays_m, only: real_array_c, int_array_c
    !> Import vector norm function
    use vectors_m, only: inf_norm_vec_real
    implicit none !< Enforce explicit variable declarations
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: root
    class(real_array_c), intent(in) :: mixing_ratios_conc !> mixing ratios concentration matrix
    class(real_array_c), intent(in) :: mixing_ratios_R !> mixing ratios reaction rates matrix
    class(int_array_c), intent(in) :: mix_conc_indices !> matrix that contains indices of target waters that mix with each target water
    class(time_discr_c), intent(in) :: time_discr !> time discretisation object (used to solve transport)
    integer(kind=4), intent(in) :: int_method_chem_reacts !> integration method for chemical reactions
!> @par Local Variables:
    !> @var i Loop counter for target waters
    integer(kind=4) :: i !< Target water loop counter [-]
    !> @var j Loop counter for species and mixing waters
    integer(kind=4) :: j !< Species/mixing water loop counter [-]
    !> @var l Loop counter for reactive zones
    integer(kind=4) :: l !< Reactive zone loop counter [-]
    !> @var k Counter for main time loop (transport time steps)
    integer(kind=4) :: k !< Time step counter [1 to Num_time]
    !> @var kk Counter for output time steps (subset of k)
    integer(kind=4) :: kk !< Output time step counter [-]
    !> @var ii Counter for output target waters (subset of all target waters)
    integer(kind=4) :: ii !< Output target water counter [-]
    !> @var num_tar_wat Number of target waters in the system
    integer(kind=4) :: num_tar_wat !< Total number of target waters [-]
    !> @var num_tar_sol Number of target solids in the system
    integer(kind=4) :: num_tar_sol !< Total number of target solids [-]
    !> @var n_p Number of primary species (components)
    integer(kind=4) :: n_p !< Number of primary species [-]
    !> @var n_v Number of variable activity species (aqueous + solid)
    integer(kind=4) :: n_v !< Total number of variable activity species [-]
    !> @var n_v_aq Number of aqueous variable activity species
    integer(kind=4) :: n_v_aq !< Number of aqueous variable activity species [-]
    !> @var n_v_aq_2 Number of aqueous secondary variable activity species
    integer(kind=4) :: n_v_aq_2 !< Number of secondary aqueous variable activity species [-]
    !> @var unit Output file unit number
    integer(kind=4) :: unit !< File unit number for output file [-]
    !> @var tar_gas_indices Indices of target gases in each reactive zone
    integer(kind=4), allocatable :: tar_gas_indices(:) !< Target gas indices [-]
    !> @var tar_sol_indices Indices of target solids in each reactive zone
    integer(kind=4), allocatable :: tar_sol_indices(:) !< Target solid indices [-]
    !> @var tar_wat_indices Indices of target waters in each reactive zone
    integer(kind=4), allocatable :: tar_wat_indices(:) !< Target water indices [-]
    !> @var perm Permutation vector for aqueous concentrations
    integer(kind=4), allocatable :: perm(:) !< Permutation array for reordering species [-]
    !> @var time Current simulation time
    REAL(KIND=8) :: time !< Current time in simulation [s]
    !> @var Delta_t Time step size
    REAL(KIND=8) :: Delta_t !< Time step size Delta_t [s]
    !> @var theta Time weighting factor (0=EE, 1=IE, 0.5=CN)
    real(kind=8) :: theta !< Time integration weighting factor [-]
    !> @var c_mix Variable activity species concentrations after mixing
    REAL(KIND=8), allocatable :: c_mix(:) !< Mixed concentrations before reaction [mol/L or mol/kg]
    !> @var conc_nc Concentrations of variable activity species
    REAL(KIND=8), allocatable :: conc_nc(:) !< Variable activity species concentrations [mol/L or mol/kg]
    !> @var conc_comp Concentrations of components
    REAL(KIND=8), allocatable :: conc_comp(:) !< Components concentrations [mol/L or mol/kg]
    !> @var target_waters_new Target waters at time step k+1
    type(aqueous_chemistry_c), allocatable :: target_waters_new(:) !< Target waters at new time level [-]
    !> @var target_waters_old Target waters at time step k
    type(aqueous_chemistry_c), allocatable :: target_waters_old(:) !< Target waters at current time level [-]
    !> @var target_waters_old_old Target waters at time step k-1
    type(aqueous_chemistry_c), allocatable :: target_waters_old_old(:) !< Target waters at previous time level [-]
    integer(kind=4) :: j_mix !< loop index for mixing waters
    integer(kind=4) :: num_mix_loc !< number of mixing waters
    real(kind=8), allocatable :: conc_old_mix(:,:) !< conc_old from each mixing water
    integer(kind=4), allocatable :: ind_aq_sp_mix(:,:) !< indices_aq_species from each mixing water

!> @par Procedure Pointers:
    !> @var p_solver Procedure pointer for reactive mixing subroutine
    procedure(mixing_iter_comp), pointer :: p_solver=>null() !< Solver procedure pointer [-]
!> @par Initialization Section:
!> Initialize target waters for time stepping algorithm
    !> Copy initial target waters to old state (time step k)
    target_waters_old=this%waters !< State at time k [mol/L]
    !> Copy to old-old state (time step k-1) for multi-step methods
    target_waters_old_old=target_waters_old !< State at time k-1 [mol/L]
    !> Initialize new state (time step k+1) as copy of current state
    target_waters_new=target_waters_old !< State at time k+1 (to be computed) [mol/L]

    !> Initialize simulation time to zero
    time=0d0 !< Current simulation time [s]
    !> Initialize output time step counter (starts at 2 since 1 is initial condition)
    kk=2 !< Next output time step index [-]
    !> Initialize output target water counter
    ii=1 !< Next output target water index [-]
    !> Set last output time step to final simulation time
    this%chem_out_options%time_steps(this%chem_out_options%num_time_steps)=time_discr%Num_time !< Final time step for output [-]
            
    !> Set file unit number for output
    unit=7 !< File unit number for main output file [-]
    !> Open formatted output file for time-series concentration data
    open(unit,file=root//'.out',form="formatted",access="sequential",status="unknown") !< Open output file
!> @par Solver Selection Based on Reaction Types:
!> Select reactive mixing subroutine depending on the nature of the chemical system
!> and the methods to compute Jacobians and integrate in time.
!> Three main cases: (1) equilibrium + kinetics, (2) equilibrium only, (3) kinetics only
    if (this%chem_syst%num_kin_reacts>0 .and. this%chem_syst%speciation_alg%num_eq_reactions>0) then !< Case 1: Both equilibrium and kinetic reactions
        if (int_method_chem_reacts==1) then !< Explicit Euler (EE) method for chemical reactions
            !> Set time weighting factor theta=0 for fully explicit treatment
            theta=0d0 !< theta=0 means explicit Euler [-]
        else if (this%Jac_opt==1) then !< Use analytical Jacobian for Newton-Raphson iterations
            if (int_method_chem_reacts==2) then !< Implicit Euler (IE)
                !> Set time weighting factor theta=1 for fully implicit treatment
                theta=1d0 !< theta=1 means fully implicit Euler [-]
            else if (int_method_chem_reacts==3) then !< Crank-Nicolson (CN)
                !> Set time weighting factor theta=0.5 for Crank-Nicolson (2nd order accurate)
                theta=5d-1 !< theta=0.5 means Crank-Nicolson [-]
            else
                !> Invalid integration method specified
                error stop "Integration method for chemical reactions not implemented yet"
            end if
            !> Check kinetic reaction rate averaging option
            if (this%rk_avg_opt==1) then !< Option 1: not implemented
                error stop "rk average option 1 not implemented yet"
            else if (this%rk_avg_opt==2) then !< Option 2: not implemented
                error stop "rk average option 2 not implemented yet"
            else
                !> Invalid averaging option specified
                error stop "rk average option not implemented yet"
            end if
        else
            !> Invalid integration method specified
            error stop "Integration method for chemical reactions not implemented yet"
        end if
    else if (this%chem_syst%speciation_alg%num_eq_reactions>0) then !< Case 2: Only equilibrium reactions (no kinetics)
        !> Assign procedure pointer to equilibrium-only solver with lumping
        p_solver=>mixing_iter_comp !< Equilibrium-only solver with lumping optimization [-]
    else !< Case 3: Only kinetic reactions (no equilibrium)
        if (int_method_chem_reacts==1) then !< Explicit Euler for kinetics only
            !> Set time weighting factor for explicit treatment
            theta=0d0 !< theta=0 means explicit Euler [-]
        else if (this%Jac_opt==1) then !< Use analytical Jacobian for kinetics-only case
            if (int_method_chem_reacts==2) then !< Implicit Euler
                !> Set time weighting factor for fully implicit treatment
                theta=1d0 !< theta=1 means fully implicit Euler [-]
            else if (int_method_chem_reacts==3) then !< Crank-Nicolson
                !> Set time weighting factor for Crank-Nicolson
                theta=5d-1 !< theta=0.5 means Crank-Nicolson [-]
            else
                !> Invalid integration method for kinetics-only case
                error stop "Integration method for chemical reactions not implemented yet"
            end if
        else
            !> Invalid Jacobian option specified
            error stop "Integration method for chemical reactions not implemented yet"
        end if
    end if     
!> @par Main Time Loop:
!> Iterate over all transport time steps to solve coupled mixing-reaction problem
!> At each time step: (1) mix waters, (2) react, (3) update minerals/gases, (4) output results
        !> Main time loop: iterate from k=1 to final time step
        do k=1,time_discr%Num_time !< Loop over all transport time steps [-]
            !> Get time step size for current iteration (may vary for adaptive stepping)
            Delta_t=time_discr%get_Delta_t(k) !< Time step size Delta_t [s]
            !> Advance simulation time by Delta_t
            time=time+Delta_t !< Current time = previous time + Delta_t [s]
            !> Check if current time step matches an output time step
            if (k==this%chem_out_options%time_steps(kk)) then !< Check if this is an output time step
                !> Write current time in scientific notation (ES15.5 format)
                write(unit,"(2x,'t = ',*(ES15.5))") time !< Output: "  t = 1.23456E+02" [s]
                !> Write column headers with aqueous species names (15 characters each)
                write(unit,"(A10,*(A15))") 'Water', (this%chem_syst%aq_phase%aq_species(this%chem_out_options%ind_aq_species(j))%name, &
                j=1,this%chem_out_options%num_aq_species) !< Species names from output options
            end if
        !> @par Inner Loop Over Domain Waters:
        !> For each domain water, solve the mixing-reaction problem
        !> This is the core computational loop of the reactive transport solver
            !> Inner loop: iterate over domain target waters only
            do i=1,this%num_target_waters !< Loop over domain target waters [-]
                !> @par System Size Extraction:
                !> Extract dimensions of chemical system for current target water
                !> These determine array sizes for reactive mixing solver
                
                !> Get number of primary species (aqueous components)
                n_p=this%waters(this%tar_wat_indices(i))%solid_chemistry%reactive_zone%speciation_alg%num_prim_species !< n_p [-]
                !> Get total number of variable activity species (aqueous + solid)
                n_v=this%waters(this%tar_wat_indices(i)&
                    )%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species !< n_v [-]
                !> Get number of aqueous variable activity species
                n_v_aq=this%waters(this%tar_wat_indices(i)&
                    )%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species !< n_v_aq [-]
                !> Get number of secondary aqueous variable activity species
                n_v_aq_2=this%waters(this%tar_wat_indices(i)&
                    )%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species !< n_v_aq_2 [-]
                !> @par Lumping Criterion Check:
                !> Check if this target water can be "lumped" (skipped) based on mixing ratios
                !> If mixing_ratio[1] ≈ 1 and all other ratios ≈ 0, then no mixing occurs
                if (abs(mixing_ratios_conc%cols(i)%col_1(1)-1d0)<&
                this%waters(this%tar_wat_indices(i))%solid_chemistry%reactive_zone%CV_params%abs_tol .and. &
                    inf_norm_vec_real(mixing_ratios_conc%cols(i)%col_1(2:mixing_ratios_conc%cols(i)%dim))<&
                this%waters(this%tar_wat_indices(i))%solid_chemistry%reactive_zone%CV_params%abs_tol) then !< Check lumping criterion
                    !> Skip this water: no significant mixing, so no computation needed
                    continue !< Move to next target water [-]
                else !< Significant mixing occurs, need to solve reactive mixing
                    !> @par Solver Selection for Cation Exchange:
                    !> Check if cation exchange reactions are present in reactive zone
                    !if (this%waters(this%tar_wat_indices(i))%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl&
                    !    >0) then !< Check for surface complexation (cation exchange)
                    !    !> Assign procedure pointer to cation exchange solver with lumping
                    !    !p_solver=>mixing_iter_comp_exch_lump !< Equilibrium solver with cation exchange and lumping [-]
                    !end if
                    !> @par Array Allocation:
                    !> Allocate working array for concentrations of variable activity species
                    allocate(conc_nc(n_v)) !< Allocate conc_nc array [n_v] [mol/L or mol/kg]
                    allocate(c_mix(n_v)) !< Allocate c_mix array [n_v] [mol/L or mol/kg]
                    allocate(conc_comp(n_p)) !< Allocate conc_comp array [num_comp_species] [mol/L or mol/kg]
                    !> @par Index Validation:
                    !> Verify that target water indices match expected order
                    !> This ensures mix_conc_indices and tar_wat_indices are consistent
                    if (mix_conc_indices%cols(i)%col_1(1)/=this%tar_wat_indices(i)) then !< Check index consistency
                        !> Error: inconsistent target water ordering
                        print *, mix_conc_indices%cols(i)%col_1(1), this%tar_wat_indices(i) !< Print mismatched indices
                        error stop "Target waters not in the right order" !< Halt execution
                    end if
                !> @par Advective Mixing Computation:
                !> Compute mixed concentrations c_mix from upstream waters
                !> Initialize with current target water concentrations
                    !c_mix=target_waters_old(this%tar_wat_indices(i))%get_conc_nc() !< Get variable activity species concentrations [mol/L or mol/kg]
                    num_mix_loc=size(mix_conc_indices%cols(i)%col_1)-1
                    allocate(conc_old_mix(size(target_waters_old(mix_conc_indices%cols(i)%col_1(1))%conc_old), num_mix_loc))
                    allocate(ind_aq_sp_mix(size(target_waters_old(mix_conc_indices%cols(i)%col_1(1))%indices_aq_species), num_mix_loc))
                    do j_mix=1,num_mix_loc
                        conc_old_mix(:,j_mix)=target_waters_old(mix_conc_indices%cols(i)%col_1(j_mix+1))%conc_old
                        ind_aq_sp_mix(:,j_mix)=target_waters_old(mix_conc_indices%cols(i)%col_1(j_mix+1))%indices_aq_species
                    end do
                    call compute_c_mix(target_waters_old(mix_conc_indices%cols(i)%col_1(1)),&
                        conc_old_mix,ind_aq_sp_mix,&
                        mixing_ratios_conc%cols(i)%col_1,c_mix)
                    deallocate(conc_old_mix,ind_aq_sp_mix)
                !> @par Reactive Mixing Solver Call:
                !> Call procedure pointer to selected solver
                !> This solves the coupled system of reaction equations after mixing
                    call p_solver(target_waters_new(this%tar_wat_indices(i)),& !< [inout] Target water to solve (updated)
                        target_waters_old_old(this%tar_wat_indices(i))%get_c1(),& !< [in] Primary species at k-1 [mol/L]
                        c_mix,& !< [in] Mixed concentrations from advection
                        mixing_ratios_R%cols(i)%col_1(1),& !< [in] Mixing ratios for reaction rates
                        Delta_t,& !< [in] Time step size
                        theta,& !< [in] Time weighting factor [-]
                        conc_nc,& !< [out] Variable activity species concentrations at k+1 [mol/L or mol/kg]
                        conc_comp) !< [out] Components concentrations at k+1 [mol/L or mol/kg]
                        
                !> @par Equilibrium Reaction Rate Computation:
                !> Compute equilibrium reaction rates Re from mass balance equation
                !> Re is calculated post hoc from the change in primary species concentrations
                    call target_waters_new(this%tar_wat_indices(i))%compute_Re(&
                        c_mix(1:n_v),& !< [in] Secondary species after mixing
                        Delta_t,& !< [in] Time step size
                        mixing_ratios_R%cols(i)%col_1(1)) !> Compute equilibrium reaction rates using lumped WMA
            !> @par Mineral Chemistry Update:
            !> Update solid chemistry state variables (mineral concentrations and volumetric fractions)
            !> This applies kinetic mineral reactions to compute new mineral masses
                    !> Check if minerals are present in this reactive zone
                    if (associated(target_waters_new(this%tar_wat_indices(i))%solid_chemistry%mineral_zone)) then !< Check for mineral zone
                    !> Compute mass balance for minerals (dissolution/precipitation)
                    !> Updates mineral volumetric fractions based on kinetic reaction rates
                        call target_waters_new(this%tar_wat_indices(i))%solid_chemistry%compute_mass_bal_mins(Delta_t) !< [in] Time step [s]
                    !> Compute updated mineral concentrations for current iteration
                    !> Ensures minerals remain non-negative and updates volumetric properties
                        call target_waters_new(this%tar_wat_indices(i))%solid_chemistry%compute_conc_minerals_iter(Delta_t) !< [in] Time step [s]
                    end if !< End mineral zone check
                    
                !> @par Gas Chemistry Update:
                !> Update gas chemistry state variables (gas phase concentrations and partial pressures)
                !> This applies equilibrium/kinetic gas reactions
                    !> Check if gas chemistry object is associated with this target water
                    if (associated(target_waters_new(this%tar_wat_indices(i))%gas_chemistry)) then !< Check for gas chemistry
                    !> Compute updated gas concentrations for current iteration
                    !> Uses equilibrium constants and mean reaction rates
                        call target_waters_new(this%tar_wat_indices(i))%gas_chemistry%compute_conc_gases_iter(&
                            Delta_t,& !< [in] Time step size [s]
                            target_waters_new(i)%volume,& !< [in] Water volume [L]
                            [target_waters_new(i)%re_mean,target_waters_new(i)%rk]) !< [in] Mean reaction rates [mol/L/s]
                    !> Compute volumetric gas concentrations [mol/m^3] from molar amounts
                        call target_waters_new(this%tar_wat_indices(i))%gas_chemistry%compute_vol_gas_species_conc() !< Compute volume [m^3]
                    !> Compute logarithms of gas phase activity coefficients (for fugacity corrections)
                        call target_waters_new(this%tar_wat_indices(i))%gas_chemistry%compute_log_act_coeffs_gases() !< Compute log(gamma) [-]
                    end if !< End gas chemistry check
                !> @par Conditional Output Writing:
                !> Write concentrations to output file if this is an output time step and target water
                !> Output format: target_water_index followed by selected aqueous species concentrations
                    if (k==this%chem_out_options%time_steps(kk) .and. i==this%chem_out_options%ind_waters(ii)) then !< Check output criteria
                        !> Write target water index and aqueous species concentrations
                        !> Format: I10 for index, ES15.5 for each concentration (scientific notation)
                        write(unit,"(I10,*(ES15.5))") i, (conc_nc(j), j=1,this%chem_out_options%num_aq_species) !< Write data to file
                        
                        !> @par Output Counter Updates:
                        !> Increment output counters after writing
                        if (ii<this%chem_out_options%num_waters) then !< More target waters at this time step
                            ii=ii+1 !< Move to next output target water
                        else if (kk<this%chem_out_options%num_time_steps) then !< More time steps remaining
                            kk=kk+1 !< Move to next output time step
                            ii=1 !< Reset target water counter to first
                        else !< All output complete
                            exit !< Exit domain water loop (all output written)
                        end if
                    end if !< End output check
                    
                !> @par Memory Cleanup:
                !> Deallocate temporary arrays for this target water
                    deallocate(c_mix,conc_nc) !< Free c_mix and conc_nc memory
                end if !< End lumping criterion check
            end do !< End inner loop over domain target waters
            
        !> @par Time Step Finalization:
        !> Shift time levels for multi-step time integration methods
        !> This prepares state variables for the next time step k+1
            !> Copy time level k to k-1 (for Crank-Nicolson and other multi-step methods)
            target_waters_old_old=target_waters_old !< c^(k-1) ← c^k
            !> Copy time level k+1 to k (advance solution to next time step)
            target_waters_old=target_waters_new !< c^k ← c^(k+1)
        end do !< End main time loop over all transport time steps
        
!> @par Final State Assignment:
!> Copy final computed target waters back to chemistry object
!> This updates the persistent state with the solution at final time
    this%waters=target_waters_new !< Transfer final solution to chemistry object
    
!> @par Output File Cleanup:
!> Close output file after all time steps are written
    close(unit) !< Close file unit 7 (root.out)
end subroutine !< End of solve_reactive_mixing_lump